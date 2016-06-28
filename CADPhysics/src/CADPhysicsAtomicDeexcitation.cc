// CADPhysicsAtomicDeexcitation.cc
//
// A general note on shell indices versus shell IDs. The shell IDs are those used in the
// Lawrence Livermore Evaluated Atomic Data Libraries, EADL, to designate subshells of atoms.
// These IDs are described in table VI, on page 7, of the document
// 'ENDL Type Formats for the Livermore Evaluated Atomic Data Library, EADL' which can be found 
// (a.o.) at http://www.nea.fr/html/dbdata/data/nds_eval_eadl.pdf.
// Shell indices, on the other hand, are an internal Geant4 designation of shells per element. The
// shell index starts at zero, and is for each element just a counter over all shells that appear in the file
// G4EMLOW/fluor/binding.dat . In general, this leads to the following designations:
//
// Name  ID  index
// K     1   0
// L1    3   1
// L2    5   2
// L3    6   3
// M1    8   4
// M2    10  5
// (...)
// M5    14  8
// N1    16  9
// etcetera.
// 
// The highest valid number for the index of course depends on the element (Z number). For instance, for hydrogen the
// highest index is 0, while for silicon (Z=14) it is 6 (shells defined up to the M3 shell, ID=11) and for gold (Z=79) it's 
// 21 (highest shell is P1, ID=41).

#include "CADPhysicsAtomicDeexcitation.hh"
#include "Randomize.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4FluoTransition.hh"

CADPhysicsAtomicDeexcitation::CADPhysicsAtomicDeexcitation():
minGammaEnergy(100.*eV),
minElectronEnergy(100.*eV),
fAuger(false),
transitionManager(G4AtomicTransitionManager::Instance())
{}

CADPhysicsAtomicDeexcitation::~CADPhysicsAtomicDeexcitation()
{}

std::vector<G4DynamicParticle*>* CADPhysicsAtomicDeexcitation::GenerateParticles(G4int Z,G4int givenShellId)
{
	// Some initialization
	VacancyEnergies.clear();
	ShellIdstobeprocessed.clear();
	std::vector<G4DynamicParticle*>* vectorOfParticles = new std::vector<G4DynamicParticle*>;
	ShellIdstobeprocessed.push_back(givenShellId);
	G4DynamicParticle* aParticle;
	G4int currentShellId = 0;// ID of the current shell under investigation (in which a vacancy is present)
	G4int provShellId = 0;// ID of a higher shell, i.e. the 'starting point' for the transition

	do {// Repeat the Fluorescence/Auger process as long as there are any shell vacancies to be processed
		currentShellId = ShellIdstobeprocessed.back();
		ShellIdstobeprocessed.pop_back();
		provShellId = SelectTypeOfTransition(Z, currentShellId);// First find out if a radiative transition
		// should take place and if so, determine the ID of the upper shell ('starting point' of the transition)
		if  (provShellId >0)// Radiative transition between the already-known shells
		{
			aParticle = GenerateFluorescence(Z,currentShellId,provShellId);  
		}
		else if (provShellId == -1)// Non-radiative transition
		{
			aParticle = GenerateAuger(Z, currentShellId);
		}
		else
		{
			G4Exception("CADPhysicsAtomicDeexcitation::GenerateParticles()","errorCADPhysicsAtomicDeexcitation01",FatalException,"CADPhysicsAtomicDeexcitation: starting shell uncorrect: check it");
		}
		if (aParticle != 0) {vectorOfParticles->push_back(aParticle);}
	}
	while (ShellIdstobeprocessed.size()>0); 

	return vectorOfParticles;
}

G4int CADPhysicsAtomicDeexcitation::SelectTypeOfTransition(G4int Z, G4int shellId)
{
	if (shellId <=0 ) G4Exception("CADPhysicsAtomicDeexcitation::SelectTypeOfTransition()","errorCADPhysicsAtomicDeexcitation02",FatalException,"CADPhysicsAtomicDeexcitation: zero or negative shellId");

	G4bool fluoTransitionFoundFlag = false;

	G4int provShellId = -1;// ID of the upper shell
	G4int shellNum = 0;// INDEX corresponding to shellId (lower shell) in a list of shells that can serve as
	// lower shells in a fluorescent transition
	G4int maxNumOfShells = transitionManager->NumberOfReachableShells(Z);  

	const G4FluoTransition* refShell = transitionManager->ReachableShell(Z,maxNumOfShells-1);
	// refShell corresponds to the highest shell that can still serve as a final state for a fluorescent transition.

	if ( shellId <= refShell->FinalShellId())// If shellId is larger than the ID of refShell,
		// then it means that there can be NO fluorescent transition that has shellId as the lower shell.
		// Hence, we can save us the trouble looking for one.
	{
		// Now find the index of shellId in the list of shells that can get a fluorescent transition.
		while (shellId != transitionManager->ReachableShell(Z,shellNum)->FinalShellId())
		{
			if(shellNum==maxNumOfShells-1) break;
			shellNum++;
		}
		G4int transProb = 1;

		G4double partialProb = G4UniformRand();      
		G4double partSum = 0;
		const G4FluoTransition* aShell = transitionManager->ReachableShell(Z,shellNum);    
		G4int trSize =  (aShell->TransitionProbabilities()).size();

		// Loop over the shells wich can provide an electron for a radiative transition towards shellId:
		// in every run of the loop the partial sum of probabilities of the first transProb shells
		// is calculated and compared with a random number [0,1].
		// If the partial sum is greater, the shell whose index is transProb
		// is chosen as the starting shell for a radiative transition
		// and its ID is returned.
		// Note that the sum of all radiative probabilities can be lower than unity, due to competition 
		// with Auger transitions. Hence, the loop may be terminated without finding
		// a radiative transition, and in that case -1 is returned and control is handed over to
		// the GenerateAuger routine.
		while(transProb < trSize){

			partSum += aShell->TransitionProbability(transProb);

			if(partialProb <= partSum)
			{
				provShellId = aShell->OriginatingShellId(transProb);
				fluoTransitionFoundFlag = true;
				break;
			}
			transProb++;
		}
	}
	// If a radiative transition was found, the ID of the 'upper shell' is returned. If not, an Auger transition may take place.
	return provShellId;
}

G4DynamicParticle* CADPhysicsAtomicDeexcitation::GenerateFluorescence(G4int Z, 
																	  G4int shellId,
																	  G4int provShellId )
{ 
	// Isotropic angular distribution for the photon
	G4double newcosTh = 1.-2.*G4UniformRand();
	G4double  newsinTh = std::sqrt(1.-newcosTh*newcosTh);
	G4double newPhi = twopi*G4UniformRand();
	G4double xDir =  newsinTh*std::sin(newPhi);
	G4double yDir = newsinTh*std::cos(newPhi);
	G4double zDir = newcosTh;
	G4ThreeVector newGammaDirection(xDir,yDir,zDir);

	// Again, find the index of shellId
	G4int shellNum = 0;
	G4int maxNumOfShells = transitionManager->NumberOfReachableShells(Z);
	while (shellId != transitionManager->
		ReachableShell(Z,shellNum)->FinalShellId())
	{
		if(shellNum == maxNumOfShells-1)
		{
			break;
		}
		shellNum++;
	}

	// Find the index of provShellId
	const G4FluoTransition* aShell = transitionManager->ReachableShell(Z,shellNum);     
	size_t transitionSize = aShell->OriginatingShellIds().size();
	size_t index = 0;
	while (provShellId != aShell->OriginatingShellId(index))
	{
		if(index ==  transitionSize-1)
		{
			break;
		}
		index++;
	}

	// Energy of the gamma leaving provShellId for shellId
	G4double transitionEnergy = aShell->TransitionEnergy(index);

	// This is the shell where the new vacancy is: it is the same
	// shell where the electron came from. It is stored to the
	// vector of shells to be processed.
	ShellIdstobeprocessed.push_back(aShell->OriginatingShellId(index));

	// Return the photon as a new particle
	G4DynamicParticle* newPart = new G4DynamicParticle(G4Gamma::Gamma(), 
		newGammaDirection,
		transitionEnergy);
	return newPart;
}

G4DynamicParticle* CADPhysicsAtomicDeexcitation::GenerateAuger(G4int Z, G4int shellId)
{
	// If Auger production was disabled, then simply push the binding energy of the 
	// current shell to the VacancyEnergies vector and return
	if(!fAuger) {
		PushEnergy(Z,shellId);
		return 0;
	}

	if (shellId <=0 ) G4Exception("CADPhysicsAtomicDeexcitation::GenerateAuger()","errorCADPhysicsAtomicDeexcitation03",FatalException,"CADPhysicsAtomicDeexcitation: zero or negative shellId");

	G4int maxNumOfShells = transitionManager->NumberOfReachableAugerShells(Z);  
	const G4AugerTransition* refAugerTransition = 
		transitionManager->ReachableAugerShell(Z,maxNumOfShells-1);
	// refAugerTransition corresponds to the highest shell that can still serve as a final state for an Auger transition.

	G4int shellNum = 0;// index corresponding to shell ID
	if ( shellId <= refAugerTransition->FinalShellId() )// If shellId is larger than the ID of refAugerTransition,
		// then it means that there can be NO Auger transition that has shellId as the lower shell.
		// Hence, we can save us the trouble looking for one.
	{
		// Now find the index of shellId in the list of shells that can get an Auger transition.
		G4int pippo = transitionManager->ReachableAugerShell(Z,shellNum)->FinalShellId();
		if (shellId  != pippo )
		{
			do { 
				shellNum++;
				if(shellNum == maxNumOfShells)// No Auger transitions for shellId!
				{
					PushEnergy(Z,shellId);// Add this shell's energy to the VacancyEnergies vector
					return 0;
				}
			}
			while (shellId != (transitionManager->ReachableAugerShell(Z,shellNum)->FinalShellId()) ) ;
		}

		// Find the Auger transition
		G4int transitionLoopShellIndex = 0;      
		G4double partSum = 0;
		const G4AugerTransition* anAugerTransition = 
			transitionManager->ReachableAugerShell(Z,shellNum);// anAugerTransition contains ALL possible
		// Auger transitions that have shellId as the lower shell

		// First sum up the transition probabilities over all possible Auger transitions that have shellId as the lower shell
		G4int transitionSize = 
			(anAugerTransition->TransitionOriginatingShellIds())->size();
		while (transitionLoopShellIndex < transitionSize)
		{
			std::vector<G4int>::const_iterator pos = 
				anAugerTransition->TransitionOriginatingShellIds()->begin();

			G4int transitionLoopShellId = *(pos+transitionLoopShellIndex);
			G4int numberOfPossibleAuger = 
				(anAugerTransition->AugerTransitionProbabilities(transitionLoopShellId))->size();
			G4int augerIndex = 0;

			if (augerIndex < numberOfPossibleAuger) {
				do 	{
					G4double thisProb = anAugerTransition->AugerTransitionProbability(augerIndex, 
						transitionLoopShellId);
					partSum += thisProb;
					augerIndex++;
				} while (augerIndex < numberOfPossibleAuger);
			}
			transitionLoopShellIndex++;
		}
		G4double totalVacancyAugerProbability = partSum;

		// Next, select a specific Auger process based on partial probabilities.
		// transitionRandomShellId corresponds to the shell from which an electron goes to shellId
		// augerId corresponds to the shell from which an electron is emitted (into vacuum)
		G4int transitionRandomShellIndex = 0;
		G4int transitionRandomShellId = 1;
		G4int augerIndex = 0;
		G4int augerId = 1;
		partSum = 0; 
		G4double partialProb = G4UniformRand();


		G4int numberOfPossibleAuger = 0;
		G4bool foundFlag = false;

		while (transitionRandomShellIndex < transitionSize)// Loop over the transitionRandomShellIds
		{
			std::vector<G4int>::const_iterator pos = 
				anAugerTransition->TransitionOriginatingShellIds()->begin();
			transitionRandomShellId = *(pos+transitionRandomShellIndex);// Find transitionRandomShellId
			//corresponding to transitionRandomShellIndex

			augerIndex = 0;
			numberOfPossibleAuger = (anAugerTransition-> 
				AugerTransitionProbabilities(transitionRandomShellId))->size();
			while (augerIndex < numberOfPossibleAuger)// Loop over the augerIds for the given transitionRandomShellId
			{
				G4double thisProb =anAugerTransition->AugerTransitionProbability(augerIndex, 
					transitionRandomShellId);// Partial probability of this specific Auger process
				partSum += thisProb;
				// EB: original code: if (partSum >= (partialProb/totalVacancyAugerProbability) )// This specific transition is selected
                                // EB: totalVacancyAugerProbability is the total probability of generating an Auger 
                                // EB: partSum should, in case we reach the end of the loop equal totalVacancyAugerProbability i.e. 100% of that probability
                                // EB: partialProb is the random number we have drawn to determine which of these transitions is selected.
                                // EB: This means that only a multiplication here makes sense. Note that in G4 9.5 this is indeed repaired with the
                                // EB: famous comment :   // was /
				if (partSum >= (partialProb*totalVacancyAugerProbability) )// This specific transition is selected
				{
					foundFlag = true;
					augerId = anAugerTransition->AugerOriginatingShellId(augerIndex,transitionRandomShellId);// Find augerId 
					// corresponding to augerIndex 
					break;
				}
				augerIndex++;
			}
			if (foundFlag) break;
			transitionRandomShellIndex++;
		}

		if (!foundFlag) {// After looping over all possible transitions, no suitable transition has been found. Note that in principle
			// this should never occur.
			PushEnergy(Z,shellId);// Simply add the binding energy corresponding to shellId to the VacancyEnergies vector
			return 0;
		}

		// Deal with the details of the newly emitted electron:
		// Isotropic angular distribution
		G4double newcosTh = 1.-2.*G4UniformRand();
		G4double newsinTh = std::sqrt(1.-newcosTh*newcosTh);
		G4double newPhi = twopi*G4UniformRand();
		G4double xDir =  newsinTh*std::sin(newPhi);
		G4double yDir = newsinTh*std::cos(newPhi);
		G4double zDir = newcosTh;
		G4ThreeVector newElectronDirection(xDir,yDir,zDir);

		// Energy of the Auger electron emitted
		G4double transitionEnergy = anAugerTransition->AugerTransitionEnergy(augerIndex, transitionRandomShellId);
    	//	G4cout << "Auger electron energy: " << transitionEnergy
		//       << " MeV for element Z=" << Z;
		//  G4cout << G4endl;

		// As a result of the Auger process, we have TWO new shells with vacancies. Both are added
		// to the vector of shells to be processed.
		ShellIdstobeprocessed.push_back(transitionRandomShellId);
		ShellIdstobeprocessed.push_back(augerId);

		G4DynamicParticle* newPart = new G4DynamicParticle(G4Electron::Electron(), 
			newElectronDirection,
			transitionEnergy);
		return newPart;
	} else {// No Auger transitions for this shell
		PushEnergy(Z,shellId);
		return 0;
	}  
}

void CADPhysicsAtomicDeexcitation::PushEnergy(G4int Z, G4int ShellId)
{
	// Loop over all shells of Z to find the shell index corresponding to ShellId.
	// Once found, immediately store the binding energy corresponding to that shell to VacancyEnergies.
	G4int maxNumOfShells = transitionManager->NumberOfShells(Z);
	G4int index = 0;
	while (index < maxNumOfShells)
	{
		if(ShellId == transitionManager->Shell(Z,index)->ShellId()) {
			VacancyEnergies.push_back(transitionManager->Shell(Z,index)->BindingEnergy());
			break;
		}
		index++;
	}
	return;
}
