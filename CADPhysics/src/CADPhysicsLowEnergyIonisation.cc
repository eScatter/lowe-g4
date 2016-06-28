// CADPhysicsLowEnergyIonisation.cc

// 2007-01-04: introducing ions as secondaries for tracing

#include "CADPhysicsLowEnergyIonisation.hh"
#include "G4eIonisationSpectrum.hh"
#include "G4eIonisationCrossSectionHandler.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4AtomicShell.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4SemiLogInterpolation.hh"
#include "G4LogLogInterpolation.hh"
#include "G4EMDataSet.hh"
#include "G4VEMDataSet.hh"
#include "G4CompositeEMDataSet.hh"
#include "G4ShellVacancy.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4ProductionCutsTable.hh"

//EK, for creation of ions
//#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4Transportation.hh"
#include "G4VProcess.hh"
#include "CADPhysicsIonKill.hh"

#include "CADPhysicsLowEnIoniMessenger.hh"

CADPhysicsLowEnergyIonisation::CADPhysicsLowEnergyIonisation(const G4String& nam)
: G4eLowEnergyLoss(nam), 
crossSectionHandler(0),
theMeanFreePath(0),
energySpectrum(0),
shellVacancy(0)
{
   MinKineticEnergy = DBL_MIN;// EK, AlongStepDoIt will kill all electrons below this energy value!
   cutForPhotons = 250.0*eV;
   cutForElectrons = 250.0*eV;
   useCutForElectrons = false;
   generateIons = false;
   verboseLevel = 0;
   electron_mass_over_amu = electron_mass_c2/amu_c2;
   messenger = new CADPhysicsLowEnIoniMessenger(this);
}


CADPhysicsLowEnergyIonisation::~CADPhysicsLowEnergyIonisation()
{
   delete crossSectionHandler;
   delete energySpectrum;
   delete theMeanFreePath;
   delete shellVacancy;
}


void CADPhysicsLowEnergyIonisation::BuildPhysicsTable(const G4ParticleDefinition& aParticleType)
{ 
   if(verboseLevel > 0) {
   G4cout << "CADPhysicsLowEnergyIonisation::BuildPhysicsTable start"
      << G4endl;
   }

   cutForDelta.clear();

   // Create and fill IonisationParameters once
   if( energySpectrum != 0 ) delete energySpectrum;
   energySpectrum = new G4eIonisationSpectrum();

   if(verboseLevel > 0) {
      G4cout << "G4VEnergySpectrum is initialized"
         << G4endl;
   }

   // Create and fill G4CrossSectionHandler once

   if ( crossSectionHandler != 0 ) delete crossSectionHandler;
   G4VDataSetAlgorithm* interpolation = new G4SemiLogInterpolation();
   G4double lowKineticEnergy  = GetLowerBoundEloss();
   G4double highKineticEnergy = GetUpperBoundEloss();
   G4int    totBin = GetNbinEloss();
   crossSectionHandler = new G4eIonisationCrossSectionHandler(energySpectrum,
      interpolation,
      lowKineticEnergy,
      highKineticEnergy,
      totBin);
   crossSectionHandler->LoadShellData("ioni/ion-ss-cs-");

   if (verboseLevel > 0) {
      G4cout << GetProcessName()
         << " is created; Cross section data: "
         << G4endl;
      crossSectionHandler->PrintData();
      G4cout << "Parameters: "
         << G4endl;
      energySpectrum->PrintData();
   }

   // Build loss table for IonisationIV

   BuildLossTable(aParticleType);

   if(verboseLevel > 0) {
      G4cout << "The loss table is built"
         << G4endl;
   }

   if (&aParticleType==G4Electron::Electron()) {

      RecorderOfElectronProcess[CounterOfElectronProcess] = (*this).theLossTable;
      CounterOfElectronProcess++;
      PrintInfoDefinition();  

   } else {

      RecorderOfPositronProcess[CounterOfPositronProcess] = (*this).theLossTable;
      CounterOfPositronProcess++;
   }

   // Build mean free path data using cut values

   if( theMeanFreePath ) delete theMeanFreePath;
   theMeanFreePath = crossSectionHandler->
      BuildMeanFreePathForMaterials(&cutForDelta);

   if(verboseLevel > 0) {
      G4cout << "The MeanFreePath table is built"
         << G4endl;
      if(verboseLevel > 1) theMeanFreePath->PrintData();
   }

   // Build common DEDX table for all ionisation processes

   BuildDEDXTable(aParticleType);

   if (verboseLevel > 0) {
      G4cout << "CADPhysicsLowEnergyIonisation::BuildPhysicsTable end"
         << G4endl;
   }
}


void CADPhysicsLowEnergyIonisation::BuildLossTable(const G4ParticleDefinition& )
{
   // Build table for energy loss due to soft brems
   // the tables are built for *MATERIALS* binning is taken from LowEnergyLoss

   G4double lowKineticEnergy  = GetLowerBoundEloss();
   G4double highKineticEnergy = GetUpperBoundEloss();
   size_t   totBin = GetNbinEloss();

   //  create table

   if (theLossTable) { 
      theLossTable->clearAndDestroy();
      delete theLossTable;
   }
   const G4ProductionCutsTable* theCoupleTable=
      G4ProductionCutsTable::GetProductionCutsTable();
   size_t numOfCouples = theCoupleTable->GetTableSize();
   theLossTable = new G4PhysicsTable(numOfCouples);

   if (shellVacancy != 0) delete shellVacancy;
   shellVacancy = new G4ShellVacancy();
   G4DataVector* ksi = 0;
   G4DataVector* energy = 0;
   size_t binForFluo = totBin/10;

   G4PhysicsLogVector* bVector = new G4PhysicsLogVector(lowKineticEnergy,
      highKineticEnergy,
      binForFluo);
   const G4AtomicTransitionManager* transitionManager = G4AtomicTransitionManager::Instance();

   // Clean up the vector of cuts

   cutForDelta.clear();

   // Loop for materials

   ionTable = G4IonTable::GetIonTable();
   //particleTable = G4ParticleTable::GetParticleTable(); // EK, used below for creation of ions
   G4ParticleDefinition* ion;//idem
   G4Element* theelement;
   G4Isotope* isotope;//idem

   for (size_t m=0; m<numOfCouples; m++) {

      // create physics vector and fill it
      G4PhysicsLogVector* aVector = new G4PhysicsLogVector(lowKineticEnergy,
         highKineticEnergy,
         totBin);

      // get material parameters needed for the energy loss calculation
      const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(m);
      const G4Material* material= couple->GetMaterial();

      G4bool isgas = (material->GetState()==kStateGas);// EK, true if the material is a gas

      // the cut cannot be below lowest limit
      G4double tCut = (*(theCoupleTable->GetEnergyCutsVector(1)))[m];
      if(tCut > highKineticEnergy) tCut = highKineticEnergy;
      cutForDelta.push_back(tCut);
      const G4ElementVector* theElementVector = material->GetElementVector();
      size_t NumberOfElements = material->GetNumberOfElements() ;
      const G4double* theAtomicNumDensityVector =
         material->GetAtomicNumDensityVector();
      if(verboseLevel > 0) {
         G4cout << "Energy loss for material # " << m
            << " tCut(keV)= " << tCut/keV
            << G4endl;
      }

      // now comes the loop for the kinetic energy values
      for (size_t i = 0; i<totBin; i++) {

         G4double lowEdgeEnergy = aVector->GetLowEdgeEnergy(i);
         G4double ionloss = 0.;

         if (isgas) {// EK, Add energy loss only for gases
            // loop for elements in the material
            for (size_t iel=0; iel<NumberOfElements; iel++ ) {

               G4int Z = (G4int)((*theElementVector)[iel]->GetZ());

               G4int nShells = transitionManager->NumberOfShells(Z);

               for (G4int n=0; n<nShells; n++) {

                  G4double e = energySpectrum->AverageEnergy(Z, 0.0, tCut,
                     lowEdgeEnergy, n);
                  G4double cs= crossSectionHandler->FindValue(Z, lowEdgeEnergy, n);
                  ionloss += e * cs * theAtomicNumDensityVector[iel];

                  if(verboseLevel > 1) {
                     G4cout << "Z= " << Z
                        << " shell= " << n
                        << " E(keV)= " << lowEdgeEnergy/keV
                        << " Eav(keV)= " << e/keV
                        << " cs= " << cs
                        << " loss= " << ionloss
                        << " rho= " << theAtomicNumDensityVector[iel]
                     << G4endl;
                  }
               }
               G4double esp = energySpectrum->Excitation(Z, lowEdgeEnergy);
               //G4cout << Z << "\t" << esp*theAtomicNumDensityVector[iel]*mm/eV << "\t";
               ionloss   += esp * theAtomicNumDensityVector[iel];

            }
         }
         // EK, the following check has been added:
         if(ionloss==0.) ionloss = DBL_MIN;
         // EK: zero ionloss leads to 'inf' range if no other processes (CADPhysicsLowEnergyBremsstrahlung) are active.
         // An 'inf' range, in turn, leads to unexpected behaviour of the AlongStepDoIt method.

         if(verboseLevel > 1) {
            G4cout << "Sum: "
               << " E(keV)= " << lowEdgeEnergy/keV
               << " loss(MeV/mm)= " << ionloss*mm/MeV
               << G4endl;
         }

         aVector->PutValue(i,ionloss);
      }
      theLossTable->insert(aVector);

      // fill data for fluorescence

      G4VDataSetAlgorithm* interp = new G4LogLogInterpolation();
      G4VEMDataSet* xsis = new G4CompositeEMDataSet(interp, 1., 1.);
      for (size_t iel=0; iel<NumberOfElements; iel++ ) {

         G4int Z = (G4int)((*theElementVector)[iel]->GetZ());
         energy = new G4DataVector();
         ksi    = new G4DataVector();

         for (size_t j = 0; j<binForFluo; j++) {

            G4double lowEdgeEnergy = bVector->GetLowEdgeEnergy(j);
            G4double cross   = 0.;
            G4double eAverage= 0.;
            G4int nShells = transitionManager->NumberOfShells(Z);

            for (G4int n=0; n<nShells; n++) {

               G4double e = energySpectrum->AverageEnergy(Z, 0.0, tCut,
                  lowEdgeEnergy, n);
               G4double pro = energySpectrum->Probability(Z, 0.0, tCut,
                  lowEdgeEnergy, n);
               G4double cs= crossSectionHandler->FindValue(Z, lowEdgeEnergy, n);
               eAverage   += e * cs * theAtomicNumDensityVector[iel];
               cross      += cs * pro * theAtomicNumDensityVector[iel];
               if(verboseLevel > 1) {
                  G4cout << "Z= " << Z
                     << " shell= " << n
                     << " E(keV)= " << lowEdgeEnergy/keV
                     << " Eav(keV)= " << e/keV
                     << " pro= " << pro
                     << " cs= " << cs
                     << G4endl;
               }//if verboselevel
            }//for n<nShells

            G4double coeff = 0.0;
            if(eAverage > 0.) {
               coeff = cross/eAverage;
               eAverage /= cross;
            }//if eAverage

            if(verboseLevel > 1) {
               G4cout << "Ksi Coefficient for Z= " << Z
                  << " E(keV)= " << lowEdgeEnergy/keV
                  << " Eav(keV)= " << eAverage/keV
                  << " coeff= " << coeff
                  << G4endl;
            }//if verboselevel

            energy->push_back(lowEdgeEnergy);
            ksi->push_back(coeff);
         }//for j<binForFluo
         interp = new G4LogLogInterpolation();
         G4VEMDataSet* set = new G4EMDataSet(Z,energy,ksi,interp,1.,1.);
         xsis->AddComponent(set);
      }//for iel<NumberOfElements
      if(verboseLevel) xsis->PrintData();
      shellVacancy->AddXsiTable(xsis);

      //EK, new: predefine all relevant ions (i.e., all single-charged positive ions of all elements in the material)
      if (isgas) {
         for (size_t iel=0; iel<NumberOfElements; iel++ ) {
            theelement = (*theElementVector)[iel];
            size_t niso = theelement->GetNumberOfIsotopes();
            if (niso==0) {
               G4int fAtomicNumber = (int)(theelement->GetZ()+0.5);// EK, rounding to the nearest integer.
               //GetZ should return (nearly) integer value anyway
               G4int fAtomicMass = (int)(theelement->GetN()+0.5);
               ion =  ionTable->FindIon( fAtomicNumber, fAtomicMass, 0);// EK, this should suffice to create the ion...
               if (ion) {// EK, now we add Transportation and a kill process
                  G4ProcessManager* pmanager = ion->GetProcessManager();
                  if (pmanager && pmanager->GetProcessListLength()==0) {// EK, Only if no processes have been added already...
                     G4Transportation* theTransportationProcess = new G4Transportation();
                     G4VProcess* ionkill                        = new CADPhysicsIonKill();
                     // EK, add transportation with ordering set to first
                     pmanager ->AddProcess(theTransportationProcess);
                     pmanager ->AddDiscreteProcess(ionkill);
                     pmanager ->SetProcessOrdering(theTransportationProcess, idxAlongStep,0);
                     pmanager ->SetProcessOrdering(theTransportationProcess, idxPostStep,0);
                     G4cerr << "CADPhysicsLowEnergyIonisation: Calling SetProcessOrdering" << G4endl;
                     pmanager ->SetProcessOrdering(ionkill,                  idxPostStep,1);
                     if (verboseLevel>1) G4cout << "Succesfully created ion with Z=" << fAtomicNumber << " and A=" << fAtomicMass << G4endl;
                  } else {
                     if (verboseLevel>1) G4cout << "Ion with Z=" << fAtomicNumber << " and A=" << fAtomicMass << " already existed " << G4endl;
                  }
                  G4cerr << "CADPhysicsLowEnergyIonisation: After call to GetProcessmanager" << G4endl;
               } else {
                  G4cerr << "Failed to create ion with Z=" << fAtomicNumber << " and A=" << fAtomicMass << G4endl;
               }
            } else {// EK, niso>0, the ion has multiple isotopes
               G4IsotopeVector* isovec = theelement->GetIsotopeVector();
               for (size_t iiso=0; iiso<niso; iiso++) {
                  isotope = (*isovec)[iiso];
                  G4int fAtomicNumber = isotope->GetZ();
                  G4int fAtomicMass = isotope->GetN();
                  ion =  ionTable->FindIon( fAtomicNumber, fAtomicMass, 0);// EK, this should suffice to create the ion...
                  if (ion) {// EK, now we add Transportation and a kill process
                     G4ProcessManager* pmanager = ion->GetProcessManager();
                     if (pmanager && pmanager->GetProcessListLength()==0) {// EK, Only if no processes have been added already...
                        G4Transportation*   theTransportationProcess   = new G4Transportation();
                        G4VProcess*         ionkill                  = new CADPhysicsIonKill();
                        // EK, add transportation with ordering = ( -1, "first", "first" )
                        pmanager ->AddProcess(theTransportationProcess);
                        pmanager ->AddDiscreteProcess(ionkill);
                        pmanager ->SetProcessOrderingToFirst(theTransportationProcess, idxAlongStep);
                        pmanager ->SetProcessOrderingToFirst(theTransportationProcess, idxPostStep);
                        pmanager ->SetProcessOrdering(ionkill,idxPostStep,2);
                        if (verboseLevel>1) G4cout << "Succesfully created isotope ion with Z=" << fAtomicNumber << " and A=" << fAtomicMass << G4endl;
                     } else {
                        if (verboseLevel>1) G4cout << "Isotope ion with Z=" << fAtomicNumber << " and A=" << fAtomicMass << " already existed " << G4endl;
                     }
                  } else {
                     G4cerr << "Failed to create isotope ion with Z=" << fAtomicNumber << " and A=" << fAtomicMass << G4endl;
                  }
               }
            }
         }//for iel<NumberOfElements
      }//if (isgas)
      // EK, end of creation of all ions

   }//for m<numCouples
   delete bVector;
}

G4VParticleChange* CADPhysicsLowEnergyIonisation::AlongStepDoIt(
   const G4Track& track,const G4Step&)//track,step
{
   // Not implemented for this process
   aParticleChange.Initialize(track);
   return &aParticleChange;
}

G4VParticleChange* CADPhysicsLowEnergyIonisation::PostStepDoIt(const G4Track& track,
                                                const G4Step&  step)
{
   // Delta electron production mechanism on base of the model
   // J. Stepanek " A program to determine the radiation spectra due
   // to a single atomic subshell ionisation by a particle or due to
   // deexcitation or decay of radionuclides",
   // Comp. Phys. Comm. 1206 pp 1-19 (1997)

   aParticleChange.Initialize(track);

   const G4MaterialCutsCouple* couple = track.GetMaterialCutsCouple();
   G4double kineticEnergy = track.GetKineticEnergy();

   // Select atom and shell

   G4int Z = crossSectionHandler->SelectRandomAtom(couple, kineticEnergy);
   G4int shell = crossSectionHandler->SelectRandomShell(Z, kineticEnergy);
   const G4AtomicShell* atomicShell =
      (G4AtomicTransitionManager::Instance())->Shell(Z, shell);
   G4double bindingEnergy = atomicShell->BindingEnergy();
   G4int shellId = atomicShell->ShellId();

   // Sample delta energy

   G4int    index  = couple->GetIndex();
   G4double tCut   = cutForDelta[index];
   G4double tmax   = energySpectrum->MaxEnergyOfSecondaries(kineticEnergy);
   G4double tDelta = energySpectrum->SampleEnergy(Z, tCut, tmax,
      kineticEnergy, shell);

   if(tDelta == 0.0)
      return G4VContinuousDiscreteProcess::PostStepDoIt(track, step);

   // Transform to shell potential
   G4double deltaKinE = tDelta + 2.0*bindingEnergy;
   G4double primaryKinE = kineticEnergy + 2.0*bindingEnergy;

   // sampling of scattering angle neglecting atomic motion
   G4double deltaMom = std::sqrt(deltaKinE*(deltaKinE + 2.0*electron_mass_c2));
   G4double primaryMom = std::sqrt(primaryKinE*(primaryKinE + 2.0*electron_mass_c2));

   G4double cost = deltaKinE * (primaryKinE + 2.0*electron_mass_c2)
      / (deltaMom * primaryMom);

   if (cost > 1.) cost = 1.;
   G4double sint = std::sqrt(1. - cost*cost);
   G4double phi  = twopi * G4UniformRand();
   G4double dirx = sint * std::cos(phi);
   G4double diry = sint * std::sin(phi);
   G4double dirz = cost;

   // Rotate to incident electron direction
   G4ThreeVector primaryDirection = track.GetMomentumDirection();
   G4ThreeVector deltaDir(dirx,diry,dirz);
   deltaDir.rotateUz(primaryDirection);
   dirx = deltaDir.x();
   diry = deltaDir.y();
   dirz = deltaDir.z();

   // Take into account atomic motion del is relative momentum of the motion
   // kinetic energy of the motion == bindingEnergy in V.Ivanchenko model
   cost = 2.0*G4UniformRand() - 1.0;
   sint = std::sqrt(1. - cost*cost);
   phi  = twopi * G4UniformRand();
   G4double del = std::sqrt(bindingEnergy *(bindingEnergy + 2.0*electron_mass_c2))/deltaMom;
   dirx += del * sint * std::cos(phi);
   diry += del * sint * std::sin(phi);
   dirz += del * cost;

   // Construct a random ion direction
   cost = 2.0*G4UniformRand() - 1.0;
   sint = std::sqrt(1. - cost*cost);
   phi  = twopi * G4UniformRand();
   G4double iondirx = sint * std::cos(phi);
   G4double iondiry = sint * std::sin(phi);
   G4double iondirz = cost;

   // Find out new primary electron direction
   G4double finalPx = primaryMom*primaryDirection.x() - deltaMom*dirx;
   G4double finalPy = primaryMom*primaryDirection.y() - deltaMom*diry;
   G4double finalPz = primaryMom*primaryDirection.z() - deltaMom*dirz;

   // create G4DynamicParticle object for delta ray
   G4DynamicParticle* theDeltaRay = new G4DynamicParticle();
   theDeltaRay->SetKineticEnergy(tDelta);
   G4double norm = 1.0/std::sqrt(dirx*dirx + diry*diry + dirz*dirz);
   dirx *= norm;
   diry *= norm;
   dirz *= norm;
   theDeltaRay->SetMomentumDirection(dirx, diry, dirz);
   theDeltaRay->SetDefinition(G4Electron::Electron());

   G4double theEnergyDeposit = bindingEnergy;

   // fill ParticleChange
   // changed energy and momentum of the actual particle

   G4double finalKinEnergy = kineticEnergy - tDelta - theEnergyDeposit;
   if(finalKinEnergy < 0.0) {
      theEnergyDeposit += finalKinEnergy;
      finalKinEnergy    = 0.0;
      aParticleChange.ProposeTrackStatus(fStopAndKill);

   } else {

      G4double norm = 1.0/std::sqrt(finalPx*finalPx+finalPy*finalPy+finalPz*finalPz);
      finalPx *= norm;
      finalPy *= norm;
      finalPz *= norm;
      aParticleChange.ProposeMomentumDirection(finalPx, finalPy, finalPz);
   }

   // EK, check on the final energy for the primary electron
   if (useCutForElectrons && finalKinEnergy<cutForElectrons) {
      aParticleChange.ProposeTrackStatus(fStopAndKill);
   } else {
      aParticleChange.ProposeEnergy(finalKinEnergy);
   }

   // Generation of Fluorescence and Auger
   size_t nSecondaries = 0;
   size_t totalNumber  = 1;

   std::vector<G4DynamicParticle*>* secondaryVector = 0;
   G4DynamicParticle* aSecondary = 0;
   G4ParticleDefinition* type = 0;

   // Fluorescence data start from element 6

   if (Fluorescence() && Z > 5 && (bindingEnergy >= cutForPhotons
      ||  bindingEnergy >= cutForElectrons)) {

         secondaryVector = deexcitationManager.GenerateParticles(Z, shellId);

         if (secondaryVector != 0) {

            nSecondaries = secondaryVector->size();
            for (size_t i = 0; i<nSecondaries; i++) {

               aSecondary = (*secondaryVector)[i];
               if (aSecondary) {

                  G4double e = aSecondary->GetKineticEnergy();
                  type = aSecondary->GetDefinition();
                  if (e < theEnergyDeposit &&
                     ((type == G4Gamma::Gamma() && e > cutForPhotons ) ||
                     (type == G4Electron::Electron() && e > cutForElectrons ))) {

                        theEnergyDeposit -= e;
                        totalNumber++;

                  } else {

                     delete aSecondary;
                     (*secondaryVector)[i] = 0;
                  }
               }
            }
         }
   }

   // Save delta-electrons
   // EK, check on the energy for the secondary electron
   G4bool rejectedtheDeltaRay=false;
   if (useCutForElectrons && tDelta<cutForElectrons) {
      //G4cout << "tDelta rejected" << G4endl;
      rejectedtheDeltaRay = true;
      totalNumber--;
      aParticleChange.SetNumberOfSecondaries(totalNumber);
      delete theDeltaRay;
   } else {
      aParticleChange.SetNumberOfSecondaries(totalNumber);
      //aParticleChange.AddSecondary(theDeltaRay);
      //Actually adding theDeltaRay postponed till later - first we need to find out if the totalNumber is affected by generating an ion.
   }

   // EK, generate the secondary ion
   if(generateIons) {
      // EK, first find the corresponding element
      const G4Material* material= couple->GetMaterial();
      size_t NumberOfElements = material->GetNumberOfElements();
      G4Element* theelement = 0;
      const G4ElementVector* theElementVector = material->GetElementVector();
      size_t iel;
      for (iel=0;iel<NumberOfElements;iel++) {
         theelement = (*theElementVector)[iel];
         G4int elZ = (G4int)(theelement->GetZ());
         if (elZ==Z) break;
      }
      G4DynamicParticle* theIon = new G4DynamicParticle();
      if (iel==NumberOfElements) {
         delete theIon;
      } else {
         // EK, select the isotope (if there's more than one)
         G4IsotopeVector* isv = theelement->GetIsotopeVector();
         G4int ni = 0;
         G4int A = 0;//The atomic mass in amu
         if(!isv) {
            A = (int)(theelement->GetN()+0.5);
         } else {
            ni = isv->size();
            if(ni == 1) {
               A = theelement->GetIsotope(0)->GetN();
            } else if(ni > 1) {
               G4double* ab = theelement->GetRelativeAbundanceVector();
               G4double y = G4UniformRand();
               G4int j = -1;
               ni--;
               do {
                  j++;
                  y -= ab[j];
               } while (y > 0.0 && j < ni);
               A = theelement->GetIsotope(j)->GetN();
            }
         }
         // EK, create new dynamic particle
         // EK, simply assume here that momentum transfer to the ion is negligible compared to thermal motion
         G4double energy = -std::log(G4UniformRand())* 1.5 * 0.025*eV;// EK, roughly room temperature
         if (verboseLevel>2) G4cout << "Generating ion with Z=" << Z << ", A=" << A << ", energy=" << energy/eV << " eV." << G4endl;
         theIon->SetKineticEnergy(energy);
         theIon->SetMomentumDirection(iondirx,iondiry,iondirz);
         theIon->SetDefinition(ionTable->FindIon( Z, A, 0.));
         theIon->SetCharge(1);
         totalNumber++;
         aParticleChange.SetNumberOfSecondaries(totalNumber);
         aParticleChange.AddSecondary(theIon);
      }
   }
   // EK, end of the ion generation part

   if(!rejectedtheDeltaRay) {
      aParticleChange.AddSecondary(theDeltaRay);
   }

   // Save Fluorescence and Auger
   if (secondaryVector) {
      for (size_t l = 0; l < nSecondaries; l++) {
         aSecondary = (*secondaryVector)[l];
         if(aSecondary) {
            aParticleChange.AddSecondary(aSecondary);
         }
      }
      delete secondaryVector;
   }

   if(theEnergyDeposit < 0.) {
      G4cout << "CADPhysicsLowEnergyIonisation: Negative energy deposit: "
         << theEnergyDeposit/eV << " eV" << G4endl;
      theEnergyDeposit = 0.0;
   }
   aParticleChange.ProposeLocalEnergyDeposit(theEnergyDeposit);

   return G4VContinuousDiscreteProcess::PostStepDoIt(track, step);
}


void CADPhysicsLowEnergyIonisation::PrintInfoDefinition()
{
   G4String comments = "Total cross sections from EEDL database.";
   comments += "\n      Gamma energy sampled from a parametrised formula.";
   comments += "\n      Implementation of the continuous dE/dx part.";
   comments += "\n      At present it can be used for electrons ";
   comments += "in the energy range [250eV,100GeV].";
   comments += "\n      The process must work with G4LowEnergyBremsstrahlung.";

   G4cout << G4endl << GetProcessName() << ":  " << comments << G4endl;
}

G4bool CADPhysicsLowEnergyIonisation::IsApplicable(const G4ParticleDefinition& particle)
{
   return ( (&particle == G4Electron::Electron()) );
}

std::vector<G4DynamicParticle*>*
CADPhysicsLowEnergyIonisation::DeexciteAtom(const G4MaterialCutsCouple* couple,
                                 G4double incidentEnergy,
                                 G4double eLoss)
{
   // create vector of secondary particles
   const G4Material* material = couple->GetMaterial();

   std::vector<G4DynamicParticle*>* partVector =
      new std::vector<G4DynamicParticle*>;

   if(eLoss > cutForPhotons && eLoss > cutForElectrons) {

      const G4AtomicTransitionManager* transitionManager =
         G4AtomicTransitionManager::Instance();

      size_t nElements = material->GetNumberOfElements();
      const G4ElementVector* theElementVector = material->GetElementVector();

      std::vector<G4DynamicParticle*>* secVector = 0;
      G4DynamicParticle* aSecondary = 0;
      G4ParticleDefinition* type = 0;
      G4double e;
      G4ThreeVector position;
      G4int shell, shellId;

      // sample secondaries

      G4double eTot = 0.0;
      std::vector<G4int> n =
         shellVacancy->GenerateNumberOfIonisations(couple,
         incidentEnergy,eLoss);
      for (size_t i=0; i<nElements; i++) {

         G4int Z = (G4int)((*theElementVector)[i]->GetZ());
         size_t nVacancies = n[i];

         G4double maxE = transitionManager->Shell(Z, 0)->BindingEnergy();

         if (nVacancies && Z > 5 && (maxE>cutForPhotons || maxE>cutForElectrons)) {

            for (size_t j=0; j<nVacancies; j++) {

               shell = crossSectionHandler->SelectRandomShell(Z, incidentEnergy);
               shellId = transitionManager->Shell(Z, shell)->ShellId();
               G4double maxEShell =
                  transitionManager->Shell(Z, shell)->BindingEnergy();

               if (maxEShell>cutForPhotons || maxEShell>cutForElectrons ) {

                  secVector = deexcitationManager.GenerateParticles(Z, shellId);

                  if (secVector != 0) {

                     for (size_t l = 0; l<secVector->size(); l++) {

                        aSecondary = (*secVector)[l];
                        if (aSecondary != 0) {

                           e = aSecondary->GetKineticEnergy();
                           type = aSecondary->GetDefinition();
                           if ( eTot + e <= eLoss &&
                              ((type == G4Gamma::Gamma() && e>cutForPhotons ) ||
                              (type == G4Electron::Electron() && e>cutForElectrons))) {

                                 eTot += e;
                                 partVector->push_back(aSecondary);

                           } else {

                              delete aSecondary;

                           }
                        }
                     }
                     delete secVector;
                  }
               }
            }
         }
      }
   }
   return partVector;
}

G4double CADPhysicsLowEnergyIonisation::GetMeanFreePath(const G4Track& track,
                                          G4double , // previousStepSize
                                          G4ForceCondition* cond)
{
   *cond = NotForced;

   // EK, act only on material in the gas state
   const G4Material* material = (track.GetMaterialCutsCouple())->GetMaterial();
   if (material->GetState()!=kStateGas) return DBL_MAX;

   G4int index = (track.GetMaterialCutsCouple())->GetIndex();
   const G4VEMDataSet* data = theMeanFreePath->GetComponent(index);
   G4double meanFreePath = data->FindValue(track.GetKineticEnergy());
   return meanFreePath;
}

void CADPhysicsLowEnergyIonisation::SetCutForLowEnSecPhotons(G4double cut)
{
   cutForPhotons = cut;
   deexcitationManager.SetCutForSecondaryPhotons(cut);
}

void CADPhysicsLowEnergyIonisation::SetCutForLowEnSecElectrons(G4double cut)
{
   cutForElectrons = cut;
   deexcitationManager.SetCutForAugerElectrons(cut);
}

G4double CADPhysicsLowEnergyIonisation::GetCutForLowEnSecElectrons()
{
   return cutForElectrons;
}

void CADPhysicsLowEnergyIonisation::SetUseCut(G4bool usecut)
{
   useCutForElectrons = usecut;
}

G4bool CADPhysicsLowEnergyIonisation::GetUseCut()
{
   return useCutForElectrons;
}

void CADPhysicsLowEnergyIonisation::SetGenerateIons(G4bool gi)
{
   generateIons = gi;
}

G4bool CADPhysicsLowEnergyIonisation::GetGenerateIons()
{
   return generateIons;
}

void CADPhysicsLowEnergyIonisation::ActivateAuger(G4bool val)
{
   deexcitationManager.ActivateAugerElectronProduction(val);
}

