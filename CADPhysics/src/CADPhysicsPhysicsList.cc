// CADPhysicsPhysicsList.cc
//
// See the header file for a general description 
//
// History:
// 2010-09-08: Added Rutherford elastic scattering model for energies above 30keV
// 2007-08-16: Added Penelope 'low-energy' models for gammas
//

#include "CADPhysicsPhysicsList.hh"
#include "globals.hh"
#include "G4ProcessManager.hh"
#include <iomanip>

CADPhysicsPhysicsList::CADPhysicsPhysicsList(): G4VUserPhysicsList()
{
	// Some initialization
	defaultCutValue = 0.;// Modifying the default G4VUserPhysicsList value of 1 mm. Deep inside Geant4, this 'minimum range' is converted
	// to a minimum energy for e+, e-, and gamma based some continuous energy loss model.
	// We don't want to be affected by this behavior and hence simply set this minimum range value to zero.
	fCutsTable->SetEnergyRange(0.10001*eV, 1.*MeV);// Modifying the default G4VUserPhysicsList values of 0.99*keV and 100*TeV. The higher limit
	// is not used anywhere in Geant4, but the lower limit is a lower bound to the 'minimum energy' as derived from the 'minimum range'.
	// The resulting number is used in certain processes - including CADPhysicsLowEnergyIonisation, G4Bremsstrahlung, and the G4Penelope processes - 
	// as the lowest kinetic energy for secondary particles.
	// Hence the abovementioned processes will not generate particles with energy lower than 0.10001 eV.
	SetVerboseLevel(1);// Setting the verbosity of this class and fCutsTable to 1.
        G4cout << "Calling CADPhysicsPhysicsList. Version December 2007." << G4endl;
        G4cout << "- Includes Rutherford scattering model above 30keV. September 2010." << G4endl;
}

CADPhysicsPhysicsList::~CADPhysicsPhysicsList() {}

#include "G4Electron.hh"
#include "G4Geantino.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"

void CADPhysicsPhysicsList::ConstructParticle()
{
	// Define all particles
	G4Electron::ElectronDefinition();
	G4Geantino::GeantinoDefinition();
	G4Positron::PositronDefinition();
	G4Gamma::GammaDefinition();
}

#include "G4Transportation.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4PenelopePhotoElectricModel.hh"
//#include "G4LivermorePhotoElectricModel.hh"
#include "G4ComptonScattering.hh"
#include "G4PenelopeComptonModel.hh"
//#include "G4LivermoreComptonModel.hh"
#include "G4GammaConversion.hh"
#include "G4PenelopeGammaConversionModel.hh"
//#include "G4LivermoreGammaConversionModel.hh"
#include "G4RayleighScattering.hh"
#include "G4PenelopeRayleighModel.hh"
//#include "G4LivermoreRayleighModel.hh"
#include "CADPhysicsTransportation.hh"
#include "CADPhysicsSingleScattering.hh"
#include "CADPhysicsDI.hh"
#include "CADPhysicseBoundary.hh"
#include "CADPhysicsLowEnergyIonisation.hh"
#include "G4LowEnergyBremsstrahlung.hh"
#include "CADPhysicsKill.hh"
#include "G4EmProcessOptions.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"
#include "G4PhysicsListHelper.hh"
#include "CADPhysicsAtomicDeexcitation.hh"

void CADPhysicsPhysicsList::ConstructProcess()
{
	// Initialization
	G4ProcessManager * pManager = 0;
	G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

	// Define the 'standard' transportation process, used for all particles except electrons
	G4Transportation* thegeneralTransportation = new G4Transportation();

	//////////////
	// Geantino //
	//////////////
	pManager = G4Geantino::Geantino()->GetProcessManager();
	pManager ->AddProcess(thegeneralTransportation);
	pManager ->SetProcessOrderingToFirst(thegeneralTransportation, idxAlongStep);
	pManager ->SetProcessOrderingToFirst(thegeneralTransportation, idxPostStep);

	///////////
	// Gamma //
	///////////
	pManager = G4Gamma::Gamma()->GetProcessManager();
	G4double PenelopeHighEnergyLimit = 1.0*GeV;

	// Define processes
	G4PhotoElectricEffect* gammaPhotoElectric = new G4PhotoElectricEffect();
	G4PenelopePhotoElectricModel* thePEPenelopeModel = new G4PenelopePhotoElectricModel();   
	thePEPenelopeModel->SetHighEnergyLimit(PenelopeHighEnergyLimit);
	gammaPhotoElectric->SetEmModel(thePEPenelopeModel, 1);
//	G4LivermorePhotoElectricModel* thePELivermoreModel = new G4LivermorePhotoElectricModel();
//	thePELivermoreModel->SetHighEnergyLimit(PenelopeHighEnergyLimit);
//	gammaPhotoElectric->SetEmModel(thePELivermoreModel,1);
	
	G4ComptonScattering* gammaCompton       = new G4ComptonScattering();
	G4PenelopeComptonModel* theComptonPenelopeModel = new G4PenelopeComptonModel();
	theComptonPenelopeModel->SetHighEnergyLimit(PenelopeHighEnergyLimit);
	gammaCompton->SetEmModel(theComptonPenelopeModel, 1);
//	G4LivermoreComptonModel* theComptonLivermoreModel = new G4LivermoreComptonModel();
//	theComptonLivermoreModel->SetHighEnergyLimit(PenelopeHighEnergyLimit);
//	gammaCompton->SetEmModel(theComptonLivermoreModel,1);

	G4GammaConversion* gammaConversion    = new G4GammaConversion();
	G4PenelopeGammaConversionModel* theGCPenelopeModel = new G4PenelopeGammaConversionModel();
	gammaConversion->SetEmModel(theGCPenelopeModel,1);
//	G4LivermoreGammaConversionModel* theGCLivermoreModel = new G4LivermoreGammaConversionModel();
//	gammaConversion->SetEmModel(theGCLivermoreModel,1);
	
	G4RayleighScattering* gammaRayleigh      = new G4RayleighScattering();
	G4PenelopeRayleighModel* theRayleighPenelopeModel = new G4PenelopeRayleighModel();
	//theRayleighPenelopeModel->SetHighEnergyLimit(PenelopeHighEnergyLimit);
	gammaRayleigh->SetEmModel(theRayleighPenelopeModel, 1);
//	G4LivermoreRayleighModel* theRayleighLivermoreModel = new G4LivermoreRayleighModel();
//	gammaRayleigh->SetEmModel(theRayleighLivermoreModel,1);

	// Subscribe processes to the process manager
//	pManager->AddProcess(thegeneralTransportation);
//	pManager->AddDiscreteProcess(gammaPhotoElectric);
//	pManager->AddDiscreteProcess(gammaCompton);
//	pManager->AddDiscreteProcess(gammaConversion);
//	pManager->AddDiscreteProcess(gammaRayleigh);

	ph->RegisterProcess(thegeneralTransportation,G4Gamma::Gamma());
	ph->RegisterProcess(gammaPhotoElectric,G4Gamma::Gamma());
	ph->RegisterProcess(gammaCompton,G4Gamma::Gamma());
	ph->RegisterProcess(gammaConversion,G4Gamma::Gamma());
	ph->RegisterProcess(gammaRayleigh,G4Gamma::Gamma());

	// Set ordering for AlongStepDoIt (only one process)
//	pManager->SetProcessOrdering(thegeneralTransportation, idxAlongStep,0);

	// Set ordering for PostStepDoIt
//	pManager->SetProcessOrdering(thegeneralTransportation, idxPostStep,0);
//	pManager->SetProcessOrdering(gammaPhotoElectric,       idxPostStep,1);
//	pManager->SetProcessOrdering(gammaCompton,             idxPostStep,2);
//	pManager->SetProcessOrdering(gammaConversion,          idxPostStep,3);
//	pManager->SetProcessOrdering(gammaRayleigh,            idxPostStep,4);


	//////////////
	// Positron //
	//////////////
	pManager = G4Positron::Positron()->GetProcessManager();
	pManager ->AddProcess(thegeneralTransportation);
	pManager ->SetProcessOrderingToFirst(thegeneralTransportation, idxAlongStep);
	pManager ->SetProcessOrderingToFirst(thegeneralTransportation, idxPostStep);


	//////////////
	// Electron //
	//////////////
	// The G4LEDATA environment variable is required to be set BEFORE the constructor of CADPhysicsDI is called. Hence we set it here.
	char* dummy = getenv("G4LEDATA");
	if (!dummy) {
		G4cout << "Environment variable G4LEDATA undefined. Setting it to default, ../G4EMLOW" << G4endl;
                sprintf(dummy,"G4LEDATA=../G4EMLOW");
		putenv(dummy);// G4LEDATA is set only if not defined previously
	}

	pManager = G4Electron::Electron()->GetProcessManager();

	// Define processes
	G4VProcess* theeminusTransportation     = CADPhysicsTransportation::GetInstance();// This process has a 'singleton' structure
	G4VProcess* theeminusElasticScattering  = CADPhysicsSingleScattering::GetInstance();// This process has a 'singleton' structure
	G4VProcess* theeminusIonisation         = new CADPhysicsDI();
	G4VProcess* theeminusBoundary           = new CADPhysicseBoundary();
	G4VProcess* theeminusLoweIonisation     = new CADPhysicsLowEnergyIonisation();
	G4VProcess* theeminusLoweBremsstrahlung = new G4LowEnergyBremsstrahlung();
	G4VProcess* theeminusKill               = new CADPhysicsKill();

	// Subscribe processes to the process manager
	pManager->AddProcess(theeminusTransportation);
	pManager->AddDiscreteProcess(theeminusElasticScattering);
	pManager->AddDiscreteProcess(theeminusIonisation);
	pManager->AddDiscreteProcess(theeminusBoundary);
	pManager->AddProcess(theeminusLoweIonisation);
	pManager->AddProcess(theeminusLoweBremsstrahlung);
	pManager->AddProcess(theeminusKill);

	// Set ordering for AlongStepDoIt
	pManager->SetProcessOrdering(theeminusTransportation,     idxAlongStep,0);
	pManager->SetProcessOrdering(theeminusElasticScattering,  idxAlongStep,1);
	pManager->SetProcessOrdering(theeminusLoweIonisation,     idxAlongStep,2);

	// Set ordering for PostStepDoIt
	pManager->SetProcessOrdering(theeminusTransportation,     idxPostStep,0);
	pManager->SetProcessOrdering(theeminusElasticScattering,  idxPostStep,1);
	pManager->SetProcessOrdering(theeminusIonisation,         idxPostStep,2);
	pManager->SetProcessOrdering(theeminusBoundary,           idxPostStep,3);
	pManager->SetProcessOrdering(theeminusLoweIonisation,     idxPostStep,4);
	pManager->SetProcessOrdering(theeminusLoweBremsstrahlung, idxPostStep,5);
	pManager->SetProcessOrdering(theeminusKill,               idxPostStep,6);

	// Set individual process parameter(s)
	G4EmProcessOptions emOptions;
	emOptions.SetMinEnergy(100*eV);
	emOptions.SetFluo(true);
	emOptions.SetAuger(true);

	// Deexcitation
	G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
	G4LossTableManager::Instance()->SetAtomDeexcitation(de);
	de->SetDeexcitationActiveRegion("World",true,true,true);
	de->SetFluo(true);
	de->SetAuger(true);
	de->InitialiseForNewRun();
}

void CADPhysicsPhysicsList::SetCuts()
{
	SetCutsWithDefault();// Calling the G4VUserPhysicsList member method. Assigns the defaultCutValue to e-, e+ and gamma.
	G4double lowLimit = 250. * eV;
	G4double highLimit = 100. * GeV;
	G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowLimit, highLimit);
}
