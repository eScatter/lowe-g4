// CADPhysicsLowEnergyIonisation.hh
//
// Process based on G4LowEnergyIonisation, Geant4 v. 8.0p1.
// This class inherits from G4eLowEnergyLoss. CADPhysicsLowEnergyIonisation is used 
// as an alternative to the DI process in the case of gas state materials.
// Accurate optical data for gases seem a bit hard to come by, and since there are no solid state effects to worry about anyway,
// this 'standard' Geant4 process seems to be a good alternative.
//
// The main modification from the standard code is that the process checks the physical state that the material is in,
// (both in BuildLossTable - for the continuous part - and in GetMeanFreePath - for the discrete part) and acts only if it is a gas.
// 
// Furthermore, some code has been added to generate ions next to electrons as secondary particles.
// A switch has been introduced so that, via a messenger, the user can turn generation of ions on and off.
// Another new option is to set an energy cut for primary/secondary electrons.
// A new messenger (CADPhysicsLowEnIoniMessenger) has been created to give the user control over these variables.
//
// Some final remarks on how to use this process: 
// - The original version was written such that it has to coexist with G4LowEnergyBremsstrahlung. Hence if the current process is
//   defined, then G4LowEnergyBremsstrahlung needs to be defined as well.
// - However, only ONE of the two should be added as a CONTINUOUS process using a call like the following:
//     pManager->SetProcessOrdering(theeminusLoweIonisation,	 idxAlongStep,2);
//   It is advised to use the current process rather than G4LowEnergyBremsstrahlung, since the latter will kill all electrons 
//   with energy below 1 eV, and this is not always desired.
// 
// Comments that mark modifications from the original Geant4 code are preceded by 'EK'.
// 
// 2007-11-28: Replaced use of G4AtomicDeexcitation by CADPhysicsAtomicDeexcitation
//             for consistency with CADPhysicsDI.
// 2007-11-30: Inheriting from G4eLowEnergyLoss instead of special CADPhysicseLowEnergyLoss.
//             Accordingly, setting MinKineticEnergy = DBL_MIN in the constructor.

#ifndef CADPhysicsLowENERGYIONISATION_HH
#define CADPhysicsLowENERGYIONISATION_HH 1

#include "G4eLowEnergyLoss.hh"
#include "CADPhysicsAtomicDeexcitation.hh"
#include "G4IonTable.hh"

class G4Track;
class G4Step;
class G4ParticleDefinition;
class G4VParticleChange;
class G4VDataSetAlgorithm;
class G4ParticleChange;
class G4VEnergySpectrum;
class G4VCrossSectionHandler;
class G4ShellVacancy;
class G4VEMDataSet;
class G4ParticleTable;// EK, for creation of ions

class CADPhysicsLowEnIoniMessenger;

class CADPhysicsLowEnergyIonisation : public G4eLowEnergyLoss
{ 
public:
 
  CADPhysicsLowEnergyIonisation(const G4String& processName = "LowEnergyIoni");
  
  ~CADPhysicsLowEnergyIonisation();
  
  G4bool IsApplicable(const G4ParticleDefinition&);
  
  void PrintInfoDefinition();
  
  void BuildPhysicsTable(const G4ParticleDefinition& ParticleType);

  G4VParticleChange* AlongStepDoIt(const G4Track& aTrack,const G4Step& aStep);
  // Perform the continuous part of the process
  
  G4VParticleChange* PostStepDoIt(const G4Track& track,         
				  const G4Step& step);                 
 
  void SetCutForLowEnSecPhotons(G4double cut);

  void SetCutForLowEnSecElectrons(G4double cut);
  G4double GetCutForLowEnSecElectrons();// EK, 'get' method for the corresponding 'set' method

  // EK, new methods for communication with the messenger
  void SetUseCut(G4bool usecut);
  G4bool GetUseCut();
  void SetGenerateIons(G4bool generateions);
  G4bool GetGenerateIons();
  
  void ActivateAuger(G4bool val);
    
protected:
 
  G4double GetMeanFreePath(const G4Track& track,
			   G4double previousStepSize,
			   G4ForceCondition* condition );

protected:

  virtual std::vector<G4DynamicParticle*>* DeexciteAtom(const G4MaterialCutsCouple* couple,
							  G4double incidentEnergy,
							  G4double eLoss);

private:

  // Hide copy constructor and assignment operator as private 
  CADPhysicsLowEnergyIonisation(const CADPhysicsLowEnergyIonisation& );
  CADPhysicsLowEnergyIonisation& operator = (const CADPhysicsLowEnergyIonisation& right);
  
  CADPhysicsLowEnIoniMessenger* messenger;

  void BuildLossTable(const G4ParticleDefinition& ParticleType);
  
  G4VCrossSectionHandler* crossSectionHandler;
  G4VEMDataSet* theMeanFreePath;
  G4VEnergySpectrum* energySpectrum;

  // Lower limit for generation of gamma in this model
  G4DataVector cutForDelta;
  G4double cutForPhotons;
  G4double cutForElectrons;
  G4bool useCutForElectrons; //EK, boolean to switch off/on the application of the energy cut for electrons
  G4bool generateIons; //EK, boolean to switch off/on the generation of an ion whenever a secondary electron is generated
  CADPhysicsAtomicDeexcitation deexcitationManager;
  G4ShellVacancy* shellVacancy;

  // EK, for creation of ions
  //G4ParticleTable* particleTable;
  G4IonTable* ionTable;
  G4double electron_mass_over_amu;
  
};

#endif
 










