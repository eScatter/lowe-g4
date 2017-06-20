// CADPhysicsDI.hh
//
// Process for inelastic scattering in (non-gaseous) materials. "DI" stands for the Dielectric
// Ionisation formalism that is applied.
// Several authors have used versions of this formalism. The present implementation of this class
// is most closely based on J.C. Ashley's description in
// J. Electr. Spectrosc. Rel. Phenom. 46: 199-214 (1988).
// This method first calculates total cross sections and differential cross sections for omega
// prime. In PostStepDoIt, omega is probed according to a function F(E,omega prime, omega)
// (eq. 7 in Ashley). See the detailed comments for more information.
// Note: this process is only active on non-gaseous materials.
//

#ifndef CADPhysicsDI_h
#define CADPhysicsDI_h 1

#include "G4VContinuousDiscreteProcess.hh"
#include "G4GPILSelection.hh"
#include "G4Electron.hh"
#include "G4DataVector.hh"
#include "CADPhysicsAtomicDeexcitation.hh"

// for read in hdf5
#include <csread/material.h>
#include <csread/units/unit_system.h>

typedef std::vector<G4DataVector*> CADPhysicsDataTable;

class CADPhysicsDIMessenger;
class G4CrossSectionHandler;

class CADPhysicsDI : public G4VContinuousDiscreteProcess

{
public:

   CADPhysicsDI(const G4String& processName="DielIoni");

   ~CADPhysicsDI();

   G4bool IsApplicable ( const G4ParticleDefinition& );
      // Return whether the current process is relevant to the given particle

   void load_material(std::string const & filename, double high_energy, size_t N_K, size_t N_P);
   void load_vacuum();

   void BuildPhysicsTable(const G4ParticleDefinition& aParticleType);
      // This function overloads the corresponding virtual function
      // of the base class G4VContinuousDiscreteProcess.
      // It is an initialization routine invoked from G4VUserPhysicsList::BuildPhysicsTable()
      // at the start of the first run (and after geometry modification)

   void PrintInfoDefinition();
      // Print information about the process, invoked by BuildPhysicsTable()

   G4double PostStepGetPhysicalInteractionLength( const G4Track& track,
      G4double   previousStepSize,
      G4ForceCondition* condition
      );
      // Return the actual step size to be taken before the next discrete event of this process,
      // and set the 'G4ForceCondition', i.e. determine whether a call to the PostStepDoIt function
      // should be forced for this step. Overloading the base class with the sole purpose of
      // adjusting the verbosity levels and corresponding output. On the meaning of various values
      // of verboseLevel:
      // 0  - Silent except for essential warning messages
      // 1  - Basic information during initialization
      // 2  - More extensive information during initialization
      // 3  - For output of e.g. cross section tables to screen during initialization
      // >3 - For (debugging) output during simulations, e.g. in
      //      PostStepGetPhysicalInteractionLength, or perhaps in PostStepDoIt

   G4double AlongStepGetPhysicalInteractionLength(const G4Track&,
      G4double  previousStepSize,
      G4double  currentMinimumStep,
      G4double& currentSafety,
      G4GPILSelection* selection);
      // Limit the step size for the continuous part of the process.
      // The function overloads the corresponding function of the base
      // class. It invokes the method GetContinuousStepLimit at every step. The overloaded version
      // follows the default behavior.

   G4double GetContinuousStepLimit(const G4Track& aTrack,
      G4double previousStepSize,
      G4double currentMinimumStep,
      G4double& currentSafety);
      // Returns the step size limit for the continuous part of the process.
      // Invoked by the AlongStepGetPhysicalInteractionLength method.

   G4double GetMeanFreePath(const G4Track& aTrack,
      G4double previousStepSize,
      G4ForceCondition* condition);
      // This function overloads a virtual function of the base class, and is invoked by
      // PostStepGetPhysicalInteractionLength()
      // Returns the mean free path of the process for this step.

   G4VParticleChange* AlongStepDoIt(const G4Track& aTrack,const G4Step& aStep);
      // Perform the continuous part of the process

   G4VParticleChange* PostStepDoIt(const G4Track& aTrack,const G4Step& aStep);
      // Perform the discrete part of the process: determine the final state of the particle
      // after the step, and generate secondary particles (if any)

   // Public inline methods
   inline void SetGenerateSecondaries(G4bool gensec) {
      // Method to determine whether secondary particles should be generated.
      generateSecondaries = gensec;
   }
   inline G4bool GetGenerateSecondaries() {
      // Corresponding 'get' method (invoked by the Messenger)
      return generateSecondaries;
   }


   inline void SetGenerateXrays(G4bool genxray) {
      // Method to determine whether X ray photons should be generated.
      generateXrays = genxray;
   }

   inline G4bool GetGenerateXrays() {
      // Corresponding 'get' method (invoked by the Messenger)
      return generateXrays;
   }

   inline void ResetCounter() {
      // Method called from the Messenger for counting generated secondary electrons
      // (for debugging purposes)
      G4cout << "DI generated total of " << pairsgenerated << " pairs." << G4endl;
      pairsgenerated = 0;
   }

   inline void SetOutput(G4bool kl) {DI_output=kl;}// Method (called by the messenger) to set DI_output

private:
   // Hiding the assignment operator and the copy constructor as private,
   // since it doesn't make sense to have multiple instances of this class
   CADPhysicsDI & operator = (const CADPhysicsDI &right);
   CADPhysicsDI ( const CADPhysicsDI &);

   // Private methods:
   void BuildEnrange();
      // Fill the enrange, Psecenergies and Psecvalues vectors. Invoked at initialization
      // (by the BuildPhysicsTable method).

   G4VParticleChange* Phononloss(const G4Track& aTrack,
      const G4Step& aStep,
      G4double kinenergy,
      G4double omegaprime);
      // Process the event in case of sub-bandgap energy loss in semiconductors and insulators.

   // Data members:
   CADPhysicsDIMessenger* messenger;

   // Material-independent data
   G4DataVector Psecenergies;
      // Energy fraction values for tabulated data used for excitation of Fermi sea electron
   G4DataVector Psecvalues;
      // Tabulated function used for excitation of Fermi sea electron

   // Vectors that contain data for all materials (used during simulations)
   //CADPhysicsDataTable* outershelltable;
      // Contains for each material a list of outer shell energies
   std::vector<G4bool> outershells;
      // Per material, a boolean to state whether it has any outer shell energies
   std::vector<G4double> vec_fermieff;
      // Fermi energy per material
   std::vector<G4double> vec_bandgap;
      // Bandgap per material
   std::vector<G4double> vec_barrier;
      // Vacuum potential barrier per material
   std::vector<typename material::conductor_type_t> vec_conductortype;
      // Conductortype per material (0=metal, 1=semiconductor, 2=insulator)

   std::vector<imfp_table<float>> inelastic_imfp_vector;
   std::vector<icdf_table<float>> inelastic_icdf_vector;
   std::vector<imfp_table<float>> elastic_imfp_vector;
   std::vector<icdf_table<float>> elastic_icdf_vector;
   std::vector<ionization_table<float>> ionization_icdf_vector;
   std::vector<std::vector<float>> outershelltable;

   G4int    findex;
   G4double ffermienergy;
      // stored values of index and energy for use in GetMeanFreePath
   G4double fbarrier;
   //G4String    fconductortype;

   // Variables controlled by the Messenger (via their Set methods)
   G4bool generateSecondaries;
      // Determines whether secondaries are generated.
   G4bool generateXrays;
      // Determines whether X ray photons are generated.
   G4int pairsgenerated;
      // Counter for the number of secondary electrons created.

   G4bool killthisone;
      // Boolean set by GetMeanFreePath that tells PostStepDoIt to kill the
      // particle in case its energy has fallen below the energy limit.

   G4bool DI_output; // Boolean to set whether all absorbed electrons are
      // written to file

   // Specific for Auger/X-ray generation
   G4CrossSectionHandler* crossSectionHandler;
   CADPhysicsAtomicDeexcitation deexcitation;

   // for the output of the absorbed electrons
   std::ofstream output;

};

inline G4bool CADPhysicsDI::IsApplicable(const G4ParticleDefinition& particle)
// Applies to electrons
{
   return ( &particle == G4Electron::Electron() );
}

inline G4double CADPhysicsDI::PostStepGetPhysicalInteractionLength(
   const G4Track& track,
   G4double   previousStepSize,
   G4ForceCondition* condition
   )
{
   if ( (previousStepSize <=0.0) || (theNumberOfInteractionLengthLeft<=0.0)) {
      // We're at the beggining of tracking for this particle,
      // or just after the previous PostStepDoIt of this process.
      // The NumberOfInteractionLengthLeft needs to be (re)set.
      ResetNumberOfInteractionLengthLeft();
   } else {
      // Update NumberOfInteractionLengthLeft for the size of the previous step taken
      SubtractNumberOfInteractionLengthLeft(previousStepSize);
      if(theNumberOfInteractionLengthLeft<0.)
         theNumberOfInteractionLengthLeft=perMillion;
   }

   // Condition is set to "Not Forced", i.e. a call to PostStepDoIt is not normally enforced
   // for this step
   *condition = NotForced;

   // Get mean free path
   currentInteractionLength = GetMeanFreePath(track, previousStepSize, condition);

#ifdef G4VERBOSE
   if ((currentInteractionLength <=0.0) || (verboseLevel>4)){
      G4cout << "CADPhysicsDI::PostStepGetPhysicalInteractionLength ";
      G4cout << "[ " << GetProcessName() << "]" <<G4endl;
      track.GetDynamicParticle()->DumpInfo();
      G4cout << " in Material  " << track.GetMaterial()->GetName() <<G4endl;
      G4cout << "MeanFreePath = " << currentInteractionLength/cm << "[cm]" <<G4endl;
   }
#endif

   // Determine the actual step length for this step
   G4double value;
   if (currentInteractionLength <DBL_MAX) {
      value = theNumberOfInteractionLengthLeft * currentInteractionLength;
   } else {
      value = DBL_MAX;
   }
#ifdef G4VERBOSE
   if (verboseLevel>3){
      G4cout << "CADPhysicsDI::PostStepGetPhysicalInteractionLength ";
      G4cout << "[ " << GetProcessName() << "]" <<G4endl;
      track.GetDynamicParticle()->DumpInfo();
      G4cout << " in Material  " <<  track.GetMaterial()->GetName() <<G4endl;
      G4cout << "InteractionLength= " << value/cm <<"[cm] " <<G4endl;
   }
#endif
   return value;
}

inline G4double CADPhysicsDI::AlongStepGetPhysicalInteractionLength(
   const G4Track& track,
   G4double previousStepSize,
   G4double currentMinimumStep,
   G4double& currentSafety,
   G4GPILSelection* selection)
{
   // Get the step limit proposed by the continuous part of the process. N.B., for this process
   // no such step limit has been implemented (hence, returns DBL_MAX)
   G4double steplength = GetContinuousStepLimit(track,previousStepSize,
      currentMinimumStep,currentSafety);

   // Set return value for G4GPILSelection to the default value:
   *selection = NotCandidateForSelection;
   return  steplength;
}

#endif
