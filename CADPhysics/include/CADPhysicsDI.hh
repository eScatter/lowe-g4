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

   inline void SetGenerateAugers(G4bool genauger) {
      // Method to determine whether Auger electrons should be generated.
      generateAugers = genauger;
      deexcitation.ActivateAugerElectronProduction(generateAugers);
   }

   inline G4bool GetGenerateAugers() {
      // Corresponding 'get' method (invoked by the Messenger)
      return generateAugers;
   }

   inline void SetRangeCut(G4bool rc) {
      // Method to determine whether 'range cut' simulation speedup method should be applied
      rangecut = rc;
   }
   inline G4bool GetRangeCut() {
      // Corresponding 'get' method (invoked by the Messenger)
      return rangecut;
   }

   inline void SetEnergyCut(G4double ec){
      // Method to set minimum energy for scattered and secondary electrons
      cutenergy = ec;
      if (ec>0.) energycut = true; else energycut = false;
   }

   inline G4double GetEnergyCut() {
      // Corresponding 'get' method (invoked by the Messenger)
      return cutenergy;
   }

   inline void ResetCounter() {
      // Method called from the Messenger for counting generated secondary electrons
      // (for debugging purposes)
      G4cout << "DI generated total of " << pairsgenerated << " pairs." << G4endl;
      pairsgenerated = 0;
   }

private:
   // Hiding the assignment operator and the copy constructor as private,
   // since it doesn't make sense to have multiple instances of this class
   CADPhysicsDI & operator = (const CADPhysicsDI &right);
   CADPhysicsDI ( const CADPhysicsDI &);

   // Private methods:
   void BuildEnrange();
      // Fill the enrange, Psecenergies and Psecvalues vectors. Invoked at initialization
      // (by the BuildPhysicsTable method).

   G4String GetdfFileForMaterial(const G4String &materialName);
      // Check existence of a df file with name "../CADPhysics/df/df_"+materialName. i
      // Invoked by ReadEpsilonFile.

   G4bool ExistsFileForMaterial(const G4String &materialName);
      // Check existence of a file with name materialName. Invoked by ExistsdfFileForMaterial.

   void ReadEpsilonFile(const G4String& materialName,
      const G4String& materialFormula,
      G4bool isgas,
      G4double barrier);
      // Read a material's inner shell energies and energy loss function from file

   void InterpolateEpsilon();
      // Interpolate the energy loss function to the 'standard' energy values stored in enrange

   void CalculateCDCS (CADPhysicsDataTable* cdcs,
      G4DataVector* csv,
      G4double fermiEnergy,
      G4double work);
      // Calculate the Cumulative Differential Cross Section and the electron's range as a function
      // of energy

   G4double GetRange(G4int index, G4double kineticEnergy);
      // Get the electron's range from the tabulated data

   void TabulateL();
      // Tabulate values of the 'L' and 'Lprim' functions that are used by CalculateCDCS
      // See Ashley eq. (20) the non-exchange corrected definition of L, and also see further
      // notes in the code.

   G4int FindIntLogen(G4double kineticEnergy);
      // Return the number of the energy bin that kineticEnergy is in

   G4VParticleChange* Phononloss(const G4Track& aTrack,
      const G4Step& aStep,
      G4double kinenergy,
      G4double omegaprime);
      // Process the event in case of sub-bandgap energy loss in semiconductors and insulators.

   // Data members:
   CADPhysicsDIMessenger* messenger;

   G4double LowestKineticEnergy;
      // Lowest energy value in enrange
   G4double HighestKineticEnergy;
      // Highest energy value in enrange
   G4int    TotBin;
      // Number of values in enrange

   // Material-independent data
   G4DataVector enrange;
      // Energy values for internally generated tabulated data
   G4DataVector Psecenergies;
      // Energy fraction values for tabulated data used for excitation of Fermi sea electron
   G4DataVector Psecvalues;
      // Tabulated function used for excitation of Fermi sea electron
   G4DataVector L;
      // 'L' function
   G4DataVector LPrim;
      // Idem, for inner-shell ionization events

   // Data for a specific material (used only during initialization)
   G4DataVector readepsenergies;
      // Energy values of the energy loss function file
   G4DataVector readepsdata;
      // Corresponding energy loss function values
   G4DataVector tmfpenergies;
      // Energy values of the 'tmfp' file
   G4DataVector tmfpdata;
      // Corresponding data values
   G4DataVector* imeps;
      // The energy loss function

   // Vectors that contain data for all materials (used during simulations)
   CADPhysicsDataTable* rangetable;
      // Contains for each material the electron range as a function of energy
   CADPhysicsDataTable* innershelltable;
      // Contains for each material a list of inner shell energies
   std::vector<G4bool> innershells;
      // Per material, a boolean to state whether it has any inner shell energies
   std::vector<G4double> vec_fermieff;
      // Fermi energy per material
   std::vector<G4double> vec_minimumlimit;
      // Lowest kinetic energy limit allowable for the material
   std::vector<G4double> vec_bandgap;
      // Bandgap per material
   std::vector<G4double> vec_barrier;
      // Vacuum potential barrier per material
   std::vector<G4int> vec_conductortype;
      // Conductortype per material (0=metal, 1=semiconductor, 2=insulator)
   CADPhysicsDataTable* lambdatable;
      // Table containing inverse mfp values
   std::vector<CADPhysicsDataTable*>* difflambdatable;
      // Differential (for omega prime) inverse mfp values

   G4double cross;
      // Cross section as calculated in GetMeanFreePath and reused in PostStepDoIt

   G4int    findex;
   G4double fenergy;
   G4double ffermienergy;
      // stored values of index, energy and lambda for use in GetMeanFreePath
   G4double fenergylimit;
      // Local value of the energy limit corrected for the material's barrier energy wrt vacuum
   G4double flambda;
   G4double fbarrier;
   G4int    fconductortype;

   // Variables controlled by the Messenger (via their Set methods)
   G4bool generateSecondaries;
      // Determines whether secondaries are generated.
   G4bool generateXrays;
      // Determines whether X ray photons are generated.
   G4bool generateAugers;
      // Determines whether Auger electrons are generated.
   G4bool rangecut;
      // Determines whether to apply the 'range cut' to speed up simulations.
   G4bool energycut;
      // Determines whether to apply an energy limit to the primary and secondary electrons.
   G4double cutenergy;
      // The corresponding value of the energy cut (in vacuum)
   G4int pairsgenerated;
      // Counter for the number of secondary electrons created.

   G4bool killthisone;
      // Boolean set by GetMeanFreePath that tells PostStepDoIt to kill the
      // particle in case its energy has fallen below the energy limit.

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
