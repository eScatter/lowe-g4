// CADPhysicsTransportation.hh
//
// Process based on G4Transportation, Geant4 v. 8.0p1. Original Geant4 disclaimer follows:
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// This process has been modified to let it work together with the multistep and diffusionstep methods
// in the CADPhysicsSingleScattering process. To this end, the original AlongStepGetPhysicalInteractionLength
// method has been split into two: a new AlongStepGetPhysicalInteractionLength and a method RealASGPIL. The latter
// does the 'real' work, and is called by CADPhysicsSingleScattering *before* AlongStepGetPhysicalInteractionLength is called.
// AlongStepDoIt has also been modified accordingly.
//
// The reason for this way of working is that normally, G4Transportation is the last process to execute its
// AlongStepGetPhysicalInteractionLength method, but once it does so, it immediately proceeds to calculate the post-step parameters
// (position, direction, time, etc.) of the particle. While CADPhysicsSingleScattering needs to know the length of the geometry-limited step
// and the current safety, it also needs to be able to 'overrule' the geometry-limited step size of G4Transportation in the case of application of the
// 'multistep' or 'diffusionstep' model. Hence, its AlongStepPhysicalInteractionLength() is called by the SteppingManager *before* that of
// the Transportation process, but 'under water' CADPhysicsSingleScattering gets all the information from the Transportation process that its needs
// via the RealASGPIL method.
//
// See further comments marked 'EK' in this file and in the abovementioned methods in CADPhysicsTransportation.cc
//
#ifndef CADPhysicsTransportation_hh
#define CADPhysicsTransportation_hh 1

#include "G4VProcess.hh"
#include "G4FieldManager.hh"

#include "CADPhysicsUnits.hh"

#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4PropagatorInField.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleChangeForTransport.hh"
class G4SafetyHelper;

class CADPhysicsTransportation : public G4VProcess
{
   // Concrete class that does the geometrical transport

public:  // with description
   //constructor method hidden as private. GetInstance() method public instead:
   static CADPhysicsTransportation* GetInstance();

   ~CADPhysicsTransportation();

   G4double AlongStepGetPhysicalInteractionLength(
      const G4Track& track,
      G4double  previousStepSize,
      G4double  currentMinimumStep,
      G4double& currentSafety,
      G4GPILSelection* selection
      );

        // EK: Added
   G4double RealASGPIL(
      const G4Track& track,
      G4double  previousStepSize,
      G4double  currentMinimumStep,
      G4double& currentSafety,
      G4GPILSelection* selection
      );
        // EK: End of addition

   G4VParticleChange* AlongStepDoIt(
      const G4Track& track,
      const G4Step& stepData
      );

   G4VParticleChange* PostStepDoIt(
      const G4Track& track,
      const G4Step&  stepData
      );
   // Responsible for the relocation.

   G4double PostStepGetPhysicalInteractionLength(
      const G4Track& ,
      G4double   previousStepSize,
      G4ForceCondition* pForceCond
      );
   // Forces the PostStepDoIt action to be called,
   // but does not limit the step.

   G4PropagatorInField* GetPropagatorInField();
   void SetPropagatorInField( G4PropagatorInField* pFieldPropagator);
   // Access/set the assistant class that Propagate in a Field.

   inline void   SetVerboseLevel( G4int verboseLevel );
   inline G4int  GetVerboseLevel() const;
   // Level of warnings regarding eg energy conservation
   // in field integration.

   inline G4double GetThresholdWarningEnergy() const;
   inline G4double GetThresholdImportantEnergy() const;
   inline G4int GetThresholdTrials() const;

   inline void SetThresholdWarningEnergy( G4double newEnWarn );
   inline void SetThresholdImportantEnergy( G4double newEnImp );
   inline void SetThresholdTrials(G4int newMaxTrials );

   // Get/Set parameters for killing loopers:
   //   Above 'important' energy a 'looping' particle in field will
   //   *NOT* be abandoned, except after fThresholdTrials attempts.
   // Below Warning energy, no verbosity for looping particles is issued

   inline G4double GetMaxEnergyKilled() const;
   inline G4double GetSumEnergyKilled() const;
   inline void ResetKilledStatistics( G4int report = 1);
   // Statistics for tracks killed (currently due to looping in field)

   inline void EnableShortStepOptimisation(G4bool optimise=true);
   // Whether short steps < safety will avoid to call Navigator (if field=0)


public:  // without description

   G4double AtRestGetPhysicalInteractionLength(
      const G4Track& ,
      G4ForceCondition*
      ) { return -1.0; };
   // No operation in  AtRestDoIt.

   G4VParticleChange* AtRestDoIt(
      const G4Track& ,
      const G4Step&
      ) {return 0;};
   // No operation in  AtRestDoIt.

   void StartTracking(G4Track* aTrack);
   // Reset state for new (potentially resumed) track

protected:

   CADPhysicsTransportation( G4int verbosityLevel= 1);

   G4bool               DoesGlobalFieldExist();
   // Checks whether a field exists for the "global" field manager.

private:
   static CADPhysicsTransportation* fTransportationInstance;//the instance

   G4Navigator*         fLinearNavigator;
   G4PropagatorInField* fFieldPropagator;
   // The Propagators used to transport the particle

   // G4FieldManager*      fGlobalFieldMgr;     // Used MagneticField CC
   // Field Manager for the whole Detector

   G4ThreeVector        fTransportEndPosition;
   G4ThreeVector        fTransportEndMomentumDir;
   G4double             fTransportEndKineticEnergy;
   G4ThreeVector        fTransportEndSpin;
   G4bool               fMomentumChanged;
   G4bool               fEnergyChanged;
   G4bool               fEndGlobalTimeComputed;
   G4double             fCandidateEndGlobalTime;
   // The particle's state after this Step, Store for DoIt

   G4bool               fParticleIsLooping;

   G4TouchableHandle    fCurrentTouchableHandle;

   // G4bool         fFieldExists;
   // Whether a magnetic field exists ...
   // A data member for this is problematic: it is useful only if it
   // can be initialised and updated -- and a scheme is not yet possible.

   G4bool fGeometryLimitedStep;
   // Flag to determine whether a boundary was reached.

   G4ThreeVector  fPreviousSftOrigin;
   G4double       fPreviousSafety;
   // Remember last safety origin & value.

   G4ParticleChangeForTransport fParticleChange;
   // New ParticleChange

   G4double endpointDistance;

   // Thresholds for looping particles:
   //
   G4double fThreshold_Warning_Energy;     //  Warn above this energy
   G4double fThreshold_Important_Energy;   //  Hesitate above this
   G4int    fThresholdTrials;              //    for this no of trials
   // Above 'important' energy a 'looping' particle in field will
   //   *NOT* be abandoned, except after fThresholdTrials attempts.
   G4double fUnimportant_Energy;
   //  Below this energy, no verbosity for looping particles is issued

   // Counter for steps in which particle reports 'looping',
   //   if it is above 'Important' Energy
   G4int    fNoLooperTrials;
   // Statistics for tracks abandoned
   G4double fSumEnergyKilled;
   G4double fMaxEnergyKilled;

        // Whether to avoid calling G4Navigator for short step ( < safety)
        //   If using it, the safety estimate for endpoint will likely be smaller.
        G4bool   fShortStepOptimisation;

        G4SafetyHelper* fpSafetyHelper;  // To pass it the safety value obtained

   // Verbosity
   G4int    fVerboseLevel;
   // Verbosity level for warnings
   // eg about energy non-conservation in magnetic field.

   // EK: new variables for communication with CADPhysicsSingleScattering:
   G4bool fDoconventionalstep; // Boolean set by AlongStepGetInteractionLength().
      // Signals AlongStepDoIt() whether it should do a conventional step.

   // EK: the following are parameters set by AlongStepGetInteractionLength()
   //     that are needed by RealASGPIL()
   G4double geometryStepLength;// Physical step length as limited by the geometry
   G4double fSafety;
   G4GPILSelection* fselection;

   G4bool asgpildone;
       // This boolean tells AlongStepGetInteractionLength()
       // whether RealASGPIL() has been executed for this step. Since
       // RealASGPIL() is normally called by CADPhysicsSingleScattering for
       // every step, if it hasn't been executed, this means that apparently
       // CADPhysicsSingleScattering is not present or inactive. In this case,
       // the current run is aborted so that the user can fix the matter.

   // EK, end of addition

   // for the output of the absorbed electrons
   std::ofstream output;
};

#include "CADPhysicsTransportation.icc"

#endif
