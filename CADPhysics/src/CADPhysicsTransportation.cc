// CADPhysicsTransportation.cc
//
// See CADPhysicsTransportation.hh for class info and the standard Geant4 disclaimer.

#include "CADPhysicsTransportation.hh"
#include "CADPhysicsTransportationMessenger.hh"
#include "CADPhysicsSingleScattering.hh"
#include "CADPhysicsOutput.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ParticleTable.hh"
#include "G4ChordFinder.hh"
#include "G4SafetyHelper.hh"
#include "G4FieldManagerStore.hh"
#include "G4RunManager.hh"
#include "G4EquationOfMotion.hh"
class G4VSensitiveDetector;

CADPhysicsTransportation* CADPhysicsTransportation::fTransportationInstance = 0;

//////////////////////////////////////////////////////////////////////////
//
// Constructor

CADPhysicsTransportation::CADPhysicsTransportation( G4int verboseLevel )
  : G4VProcess( G4String("CADTrans"), fTransportation ),
    fParticleIsLooping( false ),
    fPreviousSftOrigin (0.,0.,0.),
    fPreviousSafety    ( 0.0 ),
    fThreshold_Warning_Energy( 100 * MeV ),
    fThreshold_Important_Energy( 250 * MeV ),
    fThresholdTrials( 10 ),
    fUnimportant_Energy( 1 * MeV ),
    fNoLooperTrials(0),
    fSumEnergyKilled( 0.0 ), fMaxEnergyKilled( 0.0 ),
    fShortStepOptimisation(false),    // Old default: true (=fast short steps)
    fVerboseLevel( verboseLevel )
{
  G4TransportationManager* transportMgr ;

  transportMgr = G4TransportationManager::GetTransportationManager() ;

  fLinearNavigator = transportMgr->GetNavigatorForTracking() ;

  // fGlobalFieldMgr = transportMgr->GetFieldManager() ;

  fFieldPropagator = transportMgr->GetPropagatorInField() ;

  fpSafetyHelper =   transportMgr->GetSafetyHelper();  // New

  // Cannot determine whether a field exists here,
  //  because it would only work if the field manager has informed
  //  about the detector's field before this transportation process
  //  is constructed.
  // Instead later the method DoesGlobalFieldExist() is called

  static G4TouchableHandle nullTouchableHandle;  // Points to (G4VTouchable*) 0
  fCurrentTouchableHandle = nullTouchableHandle;

  fEndGlobalTimeComputed  = false;
  fCandidateEndGlobalTime = 0;

  // EK: For communication with CADPhysicsSingleScattering (and variants)
  fDoconventionalstep = true;
  geometryStepLength = 0.;
  fSafety = 0.;
  fselection = new G4GPILSelection();
  asgpildone = false;
  // EK: end of addition
  transportation_output = false;
  messenger = new CADPhysicsTransportationMessenger(this);
}

//////////////////////////////////////////////////////////////////////////

CADPhysicsTransportation* CADPhysicsTransportation::GetInstance()
{
   static CADPhysicsTransportation theTransportation;
   if (!fTransportationInstance)
   {
      fTransportationInstance = &theTransportation;
   }
   return fTransportationInstance;
}

//////////////////////////////////////////////////////////////////////////

CADPhysicsTransportation::~CADPhysicsTransportation()
{
   if( (fVerboseLevel > 0) && (fSumEnergyKilled > 0.0 ) ){
      G4cout << " CADPhysicsTransportation: Statistics for looping particles " << G4endl;
      G4cout << "   Sum of energy of loopers killed: " <<  fSumEnergyKilled << G4endl;
      G4cout << "   Max energy of loopers killed: " <<  fMaxEnergyKilled << G4endl;
   }
   output.close();
}

//////////////////////////////////////////////////////////////////////////
//
// Responsibilities:
//    Find whether the geometry limits the Step, and to what length
//    Calculate the new value of the safety and return it.
//    Store the final time, position and momentum.

G4double CADPhysicsTransportation::
AlongStepGetPhysicalInteractionLength( const G4Track&,//  track,
                             G4double, //  previousStepSize
                             G4double,//  currentMinimumStep,
                             G4double& currentSafety,
                             G4GPILSelection* selection )
{
   // EK, AlongStepGetPhysicalInteractionLength() has been reduced to a routine
   // that passes the parameters set by RealASGPIL to the SteppingManager.
   // This routine assumes that RealASGPIL has already been executed for the
   // current step. If not, it issues a warning and aborts the run.
   if (!asgpildone) {
      G4cerr << "CADPhysicsTransportation: process 'ssc' (defined by CADPhysicsSingleScattering) appears to be inactive!" << G4endl;
      G4cerr << "You may correct this by executing the following command:" << G4endl;
      G4cerr << " /process/activate ssc" << G4endl;
      G4cerr << "Cannot continue now, current run will be aborted." << G4endl;
      G4RunManager::GetRunManager()->AbortRun(false);
   }
   asgpildone = false;
   currentSafety = fSafety;
   *selection = *fselection;
   return geometryStepLength;
}

G4double CADPhysicsTransportation::RealASGPIL( const G4Track&  track,
                                   G4double, //  previousStepSize
                                   G4double  currentMinimumStep,
                                   G4double& currentSafety,
                                   G4GPILSelection* selection )
{
   asgpildone = true;
   //G4double geometryStepLength, newSafety ;
   G4double newSafety;
   fParticleIsLooping = false ;

   // Initial actions moved to  StartTrack()
   // --------------------------------------
   // Note: in case another process changes touchable handle
   //    it will be necessary to add here (for all steps)
   // fCurrentTouchableHandle = aTrack->GetTouchableHandle();

   // GPILSelection is set to defaule value of CandidateForSelection
   // It is a return value
   //
   *selection = CandidateForSelection ;

   // Get initial Energy/Momentum of the track
   //
   const G4DynamicParticle*    pParticle  = track.GetDynamicParticle() ;
   const G4ParticleDefinition* pParticleDef   = pParticle->GetDefinition() ;
   G4ThreeVector startMomentumDir       = pParticle->GetMomentumDirection() ;
   G4ThreeVector startPosition          = track.GetPosition() ;

   // G4double   theTime        = track.GetGlobalTime() ;

   // The Step Point safety can be limited by other geometries and/or the
   // assumptions of any process - it's not always the geometrical safety.
   // We calculate the starting point's isotropic safety here.
   //
   G4ThreeVector OriginShift = startPosition - fPreviousSftOrigin ;
   G4double      MagSqShift  = OriginShift.mag2() ;
   if( MagSqShift >= sqr(fPreviousSafety) )
   {
      currentSafety = 0.0 ;
   }
   else
   {
      currentSafety = fPreviousSafety - std::sqrt(MagSqShift) ;
   }

   // EK, next is the part that has been modified form standard G4Transportation.
   fDoconventionalstep = true;// Default is that a 'conventional' step will be done

   G4bool multistepping = false;
   G4double elsl = CADPhysicsSingleScattering::GetInstance()->GetFirstsinglestep(multistepping);
   // The above call does two things:
   // 1) CADPhysicsSingleScattering returns the size of the step before
   //    the *first* elastic scatter event
   // 2) The boolean multistepping is set to the value that reflects whether
   //    CADPhysicsSingleScattering has been set up to do multistepping. If
   //    not, then we are certain to proceed with a conventional step.
   if (elsl < currentMinimumStep) {
         // The elastic scattering interaction length is
         // smaller than any length defined by another process. If not,
         // then another process wins and we proceed to do a conventional step.
      currentMinimumStep = elsl;// Update the current minimum step for the interaction
         // length indicated by SingleScattering

      if (multistepping) {// CADPhysicsSingleScattering was told by the user that it can do multistepping
         if (elsl<currentSafety) {// Its interaction length is smaller than the
               // current safety, so Transportation knows that it should not do
               // a conventional step, and SingleScattering will perform either
               // a multi- or a diffusionstep.
            *selection = NotCandidateForSelection ;
            *fselection = *selection;
            fDoconventionalstep = false;
            fSafety = currentSafety;
            geometryStepLength = DBL_MAX;
            return DBL_MAX;
         }

         // The currently known safety is smaller than the next step length.
         // However, we can recalculate the safety and make a second attempt.
         G4double linearnewsafety = fLinearNavigator->ComputeSafety(startPosition);
         if (elsl<linearnewsafety) {// See above
            currentSafety = linearnewsafety ;
               // The newly calculated safety is stored to currentSafety for
               // the benefit of SingleScattering
            *selection = NotCandidateForSelection ;
            *fselection = *selection;
            fDoconventionalstep = false;
            fSafety = currentSafety;
            geometryStepLength = DBL_MAX;
            return DBL_MAX;
         }
      }
      // If we get here, then SingleScattering was told not to do multistepping and we proceed to do a conventional step.
   }
   CADPhysicsSingleScattering::GetInstance()->NoRaytrace();
      // Since Transportation does a conventional step, CADPhysicsSingleScattering
      // is told not to do any raytracing (i.e. no multi- or diffusionstep). It may
      // still do a single elastic scatter event, depending on its interaction
      // length compared to other processes.

   // EK, end of this modified part

   // Is the particle charged ?
   //
   G4double              particleCharge = pParticle->GetCharge() ;
   G4double              magneticMoment = pParticle->GetMagneticMoment() ;

   fGeometryLimitedStep = false ;
   // fEndGlobalTimeComputed = false ;

   // There is no need to locate the current volume. It is Done elsewhere:
   //   On track construction
   //   By the tracking, after all AlongStepDoIts, in "Relocation"

   // Check whether the particle have an (EM) field force exerting upon it
   //
   G4FieldManager* fieldMgr=0;
   G4bool          fieldExertsForce = false ;
   if( (particleCharge != 0.0) )
   {
      fieldMgr= fFieldPropagator->FindAndSetFieldManager( track.GetVolume() );
      if (fieldMgr != 0) {
        // Message the field Manager, to configure it for this track
        fieldMgr->ConfigureForTrack( &track );
        // Moved here, in order to allow a transition
        //   from a zero-field  status (with fieldMgr->(field)0
        //   to a finite field  status

         // If the field manager has no field, there is no field !
         fieldExertsForce = (fieldMgr->GetDetectorField() != 0);
      }
   }

   // Choose the calculation of the transportation: Field or not
   //
   if( !fieldExertsForce )
   {
      G4double linearStepLength ;
      if( fShortStepOptimisation && (currentMinimumStep <= currentSafety) )
      {
         // The Step is guaranteed to be taken
         //
         geometryStepLength   = currentMinimumStep ;
         fGeometryLimitedStep = false ;
      }
      else
      {
         //  Find whether the straight path intersects a volume
         //
         linearStepLength = fLinearNavigator->ComputeStep( startPosition,
                                                           startMomentumDir,
                                                           currentMinimumStep,
                                                           newSafety) ;
         // Remember last safety origin & value.
         //
         fPreviousSftOrigin = startPosition ;
         fPreviousSafety    = newSafety ;

         // The safety at the initial point has been re-calculated:
         //
         currentSafety = newSafety ;

         fGeometryLimitedStep= (linearStepLength <= currentMinimumStep);
         if( fGeometryLimitedStep )
         {
            // The geometry limits the Step size (an intersection was found.)
            //
            geometryStepLength   = linearStepLength ;
         }
         else
         {
            // The full Step is taken.
            //
            geometryStepLength   = currentMinimumStep ;
         }
      }
      endpointDistance = geometryStepLength ;

      // Calculate final position
      //
      fTransportEndPosition = startPosition+geometryStepLength*startMomentumDir ;

      // Momentum direction, energy and polarisation are unchanged by transport
      //
      fTransportEndMomentumDir   = startMomentumDir ;
      fTransportEndKineticEnergy = track.GetKineticEnergy() ;
      fTransportEndSpin          = track.GetPolarization();
      fParticleIsLooping         = false ;
      fMomentumChanged           = false ;
      fEndGlobalTimeComputed     = false ;
   }
   else   //  A field exerts force
   {
      G4double       momentumMagnitude = pParticle->GetTotalMomentum() ;
      G4ThreeVector  EndUnitMomentum ;
      G4double       lengthAlongCurve ;
      G4double       restMass = pParticleDef->GetPDGMass() ;

      G4ChargeState chargeState(particleCharge,             // The charge can change (dynamic)
                                magneticMoment,
                                pParticleDef->GetPDGSpin() );

      G4EquationOfMotion* equationOfMotion =
      (fFieldPropagator->GetChordFinder()->GetIntegrationDriver()->GetStepper())
      ->GetEquationOfMotion();

      equationOfMotion->SetChargeMomentumMass( chargeState,    // in e+ units
                                               momentumMagnitude, // in Mev/c
                                               restMass           ) ;

      G4ThreeVector spin        = track.GetPolarization() ;
      G4FieldTrack  aFieldTrack = G4FieldTrack( startPosition,
                                                track.GetMomentumDirection(),
                                                0.0,
                                                track.GetKineticEnergy(),
                                                restMass,
                                                track.GetVelocity(),
                                                track.GetGlobalTime(), // Lab.
                                                track.GetProperTime(), // Part.
                                                &spin                  ) ;
      if( currentMinimumStep > 0 )
      {
         // Do the Transport in the field (non recti-linear)
         //
         lengthAlongCurve = fFieldPropagator->ComputeStep( aFieldTrack,
                                                           currentMinimumStep,
                                                           currentSafety,
                                                           track.GetVolume() ) ;
         fGeometryLimitedStep= lengthAlongCurve < currentMinimumStep;
         if( fGeometryLimitedStep ) {
            geometryStepLength   = lengthAlongCurve ;
         } else {
            geometryStepLength   = currentMinimumStep ;
         }
      }
      else
      {
         geometryStepLength   = lengthAlongCurve= 0.0 ;
         fGeometryLimitedStep = false ;
      }

      // Remember last safety origin & value.
      //
      fPreviousSftOrigin = startPosition ;
      fPreviousSafety    = currentSafety ;
      // fpSafetyHelper->SetCurrentSafety( newSafety, startPosition);

      // Get the End-Position and End-Momentum (Dir-ection)
      //
      fTransportEndPosition = aFieldTrack.GetPosition() ;

      // Momentum:  Magnitude and direction can be changed too now ...
      //
      fMomentumChanged         = true ;
      fTransportEndMomentumDir = aFieldTrack.GetMomentumDir() ;

      fTransportEndKineticEnergy  = aFieldTrack.GetKineticEnergy() ;

      if( fFieldPropagator->GetCurrentFieldManager()->DoesFieldChangeEnergy() )
      {
         // If the field can change energy, then the time must be integrated
         //    - so this should have been updated
         //
         fCandidateEndGlobalTime   = aFieldTrack.GetLabTimeOfFlight();
         fEndGlobalTimeComputed    = true;

         // was ( fCandidateEndGlobalTime != track.GetGlobalTime() );
         // a cleaner way is to have FieldTrack knowing whether time is updated.
      }
      else
      {
         // The energy should be unchanged by field transport,
         //    - so the time changed will be calculated elsewhere
         //
         fEndGlobalTimeComputed = false;

         // Check that the integration preserved the energy
         //     -  and if not correct this!
         G4double  startEnergy= track.GetKineticEnergy();
         G4double  endEnergy= fTransportEndKineticEnergy;

         static G4int no_inexact_steps=0, no_large_ediff;
         G4double absEdiff = std::fabs(startEnergy- endEnergy);
         if( absEdiff > perMillion * endEnergy )
         {
            no_inexact_steps++;
            // Possible statistics keeping here ...
         }
         if( fVerboseLevel > 1 )
         {
            if( std::fabs(startEnergy- endEnergy) > perThousand * endEnergy )
            {
               static G4int no_warnings= 0, warnModulo=1,  moduloFactor= 10;
               no_large_ediff ++;
               if( (no_large_ediff% warnModulo) == 0 )
               {
                  no_warnings++;
                  G4cout << "WARNING - CADPhysicsTransportation::RealASGPIL() "
                         << "   Energy change in Step is above 1^-3 relative value. " << G4endl
                         << "   Relative change in 'tracking' step = "
                         << std::setw(15) << (endEnergy-startEnergy)/startEnergy << G4endl
                         << "     Starting E= " << std::setw(12) << startEnergy / MeV << " MeV " << G4endl
                         << "     Ending   E= " << std::setw(12) << endEnergy   / MeV << " MeV " << G4endl;
                  G4cout << " Energy has been corrected -- however, review"
                         << " field propagation parameters for accuracy."  << G4endl;
                  if( (fVerboseLevel > 2 ) || (no_warnings<4) || (no_large_ediff == warnModulo * moduloFactor) ){
                    G4cout << " These include EpsilonStepMax(/Min) in G4FieldManager "
                           << " which determine fractional error per step for integrated quantities. " << G4endl
                           << " Note also the influence of the permitted number of integration steps."
                           << G4endl;
                  }
                  G4cerr << "ERROR - G4Transportation::AlongStepGetPIL()" << G4endl
                         << "        Bad 'endpoint'. Energy change detected"
                         << " and corrected. "
                         << " Has occurred already "
                         << no_large_ediff << " times." << G4endl;
                  if( no_large_ediff == warnModulo * moduloFactor )
                  {
                     warnModulo *= moduloFactor;
                  }
               }
            }
         }  // end of if (fVerboseLevel)

         // Correct the energy for fields that conserve it
         //  This - hides the integration error
         //       - but gives a better physical answer
         fTransportEndKineticEnergy= track.GetKineticEnergy();
      }

      fTransportEndSpin = aFieldTrack.GetSpin();
      fParticleIsLooping = fFieldPropagator->IsParticleLooping() ;
      endpointDistance   = (fTransportEndPosition - startPosition).mag() ;
   }

   // If we are asked to go a step length of 0, and we are on a boundary
   // then a boundary will also limit the step -> we must flag this.
   //
   if( currentMinimumStep == 0.0 )
   {
      if( currentSafety == 0.0 )  fGeometryLimitedStep = true ;
   }

   // Update the safety starting from the end-point,
   // if it will become negative at the end-point.
   //
   if( currentSafety < endpointDistance )
   {
      // if( particleCharge == 0.0 )
      //    G4cout  << "  Avoiding call to ComputeSafety : charge = 0.0 " << G4endl;

      if( particleCharge != 0.0 ) {

         G4double endSafety =
            fLinearNavigator->ComputeSafety( fTransportEndPosition) ;

         // EK, we have a new safety based upon the conventional step taken. This new
         // safety should *not* be passed on to SingleScattering.
         // Therefore the following line has been commented out:
         //currentSafety      = endSafety ;
         fPreviousSftOrigin = fTransportEndPosition ;
         //fPreviousSafety    = currentSafety ;
         fPreviousSafety = endSafety;// EK, effectively does the same thing
         fpSafetyHelper->SetCurrentSafety( endSafety, fTransportEndPosition);
            // EB: Original for 4.9.2:
            // EB: fpSafetyHelper->SetCurrentSafety( currentSafety, fTransportEndPosition);
            // EB: It is meant to optimize calculation of safety.
            // EB: As we are not copying endSafety into currentSafety the call is changed

         // Because the Stepping Manager assumes it is from the start point,
         //  add the StepLength
         //
         // EK, also commented out:
         //currentSafety     += endpointDistance ;

#ifdef G4DEBUG_TRANSPORT
         G4cout.precision(12) ;
         G4cout << "***CADPhysicsTransportation::RealASGPIL ** " << G4endl  ;
         G4cout << "  Called Navigator->ComputeSafety at " << fTransportEndPosition
                << "    and it returned safety= " << endSafety << G4endl ;
         G4cout << "  Adding endpoint distance " << endpointDistance
                << "   to obtain pseudo-safety= " << currentSafety << G4endl ;
#endif
      }
   }

   fParticleChange.ProposeTrueStepLength(geometryStepLength) ;

   // EK: communication with CADPhysicsSingleStepping (and variants)A
   // EB: one would expect that G4SafetyHelper could be used but it has no
   // EB: GetCurrentSafety routine, only a ComputeSafety which may do too much.
   fSafety = currentSafety;
   *fselection = *selection;
   // EK: end of additions
   return geometryStepLength ;
}

//////////////////////////////////////////////////////////////////////////
//
//   Initialize ParticleChange  (by setting all its members equal
//                               to corresponding members in G4Track)

G4VParticleChange* CADPhysicsTransportation::AlongStepDoIt( const G4Track& track,
                                             const G4Step&  stepData )
{
   static G4int noCalls=0;
   static const G4ParticleDefinition* fOpticalPhoton =
            G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");

   noCalls++;

   fParticleChange.Initialize(track) ;

   // EK, check if we're doing a 'conventional' step:
   if (!fDoconventionalstep)
   {
      // No conventional step. However, because this process is the first in line for the AlongStepDoIt calls,
      // the ParticleChange does need to be initalized. We simply enter all the particle properties at the beginning of the step.
      const G4DynamicParticle*  pParticle = track.GetDynamicParticle();
      fParticleChange.ProposeEnergy(pParticle->GetKineticEnergy()) ;
      fParticleChange.ProposePosition(track.GetPosition()) ;
      fParticleChange.ProposeMomentumDirection(pParticle->GetMomentumDirection()) ;
      fParticleChange.SetMomentumChanged(false) ;
      fParticleChange.ProposeProperTime(track.GetProperTime()) ;
      fParticleChange.ProposeGlobalTime(track.GetGlobalTime());
      return &fParticleChange;
   }
   // EK, Conventional step: proceed with proposing the usual step parameters

   //  Code for specific process
   //
   fParticleChange.ProposePosition(fTransportEndPosition) ;
   fParticleChange.ProposeMomentumDirection(fTransportEndMomentumDir) ;
   fParticleChange.ProposeEnergy(fTransportEndKineticEnergy) ;
   fParticleChange.SetMomentumChanged(fMomentumChanged) ;

   fParticleChange.ProposePolarization(fTransportEndSpin);

   G4double deltaTime = 0.0 ;

   // Calculate  Lab Time of Flight (ONLY if field Equations used it!)
   // G4double endTime   = fCandidateEndGlobalTime;
   // G4double delta_time = endTime - startTime;

   G4double startTime = track.GetGlobalTime() ;

   if (!fEndGlobalTimeComputed)
   {
      // The time was not integrated .. make the best estimate possible
      //
      G4double finalVelocity   = track.GetVelocity() ;
      G4double initialVelocity = stepData.GetPreStepPoint()->GetVelocity() ;
      G4double stepLength      = track.GetStepLength() ;

      deltaTime= 0.0; // in case initialVelocity = 0
      const G4DynamicParticle* fpDynamicParticle = track.GetDynamicParticle();
      if (fpDynamicParticle->GetDefinition()== fOpticalPhoton)
      {
         //  A photon is in the medium of the final point
         //  during the step, so it has the final velocity.
         deltaTime = stepLength/finalVelocity ;
      }
      else if (finalVelocity > 0.0)
      {
         G4double meanInverseVelocity ;
         // deltaTime = stepLength/finalVelocity ;
         meanInverseVelocity = 0.5
            * ( 1.0 / initialVelocity + 1.0 / finalVelocity ) ;
         deltaTime = stepLength * meanInverseVelocity ;
      }
      else if (initialVelocity > 0.0 )
      {
         deltaTime = stepLength/initialVelocity ;
      }
      fCandidateEndGlobalTime   = startTime + deltaTime ;
   }
   else
   {
      deltaTime = fCandidateEndGlobalTime - startTime ;
   }

   fParticleChange.ProposeGlobalTime( fCandidateEndGlobalTime ) ;

   // Now Correct by Lorentz factor to get "proper" deltaTime

   G4double  restMass       = track.GetDynamicParticle()->GetMass() ;
   G4double deltaProperTime = deltaTime*( restMass/track.GetTotalEnergy() ) ;

   fParticleChange.ProposeProperTime(track.GetProperTime() + deltaProperTime) ;
   //fParticleChange. ProposeTrueStepLength( track.GetStepLength() ) ;

   // If the particle is caught looping or is stuck (in very difficult
   // boundaries) in a magnetic field (doing many steps)
   //   THEN this kills it ...
   //
   if ( fParticleIsLooping )
   {
      G4double endEnergy= fTransportEndKineticEnergy;

      if( (endEnergy < fThreshold_Important_Energy)
         || (fNoLooperTrials >= fThresholdTrials ) ){
         // Kill the looping particle
         //
         if(transportation_output) {
           if(!output.is_open()) {
             output.open ("Transportationout.dat");
             output << std::setprecision(10);
             output << "EventID\tTrackID\tEtotal (eV)\tEkin (eV)"
                    << "\tx (nm)\ty (nm)\tz (nm)\txOrigin (nm)\tyOrigin (nm)\tzOrigin (nm)\tx-direction (nm)\ty-direction (nm)\tz-direction (nm)"
                    << "\tMax Depth (nm)\tMax Radius (nm)" << G4endl;
           }
           output << CADPhysicsOutput::CADPhysicsOutput(track, stepData) << G4endl;
         }
         fParticleChange.ProposeTrackStatus( fStopAndKill )  ;

         // 'Bare' statistics
         fSumEnergyKilled += endEnergy;
         if( endEnergy > fMaxEnergyKilled) { fMaxEnergyKilled= endEnergy; }

#ifdef G4VERBOSE
         if( (fVerboseLevel > 1) ||
            ( endEnergy > fThreshold_Warning_Energy )  ) {
               G4cout << " CADPhysicsTransportation is killing track that is looping or stuck "
                  << G4endl
                  << "   This track has " << track.GetKineticEnergy() / MeV
                  << " MeV energy." << G4endl;
               G4cout << "   Number of trials = " << fNoLooperTrials
              << "   No of calls to AlongStepDoIt = " << noCalls
              << G4endl;
         }
#endif
         fNoLooperTrials=0;
      }
      else{
         fNoLooperTrials ++;
#ifdef G4VERBOSE
        if( (fVerboseLevel > 2) ){
          G4cout << "   CADPhysicsTransportation::AlongStepDoIt(): Particle looping -  "
                 << "   Number of trials = " << fNoLooperTrials
                 << "   No of calls to  = " << noCalls
                 << G4endl;
        }
#endif
      }
   }else{
      fNoLooperTrials=0;
   }

   // Another (sometimes better way) is to use a user-limit maximum Step size
   // to alleviate this problem ..

   // Introduce smooth curved trajectories to particle-change
   //
   fParticleChange.SetPointerToVectorOfAuxiliaryPoints
      (fFieldPropagator->GimmeTrajectoryVectorAndForgetIt() );

   return &fParticleChange ;
}

//////////////////////////////////////////////////////////////////////////
//
//  This ensures that the PostStep action is always called,
//  so that it can do the relocation if it is needed.
//

G4double CADPhysicsTransportation::
PostStepGetPhysicalInteractionLength( const G4Track&,
                            G4double, // previousStepSize
                            G4ForceCondition* pForceCond )
{
   *pForceCond = Forced ;
   return DBL_MAX ;  // was kInfinity ; but convention now is DBL_MAX
}

/////////////////////////////////////////////////////////////////////////////
//

G4VParticleChange* CADPhysicsTransportation::PostStepDoIt( const G4Track& track,
                                            const G4Step& step)
{
   G4TouchableHandle retCurrentTouchable ;   // The one to return

   // Initialize ParticleChange  (by setting all its members equal
   //                             to corresponding members in G4Track)
   // fParticleChange.Initialize(track) ;  // To initialise TouchableChange

   fParticleChange.ProposeTrackStatus(track.GetTrackStatus()) ;

   // If the Step was determined by the volume boundary,
   // logically relocate the particle

   if(fGeometryLimitedStep)
   {
      // fCurrentTouchable will now become the previous touchable,
      // and what was the previous will be freed.
      // (Needed because the preStepPoint can point to the previous touchable)

      fLinearNavigator->SetGeometricallyLimitedStep() ;
      fLinearNavigator->
         LocateGlobalPointAndUpdateTouchableHandle( track.GetPosition(),
         track.GetMomentumDirection(),
         fCurrentTouchableHandle,
         true                      ) ;
      // Check whether the particle is out of the world volume
      // If so it has exited and must be killed.
      //
      if( fCurrentTouchableHandle->GetVolume() == 0 )
      {
         if(transportation_output) {
           if(!output.is_open()) {
             output.open ("Transportationout.dat");
             output << std::setprecision(10);
             output << "EventID\tTrackID\tEtotal (eV)\tEkin (eV)"
                    << "\tx (nm)\ty (nm)\tz (nm)\txOrigin (nm)\tyOrigin (nm)\tzOrigin (nm)\tx-direction (nm)\ty-direction (nm)\tz-direction (nm)"
                    << "\tMax Depth (nm)\tMax Radius (nm)" << G4endl;
           }
           output << CADPhysicsOutput::CADPhysicsOutput(track, step) << G4endl;
         }
         fParticleChange.ProposeTrackStatus( fStopAndKill ) ;
      }
      retCurrentTouchable = fCurrentTouchableHandle ;
      fParticleChange.SetTouchableHandle( fCurrentTouchableHandle ) ;
   }
   else                 // fGeometryLimitedStep  is false
   {
      // This serves only to move the Navigator's location
      //
      fLinearNavigator->LocateGlobalPointWithinVolume( track.GetPosition() ) ;

      // The value of the track's current Touchable is retained.
      // (and it must be correct because we must use it below to
      // overwrite the (unset) one in particle change)
      // Although in general this is fCurrentTouchable, at the start of
      // a step it could be different ... ??
      //
      fParticleChange.SetTouchableHandle( track.GetTouchableHandle() ) ;
      retCurrentTouchable = track.GetTouchableHandle() ;
   }         // endif ( fGeometryLimitedStep )

   const G4VPhysicalVolume* pNewVol = retCurrentTouchable->GetVolume() ;
   const G4Material* pNewMaterial   = 0 ;
   const G4VSensitiveDetector* pNewSensitiveDetector   = 0 ;

   if( pNewVol != 0 )
   {
      pNewMaterial= pNewVol->GetLogicalVolume()->GetMaterial();
      pNewSensitiveDetector= pNewVol->GetLogicalVolume()->GetSensitiveDetector();
   }

   // ( <const_cast> pNewMaterial ) ;
   // ( <const_cast> pNewSensitiveDetector) ;

   fParticleChange.SetMaterialInTouchable( (G4Material *) pNewMaterial ) ;
   fParticleChange.SetSensitiveDetectorInTouchable( (G4VSensitiveDetector *) pNewSensitiveDetector ) ;

   const G4MaterialCutsCouple* pNewMaterialCutsCouple = 0;
   if( pNewVol != 0 )
   {
      pNewMaterialCutsCouple=pNewVol->GetLogicalVolume()->GetMaterialCutsCouple();
   }

   if( pNewVol!=0 && pNewMaterialCutsCouple!=0 && pNewMaterialCutsCouple->GetMaterial()!=pNewMaterial )
   {
      // for parametrized volume
      //
      pNewMaterialCutsCouple =
         G4ProductionCutsTable::GetProductionCutsTable()
         ->GetMaterialCutsCouple(pNewMaterial,
         pNewMaterialCutsCouple->GetProductionCuts());
   }
   fParticleChange.SetMaterialCutsCoupleInTouchable( pNewMaterialCutsCouple );

   // temporarily until Get/Set Material of ParticleChange,
   // and StepPoint can be made const.
   // Set the touchable in ParticleChange
   // this must always be done because the particle change always
   // uses this value to overwrite the current touchable pointer.
   //
   fParticleChange.SetTouchableHandle(retCurrentTouchable) ;


   return &fParticleChange ;
}

// New method takes over the responsibility to reset the state of CADPhysicsTransportation
//   object at the start of a new track or the resumption of a suspended track.

void
CADPhysicsTransportation::StartTracking(G4Track* aTrack)
{
   G4VProcess::StartTracking(aTrack);

   // The actions here are those that were taken in AlongStepGPIL
   //   when track.GetCurrentStepNumber()==1

   // reset safety value and center
   //
   fPreviousSafety    = 0.0 ;
   fPreviousSftOrigin = G4ThreeVector(0.,0.,0.) ;

   // reset looping counter -- for motion in field
   fNoLooperTrials= 0;
   // Must clear this state .. else it depends on last track's value
   //  --> a better solution would set this from state of suspended track TODO ?
   // Was if( aTrack->GetCurrentStepNumber()==1 ) { .. }

   // ChordFinder reset internal state
   //
   if( DoesGlobalFieldExist() ) {
      fFieldPropagator->ClearPropagatorState();
        // Resets all state of field propagator class (ONLY)
        //  including safety values (in case of overlaps and to wipe for first track).

      // G4ChordFinder* chordF= fFieldPropagator->GetChordFinder();
      // if( chordF ) chordF->ResetStepEstimate();
   }

   // Make sure to clear the chord finders of all fields (ie managers)
   static G4FieldManagerStore* fieldMgrStore= G4FieldManagerStore::GetInstance();
   fieldMgrStore->ClearAllChordFindersState();

   // Update the current touchable handle  (from the track's)
   //
   fCurrentTouchableHandle = aTrack->GetTouchableHandle();
}
