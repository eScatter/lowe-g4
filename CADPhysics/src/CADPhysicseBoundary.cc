// CADPhysicseBoundary.cc, loosely based on G4OpBoundaryProcess.cc
//
// 2006-07-11: sign of deltaphi2 corrected in definition of normalEnergy in DielectricDielectric
// 2007-07-19: added energy dependent check for absorption vs. reflection at interfaces
// 2007-11-15: correct bookkeeping of deposited energy for particles absorbed at boundary
// 2012-03-08: modification added. See notes in CADPhysicseBoundary.cc.120308 and notebook dd april 6 2011 EB

#include "CADPhysicseBoundary.hh"
#include "Randomize.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4TransportationManager.hh"


CADPhysicseBoundary::CADPhysicseBoundary(const G4String& processName,
                               G4ProcessType type)
                               : G4VDiscreteProcess(processName, type)
{
   if ( verboseLevel > 0) {
      G4cout << GetProcessName() << " is created " << G4endl;
   }
   theStatus = Undefined;
   lasttrackID = -1;
}

CADPhysicseBoundary::~CADPhysicseBoundary(){}

G4VParticleChange* CADPhysicseBoundary::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
   aParticleChange.Initialize(aTrack);
   G4int trackID = aTrack.GetTrackID();

   G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();

   if (pPostStepPoint->GetStepStatus() != fGeomBoundary) {
      // The particle is not at a boundary so nothing to be done
      theStatus = NotAtBoundary;
      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
   }
   // EB Added for debugging
   if (verboseLevel>1) {
      G4StepPoint* tmpPreStepPoint = aStep.GetPreStepPoint(); // Prevent unnecessary evaluation if verboseLevel=0
      G4cerr << "CADPhysicseBoundary::PostStepDoIt          " << tmpPreStepPoint->GetPhysicalVolume()->GetName() << " " << pPostStepPoint->GetPhysicalVolume()->GetName() << G4endl;
      if (verboseLevel>3) {
         Material1 = tmpPreStepPoint  -> GetMaterial();
         Material2 = pPostStepPoint -> GetMaterial();
         G4cout << "Pre-step  point =  " << tmpPreStepPoint->GetPosition() << " Material = " << Material1->GetName() << G4endl;
         G4cout << "Post-step point =  " << pPostStepPoint->GetPosition() << " Material = " << Material2->GetName() << G4endl;
         G4cout << " Status from previous step: ";
         if ( theStatus == Undefined )               G4cout << " Undefined " << G4endl;
         if ( theStatus == Refraction )              G4cout << " Refraction " << G4endl;
         if ( theStatus == TotalInternalReflection ) G4cout << " TotalInternalReflection " << G4endl;
         if ( theStatus == NotAtBoundary )           G4cout << " NotAtBoundary " << G4endl;
         if ( theStatus == SameMaterial )            G4cout << " SameMaterial " << G4endl;
         if ( theStatus == StepTooSmall )            G4cout << " StepTooSmall " << G4endl;
         if ( theStatus == ParticleKilled )          G4cout << " ParticleKilled " << G4endl;
      }
   }
   // EB End of addition

   // If the last step was a TotalInternalReflection, the momentum direction has changed. The consequence is that
   // it is now pointing back into the volume we came from before this step. The normal has now changed to point
   // into the new volume (i.e. after the first step). This means that the particle is "moving OUT" while it should
   // be moving in. This causes an oscillation that will be detected by the G4Navigator costing a lot of time.
   // Instead, we detect that the last track was TotalInternalReflection and return immediately.
   if ( theStatus == TotalInternalReflection && trackID == lasttrackID)
      // Prevent the particle from oscillating back and forth at the geom. boundary
      // If the particle was reflected at a boundary in the previous step, it will not be reflected again.
   {
      theStatus = NotAtBoundary;
      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
   }


   // Detect whether we are moving into or out of the volume. In the latter case we can return immediately
   G4bool        valid;
   G4Navigator*  theNavigator = G4TransportationManager::GetTransportationManager()-> GetNavigatorForTracking();
   G4ThreeVector theGlobalPoint = pPostStepPoint->GetPosition();
   G4ThreeVector theLocalPoint = theNavigator->GetGlobalToLocalTransform().TransformPoint(theGlobalPoint);
   // ********************************************************************
   // GetLocalExitNormal
   //
   // Obtains the Normal vector to a surface (in local coordinates)
   // pointing out of previous volume and into current volume
   // ********************************************************************
   G4ThreeVector theLocalNormal = theNavigator->GetLocalExitNormal(&valid);
   if (valid) {
      // Make sure the normal points back into volume
      theLocalNormal = -theLocalNormal;
   } else {
      // This happens if the particle is not at a boundary so nothing to be done
      // Observation: this happens if the track is curved due to fields
      if (verboseLevel > 1) {
         G4cerr << " CADPhysicseBoundary/PostStepDoIt(): " << " The Navigator reports that it returned an invalid normal" << G4endl;
         if (verboseLevel > 3) {
            theGlobalNormal = theNavigator->GetLocalToGlobalTransform().TransformAxis(theLocalNormal);
            G4cout << "theLocalPoint=" << theLocalPoint << G4endl;
            G4cout << "theGlobalPoint=" << theGlobalPoint << G4endl;
            G4cout << "theLocalNormal=" << theLocalNormal << G4endl;
            G4cout << "theGlobalNormal=" << theGlobalNormal << G4endl;
         }
      }
      theStatus = NotAtBoundary;
      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
   }
   theGlobalNormal = theNavigator->GetLocalToGlobalTransform().TransformAxis(theLocalNormal);
   if ( verboseLevel > 3 ) {
      G4cout << "theLocalPoint=" << theLocalPoint << G4endl;
      G4cout << "theGlobalPoint=" << theGlobalPoint << G4endl;
      G4cout << "theLocalNormal=" << theLocalNormal << G4endl;
      G4cout << "theGlobalNormal=" << theGlobalNormal << G4endl;
   }

   const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
   OldMomentum = aParticle->GetMomentumDirection();

   if (OldMomentum * theGlobalNormal > 0.0) {
      theStatus = TotalInternalReflection;
      if (verboseLevel>3) {
         G4cout << "Particle moving INTO volume, boundary process skipped." << G4endl;
      }
      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
   }

   G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
   theStatus = Undefined;
   lasttrackID = trackID;

   Material1 = pPreStepPoint  -> GetMaterial();
   Material2 = pPostStepPoint -> GetMaterial();

   if ( verboseLevel > 1 ) {
      G4cout << "Pre-step  point =  " << pPreStepPoint->GetPosition() << " Material = " << Material1->GetName() << G4endl;
      G4cout << "Post-step point =  " << pPostStepPoint->GetPosition() << " Material = " << Material2->GetName() << G4endl;
   }

   if (Material1 == Material2)// If both sides of the interface are of the same material, then physically speaking there is no interface
   {
      theStatus = SameMaterial;
      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
   }

   OldEnergy = aParticle->GetKineticEnergy();

   G4MaterialPropertiesTable* aMPT;
   G4int cond1 = 0;
   workFunction1 = 0;
   affinity1 = 0;
   fermiEnergy1 = 0;
   bandbending1 = 0;
   deltaphi1 = 0;
   G4double bandgap1 = 0;
   aMPT = Material1->GetMaterialPropertiesTable();
   // Get material properties. An overview:
   // Work function - energy between Fermi level and vacuum
   // Fermi energy - kinetic energy of the electron at the Fermi level (often zero for semiconductors and insulators)
   // Affinity - energy between bottom of the conduction band and vacuum. Used instead of work function for insulators.
   // Bandbending - bending of energy bands near the (vacuum) interface.
   // Deltaphi - additional potential of the whole sample (e.g. due to conductive contact with another material)

   // A note on the meaning of the latter two parameters: at an interface between two materials, bandbending is ignored while the deltaphis
   // lead to an alteration of the potential step; in other words, the potential step is assumed to be always fully localized.
   // At a vacuum interface however, the localized potential step is NOT altered by deltaphi; instead, it is assumed to be present
   // in the form of a patch field that affects the kinetic energy of the particle *before* it reaches the surface. Something similar
   // holds true for the bandbending (where the field due to deltaphi is always in the vacuum while the effect of bandbending is
   // noticeable on both sides of the surface). See also below in DielectricDielectric().

   if (aMPT) {
      if(aMPT->ConstPropertyExists("CONDUCTORTYPE")) cond1 = (int)(aMPT->GetConstProperty("CONDUCTORTYPE")); // cond==0: metal, cond==1: semi, cond==2 insulator
      if(aMPT->ConstPropertyExists("WORKFUNCTION")) workFunction1 = aMPT->GetConstProperty("WORKFUNCTION");
      if(aMPT->ConstPropertyExists("AFFINITY")) affinity1 = aMPT->GetConstProperty("AFFINITY");
      if(aMPT->ConstPropertyExists("FERMIENERGY")) fermiEnergy1 = aMPT->GetConstProperty("FERMIENERGY");
      if(aMPT->ConstPropertyExists("BANDBENDING")) bandbending1 = aMPT->GetConstProperty("BANDBENDING");
      if(aMPT->ConstPropertyExists("DELTAPHI")) deltaphi1 = aMPT->GetConstProperty("DELTAPHI");
      if(aMPT->ConstPropertyExists("BANDGAP")) bandgap1 = aMPT->GetConstProperty("BANDGAP");
      if(!workFunction1) workFunction1 = affinity1;
   }
   if(cond1!=0 && fermiEnergy1==0) {
      U1 = workFunction1 + bandgap1;
   } else if(cond1!=0) {
     U1 = workFunction1 + fermiEnergy1 + bandgap1/2.0;
   } else {
     U1 = workFunction1 + fermiEnergy1;
   } // AT: Redefined the potential step for insulators and semiconductors to take the band gap into account improved for the case when fermiEnergy is zero

   G4int cond2 = 0;
   workFunction2 = 0;
   affinity2 = 0;
   fermiEnergy2 = 0;
   bandbending2 = 0;
   deltaphi2 = 0;
   G4double bandgap2 = 0;
   aMPT = Material2->GetMaterialPropertiesTable();
   // Same for second material.
   if (aMPT) {
      if(aMPT->ConstPropertyExists("CONDUCTORTYPE")) cond2 = (int)(aMPT->GetConstProperty("CONDUCTORTYPE"));
      if(aMPT->ConstPropertyExists("WORKFUNCTION")) workFunction2 = aMPT->GetConstProperty("WORKFUNCTION");
      if(aMPT->ConstPropertyExists("AFFINITY")) affinity2 = aMPT->GetConstProperty("AFFINITY");
      if(aMPT->ConstPropertyExists("FERMIENERGY")) fermiEnergy2 = aMPT->GetConstProperty("FERMIENERGY");
      if(aMPT->ConstPropertyExists("BANDBENDING")) bandbending2 = aMPT->GetConstProperty("BANDBENDING");
      if(aMPT->ConstPropertyExists("DELTAPHI")) deltaphi2 = aMPT->GetConstProperty("DELTAPHI");
      if(aMPT->ConstPropertyExists("BANDGAP")) bandgap2 = aMPT->GetConstProperty("BANDGAP");
      if(!workFunction2) workFunction2 = affinity2;
   }
   if(cond2!=0 && fermiEnergy2==0) {
      U2 = workFunction2 + bandgap2;
   } else if(cond2!=0) {
     U2 = workFunction2 + fermiEnergy2 + bandgap2/2.0;
   } else {
     U2 = workFunction2 + fermiEnergy2;
   } // AT: Redefined the potential step for insulators and semiconductors to take the band gap into account improved for the case when fermiEnergy is zero

   Ugain = U2 + deltaphi2 - (U1 + deltaphi1);// Total kinetic energy gained when going from material 1 to material 2
   deltaU = Ugain;// Maximum height of the potential barrier to be overcome; criterion for transmission or reflection

   vacuuminterface = false;
   if ( U1==0 || U2==0 )// We're at a vacuum interface, either from inside or outside
   {
      vacuuminterface = true;
      deltaU = U2 - U1;// The height of the potential barrier at a vacuum interface does not include deltaphi,
      // since this is expressed as a patch field on the vacuum side rather than an localized potential jump
   }

   if ( verboseLevel > 3 ) {
      G4cout << " Electron at Boundary! " << G4endl;
      G4cout << " Old Momentum Direction: " << OldMomentum     << G4endl;
   }

   DielectricDielectric();// Perform the actual reflection/transmission test

   if (theStatus == TotalInternalReflection && deltaU<0.)
      // First criterion for killing the particle after reflection:
      // it may be killed if going from (e.g.) a material to vacuum, but it's always simply reflected if going from vacuum to a material
   {
      // A second criterion is about the particle's energy relative to the potential step:
      if ( exp(1.+OldEnergy/(2.*deltaU)) > G4UniformRand() ) {// Just a guess - 'absorption' is most likely at low energies but
         // reflection should be more probable at very high energies (e.g. grazing incidence at interfaces in STEM simulations)
         NewEnergy = 0.;
         aParticleChange.ProposeTrackStatus(fStopAndKill);
         aParticleChange.ProposeLocalEnergyDeposit(OldEnergy);
         theStatus = ParticleKilled;
      }
   }

   NewMomentum = NewMomentum.unit();

   if ( verboseLevel > 3) {
      G4cout << " New Momentum Direction: " << NewMomentum     << G4endl;
      if ( theStatus == Undefined ) G4cout << " *** Undefined *** " << G4endl;
      if ( theStatus == Refraction ) G4cout << " *** Refraction *** " << G4endl;
      if ( theStatus == TotalInternalReflection ) G4cout << " *** TotalInternalReflection *** " << G4endl;
      if ( theStatus == NotAtBoundary ) G4cout << " *** NotAtBoundary *** " << G4endl;
      if ( theStatus == SameMaterial ) G4cout << " *** SameMaterial *** " << G4endl;
      if ( theStatus == StepTooSmall ) G4cout << " *** StepTooSmall *** " << G4endl;
      if ( theStatus == ParticleKilled ) G4cout << " *** ParticleKilled *** " << G4endl;
   }

   aParticleChange.ProposeMomentumDirection(NewMomentum);
   aParticleChange.ProposeEnergy(NewEnergy);

   return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}


void CADPhysicseBoundary::DielectricDielectric()
// Perform the actual reflection/transmission test
// After invoking this method, NewEnergy and NewMomentum will be set to their proper final values
// and theStatus will reflect the kind of process performed:
// theStatus == TotalInternalReflection: reflection at the interface
// theStatus == Refraction: the electron gets transmitted (with adjusted energy)
{
   theFacetNormal = theGlobalNormal;

   G4double PdotN = OldMomentum * theFacetNormal;
   G4double normalEnergy = OldEnergy * PdotN * PdotN;// component of the kinetic energy normal to the surface
   if(vacuuminterface) normalEnergy += bandbending1 + bandbending2 + deltaphi2;
   // two cases at once:
   // material->vacuum: kinetic energy increases by bandbending1
   // vacuum->material: kinetic energy increases by deltaphi2 + bandbending2
   // (which could be either positive or negative)
   // normalEnergy is the component of the kin. energy normal to the interface AT the interface
   G4double newNormalEnergy = normalEnergy + deltaU - std::max(0.0,deltaphi1);
   // assuming the work functions are defined as positive numbers
   // if deltaphi1 is positive, it means the material is at a positive potential, and
   // the electron will need to overcome an extra energy to reach 'vacuum', otherwise it will be reflected
   // by the external patch field.
   // newNormalEnergy is the component of the kin. energy normal to the interface BEYOND the interface.
   // Note: in a more sophisticated model, all patch fields etc. should be modeled explicitly, instead of by
   // using the material parameters as implemented here.

   if ( verboseLevel > 3) {
      G4cout << "CADPhysicseBoundary: normalEnergy=" << normalEnergy/eV << " eV; " << G4endl;
      G4cout << "              and newNormalEnergy=" << newNormalEnergy/eV << " eV." << G4endl;
   }

   if (0. >= newNormalEnergy || 0. > normalEnergy )// Either energy is below zero, the electron cannot be transmitted
   {
      // Simulate total internal reflection
      theStatus = TotalInternalReflection;
      NewMomentum = OldMomentum - (2.*PdotN)*theFacetNormal;
      NewEnergy = OldEnergy;
   } else {
      // Do quantum mechanical calculation of transmission probability
      // Note: both a positive and a negative potential step deltaU may lead to reflection!
      G4double x = -1. * deltaU / normalEnergy;
      G4double y = sqrt(1-x);
      G4double T = 4.*y/((1.+y)*(1.+y));
      G4double randT = G4UniformRand();
      if(randT<T)
      {
         // Simulate transmission/refraction
         theStatus = Refraction;

         NewEnergy = OldEnergy + Ugain;
         G4ThreeVector perpMomentum = OldMomentum - PdotN * theFacetNormal;
         G4double alpha = NewEnergy/OldEnergy;
         G4double tmp = alpha - perpMomentum * perpMomentum;
         if (tmp<0) tmp=0;
         G4ThreeVector newnormMomentum = - sqrt(tmp) * theFacetNormal;
         NewMomentum = perpMomentum + newnormMomentum;
         NewMomentum = NewMomentum.unit();
      } else {
         // Simulate total internal reflection
         theStatus = TotalInternalReflection;
         NewMomentum = OldMomentum - (2.*PdotN)*theFacetNormal;
         NewEnergy = OldEnergy;
      }
   }
}

G4double CADPhysicseBoundary::GetMeanFreePath(const G4Track& ,
                                   G4double ,
                                   G4ForceCondition* condition)
{
   *condition = Forced;
   return DBL_MAX;
}
