// CADPhysicsSingleScattering.cc
// See CADPhysicsSingleScattering.hh for a general description of this class.

// 2006-07-07: introduced energy loss due to atom recoil, for electron energies > 50 eV
// 2006-07-10: changed folder for id-files to current working directory; always overwrite existing file
// 2006-07-17: in PostStepDoMultiStep: code for the momentum change written out (instead of using RotateUz) to save calculation time
// 2006-07-20: new implementation of cross section calculations for low energies. Now using acoustic phonon scattering cross sections,
//            derived from resistivity (for metals) and a model by Fitting and Schreiber (for semiconductors; see the header file for details).
//            As part of the change, the number of energy bins has been increased.
// 2006-07-24: changed storage of cumulative diff. cross section for angle for faster calculation in PostStepDoIt;
//            in PostStepDoMultiStep, for energybin<42, the calculation of xperel has been pulled outside the h loop
// 2006-08-23: changed selection of the model (metal/semiconductor). Instead of just reading the 'conductortype' property from the
//            MaterialPropertiesTable, the process now checks if it can find the parameters 'soundvelocity', 'defpotential' and 'lattice'
//            that it needs to apply the acoustic phonon scattering description in DInvMFPTableforSemi. If not, it will assume that the material
//            is a metal, and if this also doesn't work (and the material is not a gas, for which different processes may apply) the
//            BuildPhysicsTable method will give a hard exit.
// 2006-09-15: set fermienergy to zero for gaseous materials to prevent 'hanging' CADPhysicsDI
// 2007-07-24: (Temporarily) disabled for gases to test Cascade model
// 2007-07-30: Removed 2nd energy bin at 30 keV - instead created a more smooth transition between Mott and Browning data
// 2007-08-03: Added a switch for scattering in gases. The switch is currently set by the PhysicsList.
// 2007-11-21: Fixed bug and implemented more accurate version of momentum direction distribution in DoDiffusionStep
// 2007-11-26: Changed naming and implementation of 'id' to 'tmfp' to get more consistent terminology
// 2007-11-27: Added additional comments
// 2007-11-27: Cleared up nomenclature: replaced 'CS'/'cross section' with 'InvMFP'/'inverse mean free path' wherever the latter was meant.
// 2010-06-14: Replaced Browning cross sections with screened relativistic Rutherford cross section.
// 2012-07-27: Refined angular interpolation (for higher energies) and corrected Rutherford cross section (Reimer implementation).
// 2012-07-30: Updated angular interpolation/integration also for metals at 'Mott' energies
// 2013-01-21: Fixed error in element selection in DoSingleStep that was introduced in the previous version.

#include "CADPhysicsSingleScattering.hh"
#include "CADPhysicsSSMessenger.hh"
#include "CADPhysicsTransportation.hh"
#include "Randomize.hh"
#include "G4ProductionCutsTable.hh"
#include <sstream>
#include "G4MaterialPropertiesTable.hh"
#include "G4GPILSelection.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Material.hh"

#include "CADPhysicsVacuum.hh"

using namespace std;

CADPhysicsSingleScattering* CADPhysicsSingleScattering::fSingleScattering = 0;

CADPhysicsSingleScattering::CADPhysicsSingleScattering(const G4String& processName)
: G4VContinuousDiscreteProcess(processName),
currentCouple(0),
currentMaterial(0),
enlossconst (2.*electron_mass_c2/c_squared),
previousenergy(0.),
previousMFP(0.),
domultistep(false)
{
   SetVerboseLevel(1);
   messenger = new CADPhysicsSSMessenger(this);
   singlestep = false;
   multistep = false;
   noraytrace = false;
   safetycheck = false;
   fphononloss = 0.;
}

CADPhysicsSingleScattering* CADPhysicsSingleScattering::GetInstance()
{
   static CADPhysicsSingleScattering theSingleScattering;
   if (!fSingleScattering)
   {
      fSingleScattering = &theSingleScattering;
   }
   return fSingleScattering;
}

CADPhysicsSingleScattering::~CADPhysicsSingleScattering()
{}

void CADPhysicsSingleScattering::BuildPhysicsTable(
   const G4ParticleDefinition& part)
{
   if(verboseLevel) {
      G4cout << "CADPhysicsSingleScattering::BuildPhysicsTable() for "
         << GetProcessName()
         << " and particle " << part.GetParticleName()
         << G4endl;
   }

   // Access to materials
   const G4ProductionCutsTable* theCoupleTable=
      G4ProductionCutsTable::GetProductionCutsTable();

   // Create new instances of theInvMFPTable and theDInvMFPTable. For each material,
   // theInvMFPTable contains a pointer to a data vector that contains the total inverse mean free path
   // as a function of energy.
   // theDInvMFPTable is a four-dimensional data structure. For each element in each material,
   // and for each energy in the enrange list, it contains the cumulative differential inverse mean free path
   // as a function of scattering angle. It is used in PostStepDoIt to determine the scattering angle.
   // theSumDInvMFPTable holds similar information as theDInvMFPTable, except that the differential inverse MFPs
   // are summed over all elements in the material. This one is used if the individual element responsible
   // for scattering does not have to or cannot be known.
   vec_phononloss.clear();
   vec_effective_A.clear();

   size_t numOfCouples = theCoupleTable->GetTableSize();
   for(size_t i=0; i<numOfCouples; i++)// Loop over all materials
   {
      const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(i);
      const G4Material* material = couple->GetMaterial();
      G4String mname = material->GetName();
      std::ostringstream ost;
      ost << mname << ".mat.hdf5";
      std::string const & matfile = ost.str();
      G4cout << matfile << G4endl;
      if(verboseLevel) G4cout << "CADPhysicsSingleScattering: Read in Material " << mname << "\n";
      if (mname != "Galactic") {
        load_material(matfile, 1.e4, 128, 1024);
      } else {
        load_vacuum();
      }
   }

   if(1 < verboseLevel) {
      PrintInfoDefinition();
      G4cout << "CADPhysicsSingleScattering::BuildPhysicsTable() done for "
         << GetProcessName()
         << " and particle " << part.GetParticleName()
         << G4endl;
   }
}

void CADPhysicsSingleScattering::load_material(std::string const & filename, double high_energy, size_t N_K, size_t N_P)
{
  material mat(filename);

  // Get desired lower energy limit for use in the simulation
  const double low_energy = mat.get_barrier().value;
  // Get the energy ranges stored in the hdf5 file.
  auto el_range = mat.get_elastic_energy_range();

  // Print diagnostic info
  G4cout << "Material: " << mat.get_name() << G4endl
    << "  Phonon loss     = " << (mat.get_phonon_loss() / units::eV).value << " eV" << G4endl
    << "  Density         = " << (mat.get_density() * (units::cm * units::cm * units::cm)).value << " cm^-3" << G4endl
    << "  Elastic range   = [" << el_range.first << ", " << el_range.second << "] eV" << G4endl;

    vec_phononloss.push_back((mat.get_phonon_loss() / units::eV).value*eV);
    vec_effective_A.push_back((mat.get_effective_A() / units::g).value*g);

  // Warn if tabulated data has too narrow range
  if (el_range.first > low_energy)
    std::cerr << "WARNING: extrapolating elastic cross sections between " << el_range.first << " and " << low_energy << " eV" << std::endl;
  if (el_range.second < high_energy)
    std::cerr << "WARNING: extrapolating elastic cross sections between " << el_range.second << " and " << high_energy << " eV" << std::endl;

  /*
    * Now build the tables for use in the simulation. If the hdf5 tables
    * contain insufficient data, we store only this insufficient range
    * for use in the simulation.
    * The simulation tables extrapolate anyway, and they extrapolate in
    * the correct space (log space for imfp tables, while the material
    * class stores imfp tables in linear space).
    */
  const double elastic_low = std::max(low_energy, el_range.first);
  const double elastic_high = std::min(high_energy, el_range.second);

  elastic_imfp_vector.push_back(mat.get_elastic_imfp(elastic_low, elastic_high, N_K));
  elastic_icdf_vector.push_back(mat.get_elastic_angle_icdf(elastic_low, elastic_high, N_K, N_P));
}

void CADPhysicsSingleScattering::load_vacuum()
{
    CADPhysicsVacuum vac;
    vec_phononloss.push_back(0.);
    vec_effective_A.push_back(0.);

  /*
    * Now build the tables for use in the simulation. If the hdf5 tables
    * contain insufficient data, we store only this insufficient range
    * for use in the simulation.
    * The simulation tables extrapolate anyway, and they extrapolate in
    * the correct space (log space for imfp tables, while the material
    * class stores imfp tables in linear space).
    */

  elastic_imfp_vector.push_back(vac.get_vacuum_imfp()); // this is just to write something, so that the material index is still correct
  elastic_icdf_vector.push_back(vac.get_vacuum_icdf());
}

G4VParticleChange* CADPhysicsSingleScattering::AlongStepDoIt(
   const G4Track& track,const G4Step& /*step*/)//track,step
{
   // In the case of multistep, this is the place to propose energy loss, new position, and time spent via the ParticleChange
   aParticleChange.Initialize(track);
   if (multistep) {
      aParticleChange.ProposeLocalEnergyDeposit(energyloss);
      aParticleChange.ProposePosition(proposeposition);

      // Simple linearized approximation of the time spent in this step. Time calculations are normally done in the Transportation process,
      // but since we're bypassing it, this needs to be done here.
      G4double deltaTime = track.GetStepLength() / track.GetVelocity();
      aParticleChange.ProposeGlobalTime(track.GetGlobalTime() + deltaTime);
      aParticleChange.ProposeProperTime(track.GetProperTime() + deltaTime);
   }
   return &aParticleChange;
}

G4VParticleChange* CADPhysicsSingleScattering::PostStepDoIt(
   const G4Track& track,
   const G4Step& step)
{
   aParticleChange.Initialize(track);
   if (donothing) {
      return pParticleChange;// Do nothing, returns pointer to aParticleChange
   }

   if (singlestep) {
      return DoSingleStep(track,step);// Do a 'conventional' single scatter event
   }

   // In the case of multistep, this is the place to propose new direction and final kinetic energy via the ParticleChange
   aParticleChange.ProposeMomentumDirection(elDirectionnew.unit());
   aParticleChange.ProposeEnergy(finalT);
   return G4VContinuousDiscreteProcess::PostStepDoIt(track, step);
}

G4VParticleChange* CADPhysicsSingleScattering::DoSingleStep(
   const G4Track& track,
   const G4Step& step)
{
   // Some initialization
   aParticleChange.Initialize(track);
   finalT = track.GetKineticEnergy();
   energyloss = 0.;
   size_t i = currentMaterialIndex;
   size_t nelm = currentMaterial->GetNumberOfElements() ;
   size_t ii=0;
   G4double cosTheta = std::cos(elastic_icdf_vector[i].get(finalT/eV, G4UniformRand()));

   if(finalT<200*eV)// Below 200 eV, we use the phonon approximation to calculate the energy loss.
   {
      // Energy loss in the scatter event is taken as the fixed average phonon loss value for this material.
      energyloss = fphononloss;
      finalT -= energyloss;
      if (finalT>0.) aParticleChange.ProposeEnergy(finalT);
      aParticleChange.ProposeLocalEnergyDeposit(energyloss);
   } else {// Above 200 eV, the procedure is different since we also need to select the right element (in a compound)
      // Calculate the "elastic" energy loss due to atom recoil
      // N.B.: this energy loss is extremely small compared to the kinetic energy - the place where you are most likely to notice it
      // is in the elastic backscatter peak, for which no other loss processes are present.
      G4double effA  = vec_effective_A[i]; // This is the effective atom mass in gram
      energyloss = enlossconst * (1.-cosTheta) * finalT / effA;
      finalT -= energyloss;
      aParticleChange.ProposeEnergy(finalT);
      aParticleChange.ProposeLocalEnergyDeposit(energyloss);
   }

   // Randomly select the phi angle from a uniform distribution, and determine the new momentum vector
   G4double Phi     = twopi * G4UniformRand();
   G4double sinTheta = sqrt ( 1. - cosTheta*cosTheta );
   G4double dirx = sinTheta*cos(Phi), diry = sinTheta*sin(Phi), dirz = cosTheta;
   const G4DynamicParticle* dp = track.GetDynamicParticle();
   G4ThreeVector elDirection = dp->GetMomentumDirection();
   elDirectionnew = G4ThreeVector( dirx,diry,dirz );
   elDirectionnew.rotateUz(elDirection);

   // update G4VParticleChange for the scattered electron direction
   aParticleChange.ProposeMomentumDirection(elDirectionnew);

   return G4VContinuousDiscreteProcess::PostStepDoIt(track, step);
}

G4double CADPhysicsSingleScattering::DoMultiStep(
   const G4Track& track,
   G4double othersteplength,
   G4double safety)
   // Multistep version of the PostStepDoIt method. It may be combined with the
   // PostStepDoIt methods of other processes, but always needs to be performed first, to make sure
   // that the other processes happen at the correct position and momentum. Therefore, because it needs to propose a *position*
   // (unlike other PostStepDoIt methods), and because its total step length is not a priori known,
   // DoMultiStep is invoked from the AlongStepGetPhysicalInteractionLength method!
{
   G4bool finalstep = false;
   G4bool testforsafety = (othersteplength>safety);// othersteplength is the limiting step from the other processes. If it is less than
   // the current safety value, then it is necessarily limiting and there is no need to check against the safety after each individual step.
   // This can save us some time, and hence we define a special boolean for it.
   G4double nextstep = lastgpilvalue;// This is the step length taken *before* the first scatter event
   G4double steptaken = 0.;
   G4double safety2 = safety * safety;
   G4ThreeVector initialposition = track.GetPosition();
   G4ThreeVector momdir = track.GetMomentumDirection();

   G4double kinen = track.GetKineticEnergy();
   finalT = kinen;

   size_t i = currentMaterialIndex;
   size_t nelm = currentMaterial->GetNumberOfElements();

   G4double costheta1,sintheta1,cosphi1,sinphi1;
   G4double tempx,tempy,tempz;
   G4ThreeVector relativeglobalposition = G4ThreeVector(0.,0.,0.);

   if (kinen<200*eV) {// Below 200 eV, we use the phonon approximation to calculate the energy loss.
      // Hence, there is no need to determine on which individual element the electron was scattered
      // and we can use the summed differential inverse mean free paths per material instead of the individual
      // differential inverse mean free paths per element - this makes the simulation a little bit faster.
      energyloss = 0.;// Initialize the total energy loss to zero

      do// Inside this loop, all the processing for one single step/scatter event takes place
      {
         relativeglobalposition += momdir * nextstep;// Take the step with the previously defined step length

         // Test for safety
         if (testforsafety && relativeglobalposition.mag2()>safety2) {// Oops! We have travelled outside the sphere defined by the 'safety' parameter
            relativeglobalposition -= momdir * nextstep;// Reverse the last step
            safetycheck = true;// Signal that the current instance of DoMultiStep was limited by the safety
            if (!finalstep) lastgpilvalue = nextstep;// Store the last steplength for the next step
            break;// Exit the loop
         }

         // We passed the safety check, so the latest step length can be added to the total step taken
         steptaken += nextstep;
         if (finalstep) break;// Exit the loop if this was determined to be the final step

         // Now determine cosTheta for the next scatter event
         G4double cosTheta = -1.;// Initialize
         cosTheta = std::cos(elastic_icdf_vector[i].get(kinen/eV, G4UniformRand()));

         // Randomly select the phi angle from a uniform distribution
         G4double Phi = twopi * G4UniformRand();
         G4double cosPhi = cos(Phi);
         G4double sinPhi = sin(Phi);
         G4double sinTheta = sqrt ( 1. - cosTheta*cosTheta );

         // Determine the new momentum vector
         // Note: the vector2.rotateUZ(vector1) method could be used for this, but writing this out was found to be somewhat faster
         costheta1 = momdir.z();
         sintheta1 = std::sqrt(momdir.x()*momdir.x()+momdir.y()*momdir.y());
         cosphi1 = 1.;
         sinphi1 = 0.;
         if (sintheta1>0.) {
            cosphi1 = momdir.x()/sintheta1;
            sinphi1 = momdir.y()/sintheta1;
         }
         tempx = cosPhi * cosphi1 * costheta1 - sinPhi * sinphi1;
         tempy = cosPhi * sinphi1 * costheta1 + sinPhi * cosphi1;
         tempz = -cosPhi * sintheta1;
         tempx *= sinTheta;
         tempy *= sinTheta;
         tempz *= sinTheta;
         momdir = cosTheta * momdir + G4ThreeVector(tempx,tempy,tempz);

         // Finally add the (average) energy loss due to phonon scattering to the total energy loss
         G4double locenloss = fphononloss;
         energyloss += locenloss;
         if (energyloss>0.5*eV) break;// A limit to the total energy loss is implemented here, to make sure that all MFPs are still valid

         // Finally determine the next step size
         nextstep = -std::log(G4UniformRand())*currentInteractionLength;
         if (steptaken + nextstep > othersteplength)// If the total step length exceeds the step limit imposed by other processes,
            // then just move the particle over the remaining step length and declare this to be the final step
         {
            lastgpilvalue = nextstep;
            nextstep = othersteplength - steptaken;
            finalstep = true;
         }
      } while (1);
   } else {// Above 200 eV, we use the Mott inverse MFPs and the atomic scattering parameters for the energy loss.
      // Hence we need to select the scattering element for each individual scattering event.
      G4double sumforenloss = 0.;// Initialization for energy loss

      do{// Inside this loop, all the processing for one single step/scatter event takes place
         relativeglobalposition += momdir * nextstep;// Take the step with the previously defined step length

         // Test for safety
         if (testforsafety && relativeglobalposition.mag2()>safety2) {// Oops! We have travelled outside the sphere defined by the 'safety' parameter
            relativeglobalposition -= momdir * nextstep;// Reverse the last step
            safetycheck = true;// Signal that the current instance of DoMultiStep was limited by the safety
            if (!finalstep) lastgpilvalue = nextstep;// Store the last steplength for the next step
            break;// Exit the loop
         }

         // We passed the safety check, so the latest step length can be added to the total step taken
         steptaken += nextstep;
         if (finalstep) break;// Exit the loop if this was determined to be the final step

         // Select the element and select the scattering angle (theta) for the next scattering event
         G4double cosTheta = -1.;// Initialize
         cosTheta = std::cos(elastic_icdf_vector[i].get(kinen/eV, G4UniformRand()));

         G4double effA  = vec_effective_A[i];
         sumforenloss += (1 - cosTheta)/effA;// The energy loss is proportional to this value. Only this factor is summed
         // for every step; the actual energy loss is calculated at the end.

         // Randomly select the phi angle from a uniform distribution
         G4double Phi     = twopi * G4UniformRand();
         G4double cosPhi = cos(Phi);
         G4double sinPhi = sin(Phi);
         G4double sinTheta = sqrt ( 1. - cosTheta*cosTheta );

         // Determine the new momentum vector
         costheta1 = momdir.z();
         sintheta1 = std::sqrt(momdir.x()*momdir.x()+momdir.y()*momdir.y());
         cosphi1 = 1.;
         sinphi1 = 0.;
         if (sintheta1>0.) {
            cosphi1 = momdir.x()/sintheta1;
            sinphi1 = momdir.y()/sintheta1;
         }
         tempx = cosPhi * cosphi1 * costheta1 - sinPhi * sinphi1;
         tempy = cosPhi * sinphi1 * costheta1 + sinPhi * cosphi1;
         tempz = -cosPhi * sintheta1;
         tempx *= sinTheta;
         tempy *= sinTheta;
         tempz *= sinTheta;
         momdir = cosTheta * momdir + G4ThreeVector(tempx,tempy,tempz);

         // Finally determine the next step size
         nextstep = -std::log(G4UniformRand())*currentInteractionLength;
         if (steptaken + nextstep > othersteplength)// If the total step length exceeds the step limit imposed by other processes,
            // then just move the particle over the remaining step length and declare this to be the final step
         {
            lastgpilvalue = nextstep;
            nextstep = othersteplength - steptaken;
            finalstep = true;
         }
      } while (1);
      // Calculate the total energy loss
      energyloss = enlossconst * sumforenloss * finalT;

   }//end if(kinen)
   // Update position, energy and direction for later registration via ParticleChange.
   proposeposition = initialposition + relativeglobalposition;
   if(energyloss<finalT) finalT -= energyloss;// Extra check: energy loss cannot be larger than total available energy
   elDirectionnew = momdir;
   return steptaken;
}

G4double CADPhysicsSingleScattering::GetMeanFreePath(
   const G4Track& track,
   G4double,
   G4ForceCondition*)
{
   //G4double preStepMFP;
   G4double preStepKinEnergy = track.GetKineticEnergy();
   const G4MaterialCutsCouple* couple = track.GetMaterialCutsCouple();
   currentCouple   = couple;
   currentMaterial = couple->GetMaterial();
   currentMaterialIndex = couple->GetIndex();
   G4String mname = currentMaterial->GetName();
   G4double mfp = DBL_MAX;
   previousMFP = mfp;
   fphononloss = vec_phononloss[currentMaterialIndex];
   previousenergy = preStepKinEnergy;

   if (mname == "Galactic") { // prevent the read out of elastic_imfp_vector for Galactic
     return mfp;
  }

  G4double cross = elastic_imfp_vector[currentMaterialIndex].get(preStepKinEnergy/eV); // nm^-1
  if (cross > 0.) {
     mfp = 1./cross; // nm
     mfp = mfp * 1.e-6; // mm
     previousMFP = mfp;
  }
  return mfp;
}

G4double CADPhysicsSingleScattering::AlongStepGetPhysicalInteractionLength(
   const G4Track& track,
   G4double previousStepSize,
   G4double PhysicalStep,
   G4double& currentSafety,
   G4GPILSelection* selection)
{
   *selection = CandidateForSelection;

   G4double value = lastgpilvalue;// This is the current physical interaction length
   singlestep = false;// Initialization
   multistep = false;
   donothing = false;

   G4GPILSelection* tselection = new G4GPILSelection();// dummy
   G4double transvalue = CADPhysicsTransportation::GetInstance()->RealASGPIL(track, previousStepSize, PhysicalStep, currentSafety, tselection);
   // This does several things:
   // - it updates the currentSafety value
   // - it compares the current minimum step (including the one from this process) to the current safety.
   // If the minimum step is larger, it sets the noraytrace switch in this process to true and proceeds to do the raytracing.
   // If not, it returns a DBL_MAX value and leaves the raytracing to the current process.

   delete tselection;
   PhysicalStep = min(transvalue,PhysicalStep);

   if (noraytrace || (!domultistep)) {// Do no 'ray tracing'
   //AT: add switch. if domultistep is false, a multistep is never performed
      noraytrace = false;
      if (value > PhysicalStep) {// Another process wins, so do nothing
         donothing = true;
         return DBL_MAX;
      } else {// This process has the shortest physical interaction length, so do a conventional single step later
         singlestep = true;
         return value;
      }
   } else {// Do the 'raytracing' - N.B.: raytracing here means just taking linear steps. Use SSforCharging instead
      // if the local E-field needs to be taken into account
      G4double mysafety = min(currentSafety,PhysicalStep);
      multistep = true;// The requirements for a diffusion step were not met, hence do a 'normal' multistep
      return DoMultiStep(track,PhysicalStep,currentSafety);
   }
   return value;// As a backup, we should not actually get here
}

G4double CADPhysicsSingleScattering::PostStepGetPhysicalInteractionLength(
   const G4Track& track,
   G4double previousStepSize,
   G4ForceCondition* condition
   )
{
   if (safetycheck && currentInteractionLength > 0.) {// The previous step was a multistep, and it ended
      // on a safety check. This means that the next step size was already determined but not yet executed.
      // Here we update theNumberOfInteractionLengthLeft to account for this fact.
      safetycheck = false;
      theNumberOfInteractionLengthLeft = lastgpilvalue / currentInteractionLength;
   } else {
      safetycheck = false;
      if ( (previousStepSize <=0.0) || (theNumberOfInteractionLengthLeft<=0.0)) {
         // We're at the beginning of tracking, or just after a normal single/multistep of this process
         // In the latter case, G4VContinuousDiscreteProcess::PostStepDoIt has called the
         // G4VProcess::ClearNumberOfInteractionLengthLeft() method which by default sets
         // theNumberOfInteractionLengthLeft to -1.
         ResetNumberOfInteractionLengthLeft();// Draw a new random-distributed theNumberOfInteractionLengthLeft.
      } else {
         // Only after 'donothing':
         // subtract previousStepSize/currentInteractionLength from NumberOfInteractionLengthLeft
         SubtractNumberOfInteractionLengthLeft(previousStepSize);
         if(theNumberOfInteractionLengthLeft<0.) {
            theNumberOfInteractionLengthLeft=perMillion;
         }
      }
   }

   // Get the mean free path
   currentInteractionLength = GetMeanFreePath(track, previousStepSize);//NO G4ForceCondition* passed on

#ifdef G4VERBOSE
   if ((currentInteractionLength <=0.0) || (verboseLevel>4)){
      G4cout << "CADPhysicsSingleScattering::PostStepGetPhysicalInteractionLength ";
      G4cout << "[ " << GetProcessName() << "]" <<G4endl;
      track.GetDynamicParticle()->DumpInfo();
      G4cout << " in Material  " << track.GetMaterial()->GetName() <<G4endl;
      G4cout << "MeanFreePath = " << currentInteractionLength/cm << "[cm]" <<G4endl;
   }
#endif

   if (currentInteractionLength <DBL_MAX) {
      lastgpilvalue = theNumberOfInteractionLengthLeft * currentInteractionLength;
   } else {
      lastgpilvalue = DBL_MAX;
   }
   *condition = Forced;// Always call PostStepDoIt
   return DBL_MAX;
}


void CADPhysicsSingleScattering::PrintInfoDefinition()
{
   G4cout << "CADPhysicsSingleScattering: elastic scattering\n"
      << " based on interpolated tabulated Mott cross sections\n"
      << " from Czyzewski et al., J.Appl.Phys. 68 (7), 3066 (1990),\n"
      << " screened relativistic Rutherford cross sections (above 30keV)\n"
      << " Williams & Carter, Transmission Electron Microscopy, Part-I basics, (1996) 39-40\n"
      << " and acoustic phonon scattering based on the description in\n"
      << " H.-J.Fitting et al., J.Electron Spectrosc.Rel.Phenom. 119 (2001) 35-47, and\n"
      << " E.Schreiber & H.-J.Fitting, J.Electron Spectrosc.Rel.Phenom. 124 (2002) 25-37." << G4endl;
}
