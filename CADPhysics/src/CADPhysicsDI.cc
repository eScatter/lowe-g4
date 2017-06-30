// CADPhysicsDI.cc
// See CADPhysicsDI.hh for a general description of the CADPhysicsDI class.

// History:
// 2006-07-10: changed folder for id-files to current working directory
//             and removed several redundant warning messages in PostStepDoIt
// 2006-07-27: fixed a bug in ReadEpsilonFile. If the df file does not exist, we still need to put
//             entries into outershells and outershelltable; otherwise we get the material couple
//             indices wrong.
// 2006-08-08: improved implementation of the 'id' values exported to file
// 2006-09-07: implementation of Auger electron generation through G4 atomic decay processes
// 2006-12-08: now adding fermi energy also to Auger electron kinetic energies
// 2007-01-05: ReadEpsilonFile only exits if df file does not exist AND material is not a gas
// 2007-04-20: Extra check in GetMeanFreePath: forced kill only if energy<barrier and particle NOT
//             at a geom.boundary
// 2007-05-02: Introduced some relativistic corrections to the total CS and angular distributions.
// 2007-05-02: Implementation of a boolean 'rangecut' for killing below-range electrons and a
//             counter for e-h pairs
// 2007-07-02: Energy check modified by local potential
// 2007-07-18: Changed calculation of Lvalue in TabulateL to reduce numerical errors
// 2007-07-18: Fixed bug in calculation of outershelldifflambdavalue in CalculateCDCS
// 2007-07-18: Fixed bug: added extra energy check for 'vary omegaprime' in calculation of 'id'
//             values in CalculateCDCS; resolves buffer overrun runtime error with Microsoft
//             Visual C++ 8 compiler
// 2007-07-20: Changed default location of low-energy data to make it version independent
// 2007-08-01: Added user-controllable 'energycut' variable that sets minimum total electron energy
//             wrt vacuum level
// 2007-08-03: Moved definition of env.var. G4LEDATA to PhysicsList, since it is required by the
//             constructor of the TransitionManager
// 2007-08-15: Set cutForPhotons to 250 eV to match validity of Penelope gamma processes
// 2007-08-17: Minus sign error (introduced on 2007-07-18) in Lvalue fixed
// 2007-11-26: Changed naming and implementation of 'id' to 'tmfp' to get more consistent
//             terminology
// 2012-03-08: Count SE's that are not followed as deposited energy (around line 1190) (EB)
// 2012-07-25: Added functionality so that user can have 'private' df files along with the global
//             ones, which are supposedly 'validated'
// 2013-06-06: Reformatted to width of 100 characters, some style changes
// 2017-02-03: Modified primary el. momentum change for non-outer-shell excitations (i.e. plasmon)

#include "CADPhysicsDI.hh"
#include "CADPhysicsDIMessenger.hh"
#include "CADPhysicsOutput.hh"
#include "Randomize.hh"
#include "G4ProductionCutsTable.hh"
#include <sstream>
#include <algorithm>
#include "G4MaterialPropertiesTable.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Gamma.hh" //for X-ray fluorescence
#include "G4Material.hh"
#include "G4CrossSectionHandler.hh"

#include "CADPhysicsVacuum.hh"

#ifndef CADPHYSICS_RELATIVISTIC
#define CADPHYSICS_RELATIVISTIC 1
#endif

using namespace std;

CADPhysicsDI::CADPhysicsDI(const G4String& processName)
: G4VContinuousDiscreteProcess(processName),
outershelltable(0),
//fconductortype("metal"),
generateSecondaries(true),
generateXrays(false),
rangecut(true),
pairsgenerated(0),
crossSectionHandler(0)
{
   SetVerboseLevel(1);

   // Create the messenger
   messenger = new CADPhysicsDIMessenger(this);

   // Tell the deexcitation manager to also generate Auger electrons.
   deexcitation.ActivateAugerElectronProduction(true);

   killthisone = false;
   DI_output = false;
}

CADPhysicsDI::~CADPhysicsDI()
{
   // Free up memory in the data structures (perhaps redundant, since this method is only called
   // when the program exits anyway)
   output.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..........oooOO0OOooo....

void CADPhysicsDI::BuildPhysicsTable(const G4ParticleDefinition& /*aParticleType*/)
{
   // Does the initialization of the process: fill the required data tables.

   if(verboseLevel) G4cout << "CADPhysicsDI: Start BuildPhysicsTable" << G4endl;

   // First clean up any existing data structures. This is necessary because this method may be
   // called more than once, in case the geometry has been modified between simulations.
   //if( outershelltable!=0 ) delete outershelltable;
   //outershelltable = new CADPhysicsDataTable;
   outershells.clear();
   vec_fermieff.clear();
   vec_bandgap.clear();
   vec_conductortype.clear();
   vec_barrier.clear();
   inelastic_imfp_vector.clear();
   inelastic_icdf_vector.clear();
   ionization_icdf_vector.clear();
   outershelltable.clear();
   outershells.clear();
   electron_range_vector.clear();
   rangecut_vector.clear();


   // Create material-independent tabulated data
   BuildEnrange(); //necessary for Psecenergies and Psecvalues

   // Initialization of  the Auger/X-ray fluorescence part
   // Create and fill G4CrossSectionHandler
   if ( crossSectionHandler ) delete crossSectionHandler;
   crossSectionHandler = new G4CrossSectionHandler();
   G4String shellCrossSectionFile = "phot/pe-ss-cs-";
   crossSectionHandler->LoadShellData(shellCrossSectionFile);
   crossSectionHandler->BuildMeanFreePathForMaterials();
   // End of the Auger/X-ray part

   // Create material-dependent data
   const G4ProductionCutsTable* theCoupleTable=
      G4ProductionCutsTable::GetProductionCutsTable();
   size_t numOfCouples = theCoupleTable->GetTableSize();
   // Loop for each material in the table
   for (size_t i=0;i<numOfCouples;i++) {
      const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(i);
      const G4Material* material = couple->GetMaterial();
      G4String mname = material->GetName();
      std::ostringstream ost;
      ost << mname << ".mat.hdf5";
      std::string const & matfile = ost.str();
      G4cout << matfile << G4endl;
      if(verboseLevel) G4cout << "CADPhysicsDI: Read in Material " << mname << "\n";
      if (mname != "Galactic") {
        load_material(matfile, 1.e4, 128, 1024);
      } else {
        load_vacuum();
      }
   }

   if(verboseLevel>1) PrintInfoDefinition();
}

void CADPhysicsDI::load_material(std::string const & filename, double high_energy, size_t N_K, size_t N_P)
{
  material mat(filename);

  // Get desired lower energy limit for use in the simulation
  const double low_energy = mat.get_barrier().value;
  // Get the energy ranges stored in the hdf5 file.
  auto in_range = mat.get_inelastic_energy_range();
  auto io_range = mat.get_ionization_energy_range();
  auto ran_range = mat.get_electron_range_energy_range();

  // Print diagnostic info
  std::cout << "Material: " << mat.get_name() << std::endl
    << "  Fermi           = " << (mat.get_fermi() / units::eV).value << " eV" << std::endl
    << "  Barrier         = " << (mat.get_barrier() / units::eV).value << " eV" << std::endl
    << "  Band gap        = " << (mat.get_band_gap() / units::eV).value << " eV" << std::endl
    << "  Phonon loss     = " << (mat.get_phonon_loss() / units::eV).value << " eV" << std::endl
    << "  Density         = " << (mat.get_density() * (units::cm * units::cm * units::cm)).value << " cm^-3" << std::endl
    << "  Inelastic range = [" << in_range.first << ", " << in_range.second << "] eV" << std::endl;

    vec_fermieff.push_back(G4double((mat.get_fermi() / units::eV).value)*eV);
    vec_bandgap.push_back((mat.get_band_gap() / units::eV).value*eV);
    vec_conductortype.push_back(mat.get_conductor_type());
    vec_barrier.push_back((mat.get_barrier() / units::eV).value*eV);

  // Warn if tabulated data has too narrow range
  if (in_range.first > low_energy)
    std::cerr << "WARNING: extrapolating inelastic cross sections between " << in_range.first << " and " << low_energy << " eV" << std::endl;
  if (in_range.second < high_energy)
    std::cerr << "WARNING: extrapolating inelastic cross sections between " << in_range.second << " and " << high_energy << " eV" << std::endl;
  if (ran_range.second < 1.*keV)
    std::cerr << "WARNING: the electron range is only known up to " << ran_range.second << " eV, but it should be known up to 1 keV." << std::endl;
  // No need to warn for the ionization table. Below the lowest tabulated value we should use outer shell energies, and otherwise the band gap.

  /*
    * Now build the tables for use in the simulation. If the hdf5 tables
    * contain insufficient data, we store only this insufficient range
    * for use in the simulation.
    * The simulation tables extrapolate anyway, and they extrapolate in
    * the correct space (log space for imfp tables, while the material
    * class stores imfp tables in linear space).
    */
  const double inelastic_low = std::max(low_energy, in_range.first);
  const double inelastic_high = std::min(high_energy, in_range.second);
  const double ionization_low = std::max(low_energy, io_range.first);
  const double ionization_high = std::min(high_energy, io_range.second);
  const double range_low = std::max(low_energy, ran_range.first);
  const double range_high = ran_range.second;

  inelastic_imfp_vector.push_back(mat.get_inelastic_imfp(inelastic_low, inelastic_high, N_K));
  inelastic_icdf_vector.push_back(mat.get_inelastic_w0_icdf(inelastic_low, inelastic_high, N_K, N_P));
  ionization_icdf_vector.push_back(mat.get_ionization_icdf(ionization_low, ionization_high, N_K, N_P));
  std::vector<float> tmp_outershells = mat.get_outer_shells();
  std::sort(tmp_outershells.begin(),tmp_outershells.end());
  outershelltable.push_back(tmp_outershells);
  if (!outershelltable.back().empty()) {
    outershells.push_back(true);
    for (float& f: outershelltable.back()) {
      G4cout << "outershells: " << f << G4endl;
      f *= eV;
    }
  } else {
    outershells.push_back(false);
  }
  electron_range_vector.push_back(mat.get_electron_range(range_low, range_high, N_K));
  rangecut_vector.push_back(rangecut);
}

void CADPhysicsDI::load_vacuum()
{
    CADPhysicsVacuum vac;
    vec_fermieff.push_back(0.);
    vec_bandgap.push_back(0.);
    vec_conductortype.push_back(material::CND_METAL);
    vec_barrier.push_back(0.);

  /*
    * Now build the tables for use in the simulation. If the hdf5 tables
    * contain insufficient data, we store only this insufficient range
    * for use in the simulation.
    * The simulation tables extrapolate anyway, and they extrapolate in
    * the correct space (log space for imfp tables, while the material
    * class stores imfp tables in linear space).
    */

  inelastic_imfp_vector.push_back(vac.get_vacuum_imfp()); // this is just to write something, so that the material index is still correct
  inelastic_icdf_vector.push_back(vac.get_vacuum_icdf());
  ionization_icdf_vector.push_back(vac.get_vacuum_ionization());
  std::vector<float> tmp;
  outershelltable.push_back(tmp);
  outershells.push_back(false);
  electron_range_vector.push_back(vac.get_vacuum_range());
  rangecut_vector.push_back(false); // never perform the rangecut for vacuum
}


void CADPhysicsDI::BuildEnrange()
{

   // Fill the vectors used for Fermi sea electron excitation; see also PostStepDoIt below
   G4double tempPsecenergies[21] = { 0.001,0.001584893,0.002511886,0.003981072,0.006309573,0.01,
      0.015848932,0.025118864,0.039810717,0.063095734,0.1,
      0.158489319,0.251188643,0.398107171,0.630957344,1.,
      1.584893192,2.511886432,3.981071706,6.309573445,10. };
   for(size_t j=0;j<21;j++) Psecenergies.push_back(tempPsecenergies[j]);
   G4double tempPsecvalues[21] = { 2.10882E-05,4.20838E-05,8.39916E-05,0.000167659,0.000334757,
      0.000668663, 0.001336482,0.002673959,0.005358324,0.010763755,0.021703431, 0.04401045,
      0.089993836,0.186221015,0.39157919,0.840316775, 1.846675208,4.162555928,9.616927383,
      22.70965162,54.59532675 };
   for(size_t j=0;j<21;j++) Psecvalues.push_back(tempPsecvalues[j]);
}

G4double CADPhysicsDI::GetContinuousStepLimit(const G4Track& /* track */,
   G4double, G4double /* currentMinimumStep */, G4double&)
{
   // Not implemented for this process, hence return the 'maximum' value
   return DBL_MAX;
}

G4VParticleChange* CADPhysicsDI::AlongStepDoIt(const G4Track& track,const G4Step& /*step*/)
{
   // Not implemented for this process
   aParticleChange.Initialize(track);
   return &aParticleChange;
}

G4VParticleChange* CADPhysicsDI::Phononloss( const G4Track& track,
   const G4Step& step,
   G4double kineticEnergy,
   G4double omegaprime)
{
   // Handle the discrete event in the case of sub-bandgap energy loss in semiconductors and
   // insulators; energy loss due to LO phonon excitation is assumed.
   G4double finalKinEnergy = kineticEnergy-omegaprime;
   G4double theEnergyDeposit = omegaprime;

   // Perform checks for kinetic energy vs. energy limit and electron range vs. safety
   G4StepPoint* steppoint = step.GetPostStepPoint();
   G4double pstepsafety = steppoint->GetSafety();
   G4double electronrange = electron_range_vector[findex].get(finalKinEnergy/eV) * nm;
   if ( finalKinEnergy < max(ffermienergy, fbarrier) ||
         ((rangecut_vector[findex] && finalKinEnergy < 1.*keV) && pstepsafety > 5.*electronrange) ) {
      // Too low remaining energy, kill the particle.
      //theEnergyDeposit += finalKinEnergy - ffermienergy;
      if(DI_output) {
        if(!output.is_open()) {
          output.open ("DIout.dat");
          output << std::setprecision(10);
          output << "EventID\tTrackID\tEtotal (eV)\tEkin (eV)"
                 << "\tx (nm)\ty (nm)\tz (nm)\txOrigin (nm)\tyOrigin (nm)\tzOrigin (nm)\tx-direction (nm)\ty-direction (nm)\tz-direction (nm)"
                 << "\tMax Depth (nm)\tMax Radius (nm)" << G4endl;
          }
          output << CADPhysicsOutput::CADPhysicsOutput(track, step) << G4endl;
      }
      theEnergyDeposit += finalKinEnergy;
      // EB: commented out this line as small negative deposits can happen
      //if(theEnergyDeposit<0.) theEnergyDeposit = 0.;
      finalKinEnergy    = 0.0;
      aParticleChange.ProposeTrackStatus(fStopAndKill);
   }

   aParticleChange.ProposeLocalEnergyDeposit(theEnergyDeposit);
   aParticleChange.ProposeEnergy(finalKinEnergy);
   return G4VContinuousDiscreteProcess::PostStepDoIt(track, step);
}

G4VParticleChange* CADPhysicsDI::PostStepDoIt( const G4Track& track, const G4Step& step)
{
   // Initialization
   aParticleChange.Initialize(track);
   G4double kineticEnergy = track.GetKineticEnergy();

   const G4MaterialCutsCouple* couple = track.GetMaterialCutsCouple();

   // Kill the particle if flagged by GetMeanFreePath (through boolean killthisone)
   // and also redo the energy check (see GetMeanFreePath for details)
   if (killthisone || kineticEnergy<max(ffermienergy,fbarrier)) {
      if(DI_output) {
        if(!output.is_open()) {
          output.open ("DIout.dat");
          output << std::setprecision(10);
          output << "EventID\tTrackID\tEtotal (eV)\tEkin (eV)"
                 << "\tx (nm)\ty (nm)\tz (nm)\txOrigin (nm)\tyOrigin (nm)\tzOrigin (nm)\tx-direction (nm)\ty-direction (nm)\tz-direction (nm)"
                 << "\tMax Depth (nm)\tMax Radius (nm)" << G4endl;
          }
          output << CADPhysicsOutput::CADPhysicsOutput(track, step) << G4endl;
      }
      aParticleChange.ProposeTrackStatus(fStopAndKill);// Change the particle status
      aParticleChange.ProposeEnergy(0.);// Set energy to zero
      // Treat the remaining kinetic energy as 'deposited locally'
      //aParticleChange.ProposeLocalEnergyDeposit(kineticEnergy - ffermienergy);
      aParticleChange.ProposeLocalEnergyDeposit(kineticEnergy);
      killthisone = false;// Reset the flag
      // Return the default PostStepDoIt function. Doesn't do much extra,
      // except clearing the remaining NumberOfInteractionLengthLeft.
      return G4VContinuousDiscreteProcess::PostStepDoIt(track, step);
   }

   // Sample omegaprime
   // The total cross section is multiplied with a random number. The value for omegaprime is
   // derived from the cumulative differential cross section table.
   G4double omegaprime = inelastic_icdf_vector[findex].get(kineticEnergy/eV,G4UniformRand())*eV;

   // Check for ionization of an outer shell
   G4double Ebind = 0.;// Initialization
   G4int Z=0;
   G4int shellId = 0;
   if (omegaprime>100.*eV) {
      // For large enough omegaprime, first see if the 'standard' Geant4 low energy utilities can
      // come up with a suitable outer shell to be ionized. This will be accepted only if the
      // binding energy is larger than 50 eV. The advantage of this method is that it knows about
      // the relative ionization cross sections of the different outer shells (at a certain
      // omegaprime energy) from which it can randomly select one. A disadvantage is that it knows
      // only about individual atoms, not about solid state electronic structure.
      // Hence the energy limits introduced here.
      G4double bindingEnergy = 0.;
      // Create a margin to catch discrepancies in outer shell energies between G4 data and the
      // df files
      G4double margin = 10.*eV;
      //const G4MaterialCutsCouple* couple = track.GetMaterialCutsCouple();
      // Select an atom (=element) from the mixture
      // G4cout << "DI particle: " << track.GetDefinition()->GetParticleName() << G4endl;
      Z = crossSectionHandler->SelectRandomAtom(couple, omegaprime + margin);
      // Select a shell index from the atom
      G4int shell = crossSectionHandler->SelectRandomShell(Z, omegaprime + margin);
      // N.B.: see CADPhysicsAtomicDeexcitation for the fine distinction between 'shellId's and
      //       shell indices!
      // Get the actual shell given the index...
      const G4AtomicShell* atomicShell = (G4AtomicTransitionManager::Instance())->Shell(Z, shell);
      // ...and get the binding energy
      bindingEnergy = atomicShell->BindingEnergy();
      if(bindingEnergy>50.*eV) {// Accept the binding energy only if it's above 50 eV
         shellId = atomicShell->ShellId();
         Ebind = bindingEnergy;
      }
   }
   // If we haven't found an outer shell (yet), check the data from the df file
   if(Ebind==0 && outershells[findex]) {
      // This part is a bit less advanced. Simply find the largest outer shell ionization energy
      // in the list that is still smaller than omegaprime.
      std::vector<float> itsoutershells = (outershelltable)[findex];
      size_t listsize = itsoutershells.size();
      size_t shellnr=0;
      G4bool toolarge = false;
      do {
         if(omegaprime>itsoutershells[shellnr]) {
            // Binding energy set to this outer shell's energy
            Ebind = itsoutershells[shellnr];
         }
         else {
            // Stop looking for next outer shell
            toolarge = true;
         }
         shellnr++;
      } while (!toolarge && shellnr<listsize);
   }
   // Finished looking for outer shells.

   // Sample omega
   // First set the limits. The maximum energy loss for the primary electron is such that its
   // kinetic energy after the event remains larger than omega-omegaprime, i.e. the *additional*
   // kinetic energy transfer to the secondary electron. omegamax equals the upper integration
   // limit of eq. 9 in Ashley, but corrected for the Fermi energy. The lower limit is set equal
   // to omegaprime (for outer shells) or to eq. 10 in Ashley (otherwise), where the latter
   // takes momentum conservation into account.
   G4double omegamax = 0.5*(kineticEnergy + omegaprime - ffermienergy);
   G4double omegamin = omegaprime;
   G4double Ebindprime = 0.;
   // EK: 2013-06-12 Modified if statement
   // Previous version:
   //    if(Ebind < 50.*eV)
   //      // No 'outer shell' ionization, use Ashley's limit
   // New version:
   if(kineticEnergy > 2.*omegaprime) {
      // Use Ashley's limit, as long as momentum conservation is possible
      omegamin = 0.5*(kineticEnergy+omegaprime-sqrt(kineticEnergy*(kineticEnergy-2.*omegaprime)));
      Ebindprime = omegaprime;
   } else {
      // See below for the meaning of Ebindprime
      Ebindprime = Ebind - ffermienergy;
      if(Ebindprime>omegaprime-1.*eV) Ebindprime=omegaprime-1.*eV;
         // Upper limit necessary because of the margin introduced above
   }

   G4double omega = 0;
   if(omegamin < omegamax && Ebindprime>0.) {
      // For nonzero binding energy, sample omega according to eq. 7 in Ashley,
      // using the lower and upper limits as defined above.
      // For outer-shell ionization (Ebind > 50 eV) we substitute the Fermi-energy corrected
      // binding energy for omegaprime (so that the differential cross section becomes inversely
      // proportional to both the total energy transfer and the kinetic energy of the secondary
      // electron).
      G4double fmin = 1/Ebindprime * log ((omegamin - Ebindprime)/omegamin);
      G4double fmax = 1/Ebindprime * log ((omegamax - Ebindprime)/omegamax);
      G4double B = 1/(fmax - fmin);
      G4double C = -fmin / (fmax - fmin);
      omega = Ebindprime / ( 1 - exp ( Ebindprime * ( G4UniformRand() - C ) / B) );
   } else {
      // In some cases (typically only occuring for Ebind<50 eV) we get omegamin > omegamax.
      // This is due to our fermi Energy correction in the definition of omegamax. Physically, this
      // means that momentum cannot be conserved because the primary electron cannot have a final
      // kinetic energy that is lower than the Fermi energy. In this (relatively rare) case we have
      // to ignore momentum conservation and probe omega according to a 1/(omega)^2 distribution
      // with omegaprime and omegamax as lower and upper limits, respectively.
      omega = omegaprime/(1.-G4UniformRand() *(1.-omegaprime/omegamax));
   }

   // Set the final kinetic energy of the primary particle
   G4double finalKinEnergy = kineticEnergy - omega;

   // Some special work needs to be done for 'zero' binding energy. We distinguish three cases:
   // 1) Metals: excitation of a Fermi sea electron
   // 2) Semiconductors and insulators, and omegaprime larger than the bandgap: cross-bandgap
   //    excitation
   // 3) Insulators: Phonon excitation, with energy loss but no creation of a secondary electron.
   // Now also for semiconductors!
   // Otherwise we simply have zero binding energy.
   if(Ebind==0) {
      if(vec_conductortype[findex]==material::CND_METAL) {
        // Case 1. Metal, assume Fermi sea electron.
         // The probability distribution for the initial energy of the SE is proportional to the
         // density of states of both the initial and the final states of the SE,
         //    p~sqrt(E')*sqrt(E'+omega)
         // where E' is the initial energy above the bottom of the conduction band and omega is the
         // transferred energy.
         // By using this expression we assume a free-electron like behaviour of the electrons in
         // the conduction band (i.e. perfectly quadratic shape of the conduction band).
         // See, e.g., eq. 7 in Y.T. Yue et al., J.Phys.D 38 (2005) 1966-1977
         // (where Delta E is used instead of omega).
         //
         // A normalized cumulative distribution function for the dimensionless variable
         // y=E'/omega is stored in the vector Psecvalues.
         //
         // First limit omega to more than 0.1 and less than 1000 times the Fermi energy.
         // If omega is outside the given range, then Ebind simply remains equal to zero, which
         // effectively means that the SE was assumed to have kinetic energy equal to the Fermi
         // energy. In both extremes this is a reasonable approximation:
         // if omega<0.1*ffermienergy, the error is less than 10% of the final SE's kinetic energy
         // (and in a typical case the SE will have too low energy to escape anyway);
         // if omega>1000*ffermienergy, the error is less than 0.1% of the final SEs kinetic energy.
         if(omega>0.1*ffermienergy && omega<1000.*ffermienergy) {
            // Next, determine the limiting values that the cumulative distribution function can
            // take. These are a function of the following dimensionless variable x:
            G4double x = ffermienergy/omega;
            G4int psecnr=0;
            G4double P0=0;// lower limit
            G4double P1=0;// upper limit
            if(x>1.) {
               // omega is smaller than ffermienergy; the lower limit for E' is
               // ffermienergy - omega. Hence, y is at least x - 1.
               // Here, we find the value of the cumulative distrubtion function that matches this
               // energy value.
               do {
                  if(Psecenergies[psecnr+1]>x-1.) {
                     // We are in the right 'energy bin', do a log-log interpolation to find P0.
                     P0 = Psecvalues[psecnr] * exp(log(x-1.)-log(Psecenergies[psecnr])/
                        (log(Psecenergies[psecnr+1])-log(Psecenergies[psecnr]))
                        *(log(Psecvalues[psecnr+1])-log(Psecvalues[psecnr])));
                        //note: also works fine if x-1<0.001
                     break;
                  }
                  psecnr++;
               } while (psecnr+1<21);
               psecnr--;
            }

            // For any value of x, now find the higher limit of the cumulative distribution function
            // of y.
            // This corresponds to E' equal to the Fermi energy, or y = ffermienergy / omega = x.
            size_t j=0;
            do {
               if(Psecenergies[j+1]>x) {
                  // Do a log-log interpolation to find P1.
                  P1 = Psecvalues[j] * exp(log(x)-log(Psecenergies[j])/
                     (log(Psecenergies[j+1])-log(Psecenergies[j]))
                     *(log(Psecvalues[j+1])-log(Psecvalues[j])));
                  break;
               }
               j++;
            } while (j+1<21);
            G4double P2 = P0 + (P1-P0)*G4UniformRand();
            // Find the y value that matches P2. We can start looking in the bin that corresponds to
            // the lower limit P0.
            G4int k = psecnr;
            do {
               if(Psecvalues[k+1]>P2) {
                  // log-log interpolation of the tabulated cumulative distribution function
                  // E' = omega * y and finally the 'binding' energy is defined as
                  // ffermienergy - E', i.e. the distance of the SE's initial kinetic energy below
                  // the Fermi level.
                  Ebind = ffermienergy - omega*(Psecenergies[k] *
                     exp(log(P2)-log(Psecvalues[k])/(log(Psecvalues[k+1])-log(Psecvalues[k]))
                     *(log(Psecenergies[k+1])-log(Psecenergies[k]))));
                  break;
               }
               k++;
            } while (k+1<21);
         }
      } else {
         if (omegaprime > vec_bandgap[findex]) Ebind = vec_bandgap[findex];
            // Case 2, electron excitation across the bandgap
         else return Phononloss(track,step,kineticEnergy,omegaprime);
      }
   }// End if(Ebind==0)

   // SE kinetic energy is equal to the transferred minus the binding energy -- but counting from
   // the Fermi level, so add the Fermi energy.
   G4double secondaryKineticEnergy = omega - Ebind + ffermienergy;
   G4double theEnergyDeposit = 0;

   // EB: Subtract fermi energy from material for bookkeeping
   theEnergyDeposit += Ebind - ffermienergy;

   // Finally, sample the direction of the secondary electron
   // This part of the code has been derived from that of
   // G4LowEnergyIonisation.
   // 2017-02-03: modification for non-outer-shell secondaries. Momentum transfer from the primary
   // for plasmon excitations is only (by definition) equivalent to omega-omegaprime. deltaKinE
   // modified accordingly. N.B. - to be done: look also into mom. dist. of the secondary. For now,
   // there is no theoretical basis to assume a certain distribution, so we just leave as is.

   // First transform to the shell potential - from both the (initial) primary kinetic energy and
   // the SE kinetic energy first subtract the Fermi energy, then add the binding energy *twice*.
   // The classical physical picture is that the electric potential near the given (outer shell)
   // orbit is twice the binding energy of that shell. The bound electron has a kinetic energy
   // equal to Ebind so that it's total energy is E_pot + E_kin = -2*Ebind + Ebind = -Ebind, as it
   // is expected to be.
   G4double deltaKinE = omega + Ebind;
   if (Ebind<50.*eV) deltaKinE = omega - omegaprime;//2017-02-03
   G4double primaryKinE = kineticEnergy - ffermienergy + 2.0*Ebind;

   // After the transformation we assume that the collision behaves as if it occurs between free
   // electrons of which the SE had no net momentum prior to the collision.
#ifdef CADPHYSICS_RELATIVISTIC
   // Sample the scattering angle neglecting atomic motion, RELATIVISTIC. This is default behavior.
   G4double deltaMom = std::sqrt(deltaKinE*(deltaKinE + 2.0*electron_mass_c2));
   G4double primaryMom = std::sqrt(primaryKinE*(primaryKinE + 2.0*electron_mass_c2));
   G4double cost = deltaKinE * (primaryKinE + 2.0*electron_mass_c2)
      / (deltaMom * primaryMom);
   if (cost > 1.) cost = 1.;
   G4double sint = std::sqrt(1. - cost*cost);
   G4double dirz = cost;
#else
   // Sample the scattering angle neglecting atomic motion, NON-RELATIVISTIC but slightly faster.
   // Here we have implemented a nonrelativistic version for speed of simulation. This leads to
   // max. 1.5% relative error in the cosine of theta at 30 keV primary energy, while the absolute
   // error is even less.
   G4double cost2 = deltaKinE / primaryKinE;
   if (cost2 > 1.) cost2 = 1.;
   G4double cost = sqrt(cost2);
   G4double sint = std::sqrt(1. - cost2);
   G4double dirz = cost;
#endif

   G4double phi  = twopi * G4UniformRand();
   G4double dirx = sint * std::cos(phi);
   G4double diry = sint * std::sin(phi);

   // Rotate to align with the incident electron direction
   G4ThreeVector primaryDirection = track.GetMomentumDirection();
   G4ThreeVector deltaDir(dirx,diry,dirz);
   deltaDir.rotateUz(primaryDirection);
   dirx = deltaDir.x();
   diry = deltaDir.y();
   dirz = deltaDir.z();

   // Finally add a component to the SE's momentum that corresponds to its instantaneous momentum
   // prior to the collision (with random numbers that determine its direction).
   // See G4LowEnergyIonisation for (some) details.
   G4double delcost = 2.0*G4UniformRand() - 1.0;
   sint = std::sqrt(1. - delcost*delcost);
   phi  = twopi * G4UniformRand();
#ifdef CADPHYSICS_RELATIVISTIC
   G4double del = std::sqrt(Ebind *(Ebind + 2.0*electron_mass_c2))// RELATIVISTIC
      / deltaMom;
#else
   G4double del = sqrt(Ebind / deltaKinE);// Again NON-RELATIVISTIC
#endif
   if (Ebind<0.) abort();
   dirx += del* sint * std::cos(phi);
   diry += del* sint * std::sin(phi);
   dirz += del* delcost;
   G4ThreeVector newdir = G4ThreeVector ( dirx, diry, dirz );

   // Get the local safety for the electon range check
   G4StepPoint* steppoint = step.GetPostStepPoint();
   G4double pstepsafety = steppoint->GetSafety();

   // Generate secondary particles

   // Initialization
   std::vector<G4DynamicParticle*>* secondaryVector = 0;
      // Vector to be filled with the secondary particles
   G4DynamicParticle* aSecondary = 0;
   G4ParticleDefinition* type = 0;
   G4double cutForPhotons = 250.*eV;
      // Lower energy limit for X-ray photons
   size_t nSecondaries = 0;
      // Total number of secondary particles BEFORE filtering
   size_t totalNumber  = 0;
      // Total number of secondary particles AFTER filtering

   std::vector<G4double> VacancyEnergies;
      // for hole energies
   newdir = newdir.unit();
   if (generateSecondaries) {
      // X-ray photons and Auger electrons from atom deexcitation
      if (shellId>0 && Z>5) {
         // Data are available only for Z>5 and shellId>0 only if a suitable outer shell was
         // found previously
         secondaryVector = deexcitation.GenerateParticles(Z, shellId);
      } else {
        // If no deexcitation, start with an empty vector
        secondaryVector = new  std::vector<G4DynamicParticle*>;
      }

      // Add the 'normal' secondary electron to the secondary particle vector
      G4DynamicParticle* theDeltaRay = new G4DynamicParticle();
      theDeltaRay->SetKineticEnergy(secondaryKineticEnergy);
      theDeltaRay->SetMomentumDirection(newdir);
      theDeltaRay->SetDefinition(G4Electron::Electron());
      secondaryVector->push_back(theDeltaRay);

      // Loop over all secondary particles (if any) and filter for their energy
      if (secondaryVector != 0) {
         nSecondaries = secondaryVector->size();
         pairsgenerated += nSecondaries;
            // Count all secondaries generated, also those that have too low energy to be tracked
         for (size_t i = 0; i<nSecondaries; i++) {
            aSecondary = (*secondaryVector)[i];
            if (aSecondary) {
               G4double e = aSecondary->GetKineticEnergy();
               type = aSecondary->GetDefinition();
               if (type==G4Electron::Electron()) {
                  if (i<nSecondaries-1) {
                     // Add the Fermi energy to the kinetic energies of all Auger electrons
                     // (i.e. all electrons in the vector but the last)
                     e += ffermienergy;
                     aSecondary->SetKineticEnergy(e);
                     theEnergyDeposit -= e;
                        // Taken relative to the Fermi level. The particle contributes a
                  }
                  // Get the range for the given kinetic energy (which is different for each
                  // particle)
                  G4double electronrange = electron_range_vector[findex].get(e/eV) * nm;
                  // Energy checks.
                  // The electron is accepted as a new particle if either
                  // - its kinetic energy is larger than 1 keV
                  // - or the local safety is less than five times its range (so that it has a
                  //   realistic chance of reaching a surface)
                  // - or rangecut is false
                  // - or the material is an insulator
                  // and its kinetic energy equals at least the vacuum boundary level plus
                  // the user-determined cut energy (if any).
                  G4bool passenergy = false;
                  if (!passenergy) {
                     passenergy = (e>1.*keV);
                     if (!passenergy) {
                        passenergy = (e>max(ffermienergy,fbarrier) &&
                                       (!rangecut_vector[findex] || pstepsafety < 5.*electronrange));
                     }
                  }
                  if (passenergy) {
                     //Only check for zero energy, to avoid G4VParticleChange errors
                     //theEnergyDeposit -= e;
                     //   // Taken relative to the Fermi level. The particle contributes a
                     //   // *negative* energy deposit because it takes a share of the binding
                     //   // energy that is not deposited into the material.
                     totalNumber++;
                  } else {
                     // The particle did not pass the filter, so delete it and remove its pointer
                     // from the vector.
                     delete aSecondary;
                     (*secondaryVector)[i] = 0;
                     theEnergyDeposit += e; //EB: The energy is returned to the material
                  }
               } else if (type == G4Gamma::Gamma() && e > cutForPhotons && generateXrays) {
                  // Accept an X-ray photon...
                  // if the user switch is on and its energy is larger than cutForPhotons.
                  // This number is fixed to 250 eV since this is the lower limit that Geant4's
                  // LowEnergy gamma processes are still valid for.
                  theEnergyDeposit -= e;
                  totalNumber++;
               } else {
                  // Delete the particle in all other cases
                  delete aSecondary;
                  (*secondaryVector)[i] = 0;
               }
            }
         }
      }
   } else {
      // EB: If we are not creating secondaries the energy is considered to be deposited here
      theEnergyDeposit += secondaryKineticEnergy;
   }

   // Register the number of secondary particles that remain after filtering
   aParticleChange.SetNumberOfSecondaries(totalNumber);

   // Register the actual secondary particles that remain in the vector after filtering
   G4double SE_energy = 0.0;
   G4int    SE_number = 0;
   if (secondaryVector) {
      for (size_t l = 0; l < nSecondaries; l++) {
         aSecondary = (*secondaryVector)[l];
         if(aSecondary) {
            aParticleChange.AddSecondary(aSecondary);
            SE_energy += aSecondary->GetKineticEnergy();
            SE_number++;
         }
      }
      delete secondaryVector;
   }
   // End of generating secondary particles

   // Filter for the primary electron energy (similar to the SE filter above)
   G4double electronrange = electron_range_vector[findex].get(finalKinEnergy/eV) * nm;
   if ( finalKinEnergy > max(ffermienergy,fbarrier) &&
         (!rangecut_vector[findex] || pstepsafety < 5.*electronrange || finalKinEnergy > 1.*keV)) {
      // The particle survives
      // Update the primary electron direction assuming conservation of momentum
      G4ThreeVector finalP = primaryDirection - cost * newdir;
      finalP = finalP.unit();
      aParticleChange.ProposeMomentumDirection(finalP);
   } else {
      // The particle does not survive, prepare to kill it
      theEnergyDeposit += finalKinEnergy;
         // EB: No longer Taken relative to the Fermi level.
      finalKinEnergy    = 0.0;
      if(DI_output) {
        if(!output.is_open()) {
          output.open ("DIout.dat");
          output << std::setprecision(10);
          output << "EventID\tTrackID\tEtotal (eV)\tEkin (eV)"
                 << "\tx (nm)\ty (nm)\tz (nm)\txOrigin (nm)\tyOrigin (nm)\tzOrigin (nm)\tx-direction (nm)\ty-direction (nm)\tz-direction (nm)"
                 << "\tMax Depth (nm)\tMax Radius (nm)" << G4endl;
        }
        output << CADPhysicsOutput::CADPhysicsOutput(track, step) << G4endl;
      }
      aParticleChange.ProposeTrackStatus(fStopAndKill);
   }

   // EB: Commented out the following 4 lines. The deposit must be allowed to be negative in order
   //     to make the balance correct
   // if(theEnergyDeposit < 0.) {// Sanity check. We can only deposit net energy, not extract it.
   //    theEnergyDeposit = 0.0;
   // }
   // G4cout << "---------------------------------------------------------------------------"
   //        << G4endl;
   // G4cout << setprecision(9) << "DI: Start Kinetic energy =          " << kineticEnergy/eV
   //        << G4endl;
   // G4cout << setprecision(9) << "DI: Final Kinetic energy =          " << finalKinEnergy/eV
   //        << G4endl;
   // G4cout << setprecision(9) << "    Kinetic energy of secondaries = " << SE_energy/eV
   //        << G4endl;
   // G4cout << setprecision(9) << "    Energy deposit =                " << theEnergyDeposit/eV
   //        << G4endl;
   // G4cout << setprecision(9) << "    Total =                         "
   //        << (finalKinEnergy+SE_energy+theEnergyDeposit)/eV << G4endl ;
   // G4cout << "---------------------------------------------------------------------------"
   //        << G4endl;
   aParticleChange.ProposeLocalEnergyDeposit(theEnergyDeposit);
   aParticleChange.ProposeEnergy(finalKinEnergy);

   return G4VContinuousDiscreteProcess::PostStepDoIt(track, step);
}

G4double CADPhysicsDI::GetMeanFreePath( const G4Track& aTrack,
   G4double, G4ForceCondition* condition)
{
   // Initialization
   *condition = NotForced;
   const G4MaterialCutsCouple* couple = aTrack.GetMaterialCutsCouple();
   G4double kineticEnergy = aTrack.GetKineticEnergy();
   G4int index = couple->GetIndex();
   G4String mname = couple->GetMaterial()->GetName();

   G4int stepnr = aTrack.GetCurrentStepNumber();
   //// Temporarily enlarged to a factor 10 higher in order to follow particles longer
   //if (stepnr>=100000 || kineticEnergy > 1.)
   if (stepnr>=10000 || kineticEnergy > 1.) {
      *condition = Forced;
      killthisone = true;
      return DBL_MIN;
   }

   findex = index;

   ffermienergy = vec_fermieff[findex];
   fbarrier = vec_barrier[findex];
   //fconductortype = vec_conductortype[findex];

   G4double mfp = DBL_MAX;

   if (mname == "Galactic") { // prevent the read out of inelastic_imfp_vector for Galactic
      return mfp;
   }

   G4double cross = inelastic_imfp_vector[findex].get(kineticEnergy/eV); // nm^-1
   if (cross > 0.) {
      mfp = 1./cross; // nm
      mfp = mfp * 1.e-6; // mm
   }
   return mfp;
}

void CADPhysicsDI::PrintInfoDefinition()
{
   G4cout << "CADPhysicsDI: Process for inelastic scattering in (non-gaseous) materials,\n"
      << " based on the Dielectric Ionisation formalism. The implementation is loosely based on\n"
      << " the description in J.C. Ashley, J. Electr. Spectrosc. Rel. Phenom. 46: 199-214 (1988),\n"
      << " but with adjusted cross sections for ionization of outer shells. Furthermore, Auger\n"
      << " decay of outer-shell ionized atoms has been included.\n"
      << " The process is valid for energies up to 30 keV." << G4endl;
}
