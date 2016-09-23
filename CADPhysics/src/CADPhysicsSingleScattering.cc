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

using namespace std;

CADPhysicsSingleScattering* CADPhysicsSingleScattering::fSingleScattering = 0;

CADPhysicsSingleScattering::CADPhysicsSingleScattering(const G4String& processName)
: G4VContinuousDiscreteProcess(processName),
theInvMFPTable(0),
theTMFPTable(0),
theDInvMFPTable(0),
theSumDInvMFPTable(0),
currentCouple(0),
currentMaterial(0),
enlossconst (2.*electron_mass_c2*Avogadro/c_squared),
previousenergy(0.),
previousMFP(0.),
preStepInvMFP(0.),
energybin(0),
PSDInvMFPforEl(0),
PSDInvMFPforEnergy1(0),
PSDInvMFPforEnergy2(0),
InvMFPforMat(0),
fTMFPforMat(0),
domultistep(true),
dodiffusionstep(true),
workforgases(true),
numangles(101),
numMottangles(96)
{
   SetVerboseLevel(1);
   messenger = new CADPhysicsSSMessenger(this);
   singlestep = false;
   multistep = false;
   diffusionstep = false;
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
{
   if(theDInvMFPTable)
   {
      CADPhysicsDataCube* a=0;
      CADPhysicsDataTable* b=0;
      G4DataVector* c=0;
      while (theDInvMFPTable->size()>0)
      {
         a = theDInvMFPTable->back();
         theDInvMFPTable->pop_back();
         while (a->size()>0)
         {
            b = a->back();
            a->pop_back();
            while (b->size()>0)
            {
               c = b->back();
               b->pop_back();
               if ( c ) {delete c;}
            }
            delete b;
         }
         delete a;
      }
      delete theDInvMFPTable;
   }
   if (theInvMFPTable)
   {
      G4DataVector* a=0;
      while (theInvMFPTable->size()>0)
      {
         a = theInvMFPTable->back();
         theInvMFPTable->pop_back();
         if ( a ) {delete a;}
      }
      delete theInvMFPTable;
   }
   if (theTMFPTable)
   {
      G4DataVector* a=0;
      while (theTMFPTable->size()>0)
      {
         a = theTMFPTable->back();
         theTMFPTable->pop_back();
         if ( a ) {delete a;}
      }
      delete theTMFPTable;
   }
}

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
   if( theInvMFPTable != 0 ) delete theInvMFPTable;
   theInvMFPTable = new CADPhysicsDataTable;
   if( theTMFPTable != 0 ) delete theTMFPTable;
   theTMFPTable = new CADPhysicsDataTable;
   if( theDInvMFPTable != 0 ) delete theDInvMFPTable;
   theDInvMFPTable = new std::vector<CADPhysicsDataCube*>;
   if( theSumDInvMFPTable != 0 ) delete theSumDInvMFPTable;
   theSumDInvMFPTable = new CADPhysicsDataCube;
   vec_phononloss.clear();

   BuildRanges();// Fill the enrange and anglerange vectors

   size_t numOfCouples = theCoupleTable->GetTableSize();
   for(size_t i=0; i<numOfCouples; i++)// Loop over all materials
   {
      const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(i);
      G4int cond=0;// Conductor type (see below)
      const G4Material* material = couple->GetMaterial();
      G4MaterialPropertiesTable* aMPT = material->GetMaterialPropertiesTable();
      // Get material properties. An overview of those properties that are used here:
      // Conductortype - type of material. 0 stands for metals, 1 for semiconductors and 2 for insulators.
      // For METALS:
      //      Fermi energy - kinetic energy of the electron at the Fermi level (in Geant4 internal units)
      //      Resistivity - resistivity of the material (in units of Ohm meter)
      // For SEMICONDUCTORS/INSULATORS:
      //      Sound velocity (in units of meters/second)
      //      DEFPOTENTIAL = Acoustic deformation potential (in Geant4 internal units)
      //      LATTICE = lattice constant (in units of nm)
      if (aMPT) {
         if(aMPT->ConstPropertyExists("CONDUCTORTYPE")) {
            cond = (int)(aMPT->GetConstProperty("CONDUCTORTYPE"));
         }
      }

      // Define and initialize material parameters
      G4bool issemi = false;
      G4bool ismetal = false;
      G4double soundvelocity = 0.;// Sound velocity
      G4double defpotential = 0.;// Acoustic deformation potential
      G4double lattice = 0.;// Lattice constant
      G4double fermienergy = 0.;// Fermi energy (i.e. kinetic energy at the Fermi level)
      G4double resistivity = 0.;// Resistivity
      fphononloss = 0.00108*eV;// Average energy loss in case of an acoustic phonon scattering event. Default is the value for silicon.
      if (aMPT) {
         // Check if we have sufficient input parameters for an explicit description of the acoustic phonon scattering.
         // The approach here is: we treat the material as a semiconductor simply if we have enough information to do so -
         // regardless of whether the material was actually defined to be a semiconductor (or insulator).
         // The required parameters are, as stated above, the sound velocity, deformation potential, and lattice constant.
         if(aMPT->ConstPropertyExists("SOUNDVELOCITY")) {
            soundvelocity = aMPT->GetConstProperty("SOUNDVELOCITY");
            if(aMPT->ConstPropertyExists("LATTICE")) {
               lattice = aMPT->GetConstProperty("LATTICE");
               if (lattice>0.) {// Enough data to evaluate the phonon energy loss.
                  fphononloss = 6.463e-8*eV*soundvelocity/lattice;
                  if (fphononloss>0.05*eV) {
                     fphononloss = 0.05*eV;
                     G4cout << "CADPhysicsSingleScattering: Warning: " << G4endl;
                     G4cout << "Limiting acoustic phonon energy constant for material ";
                     G4cout << material->GetName() << " to 0.05 eV, to avoid unphysical results." << G4endl;
                  }
               }
               if(aMPT->ConstPropertyExists("DEFPOTENTIAL")) {
                  defpotential = aMPT->GetConstProperty("DEFPOTENTIAL");
                  issemi = true;// All parameters found.
               }
            }
         }
      }
      vec_phononloss.push_back(fphononloss);

      if (!issemi && aMPT) {
         // Second option (for metals): derive the acoustic phonon scattering inverse mean free path (at the low-energy limit)
         // from the metal's resistivity
         if(aMPT->ConstPropertyExists("FERMIENERGY")) {
            fermienergy = aMPT->GetConstProperty("FERMIENERGY");
            if(aMPT->ConstPropertyExists("RESISTIVITY")) {
               resistivity = aMPT->GetConstProperty("RESISTIVITY");
               ismetal = true;
            }
         }
      }
      if (!issemi && !ismetal && material->GetState() != kStateGas && material->GetName() != "Water") {
         // Both descriptions fail, and the material is not defined as a gas (or Water). Then a warning message is printed, but we continue
         // using default parameters (for silicon).
         // This provides a way to avoid defining those parameters for simulations where one is not interested in any low-energy behavior.
         G4cout << "CADPhysicsSingleScattering: found insufficient GDML data to calculate cross sections" << G4endl;
         G4cout << "for material " << material->GetName() << ". Insert either fermienergy and resistivity" << G4endl;
         G4cout << "or soundvelocity, defpotential and lattice." << G4endl;
         G4cout << "The program continues with standard values for silicon - but update your GDML for realistic output!" << G4endl;
         issemi = true;
         soundvelocity = 9040.;
         defpotential = 9.2*eV;
         lattice = 0.543;
      } else if (material->GetState() == kStateGas || material->GetName() == "Water") {
         // In the case of gases, it does not make sense to implement an 'acoustic phonon' type of scattering.
         // Instead, we just extrapolate the Mott scattering cross sections to zero energy.
         fermienergy = 0.;
         if(verboseLevel>1 ) G4cout << "CADPhysicsSingleScattering: " << material->GetName() << " is treated as a gas, " <<
            "with Mott cross sections starting from zero energy." << G4endl;
      } else if(verboseLevel>1) {
         if (issemi)
         {
            G4cout << "CADPhysicsSingleScattering: " << material->GetName() << " has" << G4endl;
            G4cout << "\tsound velocity =            " << soundvelocity << " m/s;" << G4endl;
            G4cout << "\tac. deformation potential = " << defpotential/eV << " eV;" << G4endl;
            G4cout << "\tlattice constant =          " << lattice << " nm." << G4endl;
         } else {
            G4cout << "CADPhysicsSingleScattering: " << material->GetName() << " has resistivity = " << resistivity << " Ohm meter." << G4endl;
         }
      }

      if (issemi) {
         if(1 < verboseLevel) G4cout << "Run DInvMFPTableforSemi" << G4endl;
         CADPhysicsDataCube* DInvMFPMat = DInvMFPTableforSemi(couple,soundvelocity,defpotential,lattice);// Calculate differential inverse MFPs
         theDInvMFPTable -> push_back(DInvMFPMat);
      } else {
         if(1 < verboseLevel) G4cout << "Run DInvMFPTableforMetal" << G4endl;
         CADPhysicsDataCube* DInvMFPMat = DInvMFPTableforMetal(couple,fermienergy,resistivity);// Calculate differential inverse MFPs
         theDInvMFPTable -> push_back(DInvMFPMat);
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

void CADPhysicsSingleScattering::BuildRanges()// Fill the enrange and anglerange vectors
{
   if(verboseLevel>1) G4cout << "CADPhysicsSingleScattering: Start BuildRanges\n";
   enrange.clear();
   G4double tempenrange[38] = {// Energies from 200 eV up to 30 keV match those of the tabulated Mott cross sections
      200.*eV, 300.*eV, 400.*eV, 500.*eV,
      600.*eV, 700.*eV, 800.*eV, 900.*eV,
      1.*keV, 2.*keV, 3.*keV, 4.*keV, 5.*keV,
      6.*keV, 7.*keV, 8.*keV, 9.*keV, 10.*keV,
      15.*keV, 20.*keV, 25.*keV, 30.*keV,
      38168.*eV, 48559.*eV, 61780.*eV,
      78600.*eV, 100.*keV, 158489.*eV, 251189.*eV,
      398107.*eV, 630957.*eV, 1.*MeV, 1584893.*eV,
      2511886.*eV, 3981072.*eV, 6309573.*eV, 10.*MeV
   };
   G4double logstep = pow(10.,0.1);
   G4double en = 0.01*eV;
   for (size_t i=0;i<41;i++) {
      enrange.push_back(en);
      en *= logstep;
   }
   for(size_t j=0;j<38;j++) enrange.push_back(tempenrange[j]);
   const G4double piover180 = pi / 180.;
   G4double tempanglerange[101] = {// These angles match the tabulated values for the Mott differential cross sections, with five additions
      0.05, 0.3, 0.55, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5,
      11., 13., 15., 17., 19., 21., 23., 25., 27., 29., 31., 33., 35., 37., 39.,
      41., 43., 45., 47., 49., 51., 53., 55., 57., 59., 61., 63., 65., 67., 69.,
      71., 73., 75., 77., 79., 81., 83., 85., 87., 89., 91., 93., 95., 97., 99.,
      101., 103., 105., 107., 109., 111., 113., 115., 117., 119.,
      121., 123., 125., 127., 129., 131., 133., 135., 137., 139.,
      141., 143., 145., 147., 149., 151., 153., 155., 157., 159.,
      161., 163., 165., 167., 169., 171., 173., 175., 177., 179., 180.};
   for(size_t jj=0;jj<numangles;jj++)
   {
      // Note: we store the cosines of the angles rather than the angles themselves!
      tempanglerange[jj] = std::cos(tempanglerange[jj]*piover180);
      anglerange.push_back(tempanglerange[jj]);
   }
}

CADPhysicsDataCube* CADPhysicsSingleScattering::DInvMFPTableforSemi(const G4MaterialCutsCouple* couple, G4double soundvelocity,
                                                G4double defpotential, G4double lattice)
                                                // Calculate (differential) inverse MFPs for semiconductors and insulators
{
   G4double e_mass_kg = electron_mass_c2/(c_squared*kg);
   G4double k_bz = pi/(lattice*nanometer/meter);
   G4double hbar_k_bz = (hbar_Planck/(joule*s))*k_bz;
   G4double w_bz = k_bz*soundvelocity;
   G4double EBZ = (eV/MeV)*(1.0/e_SI)*pow(hbar_k_bz,2)/(2*e_mass_kg);
   G4double screening = 5. * EBZ;// Screening parameter (see Fitting et al.)

   const G4Material* material = couple->GetMaterial();
   const size_t numel = material->GetNumberOfElements();

   G4double totaldensity = material->GetDensity();

   // Note: for calculation of the following prefactors, all parameters need to be given in the units as stated above!
   G4double preconst1 = (e_mass_kg*e_mass_kg)*(k_Boltzmann/(joule/kelvin))*(STP_Temperature/kelvin+25.)*pow((defpotential/joule),2.) /
      (pi*pow((hbar_Planck/(joule*s)),4.)*pow(soundvelocity,2.)*totaldensity*meter3/kilogram); // SI units
   preconst1 = 1.e-6*(millimeter/meter)*preconst1; // G4 internal units
   preconst1 = preconst1 * pi; // Multiply the preconst1 with pi to get the same result as in Fitting 2001
   // Constant prefactor for inverse mean free path at energy below EBZ/4

   G4double hbar_5 = pow((hbar_Planck/(joule*s)),5.);
   G4double e_mass_3 = pow(e_mass_kg,3.);
   G4double nB = 1./(exp((hbar_Planck/(joule*s)*w_bz)/((k_Boltzmann/(joule/kelvin))*(STP_Temperature/kelvin+25.)))-1.); // number density
   G4double preconst2 = 4.*(2.*nB+1.)*e_mass_3*pow((defpotential/joule),2.) /
      (pi*hbar_5*w_bz*totaldensity*(meter3/kilogram)); // SI units
   preconst2 = 1.e-6*(millimeter/meter)*(eV/joule)*preconst2; // G4 internal units and taking into account that the energy will be used in eV further down
   preconst2 = preconst2 * pi; // Multiply the preconst2 with pi to get the same result as in Fitting 2001
   // Constant prefactor for inverse mean free path at energy above EBZ

   G4String mname = material->GetName();

   const G4ElementVector* theElementVector = material->GetElementVector();
   const G4double* theAtomNumDensityVector = material->GetVecNbOfAtomsPerVolume();
   CADPhysicsDataCube* DInvMFPforMat = new CADPhysicsDataCube;
   CADPhysicsDataTable* SigmaThetaforMat = new CADPhysicsDataTable;// Used later to fill the tmfp file;
   // SigmaTheta is the inverse mean free path weighted with 1-cos(theta) where theta is the scattering angle. It is used to
   // calculate the average distance an electron has to travel for 'isotropization' of its direction by elastic scattering.

   G4int integralZ;
   G4bool fileExists = false;
   for (size_t ii=0; ii<numel; ii++) {// Loop over all elements in the material
      CADPhysicsDataTable* DInvMFPforEl = new CADPhysicsDataTable;
      G4DataVector* SigmaThetaforEl = new G4DataVector;

      // First the phonon part:
      // N.B.: this implementation is not yet very elegant. We do here for each element what needs to be done only once
      // for the material... but it works.
      size_t j=0;
      for (j=0; j<41; j++) {// Energies up to 100 eV: do pure phonon scattering
         G4double en = enrange[j];

         G4double interpolfac = 0.;// Factor for interpolation between the 'low-energy' and 'higher-energy' inverse MFPs, cf
         // Fitting et al., eq. (3a) and (3b). For E<E_BZ/4, the low-energy InvMFP is used...
         if (en > EBZ) interpolfac = 1.; else // Above E>E_BZ, we use the higher-energy InvMFP
            if (en > EBZ/4.) interpolfac = (en-EBZ/4.)/(EBZ-EBZ/4.); // And we interpolate for intermediate energies.

         // Low-energy inverse mean free path, cf eq (3a) in Fitting et al., assuming a parabolic energy band
         G4DataVector* DInvMFPforEnergy = new G4DataVector;
         G4double lowenInvMFP = preconst1 / ( 1. + en/screening);
         if (mname=="Galactic") lowenInvMFP=DBL_MIN;
         lowenInvMFP /= numel*nanometer;// Distribute the total InvMFP evenly over all elements in the compound - see above remark
         G4double preconst3 = en + screening;// A factor for renormalization of the angular distribution
         lowenInvMFP *= preconst3;

         // Higher-energy inverse mean free path, cf eq (3b) in Fitting et al.
         // assuming a parabolic energy band and m*=m_e
         G4double highenInvMFP = preconst2 * screening * screening / en *
            ( - en/(screening * (1. + en/screening)) + log(1. + en/screening));
         G4double norm = - en / (en + screening) + log ( 1. + en/screening);// Normalization factor for the angular distribution
         // of the inverse mean free path at higher energies
         if (mname=="Galactic") highenInvMFP=DBL_MIN;
         highenInvMFP /= numel*nanometer*eV;// Again, distribute the total InvMFP evenly over all elements in the compound
         highenInvMFP /= norm;

         // Calculate cumulative differential inverse mean free path for angle
         G4double lastvalue = 0.;
         G4double cumSigmaThetaInvMFP = 0.;
         for(size_t jj=0;jj<numangles;jj++)// Loop over all predefined angles
         {
            G4double fac = 2. * screening / ( 2. * screening + en - en * anglerange[jj]);
            G4double cumDInvMFP = (1-interpolfac) * lowenInvMFP / ( en + screening / ( 0.5 - 0.5 * anglerange[jj] ))// Low-energy approx.
               + interpolfac * highenInvMFP * ( -1. + fac - log(fac) );// Higher-energy approx.
            DInvMFPforEnergy->push_back(cumDInvMFP);
            if (jj>0) {
               cumSigmaThetaInvMFP += (cumDInvMFP - lastvalue) * (1.-(anglerange[jj] + anglerange[jj-1])/2);
            } else cumSigmaThetaInvMFP += cumDInvMFP * (1.-anglerange[jj]);
            lastvalue = cumDInvMFP;
         }
         DInvMFPforEl->push_back(DInvMFPforEnergy);
         SigmaThetaforEl->push_back(cumSigmaThetaInvMFP);
      }

      const G4Element* elm = (*theElementVector)[ii];
      G4double density = theAtomNumDensityVector[ii];
      integralZ = (int)std::floor(elm->GetZ()+0.5);// Convert the element's 'effective' Z value to an integer number

      // For energies above 100 eV, we first need to open the file containing the Mott cross section data
      // Note: there is a 'clear cut' in the tabulated data: up to 100 eV we have phonon scattering, and we use Mott cross sections
      // from 200 eV upward. During simulations, the (angular) inverse mean free paths are linearly interpolated between 100 and 200 eV
      std::ostringstream ost;
                char *path=getenv("CADPHYSICS_BASE");
                if (!path) {
         if(integralZ<10) ost << "/usr/local/share/CAD/mott/0" << integralZ << ".dat";
                   else ost << "/usr/local/share/CAD/mott/" << integralZ << ".dat";
                } else {
         if(integralZ<10) ost << path << "/mott/0" << integralZ << ".dat";
                   else ost << path << "/mott/" << integralZ << ".dat";
                }
      G4String name = ost.str();
      if(verboseLevel>1) G4cout << "Trying to open file for element " << integralZ << "... ";
      std::ifstream file(name);
      std::filebuf* lsdp = file.rdbuf();
      if (lsdp->is_open())
      {
         fileExists = true;
         if(verboseLevel>1) G4cout << "Opened Mott CS file for element " << integralZ << "." << G4endl;
      } else G4cout << "CADPhysicsSingleScattering: opening Mott CS file for element " << integralZ << " failed! Using zero CS instead!" << G4endl;
      G4double a = 0.;

      if(fileExists)
      {
         for (size_t dummy = 0; dummy<4; dummy++)// Skip entries for first four energy values (25, 50, 75, and 100 eV)...
         {
            for(size_t jj=0;jj<numMottangles;jj++) file >> a;
            file >> a;
            file >> a;
         }
         for(j=41; j<78; j++)// Read data for energies from 200 eV up to 30 keV... Rutherford formula for the higher energies
         {
            G4DataVector* DInvMFPforEnergy = new G4DataVector;
            G4double cumInvMFP = 0.;
            G4double InvMFP = 0.;
            G4double cumSigmaThetaInvMFP = 0.;
			G4double tmpa = 0.;
			G4double preva = 0.;
            for(size_t jj=0;jj<numangles;jj++)// ...and for all angles in the range.
            {
               if (j<63) {
				   if (jj>8) file >> a;
				   else {
					   if (0==jj || 3==jj || 5==jj || 7==jj) file >> tmpa;// Read a new value from file...
					   if (0==jj || 4==jj || 6==jj || 8==jj) a = tmpa;
					   else if (3==jj || 5==jj || 7==jj)  a = (a + tmpa)/2.;
				   }
               } else {
                  a = ComputeRutherfordCrossSectionPerAtom(enrange[j],integralZ,acos(anglerange[jj])); // Compute using screened relativistic Rutherford crossection
               }
               if(jj==0) {
                  InvMFP = a * twopi * (1. - anglerange[jj]) * density * 1.e-14;
                  // Converting from the CS per atom given in sq. Angstrom to an inverse MFP (in mm^-1)
                  cumInvMFP += InvMFP;// Build the cumulative angular-differential inverse MFP vector
                  cumSigmaThetaInvMFP += InvMFP * (1. - anglerange[jj]);
               }
               else {
                  InvMFP = (a + preva)/2. * twopi * (anglerange[jj-1] - anglerange[jj]) * density * 1.e-14;
                  // Converting from the CS per atom given in sq. Angstrom to an inverse MFP (in mm^-1)
                  cumInvMFP += InvMFP;// Build the cumulative angular-differential inverse MFP vector
                  cumSigmaThetaInvMFP += InvMFP * (1. - (anglerange[jj]+ anglerange[jj-1])/2.);
               }
               DInvMFPforEnergy->push_back(cumInvMFP);
			   preva = a;
            }
            DInvMFPforEl->push_back(DInvMFPforEnergy);
            SigmaThetaforEl->push_back(cumSigmaThetaInvMFP);
            if (j<63) {
               file >> a;
               file >> a;// Skip the cumulative cross section and mean free path for each energy
               // Instead, we recalculate the total cross section/inverse mean free path from the differential function (see below),
               // to avoid rounding errors
            }
         }
      } else {// The Mott file does not exist
         // Seems we have a problem, since the Mott files should be present for all elements up to number 94.
         G4cerr << "CADPhysicsSingleScattering: Mott file for element with Z=" << integralZ << " not found. Exiting now." << G4endl;
         exit(1);
         // Alternatively, one could just use zero cross sections but presumably this is not what the user really wants...
      }
      file.close();
      DInvMFPforMat->push_back(DInvMFPforEl);
      SigmaThetaforMat->push_back(SigmaThetaforEl);
   }
   // Calculate the total inverse MFP for each energy (sum up the highest cumulative inverse MFPs for all elements in the material)
   G4DataVector* InvMFPforMat = new G4DataVector;
   CADPhysicsDataTable* SumDInvMFPforMat = new CADPhysicsDataTable;

   if (verboseLevel>2) {
      G4cout << "Total mean free paths:" << G4endl;
      G4cout << "\tEnergy(eV)\tEMFP(nm)" << G4endl;
   }

   for(size_t j=0;j<78;j++)
   {
      G4double totinvmfp = 0.;
      for(size_t ii=0;ii<numel;ii++)
      {
         CADPhysicsDataTable* DInvMFPforEl = (*DInvMFPforMat)[ii];
         G4DataVector* DInvMFPforEnergy = (*DInvMFPforEl)[j];
         size_t vectorsize = DInvMFPforEnergy->size();
         totinvmfp += (*DInvMFPforEnergy)[vectorsize-1];
      }
      if (verboseLevel>2) G4cout << mname << "\t" << enrange[j]/eV << "\t" << 1/(totinvmfp*nanometer) << G4endl;
      InvMFPforMat->push_back(totinvmfp);
   }

   // Now also calculate the total cumulative differential inverse MFP per material.
   for(size_t j=0;j<78;j++)
   {
      G4DataVector* SumDInvMFPforEnergy = new G4DataVector;
      for(size_t jj=0;jj<numangles;jj++)
      {
         G4double SumDInvMFP = 0;
         for (size_t ii=0;ii<numel;ii++)
         {
            SumDInvMFP += (*(*(*DInvMFPforMat)[ii])[j])[jj];
         }
         SumDInvMFPforEnergy->push_back(SumDInvMFP);
      }
      SumDInvMFPforMat->push_back(SumDInvMFPforEnergy);
   }
   theSumDInvMFPTable->push_back(SumDInvMFPforMat);
   theInvMFPTable->push_back(InvMFPforMat);

   // Now calculate transport mean free path values and export to file...
   // These are the average distances the electron has to travel before 'isotropization' of its direction due to elastic scattering.
   // These values are used by the DI process to determine the effective distance electrons can travel inside the material
   // and hence decide if a certain electron has any chance to leave the current material, and therefore whether or not
   // it needs to be created or killed.
   // Also, the transport mean free path is used in the DoDiffusionStep method to determine the position and momentum distribution at
   // the end of a step.

   std::ostringstream ost;
   ost << "tmfp_" << mname << ".dat";
   G4String name = ost.str();

   G4DataVector* TMFPforMat = new G4DataVector;
   std::ofstream tmfpfile;
   tmfpfile.open (name.c_str());
   tmfpfile << setprecision(6);
   for(size_t j=0; j<78; j++) {
      G4double invtmfp = 0.;
      for(size_t ii=0;ii<numel; ii++) invtmfp += (*(*SigmaThetaforMat)[ii])[j];
      G4double tmfp = 1./invtmfp;
      tmfpfile << enrange[j] << "\t" << tmfp << G4endl;
      TMFPforMat->push_back(tmfp);
   }
   tmfpfile << "-1\t-1\n";
   tmfpfile.close();
   theTMFPTable->push_back(TMFPforMat);

   return DInvMFPforMat;
}

CADPhysicsDataCube* CADPhysicsSingleScattering::DInvMFPTableforMetal(const G4MaterialCutsCouple* couple, G4double fermienergy,
                                                 G4double resistivity)
                                                 // Calculate (differential) inverse MFPs for metals (and gases)
                                                 // Note: the Mott/Browning part is identical to that of DInvMFPTableforSemi,
                                                 // but the implementation for acoustic phonon scattering is different.
{
   G4double screening = 25.*eV;// Screening parameter (see Fitting et al.). Here we do not have information on the lattice
   // constant (and hence E_BZ) available, therefore we just define it as a constant number.

   if(verboseLevel>1) G4cout << "CADPhysicsSingleScattering: calculate DInvMFPTableforMetal ";

   // Find the energy bin for the Fermi energy
   G4int fermibin=0;
   while (enrange[fermibin]<fermienergy && fermibin<41) fermibin++;// The Fermi energy is limited to 100 eV, for safety

   const G4Material* material = couple->GetMaterial();
   G4String mname = material->GetName();
   if(verboseLevel>1) G4cout << "for material " << mname << G4endl;
   const size_t numel = material->GetNumberOfElements();
   const G4ElementVector* theElementVector = material->GetElementVector();
   const G4double* theAtomNumDensityVector = material->GetVecNbOfAtomsPerVolume();
   CADPhysicsDataCube* DInvMFPforMat = new CADPhysicsDataCube;
   CADPhysicsDataTable* SigmaThetaforMat = new CADPhysicsDataTable;// Used later to fill the tmfp file;
   // SigmaTheta is the inverse mean free path weighted with 1-cos(theta) where theta is the scattering angle. It is used to
   // calculate the average distance an electron has to travel for 'isotropization' of its direction by elastic scattering.

   G4int integralZ;
   G4bool fileExists = false;
   for (size_t ii=0; ii<numel; ii++) {// Loop over all elements in the material
      CADPhysicsDataTable* DInvMFPforEl = new CADPhysicsDataTable;
      G4DataVector* SigmaThetaforEl = new G4DataVector;

      // For metals, we use the phonon inverse MFPs only at (or up to) the Fermi energy, and interpolate linearly from this value to
      // the full Mott inverse MFPs at 100 eV.
      // Hence, we need to start with reading the Mott file entirely.
      CADPhysicsDataTable* MottTable = new CADPhysicsDataTable;
      const G4Element* elm = (*theElementVector)[ii];
      integralZ = (int)std::floor(elm->GetZ()+0.5);

      // Open the Mott cross section file for this element
      std::ostringstream ost;
                char *path=getenv("CADPHYSICS_BASE");
                if (!path) {
         if(integralZ<10) ost << "/usr/local/share/CAD/mott/0" << integralZ << ".dat";
                   else ost << "/usr/local/share/CAD/mott/" << integralZ << ".dat";
                } else {
         if(integralZ<10) ost << path << "/mott/0" << integralZ << ".dat";
                   else ost << path << "/mott/" << integralZ << ".dat";
                }
      G4String name = ost.str();
      if(verboseLevel>1) G4cout << "Trying to open file for element " << integralZ << "... ";
      std::ifstream file(name);
      std::filebuf* lsdp = file.rdbuf();
      if (lsdp->is_open())
      {
         fileExists = true;
         if(verboseLevel>1) G4cout << "Opened Mott CS file for element " << integralZ << "." << G4endl;
      } else G4cout << "CADPhysicsSingleScattering: opening Mott CS file for element " << integralZ << " failed! Using zero CS instead!" << G4endl;

      G4double a = 0.;

      if (fileExists)
      {
         for (size_t dummy = 0; dummy<26; dummy++)
         {
            G4DataVector* MottVector = new G4DataVector;
			G4double tmpa = 0.;
            for(size_t jj=0;jj<numangles;jj++)
            {
				if (jj>8) file >> a;
				else {
					if (0==jj || 3==jj || 5==jj || 7==jj) file >> tmpa;// Read a new value from file...
					if (0==jj || 4==jj || 6==jj || 8==jj) a = tmpa;
					else if (3==jj || 5==jj || 7==jj)  a = (a + tmpa)/2.;
				}
               MottVector->push_back(a);
            }
            MottTable->push_back(MottVector);
            file >> a;
            file >> a;// Skip the cumulative cross section and mean free path for each energy
            // Instead, we recalculate the total cross section/inverse mean free path from the differential function (see below),
            // to avoid rounding errors
         }
      } else {// The Mott file does not exist
         // It seems that we have a problem, since the Mott files should be present for all elements up to number 94.
         G4cerr << "CADPhysicsSingleScattering: Mott file for element with Z=" << integralZ << " not found. Exiting now." << G4endl;
         exit(1);
         // Alternatively, one could just use zero cross sections but presumably this is not what the user really wants...
         // - see also DInvMFPTableforSemi above.
      }
      file.close();

      // Now start calculating the actual inverse mean free paths.
      G4double phononInvMFP = 1./numel * fermienergy/eV * resistivity / (4.67e-15*meter);// The phonon scattering InvMFP.
      // Note: since the total scattering InvMFP is a sum over all elements in the compound, for each individual
      // element a factor 1./numel is added to the phonon InvMFP.

      // First the part UP TO the Fermi energy. In principle, the particle should never get into this regime, but the data
      // are added anyway for the sake of completeness.
      // Note that in this range the InvMFP is assumed to be NOT a function of energy - similar to the low-energy limit of
      // eq (3a) in Fitting et al.
      G4int j=0;
      for (j=0; j<fermibin; j++) {
         G4DataVector* DInvMFPforEnergy = new G4DataVector;
         G4double lastvalue = 0;
         G4double cumSigmaThetaInvMFP = 0;
         for(size_t jj=0;jj<numangles;jj++)// Loop over all angles in the range
         {
            G4double cumDInvMFP = 0.5 * phononInvMFP * (1.-anglerange[jj]);// Isotropic angular distribution
            // (as in the low-energy limit of eq. (13a) in Schreiber et al.)
            DInvMFPforEnergy->push_back(cumDInvMFP);
            if (jj>0) {
               cumSigmaThetaInvMFP += (cumDInvMFP - lastvalue) * ( 1. - (anglerange[jj] + anglerange[jj-1])/2.);
            } else cumSigmaThetaInvMFP += cumDInvMFP * (1. - anglerange[jj]);
            lastvalue = cumDInvMFP;
         }
         DInvMFPforEl->push_back(DInvMFPforEnergy);
         SigmaThetaforEl->push_back(cumSigmaThetaInvMFP);
      }

      G4double density = theAtomNumDensityVector[ii];

      // Between the Fermi energy and 100 eV, we interpolate between the phonon and Mott values
      // N.B.: if the resistivity is not given (e.g. for gases), simply the full Mott values are used.

      // Find the correct Mott data between which to interpolate
      // Mottbin is the number of the Mott data energy bin, mottfac2 represents the relative position within the bin.
      for (j=fermibin; j<41; j++) {
         G4int Mottbin = 0;
         G4double mottfac2 = 0.;
         if (enrange[j]>20.*eV) {
            Mottbin = 1;
            mottfac2 = (enrange[j]-20.*eV)/(30.*eV);
            if (enrange[j]>50.*eV) {
               Mottbin=2;
               mottfac2 = (enrange[j]-50.*eV)/(25.*eV);
               if(enrange[j]>75.*eV) {
                  Mottbin=3;
                  mottfac2 = (enrange[j]-75.*eV)/(25.*eV);
                  if (enrange[j]>100.*eV) {
                     Mottbin=4;
                     mottfac2 = (enrange[j]-100.*eV)/(100.*eV);
                  }
               }
            }
         }
         // mottfac represents the relative weight of the Mott InvMFP in the total InvMFP.
         // This relative weight goes from zero at the Fermi energy (we just use the acoustic phonon InvMFP)
         // to unity at 100 eV (actually 100.01 eV, for numerical reasons - fermienergy is limited to 100.*eV)
         // Exception: if resistivity is equal to zero (for gases), mottfac equals unity - in other words: we just use the Mott values.
         G4double mottfac = (enrange[j]-fermienergy)/(100.01*eV-fermienergy);
         if (resistivity<=0.) mottfac=1.;

         G4DataVector* DInvMFPforEnergy = new G4DataVector;

         G4double preconst3 = enrange[j] + screening;
         G4double InvMFP = phononInvMFP * preconst3;
         G4double lastvalue = 0.;
         G4double cumSigmaThetaInvMFP = 0.;
         G4double cumInvMFP = 0.;
         G4double mottInvMFP = 0.;
         G4double phononCumInvMFP = 0.;
		 G4double prevmottvalue = 0.;
		 for(size_t jj=0;jj<numangles;jj++)// Loop over all angles
         {
            G4double mottvalue = 0.;
            if (Mottbin==0) mottvalue = (*(*MottTable)[0])[jj];
            else mottvalue = (1-mottfac2) *(*(*MottTable)[Mottbin-1])[jj] +
               (mottfac2)*(*(*MottTable)[Mottbin])[jj];// The interpolated Mott cumulative differential CS
            if(jj==0) {
               mottInvMFP = mottvalue * twopi * (1. - anglerange[jj]) * density * 1.e-14;
               // Converting from the CS per atom given in sq. Angstrom to an inverse MFP (in mm^-1)

               phononCumInvMFP = InvMFP / (enrange[j] + screening/(0.5-0.5*anglerange[jj]));// Phonon cumulative differential InvMFP
               // Note: the angular distribution follows the low-energy expression of eq. (13a) in Schreiber et al.

               cumInvMFP += mottfac*mottInvMFP + (1.-mottfac)*phononCumInvMFP;// Weighted interpolation between Mott and phonon parts
               cumSigmaThetaInvMFP += (mottfac*mottInvMFP + (1.-mottfac)*phononCumInvMFP) * (1.-anglerange[jj]);
            }
            else {
               mottInvMFP = 0.5*(mottvalue + prevmottvalue)* twopi * (anglerange[jj-1] - anglerange[jj]) * density * 1.e-14;
               // Converting from the CS per atom given in sq. Angstrom to an inverse MFP (in mm^-1)

               phononCumInvMFP = InvMFP / (enrange[j] + screening/(0.5-0.5*anglerange[jj]));// Phonon cumulative differential InvMFP

               cumInvMFP += mottfac*mottInvMFP + (1.-mottfac)*(phononCumInvMFP-lastvalue);// Weighted interpolation between Mott and phonon parts
               cumSigmaThetaInvMFP += (mottfac*mottInvMFP + (1.-mottfac)*(phononCumInvMFP - lastvalue))
                  * (1.-(anglerange[jj] + anglerange[jj-1])/2.);
            }
			prevmottvalue = mottvalue;
            lastvalue = phononCumInvMFP;
            DInvMFPforEnergy->push_back(cumInvMFP);
         }
         DInvMFPforEl->push_back(DInvMFPforEnergy);
         SigmaThetaforEl->push_back(cumSigmaThetaInvMFP);
      }

      for(j=41; j<63; j++)// From 200 eV up to 30 keV, just use the Mott cross section data
      {
         G4DataVector* DInvMFPforEnergy = new G4DataVector;
         G4double cumInvMFP = 0.;
         G4double InvMFP = 0.;
         G4double cumSigmaThetaInvMFP = 0.;
         for(size_t jj=0;jj<numangles;jj++)// Loop over all angles
         {
            if(jj==0) {
               InvMFP = (*(*MottTable)[j-37])[0] * twopi * (1. - anglerange[jj]) * density * 1.e-14;
               // Converting from the CS per atom given in sq. Angstrom to an inverse MFP (in mm^-1)
               cumInvMFP += InvMFP;// Build the cumulative angular-differential inverse MFP vector
               cumSigmaThetaInvMFP += InvMFP * (1.-anglerange[jj]);
            }
            else {
               InvMFP = 0.5*((*(*MottTable)[j-37])[jj-1] + (*(*MottTable)[j-37])[jj])* twopi * (anglerange[jj-1] - anglerange[jj]) * density * 1.e-14;
               // Converting from the CS per atom given in sq. Angstrom to an inverse MFP (in mm^-1)
               cumInvMFP += InvMFP;// Build the cumulative angular-differential inverse MFP vector
               cumSigmaThetaInvMFP += InvMFP * (1.-(anglerange[jj] + anglerange[jj-1])/2.);
            }
            DInvMFPforEnergy->push_back(cumInvMFP);
         }
         DInvMFPforEl->push_back(DInvMFPforEnergy);
         SigmaThetaforEl->push_back(cumSigmaThetaInvMFP);
      }

      for(j=63; j<78; j++)// For energies above 30 keV, we use the screened relativistic Rutherford cross-section
      {
            G4DataVector* DInvMFPforEnergy = new G4DataVector;
            G4double cumInvMFP = 0.;
            G4double InvMFP = 0.;
            G4double cumSigmaThetaInvMFP = 0.;
			G4double preva = 0.;
            for(size_t jj=0;jj<numangles;jj++)// ...and for all angles in the range.
            {
               a = ComputeRutherfordCrossSectionPerAtom(enrange[j],integralZ,acos(anglerange[jj])); // Compute using screened relativistic Rutherford crossection
               if(jj==0) {
                  InvMFP = a * twopi * (1. - anglerange[jj]) * density * 1.e-14;
                  // Converting from the CS per atom given in sq. Angstrom to an inverse MFP (in mm^-1)
                  cumInvMFP += InvMFP;// Build the cumulative angular-differential inverse MFP vector
                  cumSigmaThetaInvMFP += InvMFP * (1. - anglerange[jj]);
               }
               else {
                  InvMFP = (a + preva) /2. * twopi * (anglerange[jj-1] - anglerange[jj]) * density * 1.e-14;
                  // Converting from the CS per atom given in sq. Angstrom to an inverse MFP (in mm^-1)
                  cumInvMFP += InvMFP;// Build the cumulative angular-differential inverse MFP vector
                  cumSigmaThetaInvMFP += InvMFP * (1. - (anglerange[jj]+ anglerange[jj-1])/2.);
               }
               DInvMFPforEnergy->push_back(cumInvMFP);
			   preva = a;
            }
            DInvMFPforEl->push_back(DInvMFPforEnergy);
            SigmaThetaforEl->push_back(cumSigmaThetaInvMFP);
            if (j<63) {
               file >> a;
               file >> a;// Skip the cumulative cross section and mean free path for each energy
               // Instead, we recalculate the total cross section/inverse mean free path from the differential function (see below),
               // to avoid rounding errors
            }
        }
//      for(j=63; j<78; j++)// For energies above 30 keV, no Mott data are available, and hence we use the Browning approximation
//         // (even though this approximation is also not strictly valid above 30 keV!)
//      {
//         G4DataVector* DInvMFPforEnergy = new G4DataVector;
//         G4double energy = enrange[j];
//         G4double invmfp = ComputeBrowningCrossSectionPerAtom(energy, integralZ) * density;// Convert from cross section to Inv. MFP
//         DInvMFPforEnergy->push_back(invmfp);
//         DInvMFPforEl->push_back(DInvMFPforEnergy);
//      }
		  DInvMFPforMat->push_back(DInvMFPforEl);
      SigmaThetaforMat->push_back(SigmaThetaforEl);
   }
   // Calculate the total InvMFP for each energy (sum up the highest cumulative InvMFPs for all elements in the material)
   G4DataVector* InvMFPforMat = new G4DataVector;
   CADPhysicsDataTable* SumDInvMFPforMat = new CADPhysicsDataTable;

   if (verboseLevel>2) {
      G4cout << "Total mean free paths:" << G4endl;
      G4cout << "\tEnergy(eV)\tEMFP(nm)" << G4endl;
   }

   for(size_t j=0;j<78;j++)
   {
      G4double totinvmfp = 0;
      for(size_t ii=0;ii<numel;ii++)
      {
         CADPhysicsDataTable* DInvMFPforEl = (*DInvMFPforMat)[ii];
         G4DataVector* DInvMFPforEnergy = (*DInvMFPforEl)[j];
         size_t vectorsize = DInvMFPforEnergy->size();
         totinvmfp += (*DInvMFPforEnergy)[vectorsize-1];

      }
      if (verboseLevel>2) G4cout << mname << "\t" << enrange[j]/eV << "\t" << 1/(totinvmfp*nanometer) << G4endl;
      InvMFPforMat->push_back(totinvmfp);
   }

   // Now also calculate the total cumulative differential inverse MFP,
   for(size_t j=0;j<76;j++)
   {
      G4DataVector* SumDInvMFPforEnergy = new G4DataVector;
      for(size_t jj=0;jj<numangles;jj++)
      {
         G4double SumDInvMFP = 0.;
         for (size_t ii=0;ii<numel;ii++)
         {
            SumDInvMFP += (*(*(*DInvMFPforMat)[ii])[j])[jj];
         }
         SumDInvMFPforEnergy->push_back(SumDInvMFP);
      }
      SumDInvMFPforMat->push_back(SumDInvMFPforEnergy);
   }
   theSumDInvMFPTable->push_back(SumDInvMFPforMat);
   theInvMFPTable->push_back(InvMFPforMat);

   // Now calculate transport mean free path values and export to file...
   // These are the average distances the electron has to travel before 'isotropization' of its direction due to elastic scattering.
   // These values are used by the DI process to determine the effective distance electrons can travel inside the material
   // and hence decide if a certain electron has any chance to leave the current material, and therefore whether or not
   // it needs to be created or killed.
   // Also, the transport mean free path is used in the DoDiffusionStep method to determine the position and momentum distribution at
   // the end of a step.

   std::ostringstream ost;
   ost << "tmfp_" << mname << ".dat";
   G4String name = ost.str();

   G4DataVector* TMFPforMat = new G4DataVector;//new
   std::ofstream tmfpfile;
   tmfpfile.open (name.c_str());
   tmfpfile << setprecision(6);
   for(size_t j=0; j<78; j++) {
      G4double invtmfp = 0.;
      for(size_t ii=0;ii<numel; ii++) invtmfp += (*(*SigmaThetaforMat)[ii])[j];
      G4double tmfp = 1./invtmfp;
      tmfpfile << enrange[j] << "\t" << tmfp << G4endl;
      TMFPforMat->push_back(tmfp);
   }
   tmfpfile << "-1\t-1\n";
   tmfpfile.close();
   theTMFPTable->push_back(TMFPforMat);

   return DInvMFPforMat;
}

G4double CADPhysicsSingleScattering::ComputeRutherfordCrossSectionPerAtom(
   G4double kinEnergy,
   G4double Z,
   G4double angle)
    // Compute the screened relativistic differential Rutherford crossection
	// Eqs. (5.34), (5.37) in Reimer, Transmission Electron Microscopy, Fourth Edition
	// (note: (5.38) is valid only for small scattering angles)
	// but using the definition of R according to Bethe (i.e. including prefactor 0.885)
{

    G4double lambda = h_Planck*c_light/sqrt(2.0*electron_mass_c2*kinEnergy*(1.0+kinEnergy/(2.0*electron_mass_c2)));
    G4double a0 = h_Planck*h_Planck/(pi*electron_mass_c2*mu0*e_squared);

        // EB: Original:
	//G4double R = 0.885 * a0 * pow(Z,-0.333333333333333333);
        // EB: Corrected to match NIST c.s. data 2013-06-13
        //     It is not clear (yet) why this factor two is missing but I found that
        //     for Al,Si and Au the NIST data is matched at 200keV and 50keV.
        //     Also the BSE yield curve of Si shows a discontinuity around 30keV if this
        //     factor 2 correction is not done.
	G4double R = 2.0*0.885 * a0 * pow(Z,-0.333333333333333333);

	G4double fx = Z*lambda*lambda/(lambda*lambda+4.0*pi*pi*angle*angle*R*R);
    G4double factor = lambda*lambda*lambda*lambda/(64.0*pi*pi*pi*pi*a0*a0);
	G4double enfactor = 1.0 + kinEnergy/electron_mass_c2;

	// 2012-07-26, relativistic correction for theta0
	//G4double theta0 = lambda * pow(Z,0.333333333333333333)/(2.0*pi*a0);

    G4double stheta2 = sin(angle/2.0);

	G4double result = factor*(Z-fx)*(Z-fx)*enfactor*enfactor/(stheta2*stheta2*stheta2*stheta2);

    return (1.0e14*result); // Output in the same units as in the Mott files (barn = 1e-20 m^2 = 1e-14 mm^2)
}

G4double CADPhysicsSingleScattering::ComputeBrowningCrossSectionPerAtom(
   G4double kinEnergy,
   G4double Z)
   // Method used to calculate the scattering cross section for energies above 30 keV, where the
   // Mott cross sections are not available.
   // Ref.: R. Browning et al, J.Appl.Phys. 76 (4), 2016 (1994).
{
   G4double CrossSection = 0.;
   if (Z<0.9999) return CrossSection;
   if (kinEnergy<0.1*eV) return CrossSection;

   G4double En = (1.-pow(kinEnergy/electron_mass_c2+1.,-2.))*electron_mass_c2/2.;// 'Ad hoc' relativistic correction;
   // note however the limited validity of the Browning approximation - do not use at too high energies!
   En /= keV;
   G4double Zpow = pow(Z,1.7);
   G4double a1 = 0.00030;
   G4double a2 = 0.005;
   G4double a3 = 0.0007;
   G4double sqrtE = pow(En,0.5);

   CrossSection = a1*Zpow*nanometer*nanometer;
   CrossSection /= En + a2*Zpow*sqrtE + a3*Z*Z/sqrtE ;
   return CrossSection;
}

G4VParticleChange* CADPhysicsSingleScattering::AlongStepDoIt(
   const G4Track& track,const G4Step& /*step*/)//track,step
{
   // In the case of multi/diffusionstep, this is the place to propose energy loss, new position, and time spent via the ParticleChange
   aParticleChange.Initialize(track);
   if (diffusionstep || multistep) {
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

   // In the case of multi/diffusionstep, this is the place to propose new direction and final kinetic energy via the ParticleChange
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
   G4double cosTheta = -1.;

   if(energybin<42)// Below 200 eV, we use the phonon approximation to calculate the energy loss.
      // Hence, there is no need to determine on which individual element the electron was scattered
      // and we can use the summed differential inverse mean free paths per material instead of the individual
      // differential inverse mean free paths per element - this makes the simulation a little bit faster.
   {
      CADPhysicsDataTable* SumDInvMFPforMat = (*theSumDInvMFPTable)[i];
      if (0==energybin) {
         PSDInvMFPforEnergy1 = (*SumDInvMFPforMat)[0];
      } else PSDInvMFPforEnergy1 = (*SumDInvMFPforMat)[energybin-1];
      PSDInvMFPforEnergy2 = (*SumDInvMFPforMat)[energybin];
      G4double enfrac = 0;
      if (energybin>0) enfrac = (previousenergy-enrange[energybin-1])/(enrange[energybin]-enrange[energybin-1]);

      // Determine cosTheta from the *interpolated* cumulative distribution function
      // First a random number (xforangle) is drawn and then we scan through the tabulated data to find the matching scattering angle.
      // Because of the large number of angles (96), this 'scanning' is done in two steps - a coarse one and a finer one.
      G4double xx,oldxx=0;
      G4double xforangle = G4UniformRand()*preStepInvMFP;
      size_t jj=8;
      do
      {
         xx = (1.-enfrac)*(*PSDInvMFPforEnergy1)[jj] + enfrac* (*PSDInvMFPforEnergy2)[jj];
         if (xx>xforangle) break;
         jj+=9;
      } while (jj < numangles);
      jj-=8;
      if (jj>0) oldxx = (1.-enfrac)*(*PSDInvMFPforEnergy1)[jj-1] + enfrac* (*PSDInvMFPforEnergy2)[jj-1];
      do
      {
         xx = (1.-enfrac)*(*PSDInvMFPforEnergy1)[jj] + enfrac* (*PSDInvMFPforEnergy2)[jj];
         if (xx>xforangle)
         {
            // The correct angle 'bin' has been found; interpolate between stored values to find cosTheta
            if(jj>0) cosTheta = anglerange[jj-1] + (anglerange[jj]-anglerange[jj-1])*(xforangle-oldxx)/
               (xx-oldxx);
            else cosTheta = 1. + (anglerange[jj]-1.)*xforangle/xx;
            break;
         }
         jj++;
         oldxx=xx;
      } while (jj < numangles);

      // Energy loss in the scatter event is taken as the fixed average phonon loss value for this material.
      energyloss = fphononloss;
      finalT -= energyloss;
      if (finalT>0.) aParticleChange.ProposeEnergy(finalT);
      aParticleChange.ProposeLocalEnergyDeposit(energyloss);
   } else {// Above 200 eV, the procedure is different since we also need to select the right element (in a compound)
      G4double logenfrac = 0.;
      CADPhysicsDataCube* DInvMFPforMat = (*theDInvMFPTable)[i];

      // Select the energy value for selecting the element and scattering angle. Effectively, this gives a (log-log) interpolation of
      // all scattering probabilities between subsequent tabulated energy values.
      G4int binforcs = 0;
      if(78==energybin) binforcs = 77; else {
         logenfrac = (log(previousenergy)-log(enrange[energybin-1]))/(log(enrange[energybin])-log(enrange[energybin-1]));
         if (G4UniformRand()<logenfrac) binforcs = energybin; else binforcs = energybin - 1;
      }

      // Select the element
      G4double x = G4UniformRand()*(*InvMFPforMat)[binforcs];
      G4double xperel=0.;
      do {
         PSDInvMFPforEl = (*DInvMFPforMat)[ii];
         PSDInvMFPforEnergy1 = (*PSDInvMFPforEl)[binforcs];
         xperel = (*PSDInvMFPforEnergy1)[numangles-1];// Inverse MFP for the current element. EK 2013-01-21, corrected -1 in bin number
         x -= xperel;
         ii++;
      } while (x > 0.  && ii < nelm);
      ii--;

      // Select the scattering angle and calculate the corresponding energy loss
      const G4ElementVector* theElementVector = currentMaterial->GetElementVector();
      const G4Element* elm = (*theElementVector)[ii];
        {
         x += xperel;// Reuse the same random number x; x is now uniformly randomly distributed between 0 and xperel
         // Scan the angular differential inverse MFP in two parts (same as above)
         G4double xx,oldxx=0.;
         size_t jj=8;
         do
         {
            xx = (*PSDInvMFPforEnergy1)[jj];
            if (xx>x) break;
            jj+=9;
         } while (jj < numangles);
         jj-=8;
         if (jj>0) oldxx = (*PSDInvMFPforEnergy1)[jj-1];
         do
         {
            xx = (*PSDInvMFPforEnergy1)[jj];
            if (xx>x)
            {
               // The correct angle 'bin' has been found; interpolate between stored values to find cosTheta
               if(jj>0) cosTheta = anglerange[jj-1] + (anglerange[jj]-anglerange[jj-1])*(x-oldxx)/
                  (xx-oldxx);
               else cosTheta = 1 + (anglerange[jj]-1)*x/xx;
               break;
            }
            jj++;
            oldxx=xx;
         } while (jj < numangles);
      }

      // Calculate the "elastic" energy loss due to atom recoil
      // N.B.: this energy loss is extremely small compared to the kinetic energy - the place where you are most likely to notice it
      // is in the elastic backscatter peak, for which no other loss processes are present.
      G4double effA  = elm->GetA();
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
   CADPhysicsDataCube* DInvMFPforMat = (*theDInvMFPTable)[i];

   G4double costheta1,sintheta1,cosphi1,sinphi1;
   G4double tempx,tempy,tempz;
   G4ThreeVector relativeglobalposition = G4ThreeVector(0.,0.,0.);

   if (energybin<42) {// Below 200 eV, we use the phonon approximation to calculate the energy loss.
      // Hence, there is no need to determine on which individual element the electron was scattered
      // and we can use the summed differential inverse mean free paths per material instead of the individual
      // differential inverse mean free paths per element - this makes the simulation a little bit faster.
      energyloss = 0.;// Initialize the total energy loss to zero

      // Find the tabulated data for the current energy - it is assumed that throughout one call of the DoMultiStep method,
      // the energy of the particle changes so little that we can keep using the same inverse MFP data.
      CADPhysicsDataTable* SumDInvMFPforMat = (*theSumDInvMFPTable)[i];
      if (0==energybin) {
         PSDInvMFPforEnergy1 = (*SumDInvMFPforMat)[0];
      } else PSDInvMFPforEnergy1 = (*SumDInvMFPforMat)[energybin-1];
      PSDInvMFPforEnergy2 = (*SumDInvMFPforMat)[energybin];
      G4double enfrac = 0.;
      if (energybin>0) enfrac = (previousenergy-enrange[energybin-1])/(enrange[energybin]-enrange[energybin-1]);

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
         G4double xforangle = G4UniformRand()*preStepInvMFP;// Draw a random number
         G4double xx,oldxx=0.;
         size_t jj=8;
		 size_t kkk=90;
         // Scan the angular differential inverse MFP in two parts
         do
         {
            xx = (1.-enfrac)*(*PSDInvMFPforEnergy1)[jj] + enfrac* (*PSDInvMFPforEnergy2)[jj];
            if (xx>xforangle)
            {
               kkk=jj-8;
               jj=numangles;
            }
            jj+=9;
         } while (jj < numangles);
         if (kkk>0) oldxx = (1.-enfrac)*(*PSDInvMFPforEnergy1)[kkk-1] + enfrac* (*PSDInvMFPforEnergy2)[kkk-1];
         do
         {
            xx = (1.-enfrac)*(*PSDInvMFPforEnergy1)[kkk] + enfrac* (*PSDInvMFPforEnergy2)[kkk];
            if (xx>xforangle)
            {
               // The correct angle 'bin' has been found; interpolate between stored values to find cosTheta
               if(kkk>0) cosTheta = anglerange[kkk-1] + (anglerange[kkk]-anglerange[kkk-1])*(xforangle-oldxx)/
                  (xx-oldxx);
               else cosTheta = 1. + (anglerange[kkk]-1)*xforangle/xx;
               kkk=numangles;
            }
            kkk++;
            oldxx=xx;
         } while (kkk < numangles);

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
      G4double logenfrac = 0.;
      if (energybin<78) logenfrac = (log(previousenergy)-log(enrange[energybin-1]))/(log(enrange[energybin])-log(enrange[energybin-1]));
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
         // Select the energy value for selecting the element and scattering angle. Effectively, this gives a (log-log) interpolation of
         // all scattering probabilities between subsequent tabulated energy values.
         // Since there is a random component in this procedure, this needs to be redone for every individual step!
         G4int binforcs = 0;
         if(78==energybin) binforcs = 77; else {
            if (G4UniformRand()<logenfrac) binforcs = energybin; else binforcs = energybin - 1;
         }

         G4double x = G4UniformRand()*(*InvMFPforMat)[binforcs];
         G4double xperel=0.;
         size_t ii=0;
         // Select the element
         do {
            PSDInvMFPforEl = (*DInvMFPforMat)[ii];
            PSDInvMFPforEnergy1 = (*PSDInvMFPforEl)[binforcs];
            xperel = (*PSDInvMFPforEnergy1)[95];
            x -= xperel;
            ii++;
         } while (x > 0. && ii < nelm);
         ii--;

         const G4ElementVector* theElementVector = currentMaterial->GetElementVector();
         const G4Element* elm = (*theElementVector)[ii];
         G4double effA  = elm->GetA();
          {
            x += xperel;// Reuse the same random number x; x is now uniformly randomly distributed between 0 and xperel
            // Scan the angular differential inverse MFP in two parts
            G4double xx,oldxx=0.;
            size_t jj=8;
            do
            {
               xx = (*PSDInvMFPforEnergy1)[jj];
               if (xx>x) break;
               jj+=9;
            } while (jj < numangles);
            jj-=8;
            if (jj>0) oldxx = (*PSDInvMFPforEnergy1)[jj-1];
            do
            {
               xx = (*PSDInvMFPforEnergy1)[jj];
               if (xx>x)
               {
                  // The correct angle 'bin' has been found; interpolate between stored values to find cosTheta
                  if(jj>0) cosTheta = anglerange[jj-1] + (anglerange[jj]-anglerange[jj-1])*(x-oldxx)/(xx-oldxx);
                  else cosTheta = 1. + (anglerange[jj]-1.)*x/xx;
                  break;
               }
               jj++;
               oldxx=xx;
            } while (jj < numangles);
         }
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

   }//end if(energybin)
   // Update position, energy and direction for later registration via ParticleChange.
   proposeposition = initialposition + relativeglobalposition;
   if(energyloss<finalT) finalT -= energyloss;// Extra check: energy loss cannot be larger than total available energy
   elDirectionnew = momdir;
   return steptaken;
}

G4double CADPhysicsSingleScattering::DoDiffusionStep(
   const G4Track& track,
   G4double othersteplength,
   G4double safety)
{
   // Initialization
   G4double steplength = othersteplength;
   G4ThreeVector initialposition = track.GetPosition();
   G4ThreeVector momdir = track.GetMomentumDirection();

   // The approach is as follows: first we move the particle a distance tmfp in the initial direction to an
   // intermediate position. This serves as an offset for the final position, it is the distance the particle travels on average
   // before the initial direction is 'forgotten'.
   // According to the remaining step length, the position is drawn from a Gaussian distribution following a 'diffusion' model.

   // First step
   G4ThreeVector relativeglobalposition = momdir*tmfp;
   steplength -= tmfp;

   // Probe the position
   G4double r = 0.;
   G4ThreeVector rdp;
   G4double sigma = sqrt(steplength * tmfp * 2./3.) ;// Note, the first step already subtracted from steplength
   do {
      rdp = sigma * G4ThreeVector(G4RandGauss::shoot(),G4RandGauss::shoot(),G4RandGauss::shoot());// Gaussian-distributed 3D position
      r = rdp.mag();// Linear distance from the intermediate position
   } while (r+tmfp>safety || r>steplength);// Accept this draw if the total distance from the initial position in the worst case is less
   // than the current safety, AND if the linear distance from the intermediate to the final position is less than the allowed steplength.
   // The latter check is necessary because the Gaussian distributions (unphysically) extend all the way to infinity.
   relativeglobalposition += rdp;// Add the 'diffusion' part to the total particle displacement

   // Probe the direction
   // The final direction of the particle in the diffusion model turns out to be position dependent. It gets more anisotropic towards
   // the outside of the distribution.
   // The approach followed here is that the direction probability is proportional to the position probability at one tmfp *before* the total
   // travelled path length.
   // If r represents the final position, u is the current direction and theta is the angle between r and u, then
   // define r'=r-u*tmfp*cos(theta), which is an approximation for the particle's position at the last scatter event. Now
   // p(cos(theta))~exp(-(r')^2/(2*sigma^2)) with sigma = sqrt((steplength - tmfp)*tmfp*2/3).
   // After renormalization of the probability distribution and a small simplification (using tmfp<<r) we get the following expressions:
   G4bool accept = false;
   G4double x,y,z;
   G4double fac1 = (3./2.)*r/(steplength-tmfp);
   G4double fac2 = exp(-fac1);
   do {
      // First draw u from a uniform distribution, then apply an accept/reject method to adjust it to the abovementioned distribution.
      x = -1.+2.*G4UniformRand();
      G4double phi = twopi*G4UniformRand();
      y = sqrt(1.-x*x)*cos(phi);
      z = sqrt(1.-x*x)*sin(phi);
      G4double costheta = (x*rdp.x() + y*rdp.y() + z*rdp.z())/r;
      if (G4UniformRand()<exp(fac1*costheta)*fac2) accept = true;
   } while (!accept);
   elDirectionnew = G4ThreeVector (x,y,z);
   proposeposition = initialposition + relativeglobalposition;

   // N.B.: there is no need to implement an energy loss here; it's already done
   // in AlongStepGetPhysicalInteractionLength
   return othersteplength;
}


G4double CADPhysicsSingleScattering::GetMeanFreePath(
   const G4Track& track,
   G4double,
   G4ForceCondition*)
{
   G4double preStepMFP;
   G4double preStepKinEnergy = track.GetKineticEnergy();
   const G4MaterialCutsCouple* couple = track.GetMaterialCutsCouple();
   if(couple==currentCouple) {
      // Return the old MFP value if the kinetic energy at which it was calculated was close enough to the current value, i.e. within 1 eV
      if (abs(preStepKinEnergy-previousenergy)<1.*eV) return previousMFP;
   } else DefineMaterial(couple);

   // Check the physical state of the material. The switch 'workforgases' determines if elastic scattering is done at all for gasesous materials
   // (or if this is left to more sophisticated gas scattering models)
   if (currentMaterial->GetState()==kStateGas && !workforgases)
   {
      preStepMFP = DBL_MAX;
      previousMFP = preStepMFP;
      previousenergy = preStepKinEnergy;
      return preStepMFP;
   }

   G4double InvMFPforMFP = GetInvMFPFromTable(preStepKinEnergy);// Get the inverse MFP from tabulated data.
   if(0.<InvMFPforMFP) preStepMFP = 1./InvMFPforMFP;// Convert the inverse MFP to a mean free path
   else preStepMFP = DBL_MAX;
   previousMFP = preStepMFP;
   previousenergy = preStepKinEnergy;
   return preStepMFP;
}

void CADPhysicsSingleScattering::DefineMaterial(const G4MaterialCutsCouple* couple)// Update pointers for the current material;
// These are used in GetInvMFPFromTable, AlongStepGetPhysicalInteractionLength, and the various Do methods.
{
   currentCouple   = couple;
   currentMaterial = couple->GetMaterial();
   currentMaterialIndex = couple->GetIndex();
   InvMFPforMat = (*theInvMFPTable)[currentMaterialIndex];
   fphononloss = vec_phononloss[currentMaterialIndex];
   fTMFPforMat = (*theTMFPTable)[currentMaterialIndex];
}

G4double CADPhysicsSingleScattering::GetInvMFPFromTable(G4double e)
{
   // This method calculates *TWO* inverse MFP values. 'invmfp' is a log-log interpolation between tabulated data, and it is used to calculate
   // the total mean free path.
   // 'preStepInvMFP', on the other hand, is a linear interpolation of the same number. It is needed to normalize the cumulative distributions
   // per element and per angle (below 200 eV), since these are also defined as linear interpolations between tabulated data.
   G4double invmfp = 0.;// Initialization
   tmfp = 0.;
   if(e<enrange[0])
   {
      // If the kinetic energy is below the lowest tabulated energy (0.01 eV), just use the inverse MFP value at that lowest energy
      invmfp = (*InvMFPforMat)[0];
      preStepInvMFP = invmfp;
      tmfp = (*fTMFPforMat)[0];
      energybin = 0;// We're below the lowest tabulated energy
   }
   else
   {
      size_t j = 0;
      do {
         if(enrange[j+1]>e) {
            invmfp=(*InvMFPforMat)[j]*exp((log((*InvMFPforMat)[j+1])-log((*InvMFPforMat)[j]))*(log(e)-log(enrange[j]))/(log(enrange[j+1])-log(enrange[j])));
            if (j<41) preStepInvMFP = (*InvMFPforMat)[j]+(e-enrange[j])/(enrange[j+1]-enrange[j])
               *((*InvMFPforMat)[j+1]-(*InvMFPforMat)[j]);// Linear interpolation of the inverse mean free path for energies below 200 eV
            if(j<78) tmfp=(*fTMFPforMat)[j]+(e-enrange[j])/(enrange[j+1]-enrange[j])
               *((*fTMFPforMat)[j+1]-(*fTMFPforMat)[j]);// Calculate the transport mean free path (for use in DoDiffusionStep)
            energybin = j+1;
            break;
         }
         j++;
      } while (j+1<78);
      if(j==77)
      {
         invmfp = (*InvMFPforMat)[77];// For energies above the highest tabulated value, just take the value at that highest energy
         energybin = 78;
      }
   }
   return invmfp;
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
   diffusionstep = false;
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

   if (noraytrace) {// Do no 'ray tracing'
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
      if (mysafety>10.*tmfp && dodiffusionstep && tmfp>0.) {// Both currentSafety and PhysicalStep are large enough (i.e. more than
         // ten times the transport MFP) to do diffusionstep
         diffusionstep = true;

         // The final steplength is the lesser of PhysicalStep and an upper bound determined by the diffusion model
         value = min(PhysicalStep,currentSafety*currentSafety/(10.*tmfp));
         energyloss = 0.;
         finalT = track.GetKineticEnergy();
         // For energies below 200 eV, a limit on the total energy loss is imposed
         if (finalT<200.*eV) {
            energyloss = fphononloss * value / currentInteractionLength;
            if (energyloss>0.5*eV) {
               energyloss = 0.5*eV;
               value = currentInteractionLength * energyloss / fphononloss;// The total step length is scaled down linearly with the energy loss
               if (value < 10.*tmfp) diffusionstep = false;// ...and we need to check if it is still allowed to use the diffusion model
            }
            if (energyloss<finalT) finalT -= energyloss;
         }
         if (diffusionstep) return DoDiffusionStep(track,value,currentSafety);
      }
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
         // We're at the beginning of tracking, or just after a normal single/multi/diffusionstep of this process
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
