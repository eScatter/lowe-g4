// CADPhysicsDI.cc
// See CADPhysicsDI.hh for a general description of the CADPhysicsDI class.

// History:
// 2006-07-10: changed folder for id-files to current working directory
//             and removed several redundant warning messages in PostStepDoIt
// 2006-07-27: fixed a bug in ReadEpsilonFile. If the df file does not exist, we still need to put
//             entries into innershells and innershelltable; otherwise we get the material couple
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
// 2007-07-18: Fixed bug in calculation of innershelldifflambdavalue in CalculateCDCS
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

#include "CADPhysicsDI.hh"
#include "CADPhysicsDIMessenger.hh"
#include "Randomize.hh"
#include "G4ProductionCutsTable.hh"
#include <sstream>
#include "G4MaterialPropertiesTable.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Gamma.hh" //for X-ray fluorescence
#include "G4Material.hh"
#include "G4CrossSectionHandler.hh"

#ifndef CADPHYSICS_RELATIVISTIC
#define CADPHYSICS_RELATIVISTIC 1
#endif

using namespace std;

CADPhysicsDI::CADPhysicsDI(const G4String& processName)
: G4VContinuousDiscreteProcess(processName),
LowestKineticEnergy(0.01*eV),
HighestKineticEnergy(1.*MeV),
TotBin(1001),
L(0),
LPrim(0),
imeps(0),
rangetable(0),
innershelltable(0),
lambdatable(0),
difflambdatable(0),
fconductortype(0),
generateSecondaries(true),
generateXrays(false),
generateAugers(true),
rangecut(true),
energycut(false),
cutenergy(0.),
pairsgenerated(0),
crossSectionHandler(0)
{
   SetVerboseLevel(1);

   // Create the messenger
   messenger = new CADPhysicsDIMessenger(this);

   // Tell the deexcitation manager to also generate Auger electrons.
   // Can be overruled on the fly (by calling SetGenerateAugers(false))
   deexcitation.ActivateAugerElectronProduction(true);

   killthisone = false;
}

CADPhysicsDI::~CADPhysicsDI()
{
   // Free up memory in the data structures (perhaps redundant, since this method is only called
   // when the program exits anyway)
   if(difflambdatable) {
      CADPhysicsDataTable* a=0;
      G4DataVector* b=0;
      while (difflambdatable->size()>0) {
         a = difflambdatable->back();
         difflambdatable->pop_back();
         while (a->size()>0) {
            b = a->back();
            a->pop_back();
            if ( b ) {delete b;}
         }
         delete a;
      }
      delete difflambdatable;
   }
   if (lambdatable) {
      G4DataVector* a=0;
      while (lambdatable->size()>0) {
         a = lambdatable->back();
         lambdatable->pop_back();
         if ( a ) {delete a;}
      }
      delete lambdatable;
   }
   if (imeps) {delete imeps;}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..........oooOO0OOooo....

void CADPhysicsDI::BuildPhysicsTable(const G4ParticleDefinition& /*aParticleType*/)
{
   // Does the initialization of the process: fill the required data tables.

   if(verboseLevel) G4cout << "CADPhysicsDI: Start BuildPhysicsTable" << G4endl;

   // First clean up any existing data structures. This is necessary because this method may be
   // called more than once, in case the geometry has been modified between simulations.
   if( difflambdatable!=0 ) delete difflambdatable;
   difflambdatable = new std::vector<CADPhysicsDataTable*>;
   if( lambdatable!=0 ) delete lambdatable;
   lambdatable = new CADPhysicsDataTable;
   if( innershelltable!=0 ) delete innershelltable;
   innershelltable = new CADPhysicsDataTable;
   innershells.clear();
   vec_fermieff.clear();
   vec_minimumlimit.clear();
   vec_bandgap.clear();
   vec_conductortype.clear();
   vec_barrier.clear();
   if( rangetable != 0 ) delete rangetable;
   rangetable = new CADPhysicsDataTable;

   // Create material-independent tabulated data
   BuildEnrange();
   TabulateL();

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
      G4String mformula = material->GetChemicalFormula();
      G4bool isgas = (material->GetState()==kStateGas);
      if(verboseLevel) G4cout << "CADPhysicsDI: Start Material " << mname << "\n";

      // Get material properties from the GDML data.
      //
      // An overview:
      // Work function - energy between Fermi level and vacuum
      // Fermi energy  - kinetic energy of the electron at the Fermi level (often zero for
      //                 semiconductors and insulators)
      // Affinity      - energy between bottom of the conduction band and vacuum. Used instead of
      //                 work function for insulators.
      // Bandgap       - energy between bottom of the conduction band and top of valence band (for
      //                 semiconductors and insulators). Considered as a first 'inner shell' in
      //                 this model.
      // Bandbending   - bending of energy bands near the (vacuum) interface; has an effect on the
      //                 potential barrier at vacuum.
      // Deltaphi      - additional potential of the whole sample (e.g. due to conductive contact
      //                 with another material). In the current process, reading Deltaphi has no
      //                 other purpose than to print it to the screen.
      // Barrier       - total potential barrier towards vacuum. It is not read from file but
      //                 calculated here.
      G4int cond = 0;
      G4double work = 0.;
      G4double fermieff = 0.;
      G4double affinity = 0.;
      G4double bb = 0.;
      G4double deltaphi = 0.;
      G4double bandgap = 0.;
      G4double barrier = 0.;
      G4MaterialPropertiesTable* aMPT = material->GetMaterialPropertiesTable();
      if (aMPT) {
         if(aMPT->ConstPropertyExists("CONDUCTORTYPE")) {
            cond = (int)(aMPT->GetConstProperty("CONDUCTORTYPE"));
            switch (cond) {
               case 0:
                  if(verboseLevel>1) G4cout << "Material type: metal\n";
                  break;
               case 1:
                  if(verboseLevel>1) G4cout << "Material type: semiconductor\n";
                  break;
               case 2:
                  if(verboseLevel>1) G4cout << "Material type: insulator\n";
                  break;
               default:
                  G4cout << "CADPhysicsDI: WARNING: invalid material type!\n";
                  break;
            }
         }
         if(aMPT->ConstPropertyExists("WORKFUNCTION"))
            work = aMPT->GetConstProperty("WORKFUNCTION");
         if(aMPT->ConstPropertyExists("FERMIENERGY"))
            fermieff = aMPT->GetConstProperty("FERMIENERGY");
         if(aMPT->ConstPropertyExists("AFFINITY")) affinity = aMPT->GetConstProperty("AFFINITY");
         if(aMPT->ConstPropertyExists("BANDGAP")) bandgap = aMPT->GetConstProperty("BANDGAP");
         if(aMPT->ConstPropertyExists("BANDBENDING")) bb = aMPT->GetConstProperty("BANDBENDING");
         if(aMPT->ConstPropertyExists("DELTAPHI")) deltaphi = aMPT->GetConstProperty("DELTAPHI");
      }
      if(verboseLevel>1) {
         if (work) G4cout << "Work function = " << work << G4endl;
         if (fermieff) G4cout << "Fermi energy  = " << fermieff << G4endl;
         if (affinity) G4cout << "Affinity      = " << affinity << G4endl;
         if (bb) G4cout << "Band bending  = " << bb << G4endl;
         if (deltaphi) G4cout << "Delta phi     = " << deltaphi << G4endl;
      }

      G4double minlimit = enrange[FindIntLogen(fermieff)+1];

      if (cond==2) {// Insulator

         // If the work function is not defined separately, then use 'affinity' instead.
         if(!work) work = affinity;

         vec_fermieff.push_back(fermieff);
         vec_minimumlimit.push_back(minlimit);
         vec_bandgap.push_back(bandgap);
         vec_conductortype.push_back(2);
         barrier = fermieff + work - bb;
         vec_barrier.push_back(barrier);
      } else if (cond==1) {// Semiconductor

         // If the work function is not defined separately, then use 'affinity' instead.
         if(!work) work = affinity;

         vec_fermieff.push_back(fermieff);
         vec_minimumlimit.push_back(minlimit);
         vec_bandgap.push_back(bandgap);
         vec_conductortype.push_back(1);
         barrier = fermieff + work - bb;
         vec_barrier.push_back(barrier);
      } else {// Metal
         vec_fermieff.push_back(fermieff);
         vec_minimumlimit.push_back(minlimit);
         vec_bandgap.push_back(0);
         vec_conductortype.push_back(0);
         barrier = fermieff + work - bb;
         vec_barrier.push_back(barrier);
      }

      if(verboseLevel) G4cout << "CADPhysicsDI: Calculate cross sections for " << mname << "\n";
      // Read inner shell energies and energy loss function from file
      ReadEpsilonFile(mname,mformula,isgas,barrier);
      // Interpolate to 'internal' energy values as defined in the enrange vector
      InterpolateEpsilon();
      CADPhysicsDataTable* difflambdaformat = new CADPhysicsDataTable;
      G4DataVector* lambda = new G4DataVector;
      // Calculate the (differential) cross sections and ranges
      CalculateCDCS(difflambdaformat, lambda, fermieff, barrier);
      difflambdatable->push_back(difflambdaformat);
      lambdatable->push_back(lambda);

      if(mname != "Galactic") {
         std::ostringstream ost;
         ost << "imfp_" << mname << ".dat";
         G4String name = ost.str();
         std::ofstream imfpfile;
         imfpfile.open (name.c_str());
         imfpfile << setprecision(6);
         for(i=0; i<enrange.size(); i++) {
            if (i==0) {
               imfpfile << setw(12) << "Energy" << "\t" << setw(12) << "InelasticMFP" << G4endl;
               imfpfile << setw(12) << "(eV)"   << "\t" << setw(12) << "(nm)"       << G4endl;
            }
            if ((*lambda)[i] != 0.0) {
               imfpfile << scientific << setw(12) << enrange[i]/eV << "\t" << setw(12)
                      << (1.0/(*lambda)[i])/nanometer << G4endl;
            }
         }
         imfpfile.close();
      }
   }

   if(verboseLevel>1) PrintInfoDefinition();
}

void CADPhysicsDI::BuildEnrange()
{
   // Fill the vector enrange with TotBin+1 logarithmically spaced energy values from
   // LowestKineticEnergy to HighestKineticEnergy but double the number of entries between 1 and
   // 100 eV, since this is where most of the interesting stuff happens!
   if(verboseLevel) G4cout << "CADPhysicsDI: Start BuildEnrange\n";
   enrange.clear();
   G4double logstep=pow(10.,(log(HighestKineticEnergy/LowestKineticEnergy)/log(10.)+2)/(TotBin-1));
   G4double log2step = pow(logstep,0.5);
   G4int i;
   G4double energy = LowestKineticEnergy;
   for (i=0;i<TotBin;i++) {
      enrange.push_back(energy);
      if(0.999*eV<energy && 100.*eV>energy) {energy*=log2step;} else {energy*=logstep;}
   }

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

void CADPhysicsDI::TabulateL()
{
   if(verboseLevel>1) G4cout << "CADPhysicsDI: Start TabulateL\n";
   L.clear();
   LPrim.clear();
   G4double a = 1.e-8;
   G4double Lvalue,LPrimvalue;
   G4double log2step = pow(10.,(log(HighestKineticEnergy/LowestKineticEnergy)/log(10.)+2.)/(2.*(TotBin-1.)));
   while (a<0.5) {// Loop over a, defined as omegaprime/Ekin

      // Calculate the value of L, as used by Ashley in eqs. 17 and 19
      // Our implementation differs somewhat from Ashley's.
      // To start with, we use different expressions for omegaprime below and above 50 eV, where
      // the lower omegaprime values roughly correspond to plasmon excitation and ionization from
      // the valence band, while the higher values correspond to inner shell ionization. The values
      // for omegaprime < 50 eV are stored in L, while those for higher energies are stored in
      // Lprim.
      // For L we use an expression that is not corrected for exchange; simply because
      // this appears to generate more accurate results in practice compared to the
      // exchange-corrected version.
      // Physically, this can be justified by considering that in most cases a secondary electron
      // is not created immediately, but instead indirectly through the decay of a plasmon. The
      // expression is slightly different from the one used by Ashley (eq. 20) for two reasons:
      // 1) For the higher limit of integration over omega we use
      //       omega+ = 1/2(Ekin - E_F + omegaprime)
      //    instead of just
      //                1/2(Ekin + omegaprime),
      //    where E_F is the Fermi energy. Note that since the latter is material dependent,
      //    this correction has to be done in CalculateCDCS instead of here, so an extra term is
      //    added to L there.
      // 2) Our evaluation of L simply results in a slightly different expression than Ashley's...

      G4double s = sqrt(1.-2.*a);
      Lvalue = std::log(-1. + (2./a)*(1.+s));
      // For Lprim we use a very simple expression -log(a). Unlike Ashley's expression this results
      // in a nonzero differential cross section also for omegaprime > Ekin/2 (or, equivalently,
      // a>0.5).
      // Momentum is no longer conserved if just the primary and secondary electron are taken into
      // account; however, this expression more closely reflects the 'real' shape of individual atom
      // inner-shell ionization cross sections, which are nonzero all the way down to Ekin = Ebind.
      // N.B.: apparently the excess momentum of the primary electron is taken up by the the atom.
      // N.B.2: this expression does not take into account the limitation omega<Ekin-E_F, but for
      //        omegaprime > 50.eV this extra limitation has only a very minor influence.
      LPrimvalue = -log(a);
      L.push_back(Lvalue);
      LPrim.push_back(LPrimvalue);
      a *= log2step;// Step a
   }
   while (a<1.) {// LPrim continues up to a equal unity, as explained above.
      LPrim.push_back(-log(a));
      a *= log2step;
   }
}

G4bool CADPhysicsDI::ExistsFileForMaterial(const G4String &name)
{
   G4bool fileExists = false;

   std::ifstream test(name);
   std::filebuf* lsdp = test.rdbuf();

   if (lsdp->is_open()) {
      // Return true if the file with name 'name' has been opened succesfully
      fileExists = true;
   }

   // ...but close it again first.
   test.close();
   return fileExists;
}

G4String CADPhysicsDI::GetdfFileForMaterial(const G4String &materialName)
{
   G4String name="";
   std::ostringstream ost;
   // Build the complete string identifying the 'df' file with the data set, including its
   // relative path

   // First look in the users directory, this means that the user's df file will overrule the
   // standard df file
   char *g4workdir=getenv("G4WORKDIR");
   if (g4workdir) {
      ost << g4workdir << "/df/df_" << materialName << ".dat";
      name = ost.str();
   }
   if ((name == "") || ! ExistsFileForMaterial(name) ) {
      // Now look in the standard location given by environment variable CADPHYSICS_BASE
      char *path=getenv("CADPHYSICS_BASE");
      ost.str("");
      ost.clear();
      if (!path) {
         ost << "../CADPhysics/df/df_" << materialName << ".dat";
      } else {
         ost << path << "/df/df_" << materialName << ".dat";
      }
      name = ost.str();
      // G4cerr << "Looking for file : (" << ost.str() << ")" << G4endl;

      if ( !ExistsFileForMaterial(name)) {
         name="";
      }
   }
   return name;
}

void CADPhysicsDI::ReadEpsilonFile(const G4String& materialName,
      const G4String& materialFormula,
      G4bool isgas,
      G4double barrier)
{
   // Initialization
   if(verboseLevel>1)
      G4cout << "CADPhysicsDI: Start ReadEpsilonFile for material " << materialName << "\n";
   G4bool fillzeros = false;
   G4bool anyinnershell = false;
   G4DataVector* innershellenergies = new G4DataVector;

   // Check for the existence of the material file, otherwise just assume zero cross section values.
   // First try by name...
   G4String fullFilenameForMaterial=GetdfFileForMaterial(materialName);
   if(fullFilenameForMaterial == "") {
      // ...then by formula...
      fullFilenameForMaterial=GetdfFileForMaterial(materialFormula);
      if(fullFilenameForMaterial == "") {
         // N.B.: the rationale behind this is that materials with the same composition can still
         // have different energy loss functions.
         // This is the case, e.g., for the different allotropies of carbon: Diamond, Graphite and
         // GlassyCarbon, that all have their own df files.
         // If both fail, this might be a gas (for which this process is not supposed to do
         // anything anyway)
         fillzeros = true;
         if (!isgas) {
            // But if it's not, we have a problem and it is not possible to continue.
            G4cout << "CADPhysicsDI: Error: input file for material " << materialName
                   << " not found, exiting now!" << G4endl;
            exit(1);
         }
      }
   }

   if (fillzeros) {
      // We found no data file, so just fill the corresponding data structures with zeros as
      // place holders.
      readepsenergies.push_back(1.);
      readepsdata.push_back(0.);
      innershelltable->push_back(innershellenergies);
      innershells.push_back(false);
   } else {
      // Read the file for specific material and store to array readimeps of E, Im(-1/eps)

      // Opening file for input
      if (verboseLevel>0)
         G4cout << "CADPhysicsDI: Using df file: " << fullFilenameForMaterial << G4endl;
      std::ifstream file(fullFilenameForMaterial);
      G4double a = 0.;
      do {
         file >> a;
         if (a!=-1 && a<100.) {
            innershellenergies->push_back(a*eV);
            anyinnershell = true;
         }
         // Each file should start with a list (could be empty)
         // of inner shell energies ordered from small to large,
         // followed by -1
      } while (a!=-1);
      innershelltable->push_back(innershellenergies);
      innershells.push_back(anyinnershell);
      G4int k = 1;
      readepsenergies.clear();
      readepsdata.clear();
      do {
         file >> a;
         G4int nColumns = 2;
         // The file is organized into two columns:
         // 1st column is the energy
         // 2nd column is the corresponding value
         // The file terminates with the pattern: -1   -1

         if (a == -1) {
         } else {
            if (k%nColumns != 0) {
               G4double e = a * eV;
               readepsenergies.push_back(e);
               k++;
            } else if (k%nColumns == 0) {
               G4double value = a;
               readepsdata.push_back(value);
               k = 1;
            }//if (k%nColumns !=0)
         }//if (a == -1)
      } while (a != -1); // End of file
      file.close();
   }//if (fillzeros)

   // Now try to input the 'transport mean free path' data generated by the elastic scattering
   // process...   See CADPhysicsSingleScattering.cc for more info
   tmfpenergies.clear();
   tmfpdata.clear();
   std::ostringstream ost;
   ost << "tmfp_" << materialName << ".dat";
   G4String tmfpname = ost.str();
   if(ExistsFileForMaterial(tmfpname)){// Check if the corresponding 'tmfp' file exists
      size_t k=1;
      G4double a = 0.;
      std::ifstream file(tmfpname);
      G4double energy = 0.;
      G4double tmfpmax = 0.;
      do {
         file >> a;
         G4int nColumns = 2;

         // The file is organized into two columns:
         // 1st column is the energy
         // 2nd column is the corresponding value
         // The file terminates with the pattern: -1   -1

         if (a != -1) {
            if (k%nColumns != 0) {
               tmfpenergies.push_back(a);
               energy = a;
               k++;
            } else {
               // Make the function monotonously increasing for energies above the vacuum
               // potential barrier
               if (energy > barrier) {
                  if (a>tmfpmax) tmfpmax = a; else a = tmfpmax;
               }
               tmfpdata.push_back(a);
               k = 1;
            }//if (k%nColumns !=0)
         }//if (a != -1)
      } while (a != -1); // end of file
      file.close();
   } else {
      // If it doesn't exist, just use 'infinity' values
      G4cout << "Warning: file " << tmfpname << " doesn't exist.\n";
      tmfpenergies.push_back(1.);
      tmfpdata.push_back(DBL_MAX);
   }
}

void CADPhysicsDI::InterpolateEpsilon()
{
   // Interpolate read values of the energy loss function Im(-1/eps) to the energy values in
   // enrange and store to vector imeps
   if( imeps != 0 ) delete imeps;
   imeps = new G4DataVector;
   G4double e,tempimeps;
   G4int readepssize = readepsdata.size();
   G4int j=0;
   tempimeps = readepsdata[0];
   for (int i=0;i<TotBin;i++) {
      e = enrange[i];
      while(e > readepsenergies[j] && j<readepssize-1) {j++;}
      if (e > readepsenergies[j]) {
         // We have reached the end of tabulated data.
         // From here on we use an inverse fourth power extrapolation from the final data point.
         // This is a fairly good approximation for sufficiently high energies and as long as there
         //  are no new inner shell boundaries in the 'missing' part of the data set.
         tempimeps = readepsdata[j] * pow(readepsenergies[j]/e,4.);
      } else if (j>0) {
         // Log-log interpolation of imeps between file data. This matches a power-law behavior of
         // the underlying data (which is reasonable at least at high energies - see above)
         tempimeps =
            exp(log(readepsdata[j-1])+(log(e)-log(readepsenergies[j-1]))/
                (log(readepsenergies[j]) -log(readepsenergies[j-1]))*
                (log(readepsdata[j])-log(readepsdata[j-1])) );
      }
      imeps->push_back(tempimeps);
   }
}

void CADPhysicsDI::CalculateCDCS (CADPhysicsDataTable* difflambdaformat,
      G4DataVector* lambda,
      G4double fermiEnergy,
      G4double barrier)
{
   G4int i;// Counter over the energy values in enrange for the electron energy
   G4int ii = 0;
   for(i=0;i<TotBin;i++) {// Vary the electron energy
      G4double difflambdavalue = 0.;
      G4double innershelldifflambdavalue = 0.;
      G4DataVector* difflambda = new G4DataVector;
      G4int j=0;// Counter over the energy values in enrange for omega prime
      G4int jj=0;
#ifdef CADPHYSICS_RELATIVISTIC
      // Use an approximate relativistic correction (good up to ~1 MeV):
      G4double beta2mc2over2 = (1.-pow(enrange[i]/electron_mass_c2 + 1.,-2.))*electron_mass_c2/2.;
      G4double preconst = 1./(twopi * Bohr_radius * beta2mc2over2);
#else
      // Just the classical cross section
      G4double preconst = 1./(twopi * Bohr_radius * enrange[i]);
#endif
      G4double reduceden = enrange[i] - fermiEnergy;
      G4double omegaprime = enrange[0];
      while(omegaprime<reduceden){// Loop over omega prime values
         G4double dcs = 0.;
         // k is a measure for the ratio omega prime / Ekin, used for looking up the L and Lprim
         // values
         G4int k = jj - ii + (TotBin-1)*8/5;

         // Calculation of the actual differential cross section wrt omega prime
         // This is the integrand in eq. 19 of Ashley
         // We make a distinction between omegaprime smaller or larger than 50 eV, see TabulateL
         if(omegaprime>50.*eV) dcs = (*imeps)[j] * LPrim[k];
         else if(omegaprime<0.5*enrange[i]) {
            G4double factor = log(reduceden - omegaprime) - log(reduceden + omegaprime);
            dcs = (*imeps)[j] * (factor + L[k])*3./2.;
         }
         if (dcs<0.) dcs = 0.;// For safety

         if(j>0) {
            // Calculation of the cumulative diff. cross section
            difflambdavalue += preconst * dcs * (enrange[j+1]-enrange[j-1])/2.;
            // Next, sum up the cross sections for inner shell excitation (for omega prime above
            // 1800 eV) separately
            // - for debugging purposes only
            if(omegaprime>1800.*eV)
               innershelldifflambdavalue += preconst * dcs * (enrange[j+1]-enrange[j-1])/2.;
         } else {
            // Note: special case for j=0
            difflambdavalue  += preconst * dcs * (enrange[0]+enrange[1])/2.;
         }
         difflambda->push_back(difflambdavalue);
         // Update counters for omegaprime, and step omegaprime value
         if (j>=(TotBin-1)/5 && j<(TotBin-1)*3/5) jj++; else jj+=2;
         j++;
         omegaprime = enrange[j];
      }
      difflambdaformat->push_back(difflambda);
      // The total inverse mfp equals the last value of the cumulative differential function:
      lambda->push_back(difflambdavalue);
      if (verboseLevel>2) {
         if (i==0) {
            G4cout << setw(12) << "Energy" << " " << setw(12) << "InelasticMFP" << G4endl;
            G4cout << setw(12) << "(eV)"   << " " << setw(12) << "(nm)"       << G4endl;
         }
         if (difflambdavalue != 0.0) {
            G4cout << scientific << setw(12) << enrange[i]/eV << " " << setw(12)
                   << (1.0/difflambdavalue)/nanometer << G4endl;
         }
         //G4cout << enrange[i]/eV << "\t"  << difflambdavalue << "\t"
         //       << innershelldifflambdavalue << G4endl;
      }
      if (i>=(TotBin-1)/5 && i<(TotBin-1)*3/5) ii++; else ii+=2;
   }

   // Calculate and store the electron range as a function of energy.
   // This information is used by the 'range cut' check, that stops the tracking of an electron if
   // its distance to the nearest volume boundary is larger than a certain multiple of its range at
   // that point. Given this purpose, we do not really need a precise value for the range, but
   // rather a conservative (i.e. high) estimate suffices.
   // To this end, for each energy value we calculate the range as the sum of the IMFP and a
   // weighted average of the remaining ranges at energy E - omegaprime after one inelastic event,
   // where the weights are the differential cross sections for omegaprime. The range is somewhat
   // overestimated because the real energy loss in an event is omega rather than omegaprime.
   // Furthermore, the resulting range function is adjusted so that it gets monotonically increasing
   // as a function of energy, to account for the large mean free paths of very low-energy electrons
   // in some materials.
   //
   // On the other hand, by 'folding up' the trajectories of those low-energy electrons, elastic
   // scattering helps to reduce their effective range. This is taken into account by executing a
   // 'diffusion-like' correction rangevalue -> sqrt(2*tmfp*rangevalue) where 'tmfp' stands for the
   //  transport mean free path, or the typical distance a particle travels before information on
   // its original direction is lost. This information is generated by the
   // CADPhysicsSingleScattering class, which puts it in the 'tmfp' files, in such a way that the
   // function is monotonically increasing with energy (for energies above the vacuum level).
   // This way, the calculated 'tmfp' is always a conservative estimate of its 'effective' value
   // over the entire track of the particle.
   G4DataVector* rangeformat = new G4DataVector;
   G4DataVector* temprangeformat = new G4DataVector;
   G4double maxrangesofar = 0;
   G4double cdcsvalue = 0;
   size_t tmfpsize = tmfpenergies.size();
   for(i=0;i<TotBin;i++) {// Vary the energy
      G4double rangevalue;
      if( enrange[i] < barrier ) {
         // For energies below vacuum level, the range is (per definition) zero
         rangeformat->push_back(0.);
         temprangeformat->push_back(0.);
      } else {
         if ((*lambda)[i]==0) {
            rangeformat->push_back(DBL_MAX);
            temprangeformat->push_back(DBL_MAX);
            maxrangesofar = DBL_MAX;
         } else {
            rangevalue = 1./(*lambda)[i];// Start with the IMFP

            G4double endiff = 0.;
            G4int j=0;
            // and add a weighted average of the remaining range after one energy loss event
            while(enrange[j]<0.5*enrange[i] && enrange[j]<enrange[i]-fermiEnergy){
               // Vary omega prime
               // Find the differential cross section value
               if (j>0) cdcsvalue = (*(*difflambdaformat)[i])[j] - (*(*difflambdaformat)[i])[j-1];
               else cdcsvalue = (*(*difflambdaformat)[i])[0];
               // Find the index (intforrange) corresponding to the remaining energy after one
               // energy loss
               G4int intforrange = 0;
               for(int k=0;k<i;k++) {
                  if(enrange[k+1]>enrange[i]-enrange[j]) {
                     intforrange = k;
                     break;
                  }
               }
               if(j>0) endiff = (enrange[j+1] - enrange[j-1])/2;
               // Add the remaining range after one event, weighted by the corresponding
               // differential cross section and normalized to the total cross section
               rangevalue += (*temprangeformat)[intforrange] * cdcsvalue / (*lambda)[i] ;
               j++;
            }

            temprangeformat->push_back(rangevalue);

            // For the final result, correct the range for the elastic scattering 'folding' of the
            // trajectory

            // Find the 'tmfp' value
            G4double tmfp = DBL_MAX;
            for(size_t ii=0;ii<tmfpsize;ii++) {
               if(tmfpenergies[ii]>enrange[i]) {
                  if(ii==0) tmfp = tmfpdata[0]; else
                     tmfp = tmfpdata[ii-1] + (enrange[i]-tmfpenergies[ii-1])/
                        (tmfpenergies[ii]-tmfpenergies[ii-1])*(tmfpdata[ii]-tmfpdata[ii-1]);
                  break;
               }
            }
            // Perform the correction
            if(2.*tmfp<rangevalue) rangevalue = sqrt(2.*tmfp*rangevalue);
            // And make sure the final result is monotonically increasing with energy (so that the
            // electron can never have a smaller range than one with a lower energy)
            if(rangevalue>maxrangesofar) maxrangesofar = rangevalue;
            rangeformat->push_back(maxrangesofar);
         }
      }
      // Note that although this calculation is performed over the entire energy range, it is and
      // should be applied only for primary energies up to 1 keV. At higher energies the distances
      // between the energy values become too large for this method to be completely accurate, and
      // at the same time the additional gain in simulation speed would be rather limited.
   }
   delete temprangeformat;
   rangetable->push_back(rangeformat);
}

G4double CADPhysicsDI::GetRange(G4int index, G4double e)
{
   // Get the electron's range from the tabulated data
   G4DataVector* temprangeformat = (*rangetable)[index];
   G4double localrange = DBL_MAX;

   G4int j = FindIntLogen(e);
   if (j<=0) return 0.;
   if (j>TotBin-2) return DBL_MAX;
   // Rather than interpolating to the exact energy, we just take the value at the upper end of the
   // current energy 'bin'. This is faster than interpolation and we only need an upper limit to
   // the electron range rather than an exact value.
   localrange=(*temprangeformat)[j+1];
   return localrange;
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
   G4double electronrange = GetRange(findex,finalKinEnergy);
   if ( finalKinEnergy < fenergylimit ||
        (finalKinEnergy < 1.*keV && pstepsafety > 5.*electronrange) ) {
      // Too low remaining energy, kill the particle.
      //theEnergyDeposit += finalKinEnergy - ffermienergy;
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

   // Kill the particle if flagged by GetMeanFreePath (through boolean killthisone)
   // and also redo the energy check (see GetMeanFreePath for details)
   if (killthisone || kineticEnergy<fenergylimit) {
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
   G4double omegaprimecross = cross * G4UniformRand();

   // First find the index for the current kinetic energy
   G4int ii =  FindIntLogen(kineticEnergy);
   if (ii<16) {
      // At the very lowest energies (<0.02 eV) the differential cross section table is empty.
      // This 'pathetic case' can happen only if the potential barrier
      // (work function + fermi energy) of the material is around 0.02 eV or less.
      // N.B.: the cutoff value for ii depends on the value of TotBin!
      // Since we cannot determine omegaprime, return no particle change.
      return G4VContinuousDiscreteProcess::PostStepDoIt(track, step);
   }
   if (ii>TotBin-1) ii=TotBin-1;

   // First get the data set for the current material
   CADPhysicsDataTable* difflambdaformat = (*difflambdatable)[findex];

   // Then get the data vectors for the lower and upper limits of the current energy bin
   G4DataVector* difflambda1 = (*difflambdaformat)[ii];
   G4DataVector* difflambda2 = 0;
   if (ii==TotBin-1) difflambda2 = (*difflambdaformat)[ii];
   else difflambda2 = (*difflambdaformat)[ii+1];
   // and calculate the relative position in the energy bin so that we can interpolate the
   // differential cross sections for the exact kineticEnergy
   G4double enfrac = 0.;
   if (ii<TotBin-1) {
      enfrac = (kineticEnergy-enrange[ii])/(enrange[ii+1]-enrange[ii]);
      if (enfrac<0) enfrac = 0.;
   }

   G4int i=0;// i is the energy index for omegaprime
   G4int dlsize = difflambda1->size();
   G4double olddifflambdavalue = 0.;
   G4double difflambdavalue = (*difflambda1)[0] + enfrac*((*difflambda2)[0]-(*difflambda1)[0]);
   // Find the omegaprime energy 'bin' matching the current value of omegaprimecross
   while(omegaprimecross > difflambdavalue && i<dlsize-1) {
      i++;
      olddifflambdavalue = difflambdavalue;
      difflambdavalue = (*difflambda1)[i] + enfrac*((*difflambda2)[i]-(*difflambda1)[i]);
   }

   G4double omegaprime = enrange[0];

   // Finally, calculate omegaprime by interpolation. It can happen that we've reached the end
   // difflambda1 (the differential cross sections for the lower end of the energy bin) before
   // running out of omegaprimecross; i.e., i==dlsize-1 and still omegaprimecross > difflambdavalue.
   // In that case, interpolate up to the end of the next energy bin or up to
   // kineticEnergy - ffermienergy, whichever comes first.
   if (omegaprimecross > difflambdavalue) {
      G4double maxvalue = min(kineticEnergy - ffermienergy,enrange[i+1]);
      omegaprime  = enrange[i] +
         (omegaprimecross - difflambdavalue) / (cross - difflambdavalue) * (maxvalue  - enrange[i]);
   }
   // In other cases, just interpolate between both ends of the current energy bin for omegaprime.
   else if (i>0) {
      omegaprime = enrange[i-1] + (omegaprimecross - olddifflambdavalue)/
         (difflambdavalue - olddifflambdavalue) * (enrange[i]-enrange[i-1]);
   }
   // End of sampling omegaprime

   // Check for ionization of an inner shell
   G4double Ebind = 0.;// Initialization
   G4int Z=0;
   G4int shellId = 0;
   if (omegaprime>100.*eV) {
      // For large enough omegaprime, first see if the 'standard' Geant4 low energy utilities can
      // come up with a suitable inner shell to be ionized. This will be accepted only if the
      // binding energy is larger than 50 eV. The advantage of this method is that it knows about
      // the relative ionization cross sections of the different inner shells (at a certain
      // omegaprime energy) from which it can randomly select one. A disadvantage is that it knows
      // only about individual atoms, not about solid state electronic structure.
      // Hence the energy limits introduced here.
      G4double bindingEnergy = 0.;
      // Create a margin to catch discrepancies in inner shell energies between G4 data and the
      // df files
      G4double margin = 10.*eV;
      const G4MaterialCutsCouple* couple = track.GetMaterialCutsCouple();
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
   // If we haven't found an inner shell (yet), check the data from the df file
   if(Ebind==0 && innershells[findex]) {
      // This part is a bit less advanced. Simply find the largest inner shell ionization energy
      // in the list that is still smaller than omegaprime.
      G4DataVector* itsinnershells = (*innershelltable)[findex];
      size_t listsize = itsinnershells->size();
      size_t shellnr=0;
      G4bool toolarge = false;
      do {
         if(omegaprime>(*itsinnershells)[shellnr]) {
            // Binding energy set to this inner shell's energy
            Ebind = (*itsinnershells)[shellnr];
         }
         else {
            // Stop looking for next inner shell
            toolarge = true;
         }
         shellnr++;
      } while (!toolarge && shellnr<listsize);
   }
   // Finished looking for inner shells.

   // Sample omega
   // First set the limits. The maximum energy loss for the primary electron is such that its
   // kinetic energy after the event remains larger than omega-omegaprime, i.e. the *additional*
   // kinetic energy transfer to the secondary electron. omegamax equals the upper integration
   // limit of eq. 9 in Ashley, but corrected for the Fermi energy. The lower limit is set equal
   // to omegaprime (for inner shells) or to eq. 10 in Ashley (otherwise), where the latter
   // takes momentum conservation into account.
   G4double omegamax = 0.5*(kineticEnergy + omegaprime - ffermienergy);
   G4double omegamin = omegaprime;
   G4double Ebindprime = 0.;
   // EK: 2013-06-12 Modified if statement
   // Previous version:
   //    if(Ebind < 50.*eV) {
   //      // No 'inner shell' ionization, use Ashley's limit
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
      // For inner-shell ionization (Ebind > 50 eV) we substitute the Fermi-energy corrected
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
      if(fconductortype==0) {// Case 1. Metal, assume Fermi sea electron.
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
   // This part of the code has been derived from (and is essentially still the same as) that of
   // G4LowEnergyIonisation.

   // First transform to the shell potential - from both the (initial) primary kinetic energy and
   // the SE kinetic energy first subtract the Fermi energy, then add the binding energy *twice*.
   // The classical physical picture is that the electric potential near the given (inner shell)
   // orbit is twice the binding energy of that shell. The bound electron has a kinetic energy
   // equal to Ebind so that it's total energy is E_pot + E_kin = -2*Ebind + Ebind = -Ebind, as it
   // is expected to be.
   G4double deltaKinE = omega + Ebind;
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
   G4double dirz = cost1;
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

   // Get the local safety for the electron range check
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
   if (generateSecondaries) {
      // X-ray photons and Auger electrons from atom deexcitation
      if (shellId>0 && Z>5) {
         // Data are available only for Z>5 and shellId>0 only if a suitable inner shell was
         // found previously
         secondaryVector = deexcitation.GenerateParticles(Z, shellId);
      } else {
        // If no deexcitation, start with an empty vector
        secondaryVector = new  std::vector<G4DynamicParticle*>;
      }

      // Add the 'normal' secondary electron to the secondary particle vector
      G4DynamicParticle* theDeltaRay = new G4DynamicParticle();
      theDeltaRay->SetKineticEnergy(secondaryKineticEnergy);
      newdir = newdir.unit();
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
                  G4double electronrange = GetRange(findex,e);
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
                        passenergy = (e>fenergylimit &&
                                      (!rangecut || pstepsafety < 5.*electronrange));
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
   G4double electronrange = GetRange(findex,finalKinEnergy);
   if ( finalKinEnergy > fenergylimit &&
      (!rangecut || pstepsafety < 5.*electronrange ||
      finalKinEnergy > 1.*keV) ) {
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

   G4int stepnr = aTrack.GetCurrentStepNumber();
   //// Temporarily enlarged to a factor 10 higher in order to follow particles longer
   //if (stepnr>=100000 || kineticEnergy > 1.) {
   if (stepnr>=10000 || kineticEnergy > 1.) {
      *condition = Forced;
      killthisone = true;
      return DBL_MIN;
   }

   // Shortcut: for unchanged energy and material index, return the previous value
   if(findex == index)
      if (abs(fenergy-kineticEnergy)<min(1.*eV,0.1*kineticEnergy)) return flambda;
   findex = index;
   fenergy = kineticEnergy;

   ffermienergy = vec_fermieff[findex];
   fbarrier = vec_barrier[findex];
   fenergylimit = std::max(vec_minimumlimit[findex],fbarrier + cutenergy);
      // The energy limit is the user-defined cut energy
      // plus the energy required to get to vacuum level
   fconductortype = vec_conductortype[findex];

   // Check the physical state of the material. This process is not active in gases.
   const G4Material* material = couple->GetMaterial();
   if (material->GetState()==kStateGas)
   {
      flambda = DBL_MAX;
      return flambda;
   }

   // Look up the vector for the inverse mean free path for this material
   G4DataVector* lambda = (*lambdatable)[findex];

   // Find the energy 'bin' for this kinetic energy
   G4int j = FindIntLogen(kineticEnergy);

   if (j>=TotBin-1) {
      // 1/E extrapolation beyond the largest tabulated value
      cross = (*lambda)[TotBin-1] * enrange[TotBin-1] / kineticEnergy;
   } else {
      if (j<=0) return DBL_MAX;
      // Check against the energy limit. If the kinetic energy is too low, force PostStepDoIt to be
      // executed for this step and prepare to kill this particle.
      if (kineticEnergy<fenergylimit) {
         G4StepPoint* pPreStepPoint = aTrack.GetStep()->GetPreStepPoint();
         // Make sure that the particle is not at a boundary, to avoid the case where it is simply
         // being reflected on the outside of a surface
         if (pPreStepPoint->GetStepStatus() != fGeomBoundary) {
            killthisone = true;// Flag to kill the particle in PostStepDoIt
            *condition = Forced;
            return DBL_MAX;
         }
      }
      cross = (*lambda)[j]+(kineticEnergy-enrange[j])/(enrange[j+1]-enrange[j])
         *((*lambda)[j+1]-(*lambda)[j]);// Linear interpolation of the inverse MFP
   }

   G4double mfp = DBL_MAX;
   if (cross>0.) {// The MFP is the inverse of the inverse MFP
      mfp = 1./cross;
   }
   flambda = mfp;
   if (mfp==0.) {
      G4cout << findex << "\t" << kineticEnergy << "\t" << j << "\t"
             << cross << "\t" << flambda << "\t" << mfp << G4endl;
      abort();
   }
   return mfp;
}

G4int CADPhysicsDI::FindIntLogen(G4double kineticEnergy)
{
   // Energy on a log scale starting from zero at the lowest energy
   G4double logen = (log(kineticEnergy)/log(10.)) * (TotBin-1.)/10. + (TotBin-1.)*4./5.;
   // Correct for the 'denser' energy scale between 1 and 100 eV
   if (kineticEnergy>100.*eV) logen += (TotBin-1.)/5.; else {
      if (kineticEnergy>1.*eV) logen = logen*2. - (TotBin-1.)/5.;
   }
   // Round the result to give the index of the energy bin
   if (logen<0.) logen=0.;
   return (int)logen;
}


void CADPhysicsDI::PrintInfoDefinition()
{
   G4cout << "CADPhysicsDI: Process for inelastic scattering in (non-gaseous) materials,\n"
      << " based on the Dielectric Ionisation formalism. The implementation is loosely based on\n"
      << " the description in J.C. Ashley, J. Electr. Spectrosc. Rel. Phenom. 46: 199-214 (1988),\n"
      << " but with adjusted cross sections for ionization of inner shells. Furthermore, Auger\n"
      << " decay of inner-shell ionized atoms has been included.\n"
      << " The process is valid for energies up to 30 keV." << G4endl;
}
