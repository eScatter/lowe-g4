// Jun 14 2007: Added angular window selection for energy histogram and pdf
// may 21 2007: Added new version of angular output with energy window selection
// may 23 2007: Energy selection now via messenger for GeneralDetector


#ifndef SEMGeneralDetector_h
#define SEMGeneralDetector_h 1

#include "G4VSensitiveDetector.hh"
#include "SEMTuple.hh"
#include "SEMTupleElement.hh"
#include "SEMCountsPerCategory.hh"
#include "CADPhysicsUnits.hh"

//definition of CADDType
#include "CADDetectorData.hh"

#define DT_Undefined (0)
#define DT_SolidState (1)

enum OutputType {header,unit,data};

class SEMGeneralDetectorMessenger;

class SEMGeneralDetector : public G4VSensitiveDetector
{

  public:
      SEMGeneralDetector(G4String name);
      virtual ~SEMGeneralDetector();

      void OutputHits(string filename,SEMTuple<SEMTupleElement>& tuple);
      void OutputDeposits(string filename,SEMTuple<SEMTupleElement>& tuple,G4double mindeposit);
      void OutputEnergyPDF(string filename,SEMTuple<SEMTupleElement>& tuple,G4int PE);
      void OutputGammaPDF(string filename,SEMTuple<SEMTupleElement>& tuple,G4int PE);
      void OutputEnergyHistogram(string filename,SEMTuple<SEMTupleElement>& tuple,G4double BinSize,G4int PE);
      void OutputRhoHistogram(string filename,SEMTuple<SEMTupleElement>& tuple,G4double BinSize,G4int PE);
      void OutputRhoPDF(string filename,SEMTuple<SEMTupleElement>& tuple,G4int PE);
      void OutputAngleHistogram(string filename,SEMTuple<SEMTupleElement>& tuple,G4int PE, G4ThreeVector normal, G4int nbin);
      void OutputTimeHistogram(string filename,SEMTuple<SEMTupleElement>& tuple,G4double BinSize,G4int PE);
      void OutputCounts(string filename,SEMTuple<SEMTupleElement>& tuple,G4int PE);
      inline void SetLVname(G4String lvname) {DetectorLogVolumeName = lvname;} 
      inline G4String GetLVname() {return DetectorLogVolumeName;} 

  private:
      G4String DetectorName;
      G4String DetectorLogVolumeName;
      G4double Gauss(G4double x,G4double sigma);
      string   lastfilename;
      CADDtype fDetectorType;
      const G4double norm;
      G4double Emin,Emax;
      G4bool   Ewindow,EwindowChanged;
      G4double Amin,Amax;
      G4ThreeVector Adir;
      G4bool   Awindow,AwindowChanged;
      SEMGeneralDetectorMessenger *messenger;
      G4double EDeposit;
      G4int    fDetectorModel; // Determines what model to use to compute detector output from hits

      G4String DirName(G4String source);
      G4int    OpenOutputFile(G4String,ofstream &output);

      void   OutputEID(ofstream &output, SEMTupleElement tuple, OutputType mode);
      void   OutputTID(ofstream &output, SEMTupleElement tuple, OutputType mode);
      void   OutputPID(ofstream &output, SEMTupleElement tuple, OutputType mode);
      void   OutputEtot(ofstream &output, SEMTupleElement tuple, OutputType mode,G4int accuracy);
      void   OutputEkin(ofstream &output, SEMTupleElement tuple, OutputType mode,G4int accuracy);
      void   OutputEpot(ofstream &output, SEMTupleElement tuple, OutputType mode,G4int accuracy);
      void   OutputEdeposit(ofstream &output, SEMTupleElement tuple, OutputType mode,G4int accuracy);
      void   OutputEloss(ofstream &output, SEMTupleElement tuple, OutputType mode,G4int accuracy);
      void   OutputTime(ofstream &output, SEMTupleElement tuple, OutputType mode,G4int accuracy);
      void   OutputPosition(ofstream &output, SEMTupleElement tuple, OutputType mode,G4int accuracy);
      void   OutputVertexPosition(ofstream &output, SEMTupleElement tuple, OutputType mode,G4int accuracy);
      void   OutputDirection(ofstream &output, SEMTupleElement tuple, OutputType mode,G4int accuracy);
      void   OutputDepth(ofstream &output, SEMTupleElement tuple, OutputType mode,G4int accuracy);
      void   OutputMaxRadius(ofstream &output, SEMTupleElement tuple, OutputType mode,G4int accuracy);
      void   OutputRadius(ofstream &output, SEMTupleElement tuple, OutputType mode,G4int accuracy);
      void   OutputScatteringAngle(ofstream &output, SEMTupleElement tuple, OutputType mode,G4int accuracy);
      void   OutputLogicalVolume(ofstream &output, SEMTupleElement tuple, OutputType mode);
      void   OutputProcess(ofstream &output, SEMTupleElement tuple, OutputType mode);
      void   OutputCategory(ofstream &output, SEMTupleElement tuple, OutputType mode);
      void   OutputPotential(ofstream &output, G4double V, OutputType mode,G4int accuracy);

  public:
      inline void SetDetectorType(CADDtype t) { fDetectorType = t; }    
      inline CADDtype GetDetectorType() { return fDetectorType; }    
      inline G4String GetName() {return DetectorName; }
      inline void SetDetectorModel(int t) { fDetectorModel = t; }    
      inline G4int GetDetectorModel() { return fDetectorModel; }    

      G4bool fSelectEID,fSelectTID,fSelectPID,fSelectEtot,fSelectEkin,fSelectEdeposit,fSelectEloss,fSelectTime;
      G4bool fSelectPosition,fSelectVertexPosition,fSelectDirection;
      G4bool fSelectDepth,fSelectRadius,fSelectMaxRadius,fSelectScatteringAngle;
      G4bool fSelectLogicalVolume,fSelectProcess,fSelectCategory;
      G4bool fSelectPotential;
      inline void SetOutputEID() {fSelectEID = true; G4cout << "EID selected" << G4endl;}
      inline void SetOutputTID() {fSelectTID = true; G4cout << "TID selected" << G4endl;}
      inline void SetOutputPID() {fSelectPID = true; G4cout << "PID selected" << G4endl;}
      inline void SetOutputEtot() {fSelectEtot = true; G4cout << "total energy selected" << G4endl;}
      inline void SetOutputEkin() {fSelectEkin = true; G4cout << "kinetic energy selected" << G4endl;}
      inline void SetOutputEdeposit() {fSelectEdeposit = true; G4cout << "energy deposit selected" << G4endl;}
      inline void SetOutputEloss() {fSelectEloss = true; G4cout << "energy loss selected" << G4endl;}
      inline void SetOutputTime() {fSelectTime = true; G4cout << "hit time selected" << G4endl;}
      inline void SetOutputPosition() {fSelectPosition = true; G4cout << "position selected" << G4endl;}
      inline void SetOutputVertexPosition() {fSelectVertexPosition = true; G4cout << "vertex position selected" << G4endl;}
      inline void SetOutputDirection() {fSelectDirection = true; G4cout << "direction selected" << G4endl;}
      inline void SetOutputDepth() {fSelectDepth = true; G4cout << "depth selected" << G4endl;}
      inline void SetOutputRadius() {fSelectRadius = true; G4cout << "radius selected" << G4endl;}
      inline void SetOutputMaxRadius() {fSelectMaxRadius = true; G4cout << "maximum radius selected" << G4endl;}
      inline void SetOutputScatteringAngle() {fSelectScatteringAngle = true; G4cout << "scattering angle selected" << G4endl;}
      inline void SetOutputLogicalVolume() {fSelectLogicalVolume = true; G4cout << "logical volume selected" << G4endl;}
      inline void SetOutputProcess() {fSelectProcess = true; G4cout << "process selected" << G4endl;}
      inline void SetOutputCategory() {fSelectCategory = true; G4cout << "category selected" << G4endl;}
      inline void SetOutputPotential() {fSelectPotential = true; G4cout << "potential selected" << G4endl;}

      inline void UnsetOutputEID() {fSelectEID = false; G4cout << "EID deselected" << G4endl;}
      inline void UnsetOutputTID() {fSelectTID = false; G4cout << "TID deselected" << G4endl;}
      inline void UnsetOutputPID() {fSelectPID = false; G4cout << "PID deselected" << G4endl;}
      inline void UnsetOutputEtot() {fSelectEtot = false; G4cout << "total energy deselected" << G4endl;}
      inline void UnsetOutputEkin() {fSelectEkin = false; G4cout << "kinetic energy deselected" << G4endl;}
      inline void UnsetOutputEdeposit() {fSelectEdeposit = false; G4cout << "energy deposit deselected" << G4endl;}
      inline void UnsetOutputEloss() {fSelectEloss = false; G4cout << "energy loss deselected" << G4endl;}
      inline void UnsetOutputTime() {fSelectTime = false; G4cout << "hit time deselected" << G4endl;}
      inline void UnsetOutputPosition() {fSelectPosition = false; G4cout << "position deselected" << G4endl;}
      inline void UnsetOutputVertexPosition() {fSelectVertexPosition = false; G4cout << "vertex position deselected" << G4endl;}
      inline void UnsetOutputDirection() {fSelectDirection = false; G4cout << "direction deselected" << G4endl;}
      inline void UnsetOutputDepth() {fSelectDepth = false; G4cout << "depth deselected" << G4endl;}
      inline void UnsetOutputRadius() {fSelectRadius = false; G4cout << "radius deselected" << G4endl;}
      inline void UnsetOutputMaxRadius() {fSelectMaxRadius = false; G4cout << "maximum radius deselected" << G4endl;}
      inline void UnsetOutputScatteringAngle() {fSelectScatteringAngle = false; G4cout << "scattering angle deselected" << G4endl;}
      inline void UnsetOutputLogicalVolume() {fSelectLogicalVolume = false; G4cout << "logical volume deselected" << G4endl;}
      inline void UnsetOutputProcess() {fSelectProcess = false; G4cout << "process deselected" << G4endl;}
      inline void UnsetOutputCategory() {fSelectCategory = false; G4cout << "category deselected" << G4endl;}
      inline void UnsetOutputPotential() {fSelectPotential = false; G4cout << "potential deselected" << G4endl;}

  public:
      inline void SetEwindow(G4bool a) { Ewindow = a; EwindowChanged = true;}
      inline G4bool GetEwindow() {return Ewindow;}
      inline void SetEmin(G4double a) { 
         if (a>=0.0) 
            {Emin = a;} 
         else 
            {Emin = Emax+a;} 
         EwindowChanged = true; 
         G4cout << "Minimum energy = " << Emin/eV << " eV" << G4endl;
      }
      inline G4double GetEmin() {return Emin;}
      inline void SetEmax(G4double a) { Emax = a; EwindowChanged = true; G4cout << "Maximum energy = " << Emax/eV << " eV" << G4endl;}
      inline G4double GetEmax() {return Emax;}

      inline void SetAwindow(G4bool a) { Awindow = a; AwindowChanged = true;}
      inline G4bool GetAwindow() {return Awindow;}
      inline void SetAmin(G4double a) { Amin = a; AwindowChanged = true;}
      inline G4double GetAmin() {return Amin;}
      inline void SetAmax(G4double a) { Amax = a; AwindowChanged = true;}
      inline G4double GetAmax() {return Amax;}
      inline void SetAdir(G4ThreeVector a) { Adir = a; AwindowChanged = true;}
      inline G4ThreeVector GetAdir() {return Adir;}

      virtual void SetPE(int PE);
      virtual void ClearTuple();
      inline G4double GetEDeposit() {return EDeposit;}
      inline void SetEDeposit(G4double a) {EDeposit=a;}
      inline void AddEDeposit(G4double a) {EDeposit += a;}

      // EB: After help from Harald. The & tells the compiler to pass the argument tuple by reference, i.e.
      //     no copy is made. The const tells the compiler that we don't want to change the contents of tuple.
      //     If we try to do that, the compiler will see this at compile time.
      //     Original code: SEMCountsPerCategory GetDetectorCounts(SEMTuple<SEMTupleElement> tuple);
      SEMCountsPerCategory GetDetectorCounts(SEMTuple<SEMTupleElement>& tuple);
      G4double GetDetectorCurrent(SEMTuple<SEMTupleElement>& tuple);
};

#endif

