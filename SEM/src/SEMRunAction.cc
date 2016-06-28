#include "SEMRunAction.hh"
#include "SEMScanMessenger.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"

#include "SEMDetectorConstruction.hh"
#include "SEMPrimaryGeneratorAction.hh"
#include "SEMPlane.hh"
#include "SEMBulkDet.hh"
//#include "SEMCounter.hh"
#include "SEMGeneralDetector.hh"
#include "G4GeneralParticleSource.hh"
#include "G4SingleParticleSource.hh"
#include "G4SPSAngDistribution.hh"
#include "SEMCountsPerCategory.hh"

SEMRunAction::SEMRunAction()
{
   scanMessenger = new SEMScanMessenger(this);
   scanXdir = (G4ThreeVector)0;
   scanYdir = (G4ThreeVector)0;
   scanNx = 1;
   scanNy = 1;
   scanLx = 1.0*nanometer;
   scanLy = 1.0*nanometer;
   scanOrigin = (G4ThreeVector)0;
   validscan = false;
   validimage = false;
   maptrackend = false;
   mapRotation = G4ThreeVector(0.0,0.0,0.0);
}

SEMRunAction::~SEMRunAction()
{
   delete scanMessenger;
}

void SEMRunAction::BeginOfRunAction(const G4Run* /*aRun*/)
{
   static G4RunManager* run = G4RunManager::GetRunManager();
   SEMDetectorConstruction*  detector = (SEMDetectorConstruction*) run->GetUserDetectorConstruction();

   // Clear all data stored in the detectors before each new run
   G4int ndet = detector->GetNumberOfDetectors();
   for (int i=0; i<ndet; i++) {
       SEMGeneralDetector* GDetector = (SEMGeneralDetector*) G4SDManager::GetSDMpointer()->FindSensitiveDetector(detector->GetDetectorName(i));
       GDetector->ClearTuple();
       GDetector->SetEDeposit(0.0); // Clear the deposited energy before each run
    }
}

void SEMRunAction::EndOfRunAction(const G4Run* aRun)
{
   static G4RunManager* run = G4RunManager::GetRunManager();
   SEMPrimaryGeneratorAction* gun = (SEMPrimaryGeneratorAction*) run->GetUserPrimaryGeneratorAction(); 

   G4int nrun = aRun->GetNumberOfEvent();
   G4int ngun = gun->GetNumberOfParticles();
   PE = nrun*ngun;

   SEMDetectorConstruction*  detector = (SEMDetectorConstruction*) run->GetUserDetectorConstruction();

   // At end of run we tell the detectors how many primaries have been generated 
   G4int ndet = detector->GetNumberOfDetectors();
   for (int i=0; i<ndet; i++) {
      SEMGeneralDetector* Detector = (SEMGeneralDetector*) G4SDManager::GetSDMpointer()->FindSensitiveDetector(detector->GetDetectorName(i));
      Detector->SetPE(PE);
      if (Detector->GetEDeposit() != 0.0) {
         G4cout << setprecision(10) << "Energy deposited in detector " << detector->GetDetectorName(i) << " = " << Detector->GetEDeposit()/eV << " eV" << G4endl;
      }
   }
}

void SEMRunAction::MapOn(G4String filename,G4int ntheta, G4int nphi) {
   maptrackend = true;

   static G4RunManager* run = G4RunManager::GetRunManager();
   SEMPrimaryGeneratorAction* gun = (SEMPrimaryGeneratorAction*) run->GetUserPrimaryGeneratorAction();

   G4int nstore = gun->GetNumberOfParticles();
   gun->SetNumberOfParticles(1);
   gun->GetCurrentSource()->GetAngDist()->SetAngDistType("planar");

   G4double theta,phi,ct,st,cp,sp,s,t;
   G4double Pi = 4.0*atan(1.0);

   G4double mapx,cmx,smx, mapy,cmy,smy, mapz,cmz,smz, x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3;
 
   mapx = mapRotation.x() * Pi/180.0; cmx = cos(mapx); smx = sin(mapx);
   mapy = mapRotation.y() * Pi/180.0; cmy = cos(mapy); smy = sin(mapy);
   mapz = mapRotation.z() * Pi/180.0; cmz = cos(mapz); smz = sin(mapz);
   G4cout << "MapRotation = " << mapRotation << G4endl;

   vector < vector < G4Colour > > TrackEndMap;

   if (nphi < 0) { // Polar plot
      G4int ii,jj;

      if (ntheta % 2 == 0) ntheta++; // ntheta should be odd
      G4cerr << "Ntheta = " << ntheta << G4endl;
      TrackEndMap.resize(ntheta);
      nphi = ntheta;
      for (int i=0; i<ntheta; i++) TrackEndMap[i].resize(nphi);

      for (int i = 0; i < ntheta; i++) {
         G4cout << i << G4endl;
         ii = i - (ntheta-1)/2;
         for (int j = 0; j < nphi; j++) {
            jj = j - (nphi-1)/2;
            s = (double) sqrt((double)(ii*ii)+(double)(jj*jj))/(double) ((ntheta-1)/2);
            if (s>1.0) {
               TrackEndMap[i][j] = G4Colour(0.0,0.0,0.0); // Outside circle
            } else {
               // Inside circle
               theta =  s * Pi/2.0;
               ct = cos(theta);
               st = sin(theta);
               if ((ii==0)&&(jj==0)) {
                  TrackEndMap[i][j] = G4Colour(1.0,1.0,1.0); // Center of the plot
               } else {
                  phi = atan2((double)jj,(double)ii);
//                  G4cerr << theta*180.0/Pi  << " " << phi*180.0/Pi << " ";
                  cp = cos(phi);
                  sp = sin(phi);
                  // Perform rotations around the angles given in mapRotation
                  x = st*cp;
                  y = st*sp;
                  z = ct;
                  // x-axis:
                  x1 = x;
                  y1 = cmx*y-smx*z;
                  z1 = smx*y+cmx*z;
                  // y-axis:
                  x2 = cmy*x1+smy*z1;
                  y2 = y1;
                  z2 = -smy*x1+cmy*z1;
                  // z-axis:
                  x3 = cmz*x2-smz*y2;
                  y3 = smz*x2+cmz*y2;
                  z3 = z2;
//                  G4cout << "x0: " << x << " " << y << " " << z << G4endl;
//                  G4cout << "x1: " << x1 << " " << y1 << " " << z1 << G4endl;
//                  G4cout << "x2: " << x2 << " " << y2 << " " << z2 << G4endl;
//                  G4cout << "x3: " << x3 << " " << y3 << " " << z3 << G4endl;
         
                  //gun->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(st*cp,st*sp,ct));
                  gun->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(x3,y3,z3));
                  run->BeamOn(1); 
                  TrackEndMap[i][j] = mapcolour;
//                  G4cerr << mapcolour.GetRed() << " " << mapcolour.GetGreen() << " " << mapcolour.GetBlue() << G4endl;
               }
            }
         }
      }
   } else { // Original plot (x=phi, y=theta)
      TrackEndMap.resize(ntheta);
      for (int i=0; i<ntheta; i++) TrackEndMap[i].resize(nphi);

      for (int itheta = 0; itheta < ntheta; itheta++) {
         G4cout << itheta << G4endl;
         s = (double) itheta / (double) ntheta;
         theta =  s * Pi/2.0;
         ct = cos(theta);
         st = sin(theta);
         for (int iphi = 0; iphi < nphi; iphi++) {
//            G4cerr << itheta << " " << iphi << " ";
            t = (double) iphi / (double) nphi;
            phi = t * 2.0*Pi;
//            G4cerr << theta*180.0/Pi  << " " << phi*180.0/Pi << " ";
            cp = cos(phi);
            sp = sin(phi);
            // Perform rotations around the angles given in mapRotation
            x = st*cp;
            y = st*sp;
            z = ct;
            // x-axis:
            x1 = x;
            y1 = cmx*y-smx*z;
            z1 = smx*y+cmx*z;
            // y-axis:
            x2 = cmy*x1+smy*z1;
            y2 = y1;
            z2 = -smy*x1+cmy*z1;
            // z-axis:
            x3 = cmz*x2-smz*y2;
            y3 = smz*x2+cmz*y2;
            z3 = z2;
//            G4cout << "x0: " << x << " " << y << " " << z << G4endl;
//            G4cout << "x1: " << x1 << " " << y1 << " " << z1 << G4endl;
//            G4cout << "x2: " << x2 << " " << y2 << " " << z2 << G4endl;
//            G4cout << "x3: " << x3 << " " << y3 << " " << z3 << G4endl;
         
            //gun->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(st*cp,st*sp,ct));
            gun->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(x3,y3,z3));
            run->BeamOn(1); 
            TrackEndMap[itheta][iphi] = mapcolour;
//            G4cerr << mapcolour.GetRed() << " " << mapcolour.GetGreen() << " " << mapcolour.GetBlue() << G4endl;
         }
      }
   }

   // Output the TrackEndMap
   if (filename != "") {
      unsigned char r,g,b;
      ofstream ppmfile;
      ppmfile.open(filename);
      ppmfile << "P6";

      ppmfile.close();
      ppmfile.open(filename,ios::app|ios::binary);
      r = 10;
      ppmfile.write((char*) &r, sizeof(unsigned char));
      ppmfile.close();
      ppmfile.open(filename,ios::app);

      ppmfile << "# Creator SEM v2.0 (special)" ;

      ppmfile.close();
      ppmfile.open(filename,ios::app|ios::binary);
      ppmfile.write((char*) &r, sizeof(unsigned char));
      ppmfile.close();
      ppmfile.open(filename,ios::app);

      ppmfile << "# Image dimensions: " << ntheta << " x " << nphi ;

      ppmfile.close();
      ppmfile.open(filename,ios::app|ios::binary);
      ppmfile.write((char*) &r, sizeof(unsigned char));
      ppmfile.close();
      ppmfile.open(filename,ios::app);

      ppmfile << nphi << " " << ntheta ;

      ppmfile.close();
      ppmfile.open(filename,ios::app|ios::binary);
      ppmfile.write((char*) &r, sizeof(unsigned char));
      ppmfile.close();
      ppmfile.open(filename,ios::app);
      
      ppmfile << "255";
      ppmfile.close();
      ppmfile.open(filename,ios::app|ios::binary);
      ppmfile.write((char*) &r, sizeof(unsigned char));

      for (int i=0; i<ntheta; i++) {
         for (int j=0; j<nphi; j++) {
            r = (unsigned char) (255.0*TrackEndMap[i][j].GetRed());
            g = (unsigned char) (255.0*TrackEndMap[i][j].GetGreen());
            b = (unsigned char) (255.0*TrackEndMap[i][j].GetBlue());
            ppmfile.write((char*) &r, sizeof(unsigned char));
            ppmfile.write((char*) &g, sizeof(unsigned char));
            ppmfile.write((char*) &b, sizeof(unsigned char));
         }
      }
      ppmfile.close();
   }


   gun->SetNumberOfParticles(nstore);

   maptrackend = false;
}

void SEMRunAction::ScanOn(G4int n_event) {

   G4ThreeVector currentpos;
   G4ThreeVector scanstep = (G4ThreeVector)0;

   // Define step along the scanline
   if (scanNx>1) scanstep = scanLx*scanXdir/(scanNx-1);

   static G4RunManager* run = G4RunManager::GetRunManager();
   SEMPrimaryGeneratorAction* userPrimaryGeneratorAction = (SEMPrimaryGeneratorAction*) run->GetUserPrimaryGeneratorAction();
   SEMDetectorConstruction*  detector = (SEMDetectorConstruction*) run->GetUserDetectorConstruction();

   // Save original gun position
   G4ThreeVector tempstorepos = userPrimaryGeneratorAction->GetPosition();

   // Get number of detectors
   G4int ndet = detector->GetNumberOfDetectors();
   fLineScan.resize(ndet);

   // Scan along the line
   currentpos = scanOrigin-0.5*scanLx*scanXdir;
   for (int i=0; i<ndet; i++) fLineScan[i].resize(scanNx);
   fPosLineScan.resize(scanNx);
   SEMCountsPerCategory counts;
   for (int i=0;i<scanNx;i++) {
      fPosLineScan[i] = currentpos;
      userPrimaryGeneratorAction->SetPosition(currentpos);
//      for (int j=0; j<n_event; j++) run->BeamOn(1);   // Alternative
      run->BeamOn(n_event);
      G4cout << "Current position = " << currentpos/nanometer << " nm ";
      G4cout << "Number of hits: ";
      for (int j=0; j<ndet; j++) {
          switch(detector->GetDetectorType(j)) {
             case fPlane:
                {
                  SEMPlane* Plane = (SEMPlane*) G4SDManager::GetSDMpointer()->FindSensitiveDetector(detector->GetDetectorName(j));
                  if (scanPrefix != "") {
                     char fname[80];
                     sprintf(fname,"%s_%s_%03d.hits",scanPrefix.data(),detector->GetDetectorName(j).data(),i);
                     Plane->OutputHits(fname);
                  }
                  counts = Plane->GetCounts();
                  G4cout << detector->GetDetectorName(j) << "(" << counts.nTotal << " " << counts.nBSE << " " << counts.nSE << ") ";
                }
                break;
             case fBulk:
                {
                  SEMBulkDet* Det = (SEMBulkDet*) G4SDManager::GetSDMpointer()->FindSensitiveDetector(detector->GetDetectorName(j));
                  if (scanPrefix != "") {
                     char fname[80];
                     sprintf(fname,"%s_%s_%03d.hits",scanPrefix.data(),detector->GetDetectorName(j).data(),i);
                     Det->OutputHits(fname);
                  }
                  counts = Det->GetCounts();
                  G4cout << detector->GetDetectorName(j) << "(" << counts.nTotal << " " << counts.nBSE << " " << counts.nSE << " " << counts.nSE3 << ") ";
                }
               break;
          }
          fLineScan[j][i] = counts;
       }
       G4cout << G4endl;
      
      currentpos += scanstep;
   }
 
   // Restore original gun position
   userPrimaryGeneratorAction->SetPosition(tempstorepos);
   validscan = true;
}

void SEMRunAction::ScanOn(G4int n_event,G4String posfilename) {

   // Perform a scan over points given in the file
   // File format: 
   // Line1:  n   = number of points in the scan
   // Line2..n+1: x,y,z = position of the gun in mm

   // Read in the position file
   if (posfilename!="") {
      ifstream posfile;
      posfile.open(posfilename);
      if (!posfile) {
         G4cerr << "Position file could not be opened. Aborting linescan..." << G4endl;
         return;
      } else {
         int n;
         posfile >> n;
         SetScanNx(n); // invalidates old linescan and/or image and sets scanNx
         fPosLineScan.resize(scanNx);
         for (int i=0; i<n; i++) {
            posfile >> fPosLineScan[i];
         }
      }
   }

   G4ThreeVector currentpos;
   G4ThreeVector scanstep = (G4ThreeVector)0;

   static G4RunManager* run = G4RunManager::GetRunManager();
   SEMPrimaryGeneratorAction* userPrimaryGeneratorAction = (SEMPrimaryGeneratorAction*) run->GetUserPrimaryGeneratorAction();
   SEMDetectorConstruction*  detector = (SEMDetectorConstruction*) run->GetUserDetectorConstruction();

   // Save original gun position
   G4ThreeVector tempstorepos = userPrimaryGeneratorAction->GetPosition();

   // Get number of detectors
   G4int ndet = detector->GetNumberOfDetectors();
   fLineScan.resize(ndet);

   // Scan along the line
   for (int i=0; i<ndet; i++) fLineScan[i].resize(scanNx);
   SEMCountsPerCategory counts;
   for (int i=0;i<scanNx;i++) {
      currentpos = fPosLineScan[i];
      userPrimaryGeneratorAction->SetPosition(currentpos);
//      for (int j=0; j<n_event; j++) run->BeamOn(1);   // Alternative
      run->BeamOn(n_event);
      G4cout << "Current position = " << currentpos/nanometer << " nm ";
      G4cout << "Number of hits: ";
      for (int j=0; j<ndet; j++) {
          switch(detector->GetDetectorType(j)) {
             case fPlane:
                {
                  SEMPlane* Plane = (SEMPlane*) G4SDManager::GetSDMpointer()->FindSensitiveDetector(detector->GetDetectorName(j));
                  if (scanPrefix != "") {
                     char fname[80];
                     sprintf(fname,"%s_%s_%03d.hits",scanPrefix.data(),detector->GetDetectorName(j).data(),i);
                     Plane->OutputHits(fname);
                  }
                  counts = Plane->GetCounts();
                  G4cout << detector->GetDetectorName(j) << "(" << counts.nTotal << " " << counts.nBSE << " " << counts.nSE << ") ";
                }
                break;
             case fBulk:
                {
                  SEMBulkDet* Det = (SEMBulkDet*) G4SDManager::GetSDMpointer()->FindSensitiveDetector(detector->GetDetectorName(j));
                  if (scanPrefix != "") {
                     char fname[80];
                     sprintf(fname,"%s_%s_%03d.hits",scanPrefix.data(),detector->GetDetectorName(j).data(),i);
                     Det->OutputHits(fname);
                  }
                  counts = Det->GetCounts();
                  G4cout << detector->GetDetectorName(j) << "(" << counts.nTotal << " " << counts.nBSE << " " << counts.nSE << " " << counts.nSE3 << ") ";
                }
               break;
          }
          fLineScan[j][i] = counts;
       }
       G4cout << G4endl;
      
   }
 
   // Restore original gun position
   userPrimaryGeneratorAction->SetPosition(tempstorepos);
   validscan = true;
}

void SEMRunAction::OutputLinescan(G4String filename,G4String name)
{
   if (validscan) {
      static G4RunManager* run = G4RunManager::GetRunManager();
      SEMDetectorConstruction*  detector = (SEMDetectorConstruction*) run->GetUserDetectorConstruction();

      G4int ndet = detector->GetNumberOfDetectors();
      G4int idet = -1;
      for (int i=0; i<ndet; i++) {
         if (detector->GetDetectorName(i) == name) {
            idet = i;
         }
      }
      if (idet == -1) {
         G4cerr << "ERROR: Detector " << name << " not found. No output generated." << G4endl;
      } else if (filename != "") {
         SEMGeneralDetector* GenDetector = (SEMGeneralDetector*) G4SDManager::GetSDMpointer()->FindSensitiveDetector(detector->GetDetectorName(idet));

         ofstream scanfile;

         scanfile.open(filename);
         scanfile << "Result of the line scan for detector " << name << G4endl;
         scanfile <<    " Pos.                                  Total        SE       BSE";
         if (GenDetector->GetEwindow()) {
            scanfile << "       Energy Window";
         }
         if (GenDetector->GetAwindow()) {
            scanfile << "       Angle Window";
         }
         scanfile << G4endl;

         scanfile << " (nm,nm,nm) " << G4endl;

         for (int i=0; i<scanNx; i++) {
            scanfile << " " << scientific << setw(14) << setprecision(6) 
                     << fPosLineScan[i].x()/nanometer << " " 
                     << fPosLineScan[i].y()/nanometer << " " 
                     << fPosLineScan[i].z()/nanometer << " "
                     << setw(6) << fLineScan[idet][i].nTotal << "      " 
                     << setw(6) << fLineScan[idet][i].nSE << "      " 
                     << setw(6) << fLineScan[idet][i].nBSE << "      ";
            if (GenDetector->GetEwindow()) {
               scanfile << setw(6) << fLineScan[idet][i].nEkinWindow;
            }
            if (GenDetector->GetAwindow()) {
               scanfile << setw(6) << fLineScan[idet][i].nAngWindow;
            }
            scanfile << G4endl;
         }
         scanfile.close();
      }
   } else {
       G4cerr << "ERROR: There is either no linescan present or you changed the number of pixels along the scan. No ouput generated!" << G4endl;
   }
}

void SEMRunAction::ImageOn(G4int n_event,G4String intermediatefilename) {

   G4cout << "Generating image. Temporary images are saved during the scan." << G4endl;

   G4ThreeVector currentpos,startline;
   G4ThreeVector scanstepX = (G4ThreeVector)0;
   G4ThreeVector scanstepY = (G4ThreeVector)0;

   // Define steps in X and Y directions of the image
   if (scanNx>1) scanstepX = scanLx*scanXdir/(scanNx-1);
   if (scanNy>1) scanstepY = scanLy*scanYdir/(scanNy-1);

   static G4RunManager* run = G4RunManager::GetRunManager();
   SEMPrimaryGeneratorAction* userPrimaryGeneratorAction = (SEMPrimaryGeneratorAction*) run->GetUserPrimaryGeneratorAction();
   SEMDetectorConstruction*  detector = (SEMDetectorConstruction*) run->GetUserDetectorConstruction();

   // Save original gun position
   G4ThreeVector tempstorepos = userPrimaryGeneratorAction->GetPosition();

   // Get number of detectors
   G4int ndet = detector->GetNumberOfDetectors();
   fImage.resize(ndet);  // We store the full image for all available detectors

   // Scan along the line
   currentpos = scanOrigin-0.5*scanLx*scanXdir;
   for (int i=0; i<ndet; i++) {
      fImage[i].resize(scanNy);
      for (int j=0; j<scanNy; j++) {
         fImage[i][j].resize(scanNx);
      }
   }
   SEMCountsPerCategory counts;

   // Allow intermediate images to be saved
   validimage = true;

   // Scan the image (X direction first)
   ofstream asciifile;
   if (intermediatefilename=="") {
      intermediatefilename="Intermediate.txt";
   }
   asciifile.open(intermediatefilename);
   asciifile << "Image (ASCII) intermediate results:" << G4endl;
   startline = scanOrigin-(scanLx*scanXdir+scanLy*scanYdir)/2.0;
   
   for (int i=0;i<scanNy;i++) {
      currentpos = startline;
      for (int j=0;j<scanNx;j++) {
         userPrimaryGeneratorAction->SetPosition(currentpos);
         G4cout << i << " " << j << " " << currentpos.x()/nanometer << " " << currentpos.y()/nanometer << " " ;
         asciifile << i << " " << j << " " << currentpos.x()/nanometer << " " << currentpos.y()/nanometer << " " ;
         run->BeamOn(n_event);
         for (int k=0; k<ndet; k++) {
            switch(detector->GetDetectorType(k)) {
               case fPlane:
                  {
                    SEMPlane* Plane = (SEMPlane*) G4SDManager::GetSDMpointer()->FindSensitiveDetector(detector->GetDetectorName(k));
                    counts = Plane->GetCounts();
                    G4cout << detector->GetDetectorName(k) << "(" << counts.nTotal << ") ";
                    // If no model is attached the current will just be the total number of hits, otherwise it is the current
                    // using the detector model
                    asciifile << counts.Current << " ";
                  }
                  break;
               case fBulk:
                  {
                    SEMBulkDet* Det = (SEMBulkDet*) G4SDManager::GetSDMpointer()->FindSensitiveDetector(detector->GetDetectorName(k));
                    counts = Det->GetCounts();
                    G4cout << detector->GetDetectorName(k) << "(" << counts.nTotal << ") ";
                    asciifile << counts.nTotal << " ";
                  }
                 break;
            }
            fImage[k][i][j] = counts;
         }
         G4cout << G4endl;
         asciifile << G4endl;

         currentpos += scanstepX;
      }
      // Write images to file (after each completed line), normalized to be on the same scale (for now)
      for (int k=0; k<ndet; k++) {
         G4String filename;
         filename = detector->GetDetectorName(k)+"Image.pgm";
         WritePgm(filename,fImage[k]);
      }
      startline += scanstepY;
   }

   // Restore original gun position
   userPrimaryGeneratorAction->SetPosition(tempstorepos);
}

void SEMRunAction::WritePgm(G4String filename,vector <vector <SEMCountsPerCategory> > Image) 
{
   G4int category = 8; // The default image uses the model. If no model exists this will be the total number of hits
   G4bool autoscale = true;
   G4bool logscale = false;
   G4double min =0.0 ,max = 0.0;
   WritePgm(filename, Image, category, autoscale, logscale, &min, &max);
}

void SEMRunAction::OutputPgm(G4String filename,G4String detectorname, G4int category, G4bool autoscale, G4double min, G4double max) {

   if (validimage) {
      static G4RunManager* run = G4RunManager::GetRunManager();
      SEMDetectorConstruction*  detector = (SEMDetectorConstruction*) run->GetUserDetectorConstruction();
   
      // Find the detector first
      G4int ndet = detector->GetNumberOfDetectors();
      G4int idet = -1;
      for (int i=0; i<ndet; i++) {
          if (detector->GetDetectorName(i) == detectorname) idet = i;
      }
      if (idet == -1) {
         G4cerr << "WARNING: Detector not found, no output generated.";
      } else {
         G4bool logscale = false;
         WritePgm(filename,fImage[idet],category,autoscale,logscale,&min,&max);
      }
   } else {
       G4cerr << "ERROR: There is either no valid image present or you changed the number of pixels. No ouput generated!" << G4endl;
   }
}
void SEMRunAction::OutputASCII(G4String filename,G4String detectorname) {

   if (validimage) {
      static G4RunManager* run = G4RunManager::GetRunManager();
      SEMDetectorConstruction*  detector = (SEMDetectorConstruction*) run->GetUserDetectorConstruction();
   
      G4ThreeVector currentpos,startline;
      G4ThreeVector scanstepX = (G4ThreeVector)0;
      G4ThreeVector scanstepY = (G4ThreeVector)0;

      // Define steps in X and Y directions of the image
      if (scanNx>1) scanstepX = scanLx*scanXdir/(scanNx-1);
      if (scanNy>1) scanstepY = scanLy*scanYdir/(scanNy-1);

      startline = scanOrigin-(scanLx*scanXdir+scanLy*scanYdir)/2.0;

      // Find the detector first
      G4int ndet = detector->GetNumberOfDetectors();
      G4int idet = -1;
      for (int i=0; i<ndet; i++) {
          if (detector->GetDetectorName(i) == detectorname) idet = i;
      }
      if (idet == -1) {
         G4cerr << "WARNING: Detector not found, no output generated.";
      } else {
         ofstream ASCIIfile;
         ASCIIfile.open(filename);
         ASCIIfile << "Detected hits at each position for detector " << detectorname << G4endl;
         ASCIIfile << "NOTE: The three columns for energy and angular windows are only valid if " << G4endl;
         ASCIIfile << "      an energy or angle window has been selected before imageOn. " << G4endl;
         ASCIIfile << "      The selection criteria are not saved!" << G4endl;
         ASCIIfile << "NOTE: The detector current is only valid if a model was selected before imageOn " << G4endl;
         ASCIIfile << " x (pixels)  y (pixels)  x (nm) y (nm)  z(nm) Total     BSE      SE          Windows: Etot, Ekin, Angle    Current" << G4endl;
         ASCIIfile << G4endl;
         G4int Nx = fImage[idet][0].size();
         G4int Ny = fImage[idet].size();
         for (int i=0; i<Ny; i++) {
            currentpos = startline;
            for (int j=0; j<Nx; j++) {
               ASCIIfile << i << " " << j << " " << scientific << setprecision(5) 
                         << currentpos.x()/nanometer << " " << currentpos.y()/nanometer << " " << currentpos.z()/nanometer << " " 
                         << fImage[idet][i][j].nTotal << " "
                         << fImage[idet][i][j].nBSE << " "
                         << fImage[idet][i][j].nSE << " "
                         << fImage[idet][i][j].nEtotWindow << " "
                         << fImage[idet][i][j].nEkinWindow << " "
                         << fImage[idet][i][j].nAngWindow << " "
                         << fImage[idet][i][j].Current << " ";
               ASCIIfile << G4endl;
               currentpos += scanstepX;
            }
            startline += scanstepY;
         }
      }
   } else {
       G4cerr << "ERROR: There is either no valid image present or you changed the number of pixels. No ouput generated!" << G4endl;
   }
}

void SEMRunAction::WritePgm(G4String filename,vector <vector <SEMCountsPerCategory> > Image, G4int category, G4bool autoscale, G4bool logscale, G4double *min, G4double *max) {
   // Category:
   // 1) nTotal
   // 2) nBSE
   // 3) nSE
   // 4) nSE3
   // 5) nEtotWindow
   // 6) nEkinWindow
   // 7) nAngWindow
   // 8) Current
   // Autoscale: true if autoscaling is on, otherwise min and max are used
   // Logscale: if true, the image is plotted on a logarithmic scale, autoscaling is always used in this case.

   // Reserve memory for a copy of the image. 
   G4int Nx = Image[0].size();
   G4int Ny = Image.size();
   vector <vector<G4double > > image;
   image.resize(Ny);
   for (int i=0; i<Ny; i++) {
      image[i].resize(Nx);
   }

   G4double Min,Max;

   // Generate the requested image, determine Min,Max on the fly
   Min=Max=-1.0;
   for (int i=0; i<Ny; i++) {
      for (int j=0; j<Nx; j++) {
         switch (category) {
            case 1:
               image[i][j] = (G4double) (Image[i][j].nTotal);
               break;
            case 2:
               image[i][j] = (G4double) (Image[i][j].nBSE);
               break;
            case 3:
               image[i][j] = (G4double) (Image[i][j].nSE);
               break;
            case 4: // Only for Det type detector, otherwise this will be rubbish
               image[i][j] = (G4double) (Image[i][j].nSE3);
               break;
            case 5:
               image[i][j] = (G4double) (Image[i][j].nEtotWindow);
               break;
            case 6:
               image[i][j] = (G4double) (Image[i][j].nEkinWindow);
               break;
            case 7:
               image[i][j] = (G4double) (Image[i][j].nAngWindow);
               break;
            case 8:
               image[i][j] = Image[i][j].Current;
               break;
            default:
               G4cerr << "WARNING: Unknown category selected for pgm output. No output generated." << G4endl;
               return;
         }
         if ((image[i][j] < Min) || (Min < 0.0)) Min=image[i][j];
         if ((image[i][j] > Max) || (Max < 0.0)) Max=image[i][j];
      }
   }
  
   // Perform the autoscaling if requested
   if ((autoscale) || (logscale)){
      *min = Min;
      *max = Max;
      if (Min == Max) {
         G4cerr << "WARNING: Data range is zero in pgm output. No output generated." << G4endl;
         return;
      }
   }
   if (logscale) {
      if (Min==0.0) Min = 1.0;  // Prevent taking the log of zero later on.
   }

   unsigned char byte;
   ofstream pgmfile;
   pgmfile.open(filename);
   pgmfile << "P5" << G4endl;
   pgmfile << "# Creator SEM v2.0" << G4endl;
   pgmfile << "# Image dimensions: " << scanXdir.mag()/nanometer << " x " << scanYdir.mag()/nanometer << " nm" << G4endl;
   pgmfile << scanNx << " " << scanNy << G4endl;
   pgmfile << "255" << G4endl;

   pgmfile.close();
   pgmfile.open(filename,ios::app|ios::binary);
   for (int i=0; i<scanNy; i++) {
      for (int j=0; j<scanNx; j++) {
         if (logscale) {
            if (image[i][j]<=0.0) {
               byte = 0;
            } else {
               byte = (char) (255.0*(log10((double)image[i][j])-log10((double)Min))/(log10((double)Max)-log10((double)Min)));
            }
         } else {
            if (image[i][j] < Min) {
               byte=0;
            } else if (image[i][j] > Max) {
               byte=255;
            } else {
               byte=(char) (255.0*(double)(image[i][j] - Min)/(double)(Max - Min));
            }
         }
         pgmfile.write((char *) &byte, sizeof(char));
      }
   }
   pgmfile.close();
}

