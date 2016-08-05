// Sep 27 2010: Added printing of Energy Deposits
// Jul 31 2007: Added printing of depth (minimum z-coordinate) and maximum rho (w.r.t z-axis) to output in hits file
// Jul 18 2007: Changed setw(14) to setw(15) in OutputHits to generate proper output in Windows version
// Jun 14 2007: Added angular window selection to Energy histogram and PDF
//              Changed algorithm for PDF generation. It is now slightly simpler but still slow for large hit collections.
// May 21 2007: Added energy window selection to Angular histogram
// May 22 2007: Added output of arrival statistics to ouput of time histogram (no separate command yet)
// May 23 2007: Extra option to select an energy window (kinetic energy) which will limit the output of
//                outputHits
//                outputCounts
//                outputAngleHistogram
//                to hits that are within this energy window.


#include "SEMGeneralDetector.hh"
#include "SEMGeneralDetectorMessenger.hh"

#include "G4HCofThisEvent.hh"
#include "G4RunManager.hh"
#include "SEMPrimaryGeneratorAction.hh"
#include "SEMDetectorConstruction.hh"
#include <sys/stat.h>


using namespace std;

SEMGeneralDetector::SEMGeneralDetector(G4String name)
	:G4VSensitiveDetector(name),
	norm(sqrt(twopi))
{ 
   DetectorName = name;
   lastfilename="";
   Ewindow = false;
   Emin = 0.0;
   Emax = 1.0e30*eV;
   EwindowChanged = true;
   Awindow = false;
   Amin = -180.0*deg;
   Amax = 180.0*deg;
   AwindowChanged = true;
   fDetectorModel = DT_Undefined; // Require user to define a detector type if currents are needed

   // Standard selection for OutputHits command for backwards compatibility
   fSelectEID=true;
   fSelectTID=true;
   fSelectPID=false;
   fSelectEtot=true;
   fSelectEkin=true;
   fSelectEdeposit=false;
   fSelectEloss=false;
   fSelectTime=true;
   fSelectPosition=true;
   fSelectVertexPosition=true;
   fSelectDirection=true;
   fSelectDepth=true;
   fSelectRadius=false;
   fSelectMaxRadius=true;
   fSelectScatteringAngle=false;
   fSelectLogicalVolume=true;
   fSelectProcess=true;
   fSelectCategory=true;
   fSelectPotential=false;


   messenger = new SEMGeneralDetectorMessenger(this);
}

SEMGeneralDetector::~SEMGeneralDetector(){

   delete messenger;
}

void SEMGeneralDetector::SetPE(int /*PE*/) {
   G4Exception("Calling SEMGeneralDetector::SetPE()","SetPE()Error",JustWarning,"SEMGeneralDetector(1)");
   exit(1);
}

void SEMGeneralDetector::ClearTuple() {
}

void SEMGeneralDetector::OutputEID(ofstream &output,SEMTupleElement tupel,OutputType mode) {
   if (fSelectEID) {
      switch (mode) {
         case header:
            output << "EventID ";
            break;
         case unit:
            output << "    (-) ";
            break;
         default:
            output << setw(7) << tupel.EventID << " ";
      }
   }
}

void SEMGeneralDetector::OutputTID(ofstream &output,SEMTupleElement tupel,OutputType mode) {
   if (fSelectTID) {
      switch (mode) {
         case header:
            output << "TrackID ";
            break;
         case unit:
            output << "    (-) ";
            break;
         default:
            output << setw(7) << tupel.TrackID << " ";
      }
   }
}

void SEMGeneralDetector::OutputPID(ofstream &output,SEMTupleElement tupel,OutputType mode) {
   if (fSelectPID) {
      switch (mode) {
         case header:
            output << "ParentID ";
            break;
         case unit:
            output << "     (-) ";
            break;
         default:
            output << setw(8) << tupel.ParentID << " ";
      }
   }
}

void SEMGeneralDetector::OutputEtot(ofstream &output,SEMTupleElement tupel,OutputType mode,G4int accuracy) {
   if (fSelectEtot) {
   switch (mode) {
      case header:
         output << setw(7+accuracy) << "Etotal" << " ";
         break;
      case unit:
         output << setw(7+accuracy) << "(eV)" << " ";
         break;
      default:
         output << scientific << setprecision(accuracy) << setw(7+accuracy) << tupel.toten/eV << " ";
   }
   }
}

void SEMGeneralDetector::OutputEkin(ofstream &output,SEMTupleElement tupel,OutputType mode,G4int accuracy) {
   if (fSelectEkin) {
   switch (mode) {
      case header:
         output << setw(7+accuracy) << "Ekin" << " ";
         break;
      case unit:
         output << setw(7+accuracy) << "(eV)" << " ";
         break;
      default:
         output << scientific << setprecision(accuracy) << setw(7+accuracy) << tupel.kinen/eV << " ";
   }
   }
}

void SEMGeneralDetector::OutputEdeposit(ofstream &output,SEMTupleElement tupel,OutputType mode,G4int accuracy) {
   if (fSelectEdeposit) {
   if (accuracy < 2) accuracy=2; // Make the header fit
   switch (mode) {
      case header:
         output << setw(7+accuracy) << "Edeposit" << " ";
         break;
      case unit:
         output << setw(7+accuracy) << "(eV)" << " ";
         break;
      default:
         output << scientific << setprecision(accuracy) << setw(7+accuracy) << tupel.depositen/eV << " ";
   }
   }
}

void SEMGeneralDetector::OutputEloss(ofstream &output,SEMTupleElement tupel,OutputType mode,G4int accuracy) {
   if (fSelectEloss) {
   if (accuracy < 2) accuracy=2; // Make the header fit
   static G4RunManager* run = G4RunManager::GetRunManager();
   SEMPrimaryGeneratorAction* primgen = (SEMPrimaryGeneratorAction*) run->GetUserPrimaryGeneratorAction();
   G4double GunEnergy = primgen->GetParticleEnergy();
   switch (mode) {
      case header:
         output << setw(7+accuracy) << "Eloss" << " ";
         break;
      case unit:
         output << setw(7+accuracy) << "(eV)" << " ";
         break;
      default:
         output << scientific << setprecision(accuracy) << setw(7+accuracy) << (GunEnergy-tupel.kinen)/eV << " ";
   }
   }
}

void SEMGeneralDetector::OutputTime(ofstream &output,SEMTupleElement tupel,OutputType mode,G4int accuracy) {
   if (fSelectTime) {
   if (accuracy < 2) accuracy=2; // Make the header fit
   switch (mode) {
      case header:
         output << setw(7+accuracy) << "Hit time" << " ";
         break;
      case unit:
         output << setw(7+accuracy) << "(ns)" << " ";
         break;
      default:
         output << scientific << setprecision(accuracy) << setw(7+accuracy) << tupel.Time/ns << " ";
   }
   }
}

void SEMGeneralDetector::OutputPosition(ofstream &output,SEMTupleElement tupel,OutputType mode,G4int accuracy) {
   if (fSelectPosition) {
   switch (mode) {
      case header:
         output << setw(7+accuracy) << "x" << " ";
         output << setw(7+accuracy) << "y" << " ";
         output << setw(7+accuracy) << "z" << " ";
         break;
      case unit:
         output << setw(7+accuracy) << "(nm)" << " ";
         output << setw(7+accuracy) << "(nm)" << " ";
         output << setw(7+accuracy) << "(nm)" << " ";
         break;
      default:
         output << scientific << setprecision(accuracy) << setw(7+accuracy) << tupel.worldposition.x()/nanometer << " ";
         output << scientific << setprecision(accuracy) << setw(7+accuracy) << tupel.worldposition.y()/nanometer << " ";
         output << scientific << setprecision(accuracy) << setw(7+accuracy) << tupel.worldposition.z()/nanometer << " ";
   }
   }
}

void SEMGeneralDetector::OutputVertexPosition(ofstream &output,SEMTupleElement tupel,OutputType mode,G4int accuracy) {
   if (fSelectVertexPosition) {
   switch (mode) {
      case header:
         output << setw(7+accuracy) << "xOrigin" << " ";
         output << setw(7+accuracy) << "yOrigin" << " ";
         output << setw(7+accuracy) << "zOrigin" << " ";
         break;
      case unit:
         output << setw(7+accuracy) << "(nm)" << " ";
         output << setw(7+accuracy) << "(nm)" << " ";
         output << setw(7+accuracy) << "(nm)" << " ";
         break;
      default:
         output << scientific << setprecision(accuracy) << setw(7+accuracy) << tupel.vertexposition.x()/nanometer << " ";
         output << scientific << setprecision(accuracy) << setw(7+accuracy) << tupel.vertexposition.y()/nanometer << " ";
         output << scientific << setprecision(accuracy) << setw(7+accuracy) << tupel.vertexposition.z()/nanometer << " ";
   }
   }
}



void SEMGeneralDetector::OutputDirection(ofstream &output,SEMTupleElement tupel,OutputType mode,G4int accuracy) {
   if (fSelectDirection) {
   if (accuracy < 4) accuracy=4; // Make the header fit
   switch (mode) {
      case header:
         output << setw(7+accuracy) << "direction-x" << " ";
         output << setw(7+accuracy) << "direction-y" << " ";
         output << setw(7+accuracy) << "direction-z" << " ";
         break;
      case unit:
         output << setw(7+accuracy) << "(-)" << " ";
         output << setw(7+accuracy) << "(-)" << " ";
         output << setw(7+accuracy) << "(-)" << " ";
         break;
      default:
         output << scientific << setprecision(accuracy) << setw(7+accuracy) << tupel.momentumdirection.x() << " ";
         output << scientific << setprecision(accuracy) << setw(7+accuracy) << tupel.momentumdirection.y() << " ";
         output << scientific << setprecision(accuracy) << setw(7+accuracy) << tupel.momentumdirection.z() << " ";
   }
   }
}

void SEMGeneralDetector::OutputDepth(ofstream &output,SEMTupleElement tupel,OutputType mode,G4int accuracy) {
   if (fSelectDepth) {
   if (accuracy < 2) accuracy=2; // Make the header fit
   switch (mode) {
      case header:
         output << setw(7+accuracy) << "Max Depth" << " ";
         break;
      case unit:
         output << setw(7+accuracy) << "(nm)" << " ";
         break;
      default:
         output << scientific << setprecision(accuracy) << setw(7+accuracy) << tupel.Depth/nanometer << " ";
   }
   }
}

void SEMGeneralDetector::OutputMaxRadius(ofstream &output,SEMTupleElement tupel,OutputType mode,G4int accuracy) {
   if (fSelectMaxRadius) {
   if (accuracy < 3) accuracy=3; // Make the header fit
   switch (mode) {
      case header:
         output << setw(7+accuracy) << "Max Radius" << " ";
         break;
      case unit:
         output << setw(7+accuracy) << "(nm)" << " ";
         break;
      default:
         output << scientific << setprecision(accuracy) << setw(7+accuracy) << sqrt(tupel.Rho)/nanometer << " ";
   }
   }
}

void SEMGeneralDetector::OutputScatteringAngle(ofstream &output,SEMTupleElement tupel,OutputType mode,G4int accuracy) {
   if (fSelectScatteringAngle) {
   if (accuracy < 2) accuracy=2; // Make the header fit
   switch (mode) {
      case header:
         output << setw(7+accuracy) << "Sc.Angle" << " ";
         break;
      case unit:
         output << setw(7+accuracy) << "(rad)" << " ";
         break;
      default:
         output << scientific << setprecision(accuracy) << setw(7+accuracy) << acos(fabs(tupel.momentumdirection.z())) << " ";
   }
   }
}

void SEMGeneralDetector::OutputLogicalVolume(ofstream &output,SEMTupleElement tupel,OutputType mode) {
   if (fSelectLogicalVolume) {
   switch (mode) {
      case header:
         output <<setw(15) << "Logical Volume" << " ";
         break;
      case unit:
         output << setw(15) << " " << " ";
         break;
      default:
         output << " <" << tupel.lv << "> ";
   }
   }
}

void SEMGeneralDetector::OutputProcess(ofstream &output,SEMTupleElement tupel,OutputType mode) {
   if (fSelectProcess) {
   switch (mode) {
      case header:
         output <<setw(15) << "Process" << " ";
         break;
      case unit:
         output << setw(15) << " " << " ";
         break;
      default:
         output << " <" << tupel.process << "> ";
   }
   }
}

void SEMGeneralDetector::OutputCategory(ofstream &output,SEMTupleElement tupel,OutputType mode) {
   if (fSelectCategory) {
   switch (mode) {
      case header:
         output <<setw(15) << "Category" << " ";
         break;
      case unit:
         output << setw(15) << " " << " ";
         break;
      default:
         output << " <" << tupel.category << "> ";
   }
   }
}

void SEMGeneralDetector::OutputPotential(ofstream &output,G4double V,OutputType mode,G4int accuracy) {
   if (fSelectPotential) {
   switch (mode) {
      case header:
         output << setw(7+accuracy) << "V" << " ";
         break;
      case unit:
         output << setw(7+accuracy) << "(V)" << " ";
         break;
      default:
         output << scientific << setprecision(accuracy) << setw(7+accuracy) << V/volt << " ";
   }
   }
}



void SEMGeneralDetector::OutputHits(string filename,SEMTuple<SEMTupleElement>& tuple)
{
   if (filename!="") {
      ofstream output;
      OpenOutputFile(filename,output);
      if (!output.is_open()) {
         G4cerr << "File " << filename << " could not be opened!" << G4endl;
      } else {
      G4int length = tuple.Count();
      if (length>0) {
         SEMTupleElement tempel2;

         static G4RunManager* run = G4RunManager::GetRunManager();
         SEMDetectorConstruction* detectorConstruction = (SEMDetectorConstruction*) run->GetUserDetectorConstruction();

         OutputEID(output,tempel2,header); 
         OutputTID(output,tempel2,header); 
         OutputPID(output,tempel2,header); 
         OutputEtot(output,tempel2,header,6);
         OutputEkin(output,tempel2,header,6);
         OutputEdeposit(output,tempel2,header,6);
         OutputEloss(output,tempel2,header,6);
         OutputTime(output,tempel2,header,3);
         OutputPosition(output,tempel2,header,6);
         OutputVertexPosition(output,tempel2,header,3);
         OutputDirection(output,tempel2,header,6);
         OutputDepth(output,tempel2,header,6);
         OutputMaxRadius(output,tempel2,header,6);
         OutputScatteringAngle(output,tempel2,header,6);
         OutputLogicalVolume(output,tempel2,header);
         OutputProcess(output,tempel2,header);
         OutputCategory(output,tempel2,header);
         output << G4endl;
         OutputEID(output,tempel2,unit); 
         OutputTID(output,tempel2,unit); 
         OutputPID(output,tempel2,unit); 
         OutputEtot(output,tempel2,unit,6);
         OutputEkin(output,tempel2,unit,6);
         OutputEdeposit(output,tempel2,unit,6);
         OutputEloss(output,tempel2,unit,6);
         OutputTime(output,tempel2,unit,3);
         OutputPosition(output,tempel2,unit,6);
         OutputVertexPosition(output,tempel2,unit,3);
         OutputDirection(output,tempel2,unit,6);
         OutputDepth(output,tempel2,unit,6);
         OutputMaxRadius(output,tempel2,unit,6);
         OutputScatteringAngle(output,tempel2,unit,6);
         output << G4endl;
         for (int i=0;i<length;i++) {
            tempel2 = tuple[i];
            if ((!Ewindow) || ((tempel2.kinen >= Emin) && (tempel2.kinen <= Emax))) { // Note: boundaries of interval are INCLUSIVE
               OutputEID(output,tempel2,data);
               OutputTID(output,tempel2,data);
               OutputPID(output,tempel2,data);
               OutputEtot(output,tempel2,data,6);
               OutputEkin(output,tempel2,data,6);
               OutputEdeposit(output,tempel2,data,6);
               OutputEloss(output,tempel2,data,6);
               OutputTime(output,tempel2,data,3);
               OutputPosition(output,tempel2,data,6);
               OutputVertexPosition(output,tempel2,data,3);
               OutputDirection(output,tempel2,data,6);
               OutputDepth(output,tempel2,data,6);
               OutputMaxRadius(output,tempel2,data,6);
               OutputScatteringAngle(output,tempel2,data,6);
               OutputLogicalVolume(output,tempel2,data);
               OutputProcess(output,tempel2,data);
               OutputCategory(output,tempel2,data);
               output << G4endl;
            }
         }
      }
      output.close();
      G4cout << "Detector hits written to file " << filename << G4endl;
      }
   } else {
      G4cout << "No filename given, no output generated" << G4endl;
   }
}

void SEMGeneralDetector::OutputDeposits(string filename,SEMTuple<SEMTupleElement>& tuple,G4double mindeposit)
{
   if (filename!="") {
      ofstream output;
      OpenOutputFile(filename,output);
      if (!output.is_open()) {
         G4cerr << "File " << filename << " could not be opened!" << G4endl;
      } else {
         G4int length = tuple.Count();
         G4bool tmpSelectEID=fSelectEID;
         G4bool tmpSelectTID=fSelectTID;
         G4bool tmpSelectEdeposit=fSelectEdeposit;
         G4bool tmpSelectPosition=fSelectPosition;
         fSelectEID=true;
         fSelectTID=true;
         fSelectEdeposit=true;
         fSelectPosition=true;
         if (length>0) {
            SEMTupleElement tempel2;
            OutputEID(output,tempel2,header); 
            OutputTID(output,tempel2,header); 
            OutputEdeposit(output,tempel2,header,6); 
            OutputPosition(output,tempel2,header,6); 
            output << G4endl;
            OutputEID(output,tempel2,unit); 
            OutputTID(output,tempel2,unit); 
            OutputEdeposit(output,tempel2,unit,6); 
            OutputPosition(output,tempel2,unit,6); 
            output << G4endl;
            for (int i=0;i<length;i++) {
               tempel2 = tuple[i];
               if (fabs(tempel2.depositen) > mindeposit) {
                  OutputEID(output,tempel2,data); 
                  OutputTID(output,tempel2,data); 
                  OutputEdeposit(output,tempel2,data,6); 
                  OutputPosition(output,tempel2,data,6); 
                  output << G4endl;
               }
            }
         }
         fSelectEID=tmpSelectEID;
         fSelectTID=tmpSelectTID;
         fSelectEdeposit=tmpSelectEdeposit;
         fSelectPosition=tmpSelectPosition;
         output.close();
         G4cout << "Detector deposits written to file " << filename << G4endl;
      }
   } else {
      G4cout << "No filename given, no output generated" << G4endl;
   }
}

G4double SEMGeneralDetector::Gauss(G4double x,G4double sigma)
{
   return( exp(-1.0*x*x/(2.0*sigma*sigma))/(sigma*norm));
}

void SEMGeneralDetector::OutputEnergyPDF(string filename,SEMTuple<SEMTupleElement>& tuple,G4int PE)
{
   if (filename!="") {
      ofstream output;
      OpenOutputFile(filename,output);
      if (!output.is_open()) {
         G4cerr << "File " << filename << " could not be opened!" << G4endl;
      } else {
      output << setprecision(10);
      G4int length = tuple.Count();
      G4int minhits = 100;

     
      G4int  i,j;
      G4int    alength=0; 
      G4double Norm=0.0,inprod,momentum,angle;
      SEMTuple<SEMTupleElement> atuple;
      vector <G4int> index;
      index.resize(length);
      if (Awindow) Norm = sqrt(Adir.x()*Adir.x()+Adir.y()*Adir.y()+Adir.z()*Adir.z());
      for (i=0; i<length; i++) {
         if (Awindow) {
            inprod = tuple[i].momentumdirection.x()*Adir.x()+ 
                     tuple[i].momentumdirection.y()*Adir.y()+ 
                     tuple[i].momentumdirection.z()*Adir.z();
            momentum = sqrt(tuple[i].momentumdirection.x()*tuple[i].momentumdirection.x()+
                            tuple[i].momentumdirection.y()*tuple[i].momentumdirection.y()+
                            tuple[i].momentumdirection.z()*tuple[i].momentumdirection.z());
            angle = acos(inprod/(momentum*Norm));
            if ((angle>=Amin) && (angle <=Amax)) { // NOTE: Intervals include boundaries
               index[i] = i;
               alength ++;
            } else {
               index[i] = 0;
            }
         } else {
            index[i] = i;
         }
      }
      if (Awindow) {
         // Copy the hits that are in the right angle window to a new vector
         j=0;
         for (i=0; i<length; i++) {
            if (index[i] !=0) atuple.Add(tuple[i]);
         }
         length=alength;  // DANGER!!
      }

      if (length>minhits) {
         // Several ideas are implemented here at the same time:
         // - Use a Parzen window approach with a Gaussian as a kernel
         // - The standard deviation of the Gaussian is adapted to the local situation i.e. if the data points are very close
         //   we have a small sigma, if the data are far apart (in energy) sigma is large. A lower limit on sigma is used to
         //   to prevent unphysical behaviour (the physics is not accurate enough to allow resolving details << 1eV).
         // - The points at which the pdf is evaluated are spaced based on their distribution: more evaluations where the
         //   pdf has large values. Special care must be taken not to step over sharp features (abrupt change in sigma).

         vector <G4double> sigma;  // Store the local standard deviation at each position in the spectrum.
         G4int  nav = 10;
         G4double sum;
         sigma.resize(length);


         // Compute the local standard deviation at each hit, skip the first few at the beginning and the end of the list.
         for (i=nav; i<length-nav; i++) {
            sum=0.0;
            for (j=i-nav; j<= i+nav; j++) {
               if (Awindow) {
                  sum += atuple[j].kinen;
               } else {
                  sum += tuple[j].kinen;
               }
            }
            sum /= (2.0*nav+1.0); // This is the average energy around this hit
            sigma[i] = 0.0;
            for (j=i-nav; j<=i+nav; j++) {
               if (Awindow) {
                  sigma[i] += (atuple[j].kinen-sum)*(atuple[j].kinen-sum);
               } else {
                  sigma[i] += (tuple[j].kinen-sum)*(tuple[j].kinen-sum);
               }
            }
            sigma[i] = sqrt(sigma[i]/(2.0*nav));
            // We are not allowing too fine detail in the PDF (the physics is not that accurate anyway)
            if (sigma[i] < 0.2*eV) sigma[i] = 0.2*eV;
         }
         // Just extrapolate the sigmas at the begin and end of the list.
         for (i=0; i<nav; i++) {
            sigma[i] = sigma[nav];
            sigma[length-i-1] = sigma[length-nav-1];
         }

         // Just an attempt for an alternative approach...
//         G4double integral;
//         integral=0.0;
//         for (i=0; i<length-1; i++) {
//              integral += (tuple[i+1].kinen/eV-tuple[i].kinen/eV)*(1.0/sigma[i+1]+1.0/sigma[i])/2.0;
//         }
//         output << "Integral = " << integral << G4endl;
//         for (i=0; i<length; i++) {
//            output << tuple[i].kinen/eV << " " << 1.0/(sigma[i]*integral) << G4endl;
//         }

         // Write header to file
         output << "Kinetic energy probability density function (PDF)" << G4endl;
         output << "Lowest energy = " << tuple[0].kinen/eV << " eV " << G4endl;
         output << "Highest energy = " << tuple[length-1].kinen/eV << " eV " << G4endl;
         output << "Potential energy = " << (tuple[0].toten - tuple[0].kinen)/eV << " eV " << G4endl;
         output << "Number of primaries = " << PE << G4endl;
         output << "Number of hits = " << length << G4endl;
         if (Awindow) {
            output << "Reference axis for angular window : " << Adir << G4endl;
            output << "Angular window : " << Amin/deg << " < theta <= " << Amax/deg << " (deg) " << G4endl;
            output << "Energy (eV)       Probability (-)      Energy*Probability (eV)            Normalized on total hits " << G4endl;
         } else {
            output << "Energy (eV)       Probability (-)      Energy*Probability (eV)" << G4endl;
         }

         G4double E,S,f;
         G4bool finished = false;

         i=j=0;
         E=tuple[0].kinen;
         do {
            if (2.0*sigma[i]< 0.2*eV) {
               S = 0.2*eV;  // This is the minimal stepsize to prevent too much detail in the distribution
            } else {
               S = 2.0*sigma[i]; // This is the local sigma (determines size of next step)
            }

            // Add Gaussians with local sigmas to probability distribution at present energy E
            f = 0.0;
            for (G4int k=0; k<length; k++) {
               if (Awindow) {
                  if (fabs((atuple[k].kinen - E)/sigma[k])<10.0) {
                     f += Gauss(atuple[k].kinen - E,sigma[k]);
                  }
               } else {
                  if (fabs((tuple[k].kinen - E)/sigma[k])<10.0) {
                     f += Gauss(tuple[k].kinen - E,sigma[k]);
                  }
               }
            }

            output << scientific << setprecision(6) << setw(14) << E/eV << setw(14) << f*eV/length << setw(14) << f*E/length;
            if (Awindow) {
               output << scientific << setprecision(6) << setw(14) << f*eV/tuple.Count() << setw(14) << f*E/tuple.Count();
            } 
            output << G4endl;

            // First guess for the next evaluation point is E+S
            // Also check for sudden decrease in density of samples before accepting a next evaluation point
            // NOTE: Assumption: we start with E_i < E+S
            i++; 
            if (i>=length) {
               i=length-1;
               finished=true;
            }
            if (Awindow) {
               while ((atuple[i].kinen < (E+S)) && (atuple[i].kinen-sigma[i] < E ) && (!finished)) {
                  i++; 
                  if (i>=length) {
                     i=length-1;
                     finished=true;
                  }
               }
               E = atuple[i].kinen;
            } else {
               while ((tuple[i].kinen < (E+S)) && (tuple[i].kinen-sigma[i] < E ) && (!finished)) {
                  i++; 
                  if (i>=length) {
                     i=length-1;
                     finished=true;
                  }
               }
               E = tuple[i].kinen;
            }
         } while (!finished);

         G4cout << "Energy PDF written to file " << filename << " ( " << length << " hits)" << G4endl;
      } else {
         output << "No Energy PDF computed for less than " << minhits << " detector hits (only " << length << " hits). "<< G4endl;
         G4cout << "No Energy PDF computed for less than " << minhits << " detector hits (only " << length << " hits). "<< G4endl;
      }
      output.close();
      }
   }
}

void SEMGeneralDetector::OutputGammaPDF(string filename,SEMTuple<SEMTupleElement>& tuple,G4int PE)
{

   // Generate a spectrum for gamma photons only
   if (filename!="") {
      ofstream output;
      OpenOutputFile(filename,output);
      if (!output.is_open()) {
         G4cerr << "File " << filename << " could not be opened!" << G4endl;
      } else {
      output << setprecision(10);
      G4int length = tuple.Count();
      G4int minhits = 100;
     
      G4int  i,j;
      G4int    glength=0; 
      SEMTuple<SEMTupleElement> gtuple;
      vector <G4int> index;
      index.resize(length);
      for (i=0; i<length; i++) {
         if ( tuple[i].category=="gamma") {
            glength++;
            index[i]=i;
         } else {
            index[i]=0;
         }
      }
      // Copy the hits that are in the right angle window to a new vector
      j=0;
      for (i=0; i<length; i++) {
         if (index[i] !=0) gtuple.Add(tuple[i]);
      }
      length=glength;  // DANGER!!

      if (length>minhits) {
         // Several ideas are implemented here at the same time:
         // - Use a Parzen window approach with a Gaussian as a kernel
         // - The standard deviation of the Gaussian is adapted to the local situation i.e. if the data points are very close
         //   we have a small sigma, if the data are far apart (in energy) sigma is large. A lower limit on sigma is used to
         //   to prevent unphysical behaviour (the physics is not accurate enough to allow resolving details << 1eV).
         // - The points at which the pdf is evaluated are spaced based on their distribution: more evaluations where the
         //   pdf has large values. Special care must be taken not to step over sharp features (abrupt change in sigma).

         vector <G4double> sigma;  // Store the local standard deviation at each position in the spectrum.
         G4int  nav = 10;
         G4double sum;
         sigma.resize(length);


         // Compute the local standard deviation at each hit, skip the first few at the beginning and the end of the list.
         for (i=nav; i<length-nav; i++) {
            sum=0.0;
            for (j=i-nav; j<= i+nav; j++) {
               sum += gtuple[j].kinen;
            }
            sum /= (2.0*nav+1.0); // This is the average energy around this hit
            sigma[i] = 0.0;
            for (j=i-nav; j<=i+nav; j++) {
               sigma[i] += (gtuple[j].kinen-sum)*(gtuple[j].kinen-sum);
            }
            sigma[i] = sqrt(sigma[i]/(2.0*nav));
            // We are not allowing too fine detail in the PDF (the physics is not that accurate anyway)
            if (sigma[i] < 0.2*eV) sigma[i] = 0.2*eV;
         }
         // Just extrapolate the sigmas at the begin and end of the list.
         for (i=0; i<nav; i++) {
            sigma[i] = sigma[nav];
            sigma[length-i-1] = sigma[length-nav-1];
         }

         output << "Kinetic energy probability density function (PDF) for photons only" << G4endl;
         output << "Lowest energy = " << tuple[0].kinen/eV << " eV " << G4endl;
         output << "Highest energy = " << tuple[length-1].kinen/eV << " eV " << G4endl;
         output << "Potential energy = " << (tuple[0].toten - tuple[0].kinen)/eV << " eV " << G4endl;
         output << "Number of primaries = " << PE << G4endl;
         output << "Number of hits = " << length << G4endl;
         output << "Energy (eV)       Probability (-)      Energy*Probability (eV)" << G4endl;

         G4double E,S,f;
         G4bool finished = false;

         i=j=0;
         E=gtuple[0].kinen;
         do {
            if (2.0*sigma[i]< 0.2*eV) {
               S = 0.2*eV;  // This is the minimal stepsize to prevent too much detail in the distribution
            } else {
               S = 2.0*sigma[i]; // This is the local sigma (determines size of next step)
            }

            // Add Gaussians with local sigmas to probability distribution at present energy E
            f = 0.0;
            for (G4int k=0; k<length; k++) {
               if (fabs((gtuple[k].kinen - E)/sigma[k])<10.0) {
                     f += Gauss(gtuple[k].kinen - E,sigma[k]);
               }
            }

            output << scientific << setprecision(6) << setw(14) << E/eV << setw(14) << f*eV/length << setw(14) << f*E/length;
            output << G4endl;

            // First guess for the next evaluation point is E+S
            // Also check for sudden decrease in density of samples before accepting a next evaluation point
            // NOTE: Assumption: we start with E_i < E+S
            i++; 
            if (i>=length) {
               i=length-1;
               finished=true;
            }
            while ((gtuple[i].kinen < (E+S)) && (gtuple[i].kinen-sigma[i] < E ) && (!finished)) {
               i++; 
               if (i>=length) {
                  i=length-1;
                  finished=true;
               }
            }
            E = gtuple[i].kinen;
         } while (!finished);

         G4cout << "Gamma PDF written to file " << filename << " ( " << length << " hits)" << G4endl;
      } else {
         output << "No Gamma PDF computed for less than " << minhits << " detector hits (only " << length << " hits). "<< G4endl;
         G4cout << "No Gamma PDF computed for less than " << minhits << " detector hits (only " << length << " hits). "<< G4endl;
      }
      output.close();
      }
   }
}

void SEMGeneralDetector::OutputEnergyHistogram(string filename,SEMTuple<SEMTupleElement>& tuple,G4double BinSize,G4int PE)
{
   if (filename!="") {
      ofstream output;
      OpenOutputFile(filename,output);
      if (!output.is_open()) {
         G4cerr << "File " << filename << " could not be opened!" << G4endl;
      } else {
      output << setprecision(10);
      G4int length = tuple.Count();
      G4double Norm=0.0,inprod,momentum,angle;
      if (length>0) {
         G4double            min,max;
         SEMTupleElement tempel2;
         max=min=0;
         for (int i=0;i<length;i++) {
            tempel2 = tuple[i];
            if (i==0) {
               max=min=tempel2.kinen;
            }
            if (tempel2.kinen > max) max=tempel2.kinen;
            if (tempel2.kinen < min) min=tempel2.kinen;
         }
         output << "Kinetic energy histogram" << G4endl;
         output << "Lowest energy = " << min/eV << " eV " << G4endl;
         output << "Highest energy = " << max/eV << " eV " << G4endl;
         output << "Potential energy = " << (tuple[0].toten - tuple[0].kinen)/eV << " eV " << G4endl;
         output << "Number of primaries = " << PE << G4endl;
         output << "Number of hits = " << length << G4endl;
         if (Awindow) {
             output << "Reference axis for angular window : " << Adir << G4endl;
             output << "Angular window : " << Amin/deg << " < theta <= " << Amax/deg << " (deg) " << G4endl;
             Norm = sqrt(Adir.x()*Adir.x()+Adir.y()*Adir.y()+Adir.z()*Adir.z());
         }
         if (min<max) {
            G4int  nSE,nBSE,nSE3;
            G4int  nbin;
            vector <double> Histogram;
            nbin = (int) ((max-min)/BinSize);
            nSE=nBSE=nSE3=0;
            Histogram.resize(nbin+1);
            for (int i=0; i<length; i++) {
               G4double Ekin;
               tempel2 = tuple[i];
               Ekin = tempel2.kinen;
               if (tempel2.category=="SE3") {
                  nSE3++;
               } else if (Ekin<=50.0*eV) {
                  nSE++;
               } else {
                  nBSE++;
               }
               if (Awindow) {
                  inprod = tuple[i].momentumdirection.x()*Adir.x()+ 
                           tuple[i].momentumdirection.y()*Adir.y()+ 
                           tuple[i].momentumdirection.z()*Adir.z();
                  momentum = sqrt(tuple[i].momentumdirection.x()*tuple[i].momentumdirection.x()+
                               tuple[i].momentumdirection.y()*tuple[i].momentumdirection.y()+
                               tuple[i].momentumdirection.z()*tuple[i].momentumdirection.z());
                  angle = acos(inprod/(momentum*Norm));
                  if ((angle>=Amin) && (angle <=Amax)) { // NOTE: Intervals include boundaries
                     Histogram[(int) ((Ekin-min)/BinSize)] ++;
                  }
               } else {
                  Histogram[(int) ((Ekin-min)/BinSize)] ++;
               }
            }
            output << "Number of SEs = " << nSE << G4endl;
            output << "Number of BSEs = " << nBSE << G4endl;
            if (nSE3!=0) output << "Number of SE3s = " << nSE3 << G4endl;
            G4double cumul=0.0;
            for (int i=0; i<=nbin; i++) {
               G4double E;
               E = ((double)i+0.5)*BinSize+min;
               cumul+=Histogram[i];
               output << i << " " << E/eV << " " << Histogram[i] << " " << cumul << G4endl;
            }
         }
      }
      output.close();
      G4cout << "Energy histogram written to file " << filename << G4endl;
      }
   }
}

void SEMGeneralDetector::OutputRhoHistogram(string filename,SEMTuple<SEMTupleElement>& tuple,G4double BinSize,G4int PE)
{
   if (filename!="") {
      ofstream output;
      OpenOutputFile(filename,output);
      if (!output.is_open()) {
         G4cerr << "File " << filename << " could not be opened!" << G4endl;
      } else {
      output << setprecision(10);
      G4int length = tuple.Count();
      G4double Norm=0.0,inprod,momentum,angle;
      G4double pi = 4.0*atan(1.0);
      if (length>0) {
         G4double            min,max,rho;
         max=min=0.0;
         for (int i=0;i<length;i++) {
            rho = tuple[i].localposition.rho();
            if (i==0) {
               max=min=rho;
            }
            if (rho > max) max=rho;
            if (rho < min) min=rho;
         }
         output << "Radius histogram" << G4endl;
         output << "Smallest radius = " << min/mm << " mm (histogram will start at r=0)" << G4endl;
         output << "Largest radius  = " << max/mm << " mm " << G4endl;
         output << "Number of primaries = " << PE << G4endl;
         output << "Number of hits = " << length << G4endl;
         if (Awindow) {
             output << "Reference axis for angular window : " << Adir << G4endl;
             output << "Angular window : " << Amin/deg << " < theta <= " << Amax/deg << " (deg) " << G4endl;
             Norm = sqrt(Adir.x()*Adir.x()+Adir.y()*Adir.y()+Adir.z()*Adir.z());
         }
         if (Ewindow) {
             output << "Kinetic energy window : " << Emin/eV << " eV <= Ekin <= " << Emax/eV << " eV " << G4endl;
         }
         min=0.0;  // Start histogram at r=0
         G4int  nSE,nBSE,nSE3;
         G4int  nbin;
         vector <double> Histogram;
         nbin = (int) ((max-min)/BinSize);
         nSE=nBSE=nSE3=0;
         Histogram.resize(nbin+1);
         for (int i=0; i<length; i++) {
            rho = tuple[i].localposition.rho();
            if (Awindow || Ewindow) { // Energy and/or angle selection are active
               G4bool selected=false;
               if (Awindow) {
                  inprod = tuple[i].momentumdirection.x()*Adir.x()+ 
                           tuple[i].momentumdirection.y()*Adir.y()+ 
                           tuple[i].momentumdirection.z()*Adir.z();
                  momentum = sqrt(tuple[i].momentumdirection.x()*tuple[i].momentumdirection.x()+
                               tuple[i].momentumdirection.y()*tuple[i].momentumdirection.y()+
                               tuple[i].momentumdirection.z()*tuple[i].momentumdirection.z());
                  angle = acos(inprod/(momentum*Norm));
                  if ((angle>=Amin) && (angle <=Amax)) selected=true; // NOTE: Intervals include boundaries
               } 
               if (Ewindow) {
                  if ((tuple[i].kinen >= Emin) && (tuple[i].kinen <= Emax)) selected=true; // Note: boundaries of interval are INCLUSIVE
               }
               if (selected) {
                  Histogram[(int) ((rho-min)/BinSize)] ++;
               }
            } else { // No energy or angle selection active
               Histogram[(int) ((rho-min)/BinSize)] ++;
            }
         }
         output << " R (nm)    Number (-)     Number/Area (1/nm^2) " << G4endl;
         for (int i=0; i<=nbin; i++) {
            G4double R0,R1,R,Area;
            R0 = ((double)i*BinSize+min)/nanometer;
            R1 = (((double)i+1.0)*BinSize+min)/nanometer;
            R = sqrt((R1*R1+R0*R0)/2.0);
            Area = pi*(R1*R1-R0*R0);
            output << R << " " << Histogram[i] << " " << Histogram[i]/Area << G4endl;
         }
      }
      output.close();
      G4cout << "Radius histogram written to file " << filename << G4endl;
      }
   }
}

void SEMGeneralDetector::OutputRhoPDF(string filename,SEMTuple<SEMTupleElement>& tuple,G4int PE)
{
   if (filename!="") {
      ofstream output;
      OpenOutputFile(filename,output);
      if (!output.is_open()) {
         G4cerr << "File " << filename << " could not be opened!" << G4endl;
      } else {
      output << setprecision(10);
      G4int length = tuple.Count();
      G4int minhits = 100;

      if (length>minhits) {
         vector <G4double> sigma;  // Store the local standard deviation at each position in the histogram.
         G4int  nav = 10;
         G4int  i,j;
         G4double sum;
         sigma.resize(length);

         // Compute the local standard deviation at each hit, skip the first few at the beginning and the end of the list.
         for (i=nav; i<length-nav; i++) {
            sum=0.0;
            for (j=i-nav; j<= i+nav; j++) {
               sum += tuple[j].localposition.rho();
            }
            sum /= (2.0*nav+1.0); // This is the average radius around this hit
            sigma[i] = 0.0;
            for (j=i-nav; j<=i+nav; j++) {
               sigma[i] += (tuple[j].localposition.rho()-sum)*(tuple[j].localposition.rho()-sum);
            }
            sigma[i] = sqrt(sigma[i]/(2.0*nav));
            if (sigma[i] < tuple[i].localposition.rho()/20.0) {
                  sigma[i] = tuple[i].localposition.rho()/20.0;
            }
            if (sigma[i] == 0.0) {
                  sigma[i] = 0.01*nanometer;
            }
         }
         // Just extrapolate the sigmas at the begin and end of the list.
         for (i=0; i<nav; i++) {
            sigma[i] = sigma[nav];
            sigma[length-i-1] = sigma[length-nav-1];
         }

         // Write header to file
         output << "Radius histogram" << G4endl;
         output << "Smallest radius = " << tuple[0].localposition.rho()/nanometer << G4endl;
         output << "Largest radius = " << tuple[length-1].localposition.rho()/nanometer << G4endl;
         output << "Number of primaries = " << PE << G4endl;
         output << "Number of hits = " << length << G4endl;
         output << "Radius (nm)       Probability (1/nm)      " << G4endl;

         // Add Gaussians with local sigmas to probability distribution
         i=j=0;
         G4double R,S,f;

         G4double s0,s1,R0,f0;
         f0=0.0;
         R0=tuple[0].localposition.rho();
         s0=s1=0.0;
         R=tuple[0].localposition.rho();
         do {
            S = sigma[i]; // This is the local sigma (determines size of next step)

            f = 0.0;
            for (G4int k=0; k<length; k++) {
               if (fabs((tuple[k].localposition.rho() - R)/sigma[k])<10.0) {
                  f += Gauss(tuple[k].localposition.rho() - R,sigma[k]); 
               }
            }
            s0 += f*(R-R0)/length;
            s1 += f0*(R-R0)/length;
            R0=R;
            f0=f;

            output << R/nanometer << " " << f*nanometer/length << G4endl;

            R += S;
            G4bool finished=false;
            if (tuple[i].localposition.rho() >=R) i--;
            j=0;
            while (!finished) {
               finished = (tuple[i].localposition.rho() >= R);
               i++;

               // Limit the number of hits to be skipped in the list.
               // This deals with sudden increases in the probability (e.g. elastic peak).
               j++; if (j>sqrt((G4double)length)/2.0) {R = tuple[i].localposition.rho(); finished=true;}

               // Set up for last entry in list
               if (i==length-1) {R = tuple[i].localposition.rho(); finished=true;}

               // Signal exit of loop
               if (i==length) {R = -1.0; finished=true;}
            }
         } while (R>=0.0);
         G4cout << "Radial histogram written to file " << filename << G4endl;
      } else {
         output << "No Radial PDF computed for less than " << minhits << " detector hits." << G4endl;
         G4cout << "No Radial PDF computed for less than " << minhits << " detector hits." << G4endl;
      }
      output.close();
      }
   }
}

void SEMGeneralDetector::OutputAngleHistogram(string filename,SEMTuple<SEMTupleElement>& tuple,G4int PE, G4ThreeVector normal, G4int nbin) 
{
   // Angular histogram
   // Compute a histogram of angles w.r.t. the normal given
   // Kinetic energy is in given interval, boudaries of interval are included in the interval. If Emin < 0 no energy selection is done (which is the default)

   if (filename!="") {
      ofstream output;
      OpenOutputFile(filename,output);
      if (!output.is_open()) {
         G4cerr << "File " << filename << " could not be opened!" << G4endl;
      } else {
      output << setprecision(10);
      G4int length = tuple.Count();

      vector <G4double> AngHist;
      AngHist.resize(nbin);
      G4double rad2deg = 45.0/atan(1.0);
      G4double twopi = 8.0*atan(1.0);

      vector <G4double> theta;   
      theta.resize(nbin+1);
      theta[0] = 0.0;
      for (G4int i=1; i<nbin; i++) {
         theta[i] = acos(cos(theta[i-1])-2.0/(double)nbin);
      }
      theta[nbin] = twopi/2.0;

      G4double solidangle,thetacenter;
      G4double inprod,momentum,Norm,angle,kinen;
      G4int    ii,inwindow;

      Norm = sqrt(normal.x()*normal.x()+normal.y()*normal.y()+normal.z()*normal.z());
      inwindow=0;
      for (G4int i=0; i<length; i++) {
         kinen = tuple[i].kinen;
         if ((!Ewindow) || ((kinen >= Emin) && (kinen <= Emax))) { // Note: boundaries of interval are INCLUSIVE
            inprod = tuple[i].momentumdirection.x()*normal.x()+ 
                     tuple[i].momentumdirection.y()*normal.y()+ 
                     tuple[i].momentumdirection.z()*normal.z();
            momentum = sqrt(tuple[i].momentumdirection.x()*tuple[i].momentumdirection.x()+
                            tuple[i].momentumdirection.y()*tuple[i].momentumdirection.y()+
                            tuple[i].momentumdirection.z()*tuple[i].momentumdirection.z());
            angle = acos(inprod/(momentum*Norm));
            ii=0;
            while(angle > theta[ii+1]) ii++;
            AngHist[ii]++;
            inwindow++;
         }
      }

      output << "Angular distribution" << G4endl;
      output << "Number of primaries = " << PE << G4endl;
      output << "Number of hits = " << length << G4endl;
      if (Ewindow) output << "Number of hits in energy window from " << Emin/eV << " to " << Emax/eV << " = " << inwindow << G4endl;
      output << "Direction of normal = " << normal.x()/Norm << " " << normal.y()/Norm << " " << normal.z()/Norm << G4endl;
      output << "Angle (deg)    Number (-)   Fraction (-)    Flux (1/sr)    Angle (deg) Solid angle (sr)    Integral (-)    Fraction (-) " << G4endl;

      // The meaning of these numbers:
      // theta[i+1] : end of the angle bin
      // AngHist[i] : number of counts in the interval theta[i]..theta[i+1]
      // Fraction : fraction of the total number of counts in this interval
      // Flux : particle flux (counts/ solid angle)
      // Solid angle: total solid angle integrated from 0.. theta[i+1]
      // Sum: number of counts in the interval 0.. theta[i+1]
      // Fraction: fraction of the total number of counts in this interval
      G4double sum=0.0;
      for (G4int i=0; i<nbin; i++) {
         solidangle = twopi*(cos(theta[i])-cos(theta[i+1]));
         thetacenter = twopi*(sin(theta[i+1])-sin(theta[i])+theta[i]*cos(theta[i])-theta[i+1]*cos(theta[i+1]))/solidangle;
         sum += AngHist[i];
         output << thetacenter*rad2deg << " "
                << AngHist[i] << " " << AngHist[i]/inwindow << " " << AngHist[i]/solidangle << " " 
                << theta[i+1]*rad2deg << " " <<  twopi*(1.0-cos(theta[i+1])) << " " << sum << " " << sum/inwindow << G4endl;
      }

      output.close();
      G4cout << "Angular histogram written to file " << filename << G4endl;
      }
   }
}

void SEMGeneralDetector::OutputTimeHistogram(string filename,SEMTuple<SEMTupleElement>& tuple,G4double BinSize,G4int PE)
{
   if (filename!="") {
      ofstream output;
      OpenOutputFile(filename,output);
      if (!output.is_open()) {
         G4cerr << "File " << filename << " could not be opened!" << G4endl;
      } else {
      output << setprecision(10);
      G4int length = tuple.Count();
      if (length>0) {
         G4double            min,max;
         SEMTupleElement tempel2;
         max=min=0;
         for (int i=0;i<length;i++) {
            tempel2 = tuple[i];
            if (i==0) {
               max=min=tempel2.Time;
            }
            if (tempel2.Time > max) max=tempel2.Time;
            if (tempel2.Time < min) min=tempel2.Time;
         }
         output << "Time histogram" << G4endl;
         output << "Fastest time = " << min/ns << G4endl;
         output << "Slowest time = " << max/ns << G4endl;
         output << "Number of primaries = " << PE << G4endl;
         output << "Number of hits = " << length << G4endl;
         if (min<max) {
            G4int  nbin;
            vector <double> Histogram;
            nbin = (int) ((max-min)/BinSize);
            Histogram.resize(nbin+1);
            for (int i=0; i<length; i++) {
               G4double TOF;
               tempel2 = tuple[i];
               TOF = tempel2.Time;
               Histogram[(int) ((TOF-min)/BinSize)] ++;
            }
            for (int i=0; i<=nbin; i++) {
               G4double TOF;
               TOF = ((double)i+0.5)*BinSize+min;
               output << i << " " << TOF/ns << " " << Histogram[i] << G4endl;
            }
         }
      }
      G4cout << "Time histogram written to file " << filename << G4endl;

      // Arrival statistics...
      if (length != 0 ) {
         int EID = tuple[0].EventID;
         int n = 1;
         int nmax = 20; // Maximum number of arriving electrons for one event
         vector < int > H;
         H.resize(nmax);
         for (int i=0; i<nmax; i++) H[i] = 0;
         H[0] = PE;
         for (int i=0; i<length; i++) {
            if (EID != tuple[i].EventID) {  // Found next event, update histogram
               if (n<nmax) {
                  H[n] ++;
                  H[0] --;
               }
               EID = tuple[i].EventID;
               n=1;
            } else { // Still the same event, update counter
               n++;
            }
         }
         output << "Statistics of arrivals : " << G4endl;
         for (int i=0; i<nmax; i++) {
            output << i << " " << H[i] << G4endl;
         }
         G4cout << "Arrival statistics written to file " << filename << G4endl;
      }
      output.close();
      }
   }
}

void SEMGeneralDetector::OutputCounts(string filename,SEMTuple<SEMTupleElement>& tuple,G4int PE)
{
   if (filename!="") {
      ofstream output;
      static G4RunManager* run = G4RunManager::GetRunManager();
      if (filename != lastfilename) {
         SEMDetectorConstruction* detectorConstruction = (SEMDetectorConstruction*) run->GetUserDetectorConstruction();

         OpenOutputFile(filename,output);
         if (!output.is_open()) {
            G4cerr << "File " << filename << " could not be opened!" << G4endl;
         }
         detectorConstruction->PrintTree(output);
         // Column headers
         output << setw(13) << "Energy" 
                << setw(9)  << "n_PE" 
                << setw(9)  << "Total" 
                << setw(13) << "delta" 
                << setw(13) << "eta" 
                << setw(9)  << "n_SE3" ;
         if (Ewindow) {
            output << setw(9) << "n_Energy";
         }
         if (Awindow) {
            output << setw(9) << "n_Angle";
         }
         output << setw(17) << "E_deposit"
                << setw(17)  << "E_detected";
         if (fDetectorModel!=DT_Undefined) {
            output << setw(17) << "I_detector";
         }
         output << G4endl;

         // Units
         output << setw(13) << "(eV)" 
                << setw(9) << "(-)" 
                << setw(9) << "(-)" 
                << setw(13) << "(-)" 
                << setw(13) << "(-)" 
                << setw(9) << "(-)";
         if (Ewindow) {
            output << setw(9) << "(-)";
         }
         if (Awindow) {
            output << setw(9) << "(-)";
         }
         output << setw(17) << "(eV)"
                << setw(17) << "(eV)";
         if (fDetectorModel!=DT_Undefined) {
            output << setw(17) << "(arb.units)";
         }
         output << G4endl;
         if (Ewindow) {
            output << "Energy window: " << Emin/eV << " eV  to " << Emax/eV << " eV" << G4endl;
         }
         if (Awindow) {
            output << "Angle window: " << Amin/deg << " deg  to " << Amax/deg << " deg" << G4endl;
         }
         lastfilename=filename;
      } else {
         output.open (filename.c_str(),ios::app);
         if ((Ewindow) && (EwindowChanged)) {
            output << "Energy window: " << Emin/eV << " eV  to " << Emax/eV << " eV" << G4endl;
         }
         if ((Awindow) && (AwindowChanged)) {
            output << "Angle window: " << Amin/deg << " deg  to " << Amax/deg << " deg" << G4endl;
         }
      }
      EwindowChanged=false;
      AwindowChanged=false;

      SEMPrimaryGeneratorAction* gun = (SEMPrimaryGeneratorAction*) run->GetUserPrimaryGeneratorAction();

      output << scientific << setw(13) << setprecision(5) << gun->GetParticleEnergy()/eV << setw(9) << PE ;

      G4int length = tuple.Count();
      G4double totEkin = 0.0;
//      if (length>0) {
         SEMCountsPerCategory DetectorCounts;
         G4int  nSE,nBSE,nSE3,nEWindow,nAWindow;
         G4double Current;

         for (int i=0; i<length; i++) {
            totEkin += tuple[i].kinen;
         }

         DetectorCounts = GetDetectorCounts(tuple);
         nSE = DetectorCounts.nSE;
         nBSE = DetectorCounts.nBSE;
         nSE3 = DetectorCounts.nSE3;
         nEWindow = DetectorCounts.nEkinWindow;
         nAWindow = DetectorCounts.nAngWindow;
         Current = DetectorCounts.Current;
         output << setw(9) << nSE+nBSE+nSE3 << setw(13) << (double)nSE/(double)PE << setw(13) << (double)nBSE/(double)PE << setw(9) << nSE3;
         if (Ewindow) {
            output << setw(9) << nEWindow;
         }
         if (Awindow) {
            output << setw(9) << nAWindow;
         }
         output << setw(17) << setprecision(8) << GetEDeposit()/eV << setw(17) << totEkin/eV;
         if (fDetectorModel!=DT_Undefined) {
            output << setw(17) << Current;
         }
         output << G4endl;
//      } else {
//         output << "Nothing detected" << G4endl;
//      }
      G4cout << "Detector counts written to file " << filename << G4endl;
   } else {
      G4cout << "No filename given, no output generated " << G4endl;
   }
}

SEMCountsPerCategory SEMGeneralDetector::GetDetectorCounts(SEMTuple<SEMTupleElement>& tuple)
{
   SEMCountsPerCategory DetectorCounts;
   G4int length = tuple.Count();
   G4double totEkin = 0.0;
   if (length>0) {
      SEMTupleElement tempel2;
      G4double Norm=0;
      G4int  nSE,nBSE,nSE3,nEtotWindow,nEkinWindow,nAngWindow;
      nSE=nBSE=nSE3=nEtotWindow=nEkinWindow=nAngWindow=0;

      if (Awindow) Norm = sqrt(Adir.x()*Adir.x()+Adir.y()*Adir.y()+Adir.z()*Adir.z());
      for (int i=0; i<length; i++) {
         G4double Ekin,Etot,Angle,inprod,momentum;
         tempel2 = tuple[i];
         Ekin = tempel2.kinen;
         Etot = tempel2.toten;
         totEkin += Ekin;
         if (tempel2.category=="SE3") {
            nSE3++;
         } else if (Ekin<=50.0*eV) {
            nSE++;
         } else {
            nBSE++;
         }
         if ((Ewindow) && (Ekin>=Emin) && (Ekin <= Emax)) nEkinWindow++; // NOTE: Same window for both Etot and Ekin !!
         if ((Ewindow) && (Etot>=Emin) && (Etot <= Emax)) nEtotWindow++;
         if (Awindow) {
            inprod = tuple[i].momentumdirection.x()*Adir.x()+ 
                     tuple[i].momentumdirection.y()*Adir.y()+ 
                     tuple[i].momentumdirection.z()*Adir.z();
            momentum = sqrt(tuple[i].momentumdirection.x()*tuple[i].momentumdirection.x()+
                            tuple[i].momentumdirection.y()*tuple[i].momentumdirection.y()+
                            tuple[i].momentumdirection.z()*tuple[i].momentumdirection.z());
            Angle = acos(inprod/(momentum*Norm));
            if ((Angle>=Amin) && (Angle <=Amax)) nAngWindow++; // NOTE: Intervals include boundaries
         }
      }
      DetectorCounts.nTotal = length;
      DetectorCounts.nSE = nSE;
      DetectorCounts.nBSE = nBSE;
      DetectorCounts.nSE3 = nSE3;
      DetectorCounts.nEtotWindow = nEtotWindow;
      DetectorCounts.nEkinWindow = nEkinWindow;
      DetectorCounts.nAngWindow = nAngWindow;
      DetectorCounts.Current = GetDetectorCurrent(tuple); 
   }
   return(DetectorCounts);
}

G4double SEMGeneralDetector::GetDetectorCurrent(SEMTuple<SEMTupleElement>& tuple)
{
   G4double Current = 0.0;
   G4int length = tuple.Count();
   if (length>0) {
      switch (fDetectorModel) {
         case DT_SolidState:
            // Simulation of detector current for the Solid State (Si) detector
            // I_det = I_beam * E_beam/3.62eV *(1-eta)*R
            // Where
            // eta = backscatter coefficient of Si
            // R = detector response curve. Modeled here as: R = 1 - exp(-0.268*E_beam)
            // E_beam = energy of incoming beam
            // I_beam = current in incoming beam
            // Implementation:
            // Return N* = sum (i=1,N) E_kin(i)/3.62eV * (1-exp(-268*E_kin(i)))
            // with
            // N: the number of counts
            // E_kin(i): the kinetic energy of the i-th particle (in eV)
            for (int i=0; i<length; i++) {
               G4double Ekin;
               Ekin = tuple[i].kinen;
               Current += (Ekin/eV)*(1.0-exp(-0.268*Ekin/keV))/3.62;
            }
            break;
         default:
            // No model attached, don't know how to compute current so set it to total number of hits
            Current=(G4double)length;
      }
   }
   return(Current);
}

G4String SEMGeneralDetector::DirName(G4String source)
{
   source.erase(find(source.rbegin(), source.rend(), '/').base(), source.end());
   return source;
}


G4int SEMGeneralDetector::OpenOutputFile(G4String filename,ofstream &output) {

   // Open the file with name filename
   // If the directory does not exist, it will be created, unless it is more than one level deep
   // Return value: 
   // 0: Succes
   // 1: Directory could not be created
   // 2: File could not be opened

   G4String dirname = DirName(filename);
   struct stat sb;

   // Try creating the directory if it not already exists
   if (dirname!="") {
      if (!(stat(dirname.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))) {
         G4int status = mkdir(dirname.c_str(),S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH |S_IXOTH );
         if (status) {
            return 1;
         }
      }
   }
   // Try opening the file
   output.open(filename.c_str());
   if (output.fail()) {
      return 2;
   }
   return 0;
}
