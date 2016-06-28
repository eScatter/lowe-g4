#include "SEMSteppingVerbose.hh"

#include "G4SteppingManager.hh"
#include "G4UnitsTable.hh"

using namespace std;

SEMSteppingVerbose::SEMSteppingVerbose()
{
}

SEMSteppingVerbose::~SEMSteppingVerbose()
{}
 
void SEMSteppingVerbose::StepInfo()
{
  CopyState();
  
  G4int prec = G4cout.precision(12);

  if( verboseLevel >= 1 ){
    if( verboseLevel >= 4 ) VerboseTrack();
    if( verboseLevel >= 3 ){
      G4cout << G4endl;    
      G4cout << setw( 5) << "#Step#"     << " "
	     << setw(12) << "X"          << " "
	     << setw(12) << "Y"          << " "  
	     << setw(12) << "Z"          << " "
	     << setw(12) << "KineE"      << " "
	     << setw( 9) << "dEStep"     << " "  
	     << setw(10) << "StepLeng"     
	     << setw(10) << "TrakLeng" 
	     << setw(10) << "Volume"    << "  "
	     << setw(10) << "Process"   << G4endl;	          
    }

//	<< setw(8) << G4BestUnit(fTrack->GetPosition().x(),"Length")
//	<< setw(8) << G4BestUnit(fTrack->GetPosition().y(),"Length")
//	<< setw(8) << G4BestUnit(fTrack->GetPosition().z(),"Length")
//	<< setw(8) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")
    G4cout << setw(5) << fTrack->GetCurrentStepNumber() << " "
        << scientific << setprecision(7) 
	<< setw(14) << fTrack->GetPosition().x()/mm << " "
	<< setw(14) << fTrack->GetPosition().y()/mm << " "
	<< setw(14) << fTrack->GetPosition().z()/mm << " "
        << setprecision(5)
	<< setw(12) << fTrack->GetKineticEnergy()/eV << " "
        << setprecision(2)
	<< setw(9) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy") << " "
	<< setw(9) << G4BestUnit(fStep->GetStepLength(),"Length") << " "
	<< setw(9) << G4BestUnit(fTrack->GetTrackLength(),"Length")
	<< "  ";

    // if( fStepStatus != fWorldBoundary) 
    if( fTrack->GetNextVolume() != 0 ) { 
      G4cout << setw(10) << fTrack->GetVolume()->GetName();
    } else {
      G4cout << setw(10) << "OutOfWorld";
    }

    if(fStep->GetPostStepPoint()->GetProcessDefinedStep() != NULL){
      G4cout << "  " 
        << setw(10) << fStep->GetPostStepPoint()->GetProcessDefinedStep()
	                                ->GetProcessName();
    } else {
      G4cout << "   UserLimit";
    }

    G4cout << G4endl;

    if( verboseLevel == 2 ){
      G4int tN2ndariesTot = fN2ndariesAtRestDoIt +
	                    fN2ndariesAlongStepDoIt +
	                    fN2ndariesPostStepDoIt;
      if(tN2ndariesTot>0){
	G4cout << "    :----- List of 2ndaries - "
	       << "#SpawnInStep=" << setw(3) << tN2ndariesTot 
	       << "(Rest="  << setw(2) << fN2ndariesAtRestDoIt
	       << ",Along=" << setw(2) << fN2ndariesAlongStepDoIt
	       << ",Post="  << setw(2) << fN2ndariesPostStepDoIt
	       << "), "
	       << "#SpawnTotal=" << setw(3) << (*fSecondary).size()
	       << " ---------------"
	       << G4endl;

	for(size_t lp1=(*fSecondary).size()-tN2ndariesTot; 
                        lp1<(*fSecondary).size(); lp1++){
	  G4cout << "    : "
		 << setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().x(),"Length")
		 << setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().y(),"Length")
		 << setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().z(),"Length")
		 << setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(),"Energy")
		 << setw(10)
		 << (*fSecondary)[lp1]->GetDefinition()->GetParticleName();
	  G4cout << G4endl;
	}
              
	G4cout << "    :-----------------------------"
	       << "----------------------------------"
	       << "-- EndOf2ndaries Info ---------------"
	       << G4endl;
      }
    }
    
  }
  G4cout.precision(prec);
}

void SEMSteppingVerbose::TrackingStarted()
{

  CopyState();
G4int prec = G4cout.precision(12);
  if( verboseLevel > 0 ){

    G4cout << setw( 5) << "Step#"      << " "
           << setw(12) << "X"          << " "
	   << setw(12 ) << "Y"         << " "  
	   << setw(12) << "Z"          << " "
	   << setw(12) << "KineE"      << " "
	   << setw( 9) << "dEStep"     << " "  
	   << setw(10) << "StepLeng"  
	   << setw(10) << "TrakLeng"
	   << setw(10) << "Volume"     << "  "
	   << setw(10) << "Process"    << G4endl;	     

//	<< setw(8) << G4BestUnit(fTrack->GetPosition().x(),"Length")
//	<< setw(8) << G4BestUnit(fTrack->GetPosition().y(),"Length")
//	<< setw(8) << G4BestUnit(fTrack->GetPosition().z(),"Length")
    G4cout << setw(5) << fTrack->GetCurrentStepNumber() << " "
        << scientific << setprecision(7) 
	<< setw(14) << fTrack->GetPosition().x()/mm << " "
	<< setw(14) << fTrack->GetPosition().y()/mm << " "
	<< setw(14) << fTrack->GetPosition().z()/mm << " "
        << setprecision(5)
	<< setw(12) << fTrack->GetKineticEnergy()/eV << " "
        << setprecision(2)
	<< setw(9) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy") << " "
	<< setw(9) << G4BestUnit(fStep->GetStepLength(),"Length") << " "
	<< setw(9) << G4BestUnit(fTrack->GetTrackLength(),"Length")
	<< "  ";

    if (fTrack->GetCurrentStepNumber() >= 10000) G4Exception("User abort...","CurrentStepNumber>=10000",FatalException,"SEMSteppingVerbose(1)");

    if(fTrack->GetNextVolume()){
      G4cout << setw(10) << fTrack->GetVolume()->GetName();
    } else {
      G4cout << setw(10) << "OutOfWorld";
    }
    G4cout  << "    initStep" << G4endl;
  }
  G4cout.precision(prec);
}

