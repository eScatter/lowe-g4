// SEMVRML2FileViewer.hh, based on
// G4VRML2FileViewer.hh
// Satoshi Tanaka and Yasuhide Sawada

#ifndef SEMVRML2File_VIEWER_HH
#define SEMVRML2File_VIEWER_HH

#include <fstream>
#include "G4VViewer.hh"
#include "globals.hh"

class SEMVRML2FileSceneHandler;

class SEMVRML2FileViewer: public G4VViewer {
public:
	SEMVRML2FileViewer(SEMVRML2FileSceneHandler& scene, const G4String& name = "");
	virtual ~SEMVRML2FileViewer();
	void ClearView();
	void DrawView();
	void ShowView();
	void FinishView();
private:
	void SetView(); // Do nothing. SendViewParameters will do its job.
	void SendViewParameters ()  ;

private:
	SEMVRML2FileSceneHandler& fSceneHandler; // Reference to Graphics Scene for this view.
	std::ofstream&         fDest ;

	G4double    fViewHalfAngle ;	
	G4double    fsin_VHA       ;	

	G4double	InProduct(G4Vector3D,G4Vector3D);
	G4Vector3D	OutProduct(G4Vector3D,G4Vector3D);
	G4Vector3D	FindVectorAndRotation(G4Vector3D,G4Vector3D);//Calculates rotation vector and also updates Rotation
	G4double	Rotation;
	
};

#endif //SEMVRML2File_VIEWER_HH
