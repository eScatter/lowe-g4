// SEMVRML2File.hh, based on G4VRML2File

#ifndef SEMVRML2FILE_HH
#define SEMVRML2FILE_HH

#include "G4VGraphicsSystem.hh"

class G4VSceneHandler;

class SEMVRML2File: public G4VGraphicsSystem {
public:
	SEMVRML2File(); 
	virtual ~SEMVRML2File();
	G4VSceneHandler* CreateSceneHandler(const G4String& name = "");
	G4VViewer*  CreateViewer(G4VSceneHandler&, const G4String& name = "");

};

#endif //SEMVRML2FILE_HH
