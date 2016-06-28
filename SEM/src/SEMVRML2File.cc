// SEMVRML2File,cc, based on
// G4VRML2File.cc
// Satoshi Tanaka & Yasuhide Sawada

#include <stdio.h> // sscanf
#include <stdlib.h> // getenv

#include "G4VSceneHandler.hh"

#include "SEMVRML2File.hh"
#include "SEMVRML2FileSceneHandler.hh"
#include "SEMVRML2FileViewer.hh"

#include "G4FRClient.hh"


SEMVRML2File::SEMVRML2File() :
	G4VGraphicsSystem("VRML2FILE", "VRML2FILE", G4VGraphicsSystem::threeD)
{
}

SEMVRML2File::~SEMVRML2File()
{
}


G4VSceneHandler* SEMVRML2File::CreateSceneHandler(const G4String& name) 
{
	G4VSceneHandler *p = NULL;

	p = new SEMVRML2FileSceneHandler(*this, name);

	return p;
}

G4VViewer* SEMVRML2File::CreateViewer(G4VSceneHandler& scene, const G4String& name)
{
	G4VViewer* pView = NULL;

	SEMVRML2FileSceneHandler* pScene = (SEMVRML2FileSceneHandler*)&scene;
	pView = new SEMVRML2FileViewer(*pScene, name);

	return pView;
}
