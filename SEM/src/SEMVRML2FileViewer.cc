// SEMVRML2FileViewer.cc, based on
// G4VRML2FileViewer.cc
// Satoshi Tanaka & Yasuhide Sawada


//#define DEBUG_FR_VIEW

#include <cmath>

#include "G4Scene.hh"
#include "SEMVRML2FileViewer.hh"
#include "SEMVRML2FileSceneHandler.hh"
#include "SEMVRML2File.hh"
#include "G4ios.hh"

SEMVRML2FileViewer::SEMVRML2FileViewer(SEMVRML2FileSceneHandler& sceneHandler,
				 const G4String& name) :
 G4VViewer(sceneHandler,
	   sceneHandler.IncrementViewCount(),
	   name),
 fSceneHandler(sceneHandler),
 fDest(sceneHandler.fDest)
{
	fViewHalfAngle = 0.5 * 0.785398 ; // 0.5 * 45*deg
	fsin_VHA       = std::sin ( fViewHalfAngle ) ;	
}

SEMVRML2FileViewer::~SEMVRML2FileViewer()
{}

void SEMVRML2FileViewer::SetView()
{
#if defined DEBUG_FR_VIEW
  G4cerr << "***** SEMVRML2FileViewer::SetView(): No effects" << G4endl;
#endif

// Do nothing, since VRML a browser is running as a different process.
// SendViewParameters () will do this job instead.

}

void SEMVRML2FileViewer::DrawView()
{
#if defined DEBUG_FR_VIEW
	G4cerr << "***** SEMVRML2FileViewer::DrawView()" << G4endl;
#endif

	fSceneHandler.VRMLBeginModeling() ; 

    G4double writeunit = fSceneHandler.Setlunit();

	fDest << "#Length unit for visualization: " << writeunit/mm << " mm." << G4endl;

    // Viewpoint node
    SendViewParameters(); 

	// Here is a minimal DrawView() function.
	NeedKernelVisit();
	ProcessView();
	FinishView();
}

void SEMVRML2FileViewer::ClearView(void)
{
#if defined DEBUG_FR_VIEW
  G4cerr << "***** SEMVRML2File1View::ClearView(): No effects" << G4endl;
#endif
}

void SEMVRML2FileViewer::ShowView(void)
{
#if defined DEBUG_FR_VIEW
  G4cerr << "***** SEMVRML2FileViewer::ShowView()" << G4endl;
#endif
	fSceneHandler.VRMLEndModeling();
}

void SEMVRML2FileViewer::FinishView(void)
{
#if defined DEBUG_FR_VIEW
  G4cerr << "***** SEMVRML2FileViewer::FinishView(): No effects" << G4endl;
#endif
}

void SEMVRML2FileViewer::SendViewParameters () 
{
  // Calculates view representation based on extent of object being
  // viewed and (initial) direction of camera.  (Note: it can change
  // later due to user interaction via visualization system's GUI.)

#if defined DEBUG_FR_VIEW
      G4cerr << "***** SEMVRML2FileViewer::SendViewParameters()\n";
#endif 
    
	fViewHalfAngle = fVP.GetFieldHalfAngle();
	fsin_VHA       = std::sin ( fViewHalfAngle ) ;	

	// error recovery
	if ( fsin_VHA < 1.0e-6 ) { //return ; } 
		fViewHalfAngle = 0.5 * 0.785398 ; // 0.5 * 45*deg
		fsin_VHA       = std::sin ( fViewHalfAngle ) ;
    }

	// camera distance
	G4double extent_radius = fSceneHandler.GetScene()->GetExtent().GetExtentRadius();
	G4double camera_distance = extent_radius / fsin_VHA ;

	// camera position on Z axis
	const G4Point3D&	target_point
	  = fSceneHandler.GetScene()->GetStandardTargetPoint()
	  + fVP.GetCurrentTargetPoint();
	G4Vector3D		CD = fVP.GetViewpointDirection();
	G4Vector3D		upvector = fVP.GetUpVector();
	Rotation = 0;
	G4Vector3D		Rotvector = FindVectorAndRotation(CD,upvector);//also updates Rotation
	G4Vector3D		camera_vector = camera_distance * CD;
	G4double		E_x = target_point.x() + camera_vector.x();
	G4double		E_y = target_point.y() + camera_vector.y();
	G4double		E_z = target_point.z() + camera_vector.z();
	G4Point3D		E(E_x, E_y, E_z );
	//G4double		R_x = 1;
	//G4double		R_y = 0;
	//if (1!=CD.z()) {
	//  G4double CD_xy = sqrt(1-CD.z()*CD.z());
	//  R_x = -1*CD.y()/CD_xy;
	//  R_y = CD.x()/CD_xy;
	//  Rotation = acos(CD.z());
	//}
	G4double lunit = fSceneHandler.Getlunit();
	G4double fieldOfView = 1.5*fViewHalfAngle/fVP.GetZoomFactor();
	
	// VRML codes are generated below	
	fDest << G4endl;
        // EB: Added some extra lights
	fDest << "#---------- LIGHTS" << G4endl;
	fDest << "DirectionalLight {"         << G4endl;
	fDest << "\t" << "intensity 0.5 "           ;
	fDest << "\t" << "direction 0.0 1.0 0.0 "           ;
	fDest << "}" << G4endl;
	fDest << "DirectionalLight {"         << G4endl;
	fDest << "\t" << "intensity 0.5 "           ;
	fDest << "\t" << "direction 0.0 0.0 -1.0 "           ;
	fDest << "}" << G4endl;
	fDest << "DirectionalLight {"         << G4endl;
	fDest << "\t" << "intensity 0.5 "           ;
	fDest << "\t" << "direction 0.0 0.0  1.0 "           ;
	fDest << "}" << G4endl;
	fDest << G4endl;
        // EB: end of addition
	fDest << "#---------- CAMERA" << G4endl;
	fDest << "Viewpoint {"         << G4endl;
	fDest << "\t" << "position "           ;
	fDest                 << E.x()/lunit << " "  ;
	fDest                 << E.y()/lunit << " "  ;
	fDest                 << E.z()/lunit << G4endl ;
	fDest << "\t" << "orientation "		   ;
	fDest				  << Rotvector.x() << " "  ;
	fDest				  << Rotvector.y() << " "  ;
	fDest				  << Rotvector.z() << " "  ;
	fDest				  << Rotation << G4endl;
	fDest << "\t" << "fieldOfView ";
	fDest				  << fieldOfView << G4endl;
	fDest << "}" << G4endl;
        // EB: Added extra viewpoints
	fDest << "Viewpoint {"         << G4endl;
	fDest << "\t" << "position 0.0 -100.0 1.3 " << G4endl;
	fDest << "\t" << "fieldOfView 0.05 " << G4endl;
        fDest << "\t" << "description \"Elstar\" "  << G4endl;
	fDest << "}" << G4endl;
	fDest << "Viewpoint {"         << G4endl;
	fDest << "\t" << "position 0.0 -90.0 15.0 " << G4endl;
	fDest << "\t" << "orientation 1.0 0.0 0.0 1.570796 " << G4endl;
	fDest << "\t" << "fieldOfView 0.05 " << G4endl;
        fDest << "\t" << "description \"AngDist\" "  << G4endl;
	fDest << "}" << G4endl;
        // EB: end of addition
	fDest << G4endl;
}

G4double SEMVRML2FileViewer::InProduct(G4Vector3D vecu, G4Vector3D vecv)
{
	G4double result = vecu.x()*vecv.x()+vecu.y()*vecv.y()+vecu.z()*vecv.z();
	return result;
}

G4Vector3D SEMVRML2FileViewer::OutProduct(G4Vector3D vecu, G4Vector3D vecv)
{
	G4double resultx = vecu.y()*vecv.z()-vecu.z()*vecv.y();
    G4double resulty = vecu.z()*vecv.x()-vecu.x()*vecv.z();
	G4double resultz = vecu.x()*vecv.y()-vecu.y()*vecv.x();
	return G4Vector3D(resultx,resulty,resultz);
}

G4Vector3D SEMVRML2FileViewer::FindVectorAndRotation(G4Vector3D Cdin, G4Vector3D upvectorin) 
{
	G4double cdinup = InProduct(Cdin,upvectorin);//inner product
	G4Vector3D normup = upvectorin - cdinup*Cdin;//make normup orthogonal to Cdin
	normup = normup.unit();//normalize to unit length
	G4Vector3D t = normup - G4Vector3D(0,1,0);
	G4Vector3D u = Cdin - G4Vector3D(0,0,1);
	G4Vector3D r = OutProduct(t,u);
	if(r!=G4Vector3D(0,0,0)) {
		r=r.unit();
	} else {
		G4double r_x=0;
		G4double r_y=-t.z();
		G4double r_z=t.y();
		r = G4Vector3D(r_x,r_y,r_z);
		if(r!=G4Vector3D(0,0,0)) { r=r.unit(); } else {
			Rotation = 0;
			return G4Vector3D(1,0,0);
		}
	}
	G4Vector3D p = G4Vector3D(0,1,0);
	G4Vector3D q = normup;
	G4double k = InProduct(p,r);
	if(1==k*k) {
	  p = G4Vector3D(0,0,1);
	  q = Cdin;
	  k = InProduct(p,r);
	}
    G4double l = InProduct(OutProduct(p,q),r)/(1-k*k);
	Rotation = asin(l);
	if(InProduct(p,q)-k*k<0) Rotation = pi - Rotation;
	return r;
}
