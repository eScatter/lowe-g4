// CADPhysicsUserTrackInfo.hh
//

#ifndef CADPhysicsUserTrackInfo_H
#define CADPhysicsUserTrackInfo_H 1

#include "globals.hh"
#include "G4VUserTrackInformation.hh"

class CADPhysicsUserTrackInfo : public G4VUserTrackInformation
{
public:
	CADPhysicsUserTrackInfo(G4double ee = 0.);
	~CADPhysicsUserTrackInfo();

	inline void Print() const {};// Dummy implementation

	inline void SetExcitationEnergy(G4double ee) { excitationenergy = ee; }// Excitation energy for excited molecules (used in CADPhysicsGasScattering)
	inline G4double GetExcitationEnergy() { return excitationenergy; }
	inline void SetMinZ(G4double value) { minZ = value; }// 'Deepest' position a particle reaches anywhere along its trajectory - used in
	// SEMDet, SEMPlane, SEMSteppingAction and SEMTrackingAction.
	inline G4double GetMinZ() { return minZ; }
	inline void SetMaxRho(G4double value) { maxRho = value; }// Largest 'rho' value (actually the radius squared) a particle reaches anywhere
	// along its trajectory; maxRho = x^2+y^2. Used in the same classes as above.
	inline G4double GetMaxRho() { return maxRho; }
	inline void SetUser1(G4double value) { uservar1 = value; }// Three extra 'G4double' variables available for arbitrary purposes
	inline G4double GetUser1() { return uservar1; }
	inline void SetUser2(G4double value) { uservar2 = value; }
	inline G4double GetUser2() { return uservar2; }
	inline void SetUser3(G4double value) { uservar3 = value; }
	inline G4double GetUser3() { return uservar3; }

private:
	G4double excitationenergy,minZ,maxRho,uservar1,uservar2,uservar3;
};

#endif
