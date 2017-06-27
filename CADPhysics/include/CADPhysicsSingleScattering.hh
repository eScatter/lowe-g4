// CADPhysicsSingleScattering.hh
//
// The process that takes care of elastic scattering of the electron in any material.
// This process combines various descriptions for different energy ranges and types of material - going from high to low energy:
// - Above 30 keV, Browning's approximation is used for the cross sections and angular dependencies;
//   Ref.: R. Browning et al, J.Appl.Phys. 76 (4), 2016 (1994).
// - Between 30 keV and 100 eV (metals) or 200 eV (semiconductors and insulators), the Mott cross sections are applied;
// - For metals, between 100 eV and the Fermi energy the (differential) inverse mean free paths are a linear interpolation between the Mott values
//   and the acoustic phonon scattering inverse mean free path at the Fermi energy. The latter is derived from the material's resistivity.
// - For semiconductors and insulators, below 100 eV (differential) inverse mean free paths for phonon scattering are used. The model is a slight simplification
//   of the one described in H.-J.Fitting et al., J.Electron Spectrosc.Rel.Phenom. 119 (2001) 35-47 (scatter rates)
//   and E.Schreiber & H.-J.Fitting, J.Electron Spectrosc.Rel.Phenom. 124 (2002) 25-37 (angular dependencies).
//   Between 100 eV and 200 eV, an interpolation between the phonon and Mott data is applied.
//

#ifndef CADPhysicsSingleScattering_h
#define CADPhysicsSingleScattering_h 1
#include "G4VContinuousDiscreteProcess.hh"

#include "CADPhysicsUnits.hh"

#include "G4Electron.hh"
#include "G4DataVector.hh"

// for read in hdf5
#include <csread/material.h>
#include <csread/units/unit_system.h>

typedef std::vector<G4DataVector*> CADPhysicsDataTable;
typedef std::vector<CADPhysicsDataTable*> CADPhysicsDataCube;

class CADPhysicsSSMessenger;
class G4MaterialCutsCouple;
class G4Material;

class CADPhysicsSingleScattering : public G4VContinuousDiscreteProcess

{
public:
	// The constructor method is hidden as a 'protected' method. There is a public GetInstance() method instead.
	// The reason is that the CADPhysicsTransportation and CADPhysicsSingleScattering methods need to be able to work together
	// 'under water', and for this the Transportation method always needs to be able to find the right instance of SingleScattering.
	// The simplest way to ensure this is to make sure that only a single instance (the 'singleton') of the SingleScattering class
	// can ever be created.
	static CADPhysicsSingleScattering* GetInstance();

	~CADPhysicsSingleScattering();

	G4bool IsApplicable ( const G4ParticleDefinition& );
	// Returns true for the electron, false otherwise

  void load_material(std::string const & filename, double high_energy, size_t N_K, size_t N_P);
  void load_vacuum();

	void BuildPhysicsTable(const G4ParticleDefinition& aParticleType);
	// This function overloads the corresponding virtual function
	// of the base class G4VContinuousDiscreteProcess.
	// It is an initialization routine invoked from G4VUserPhysicsList::BuildPhysicsTable()
	// at the start of the first run (and after geometry modification)
	// It takes care of importing/calculating differential inverse mean free paths and (transport) mean free paths.

	void PrintInfoDefinition();
	// Print few lines of informations about the process: validity range,
	// origin of data, etc.
	// Invoked by BuildPhysicsTable().

	G4double PostStepGetPhysicalInteractionLength(
		const G4Track& track,
		G4double   previousStepSize,
		G4ForceCondition* condition
		);
	// Overloading the base class, with the main purpose of implementing the 'safetycheck' feature in DoMultistep.
	// In case the particle crosses the 'safety' boundary in DoMultistep, the already calculated length of the next step
	// is not executed but rather stored for the subsequent step.
	//
	// Note: also modifies the verbosity level criteria of the base class. The various CADPhysics processes do not have
	// consistent verbosity levels (yet). The following might serve as a guideline for future modifications:
	// 0 - Silent except for essential warning messages
	// 1 - Basic information during initialization
	// 2 - More extensive information during initialization
	// 3 - For output of e.g. cross section tables to screen during initialization
	// 4 and higher - For (debugging) output during simulations, e.g. in
	//                PostStepGetPhysicalInteractionLength, or perhaps in PostStepDoIt

	G4double AlongStepGetPhysicalInteractionLength(const G4Track&,
		G4double  previousStepSize,
		G4double  currentMinimumStep,
		G4double& currentSafety,
		G4GPILSelection* selection);
	// Special implementation of this method. Determines what kind of step to do (single step, multistep or
	// nothing at all). Communicates with the Transportation process to learn what the geometry-limited stepsize is and to get an
	// updated safety value. To this end, a special RealASGPIL method has been created in CADPhysicsTransportation.
	// In case of a multistep, it also invokes the corresponding DoMultiStep method so that it can
	// return the value of the actual step length taken (which can then be determined 'on the fly' by the abovenmentioned methods).

	inline G4double GetContinuousStepLimit(const G4Track&,G4double,G4double,G4double&){return DBL_MAX;}// Dummy implementation only

	G4double GetMeanFreePath(const G4Track& aTrack,
		G4double previousStepSize,
		G4ForceCondition* condition=0);
	// Returns the current mean free path. Called from PostStepGetPhysicalInteractionLength.

	G4VParticleChange* AlongStepDoIt(const G4Track& aTrack,const G4Step& aStep);

	G4VParticleChange* PostStepDoIt(const G4Track& aTrack,const G4Step& aStep);

	void SetDoMultistep(G4bool dm);// Set the value of the 'domultistep' boolean, accessible by the user via the Messenger.
	G4bool GetDoMultistep();// Returns the value of domultistep

	G4double GetFirstsinglestep(G4bool& multistepping);// For the benefit of CADPhysicsTransportation
	// Returns lastgpilvalue, while multistepping is set to domultistep.

	void NoRaytrace();// Sets the noraytrace boolean to true. Can be called from CADPhysicsTransportation.

protected:
	// 'Hidden' constructor method
	CADPhysicsSingleScattering(const G4String& processName="ssc");

private:
	//  Hidden assignment operator
	CADPhysicsSingleScattering & operator = (const CADPhysicsSingleScattering &right);
	CADPhysicsSingleScattering ( const CADPhysicsSingleScattering &);

	G4double DoMultiStep(const G4Track& track,
		G4double othersteplength,
		G4double safety);
	// Multistep version of the PostStepDoIt method. In fact, as it needs to propose a *position*,
	// and to set a flexible total step length, it is invoked from AlongStepGetPhysicalInteractionLength.
	// It returns the total step taken, and stores the proposed final position, direction, energy loss and final energy.
	// These parameters are passed on to the ParticleChange in the AlongStep/PostStepDoIt methods.

	G4VParticleChange* DoSingleStep(const G4Track& track, const G4Step& step);
	// Process a single scattering event in the conventional way. Called by PostStepDoIt.

private: // Data members
	static CADPhysicsSingleScattering* fSingleScattering;// The instance

	CADPhysicsSSMessenger* messenger;// The messenger

  // Some material parameters
  std::vector<G4double> vec_effective_A;
     // effective atomic mass in gram per material

  std::vector<imfp_table<float>> elastic_imfp_vector;
  std::vector<icdf_table<float>> elastic_icdf_vector;

	// Pointers for the current material as determined by DefineMaterial each time we're in a different material
	const G4MaterialCutsCouple* currentCouple;
	const G4Material*           currentMaterial;
	size_t                      currentMaterialIndex;

	// Some constants for energy loss
	std::vector<G4double> vec_phononloss;
	G4double fphononloss;// Phonon loss for current material
	const G4double enlossconst;// Energy loss constant for non-phonon energy loss

	// Some variables shared among various methods during stepping
	G4double previousenergy;
	G4double previousMFP;

	// Booleans that can be set by the user, to switch multistep method on or off
	G4bool	domultistep;

	// (Mean) free path
	G4double lastgpilvalue;// Actual next step to be taken as determined by PostStepGetPhysicalInteractionLength()

	// Booleans that determine the kind of step
	G4bool donothing;
	G4bool singlestep;
	G4bool multistep;
	G4bool noraytrace;// Set by CADPhysicsTransportation
	G4bool safetycheck;// Set by DoMultiStep if safety was limiting

	// Results stored by DoMultiStep that are passed on to the ParticleChange by AlongStepDoIt and PostStepDoIt
	G4ThreeVector elDirectionnew;
	G4double finalT;
	G4double energyloss;
	G4ThreeVector proposeposition;
};

inline G4bool CADPhysicsSingleScattering::IsApplicable(const G4ParticleDefinition& particle)
{
	return ( &particle == G4Electron::Electron() );
}

inline void CADPhysicsSingleScattering::SetDoMultistep (G4bool dm)// Set 'domultistep' boolean.
{
	domultistep = dm;
	G4cout << "domultistep     ";
	if (domultistep) G4cout << "true" << G4endl; else G4cout << "false" << G4endl;
}

inline G4bool CADPhysicsSingleScattering::GetDoMultistep()// Corresponding 'get' method.
{
	return domultistep;
}

inline G4double CADPhysicsSingleScattering::GetFirstsinglestep(G4bool& multistepping)
{
	multistepping = domultistep;
	return lastgpilvalue;
}

inline void CADPhysicsSingleScattering::NoRaytrace()
{
	noraytrace = true;
}

#endif
