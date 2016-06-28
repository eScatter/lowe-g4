// SEMVisExecutive.hh

#ifndef SEMVISEXECUTIVE_HH
#define SEMVISEXECUTIVE_HH

#include "G4VisExecutive.hh"

class SEMVisExecutive: public G4VisExecutive {

public: // With description

  SEMVisExecutive ();

private:

  void RegisterGraphicsSystems ();

};

#include "SEMVisExecutive.icc"

#endif

