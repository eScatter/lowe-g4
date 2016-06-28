#ifndef SEMCountsPerCategory_h
#define SEMCountsPerCategory_h 1

#include "globals.hh"

class SEMCountsPerCategory
{
   public:
      SEMCountsPerCategory();
      G4int   nTotal,nBSE,nSE,nSE3,nEtotWindow,nEkinWindow,nAngWindow;
      G4double Current;  // Intended for output of detector current based on its characteristics (see SEMGeneralDetector)
};

#endif
