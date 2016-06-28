#include "SEMBulkDetHit.hh"

#include "G4Circle.hh"
#include "G4VVisManager.hh"
//#include "G4Colour.hh"

G4Allocator<SEMBulkDetHit> SEMBulkDetHitAllocator;

SEMBulkDetHit::SEMBulkDetHit()
{
  time = 0.;
}

SEMBulkDetHit::~SEMBulkDetHit()
{;}

SEMBulkDetHit::SEMBulkDetHit(const SEMBulkDetHit &right)
    : G4VHit() {
  time               = right.time;
  localPos           = right.localPos;
  worldPos           = right.worldPos;
  momentumdirection  = right.momentumdirection;
  vertexPos          = right.vertexPos;
  kinenergy          = right.kinenergy;
  lv                 = right.lv;
  process            = right.process;
  category           = right.category;
  ParentID           = right.ParentID;
  TrackID            = right.TrackID;
  EventID            = right.EventID;
  Depth              = right.Depth;
  Rho                = right.Rho;
}

const SEMBulkDetHit& SEMBulkDetHit::operator=(const SEMBulkDetHit &right)
{
  time               = right.time;
  localPos           = right.localPos;
  worldPos           = right.worldPos;
  momentumdirection  = right.momentumdirection;
  vertexPos          = right.vertexPos;
  kinenergy          = right.kinenergy;
  lv                 = right.lv;
  process            = right.process;
  category           = right.category;
  ParentID           = right.ParentID;
  TrackID            = right.TrackID;
  EventID            = right.EventID;
  Depth              = right.Depth;
  Rho                = right.Rho;
  return *this;
}

int SEMBulkDetHit::operator==(const SEMBulkDetHit &/*right*/) const
{
  return 0;
}

void SEMBulkDetHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(worldPos);
    circle.SetScreenSize(2);
    circle.SetFillStyle(G4Circle::filled);
//    G4Colour colour(1.,1.,0.);
//    G4VisAttributes attribs(colour);
//    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

void SEMBulkDetHit::Print()
{
 // Detector gives no direct output to screen. Use output routines instead.
}


