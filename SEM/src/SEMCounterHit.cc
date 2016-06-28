#include "SEMCounterHit.hh"

#include "G4Circle.hh"
#include "G4VVisManager.hh"
//#include "G4Colour.hh"

G4Allocator<SEMCounterHit> SEMCounterHitAllocator;

SEMCounterHit::SEMCounterHit()
{
  time = 0.;
}

SEMCounterHit::~SEMCounterHit()
{;}

SEMCounterHit::SEMCounterHit(const SEMCounterHit &right)
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
  depositedenergy    = right.depositedenergy;
  ParentID           = right.ParentID;
  TrackID            = right.TrackID;
  EventID            = right.EventID;
  Depth              = right.Depth;
  Rho                = right.Rho;
}

const SEMCounterHit& SEMCounterHit::operator=(const SEMCounterHit &right)
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
  depositedenergy    = right.depositedenergy;
  ParentID           = right.ParentID;
  TrackID            = right.TrackID;
  EventID            = right.EventID;
  Depth              = right.Depth;
  Rho                = right.Rho;
  return *this;
}

int SEMCounterHit::operator==(const SEMCounterHit &/*right*/) const
{
  return 0;
}

void SEMCounterHit::Draw()
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

void SEMCounterHit::Print()
{
  G4cout << "time " << time/ns
//         << " (nsec) --- local (x,y) " << localPos.x()
         << " (nsec) --- energy " << kinenergy/eV << G4endl;
//         << ", " << localPos.y() << G4endl;
}


