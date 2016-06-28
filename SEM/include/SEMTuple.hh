// SEMTuple.hh

#ifndef SEMTuple_h
#define SEMTuple_h 1

#include "globals.hh"
#include <vector>
#include <algorithm>
using namespace std;

template <class SEMTupleElement> class SEMTuple
{
protected:
   vector<SEMTupleElement> m_items;

public:
   inline void Add(SEMTupleElement element) {
     m_items.push_back(element); 
   }

   //Return the memory address of a specific item
   inline SEMTupleElement* GetAddress(int ItemKey) {
     return &(m_items[ItemKey]);
   }

   inline void SortEnergy(void) { for (unsigned int i=0; i<m_items.size(); i++) m_items[i].sort=m_items[i].kinen; sort(m_items.begin(),m_items.end());}
   inline void SortRho(void) { for (unsigned int i=0; i<m_items.size(); i++) m_items[i].sort=m_items[i].localposition.rho(); sort(m_items.begin(),m_items.end());}

   //Clear the collection 
   inline void Clear(void) {
     m_items.clear();
   }
   //Return the number of items in collection
   inline int Count(void) {
     return m_items.size(); //One Based
   }
   //Operator Returning a reference to TBase
   inline SEMTupleElement& operator [](int ItemKey) {
     return m_items[ItemKey];
   }
};

#endif
