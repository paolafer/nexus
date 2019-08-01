#ifndef CEVNS_H
#define CEVNS_H

#include "BaseGeometry.h"

namespace nexus {

  class Cevns : public BaseGeometry
  {
  public:
    ///Constructor
    Cevns();

    ///Destructor
    ~Cevns();

    void Construct();

    /// Generates a vertex within a given region of the geometry
    G4ThreeVector GenerateVertex(const G4String& region) const;

    
  private:
    
  };
}

#endif
