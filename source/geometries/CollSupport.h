#ifndef COLL_SUPPORT_H
#define COLL_SUPPORT_H

#include "BaseGeometry.h"

namespace nexus {
  class CollSupport : public BaseGeometry {

  public:
    ///Constructor
    CollSupport();

    ///Destructor
    ~CollSupport();

    void Construct();
    G4double GetAxisCentre();
    G4double GetYDisplacement();

  private:
    G4double _axis_centre;
    G4double _y_displ;
  };
}

#endif