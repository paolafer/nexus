// ----------------------------------------------------------------------------
// nexus | Lab.h
//
// This class consists of two LXe cells/crystals placed opposite to each other.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef LAB_H
#define LAB_H

#include "BaseGeometry.h"

class G4GenericMessenger;
namespace nexus {class Pet2boxes;}
namespace nexus {class PetLYSObox;}
namespace nexus {class PetLXeCell;}
namespace nexus {class LYSOCrystal;}
namespace nexus {class LSOCeCaCrystal;}
namespace nexus {class LXeMiniCell;}
namespace nexus {class BoxPointSampler;}

namespace nexus {

  class Lab: public BaseGeometry
  {
  public:
    /// Constructor
    Lab();
    /// Destructor
    ~Lab();

    /// Return vertex within region <region> of the chamber
    virtual G4ThreeVector GenerateVertex(const G4String& region) const;
    virtual std::pair<G4ThreeVector, G4ThreeVector> GenerateVertices(const G4String& /*region*/) const;

    virtual void Construct();

  private:
    /// Messenger for the definition of control commands
    G4GenericMessenger* msg_;
    G4String mat_;
    G4double z_;

    G4double vx_, vy_, vz_;

    LYSOCrystal* lyso_module_;
    LSOCeCaCrystal* lso_module_;
    LXeMiniCell* lxe_module_;

    BoxPointSampler* generator1_;
    BoxPointSampler* generator2_;

  };

} // end namespace nexus

#endif
