// ----------------------------------------------------------------------------
// nexus | LYSOLab.h
//
// This class consists of a wrapper for an arbitrary number of LYSO crystals
// arranged in rings.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef LYSO_LAB_H
#define LYSO_LAB_H

#include "BaseGeometry.h"

class G4GenericMessenger;
namespace nexus {class LYSOCrystal;}

namespace nexus {

  class LYSOLab: public BaseGeometry
  {
  public:
    /// Constructor
    LYSOLab();
    /// Destructor
    ~LYSOLab();

    /// Return vertex within region <region> of the chamber
    virtual G4ThreeVector GenerateVertex(const G4String& region) const;

    virtual void Construct();

  private:
    /// Messenger for the definition of control commands
    G4GenericMessenger* msg_;
    LYSOCrystal* lyso_module_;

    G4int n_crystal_rows_;
    G4double crystal_pitch_;
    G4double inner_radius_, depth_, sec_size_;

  };

} // end namespace nexus

#endif
