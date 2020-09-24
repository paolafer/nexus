// ----------------------------------------------------------------------------
// nexus | PetLYSOCrystal.h
//
// Basic cell made of LYSO.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef PET_LYSO_CRYSTAL_H
#define PET_LYSO_CRYSTAL_H

#include "BaseGeometry.h"

class G4LogicalVolume;
class G4GenericMessenger;

/* namespace nexus {class SiPMpetVUV;} */
/* namespace nexus {class SiPMpetTPB;} */

namespace nexus {

  class LYSOCrystal: public BaseGeometry {
  public:
    /// Constructor
    LYSOCrystal();

    /// Destructor
    ~LYSOCrystal();

    /// Generate a vertex within a given region of the geometry
    G4ThreeVector GenerateVertex(const G4String& region) const;

    void Construct();

  private:

    // Detector dimensions
    const G4double active_size_; /// Size of the LYSO active volume
    
    // Parameters
    G4double max_step_size_;  /// Maximum Step Size

    /// Messenger for the definition of control commands
    G4GenericMessenger* msg_;

    G4double lyso_zsize_;

  };

} // end namespace nexus

#endif
