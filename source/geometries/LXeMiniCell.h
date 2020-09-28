// ----------------------------------------------------------------------------
// nexus | LXeMiniCell.h
//
// Basic cell made of LXe.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef LXE_MINI_CELL_H
#define LXE_MINI_CELL_H

#include "BaseGeometry.h"

class G4LogicalVolume;
class G4GenericMessenger;

/* namespace nexus {class SiPMpetVUV;} */
/* namespace nexus {class SiPMpetTPB;} */

namespace nexus {

  class LXeMiniCell: public BaseGeometry {
  public:
    /// Constructor
    LXeMiniCell();

    /// Destructor
    ~LXeMiniCell();

    /// Generate a vertex within a given region of the geometry
    G4ThreeVector GenerateVertex(const G4String& region) const;

    void Construct();

  private:

    // Detector dimensions
    const G4double active_size_; /// Size of the LYSO active volume
    
    // Parameters
    G4double max_step_size_;  /// Maximum Step Size

    /// Messenger for the definition of control commands
    // G4GenericMessenger* msg_;

    G4double lxe_zsize_;

  };

} // end namespace nexus

#endif
