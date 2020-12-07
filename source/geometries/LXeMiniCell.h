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

    void SetSecSize(G4double size);
    void SetDOISize(G4double size);

    G4double GetSecSize() const;
    G4double GetDOISize() const;

  private:

    // Detector dimensions
    G4double active_size_; /// Size of the LYSO active volume
    
    // Parameters
    G4double max_step_size_;  /// Maximum Step Size

    /// Messenger for the definition of control commands
    // G4GenericMessenger* msg_;

    G4double lxe_zsize_;

  };

  inline void LXeMiniCell::SetSecSize(G4double size) {active_size_ = size;}
  inline void LXeMiniCell::SetDOISize(G4double size) {lxe_zsize_   = size;}

  inline G4double LXeMiniCell::GetSecSize() const {return active_size_;}
  inline G4double LXeMiniCell::GetDOISize() const {return lxe_zsize_;}

} // end namespace nexus

#endif
