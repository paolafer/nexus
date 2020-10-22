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

//class G4LogicalVolume;
class G4GenericMessenger;

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

    void SetSecSize(G4double size);
    void SetDOISize(G4double size);

    G4double GetSecSize() const;
    G4double GetDOISize() const;

  private:

    G4double active_size_;
    G4double lyso_zsize_;
    G4double max_step_size_;

  };

  inline void LYSOCrystal::SetSecSize(G4double size) {active_size_ = size;}
  inline void LYSOCrystal::SetDOISize(G4double size) {lyso_zsize_  = size;}

  inline G4double LYSOCrystal::GetSecSize() const {return active_size_;}
  inline G4double LYSOCrystal::GetDOISize() const {return lyso_zsize_;}

} // end namespace nexus

#endif
