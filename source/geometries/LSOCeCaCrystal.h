// ----------------------------------------------------------------------------
// nexus | PetLSOCeCaCrystal.h
//
// Basic cell made of LSOCeCa.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef PET_LSOCeCa_CRYSTAL_H
#define PET_LSOCeCa_CRYSTAL_H

#include "BaseGeometry.h"

//class G4LogicalVolume;
class G4GenericMessenger;

namespace nexus {

  class LSOCeCaCrystal: public BaseGeometry {
  public:
    /// Constructor
    LSOCeCaCrystal();

    /// Destructor
    ~LSOCeCaCrystal();

    /// Generate a vertex within a given region of the geometry
    G4ThreeVector GenerateVertex(const G4String& region) const;

    void Construct();

    void SetSecSize(G4double size);
    void SetDOISize(G4double size);

    G4double GetSecSize() const;
    G4double GetDOISize() const;

    void SetReflectivity(G4double refl);

  private:

    G4double active_size_;
    G4double lso_zsize_;
    G4double max_step_size_;

    G4double reflectivity_;

  };

  inline void LSOCeCaCrystal::SetSecSize(G4double size) {active_size_ = size;}
  inline void LSOCeCaCrystal::SetDOISize(G4double size) {lso_zsize_  = size;}

  inline G4double LSOCeCaCrystal::GetSecSize() const {return active_size_;}
  inline G4double LSOCeCaCrystal::GetDOISize() const {return lso_zsize_;}

  inline void LSOCeCaCrystal::SetReflectivity(G4double refl) {reflectivity_ = refl;} 

} // end namespace nexus

#endif
