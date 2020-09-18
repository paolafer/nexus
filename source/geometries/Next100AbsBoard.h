// -----------------------------------------------------------------------------
// nexus | Next100AbsBoard.h
//
// Geometry of the NEXT-100 SiPM board, consisting of an 8x8 array of
// silicon photomultipliers (1.3x1.3 mm2 of active area) mounted on a Kapton
// board covered with a completely absorbant surface.
//
// The NEXT Collaboration
// -----------------------------------------------------------------------------

#ifndef NEXT100_ABS_BOARD_H
#define NEXT100_ABS_BOARD_H

#include "BaseGeometry.h"
#include <G4ThreeVector.hh>
#include <vector>

class G4VPhysicalVolume;
class G4GenericMessenger;

namespace nexus {

  class BoxPointSampler;
  class GenericPhotosensor;

  // Geometry of the 8x8 SiPM boards used in the tracking plane of NEXT-100

  class Next100AbsBoard: public BaseGeometry
  {
  public:
    // Default constructor
    Next100AbsBoard();
    // Destructor
    ~Next100AbsBoard();
    //
    void SetMotherPhysicalVolume(G4VPhysicalVolume*);
    //
    void Construct() override;
    //
    G4ThreeVector GenerateVertex(const G4String&) const override;

    G4double GetSize() const;
    G4double GetThickness() const;

    const std::vector<G4ThreeVector>& GetSiPMPositions() const;

  private:
    G4GenericMessenger* msg_;
    G4double size_, pitch_, margin_;
    G4double hole_diam_;
    G4double board_thickness_, mask_thickness_;
    G4double time_binning_;
    std::vector<G4ThreeVector> sipm_positions_;
    G4bool   visibility_, sipm_visibility_;
    G4VPhysicalVolume*  mpv_;
    BoxPointSampler*    vtxgen_;
    GenericPhotosensor* sipm_;
  };

  inline void Next100AbsBoard::SetMotherPhysicalVolume(G4VPhysicalVolume* p)
  { mpv_ = p;}

  inline G4double Next100AbsBoard::GetSize() const
  { return size_; }

  inline G4double Next100AbsBoard::GetThickness() const
  { return (board_thickness_ + mask_thickness_); }

  inline const std::vector<G4ThreeVector>& Next100AbsBoard::GetSiPMPositions() const
  { return sipm_positions_; }

} // namespace nexus

#endif
