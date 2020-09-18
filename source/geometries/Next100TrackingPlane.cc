// ----------------------------------------------------------------------------
// nexus | Next100TrackingPlane.cc
//
// Tracking plane of the NEXT-100 geometry.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "Next100TrackingPlane.h"

#include "Next100AbsBoard.h"
#include "CylinderPointSampler2020.h"
#include "Visibilities.h"

#include <G4GenericMessenger.hh>
#include <G4Tubs.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4NistManager.hh>
#include <Randomize.hh>
#include <G4VisAttributes.hh>

using namespace nexus;


Next100TrackingPlane::Next100TrackingPlane(G4double origin_z_coord):
  BaseGeometry(),
  z0_(origin_z_coord), // Distance between gate and inner face of copper plate
  copper_plate_diameter_  (1340.*mm),
  copper_plate_thickness_ ( 120.*mm),
  distance_board_board_   (   1.*mm),
  visibility_(true),
  sipm_board_geom_(new Next100AbsBoard),
  copper_plate_gen_(nullptr),
  mpv_(nullptr),
  msg_(nullptr)
{
  msg_ = new G4GenericMessenger(this, "/Geometry/Next100/",
                                "Control commands of the NEXT-100 geometry.");

  msg_->DeclareProperty("tracking_plane_vis", visibility_,
                        "Visibility of the tracking plane volumes.");
}


Next100TrackingPlane::~Next100TrackingPlane()
{
  delete msg_;
  delete sipm_board_geom_;
  delete copper_plate_gen_;
}


void Next100TrackingPlane::Construct()
{
  // Make sure the pointer to the mother volume is actually defined
  if (!mpv_)
    G4Exception("[Next100TrackingPlane]", "Construct()",
                FatalException, "Mother volume is a nullptr.");

  // COPPER PLATE ////////////////////////////////////////////////////

  G4String copper_plate_name = "TP_COPPER_PLATE";

  G4Tubs* copper_plate_solid_vol = new G4Tubs(copper_plate_name,
                                             0., copper_plate_diameter_/2.,
                                             copper_plate_thickness_/2.,
                                             0., 360.*deg);

  G4LogicalVolume* copper_plate_logic_vol =
    new G4LogicalVolume(copper_plate_solid_vol,
                        G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu"),
                        copper_plate_name);

  G4double zpos = GetELzCoord() - z0_ - copper_plate_thickness_/2.;

  G4VPhysicalVolume* copper_plate_phys_vol =
    new G4PVPlacement(nullptr, G4ThreeVector(0.,0.,zpos),
                      copper_plate_logic_vol, copper_plate_name, mpv_->GetLogicalVolume(),
                      false, 0, false);

  copper_plate_gen_ = new CylinderPointSampler2020(copper_plate_phys_vol);

  // SIPM BOARDS /////////////////////////////////////////////////////

  sipm_board_geom_->SetMotherPhysicalVolume(mpv_);
  sipm_board_geom_->Construct();
  G4LogicalVolume* sipm_board_logic_vol = sipm_board_geom_->GetLogicalVolume();

  zpos = GetELzCoord() - z0_ + sipm_board_geom_->GetThickness()/2.;

  // SiPM boards are positioned bottom (negative Y) to top (positive Y)
  // and left (negative X) to right (positive X).

  G4int board_index = 1;

  // Column on the far left has 5 boards.
  // It is located 3.5 boards away from the center.
  PlaceSiPMBoardColumns(5, -3.5, zpos, board_index, sipm_board_logic_vol);

  // Second column from the left has 7 boards.
  // It is located 2.5 boards away from the center.
  PlaceSiPMBoardColumns(7, -2.5, zpos, board_index, sipm_board_logic_vol);

  // Central block of 4 columns with 8 boards each
  PlaceSiPMBoardColumns(8, -1.5, zpos, board_index, sipm_board_logic_vol);
  PlaceSiPMBoardColumns(8, -0.5, zpos, board_index, sipm_board_logic_vol);
  PlaceSiPMBoardColumns(8,  0.5, zpos, board_index, sipm_board_logic_vol);
  PlaceSiPMBoardColumns(8,  1.5, zpos, board_index, sipm_board_logic_vol);

  // Second column from the right has 7 boards and
  // it's located 2.5 boards away from the center.
  PlaceSiPMBoardColumns(7,  2.5, zpos, board_index, sipm_board_logic_vol);

  // Column on the far right has 5 boards.
  // It is located 3.5 boards away from the center.
  PlaceSiPMBoardColumns(5,  3.5, zpos, board_index, sipm_board_logic_vol);

  // VISIBILITIES //////////////////////////////////////////

  if (visibility_) {
    G4VisAttributes copper_brown = CopperBrown();
    copper_plate_logic_vol->SetVisAttributes(copper_brown);
  } else {
    copper_plate_logic_vol->SetVisAttributes(G4VisAttributes::Invisible);
    sipm_board_logic_vol  ->SetVisAttributes(G4VisAttributes::Invisible);
  }

}


void Next100TrackingPlane::PlaceSiPMBoardColumns(G4int num_boards,
                                                 G4double distance_from_center,
                                                 G4double zpos,
                                                 G4int& board_index,
                                                 G4LogicalVolume* logic_vol)
{
  G4double size = sipm_board_geom_->GetSize() + distance_board_board_;

  G4double xpos = distance_from_center * size;

  for (auto i=0; i<num_boards; i++) {
    G4double ypos = (- 0.5 * (num_boards - 1) + i ) * size;
    G4ThreeVector position(xpos, ypos, zpos); board_pos_.push_back(position);
    new G4PVPlacement(nullptr, position,
                      logic_vol, logic_vol->GetName(), mpv_->GetLogicalVolume(),
                      false, board_index, false);
    board_index++;
  }
}


void Next100TrackingPlane::PrintSiPMPositions() const
{
  auto sipm_positions = sipm_board_geom_->GetSiPMPositions();

  for (unsigned int i=0; i<board_pos_.size(); ++i) {
    for (unsigned int j=0; j<sipm_positions.size(); ++j) {
      G4int id = 1000 * (i+1) + j;
      G4ThreeVector pos = sipm_positions[j] + board_pos_[i];
      G4cout << "SiPM " << id << "  " << pos << G4endl;
    }
  }
}


G4ThreeVector Next100TrackingPlane::GenerateVertex(const G4String& region) const
{
  G4ThreeVector vertex;

  if (region == "SIPM_BOARD") {
    vertex = sipm_board_geom_->GenerateVertex("");
    G4int board_num = G4RandFlat::shootInt((long) 0, board_pos_.size());
    vertex += board_pos_[board_num];
  }
  else if (region == "TP_COPPER_PLATE") {
    vertex = copper_plate_gen_->GenerateVertex("VOLUME");
  }

  return vertex;
}
