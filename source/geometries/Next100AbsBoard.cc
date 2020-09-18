// -----------------------------------------------------------------------------
// nexus | Next100AbsBoard.cc
//
// Geometry of the NEXT-100 SiPM board, consisting of an 8x8 array of
// silicon photomultipliers (1.3x1.3 mm2 of active area) mounted on a Kapton
// board covered with a TPB-coated teflon mask.
//
// The NEXT Collaboration
// -----------------------------------------------------------------------------

#include "Next100AbsBoard.h"

#include "MaterialsList.h"
#include "GenericPhotosensor.h"
#include "OpticalMaterialProperties.h"
#include "BoxPointSampler.h"
#include "Visibilities.h"

#include <G4GenericMessenger.hh>
#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4Material.hh>
#include <G4NistManager.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4LogicalBorderSurface.hh>

using namespace nexus;


Next100AbsBoard::Next100AbsBoard():
  BaseGeometry     (),
  size_            (123.40 * mm),
  pitch_           ( 15.55 * mm),
  margin_          (  7.33 * mm),
  hole_diam_       (  7.00 * mm),
  board_thickness_ (  0.5  * mm),
  //  mask_thickness_  (  2.1  * mm), // Made slightly thicker to fit SiPM
  time_binning_    (1. * microsecond),
  visibility_      (true),
  sipm_visibility_ (false),
  mpv_             (nullptr),
  vtxgen_          (nullptr),
  sipm_            (new GenericPhotosensor("SiPM", 1.3 * mm))
{
  msg_ = new G4GenericMessenger(this, "/Geometry/Next100/",
                                "Control commands of the NEXT-100 geometry.");

  msg_->DeclareProperty("sipm_board_vis", visibility_,
                        "Visibility of Next100AbsBoard.");

  msg_->DeclareProperty("sipm_vis", sipm_visibility_,
                        "Visibility of Next100 SiPMs.");

  G4GenericMessenger::Command& time_binning_cmd =
  msg_->DeclareProperty("sipm_time_binning", time_binning_,
                        "TP SiPMs time binning.");
  time_binning_cmd.SetParameterName("sipm_time_binning", false);
  time_binning_cmd.SetUnitCategory("Time");
  time_binning_cmd.SetRange("sipm_time_binning>0.");
}


Next100AbsBoard::~Next100AbsBoard()
{
  delete msg_;
  delete vtxgen_;
  delete sipm_;
}


void Next100AbsBoard::Construct()
{
  // Make sure the mother physical volume is actually valid
  if (!mpv_)
    G4Exception("[Next100AbsBoard]", "Construct()",
                FatalException, "Mother volume is a nullptr.");

  G4Material* mother_gas = mpv_->GetLogicalVolume()->GetMaterial();

  // SILICON PHOTOMULTIPLIER (SIPM) //////////////////////////////////

  // We use for now the generic photosensor until the exact features
  // of the new Hamamatsu SiPMs are known.

  G4MaterialPropertiesTable* photosensor_mpt = new G4MaterialPropertiesTable();
  G4double energy[]       = {0.2 * eV, 11.5 * eV};
  G4double reflectivity[] = {0.0     ,  0.0};
  G4double efficiency[]   = {1.0     ,  1.0};
  photosensor_mpt->AddProperty("REFLECTIVITY", energy, reflectivity, 2);
  photosensor_mpt->AddProperty("EFFICIENCY",   energy, efficiency,   2);
  sipm_->SetVisibility(sipm_visibility_);
  sipm_->SetOpticalProperties(photosensor_mpt);
  sipm_->SetWithWLSCoating(true);
  sipm_->SetTimeBinning(time_binning_);
  sipm_->SetSensorDepth(2);
  sipm_->SetMotherDepth(4);
  sipm_->SetNamingOrder(1000);
  sipm_->Construct();

  G4String board_name = "SIPM_BOARD";

  // KAPTON BOARD ////////////////////////////////////////////////////
  // Wrapper volume that contains all other elements.

  mask_thickness_ = sipm_->GetThickness();
  G4cout << "SiPM thickness = " << sipm_->GetThickness() << G4endl;

  G4Box* board_solid_vol =
    new G4Box(board_name, size_/2., size_/2., (board_thickness_ + mask_thickness_)/2.);

  G4LogicalVolume* board_logic_vol =
    new G4LogicalVolume(board_solid_vol,
                        G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON"),
                        board_name);

  BaseGeometry::SetLogicalVolume(board_logic_vol);


  // TEFLON MASK /////////////////////////////////////////////////////

  G4String mask_name = "SIPM_BOARD_MASK";
  G4double mask_zpos = board_thickness_/2.;

  G4Box* mask_solid_vol =
    new G4Box(mask_name, size_/2., size_/2., mask_thickness_/2.);

  G4LogicalVolume* mask_logic_vol =
      new G4LogicalVolume(mask_solid_vol,
                          G4NistManager::Instance()->FindOrBuildMaterial("G4_TEFLON"),
                          mask_name);

  new G4PVPlacement(nullptr, G4ThreeVector(0., 0., mask_zpos),
                    mask_logic_vol, mask_name, board_logic_vol, false, 0, false);

  G4OpticalSurface* mask_opsurf =
    new G4OpticalSurface(mask_name+"_OPSURF", unified, ground, dielectric_metal);

  G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();
  const G4int REFL_NUMENTRIES = 6;
  G4double ENERGIES[REFL_NUMENTRIES] =
    {0.2*eV, 2.*eV, 4.*eV, 6.*eV, 8.*eV, 11.5*eV};
  G4double REFLECTIVITY[REFL_NUMENTRIES] = {0., 0., 0., 0., 0., 0.};
  mpt->AddProperty("REFLECTIVITY", ENERGIES, REFLECTIVITY, REFL_NUMENTRIES);
  mask_opsurf->SetMaterialPropertiesTable(mpt);
  new G4LogicalSkinSurface(mask_name+"_OPSURF", mask_logic_vol, mask_opsurf);


  // MASK GAS HOLE ///////////////////////////////////////////////////

  G4String mask_hole_name   = "SIPM_BOARD_MASK_HOLE";
  G4double mask_hole_length = mask_thickness_;
  G4double mask_hole_zpos   = - mask_thickness_/2. + mask_hole_length/2.;

  G4Tubs* mask_hole_solid_vol =
    new G4Tubs(mask_hole_name, 0., hole_diam_/2., mask_hole_length/2., 0, 360.*deg);

  G4LogicalVolume* mask_hole_logic_vol =
    new G4LogicalVolume(mask_hole_solid_vol, mother_gas, mask_hole_name);

  // (Placement of this volume below.)



  G4double sipm_zpos = - mask_hole_length/2. + sipm_->GetThickness()/2.;

  new G4PVPlacement(nullptr, G4ThreeVector(0., 0., sipm_zpos),
                    sipm_->GetLogicalVolume(), sipm_->GetLogicalVolume()->GetName(), mask_hole_logic_vol,
                    false, 0, false);


  ////////////////////////////////////////////////////////////////////

  // Placing now 8x8 replicas of the gas hole and SiPM

  G4double zpos = board_thickness_ + sipm_->GetThickness()/2.;

  G4int counter = 0;

  for (auto i=0; i<8; i++) {

    G4double xpos = -size_/2. + margin_ + i * pitch_;

    for (auto j=0; j<8; j++) {

      G4double ypos = -size_/2. + margin_ + j * pitch_;

      G4ThreeVector sipm_position(xpos, ypos, zpos);
      sipm_positions_.push_back(sipm_position);

      // Placement of the hole+SiPM
      new G4PVPlacement(nullptr, G4ThreeVector(xpos, ypos, mask_hole_zpos),
                        mask_hole_logic_vol, mask_hole_name, mask_logic_vol,
                        false, counter, false);

      counter++;
    }
  }

  // VERTEX GENERATOR ////////////////////////////////////////////////

  vtxgen_ = new BoxPointSampler(size_, size_, board_thickness_, 0.,
                                G4ThreeVector(0., 0., -mask_thickness_/2.));

  // VISIBILITIES ////////////////////////////////////////////////////
  if (visibility_) {
    G4VisAttributes blue       = Blue();
    G4VisAttributes light_blue = LightBlue();
    board_logic_vol ->SetVisAttributes(blue);
    mask_logic_vol  ->SetVisAttributes(light_blue);
  }
  else {
    board_logic_vol ->SetVisAttributes(G4VisAttributes::Invisible);
    mask_logic_vol  ->SetVisAttributes(G4VisAttributes::Invisible);
  }
  mask_hole_logic_vol    ->SetVisAttributes(G4VisAttributes::Invisible);
}



G4ThreeVector Next100AbsBoard::GenerateVertex(const G4String&) const
{
  // Only one generation region available at the moment
  return vtxgen_->GenerateVertex("INSIDE");
}
