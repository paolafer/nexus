// ----------------------------------------------------------------------------
// nexus | LYSOLab.cc
//
// This class consists of two cells placed opposite to each other.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "LYSOLab.h"
#include "LYSOCrystal.h"
#include "MaterialsList.h"
#include "OpticalMaterialProperties.h"
#include "IonizationSD.h"

#include <G4GenericMessenger.hh>
#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4SDManager.hh>
#include <G4NistManager.hh>
#include <G4VisAttributes.hh>


namespace nexus {

  using namespace CLHEP;

  LYSOLab::LYSOLab(): BaseGeometry(),
                      msg_(0),
                      n_crystal_rows_(1),
                      crystal_pitch_(3.5 * mm),
                      inner_radius_(380. * mm),
                      sec_size_(3. * mm)

  {
    msg_ = new G4GenericMessenger(this, "/Geometry/LYSOLab/",
				  "Control commands of geometry LYSOLab.");
    msg_->DeclareProperty("crystal_rows", n_crystal_rows_,
                          "Number of crystal rows");

    G4GenericMessenger::Command& pitch_cmd =
      msg_->DeclareProperty("pitch", crystal_pitch_, "Pitch of crystals");
    pitch_cmd.SetUnitCategory("Length");
    pitch_cmd.SetParameterName("pitch", false);
    pitch_cmd.SetRange("pitch>0.");

    G4GenericMessenger::Command& inner_r_cmd =
      msg_->DeclareProperty("inner_radius", inner_radius_,
                            "Inner radius of ring");
    inner_r_cmd.SetUnitCategory("Length");
    inner_r_cmd.SetParameterName("inner_radius", false);
    inner_r_cmd.SetRange("inner_radius>0.");

    G4GenericMessenger::Command& depth_cmd =
      msg_->DeclareProperty("depth", depth_, "Dimension in DOI");
    depth_cmd.SetUnitCategory("Length");
    depth_cmd.SetParameterName("depth", false);
    depth_cmd.SetRange("depth>0.");

    G4GenericMessenger::Command& sec_cmd =
      msg_->DeclareProperty("sec_size", sec_size_, "Dimension in xy");
    sec_cmd.SetUnitCategory("Length");
    sec_cmd.SetParameterName("sec_size", false);
    sec_cmd.SetRange("sec_size>0.");
  }



  LYSOLab::~LYSOLab()
  {
    delete msg_;
  }



  void LYSOLab::Construct()
  {
    // LAB //
    // This is just a volume of air surrounding the detector so that
    // events can be generated on the outside.

    G4double lab_size (2. * m);
    G4Box* lab_solid =
      new G4Box("LAB", lab_size/2., lab_size/2., lab_size/2.);

    G4LogicalVolume* lab_logic =
      new G4LogicalVolume(lab_solid, G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"), "AIR");
    lab_logic->SetVisAttributes(G4VisAttributes::Invisible);

    // Set this volume as the wrapper for the whole geometry
    // (i.e., this is the volume that will be placed in the world)
    this->SetLogicalVolume(lab_logic);


    lyso_module_ = new LYSOCrystal();
    lyso_module_->SetSecSize(sec_size_);
    lyso_module_->SetDOISize(depth_);
    lyso_module_->Construct();
    G4LogicalVolume* crystal_logic = lyso_module_->GetLogicalVolume();

    G4ThreeVector crystal_dim = lyso_module_->GetDimensions();

    G4double offset = 0.1 * mm;
    G4int n_crystals = 2*pi*inner_radius_/crystal_pitch_;
    G4cout << "Number of crystals: " <<  n_crystals * n_crystal_rows_ << G4endl;
    G4double step = 2.*pi/n_crystals;
    G4double radius = inner_radius_ + crystal_dim.z()/2.;
    G4double axial_length = crystal_pitch_ *  n_crystal_rows_;

    G4RotationMatrix rot;
    rot.rotateX(pi/2.);

    G4int copy_no = 0;
    for (G4int j=0; j<n_crystal_rows_; j++) {
      // The first must be positioned outside the loop
      if (j!=0) rot.rotateZ(step);
      G4double z = -axial_length/2. + (j + 0.5) * crystal_pitch_;
      G4ThreeVector position(0., radius, z);
      copy_no += 1;
      G4String vol_name = "CRYSTAL_" + std::to_string(copy_no);
      new G4PVPlacement(G4Transform3D(rot, position), crystal_logic,
                        vol_name, lab_logic, false, copy_no, true);

      for (G4int i=2; i<=n_crystals; ++i) {
        G4double angle = (i-1)*step;
        rot.rotateZ(step);
        position.setX(-radius*sin(angle));
        position.setY(radius*cos(angle));
        copy_no = copy_no + 1;
        vol_name = "CRYSTAL_" + std::to_string(copy_no);
        new G4PVPlacement(G4Transform3D(rot, position), crystal_logic,
                          vol_name, lab_logic, false, copy_no, false);
      }
    }

  }



  G4ThreeVector LYSOLab::GenerateVertex(const G4String& region) const
  {
    G4ThreeVector vertex(0., 0., 0.);

    if (region == "CENTER") {
      vertex = G4ThreeVector(0.,0.,0.);
    } else {
       G4Exception("[LYSOLab]", "GenerateVertex()", FatalException,
		  "Unknown vertex generation region!");
    }

     return vertex;
  }



} // end namespace nexus
