// ----------------------------------------------------------------------------
// nexus | Lab.cc
//
// This class consists of two cells placed opposite to each other.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "Lab.h"

#include "Pet2boxes.h"
#include "PetLXeCell.h"
#include "LYSOCrystal.h"
#include "LSOCeCaCrystal.h"
#include "LXeMiniCell.h"
#include "PetLYSOCell.h"
#include "MaterialsList.h"
#include "OpticalMaterialProperties.h"
#include "IonizationSD.h"
#include "BoxPointSampler.h"

#include <G4GenericMessenger.hh>
#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4SDManager.hh>
#include <G4NistManager.hh>
#include <G4VisAttributes.hh>


namespace nexus {

  using namespace CLHEP;

  Lab::Lab():
    BaseGeometry(), msg_(0), mat_("lyso")
  {
    msg_ = new G4GenericMessenger(this, "/Geometry/Lab/",
				  "Control commands of geometry Lab.");
    msg_->DeclareProperty("material", mat_,
                          "Material of modules.");
    G4GenericMessenger::Command& z_cmd =
        msg_->DeclareProperty("interaction_z", z_,
                              "Z pos of gamma interaction");
                              z_cmd.SetUnitCategory("Length");
                              z_cmd.SetParameterName("interaction_z", false);
                              z_cmd.SetRange("interaction_z>0.");
  }



  Lab::~Lab()
  {
    delete msg_;
  }



  void Lab::Construct()
  {
    // LAB /////////////////////////////////////////////////////////////
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

    G4LogicalVolume* module_logic = nullptr;
    if (mat_ == "lyso") {
        lyso_module_ = new LYSOCrystal();
        lyso_module_->Construct();
        module_logic = lyso_module_->GetLogicalVolume();
    } else if (mat_ == "lxe"){
        lxe_module_ = new LXeMiniCell();
        lxe_module_->Construct();
        module_logic = lxe_module_->GetLogicalVolume();
    } else if (mat_ == "lso") {
        lso_module_ = new LSOCeCaCrystal();
        lso_module_->Construct();
        module_logic = lso_module_->GetLogicalVolume();

        generator1_ =
            new BoxPointSampler(lso_module_->GetSecSize(),
                                lso_module_->GetSecSize(),
                                1.*mm, 0,
                                G4ThreeVector(0., 0., -z_), 0);
        generator2_ =
            new BoxPointSampler(lso_module_->GetSecSize(),
                                lso_module_->GetSecSize(),
                                1.*mm, 0,
                                G4ThreeVector(0., 0., z_), 0);
    }

    G4double z_placement = 3.*cm;
    new G4PVPlacement(0, G4ThreeVector(0., 0., -z_placement), module_logic, "MODULE_0",
        lab_logic, false, 1, true);

    G4RotationMatrix rot;
    rot.rotateY(pi);
    new G4PVPlacement(G4Transform3D(rot, G4ThreeVector(0., 0., z_placement)), module_logic,
                      "MODULE_1", lab_logic, false, 2, true);

  }



  G4ThreeVector Lab::GenerateVertex(const G4String& region) const
  {
    G4ThreeVector vertex(0., 0., 0.);

    if (region == "CENTER") {
      vertex = G4ThreeVector(0.,0.,0.);
    } else {
       G4Exception("[Lab]", "GenerateVertex()", FatalException,
		  "Unknown vertex generation region!");
    }

     return vertex;
  }

  std::pair<G4ThreeVector, G4ThreeVector> Lab::GenerateVertices(const G4String& /*region*/) const
  {
    auto vertex1 = generator1_->GenerateVertex("INSIDE");
    auto vertex2 = generator2_->GenerateVertex("INSIDE");

    auto vertices = std::make_pair(vertex1, vertex2);

    return vertices;
  }



} // end namespace nexus
