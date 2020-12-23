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
    BaseGeometry(), msg_(0), mat_("lyso"), reflectivity_(0.95), surf_n_(1)
  {
    msg_ = new G4GenericMessenger(this, "/Geometry/Lab/",
				  "Control commands of geometry Lab.");
    msg_->DeclareProperty("material", mat_,
                          "Material of modules.");
    msg_->DeclareProperty("reflectivity", reflectivity_,
                          "Reflectivity if walls.");
    msg_->DeclareProperty("surface_refr_index", surf_n_,
                          "Refractive index of layer between surface and reflector.");
    
    G4GenericMessenger::Command& z_size_cmd =
        msg_->DeclareProperty("active_z_size", active_z_,
                              "Z dimension of crystal/cell");
                              z_size_cmd.SetUnitCategory("Length");
                              z_size_cmd.SetParameterName("active_z_size", false);
                              z_size_cmd.SetRange("active_z_size>0.");

    G4GenericMessenger::Command& z_cmd =
        msg_->DeclareProperty("interaction_z", z_,
                              "Z pos of gamma interaction");
                              z_cmd.SetUnitCategory("Length");
                              z_cmd.SetParameterName("interaction_z", false);
                              z_cmd.SetRange("interaction_z>0.");

    G4GenericMessenger::Command& vx_cmd =
        msg_->DeclareProperty("specific_vertex_X", vx_,
                              "If region is AD_HOC, x coord of primary particles");
                              vx_cmd.SetUnitCategory("Length");
                              vx_cmd.SetParameterName("specific_vertex_X", false);

    G4GenericMessenger::Command& vy_cmd =
        msg_->DeclareProperty("specific_vertex_Y", vy_,
                              "If region is AD_HOC, y coord of primary particles");
                                vy_cmd.SetUnitCategory("Length");
                                vy_cmd.SetParameterName("specific_vertex_Y", false);

    G4GenericMessenger::Command& vz_cmd =
        msg_->DeclareProperty("specific_vertex_Z", vz_,
                              "If region is AD_HOC, z coord of primary particles");
                              vz_cmd.SetUnitCategory("Length");
                              vz_cmd.SetParameterName("specific_vertex_Z", false);
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

    G4double z_dist_entry_surf = 3.*cm;

    G4LogicalVolume* module_logic = nullptr;
    G4double z_dim = 0.;
    if (mat_ == "lyso") {
        lyso_module_ = new LYSOCrystal();
        lyso_module_->Construct();
        z_dim = lyso_module_->GetDimensions().z();
        module_logic = lyso_module_->GetLogicalVolume();
    } else if (mat_ == "lxe"){
        lxe_module_ = new LXeMiniCell();
        lxe_module_->SetDOISize(active_z_);
        lxe_module_->SetReflectivity(reflectivity_);
        lxe_module_->Construct();
        z_dim = lxe_module_->GetDimensions().z();
        module_logic = lxe_module_->GetLogicalVolume();
    } else if (mat_ == "lso") {
        lso_module_ = new LSOCeCaCrystal();
        lso_module_->SetDOISize(active_z_);
        lso_module_->SetReflectivity(reflectivity_);
        lso_module_->Construct();
        z_dim = lso_module_->GetDimensions().z();
        module_logic = lso_module_->GetLogicalVolume();

        generator1_ =
            new BoxPointSampler(lso_module_->GetSecSize(),
                                lso_module_->GetSecSize(),
                                1.*mm, 0,
                                G4ThreeVector(0., 0., -z_dist_entry_surf-z_), 0);
        generator2_ =
            new BoxPointSampler(lso_module_->GetSecSize(),
                                lso_module_->GetSecSize(),
                                1.*mm, 0,
                                G4ThreeVector(0., 0., z_dist_entry_surf+z_), 0);
    }


    new G4PVPlacement(0, G4ThreeVector(0., 0., -z_dist_entry_surf-z_dim/2.), module_logic, "MODULE_0",
        lab_logic, false, 1, true);

    G4RotationMatrix rot;
    rot.rotateY(pi);
    new G4PVPlacement(G4Transform3D(rot, G4ThreeVector(0., 0., z_dist_entry_surf+z_dim/2.)), module_logic,
                      "MODULE_1", lab_logic, false, 2, true);

  }



  G4ThreeVector Lab::GenerateVertex(const G4String& region) const
  {
    G4ThreeVector vertex(0., 0., 0.);

    if (region == "CENTER") {
      vertex = G4ThreeVector(0.,0.,0.);
    } else if (region == "AD_HOC") {
      vertex = G4ThreeVector(vx_, vy_, vz_);
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
