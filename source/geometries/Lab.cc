// ----------------------------------------------------------------------------
// nexus | Lab.cc
//
// This class consists of two LXe cells placed opposite to each other.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "Lab.h"

#include "Pet2boxes.h"
#include "PetLXeCell.h"
#include "LYSOCrystal.h"
#include "LXeMiniCell.h"
#include "PetLYSOCell.h"
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
#include <G4LogicalBorderSurface.hh>
#include <G4OpticalSurface.hh>

#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>
#include <stdexcept>

#include <TTree.h>


namespace nexus {

  using namespace CLHEP;

  Lab::Lab():
    BaseGeometry(), msg_(0), lyso_(true)
  {
    msg_ = new G4GenericMessenger(this, "/Geometry/Lab/",
				  "Control commands of geometry Lab.");
    msg_->DeclareProperty("lyso", lyso_,
                          "True if lyso crystals, false if LXe mini cells.");
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

    G4double lab_size = 30*mm;
    G4Box* lab_solid =
      new G4Box("LAB", lab_size/2., lab_size/2., lab_size/2.);

    G4LogicalVolume* lab_logic =
      new G4LogicalVolume(lab_solid,
                          G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"), "LAB");
    lab_logic->SetVisAttributes(G4VisAttributes::Invisible);

    // Set this volume as the wrapper for the whole geometry
    // (i.e., this is the volume that will be placed in the world)
    this->SetLogicalVolume(lab_logic);

    G4double air_size = 29.*mm;
    G4Box* air_solid =
      new G4Box("AIR", air_size/2., air_size/2., air_size/2.);
    G4LogicalVolume* air_logic =
      new G4LogicalVolume(air_solid,
                          G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"), "AIR");
    G4PVPlacement*  air_phys =
      new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), air_logic, "AIR",
                        lab_logic, false, 0, true);

    if (lyso_) {
      lyso_module_ = new LYSOCrystal();
    } else {
      lxe_module_ = new LXeMiniCell();
    }

    G4LogicalVolume* module_logic = nullptr;
     if (lyso_) {
      lyso_module_->SetMotherPhysicalVolume(air_phys);
      lyso_module_->Construct();
      module_logic = lyso_module_->GetLogicalVolume();
    } else {
      lxe_module_->Construct();
      module_logic = lxe_module_->GetLogicalVolume();
    }

     G4PVPlacement* mod1_phys  =
       new G4PVPlacement(0, G4ThreeVector(0., 0., -1.*cm), module_logic, "MODULE_0",
                         air_logic, false, 1, true);

    G4RotationMatrix rot;
    rot.rotateY(pi);

    G4PVPlacement* mod2_phys  =
      new G4PVPlacement(G4Transform3D(rot, G4ThreeVector(0., 0., 1*cm)),
                        module_logic, "MODULE_1", air_logic, false, 2, true);

    
    // Optical surfaces around crystals
    G4OpticalSurface* lyso_refl_surf =
      new G4OpticalSurface("LYSO_AIR_OPSURF"); //, unified, groundbackpainted,
    //dielectric_dielectric, 0.);

    G4double sigma_alpha = 0.;

    lyso_refl_surf->SetType(dielectric_dielectric);
    lyso_refl_surf->SetModel(unified);
    lyso_refl_surf->SetFinish(groundbackpainted);
    lyso_refl_surf->SetSigmaAlpha(sigma_alpha);

    new G4LogicalBorderSurface("LYSO_AIR_OPSURF", mod1_phys, air_phys,
                               lyso_refl_surf);
    new G4LogicalBorderSurface("LYSO_AIR_OPSURF", mod2_phys, air_phys,
                               lyso_refl_surf);

    const G4int entries = 2;
    G4double energies[entries]      = {1.*eV, 8.*eV};
    G4double specularlobe[entries]  = {0., 0.};
    G4double specularspike[entries] = {0., 0.};
    G4double backscatter[entries]   = {0., 0.};
    G4double rindex[entries]        = {1., 1.}; // that of air.
    G4double reflectivity[entries]  = {.95, .95};
    G4double efficiency[entries]    = {0., 0.};
    

    G4MaterialPropertiesTable* smpt = new G4MaterialPropertiesTable();

    smpt->AddProperty("RINDEX", energies, rindex, entries);
    smpt->AddProperty("SPECULARLOBECONSTANT", energies, specularlobe, entries);
    smpt->AddProperty("SPECULARSPIKECONSTANT", energies, specularspike, entries);
    smpt->AddProperty("BACKSCATTERCONSTANT", energies, backscatter, entries);
    smpt->AddProperty("REFLECTIVITY", energies, reflectivity, entries);
    smpt->AddProperty("EFFICIENCY", energies, efficiency, entries);
   
    lyso_refl_surf->SetMaterialPropertiesTable(smpt);
    
    /*
    G4OpticalSurface* lyso_refl_surf =
      new G4OpticalSurface("LYSO_AIR_OPSURF", glisur, ground,
                           dielectric_metal, 0.);

    G4MaterialPropertiesTable* smpt = new G4MaterialPropertiesTable();
    const G4int entries = 4;
    G4double energies[entries]      = {1.*eV, 4*eV, 6.*eV, 8.*eV};
    G4double reflectivity[entries]  = {0.95, 0.95, 0.95, 0.95};
    smpt->AddProperty("REFLECTIVITY", energies, reflectivity, entries);
    */
    
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



} // end namespace nexus
