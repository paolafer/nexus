// ----------------------------------------------------------------------------
// nexus | LYSOCrystal.cc
//
// Basic cell made of LYSO.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "LYSOCrystal.h"
#include "MaterialsList.h"
#include "IonizationSD.h"
#include "PetKDBFixedPitch.h"
#include "PetPlainDice.h"
#include "OpticalMaterialProperties.h"
#include "BoxPointSampler.h"
#include "ToFSD.h"
#include "Visibilities.h"

#include <G4GenericMessenger.hh>
#include <G4Box.hh>
#include <G4Material.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4VisAttributes.hh>
#include <G4UserLimits.hh>
#include <G4NistManager.hh>
#include <G4SDManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4RotationMatrix.hh>


namespace nexus {

  LYSOCrystal::LYSOCrystal():
    BaseGeometry(),

    // Detector dimensions
    active_size_ (3.*mm),
    max_step_size_ (1.*mm),
    lyso_zsize_ (5.*mm)

  {
    // Messenger
    // msg_ = new G4GenericMessenger(this, "/Geometry/LYSOCrystal/",
    //                               "Control commands of geometry Petalo.");

    // // z size
    //  G4GenericMessenger::Command& zsize_cmd =
    //    msg_->DeclareProperty("z_size", lyso_zsize_, "z dimension");
    //  zsize_cmd.SetUnitCategory("Length");
    //  zsize_cmd.SetParameterName("z_size", false);
    //  zsize_cmd.SetRange("z_size>0.");

    // // Maximum Step Size
    // G4GenericMessenger::Command& step_cmd =
    //   msg_->DeclareProperty("max_step_size", max_step_size_,
    //                         "Maximum step size");
    // step_cmd.SetUnitCategory("Length");
    // step_cmd.SetParameterName("max_step_size", false);
    // step_cmd.SetRange("max_step_size>0.");
  }


  LYSOCrystal::~LYSOCrystal()
  {
  }


  void LYSOCrystal::Construct()
  {
    G4Material* lyso = MaterialsList::LYSO();
    lyso->SetMaterialPropertiesTable(OpticalMaterialProperties::LYSO());
    G4Material* kapton =
      G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");

    G4Material* resin = MaterialsList::Epoxy(); // what matters is n
    // resin->SetMaterialPropertiesTable(OpticalMaterialProperties::EpoxyFixedRefr(1.49));
    resin->SetMaterialPropertiesTable(OpticalMaterialProperties::EpoxyFixedRefr(1.54));
    G4Material* opt_gel = MaterialsList::OpticalSilicone(); // what matters is n
    opt_gel->SetMaterialPropertiesTable(OpticalMaterialProperties::EpoxyFixedRefr(1.6));

    G4double sipm_x = active_size_;
    G4double sipm_y = active_size_;
    G4double sipm_z = 1.55 * mm;
    G4double reflector_thickn = 0.1 * mm;
    G4double opt_gel_thickn = 0. * mm;

    G4double tot_xy_size = active_size_ + 2.*reflector_thickn;
    G4double tot_z_size  = lyso_zsize_ + sipm_z + opt_gel_thickn +
      reflector_thickn;

    // G4double refl_xy = tot_xy_size - 2.*container_thickn;
    // G4double refl_z = tot_z_size - 2.*container_thickn;
    // G4Box* refl_solid =
    //   new G4Box("REFLECTOR", refl_xy/2., refl_xy/2., refl_z/2.);
    // G4LogicalVolume* refl_logic =
    //   new G4LogicalVolume(refl_solid, kapton, "REFLECTOR");
    // G4PVPlacement* refl_phys =
    //   new G4PVPlacement(0, G4ThreeVector(0., 0., 0.),
    //                     refl_logic, "REFLECTOR", container_logic, false, 0, true);

    //G4double internal_z_size = lyso_zsize_ + sipm_z + opt_gel_thickn;
    G4Box* lyso_solid =
      new G4Box("LYSO", tot_xy_size/2., tot_xy_size/2., tot_z_size/2.);
      // new G4Box("LYSO", active_size_/2., active_size_/2., internal_z_size/2.);
    G4LogicalVolume* lyso_logic = new G4LogicalVolume(lyso_solid, lyso, "LYSO");
     // G4PVPlacement* lyso_phys =
     //   new G4PVPlacement(0, G4ThreeVector(0., 0., 0.),
     //                     lyso_logic, "LYSO", refl_logic, false, 0, true);
    this->SetLogicalVolume(lyso_logic);

    G4Box* refl1_solid =
      new G4Box("REFLECTOR", tot_xy_size/2., reflector_thickn/2., tot_z_size/2.);
    G4LogicalVolume* refl1_logic =
      new G4LogicalVolume(refl1_solid, kapton, "REFLECTOR");
   G4PVPlacement* refl1_phys =
     new G4PVPlacement(0, G4ThreeVector(0., tot_xy_size/2. - reflector_thickn/2., 0.),
                      refl1_logic, "REFLECTOR", lyso_logic, false, 0, true);
   G4PVPlacement* refl2_phys =
    new G4PVPlacement(0, G4ThreeVector(0., -tot_xy_size/2. + reflector_thickn/2., 0.),
                      refl1_logic, "REFLECTOR", lyso_logic, false, 1, true);

   G4Box* refl2_solid =
     new G4Box("REFLECTOR", reflector_thickn/2., active_size_/2., tot_z_size/2.);
   G4LogicalVolume* refl2_logic =
     new G4LogicalVolume(refl2_solid, kapton, "REFLECTOR");
   G4PVPlacement* refl3_phys =
     new G4PVPlacement(0, G4ThreeVector(tot_xy_size/2. - reflector_thickn/2., 0., 0.),
                       refl2_logic, "REFLECTOR", lyso_logic, false, 0, true);
   G4PVPlacement* refl4_phys =
     new G4PVPlacement(0, G4ThreeVector(-tot_xy_size/2. + reflector_thickn/2., 0., 0.),
                       refl2_logic, "REFLECTOR", lyso_logic, false, 1, true);

   G4Box* refl3_solid =
      new G4Box("REFLECTOR", active_size_/2., active_size_/2., reflector_thickn/2.);
    G4LogicalVolume* refl3_logic =
      new G4LogicalVolume(refl3_solid, kapton, "REFLECTOR");
   G4PVPlacement* refl5_phys =
     new G4PVPlacement(0, G4ThreeVector(0., 0., tot_z_size/2. - reflector_thickn/2.),
                       refl3_logic, "REFLECTOR", lyso_logic, false, 0, true);

    /*
    G4OpticalSurface* lyso_refl_surf =
           new G4OpticalSurface("LYSO_REFL_OPSURF", glisur, ground,
                                dielectric_metal, .01);
    lyso_refl_surf->SetMaterialPropertiesTable(OpticalMaterialProperties::ReflectantSurface(0.95));
    new G4LogicalSkinSurface("LYSO_REFL_OPSURF", refl1_logic, lyso_refl_surf);
    new G4LogicalSkinSurface("LYSO_REFL_OPSURF", refl2_logic, lyso_refl_surf);
    new G4LogicalSkinSurface("LYSO_REFL_OPSURF", refl3_logic, lyso_refl_surf);
    */


    // Build and place the SiPM
    G4Box* sipm_solid = new G4Box("SIPMpet", sipm_x/2., sipm_y/2., sipm_z/2);
    G4LogicalVolume* sipm_logic =
      new G4LogicalVolume(sipm_solid, resin, "SIPMpet");

    G4double wndw_depth = 0.01 * mm;
    G4double offset = 0.01 * mm;
    G4Box* wndw_solid =
      new G4Box("PHOTODIODES", active_size_/2., active_size_/2., wndw_depth/2);
    G4Material* silicon =
      G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
    G4LogicalVolume* wndw_logic =
      new G4LogicalVolume(wndw_solid, silicon, "PHOTODIODES");
    new G4PVPlacement(0, G4ThreeVector(0., 0., sipm_z/2. - wndw_depth/2. - offset),
                      wndw_logic, "PHOTODIODES", sipm_logic, false, 0, true);


    const G4int n = 4;
    G4double sipm_energies[n]     = {1.5*eV, 4*eV, 6.*eV, 8.26558*eV};
    G4double sipm_reflectivity[n] = {0., 0., 0., 0.};
    G4double sipm_efficiency[n]   = {0.45, 0.45, 0.45, 0.45};

    G4MaterialPropertiesTable* sipm_mt = new G4MaterialPropertiesTable();
    sipm_mt->AddProperty("EFFICIENCY", sipm_energies, sipm_efficiency, n);
    sipm_mt->AddProperty("REFLECTIVITY", sipm_energies, sipm_reflectivity, n);

    G4OpticalSurface* sipm_opsurf =
      new G4OpticalSurface("SIPM_OPSURF", unified, polished, dielectric_metal);
    sipm_opsurf->SetMaterialPropertiesTable(sipm_mt);
    new G4LogicalSkinSurface("SIPM_OPSURF", wndw_logic, sipm_opsurf);

    // sensitive detector
    G4String sdname = "/SIPM/SiPMpet";
    G4SDManager* sdmgr = G4SDManager::GetSDMpointer();

    if (!sdmgr->FindSensitiveDetector(sdname, false)) {
      ToFSD* sipmsd = new ToFSD(sdname);
      sipmsd->SetDetectorVolumeDepth(2); // with reflector and container: 4
      // sipmsd->SetMotherVolumeDepth(2);
      // sipmsd->SetDetectorNamingOrder(1000.);
      sipmsd->SetTimeBinning(2. * microsecond);
      G4SDManager::GetSDMpointer()->AddNewDetector(sipmsd);
      wndw_logic->SetSensitiveDetector(sipmsd);
    }

    new G4PVPlacement(0, G4ThreeVector(0., 0., -tot_z_size/2. + sipm_z/2.),
                      sipm_logic, "SIPMpet", lyso_logic, false, 0, true);


    // G4Box* opt_gel_solid =
    //   new G4Box("OPTICAL_GEL", active_size_/2., active_size_/2., opt_gel_thickn/2.);
    // G4LogicalVolume* opt_gel_logic =
    //   new G4LogicalVolume(opt_gel_solid, opt_gel, "OPTICAL_GEL");
    //    new G4PVPlacement(0, G4ThreeVector(0., 0., -internal_z_size/2. + sipm_z + opt_gel_thickn/2.),
    //                     opt_gel_logic, "OPTICAL_GEL", lyso_logic, false, 0, true);


    G4Box* active_solid =
      new G4Box("ACTIVE_LYSO", active_size_/2., active_size_/2., lyso_zsize_/2.);
    G4LogicalVolume* active_logic =
      new G4LogicalVolume(active_solid, lyso, "ACTIVE_LYSO");
    active_logic->SetVisAttributes(G4VisAttributes::Invisible);

    G4double crystal_zpos = tot_z_size/2. - reflector_thickn  - lyso_zsize_/2.;
    G4PVPlacement* active_phys =
      new G4PVPlacement(0, G4ThreeVector(0., 0., crystal_zpos), active_logic,
		      "ACTIVE_LYSO", lyso_logic, false, 0, true);

    active_logic->SetUserLimits(new G4UserLimits(max_step_size_));
    // Set the ACTIVE volume as an ionization sensitive active
    IonizationSD* ionisd = new IonizationSD("/PETALX/ACTIVE_LYSO");
    active_logic->SetSensitiveDetector(ionisd);
    G4SDManager::GetSDMpointer()->AddNewDetector(ionisd);


    G4OpticalSurface* lyso_refl_surf =
      new G4OpticalSurface("LYSO_REFL_OPSURF", unified, groundbackpainted,
                           dielectric_dielectric, 0.);

    new G4LogicalBorderSurface("LYSO_AIR_OPSURF1", active_phys, refl1_phys,
                               lyso_refl_surf);
    new G4LogicalBorderSurface("LYSO_AIR_OPSURF2", active_phys, refl2_phys,
                               lyso_refl_surf);
    new G4LogicalBorderSurface("LYSO_AIR_OPSURF3", active_phys, refl3_phys,
                               lyso_refl_surf);
    new G4LogicalBorderSurface("LYSO_AIR_OPSURF4", active_phys, refl4_phys,
                               lyso_refl_surf);
    new G4LogicalBorderSurface("LYSO_AIR_OPSURF5", active_phys, refl5_phys,
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


    G4VisAttributes sipm_col = nexus::Yellow();
    sipm_logic->SetVisAttributes(sipm_col);
    G4VisAttributes wndw_col = nexus::Red();
    wndw_col.SetForceSolid(true);
    wndw_logic->SetVisAttributes(wndw_col);

    G4VisAttributes active_col = nexus::Blue();
    active_col.SetForceSolid(true);
    active_logic->SetVisAttributes(active_col);
  }


  G4ThreeVector LYSOCrystal::GenerateVertex(const G4String& /*region*/) const
  {
    G4ThreeVector vertex(0.,0.,0.);

    return vertex;
  }


} //end namespace nexus
