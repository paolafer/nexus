// ----------------------------------------------------------------------------
// nexus | LXeMiniCell.cc
//
// Basic cell made of LXe.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "LXeMiniCell.h"
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

  LXeMiniCell::LXeMiniCell():
    BaseGeometry(),

    // Detector dimensions
    active_size_ (3.*mm),
    max_step_size_ (1.*mm),
    lxe_zsize_ (5.*mm)

  {
    // Messenger
    // msg_ = new G4GenericMessenger(this, "/Geometry/LXeMiniCell/",
    //                               "Control commands of geometry Petalo.");

    // z size
    //  G4GenericMessenger::Command& zsize_cmd =
    //    msg_->DeclareProperty("z_size", lxe_zsize_, "z dimension");
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


  LXeMiniCell::~LXeMiniCell()
  {
  }


  void LXeMiniCell::Construct()
  {
    G4Material* lxe = G4NistManager::Instance()->FindOrBuildMaterial("G4_lXe");
    lxe->SetMaterialPropertiesTable(OpticalMaterialProperties::LXe());
    G4Material* kapton =
      G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");

    G4double sipm_x = active_size_;
    G4double sipm_y = active_size_;
    G4double sipm_z = 1.55 * mm;
    G4double reflector_thickn = 0.1 * mm;

    G4double tot_xy_size = active_size_ + 2.*reflector_thickn;
    G4double tot_z_size  = lxe_zsize_ + sipm_z + reflector_thickn;

    G4Box* lxe_solid =
      new G4Box("LXE", tot_xy_size/2., tot_xy_size/2., tot_z_size/2.);
    G4LogicalVolume* lxe_logic = new G4LogicalVolume(lxe_solid, lxe, "LXE");
    this->SetLogicalVolume(lxe_logic);

    G4Box* refl1_solid =
      new G4Box("REFLECTOR", tot_xy_size/2., reflector_thickn/2., tot_z_size/2.);
    G4LogicalVolume* refl1_logic =
      new G4LogicalVolume(refl1_solid, kapton, "REFLECTOR");
    G4PVPlacement* refl1_phys =
      new G4PVPlacement(0, G4ThreeVector(0., tot_xy_size/2. - reflector_thickn/2., 0.),
                        refl1_logic, "REFLECTOR", lxe_logic, false, 0, true);
    G4PVPlacement* refl2_phys =
      new G4PVPlacement(0, G4ThreeVector(0., -tot_xy_size/2. + reflector_thickn/2., 0.),
                        refl1_logic, "REFLECTOR", lxe_logic, false, 1, true);

    G4Box* refl2_solid =
      new G4Box("REFLECTOR", reflector_thickn/2., active_size_/2., tot_z_size/2.);
    G4LogicalVolume* refl2_logic =
      new G4LogicalVolume(refl2_solid, kapton, "REFLECTOR");
    G4PVPlacement* refl3_phys =
      new G4PVPlacement(0, G4ThreeVector(tot_xy_size/2. - reflector_thickn/2., 0., 0.),
                        refl2_logic, "REFLECTOR", lxe_logic, false, 0, true);
    G4PVPlacement* refl4_phys =
      new G4PVPlacement(0, G4ThreeVector(-tot_xy_size/2. + reflector_thickn/2., 0., 0.),
                        refl2_logic, "REFLECTOR", lxe_logic, false, 1, true);

    G4Box* refl3_solid =
      new G4Box("REFLECTOR", active_size_/2., active_size_/2., reflector_thickn/2.);
    G4LogicalVolume* refl3_logic =
      new G4LogicalVolume(refl3_solid, kapton, "REFLECTOR");
    G4PVPlacement* refl5_phys =
      new G4PVPlacement(0, G4ThreeVector(0., 0., tot_z_size/2. - reflector_thickn/2.),
                        refl3_logic, "REFLECTOR", lxe_logic, false, 0, true);

    // Build and place the SiPM
    // No window (FBK style), support made of whatever
    G4Box* sipm_solid = new G4Box("SIPMpet", sipm_x/2., sipm_y/2., sipm_z/2);
    G4Material* sipm_mat = MaterialsList::FR4();
    G4LogicalVolume* sipm_logic =
      new G4LogicalVolume(sipm_solid, sipm_mat, "SIPMpet");

    G4double wndw_depth = 0.01 * mm;
    G4Box* wndw_solid =
      new G4Box("PHOTODIODES", active_size_/2., active_size_/2., wndw_depth/2);
    G4Material* silicon =
      G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
    G4LogicalVolume* wndw_logic =
      new G4LogicalVolume(wndw_solid, silicon, "PHOTODIODES");
    new G4PVPlacement(0, G4ThreeVector(0., 0., sipm_z/2. - wndw_depth/2.),
                      wndw_logic, "PHOTODIODES", sipm_logic, false, 0, true);


    const G4int n = 4;
    G4double sipm_energies[n]     = {1.5*eV, 4*eV, 6.*eV, 8.26558*eV};
    G4double sipm_reflectivity[n] = {0., 0., 0., 0.};
    G4double sipm_efficiency[n]   = {0.2, 0.2, 0.2, 0.2};

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
      sipmsd->SetDetectorVolumeDepth(3);
      // sipmsd->SetMotherVolumeDepth(2);
      // sipmsd->SetDetectorNamingOrder(1000.);
      sipmsd->SetTimeBinning(2. * microsecond);
      G4SDManager::GetSDMpointer()->AddNewDetector(sipmsd);
      wndw_logic->SetSensitiveDetector(sipmsd);
    }

    new G4PVPlacement(0, G4ThreeVector(0., 0., -tot_z_size/2. + sipm_z/2.),
                      sipm_logic, "SIPMpet", lxe_logic, false, 0, true);

    G4Box* active_solid =
      new G4Box("ACTIVE_LXE", active_size_/2., active_size_/2., lxe_zsize_/2.);
    G4LogicalVolume* active_logic =
      new G4LogicalVolume(active_solid, lxe, "ACTIVE_LXE");
    active_logic->SetVisAttributes(G4VisAttributes::Invisible);

    G4double lxe_zpos = tot_z_size/2. - reflector_thickn - lxe_zsize_/2.;
    G4PVPlacement* active_phys =
    new G4PVPlacement(0, G4ThreeVector(0., 0., lxe_zpos), active_logic,
		      "ACTIVE_LXE", lxe_logic, false, 0, true);

    active_logic->SetUserLimits(new G4UserLimits(max_step_size_));
    // Set the ACTIVE volume as an ionization sensitive active
    IonizationSD* ionisd = new IonizationSD("/PETALX/ACTIVE_LXE");
    active_logic->SetSensitiveDetector(ionisd);
    G4SDManager::GetSDMpointer()->AddNewDetector(ionisd);

    // Optical surfaces
    G4OpticalSurface* lxe_refl_surf =
      new G4OpticalSurface("LXE_REFL_OPSURF", unified, groundbackpainted,
                           dielectric_dielectric, 0.);

    new G4LogicalBorderSurface("LXE_AIR_OPSURF1", active_phys, refl1_phys,
                               lxe_refl_surf);
    new G4LogicalBorderSurface("LXE_AIR_OPSURF2", active_phys, refl2_phys,
                               lxe_refl_surf);
    new G4LogicalBorderSurface("LXE_AIR_OPSURF3", active_phys, refl3_phys,
                               lxe_refl_surf);
    new G4LogicalBorderSurface("LXE_AIR_OPSURF4", active_phys, refl4_phys,
                               lxe_refl_surf);
    new G4LogicalBorderSurface("LXE_AIR_OPSURF5", active_phys, refl5_phys,
                               lxe_refl_surf);

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

    lxe_refl_surf->SetMaterialPropertiesTable(smpt);

    G4VisAttributes sipm_col = nexus::Yellow();
    sipm_logic->SetVisAttributes(sipm_col);
    G4VisAttributes wndw_col = nexus::Red();
    wndw_col.SetForceSolid(true);
    wndw_logic->SetVisAttributes(wndw_col);

    G4VisAttributes active_col = nexus::Blue();
    active_col.SetForceSolid(true);
    active_logic->SetVisAttributes(active_col);
  }


  G4ThreeVector LXeMiniCell::GenerateVertex(const G4String& /*region*/) const
  {
    G4ThreeVector vertex(0.,0.,0.);

    return vertex;
  }


} //end namespace nexus
