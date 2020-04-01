// ----------------------------------------------------------------------------
//  Author:  <paola.ferrario@dipc.org>
//  Created: 5 March 2020
// ----------------------------------------------------------------------------

#include "SiPMpet.h"
#include "PmtSD.h"
#include "ToFSD.h"
#include "MaterialsList.h"
#include "OpticalMaterialProperties.h"
#include "Visibilities.h"

#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4VisAttributes.hh>
#include <G4PVPlacement.hh>
#include <G4Material.hh>
#include <G4NistManager.hh>
#include <G4SDManager.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4PhysicalConstants.hh>
#include <G4GenericMessenger.hh>

#include <CLHEP/Units/SystemOfUnits.h>


namespace nexus {

  using namespace CLHEP;

  SiPMpet::SiPMpet(): BaseGeometry(),
		      visibility_(0),
		      refr_index_(1.54),
		      eff_(1.),
		      time_binning_(1.*microsecond),
		      sipm_size_(3.5*mm),
                      sipm_z_(0.6 * mm),
                      offset_(0.1*mm),
                      active_depth_(0.1*mm),
                      active_size_(3.*mm),
		      mother_depth_(0),
		      naming_order_(0),
		      wls_coating_(false),
		      decay_time_(2.2*nanosecond),
		      qe_(0.8)
  {
    /// Messenger
    msg_ = new G4GenericMessenger(this, "/Geometry/SiPMpet/",
				  "Control commands of geometry.");
    msg_->DeclareProperty("visibility", visibility_, "SiPM visibility");
    msg_->DeclareProperty("refr_index", refr_index_, "Refraction index for epoxy");
    msg_->DeclareProperty("efficiency", eff_, "Detection efficiency of SiPM");

    G4GenericMessenger::Command& time_cmd =
      msg_->DeclareProperty("time_binning", time_binning_,
			    "Time binning for the sensor");
    time_cmd.SetUnitCategory("Time");
    time_cmd.SetParameterName("time_binning", false);
    time_cmd.SetRange("time_binning>0.");

    G4GenericMessenger::Command& size_cmd =
      msg_->DeclareProperty("size", sipm_size_, "Size of SiPM");
    size_cmd.SetUnitCategory("Length");
    size_cmd.SetParameterName("size", false);
    size_cmd.SetRange("size>0.");

    G4GenericMessenger::Command& active_size_cmd =
      msg_->DeclareProperty("active_size", active_size_, "Size of SiPM active area");
    active_size_cmd.SetUnitCategory("Length");
    active_size_cmd.SetParameterName("active_size", false);
    active_size_cmd.SetRange("active_size>0.");

    msg_->DeclareProperty("mother_depth", mother_depth_,
			  "Depth of sipm in its mother, being replicated");
    msg_->DeclareProperty("naming_order", naming_order_,
			  "To start numbering at different place than zero.");

    msg_->DeclareProperty("wls_coating", wls_coating_,
			  "True if a wavelength shifting coating is applied");

    G4GenericMessenger::Command& decay_time_cmd =
      msg_->DeclareProperty("decay_time", decay_time_,
			    "Decay time of WLS");
    decay_time_cmd.SetUnitCategory("Time");
    decay_time_cmd.SetParameterName("decay_time", false);
    decay_time_cmd.SetRange("decay_time>0.");

    msg_->DeclareProperty("wls_qe", qe_, "Quantum efficiency of WLS");
  }



  SiPMpet::~SiPMpet()
  {
  }



  void SiPMpet::Construct()
  {
    SetDimensions(G4ThreeVector(sipm_size_, sipm_size_, sipm_z_));

    G4String name = "SIPMpet";
    G4Box* sipm_solid = new G4Box(name, sipm_size_/2., sipm_size_/2., sipm_z_/2);

    G4Material* epoxy = MaterialsList::Epoxy();
    if (refr_index_ > 0) {
      epoxy->SetMaterialPropertiesTable(OpticalMaterialProperties::EpoxyFixedRefr(refr_index_));
    } else {
      epoxy->SetMaterialPropertiesTable(OpticalMaterialProperties::EpoxyLXeRefr());
    }

    G4LogicalVolume* sipm_logic = new G4LogicalVolume(sipm_solid, epoxy, name);

    this->SetLogicalVolume(sipm_logic);

    G4double pos_z;

    if (wls_coating_) {
      name = "WLS";
      G4double tpb_z = 0.001 * mm;
      G4Box* tpb_solid = new G4Box(name, sipm_size_/2., sipm_size_/2., tpb_z/2);
      G4Material* TPB = MaterialsList::TPB();
      TPB->SetMaterialPropertiesTable(OpticalMaterialProperties::GenericWLS(decay_time_,
									    qe_));
      G4LogicalVolume* tpb_logic = new G4LogicalVolume(tpb_solid, TPB, name);

      pos_z = (sipm_z_ - tpb_z) / 2.;
      new G4PVPlacement(0, G4ThreeVector(0.,0.,pos_z), tpb_logic,
                        name, sipm_logic, false, 0, true);

    }

    // PLASTIC SUPPORT //

    name = "SIPM_SUPPORT";
    G4double support_depth = sipm_z_ - offset_;

    G4Box* support_solid =
      new G4Box(name, sipm_size_/2., sipm_size_/2., support_depth/2);

    G4Material* plastic =
      G4NistManager::Instance()->FindOrBuildMaterial("G4_POLYCARBONATE");

    G4LogicalVolume* support_logic =
      new G4LogicalVolume(support_solid, plastic, name);

    pos_z = -sipm_z_/2. + support_depth/2.;
    new G4PVPlacement(0, G4ThreeVector(0., 0., pos_z),
		      support_logic, name, sipm_logic, false, 0, true);

    // ACTIVE REGION //

    name = "PHOTODIODES";

    G4Box* active_solid =
      new G4Box(name, active_size_/2., active_size_/2., active_depth_/2);

    G4Material* silicon =
      G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");

    G4LogicalVolume* active_logic =
      new G4LogicalVolume(active_solid, silicon, name);

    pos_z = support_depth/2. - active_depth_/2.;
    new G4PVPlacement(0, G4ThreeVector(0., 0., pos_z),
		      active_logic, name, support_logic, false, 0, true);


    // OPTICAL SURFACES //////////////////////////////////////////////

    const G4int entries            = 2;
    G4double energies[entries]     = {1.54980*eV, 8.26558*eV};
    G4double reflectivity[entries] = {0. ,0. };
    G4double efficiency[entries]   = {eff_, eff_};

    G4MaterialPropertiesTable* sipm_mpt = new G4MaterialPropertiesTable();
    sipm_mpt->AddProperty("EFFICIENCY", energies, efficiency, entries);
    sipm_mpt->AddProperty("REFLECTIVITY", energies, reflectivity, entries);

    G4OpticalSurface* sipm_opsurf =
      new G4OpticalSurface("SIPM_OPSURF", unified, polished, dielectric_metal);
    sipm_opsurf->SetMaterialPropertiesTable(sipm_mpt);

    new G4LogicalSkinSurface("SIPM_OPSURF", active_logic, sipm_opsurf);


    // SENSITIVE DETECTOR ////////////////////////////////////////////

    G4String sdname = "/SIPM/SiPMpet";
    G4SDManager* sdmgr = G4SDManager::GetSDMpointer();

    if (!sdmgr->FindSensitiveDetector(sdname, false)) {
      ToFSD* sipmsd = new ToFSD(sdname);
      sipmsd->SetDetectorNamingOrder(naming_order_);
      sipmsd->SetMotherVolumeDepth(mother_depth_);
      sipmsd->SetTimeBinning(time_binning_);

      G4SDManager::GetSDMpointer()->AddNewDetector(sipmsd);
      sipm_logic->SetSensitiveDetector(sipmsd);
    }

    // Visibilities
    if (visibility_) {
      G4VisAttributes sipm_col = nexus::Brown();
      sipm_logic->SetVisAttributes(sipm_col);
      G4VisAttributes active_col = nexus::Red();
      active_col.SetForceSolid(true);
      active_logic->SetVisAttributes(active_col);
      G4VisAttributes support_col = nexus::Yellow();
      support_logic->SetVisAttributes(support_col);
    }
    else {
      active_logic->SetVisAttributes(G4VisAttributes::Invisible);
      support_logic->SetVisAttributes(G4VisAttributes::Invisible);
    }
  }


} // end namespace nexus
