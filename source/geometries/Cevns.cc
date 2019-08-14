#include "Cevns.h"

#include <G4Box.hh>
#include <G4Tubs.hh>

#include <G4GenericMessenger.hh>
#include <G4NistManager.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <Randomize.hh>
#include <G4VisAttributes.hh>
#include <CLHEP/Units/SystemOfUnits.h>

namespace nexus {
  
  using namespace CLHEP;

  Cevns::Cevns(): _ess(1)
  {
    /// Messenger
    _msg = new G4GenericMessenger(this, "/Geometry/Cevns/",
                                  "Control commands of geometry Cevns.");

    _msg->DeclareProperty("ess", _ess,
                          "1 if ESS, 0 if SNS");

  }

  Cevns::~Cevns()
  {
  }

 void Cevns::Construct()
  {
    G4double size = 10. * m;
    G4Box* lab_solid =
      new G4Box("LAB", size/2., size/2., size/2.);

    G4Material* air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    
    G4LogicalVolume* lab_logic = new G4LogicalVolume(lab_solid, air, "LAB");
    
    lab_logic->SetVisAttributes(G4VisAttributes::Invisible);

    // Set this volume as the wrapper for the whole geometry 
    // (i.e., this is the volume that will be placed in the world)
    this->SetLogicalVolume(lab_logic);

    // First, build the moderator: a big steel cylinder.
    G4double mod_diam   = 6. * m;
    G4double mod_height = 8. * m;
    G4Tubs* moderator_solid =
      new G4Tubs("MODERATOR",  0., mod_diam/2., mod_height/2., 0, twopi);

    G4Material* steel = G4NistManager::Instance()->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    G4LogicalVolume* moderator_logic =
      new G4LogicalVolume(moderator_solid, steel, "MODERATOR");
    new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), moderator_logic,
                      "MODERATOR", lab_logic, false, 0, true);

    // Then, put air inside it
    G4double mod_thickess = 5. * cm;
    G4Tubs* inner_air_solid =
      new G4Tubs("INNER_AIR",  0., (mod_diam - 2*mod_thickess)/2., (mod_height - 2*mod_thickess)/2., 0, twopi);
    G4LogicalVolume* inner_air_logic =
      new G4LogicalVolume(inner_air_solid, air, "INNER_AIR");
    new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), inner_air_logic,
                      "INNER_AIR", moderator_logic, false, 0, true);

    // Finally, the wheel is placed inside the air.
    G4double wheel_diam = 2.5 * m;
    G4double wheel_height = 10 * cm; // to check

    G4Tubs* wheel_solid =
      new G4Tubs("WHEEL",  0., wheel_diam/2., wheel_height/2., 0, twopi);

    G4Material* wheel_mat;
    if (_ess){
      wheel_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_W");
      G4cout << "Wheel made of tungsten" << G4endl;
    } else {
      wheel_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Hg");
      G4cout << "Wheel made of mercury" << G4endl;
    }
    //G4Material* tungsten = G4NistManager::Instance()->FindOrBuildMaterial("Tungsten");
    G4LogicalVolume* wheel_logic =
      new G4LogicalVolume(wheel_solid, wheel_mat, "WHEEL");
    new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), wheel_logic,
                      "WHEEL", inner_air_logic, false, 0, true);

    
  }

  G4ThreeVector Cevns::GenerateVertex(const G4String& region) const
  {
    G4ThreeVector vertex(0., 0., 0.);
    G4double y_starting_point = 1.75 * m;
    
    if (region == "SINGLE_POINT") {
      
      vertex = G4ThreeVector(0., y_starting_point, 0.);
      
    } else if (region == "BEAM") {
      
      G4double xz_rms = 2.2 * mm;
      G4double t_rms = 3 * picosecond;
      G4double proton_speed = 2.9*1.e8 * m / s;
      G4double y_rms = t_rms * proton_speed;

      for (G4int i=0; i<3; i++) {
	if (i != 1)  { // Coordinates transverse to beam direction
	  vertex[i] = G4RandGauss::shoot(0, xz_rms);
	} else {
	  vertex[i] = G4RandGauss::shoot(y_starting_point, y_rms);
	}
      }
      
    } else {
      G4Exception("[NextNewFieldCage]", "GenerateVertex()", FatalException,
                  "Unknown vertex generation region!");
    }

    return vertex;
  }
  
}

