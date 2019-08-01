#include "Cevns.h"

#include <G4Box.hh>
#include <G4Tubs.hh>

#include <G4NistManager.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <Randomize.hh>
#include <G4VisAttributes.hh>
#include <CLHEP/Units/SystemOfUnits.h>

namespace nexus {
  
  using namespace CLHEP;

  Cevns::Cevns()
  {
  }

  Cevns::~Cevns()
  {
  }

 void Cevns::Construct()
  {
    G4double size = 4. * m;
    G4Box* air_solid = 
      new G4Box("AIR", size/2., size/2., size/2.);
    
    G4LogicalVolume* air_logic = new G4LogicalVolume(air_solid, G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"), "AIR");
    
    air_logic->SetVisAttributes(G4VisAttributes::Invisible);

    // Set this volume as the wrapper for the whole geometry 
    // (i.e., this is the volume that will be placed in the world)
    this->SetLogicalVolume(air_logic);

    G4double wheel_diam = 2.5 * m;
    G4double wheel_height = 10 * cm; // to check

    G4Tubs* wheel_solid =
      new G4Tubs("WHEEL",  0., wheel_diam/2., wheel_height/2., 0, twopi);

    G4Material* tungsten = G4NistManager::Instance()->FindOrBuildMaterial("G4_W");
    //G4Material* tungsten = G4NistManager::Instance()->FindOrBuildMaterial("Tungsten");
    G4LogicalVolume* wheel_logic =
      new G4LogicalVolume(wheel_solid, tungsten, "WHEEL");
    new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), wheel_logic,
                      "WHEEL", air_logic, false, 0, false);

    
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

