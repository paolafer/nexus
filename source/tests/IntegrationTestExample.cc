
// I'm instantiating the NexusApp class directly because I'm finding it hard to
// untangle the initialization of geometries and actions, we'll probably want to
// instantiate a simplified G4RunManager directly and not depend on macro files
// and convoluted configuration for test cases.

#include "NexusApp.h"

#include <G4UImanager.hh>
#include <catch.hpp>

using namespace nexus;

void run_macro(G4String init_macro) {
  NexusApp *app = new NexusApp(init_macro);
  app->Initialize();
  G4UImanager *UI = G4UImanager::GetUIpointer();
  app->BeamOn(1);
  delete app;
}

// TEST_CASE("macros/PETit_ring.init.mac test") {
//   run_macro("macros/PETit_ring.init.mac");
// }
//
// // Should be reentrant, the second run segfaults in my current setup at
// // G4UImanager::AddNewCommand(G4UIcommand*)
// TEST_CASE("macros/PETit_ring.init.mac test again") {
//   run_macro("macros/PETit_ring.init.mac");
// }
