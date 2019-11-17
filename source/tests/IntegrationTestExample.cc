// First attempt at creating tests for petalo, starting from integration tests
// as these are easiest to understand conceptually.
//
// I'm instantiating the NexusApp class directly because I'm finding it hard to
// untangle the initialization of geometries and actions, we'll probably want to
// instantiate a simplified G4RunManager directly and not depend on macro files
// and convoluted configuration for test cases.
//
// My initial goal is to successfully launch a few macros from the tests,
// this will at least validate the macro file syntax and tell us if we have made
// a change that makes them segfault.
//
// This might be harder than it looks for several reasons:
// * destructors not cleaning up after themselves properly
// * my lack of understanding of Geant4 setup/teardown semantics

#include "NexusApp.h"
#include <catch.hpp>

using namespace nexus;

void run_macro(G4String init_macro) {
  NexusApp *app = new NexusApp(init_macro);
  app->Initialize();
  app->BeamOn(1);
  delete app;
}

TEST_CASE("macros/PETit_ring.init.mac test") {
  run_macro("macros/PETit_ring.init.mac");
}

// Should be reentrant, the second run segfaults in my current setup at
// G4UImanager::AddNewCommand(G4UIcommand*)
// probably because I've failed to install some UI deps or failed to set flags
// to disable compiling these
TEST_CASE("macros/PETit_ring.init.mac test again") {
  run_macro("macros/PETit_ring.init.mac");
}
