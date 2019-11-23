
// I'm instantiating the NexusApp class directly because I'm finding it hard to
// untangle the initialization of geometries and actions, we'll probably want to
// instantiate a simplified G4RunManager directly and not depend on macro files
// and convoluted configuration for test cases.

// Right now this needs to be run from a specific path to find the macro file

#include "NexusApp.h"

#include <catch.hpp>

using namespace nexus;

TEST_CASE("macros/PETit_ring.init.mac test") {
  NexusApp *app = new NexusApp("macros/PETit_ring.init.mac");
  app->Initialize();
  app->BeamOn(1);
  delete app;
}
