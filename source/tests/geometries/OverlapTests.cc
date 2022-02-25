#include "NexusApp.h"

#include <G4UImanager.hh>

#include <catch.hpp>

using namespace nexus;

TEST_CASE("macros/NEW.init.mac test") {
  NexusApp *app = new NexusApp("macros/NEW.init.mac");
  app->Initialize();
  
  G4UImanager* UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/geometry/test/run");
  
  delete app;
}
