#include <G4RunManager.hh>

#include <catch.hpp>

#include <G4UImanager.hh>
#include <stdio.h>

void run_test() {
  printf("TEST START\n");

  // G4UImanager is available on first run, nil on second
  G4UImanager *UI = G4UImanager::GetUIpointer();

  printf("UI %p\n", (void*) UI);

  G4RunManager *runManager = new G4RunManager;

  delete runManager;
  printf("TEST END\n");
}

TEST_CASE("integration test") {
  run_test();
  run_test();
}
