# Overview

Attempt to add [catch2](https://github.com/catchorg/Catch2)-based tests to nexus.

First attempt at creating tests for petalo.

setup/teardown of the geant4 machinery is nontrivial (see test reentrancy below): this means that isolating test cases that involve G4RunManager from each other is nontrivial.

The approach I'm taking is creating a different executable target for each required G4RunManager instance (there's a target per directory); we might put all unit tests in the same target if they don't require different G4RunManagers, but more complex integration tests will require their own executable target:

```
tests/
|-- PETit_ring_test
|   |-- PETit_ring.cc
|   `-- test.cc
|-- reentrancy_test
|   |-- G4UImanagerTest.cc
|   `-- reentrancy_test.cc
`-- unit_test
    |-- BoxPointSamplerTests.cc
    |-- TrivialTest.cc
    `-- unit_tests.cc
```

# Usage

`catch2` is integrated into the scons build as another dependency.

1. ensure `catch.hpp` is available in your include path eg. in ubuntu:
  `sudo apt-get install catch`
  or just download the header to `/usr/include`
2. build:
  `scons`
3. run tests:
  `./unit_tests`, `./reentrancy_test`


## Test reentrancy

Found out that the G4RunManager destructor calls the G4UImanager destructor; this means the G4UImanager::GetUIpointer() call will return a valid value the first time, and null pointers thereafter, making execution of more than 1 test break.

Illustrated in `source/tests/reentrancy_test`, if run in gdb:
* `gdb reentrancy_test`
* `break G4UImanager::~G4UImanager()`
* `run`
* `backtrace`

we get

```
Breakpoint 1, G4UImanager::~G4UImanager (this=0x55ad7917c940, __in_chrg=<optimized out>) at /geant4/source/intercoms/src/G4UImanager.cc:119
119	G4UImanager::~G4UImanager()
(gdb) bt
#0  G4UImanager::~G4UImanager (this=0x55ad7917c940, __in_chrg=<optimized out>) at /geant4/source/intercoms/src/G4UImanager.cc:119
#1  0x00007f809278acb0 in G4RunManagerKernel::~G4RunManagerKernel (this=0x55ad7918b9f0, __in_chrg=<optimized out>)
    at /geant4/source/run/src/G4RunManagerKernel.cc:348
#2  0x00007f809278b409 in G4RunManagerKernel::~G4RunManagerKernel (this=0x55ad7918b9f0, __in_chrg=<optimized out>)
    at /geant4/source/run/src/G4RunManagerKernel.cc:356
#3  0x00007f8092779afb in G4RunManager::~G4RunManager (this=0x55ad7918b850, __in_chrg=<optimized out>)
    at /geant4/source/run/src/G4RunManager.cc:227
#4  0x00007f8092779e19 in G4RunManager::~G4RunManager (this=0x55ad7918b850, __in_chrg=<optimized out>)
    at /geant4/source/run/src/G4RunManager.cc:230
#5  0x000055ad77217d64 in run_test () at source/tests/IntegrationTest2.cc:17
#6  0x000055ad77217d9a in ::____C_A_T_C_H____T_E_S_T__ () at source/tests/IntegrationTest2.cc:22
#7  0x000055ad770e0008 in Catch::FreeFunctionTestCase::invoke (this=0x55ad79174a10) at /usr/include/catch.hpp:7343
#8  0x000055ad770c9c93 in Catch::TestCase::invoke (this=0x55ad7917a950) at /usr/include/catch.hpp:8341
#9  0x000055ad770dec30 in Catch::RunContext::invokeActiveTestCase (this=0x7ffcefb6ae70) at /usr/include/catch.hpp:6869
#10 0x000055ad770de99c in Catch::RunContext::runCurrentTest (this=0x7ffcefb6ae70, redirectedCout="", redirectedCerr="")
    at /usr/include/catch.hpp:6842
#11 0x000055ad770dd5b4 in Catch::RunContext::runTest (this=0x7ffcefb6ae70, testCase=...) at /usr/include/catch.hpp:6634
#12 0x000055ad770c6c81 in Catch::runTests (config=...) at /usr/include/catch.hpp:7014
#13 0x000055ad770df8ca in Catch::Session::runInternal (this=0x7ffcefb6b180) at /usr/include/catch.hpp:7186
#14 0x000055ad770df6b5 in Catch::Session::run (this=0x7ffcefb6b180) at /usr/include/catch.hpp:7145
#15 0x000055ad770df656 in Catch::Session::run (this=0x7ffcefb6b180, argc=1, argv=0x7ffcefb6b3e8) at /usr/include/catch.hpp:7110
#16 0x000055ad770d01a7 in main (argc=1, argv=0x7ffcefb6b3e8) at /usr/include/catch.hpp:11390
```
