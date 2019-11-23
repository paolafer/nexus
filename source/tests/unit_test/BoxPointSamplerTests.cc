#include <BoxPointSampler.h>
#include <catch.hpp>

// Should do something more than instantiate the class
TEST_CASE("BoxPointSampler") {
  auto sampler = nexus::BoxPointSampler(1, 1, 1, 1);

  REQUIRE(1 == 1);
}
