#include <BoxPointSampler.h>

#include <catch2/catch.hpp>

TEST_CASE("BoxPointSampler") {
  auto a = 10.;
  auto b = 10.;
  auto c = 10.;
  auto thick = 1.;
  auto sampler = nexus::BoxPointSampler(a, b, c, thick);
  auto vertex  = sampler.GenerateVertex("WHOLE_VOL");
  if (((vertex.x() >  a/2.)         & (vertex.x() <  a/2. + thick)) ||
      ((vertex.x() > -a/2. - thick) & (vertex.x() < -a/2.))  ) {
    REQUIRE(vertex.y() > -b/2 - thick);
    REQUIRE(vertex.y() <  b/2 + thick);
    REQUIRE(vertex.z() > -c/2 - thick);
    REQUIRE(vertex.z() <  c/2 + thick);
  }

  if (((vertex.y() >  b/2.)         & (vertex.y() <  b/2. + thick)) ||
      ((vertex.y() > -b/2. - thick) & (vertex.y() < -b/2.))  ) {
    REQUIRE(vertex.x() > -a/2 - thick);
    REQUIRE(vertex.x() <  a/2 + thick);
    REQUIRE(vertex.z() > -c/2 - thick);
    REQUIRE(vertex.z() <  c/2 + thick);
  }

  if (((vertex.z() >  c/2.)         & (vertex.z() <  c/2. + thick)) ||
      ((vertex.z() > -c/2. - thick) & (vertex.z() < -c/2.))  ) {
    REQUIRE(vertex.x() > -a/2 - thick);
    REQUIRE(vertex.x() <  a/2 + thick);
    REQUIRE(vertex.y() > -b/2 - thick);
    REQUIRE(vertex.y() <  b/2 + thick);
  }


}
