#include <BoxPointSampler.h>
#include <Randomize.hh>

#include <catch2/catch.hpp>

TEST_CASE("BoxPointSampler") {

  auto a = 2 + 20 * G4UniformRand();
  auto b = a + 2;
  auto c = a + 3;

  for (G4int i=0; i<20; i++) {
    auto thick = G4UniformRand();
    auto sampler = nexus::BoxPointSampler(a, b, c, thick);
    auto vertex  = sampler.GenerateVertex("WHOLE_VOL");

    REQUIRE(vertex.x() > -a/2 - thick);
    REQUIRE(vertex.x() <  a/2 + thick);
    REQUIRE(vertex.y() > -b/2 - thick);
    REQUIRE(vertex.y() <  b/2 + thick);
    REQUIRE(vertex.z() > -c/2 - thick);
    REQUIRE(vertex.z() <  c/2 + thick);
    
    if (((vertex.x() >  a/2.)         & (vertex.x() <  a/2. + thick)) ||
	((vertex.x() > -a/2. - thick) & (vertex.x() < -a/2.))) {
      REQUIRE(vertex.y() > -b/2 - thick);
      REQUIRE(vertex.y() <  b/2 + thick);
      REQUIRE(vertex.z() > -c/2 - thick);
      REQUIRE(vertex.z() <  c/2 + thick);
    } else if ((vertex.x() > -a/2) & (vertex.x() < a/2)) {

      
     
      G4cout << a << ", " << b << ", " << c << ", " << thick << G4endl;
      G4cout << vertex.x() << ", " << vertex.y() << ", " << vertex.z() << G4endl;
      REQUIRE((((vertex.z() >= -c/2 - thick) & (vertex.z() <= -c/2)) ||
	      ((vertex.z() >= c/2) & (vertex.z() <= c/2 + thick))));
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

}
