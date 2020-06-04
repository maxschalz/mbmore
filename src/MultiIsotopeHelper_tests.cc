
#include <gtest/gtest.h>

#include "composition.h"
#include "material.h"

#include "MultiIsotopeHelper.h"

namespace mbmore {

namespace multiisotopehelpertest {

cyclus::Composition::Ptr comp_natU() {
  std::map<int,double> m;
  m[922340000] = 5.5e-3;
  m[922350000] = 0.711;
  m[922380000] = 99.2835;
  return cyclus::Composition::CreateFromMass(m);
};
cyclus::Material::Ptr mat_natU() {
  cyclus::Composition::Ptr comp = comp_natU();
  double qty = 1.;
  return cyclus::Material::CreateUntracked(qty, comp);
};

} // namespace multiisotopehelpertest

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST(MultiIsotopeHelperTest, CheckFractionsComposition) {
  double expected_mass235 = 0.00711;
  double expected_atom235 = 0.00711 / 235.
                            / (5.5e-5/234. + 0.00711/235. + 0.992835/238.);
  cyclus::Composition::Ptr comp = multiisotopehelpertest::comp_natU();
  
  double atom_assay = MultiIsotopeAtomAssay(comp);
  EXPECT_EQ(expected_atom235, atom_assay);
  EXPECT_EQ(MultiIsotopeAtomFrac(comp, 235), expected_atom235);
  EXPECT_EQ(expected_mass235, MultiIsotopeMassAssay(comp));
  EXPECT_EQ(expected_mass235, MultiIsotopeMassFrac(comp, 235));
}


} // namespace mbmore

#ifndef CYCLUS_AGENT_TESTS_CONNECTED
int ConnectAgentTests();
static int cyclus_agent_tests_connected = ConnectAgentTests();
#define CYCLUS_AGENT_TESTS_CONNECTED cyclus_agent_tests_connected
#endif  // CYCLUS_AGENT_TESTS_CONNECTED

/*
cyclus::Agent* EnrichmentConstructor(cyclus::Context* ctx) {
  return new cycamore::Enrichment(ctx);
}

// required to get functionality in cyclus agent unit tests library
#ifndef CYCLUS_AGENT_TESTS_CONNECTED
int ConnectAgentTests();
static int cyclus_agent_tests_connected = ConnectAgentTests();
#define CYCLUS_AGENT_TESTS_CONNECTED cyclus_agent_tests_connected
#endif  // CYCLUS_AGENT_TESTS_CONNECTED

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INSTANTIATE_TEST_CASE_P(EnrichmentFac, FacilityTests,
                        Values(&EnrichmentConstructor));
INSTANTIATE_TEST_CASE_P(EnrichmentFac, AgentTests,
                        Values(&EnrichmentConstructor));
*/
