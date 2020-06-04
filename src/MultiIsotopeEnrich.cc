#include "MultiIsotopeEnrich.h"

namespace mbmore {

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
MultiIsotopeEnrich::MultiIsotopeEnrich(cyclus::Context* ctx)
    : cyclus::Facility(ctx),
      feed_recipe(""),
      initial_feed(),
      design_feed_flow(),
      max_centrifuges(),
      design_feed_assay(),
      design_product_assay(),
      design_tails_assay(),
      temperature(320.),
      velocity(485.),
      height(0.5),
      diameter(0.15),
      countercurrent_to_feed(2.),
      machine_feed(15.),
      x(1000.),
      eff(1.),
      feed_commod(""),
      product_commod(""),
      tails_commod(""),
      order_prefs(true) {}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
MultiIsotopeEnrich::~MultiIsotopeEnrich() {}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
std::string MultiIsotopeEnrich::str() {
  std::stringstream ss;
  ss << cyclus::Facility::str() << " with enrichment facility parameters:"
     << " * Tails assay: " << design_tails_assay 
     << " * Feed assay: " << design_feed_assay
     << " * Input cyclus::Commodity: " << feed_commod
     << " * Output cyclus::Commodity: " << product_commod
     << " * Tails cyclus::Commodity: " << tails_commod;
  return ss.str();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void MultiIsotopeEnrich::EnterNotify() {
  cyclus::Facility::EnterNotify();

  // TODO make feed_composition a class variable?
  cyclus::Composition* feed_composition = context->GetRecipe(feed_recipe);
  centrifuge = MultiIsotopeCentrifuge(velocity, height, diameter, 
                                      temperature, x, machine_feed,
                                      countercurrent_to_feed, eff);
  // create cascade object

}

} // namespace mbmore
