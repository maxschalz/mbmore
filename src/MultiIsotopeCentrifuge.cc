#include "MultiIsotopeCentrifuge.h"

namespace mbmore {

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
MultiIsotopeCentrifuge::MultiIsotopeCentrifuge() {
  // These are the standard values for a P1-type centrifuge
  velocity = 320;             // m s^-1
  height = 1.8;               // m
  diameter = 0.1;             // m
  temperature = 320;          // K
  
  x = 1000;
  machine_feed = 12.6e-6;     // kg s^-1
  countercurrent_to_feed = 2;

  eff = 1;
  delta_molar_mass = 0.003;   // kg mol^-1, M(U238) - M(U235)
  molar_mass = 0.352;         // kg mol^-1, M(UF6) with U238
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
MultiIsotopeCentrifuge::MultiIsotopeCentrifuge(
    double velocity, double height, double diameter, double temperature, 
    double x, double machine_feed, double countercurrent_to_feed, 
    double eff) {
  
  velocity = velocity;
  height = height;
  diameter = diameter;
  temperature = temperature;

  x = x;
  machine_feed = machine_feed;
  countercurrent_to_feed = countercurrent_to_feed;

  eff = eff;
  molar_mass = 0.352;           // kg mol^-1, molar mass of UF6 with U238
  delta_molar_mass = 0.003;     // kg mol^-1, molar mass difference between
                                // U235 and U238

}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
MultiIsotopeCentrifuge::MultiIsotopeCentrifuge(
    double velocity, double height, double diameter, double temperature, 
    double x, double machine_feed, double countercurrent_to_feed,
    double eff, double molar_mass = 0.352, double delta_molar_mass = 0.003) {
  
  velocity = velocity;
  height = height;
  diameter = diameter;
  temperature = temperature;

  x = x;
  machine_feed = machine_feed;
  countercurrent_to_feed = countercurrent_to_feed;

  eff = eff;
  molar_mass = molar_mass;
  delta_molar_mass = delta_molar_mass;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
MultiIsotopeCentrifuge MultiIsotopeCentrifuge::operator= (
    const MultiIsotopeCentrifuge &centrifuge) {
  MultiIsotopeCentrifuge cent;

  cent.delta_molar_mass = centrifuge.delta_molar_mass;
  cent.molar_mass = centrifuge.molar_mass;
  cent.x = centrifuge.x;
  cent.countercurrent_to_feed = centrifuge.countercurrent_to_feed;
  cent.velocity = centrifuge.velocity;
  cent.height = centrifuge.height;
  cent.diameter = centrifuge.diameter;
  cent.machine_feed = centrifuge.machine_feed;
  cent.temperature = centrifuge.temperature;
  cent.eff = centrifuge.eff;

  return cent;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MultiIsotopeCentrifuge::ComputeDeltaU(double cut) {
  // The naming of the variables partly follows Glaser 2009.
  const double radius = diameter / 2.;
  const double r2 = 0.975 * radius;  // Glaser 2009: between 0.96 and 0.99
  
  double r1_over_r2;
  if (velocity < 380) {
    r1_over_r2 = 0.534;
  } else {
    double factor_sqrt = 2. * gas_constant * temperature * std::log(x) 
                         / molar_mass / std::pow(velocity, 2.);
    r1_over_r2 = std::pow(1. - factor_sqrt, 0.5);
  }
  
  // The factor molar_mass_238 / molar_mass assumes that UF6 is composed
  // of U238 and 6 F19 atoms. It stems from the fact that only the uranium
  // part of UF6 is getting enriched. The U238 assumption is only valid at
  // low enrichment levels. For the same reason, the feed has to be scaled
  // with this factor.
  // TODO CALCULATE ATOMIC MASS FROM FEED ASSAY
  double machine_feed_U = machine_feed * molar_mass_238 / molar_mass;
  double factor_a = -2. * M_PI * D_rho * molar_mass_238 / molar_mass
                    / std::log(r1_over_r2);
  double A_p = factor_a * cut / machine_feed_U 
               / (1.+countercurrent_to_feed) 
               / (1.-cut+countercurrent_to_feed);
  double A_w = factor_a * (1.-cut) / machine_feed_U 
               / countercurrent_to_feed / (1.-cut+countercurrent_to_feed);
  
  // Calculate the optimum rectifier length, i.e., position of the feed 
  // point with respect to the position of the product scoop.
  double rectifier_length = (1-cut) * (1+countercurrent_to_feed) * height
                            / (1-cut+countercurrent_to_feed);
  
  double deltaU_1 = 0.5 * delta_molar_mass * std::pow(velocity,2.)
                    / gas_constant / temperature;
  double deltaU_line1 = 0.5 * machine_feed_U * cut * (1-cut) * std::pow(deltaU_1, 2.)
                        * std::pow(r2/radius, 4.)
                        * std::pow(1. - std::pow(r1_over_r2, 2.), 2.);
  double deltaU_line2 = (1.+countercurrent_to_feed) / cut
                        * (1. - std::exp(-A_p*rectifier_length));
  double deltaU_line3 = countercurrent_to_feed / (1.-cut)
                        * (1. - std::exp(-A_w*(height-rectifier_length)));
  double deltaU = deltaU_line1 * std::pow(deltaU_line2 + deltaU_line3, 2.);

  return deltaU;
}

} // namespace mbmore

