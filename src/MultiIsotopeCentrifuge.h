#ifndef MBMORE_SRC_MULTI_ISOTOPE_CENTRIFUGE_H_
#define MBMORE_SRC_MULTI_ISOTOPE_CENTRIFUGE_H_

#include <cmath>

namespace mbmore {

class MultiIsotopeCentrifuge {
 public:
  MultiIsotopeCentrifuge();
  MultiIsotopeCentrifuge(double velocity, double height, double diameter,
                         double machine_feed, double temperature, double eff,
                         double molar_mass, double delta_molar_mass, 
                         double x, double countercurrent_to_feed);
  
  double ComputeDeltaU(double cut);

  double delta_molar_mass;
  double molar_mass;
  double x;
  double countercurrent_to_feed;
  double velocity;
  double height;
  double diameter;
  double machine_feed;
  double temperature;
  double eff;

 private:
  const double D_rho = 2.2e-5;          // kg m^-1 s^-1
  const double gas_constant = 8.3145;   // J K^-1 mol^-1
  const double molar_mass_238 = 0.238;  // kg mol^-1

};

} // namespace mbmore

#endif  // MBMORE_SRC_MULTI_ISOTOPE_CENTRIFUGE_H_

