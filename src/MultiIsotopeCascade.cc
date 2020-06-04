#include "MultiIsotopeCascade.h"

#include <cmath>
#include <map>

namespace mbmore {

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
MultiIsotopeCascade::MultiIsotopeCascade() 
    : n_enriching(0),
      n_stripping(0),
      n_centrifuges(0),
      cascade_feed(0),
      design_product_assay(0),
      design_tails_assay(0),
      precision(1e-8) {}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
MultiIsotopeCascade::MultiIsotopeCascade(MultiIsotopeCentrifuge centrifuge,
    std::map<int,double> feed_composition, 
    double design_product_assay, double design_tails_assay, 
    double max_cascade_feed, int max_centrifuges) 
  : MultiIsotopeCascade::MultiIsotopeCascade(centrifuge, feed_composition,
      design_product_assay, design_tails_assay, max_cascade_feed, 
      max_centrifuges, 1e-8) {}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
MultiIsotopeCascade::MultiIsotopeCascade(MultiIsotopeCentrifuge centrifuge,
    std::map<int,double> feed_composition, 
    double design_product_assay, double design_tails_assay, 
    double max_cascade_feed, int max_centrifuges, double precision) {
  
  centrifuge = centrifuge;
  
  feed_composition = feed_composition;
  design_product_assay = design_product_assay;
  design_tails_assay = design_tails_assay;

  cascade_feed = max_cascade_feed;
  n_centrifuges = max_centrifuges;
  n_enriching = 0;
  n_stripping = 0;

  precision = precision;

  // BuildCascade;
  // ScaleCascade;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void MultiIsotopeCascade::CalculateNStages() {
  // TODO include code for optimisation, e.g., L-BFGS-B
  CalculateConcentrations(n_enriching, n_stripping);
  
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MultiIsotopeCascade::CalculateConcentrations(double n_enriching,
                                                    double n_stripping) {
  // Calculates gamma = alpha*beta for all isotopes.
  // Variable naming follows E. von Halle, the equation numbers also refer
  // to his article.
  std::vector<int> isotopes;
  IsotopesNucID(isotopes);
  std::map<int,double> separation_factors;
  std::map<int,double> alpha_star;
  std::map<int,double> e;
  std::map<int,double> s;
  
  double alpha_235;
  double atom_frac;
  double delta_concentration;
  double e_sum = 0;
  double s_sum = 0;
  
  // Define the above declared variables.
  // TODO: Calculate alpha_235 of centrifuge such that alpha = beta
  alpha_235 = -1e299;
  separation_factors = CalculateSeparationFactor(alpha_235);
  for (int i : isotopes) {
    // Eq. (15)
    alpha_star[i] = separation_factors[i]
                    / std::sqrt(separation_factors[IsotopeToNucID(235)]); 
    // Eq. (37)
    e[i] = 1. / alpha_star[i] 
           / (1-std::pow(alpha_star[i],-n_enriching));
    // Eq. (39)
    s[i] = 1. / alpha_star[i] 
           / (std::pow(alpha_star[i], n_stripping+1)-1);

    atom_frac = MultiIsotopeAtomFrac(feed_composition, i);
    e_sum += e[i] * atom_frac / (e[i]+s[i]);  // Eq. (48) denominator
    s_sum += s[i] * atom_frac / (e[i]+s[i]);  // Eq. (51) denominator
  }
  
  // Calculate the compositions of product and tails.
  for (int i : isotopes) {
    atom_frac = MultiIsotopeAtomFrac(feed_composition, i);
    product_composition[i] = e[i] * atom_frac / (e[i]+s[i]) / e_sum;
    tails_composition[i] = s[i] * atom_frac / (e[i]+s[i]) / s_sum;
  }

  delta_concentration = ConcentrationDifference();
  return delta_concentration;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MultiIsotopeCascade::ConcentrationDifference() {
  int nuc_id235 = IsotopeToNucID(235);
  double product_assay = product_composition[nuc_id235];
  double tails_assay = tails_composition[nuc_id235];

  double delta_product = (product_assay-design_product_assay)
                         / design_product_assay;
  double delta_tails = (tails_assay-design_tails_assay) 
                       / design_tails_assay;
  double delta = std::sqrt(std::pow(delta_product, 2.) 
                           + std::pow(delta_tails, 2.));
  return delta;
}

} // namespace mbmore

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

/*
  using cyclus::Composition
  
*/

