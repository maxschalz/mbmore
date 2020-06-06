#include "MultiIsotopeCascade.h"

#include <cmath>
#include <map>

#include "error.h"

namespace mbmore {

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
MultiIsotopeCascade::MultiIsotopeCascade() 
    : n_enriching(1),
      n_stripping(1),
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
  n_enriching = 1;
  n_stripping = 1;

  precision = precision;
  
  BuildIdealCascade();
  // ScaleCascade;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void MultiIsotopeCascade::BuildIdealCascade() {
  // Determine the number of stages in enriching and stripping section.
  CalculateNStages(n_enriching);
  CalculateNStages(n_stripping);
  
  std::map<int,MultiIsotopeStage> ideal_stages;
  
  // TODO complete initialisation of stage
  MultiIsotopeStage stage();
  double alpha, delta_U;
  
  int i = 0;
  while (i < n_enriching) {
    
    i++;
  }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void MultiIsotopeCascade::CalculateNStages(double &n_stages) {
  const double iter_max = 100;
  const double eps = 1e-4;
  double delta = 1e299;
  double previous_delta;
  
  n_stages = 0;
  do {
    n_stages++;
    previous_delta = delta;
    delta = CalculateConcentrations();
  } while (delta < previous_delta && n_stages != iter_max);
  if (n_stages == iter_max) {
    throw cyclus::Error("Unable to determine the number of stages!");
  }
  n_stages = n_stages - (1+eps);
  n_stages = previous_delta > CalculateConcentrations()
             ? n_stages + eps : n_stages + (1+eps);

  // ensure that the number of stages is an integer (stored as double)
  if (std::fmod(n_stages, 1.) < 1e-9) {
    throw cyclus::ValueError("n_stages is not a whole number!");
  }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MultiIsotopeCascade::CalculateConcentrations() {
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

