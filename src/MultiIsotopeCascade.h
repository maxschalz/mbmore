#ifndef MBMORE_SRC_MULTI_ISOTOPE_CASCADE_H_
#define MBMORE_SRC_MULTI_ISOTOPE_CASCADE_H_

#include <map>

#include "MultiIsotopeCentrifuge.h"
#include "MultiIsotopeHelper.h"
#include "MultiIsotopeStage.h"

namespace mbmore {

class MultiIsotopeCascade {
 public:
  MultiIsotopeCascade();
  MultiIsotopeCascade(MultiIsotopeCentrifuge centrifuge, 
                      std::map<int,double> feed_composition, 
                      double design_product_assay, 
                      double design_tails_assay, double max_cascade_feed,
                      int max_centrifuges);
  MultiIsotopeCascade(MultiIsotopeCentrifuge centrifuge, 
                      std::map<int,double> feed_composition, 
                      double design_product_assay, 
                      double design_tails_assay, double max_cascade_feed,
                      int max_centrifuges, double precision);

 private:
  // Material::Create(this, initial_feed, 
  //                  context()->GetRecipe(feed_recipe))
  MultiIsotopeCentrifuge centrifuge;
  std::map<int, MultiIsotopeStage> stage_config;
  double n_enriching;
  double n_stripping;

  std::map<int,double> feed_composition;
  std::map<int,double> product_composition;
  std::map<int,double> tails_composition;
  double design_product_assay;
  double design_tails_assay;

  double cascade_feed;
  int n_centrifuges;

  double precision;
  
  void BuildIdealCascade();
  void CalculateNStages(double &n_stages);
  double CalculateConcentrations();
  double ConcentrationDifference();

};

} // namespace mbmore

#endif  // MBMORE_SRC_MULTI_ISOTOPE_CASCADE_H_
