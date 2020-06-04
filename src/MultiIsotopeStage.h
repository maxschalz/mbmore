#ifndef MBMORE_SRC_MULTI_ISOTOPE_STAGE_H_
#define MBMORE_SRC_MULTI_ISOTOPE_STAGE_H_

// put #include <string> etc. here
#include "MultiIsotopeCentrifuge.h"
#include "MultiIsotopeHelper.h"

namespace mbmore {

class MultiIsotopeStage {
 public:
  MultiIsotopeStage();
  MultiIsotopeStage(MultiIsotopeCentrifuge centrifuge, double feed_assay,
                    double precision, double feed_flow);
  MultiIsotopeStage(double feed_assay, double feed_flow, double deltaU, 
                    double cut, double alpha, double precision);

  // Build an ideal, i.e., symetric stage where alpha equals beta.
  void BuildIdealStage();

  // Calculate the feed-to-product factor alpha using separative 
  // performance dU, cut and the centrifuge feed flow.
  void AlphaByDU();

  // Calculate the feed-to-tails factor beta using alpha, cut and feed assay.
  void BetaByAlphaAndCut();

  // Calculate the cut using alpha and beta.
  void CutByAlphaBeta();
  

 private:
  void CutIdealStage();

  double precision;
  double cut;
  double stage_feed;
  
  double deltaU;
  double alpha;
  double beta;

  double feed_assay;
  double product_assay;
  double tails_assay;
};

} // namespace mbmore

#endif  // MBMORE_SRC_MULTI_ISOTOPE_STAGE_H_
