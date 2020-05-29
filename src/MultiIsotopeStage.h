#ifndef MBMORE_SRC_MULTI_ISOTOPE_STAGE_H_
#define MBMORE_SRC_MULTI_ISOTOPE_STAGE_H_

// put #include <string> etc. here
#include "MultiIsotopeHelper.h"

namespace mbmore {

class MultiIsotopeCentrifuge;

class MultiIsotopeStage {
 public:
  MultiIsotopeStage() { ; }
  MultiIsotopeStage(MultiIsotopeCentrifuge centrifuge //TODO FINISH THIS)

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
  
};

} // namespace mbmore

#endif  // MBMORE_SRC_MULTI_ISOTOPE_STAGE_H_
