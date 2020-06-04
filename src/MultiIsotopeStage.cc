#include "MultiIsotopeStage.h"

namespace mbmore {

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
MultiIsotopeStage::MultiIsotopeStage() {}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
MultiIsotopeStage::MultiIsotopeStage(MultiIsotopeCentrifuge centrifuge, 
                                     double feed_assay, 
                                     double precision = 1e-8, 
                                     double feed_flow = -1) {
  centrifuge = centrifuge;
  feed_assay = feed_assay;
  feed_flow = feed_flow;
  precision = precision;

  // TODO Finish constructor
  //BuildIdealStage();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
MultiIsotopeStage::MultiIsotopeStage(double feed_assay, double feed_flow, 
                                     double deltaU, double cut, 
                                     double alpha = -1, 
                                     double precision = 1e-8) {
  // TODO Finish constructor
}


} // namespace mbmore
