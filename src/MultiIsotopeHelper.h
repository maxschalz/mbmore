#ifndef MBMORE_SRC_MULTI_ISOTOPE_HELPER_H_
#define MBMORE_SRC_MULTI_ISOTOPE_HELPER_H_

#include <map>

#include "cyclus.h"

namespace mbmore {

double MultiUraniumAssayAtom(Material::Ptr rsrc);

double MultiUraniumAssayMass(Material::Ptr rsrc);

// @brief Calculate stage separation factor for all isotopes starting from
//        the given U235 stage separation factor.
// @params stage separation factor for U235
// @returns a map containing the stage separation factors for all U isotopes
//          with the keys being the isotopes' mass
// @details This function calculates the stage separation factor for all 
//          uranium isotopes (U-232 to U-238, excluding U-237) starting from 
//          the U-235 stage separation factor that is passed as an argument.
//          Note that the stage separation factor is defined as the ratio of
//          abundance ratio in product to abundance ratio in tails. This 
//          method follows Houston G. Wood 'Effects of Separation Processes
//          on Minor Uranium Isotopes in Enrichment Cascades'. In: Science
//          and Global Security, 16:26--36 (2008). ISSN: 0892-9882. 
//          DOI: 10.1080/08929880802361796
std::map<int,double> CalculateSeparationFactor(double alpha_235);

} // namespace mbmore

#endif  // MBMORE_SRC_MULTI_ISOTOPE_HELPER_H_
