#ifndef MBMORE_SRC_MULTI_ISOTOPE_HELPER_H_
#define MBMORE_SRC_MULTI_ISOTOPE_HELPER_H_

#include <map>

#include "cyclus.h"

namespace mbmore {

double MultiIsotopeAssayAtom(Material::Ptr rsrc);
double MultiIsotopeAssayMass(Material::Ptr rsrc);
double MultiIsotopeAtomFrac(Composition::Ptr rsrc, int isotope);
double MultiIsotopeAtomFrac(Material::Ptr rsrc, int isotope);
double MultiIsotopeMassFrac(Composition::Ptr rsrc, int isotope);
double MultiIsotopeMassFrac(Material::Ptr rsrc, int isotope);

// Calculates the stage separation factor for all isotopes starting from 
// the given U235 product to feed separation factor.
//
// Returns a map containing the stage separation factors for all U isotopes
// with the keys being the isotopes' mass.
//
// Note that the stage separation factor is defined as the ratio of 
// abundance ratio in product to abundance ratio in tails. This method 
// follows Houston G. Wood 'Effects of Separation Processes on Minor 
// Uranium Isotopes in Enrichment Cascades'. In: Science and Global 
// Security, 16:26--36 (2008). ISSN: 0892-9882.
// DOI: 10.1080/08929880802361796
std::map<int,double> CalculateSeparationFactor(double alpha_235);

} // namespace mbmore

#endif  // MBMORE_SRC_MULTI_ISOTOPE_HELPER_H_
