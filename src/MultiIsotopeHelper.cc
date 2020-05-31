# include "MultiIsotopeHelper.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MultiIsotopeAssayAtom(Material::Ptr rsrc) {
  return MultiIsotopeAtomFrac(rsrc, 235);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MultiIsotopeAssayMass(Material::Ptr rsrc) {
  return MultiUraniumIsotopeMass(rsrc, 235);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MultiIsotopeAtomFrac(Composition::Ptr composition, int isotope) {
  int nuc_id = (92*1000 + isotope) * 10000; // isotope to unique NucID
  std::map<int, double> v = composition->atom();
  compmath::Normalize(&v);
  return v[nuc_id];
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MultiIsotopeAtomFrac(Material::Ptr rsrc, int isotope) {
  MatQuery mq(rsrc);
  const int isotopes[6] = {232, 233, 234, 235, 236, 238};
  double isotope_assay;
  double uranium_atom_frac = 0;
  
  // Get total uranium mole fraction, all non-uranium elements are not 
  // considered here as they are directly sent to the tails.
  for (int i : isotopes) {
    int nuc_id = (92*1000 + i) * 10000;
    uranium_atom_frac += mq.atom_frac(nuc_id);
    if (i==isotope) {
      isotope_assay = mq.atom_frac(nuc_id);
    }
  }
  isotope_assay /= uranium_atom_frac;
  
  return isotope_assay ;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MultiIsotopeMassFrac(Composition::Ptr composition, int isotope) {
  int nuc_id = (92*1000 + isotope) * 10000; // isotope to unique NucID
  std::map<int, double> v = composition->mass();
  compmath::Normalize(&v);
  return v[nuc_id];
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MultiIsotopeMassFrac(Material::Ptr rsrc, int isotope) {
  MatQuery mq(rsrc);
  const int isotopes[6] = {232, 233, 234, 235, 236, 238};
  double isotope_assay;
  double uranium_mass_frac = 0;
  
  // Get total uranium mass fraction, all non-uranium elements are not 
  // considered here as they are directly sent to the tails.
  for (int i : isotopes) {
    int nuc_id = (92*1000 + i) * 10000;
    uranium_mass_frac += mq.mass_frac(nuc_id);
    if (i==isotope) {
      isotope_assay = mq.mass_frac(nuc_id);
    }
  }
  isotope_assay /= uranium_atom_frac;
  
  return isotope_assay ;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
std::map<int,double> CalculateSeparationFactor(double alpha_235) {
  std::map<int,double> separation_factors;
  
  // Convert the product to feed separation factor to overall stage 
  // separation factor.
  alpha_235 *= alpha_235
  
  // We consider U-238 to be the key component hence the mass differences
  // are calculated with respect to this isotope.
  for (int isotope = 232; isotope < 239; isotope++) {
    if (isotope != 237) {
      double delta_mass = 238. - isotope;
      double alpha = 1. + delta_mass*(alpha_235-1.) / (238.-235.);
      separation_factors[isotope] = alpha;
    }
  }
  return separation_factors;
}
