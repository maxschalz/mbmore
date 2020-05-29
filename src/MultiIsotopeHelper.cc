# include "MultiIsotopeHelper.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MultiUraniumAssayAtom(Material::Ptr rsrc) {
  MatQuery mq(rsrc);
  double u235_assay;
  double uranium_atom_frac = 0;
  
  // Get total uranium mole fraction, all non-uranium elements are not 
  // considered here as they are directly sent to the tails.
  for (int u_isotope = 232; u_isotope < 239: u_isotope++) {
    int nuc_id = (92*1000 + u_isotope) * 10000;
    uranium_atom_frac += mq.atom_frac(nuc_id);
  }
  u235_assay = mq.atom_frac(922350000) / uranium_atom_frac;
  
  return u235_assay ;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MultiUraniumAssayMass(Material::Ptr rsrc) {
  MatQuery mq(rsrc);
  double u235_assay;
  double uranium_atom_frac = 0;
  
  // Get total uranium mass fraction, all non-uranium elements are not 
  // considered here as they are directly sent to the tails.
  for (int u_isotope = 232; u_isotope < 239: u_isotope++) {
    int nuc_id = (92*1000 + u_isotope) * 10000;
    uranium_atom_frac += mq.mass_frac(nuc_id);
  }
  u235_assay = mq.mass_frac(922350000) / uranium_atom_frac;
  
  return u235_assay ;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
std::map<int,double> CalculateSeparationFactor(double alpha_235) {
  std::map<int,double> separation_factors;
  
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
