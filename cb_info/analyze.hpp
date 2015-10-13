// cb_info allows to parse and manipulate MCell's cellblender output files
//
// (C) Markus Dittrich, 2015
// Licenses under a BSD license, see LICENSE file for details

#ifndef ANALYZE_HPP
#define ANALYZE_HPP

#include "cb_info.hpp"


// 0.01 percentile critical values for chi-squared distribution of
// the given number of degrees of freedom (n-p) wher n is the number
// of spatial sampling boxes and p = 1 (due to the constraint that
// the number of counts is equal to the total number of atoms).
constexpr double chi2_ref_999 = 1105.917;
constexpr double chi2_ref_7999 = 8296.182;

void analyze_mol_positions(const SpecMap& specs,
                           const std::vector<std::string>& cmdl);

#endif
