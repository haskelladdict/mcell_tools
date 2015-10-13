// cb_info allows to parse and manipulate MCell's cellblender output files
//
// (C) Markus Dittrich, 2015
// Licenses under a BSD license, see LICENSE file for details

#ifndef ANALYZE_HPP
#define ANALYZE_HPP

#include "cb_info.hpp"

void analyze_mol_positions(const SpecMap& specs,
                           const std::vector<std::string>& cmdl);

#endif
