// cb_info allows to parse and manipulate MCell's cellblender output files
//
// (C) Markus Dittrich, 2015
// Licenses under a BSD license, see LICENSE file for details

#include <iostream>
#include <limits>

#include "analyze.hpp"

static std::tuple<Vec3, Vec3> compute_bounds(
    const SpecMap& specMap, const std::vector<std::string>& specs);

// analyze_mol_positions tests if molecules are uniformly distributed
void analyze_mol_positions(const SpecMap& specMap,
                           const std::vector<std::string>& specs) {
  Vec3 llc;
  Vec3 urc;
  std::tie(llc, urc) = compute_bounds(specMap, specs);

  std::cout << "llc: " << llc << "\n";
  std::cout << "urc: " << urc << "\n";
}


// compute_bounds computes the system bounds based on the position of all
// molecules
std::tuple<Vec3, Vec3> compute_bounds(const SpecMap& specMap,
                                      const std::vector<std::string>& specs) {
  Vec3 llc = {std::numeric_limits<double>::max(),
              std::numeric_limits<double>::max(),
              std::numeric_limits<double>::max()};
  Vec3 urc = {-std::numeric_limits<double>::max(),
              -std::numeric_limits<double>::max(),
              -std::numeric_limits<double>::max()};

  for (const auto& s : specs) {
    for (const auto& p : specMap.at(s).pos) {
      if (p.x < llc.x) {
        llc.x = p.x;
      }
      if (p.y < llc.y) {
        llc.y = p.y;
      }
      if (p.z < llc.z) {
        llc.z = p.z;
      }
      if (p.x > urc.x) {
        urc.x = p.x;
      }
      if (p.y > urc.y) {
        urc.y = p.y;
      }
      if (p.z > urc.z) {
        urc.z = p.z;
      }
    }
  }
  return std::make_tuple(llc, urc);
}

