// cb_info allows to parse and manipulate MCell's cellblender output files
//
// (C) Markus Dittrich, 2015
// Licenses under a BSD license, see LICENSE file for details

#include <cassert>
#include <iostream>
#include <limits>

#include "analyze.hpp"

static std::tuple<long long, Vec3, Vec3> compute_bounds(
    const SpecMap& specMap, const std::vector<std::string>& specs);

// analyze_mol_positions tests if molecules are uniformly distributed
void analyze_mol_positions(const SpecMap& specMap,
                           const std::vector<std::string>& specs) {
  Vec3 llc;
  Vec3 urc;
  long long numMols;
  std::tie(numMols, llc, urc) = compute_bounds(specMap, specs);
  if (numMols == 0) {
    return;
  }

  std::cout << "# mols: " << numMols << "\n"; 
  std::cout << "llc: " << llc << "\n";
  std::cout << "urc: " << urc << "\n";

  constexpr long long N = 10;
  constexpr long long N2 = N*N;
  constexpr long long N3 = N*N*N;
  Vec3 diff = urc - llc;
  Vec3 delta = diff/N;
  std::array<long long, N3> bin{};
  for (const auto& s : specs) {
    for (const auto& p : specMap.at(s).pos) {
      auto v = p - llc;
      assert(v.x >= 0);
      size_t bin_x = static_cast<int>(v.x/delta.x);
      if (bin_x < 0) {
        bin_x = 0;
      } else if (bin_x >= N) {
        bin_x = N-1;
      }

      assert(v.y >= 0);
      size_t bin_y = static_cast<int>(v.y/delta.y);
      if (bin_y < 0) {
        bin_y = 0;
      } else if (bin_y >= N) {
        bin_y = N-1;
      }

      assert(v.z >= 0);
      size_t bin_z = static_cast<int>(v.z/delta.z);
      if (bin_z < 0) {
        bin_z = 0;
      } else if (bin_z >= N) {
        bin_z = N-1;
      }

      ++bin[bin_x + bin_y * N + bin_z*N2];
    }
  }

  double expected = numMols / static_cast<double>(N3);
  double chi2 = 0;
  for (int z = 0; z < N; ++z) {
    for (int y = 0; y < N; ++y) {
      for (int x = 0; x < N; ++x) {
        double val = bin[x + y * N + z * N2] - expected;
        chi2 += val * val;
      }
    }
  }
  chi2 /= expected;

  std::cout << "chi2 computed : " << chi2 << "\n";
  std::cout << "chi2 theor    : " << chi2_ref_999 << "\n";
}


// compute_bounds computes the system bounds based on the position of all
// molecules
std::tuple<long long, Vec3, Vec3> compute_bounds(
    const SpecMap& specMap, const std::vector<std::string>& specs) {
  Vec3 llc = {std::numeric_limits<double>::max(),
              std::numeric_limits<double>::max(),
              std::numeric_limits<double>::max()};
  Vec3 urc = {-std::numeric_limits<double>::max(),
              -std::numeric_limits<double>::max(),
              -std::numeric_limits<double>::max()};

  long long numMols = 0;
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
      ++numMols;
    }
  }
  return std::make_tuple(numMols, llc, urc);
}

