// cb_info allows to parse and manipulate MCell's cellblender output files
//
// (C) Markus Dittrich, 2015
// Licenses under a BSD license, see LICENSE file for details

#include <cassert>
#include <iostream>
#include <limits>

#include "analyze.hpp"

using binArray = std::array<long long, N3>;

// static helper functions
static std::tuple<Vec3, Vec3> compute_bounds(
    const SpecMap& specMap, const std::vector<std::string>& specs);
static std::tuple<binArray, long long> compute_bins(
    const SpecMap& specMap, const std::vector<std::string>& specs,
    const Vec3& llc, const Vec3& urc);
static void print_results(const Vec3& v1, const Vec3& v2, double chi);

// analyze_mol_positions tests if molecules are uniformly distributed
std::string analyze_mol_positions(const SpecMap& specMap,
                                  const std::vector<std::string>& specs) {
  Vec3 llc;
  Vec3 urc;
  std::tie(llc, urc) = compute_bounds(specMap, specs);

  binArray bin;
  long long numMols;
  std::tie(bin, numMols) = compute_bins(specMap, specs, llc, urc);
  if (numMols == 0) {
    return "selection contains no molecules";
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
  print_results(llc, urc, chi2);
  return "";
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

// print_results outputs the results of the chi squared analysis
void print_results(const Vec3& llc, const Vec3& urc, double chi2) {
  std::cout << "\n------ system dimensions ------------------\n"
            << "LLC: " << llc << "\n"
            << "URC: " << urc << "\n\n";

  if (chi2 < chi2_ref_999) {
    std::cout << "selected molecules are uniformly distributed (p = 0.01)\n";
  } else {
    std::cout
        << "selected molecules are *not* uniformly distributed (p = 0.01)\n";
  }
  std::cout << "CHI^2: " << chi2 << "/" << chi2_ref_999
            << " (computed/cutoff)\n";
}

// compute_bins computes the 3D binning of all selected molecules
std::tuple<binArray, long long> compute_bins(
    const SpecMap& specMap, const std::vector<std::string>& specs,
    const Vec3& llc, const Vec3& urc) {

  long long numMols = 0;
  Vec3 delta = (urc - llc) / N;
  binArray bin{};
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
      ++numMols;
    }
  }
  return std::make_tuple(bin, numMols);
}

