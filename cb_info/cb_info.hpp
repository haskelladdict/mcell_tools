// cb_info allows to parse and manipulate MCell's cellblender output files
//
// (C) Markus Dittrich, 2015
// Licenses under a BSD license, see LICENSE file for details

#ifndef CB_INFO_HPP
#define CB_INFO_HPP

#include <fstream>
#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>

// Vec3 describes a 3D vector consisting of doubles
struct Vec3 {
  double x;
  double y;
  double z;
};

// several convenience operator overloads for Vec3 class
inline std::ostream& operator<<(std::ostream& s, const Vec3& v) {
  s << v.x << " " << v.y << " " << v.z;
  return s;
}

inline Vec3 operator*(const Vec3& v1, double a) {
  return Vec3{v1.x * a, v1.y * a, v1.z * a};
}

inline Vec3 operator+(const Vec3& v1, const Vec3& v2) {
  return Vec3{v1.x + v2.x, v1.y + v2.y, v1.z + v2.z};
}

inline Vec3 operator-(const Vec3& v1, const Vec3& v2) {
  return Vec3{v1.x - v2.x, v1.y - v2.y, v1.z - v2.z};
}

inline Vec3 operator/(const Vec3& v1, double a) {
  return Vec3{v1.x / a, v1.y / a, v1.z / a};
}

// Species keeps the position and orientation info for a single
// molecular species
struct Species {
  bool isVolMol = true;    // is the molecule a volume molecule
  std::vector<Vec3> pos;
  std::vector<Vec3> orient; // orientations are only used for surface species 
};

using SpecMap = std::unordered_map<std::string, Species>;

// parse_cb does the actual work of parsing the cellblender input file
SpecMap parse_cb(const std::string& fileName);

// read_val is a helper template for reading typed data from an underlying
// binary ifstream
template <typename T>
void read_val(std::ifstream& s, T& val) {
  s.read(reinterpret_cast<char*>(&val), sizeof(T));
}

// CmdlOpts describes the user selected command line options
struct CmdlOpts {
  bool info = false;              // request basic file info
  bool addSeparator = false;      // add separator between species when
                                  // outputing positions and orientations
  bool analyzePositions = false;  // check if volume molecule positions
                                  // are uniformly distributed
  bool listMolPos = false;        // request output of molecule positions
  bool listMolOrient = false;     // request output of molecule orientation
  std::vector<std::string>
      specs;  // species to act on; empty string implies all species
  std::vector<std::string> files;  // files to parse
};

#endif
