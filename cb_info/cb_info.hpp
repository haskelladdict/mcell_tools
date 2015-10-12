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
// NOTE: MCell uses floats instead of doubles inside CellBlender format
struct Vec3 {
  float x;
  float y;
  float z;
};

// operator<< overload to print Vec3s
std::ostream& operator<<(std::ostream& s, const Vec3& v) {
  s << v.x << " " << v.y << " " << v.z;
  return s;
}

// Species keeps the position and orientation info for a single
// molecular species
struct Species {
  bool isVolMol = true;    // is the molecule a volume molecule
  std::vector<Vec3> pos;
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
  bool info = false;           // request basic file info
  bool addSeparator = false;   // add separator between species when
                               // outputing positions and orientations
  bool listMolPos = false;     // request output of molecule positions
  bool listMolOrient = false;  // request output of molecule orientation
  std::string spec;  // species to act on; empty string implies all species
  std::vector<std::string> files; // files to parse
};

#endif
