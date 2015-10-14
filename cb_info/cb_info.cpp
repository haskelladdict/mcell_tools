// cb_info allows to parse and manipulate MCell's cellblender output files
//
// (C) Markus Dittrich, 2015
// Licenses under a BSD license, see LICENSE file for details

#include <cassert>
#include <exception>
#include <getopt.h>
#include <iostream>
#include <string>

#include "analyze.hpp"
#include "cb_info.hpp"

// forward declaration of static funcs
static void usage();
static void error_and_exit(const std::string& errMsg);
static CmdlOpts parse_cmdline(int argc, char** argv);
static void read_Vec3s(std::ifstream& file, std::vector<Vec3>& vecs,
                       size_t numMols);
static std::string check_cmdline_opts(const CmdlOpts& cmdlOpts);
static std::string extract_and_check_species(const SpecMap& specMap,
                                             CmdlOpts& cmdlOpts);
static void print_positions(const SpecMap& specMap, const CmdlOpts& cmdlOpts);
static void print_orientations(const SpecMap& specMap,
                               const CmdlOpts& cmdlOpts);

// main entry point
int main(int argc, char** argv) {

  auto cmdl = parse_cmdline(argc, argv);
  std::string err = check_cmdline_opts(cmdl);
  if (err != "") {
    error_and_exit(err);
  }

  for (const auto& fileName : cmdl.files) {
    SpecMap specMap;
    try {
      specMap = parse_cb(fileName);
    } catch (std::exception& e) {
      error_and_exit("Failed to parse CellBlender file " + fileName + ": " +
                     e.what());
    }

    err = extract_and_check_species(specMap, cmdl);
    if (err != "") {
      error_and_exit(err);
    }

    if (cmdl.info) {
      for (const auto& s : specMap) {
        std::cout << s.first << "  " << s.second.pos.size()
          << "  " << (s.second.isVolMol ? "VOL" : "SURF") << "\n";
      }
    }

    if (cmdl.listMolPos) {
      print_positions(specMap, cmdl);
    }

    if (cmdl.listMolOrient) {
      for (const auto& s : cmdl.specs) {
        if (specMap.at(s).isVolMol) {
          error_and_exit("Cannot list orientations for volume mol " + s);
        }
      }
      print_orientations(specMap, cmdl);
    }

    if (cmdl.analyzePositions) {
      analyze_mol_positions(specMap, cmdl.specs);
    }
  }

  return EXIT_SUCCESS;
}

// parse_cb does the actual work of parsing the cellblender input file
SpecMap parse_cb(const std::string& fileName) {
  std::ifstream file(fileName);
  if (file.fail()) {
    throw std::runtime_error("failed to open");
  }

  SpecMap specs;
  uint32_t version;
  read_val(file, version);

  while (file) {
    unsigned char nameLength;
    read_val(file, nameLength);

    Species spec;
    char specName[nameLength + 1];
    for (int i = 0; i < nameLength; ++i) {
      read_val(file, specName[i]);
    }
    specName[nameLength] = '\0';

    unsigned char type;
    read_val(file, type);
    spec.isVolMol = (type == 0) ? true : false;

    uint32_t numMols;
    read_val(file, numMols);

    read_Vec3s(file, spec.pos, numMols);
    if (!spec.isVolMol) {
      read_Vec3s(file, spec.orient, numMols);
    }
    specs[specName] = std::move(spec);

    // check if we're done
    if (file.peek() == EOF) {
      break;
    }
  }

  return specs;
}

// read_Vec3s reads numMols worth of Vec3s from the passed ifstream
void read_Vec3s(std::ifstream& file, std::vector<Vec3>& vecs, size_t numMols) {
  vecs.reserve(numMols);
  assert(numMols % 3 == 0);
  for (size_t i = 0; i < numMols / 3; i++) {
    float x, y, z;  // NOTE: cellblender format stores floats
    read_val(file, x);
    read_val(file, y);
    read_val(file, z);
    vecs.emplace_back(std::move(Vec3{x, y, z}));
  }
}

// print_orentations prints the position info of the requested molecules
// NOTE: this function assumes that all species requested for printing
// actually exist and are surface molecule species.
void print_orientations(const SpecMap& specs, const CmdlOpts& cmdlOpts) {
  for (const auto& s : cmdlOpts.specs) {
    if (cmdlOpts.addSeparator) {
      std::cout << "--- " << s << "\n";
    }
    for (const auto& v : specs.at(s).orient) {
      std::cout << v << "\n";
    }
  }
}



// print_positions prints the position info of the requested molecules
// NOTE: this function assumes that all species requested for printing
// actually exist.
void print_positions(const SpecMap& specs, const CmdlOpts& cmdlOpts) {
  for (const auto& s : cmdlOpts.specs) {
    if (cmdlOpts.addSeparator) {
      std::cout << "--- " << s << "\n";
    }
    for (const auto& v : specs.at(s).pos) {
      std::cout << v << "\n";
    }
  }
}

// extract_and_check_species check is the requested species (for printing)
// exists. If not returns an error string. If no species were requested
// picks all available species as default. If no errors were encountered
// return an empty string.
std::string extract_and_check_species(const SpecMap& specMap,
                                      CmdlOpts& cmdlOpts) {
  for (const auto& s : cmdlOpts.specs) {
    if (specMap.find(s) == specMap.end()) {
      return "Unknown species " + s + " requested";
    }
  }

  if (cmdlOpts.specs.size() == 0) {
    for (const auto& m : specMap) {
      cmdlOpts.specs.push_back(m.first);
    }
  }

  return "";
}


// usage prints a quick usage info
void usage() {
  std::cout
      << "cb_info v0.1    (C) 2015 Markus Dittrich\n\n"
      << "usage: cb_info [options] <file1> <file2> ....\n\n"
      << "Options:"
      << "\n"
      << "\t-i, --species_info       print names of species and number of\n"
      << "\t                         available molecules\n"
      << "\t-p, --print_mol_positions <species type>\n"
      << "\t                         print the (x,y,z) positions of all\n"
      << "\t                         molecules of species type\n"
      << "\t-o, --print_mol_orientations <species type>\n"
      << "\t                         print the orientations of all\n"
      << "\t                         molecules of species type\n"
      << "\t                         (empty string \"\" implies all species)\n"
      << "\t-s, --add_separator      add separator between species in "
         "printout\n"
      << "\t-a, --analyze_positions  checks if molecules are uniformly "
         "distributed\n"
      << "\t-n, --species_name       name of species to print\n"
      << "\t-h, --help               this help message\n" << std::endl;
}


// longOptions provides the parameters used for command line parsing
static struct option long_options[] = {
    {"species_info", no_argument, NULL, 'i'},
    {"print_mol_positions", no_argument, NULL, 'p'},
    {"print_mol_orientations", no_argument, NULL, 'o'},
    {"add_separator", no_argument, NULL, 's'},
    {"species_name", required_argument, NULL, 'n'},
    {"analyze_positions", no_argument, NULL, 'a'},
    {"help", no_argument, NULL, 'h'}};

// parse_cmdline parses the user provided command line options and
// does some basic sanity checks on them
static CmdlOpts parse_cmdline(int argc, char** argv) {

  int c;
  CmdlOpts cmdlOpts;
  while ((c = getopt_long(argc, argv, "ain:spolh", long_options, NULL)) != -1) {
    switch (c) {
      case 'i':
        cmdlOpts.info = true;
        break;

      case 'p':
        cmdlOpts.listMolPos = true;
        break;

      case 'o':
        cmdlOpts.listMolOrient = true;
        break;

      case 'n':
        cmdlOpts.specs.push_back(optarg);
        break;

      case 's':
        cmdlOpts.addSeparator = true;
        break;

      case 'a':
        cmdlOpts.analyzePositions = true;
        break;

      case 'h':
      default:
        usage();
    }
  }
  while (optind < argc) {
    cmdlOpts.files.push_back(argv[optind++]);
  }

  return cmdlOpts;
}

// check_cmdline_opts does a sanity check of the provided command line
// options. If it finds a problem it returns a string with a description
// and an empty string otherwise.
std::string check_cmdline_opts(const CmdlOpts& cmdl) {
  if (cmdl.files.size() == 0) {
    return "No MCell viz files specified to operate on";
  }
  return "";
}

// error_and_exit aborts the program after printing the provided error
// message and usage info.
void error_and_exit(const std::string& errMsg) {
  std::cerr << "***** ERROR: " << errMsg << "\n\n";
  usage();
  exit(1);
}
