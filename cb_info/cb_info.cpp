// cb_info allows to parse and manipulate MCell's cellblender output files
//
// (C) Markus Dittrich, 2015
// Licenses under a BSD license, see LICENSE file for details

#include <cassert>
#include <exception>
#include <getopt.h>
#include <iostream>
#include <string>

#include "cb_info.hpp"

// forward declaration of static funcs
static void usage(); 
static CmdlOpts parse_cmdline(int argc, char** argv);

int main(int argc, char** argv) {

  if (argc <= 1) {
    usage();
    return EXIT_FAILURE;
  }
  auto cmdl = parse_cmdline(argc, argv);
  if (cmdl.files.size() == 0) {
    std::cerr << "***** ERROR: No MCell viz files specified to operate on\n\n";
    usage();
    return EXIT_FAILURE;
  }

  for (const auto& fileName : cmdl.files) {
    std::unordered_map<std::string, Species> specs;
    try {
      specs = parse_cb(fileName);
    } catch (std::exception& e) {
      std::cerr << "Failed to parse CellBlender file " << fileName << ": "
                << e.what() << std::endl;
      return EXIT_FAILURE;
    }

    if (cmdl.info) {
      for (const auto& s : specs) {
        std::cout << s.first << "  " << s.second.pos.size() << "\n";
      }
    }

    if (cmdl.listMolPos) {
      if (cmdl.spec != "" && specs.find(cmdl.spec) != specs.end()) {
        for (const auto& v : specs[cmdl.spec].pos) {
          std::cout << v << "\n";
        }
      }
      if (cmdl.spec == "") {
        for (const auto& s : specs) {
          std::cout << "--------- " << s.first << "\n";
          for (const auto& v : s.second.pos) {
            std::cout << v << "\n";
          }
        }
      }
    }
  }

  return EXIT_SUCCESS;
}

// parse_cb does the actual work of parsing the cellblender input file
std::unordered_map<std::string, Species> parse_cb(const std::string& fileName) {
  std::ifstream file(fileName);
  if (file.fail()) {
    throw std::runtime_error("failed to open");
  }

  std::unordered_map<std::string, Species> specs;
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
    if (type != 0) {
      throw std::runtime_error("non volume molecule encountered");
    }

    uint32_t numMols;
    read_val(file, numMols);
    spec.pos.reserve(numMols);
    assert(numMols % 3 == 0);
    for (int i = 0; i < numMols / 3; i++) {
      Vec3 v;
      read_val(file, v.x);
      read_val(file, v.y);
      read_val(file, v.z);
      spec.pos.emplace_back(std::move(v));
    }
    specs[specName] = std::move(spec);

    // check if we're done
    if (file.peek() == EOF) {
      break;
    }
  }

  return specs;
}


// usage prints a quick usage info
void usage() {
  std::cout
      << "cb_info v0.1    (C) 2015 Markus Dittrich\n\n"
      << "usage: cb_info [options] <file1> <file2> ....\n\n"
      << "Options:"
      << "\n"
      << "\t-s, --species_info       print names of species and number of\n"
      << "\t                         available molecules\n"
      << "\t-p, --print_mol_positions <species type>\n"
      << "\t                         print the (x,y,z) positions of all\n"
      << "\t                         molecules of species type\n"
      << "\t-o, --print_mol_orientations <species type>\n"
      << "\t                         print the orientations of all\n"
      << "\t                         molecules of species type\n"
      << "\t                         (empty string \"\" implies all species)\n"
      << "\t-h, --help               this help message\n" << std::endl;
}


// longOptions provides the parameters used for command line parsing
static struct option long_options[] = {
    {"species_info", no_argument, NULL, 'i'},
    {"print_mol_positions", no_argument, NULL, 'p'},
    {"print_mol_orientations", no_argument, NULL, 'o'},
    {"species_name", required_argument, NULL, 's'},
    {"help", no_argument, NULL, 'h'}};

// parse_cmdline parses the user provided command line options and
// does some basic sanity checks on them
static CmdlOpts parse_cmdline(int argc, char** argv) {

  int c;
  CmdlOpts cmdlOpts;
  while ((c = getopt_long(argc, argv, "is:polh", long_options, NULL)) != -1) {
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

      case 's':
        cmdlOpts.spec = optarg;
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

