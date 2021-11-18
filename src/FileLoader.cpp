#include <DD4hep/DetFactoryHelper.h>
#include <DD4hep/Factories.h>
#include <DD4hep/Printout.h>
#include <XML/Utilities.h>

#include <filesystem>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

using namespace dd4hep;

void usage(int argc, char** argv)    {
  std::cout <<
    "Usage: -plugin <name> -arg [-arg]                                                  \n"
    "     name:   factory name     FileLoader                                           \n"
    "     file:<string>            file location                                        \n"
    "     url:<string>             url location                                         \n"
    "\tArguments given: " << arguments(argc,argv) << std::endl;
  std::exit(EINVAL);
}

// Plugin to download files
long load_file(
    Detector& /* desc */,
    int argc,
    char** argv
) {
  std::string file, url;
  for (int i = 0; i < argc && argv[i]; ++i) {
    if      (0 == std::strncmp("file:",argv[i], 5)) file = (argv[i] + 5);
    else if (0 == std::strncmp("url:", argv[i], 4)) url  = (argv[i] + 4);
    else usage(argc, argv);
  }
  std::cout << "Loading " << file << " from " << url << std::endl;

  if (!fs::exists(fs::path(file))) {
    std::string parent_path = fs::path(file).parent_path();
    auto ret = std::system(("mkdir -p " + parent_path + " && "
                            "curl --retry 5 -f " + url + " -o " + file).c_str());
    if (!fs::exists(fs::path(file))) {
      std::cerr << "ERROR: file, " << file << ", does not exist\n";
      std::quick_exit(1);
    }
  }

  return 0;
}

DECLARE_APPLY(FileLoader, load_file)
