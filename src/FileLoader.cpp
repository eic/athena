#include <DD4hep/DetFactoryHelper.h>
#include <DD4hep/Primitives.h>
#include <DD4hep/Factories.h>
#include <DD4hep/Printout.h>
#include <XML/Utilities.h>

#include <fmt/core.h>

#include <filesystem>
#include <iostream>
#include <cstdlib>
#include <string>

namespace fs = std::filesystem;

using namespace dd4hep;

void usage(int argc, char** argv) {
  std::cout <<
    "Usage: -plugin <name> -arg [-arg]                                                  \n"
    "     name:   factory name     FileLoader                                           \n"
    "     cache:<string>           cache location (may be read-only)                    \n"
    "     file:<string>            file location                                        \n"
    "     url:<string>             url location                                         \n"
    "     cmd:<string>             download command with {0} for url, {1} for output    \n"
    "\tArguments given: " << arguments(argc,argv) << std::endl;
  std::exit(EINVAL);
}

// Plugin to download files
long load_file(
    Detector& /* desc */,
    int argc,
    char** argv
) {
  // argument parsing
  std::string cache, file, url;
  std::string cmd("curl --retry 5 -f {0} -o {1}");
  for (int i = 0; i < argc && argv[i]; ++i) {
    if      (0 == std::strncmp("cache:", argv[i], 6)) cache = (argv[i] + 6);
    else if (0 == std::strncmp("file:", argv[i], 5)) file = (argv[i] + 5);
    else if (0 == std::strncmp("url:", argv[i], 4)) url = (argv[i] + 4);
    else if (0 == std::strncmp("cmd:", argv[i], 4)) cmd = (argv[i] + 4);
    else usage(argc, argv);
  }
  printout(DEBUG, "FileLoader", "arg cache: " + cache);
  printout(DEBUG, "FileLoader", "arg file: " + file);
  printout(DEBUG, "FileLoader", "arg url: " + url);
  printout(DEBUG, "FileLoader", "arg cmd: " + cmd);

  // if file or url is empty, do nothing
  if (file.empty()) {
    printout(WARNING, "FileLoader", "no file specified");
    return 0;
  }
  if (url.empty()) {
    printout(WARNING, "FileLoader", "no url specified");
    return 0;
  }

  // parse cache for environment variables
  auto pos = std::string::npos;
  while ((pos = cache.find('$')) != std::string::npos) {
    auto after = cache.find_first_not_of(
      "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
      "abcdefghijklmnopqrstuvwxyz"
      "0123456789"
      "_",
      pos + 1);
    if (after == std::string::npos) after = cache.size(); // cache ends on env var
    auto env_name = cache.substr(pos + 1, after - pos - 1);
    auto env_value = std::getenv(env_name.c_str());
    if (env_value == nullptr) env_value = "";
    cache.erase(pos, after - pos);
    cache.insert(pos, env_value);
    printout(INFO, "FileLoader", "$" + env_name + " -> " + env_value);
  }

  // create file path
  fs::path file_path(file);

  // create hash from url, hex of unsigned long long
  std::string hash = fmt::format("{:016x}", dd4hep::detail::hash64(url)); // TODO: Use c++20 std::fmt

  // create file parent path, if not exists
  fs::path parent_path = file_path.parent_path();
  if (!fs::exists(parent_path)) {
    if (fs::create_directories(parent_path) == false) {
      printout(ERROR, "FileLoader", "parent path " + parent_path.string() + " cannot be created");
      printout(ERROR, "FileLoader", "check permissions and retry");
      std::quick_exit(1);
    }
  }

  // if file exists and is symlink to correct hash
  fs::path hash_path(parent_path / hash);
  if (fs::exists(file_path)
   && fs::equivalent(file_path, hash_path)) {
    printout(INFO, "FileLoader", "Link " + file + " -> hash " + hash + " already exists");
    return 0;
  }

  // if hash does not exist, we try to retrieve file from cache
  if (!fs::exists(hash_path)) {
    // recursive loop into cache directory
    fs::path cache_path(cache);
    printout(INFO, "FileLoader", "Cache " + cache_path.string());
    if (fs::exists(cache_path)) {
      for (auto const& dir_entry: fs::recursive_directory_iterator(cache_path)) {
        if (!dir_entry.is_directory()) continue;
        fs::path cache_dir_path = cache_path / dir_entry;
        printout(INFO, "FileLoader", "Checking " + cache_dir_path.string());
        fs::path cache_hash_path = cache_dir_path / hash;
        if (fs::exists(cache_hash_path)) {
          // symlink hash to cache/.../hash
          printout(INFO, "FileLoader", "File " + file + " with hash " + hash + " found in " + cache_hash_path.string());
          try {
            fs::create_symlink(cache_hash_path, hash_path);
          } catch (const fs::filesystem_error&) {
            printout(ERROR, "FileLoader", "unable to link from " + hash_path.string() + " to " + cache_hash_path.string());
            printout(ERROR, "FileLoader", "check permissions and retry");
            std::quick_exit(1);
          }
          break;
        }
      }
    }
  }

  // if hash does not exist, we try to retrieve file from url
  if (!fs::exists(hash_path)) {
    cmd = fmt::format(cmd, url, hash_path.c_str()); // TODO: Use c++20 std::fmt
    printout(INFO, "FileLoader", "Downloading " + file + " as hash " + hash + " with " + cmd);
    // run cmd
    auto ret = std::system(cmd.c_str());
    if (!fs::exists(hash_path)) {
      printout(ERROR, "FileLoader", "unable to run cmd " + cmd);
      printout(ERROR, "FileLoader", "check command and retry");
      std::quick_exit(1);
    }
  }

  // check if file already exists
  if (fs::exists(file_path)) {
    // file already exists
    if (fs::is_symlink(file_path)) {
      // file is symlink
      if (fs::equivalent(hash_path, fs::read_symlink(file_path))) {
        // link points to correct path
        return 0;
      } else {
        // link points to incorrect path 
        if (fs::remove(file_path) == false) {
          printout(ERROR, "FileLoader", "unable to remove symlink " + file_path.string());
          printout(ERROR, "FileLoader", "check permissions or remove manually");
          std::quick_exit(1);
        }
      }
    } else {
      // file exists but not symlink
      printout(ERROR, "FileLoader", "will not remove actual file " + file_path.string());
      printout(ERROR, "FileLoader", "check content, remove manually, and retry");
      std::quick_exit(1);
    }
  }
  // file_path now does not exist

  // symlink file_path to hash_path
  try {
    // use new path from hash so file link is local
    fs::create_symlink(fs::path(hash), file_path);
  } catch (const fs::filesystem_error&) {
    printout(ERROR, "FileLoader", "unable to link from " + file_path.string() + " to " + hash_path.string());
    printout(ERROR, "FileLoader", "check permissions and retry");
    std::quick_exit(1);
  }

  return 0;
}

DECLARE_APPLY(FileLoader, load_file)
