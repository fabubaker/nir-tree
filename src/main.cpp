#include <bench/randomPoints.h>
#include <globals/globals.h>
#include <iostream>
#include <map>
#include <nirtree/nirtree.h>
#include <rplustree/rplustree.h>
#include <rstartree/rstartree.h>
#include <rtree/rtree.h>
#include <string>
#include <unistd.h>

void parameters(std::map<std::string, uint64_t> &configU, std::map<std::string, double> configD, std::map<std::string, std::string> &configS) {
  std::string treeTypes[] = {"R_TREE", "R_PLUS_TREE", "R_STAR_TREE", "NIR_TREE", "QUAD_TREE", "REVISED_R_STAR_TREE"};
  std::string benchTypes[] = {
      "UNIFORM", "SKEW", "CLUSTER", "CALIFORNIA", "BIOLOGICAL", "FOREST",
      "CANADA", "GAIA", "MICROSOFTBUILDINGS", "ZIPF", "GAUSS", "NYCTAXI",
      "TWEETS"};

  std::cout << "### BENCHMARK PARAMETERS ###" << std::endl;
  std::cout << "  tree = " << treeTypes[configU["tree"]] << std::endl;
  std::cout << "  benchmark = " << benchTypes[configU["distribution"]] << std::endl;
  std::cout << "  min/max branches = " << configU["minfanout"] << "/" << configU["maxfanout"] << std::endl;
  std::cout << "  n = " << configU["size"] << std::endl;
  std::cout << "  dimensions = " << dimensions << std::endl;
  std::cout << "  seed = " << configU["seed"] << std::endl;
  std::cout << "  search rectangles = " << configU["rectanglescount"] << std::endl;
  std::cout << "  visualization = " << (configU["visualization"] ? "on" : "off") << std::endl;
  std::cout << "  zipf = " << configU["alpha"] << std::endl;
  std::cout << "  buffer pool memory = " << configU["buffer_pool_memory"] << std::endl;
  std::cout << "  points per rectangle = " << configU["points_per_rectangle"] << std::endl;
  std::cout << "  points to search = " << configU["num_points_to_search"] << std::endl;
  std::cout << "  db file name = " << configS["db_file_name"] << std::endl;
  std::cout << "### ### ### ### ### ###" << std::endl << std::endl;
}

void randomPoints(std::map<std::string, uint64_t> &configU, std::map<std::string, double> &configD, std::map<std::string, std::string> &configS) {
  switch (configU["distribution"]) {
  case UNIFORM: {
    BenchTypeClasses::Uniform::size = configU["size"];
    BenchTypeClasses::Uniform::dimensions = dimensions;
    BenchTypeClasses::Uniform::seed = configU["seed"];
    if (configU.count("precision")) {
      BenchTypeClasses::Uniform::precision = configU["precision"];
      std::cout << "Using precision " << configU["precision"] << std::endl;
    }

    PointGenerator<BenchTypeClasses::Uniform> pointGen;
    runBench(pointGen, configU, configD, configS);
    break;
  }
  case SKEW: {
    PointGenerator<BenchTypeClasses::Skew> pointGen;
    runBench(pointGen, configU, configD, configS);
    break;
  }
  case CALIFORNIA: {
    PointGenerator<BenchTypeClasses::California> pointGen;
    runBench(pointGen, configU, configD, configS);
    break;
  }
  case BIOLOGICAL: {
    PointGenerator<BenchTypeClasses::Biological> pointGen;
    runBench(pointGen, configU, configD, configS);
    break;
  }
  case FOREST: {
    PointGenerator<BenchTypeClasses::Forest> pointGen;
    runBench(pointGen, configU, configD, configS);
    break;
  }
  case CANADA: {
    PointGenerator<BenchTypeClasses::Canada> pointGen;
    runBench(pointGen, configU, configD, configS);
    break;
  }
  case GAIA: {
    PointGenerator<BenchTypeClasses::Gaia> pointGen;
    runBench(pointGen, configU, configD, configS);
    break;
  }
  case MICROSOFTBUILDINGS: {
    PointGenerator<BenchTypeClasses::MicrosoftBuildings> pointGen;
    runBench(pointGen, configU, configD, configS);
    break;
  }
  case ZIPF: {
    BenchTypeClasses::Zipf::size = configU["size"];
    BenchTypeClasses::Zipf::dimensions = dimensions;
    BenchTypeClasses::Zipf::seed = configU["seed"];
    BenchTypeClasses::Zipf::num_elements = configU["num_elements"];
    BenchTypeClasses::Zipf::alpha = configD["alpha"];
    PointGenerator<BenchTypeClasses::Zipf> pointGen;
    runBench(pointGen, configU, configD, configS);
    break;
  }
  case NYCTAXI: {
    PointGenerator<BenchTypeClasses::NYCTaxi> pointGen;
    runBench(pointGen, configU, configD, configS);
    break;
  }
  case TWEETS: {
    PointGenerator<BenchTypeClasses::Tweets> pointGen;
    runBench(pointGen, configU, configD, configS);
    break;
  }
  default: {
    std::cout << "Unknown bench type." << std::endl;
    exit(1);
  }
  }
}

int main(int argc, char *argv[]) {
  // Process command line options
  int option;

  // Benchmark default configuration
  std::map<std::string, uint64_t> configU;
  configU.emplace("tree", NIR_TREE);
  configU.emplace("minfanout", 25);
  configU.emplace("maxfanout", 50);
  configU.emplace("size", 10000);
  configU.emplace("distribution", UNIFORM);
  configU.emplace("seed", 3141);
  configU.emplace("rectanglescount", 5000);
  configU.emplace("visualization", false);

  std::map<std::string, double> configD;
  std::map<std::string, std::string> configS;

  while ((option = getopt(argc, argv, "t:m:a:b:n:s:r:v:z:g:p:B:P:S:f:")) != -1) {
    switch (option) {
      case 't': // Tree
      {
        configU["tree"] = (TreeType) std::stoull(optarg);
        break;
      }
      case 'm': // Benchmark type
      {
        configU["distribution"] = (BenchType) std::stoull(optarg);
        break;
      }
      case 'a': // Minimum fanout
      {
        configU["minfanout"] = std::stoull(optarg);
        break;
      }
      case 'b': // Maximum fanout
      {
        configU["maxfanout"] = std::stoull(optarg);
        break;
      }
      case 'n': // Benchmark size
      {
        configU["size"] = std::stoull(optarg);
        break;
      }
      case 's': // Benchmark seed
      {
        configU["seed"] = std::stoull(optarg);
        break;
      }
      case 'r': // Number of search rectangles
      {
        configU["rectanglescount"] = std::stoull(optarg);
        break;
      }
      case 'v': // Visualization
      {
        configU["visualization"] = true;
        break;
      }
      case 'z': // Zipf
      {
        configD["alpha"] = std::stod(optarg);
        break;
      }
      case 'g': //num_elements
      {
        configU["num_elements"] = std::stoull(optarg);
        break;
      }
      case 'p': //precision
      {
        configU["precision"] = std::stoull(optarg);
        break;
      }
      case 'P': // Number of points per rectangle
      {
        configU["points_per_rectangle"] = std::stoull(optarg);
        break;
      }
      case 'B': // buffer pool memory in MB
      {
        configU["buffer_pool_memory"] = std::stoull(optarg) * 1000000;
        break;
      }
      case 'S': // Number of search points
      {
        configU["num_points_to_search"] = std::stoull(optarg);
        break;
      }
      case 'f': // db file name
      {
        configS["db_file_name"] = optarg;
        break;
      }

    default: {
      std::cout << "Bad option. Usage:" << std::endl;
      std::cout << "    -t  Specifies tree type {0 = R-Tree, 1 = R+-Tree, 2 = R*-Tree, 3 = NIR-Tree, 4 = Quad-Tree, 5 = RR*-Tree}" << std::endl;
      std::cout << "    -m  Specifies benchmark type {0 = Uniform, 1 = Skew, 2 = Clustered, 3 = California, 4 = Biological, 5 = Forest, 6 = Canada, 7 = Gaia, 8 = MSBuildings}" << std::endl;
      std::cout << "    -a  Minimum fanout for nodes in the selected tree" << std::endl;
      std::cout << "    -b  Maximum fanout for nodes in the selected tree" << std::endl;
      std::cout << "    -n  Specified benchmark size if size is not constant for benchmark type" << std::endl;
      std::cout << "    -s  Specifies benchmark seed if benchmark type is randomly generated" << std::endl;
      std::cout << "    -r  Specifies number of rectangles to search in benchmark if size is not constant for benchmark type" << std::endl;
      std::cout << "    -v  Turns visualization on or off for first two dimensions of the selected tree" << std::endl;
      return 1;
    }
    }
  }

  if (configS["db_file_name"].empty()) {
    std::cout << "Need to specify database file name using -f! Exiting..." << std::endl;
    exit(1);
  }

  // Print test parameters
  parameters(configU, configD, configS);

  // Run the benchmark
  randomPoints(configU, configD, configS);
}
