#include <bulk_load.h>
#include <vector>
#include <globals/globals.h>
#include <string>

void parameters(std::map<std::string, uint64_t> &configU, std::map<std::string, double> &configD) {
  std::string treeTypes[] = {
          "R_TREE", "R_PLUS_TREE", "R_STAR_TREE",
          "NIR_TREE", "QUAD_TREE", "REVISED_R_STAR_TREE"
  };
  std::string benchTypes[] = {
          "UNIFORM", "SKEW", "CLUSTER", "CALIFORNIA", "BIOLOGICAL", "FOREST",
          "CANADA", "GAIA", "MICROSOFTBUILDINGS", "ZIPF", "GAUSS", "POIS",
          "TWEETS"};

  std::string bulkloadAlgs[] = {"STR", "QTS", "TGS"};

  std::cout << "### GEN TREE PARAMETERS ###" << std::endl;
  std::cout << "  tree = " << treeTypes[configU["tree"]] << std::endl;
  std::cout << "  benchmark = " << benchTypes[configU["distribution"]] << std::endl;
  std::cout << "  size = " << configU["size"] << std::endl;
  std::cout << "  seed = " << configU["seed"] << std::endl;
  std::cout << "  buffer pool memory = " << configU["buffer_pool_memory"] << std::endl;
  std::cout << "  bulk load percentage = " << configD["bulk_load_pct"] << std::endl;
  std::cout << "  bulk load alg = " << bulkloadAlgs[configU["bulk_load_alg"]] << std::endl;
  std::cout << "### ### ### ### ### ###" << std::endl << std::endl;
}

void generate_tree(std::map<std::string, size_t> &configU, std::map<std::string, double> &configD) {
  std::string backing_file = "bulkloaded_tree.txt";
  unlink(backing_file.c_str());

  std::vector<Point> all_points;
  std::optional<Point> next;

  if (configU["distribution"] == CALIFORNIA) {
    PointGenerator<BenchTypeClasses::California> points;
    while ((next = points.nextPoint())) {
      all_points.push_back(next.value());
    }
  } else if (configU["distribution"] == UNIFORM) {
    BenchTypeClasses::Uniform::size = configU["size"];
    BenchTypeClasses::Uniform::dimensions = dimensions;
    BenchTypeClasses::Uniform::seed = configU["seed"];
    if (configU.count("precision")) {
      BenchTypeClasses::Uniform::precision = configU["precision"];
      std::cout << "Using precision " << configU["precision"] << std::endl;
    }
    PointGenerator<BenchTypeClasses::Uniform> points;
    while ((next = points.nextPoint())) {
      all_points.push_back(next.value());
    }
  } else if (configU["distribution"] == ZIPF) {
    BenchTypeClasses::Zipf::size = configU["size"];
    BenchTypeClasses::Zipf::dimensions = dimensions;
    BenchTypeClasses::Zipf::seed = configU["seed"];
    BenchTypeClasses::Zipf::num_elements = configU["num_elements"];
    BenchTypeClasses::Zipf::alpha = configD["alpha"];

    PointGenerator<BenchTypeClasses::Zipf> points;
    while ((next = points.nextPoint())) {
      all_points.push_back(next.value());
    }
  } else if (configU["distribution"] == GAUSS) {
    BenchTypeClasses::Gauss::size = configU["size"];
    BenchTypeClasses::Gauss::dimensions = dimensions;
    BenchTypeClasses::Gauss::seed = configU["seed"];

    PointGenerator<BenchTypeClasses::Gauss> points;
    while ((next = points.nextPoint())) {
      all_points.push_back(next.value());
    }
  } else if (configU["distribution"] == POIS) {
    PointGenerator<BenchTypeClasses::Pois> points;
    while ((next = points.nextPoint())) {
      all_points.push_back(next.value());
    }
  } else if (configU["distribution"] == TWEETS) {
    PointGenerator<BenchTypeClasses::Tweets> points;
    while ((next = points.nextPoint())) {
      all_points.push_back(next.value());
    }
  }

  double bulk_load_pct = configD["bulk_load_pct"];

  uint64_t cut_off_bulk_load = std::floor(bulk_load_pct * all_points.size());
  std::cout << "Bulk loading " << cut_off_bulk_load << " points." << std::endl;
  std::cout << "Sequentially inserting " << all_points.size() - cut_off_bulk_load << " points." << std::endl;

  // Grab a reference to the buffer pool to print out stats
  buffer_pool *bufferPool;
  Index *spatialIndex;

  if (configU["tree"] == NIR_TREE) {
    nirtreedisk::NIRTreeDisk<NIR_MIN_FANOUT, NIR_MAX_FANOUT> *tree =
            new nirtreedisk::NIRTreeDisk<NIR_MIN_FANOUT, NIR_MAX_FANOUT>(
            configU["buffer_pool_memory"], backing_file, nirtreedisk::LINE_MINIMIZE_DOWN_SPLITS
    );

    spatialIndex = tree;
    bufferPool = &(tree->node_allocator_->buffer_pool_);

    // start with bulk load:
    std::cout << "Bulk Loading..." << std::endl;
    std::cout << "Creating tree with " << configU["buffer_pool_memory"] << "bytes" << std::endl;
    bulk_load_tree(
      tree, configU, all_points.begin(), all_points.begin() + cut_off_bulk_load, NIR_MAX_FANOUT,
      (nirtreedisk::LeafNode<NIR_MIN_FANOUT, NIR_MAX_FANOUT> *) nullptr,
      (nirtreedisk::BranchNode<NIR_MIN_FANOUT, NIR_MAX_FANOUT> *) nullptr
    );

    std::cout << "Buffer pool stats after bulk-loading: " << std::endl;
    bufferPool->stat();
    bufferPool->resetStat();

    // insert the rest of points:
    std::cout << "Sequential Inserting..." << std::endl;
    sequential_insert_tree(
            tree, configU,
            all_points.begin() + cut_off_bulk_load, all_points.end(),
            NIR_MAX_FANOUT
    );
    std::cout << "Created NIRTree." << std::endl;

    std::cout << "Buffer pool stats after sequential inserts: " << std::endl;
    bufferPool->stat();
    bufferPool->resetStat();
  } else if (configU["tree"] == R_STAR_TREE) {
    rstartreedisk::RStarTreeDisk<R_STAR_MIN_FANOUT, R_STAR_MAX_FANOUT> *tree = new rstartreedisk::RStarTreeDisk<R_STAR_MIN_FANOUT, R_STAR_MAX_FANOUT>(
            configU["buffer_pool_memory"], backing_file
    );
    spatialIndex = tree;
    bufferPool = &(tree->node_allocator_->buffer_pool_);

    std::cout << "Bulk Loading..." << std::endl;
    std::cout << "Creating tree with " << configU["buffer_pool_memory"] << "bytes" << std::endl;
    bulk_load_tree(
      tree, configU, all_points.begin(), all_points.begin() + cut_off_bulk_load, R_STAR_MAX_FANOUT,
      (rstartreedisk::LeafNode<R_STAR_MIN_FANOUT, R_STAR_MAX_FANOUT> *) nullptr,
      (rstartreedisk::BranchNode<R_STAR_MIN_FANOUT, R_STAR_MAX_FANOUT> *) nullptr
    );

    std::cout << "Buffer pool stats after bulk-loading: " << std::endl;
    bufferPool->stat();
    bufferPool->resetStat();

    // insert the rest of points:
    std::cout << "Sequential Inserting..." << std::endl;
    sequential_insert_tree(
            tree, configU,
            all_points.begin() + cut_off_bulk_load, all_points.end(),
            R_STAR_MAX_FANOUT
    );
    std::cout << "Created R*Tree" << std::endl;

    std::cout << "Buffer pool stats after sequential inserts: " << std::endl;
    bufferPool->stat();
    bufferPool->resetStat();
  } else if (configU["tree"] == R_PLUS_TREE) {
    rplustreedisk::RPlusTreeDisk<R_PLUS_MIN_FANOUT, R_PLUS_MAX_FANOUT> *tree = new rplustreedisk::RPlusTreeDisk<R_PLUS_MIN_FANOUT, R_PLUS_MAX_FANOUT>(
            configU["buffer_pool_memory"], backing_file
    );
    spatialIndex = tree;
    bufferPool = &(tree->node_allocator_->buffer_pool_);

    std::cout << "Bulk Loading..." << std::endl;
    std::cout << "Creating tree with " << configU["buffer_pool_memory"] << "bytes" << std::endl;
    bulk_load_tree(
            tree, configU, all_points.begin(), all_points.begin() + cut_off_bulk_load, R_STAR_MAX_FANOUT,
            (rplustreedisk::LeafNode<R_PLUS_MIN_FANOUT, R_PLUS_MAX_FANOUT> *) nullptr,
            (rplustreedisk::BranchNode<R_PLUS_MIN_FANOUT, R_PLUS_MAX_FANOUT> *) nullptr
    );

    std::cout << "Buffer pool stats after bulk-loading: " << std::endl;
    bufferPool->stat();
    bufferPool->resetStat();

    // insert the rest of points:
    std::cout << "Sequential Inserting..." << std::endl;
    sequential_insert_tree(
            tree, configU,
            all_points.begin() + cut_off_bulk_load, all_points.end(),
            R_PLUS_MAX_FANOUT
    );
    std::cout << "Created R+Tree" << std::endl;

    std::cout << "Buffer pool stats after sequential inserts: " << std::endl;
    bufferPool->stat();
    bufferPool->resetStat();
  } else {
    std::cout << "Only Supports NIR_Tree and R_STAR_TREE for gen_tree" << std::endl; 
    abort();
  }

  // Quick Test: Searching first 5000 points which are bulk loaded 
  unsigned totalSearchesLoaded = 0;
  double totalTimeSearches = 0.0;
  std::cout << "Searching for bulk-loaded points..." << std::endl;
  for (auto iter = all_points.begin(); iter < all_points.begin() + cut_off_bulk_load; iter++ ) {
    Point p = *iter; 
    std::chrono::high_resolution_clock::time_point begin = std::chrono::high_resolution_clock::now();
    std::vector<Point> out = spatialIndex->exhaustiveSearch(p);
    if (out.size() != 1) {
      int index = std::distance(all_points.begin(), iter);
      std::cout << "Could not find bulk loaded point " << p << " at index "<< index << std::endl;
      std::cout << "Output size is " << out.size() << std::endl;
      std::cout << "Total successful searches: " << totalSearchesLoaded << std::endl;
      abort();
    }
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> delta = std::chrono::duration_cast<std::chrono::duration<double>>(end - begin);
    totalTimeSearches += delta.count();
    totalSearchesLoaded += 1;
    if (totalSearchesLoaded >= 5000) {
      break;
    }
  }

  // Quick Test: Searching first 5000 points which are inserted
  std::cout << "Searching for inserted points..." << std::endl;
  unsigned totalSearchesInserted = 0;
  for (auto iter = all_points.begin() + cut_off_bulk_load; iter < all_points.end(); iter++ ) {
    Point p = *iter; 
    std::chrono::high_resolution_clock::time_point begin = std::chrono::high_resolution_clock::now();
    std::vector<Point> out = spatialIndex->search(p);
    if (out.size() != 1) {
      int index = std::distance(all_points.begin(), iter);
      std::cout << "Could not find inserted point" << p << " at index "<< index << std::endl;
      std::cout << "Output size is " << out.size() << std::endl;
      std::cout << "Total successful searches: " << totalSearchesInserted << std::endl;
      abort();
    }
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> delta = std::chrono::duration_cast<std::chrono::duration<double>>(end - begin);
    totalTimeSearches += delta.count();
    totalSearchesInserted += 1;
    if (totalSearchesInserted >= 5000) {
      break;
    }
  }
  
  spatialIndex->stat();

  std::cout << "Total time to search: " << totalTimeSearches << "s" << std::endl;
  std::cout << "Avg time to search: " << totalTimeSearches / (totalSearchesInserted + totalSearchesLoaded) << "s" << std::endl;

  return;
}

int main(int argc, char **argv) {
  int option;
  std::map<std::string, uint64_t> configU;
  std::map<std::string, double> configD;

  configU.emplace("tree", NIR_TREE);
  configU.emplace("distribution", CALIFORNIA);
  configU.emplace("seed", 0);
  configD.emplace("bulk_load_pct", 1.0);

  while ((option = getopt(argc, argv, "t:m:n:s:p:g:z:B:A:b:")) != -1) {
    switch (option) {
    case 't': {
      configU["tree"] = (TreeType)std::stoull(optarg);
      break;
    }
    case 'm': {
      configU["distribution"] = (BenchType)std::stoull(optarg);
      break;
    }
    case 'n': {
      configU["size"] = std::stoull(optarg);
      break;
    }
    case 'g': {
      configU["num_elements"] = std::stoull(optarg);
      break;
    }
    case 'z': {
      configD["alpha"] = std::stod(optarg);
      break;
    }
    case 's': {
      configU["seed"] = std::stoull(optarg);
      break;
    }
    case 'p': {
      configU["precision"] = std::stoull(optarg);
      break;
    }
    case 'B': // buffer pool memory in MB
    {
      configU["buffer_pool_memory"] = std::stoull(optarg) * 1000000;
      break;
    }
    case 'A':
    {
      configU["bulk_load_alg"] = std::stoull(optarg);
      break;
    }
    case 'b': {
      configD["bulk_load_pct"] = std::stod(optarg);
      break;
    }
    }
  }

  // Print gen_tree parameters
  parameters(configU, configD);

  // Generate the tree
  generate_tree(configU, configD);

  return 0;
}
