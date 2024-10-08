#pragma once

#include <chrono>
#include <cmath>
#include <fstream>
#include <index/index.h>
#include <iostream>
#include <map>
#include <nirtree/nirtree.h>
#include <nirtreedisk/nirtreedisk.h>
#include <optional>
#include <quadtree/quadtree.h>
#include <random>
#include <revisedrstartree/revisedrstartree.h>
#include <rplustree/rplustree.h>
#include <rplustreedisk/rplustreedisk.h>
#include <rstartree/rstartree.h>
#include <rstartreedisk/rstartreedisk.h>
#include <rtree/rtree.h>
#include <rtreedisk/rtreedisk.h>
#include <string>
#include <globals/globals.h>
#include <util/env.h>

enum BenchType { UNIFORM,
                 ZIPF,
                 GAUSS,
                 DATASET_FROM_FILE };
enum TreeType { R_TREE,
                R_PLUS_TREE,
                R_STAR_TREE,
                NIR_TREE,
                QUAD_TREE,
                REVISED_R_STAR_TREE };

// Tags defining how the benchmark is generated
namespace BenchTag {
struct DistributionGenerated {};
struct FileBackedReadAll {};
struct FileBackedReadChunksAtATime {};
struct Error {};
}; // namespace BenchTag

// Classes for each benchmark with their relevant constants
namespace BenchTypeClasses {
class Benchmark {};

class Uniform : public Benchmark {
public:
  static size_t size;
  static unsigned dimensions;
  static unsigned seed;
  static int precision;
};

class Zipf : public Benchmark {
public:
  static size_t size;
  static double alpha;
  static unsigned dimensions;
  static unsigned seed;
  static unsigned num_elements;

  static unsigned binary_search(double needle, std::vector<double> &cummulative);
  static double invert_prob(unsigned position);
};

class Gauss : public Benchmark {
public:
  static size_t size;
  static unsigned dimensions;
  static unsigned seed;
};
} // namespace BenchTypeClasses

// Mappings from each benchmark type to its tag
namespace BenchDetail {
template <typename T>
struct getBenchTag : BenchTag::Error {};

template <>
struct getBenchTag<BenchTypeClasses::Uniform> : BenchTag::DistributionGenerated {};

template <>
struct getBenchTag<BenchTypeClasses::Zipf> : BenchTag::DistributionGenerated {};

template <>
struct getBenchTag<BenchTypeClasses::Gauss> : BenchTag::DistributionGenerated {};
} // namespace BenchDetail

void load_dataset(std::vector<Point> &points, std::string dataset_file_name);

template <typename T>
class PointGenerator {
private:
  PointGenerator(BenchTag::DistributionGenerated);
  PointGenerator(BenchTag::FileBackedReadAll);
  PointGenerator(BenchTag::FileBackedReadChunksAtATime);

  void reset(BenchTag::DistributionGenerated);
  void reset(BenchTag::FileBackedReadAll);
  void reset(BenchTag::FileBackedReadChunksAtATime);
  std::optional<Point> nextPoint(BenchTag::DistributionGenerated);
  std::optional<Point> nextPoint(BenchTag::FileBackedReadAll);
  std::optional<Point> nextPoint(BenchTag::FileBackedReadChunksAtATime);

  // Class members
  size_t benchmarkSize;
  unsigned seed;
  std::fstream backingFile;
  size_t offset;

public:
  std::vector<Point> pointBuffer;

  static_assert(std::is_base_of<BenchTypeClasses::Benchmark, T>::value &&
                    !std::is_same<T, BenchTypeClasses::Benchmark>::value,
                "PointGenerator must take a Benchmark subclass");

  PointGenerator() : PointGenerator(BenchDetail::getBenchTag<T>{}) {
    if (T::dimensions != 0 && T::dimensions != DIM) {
      throw std::runtime_error("Wrong number of dimensions configured.");
    }
  }

  std::optional<Point> nextPoint();
  void reset();
  void generate();
};

static void fileGoodOrDie(std::fstream &file) {
  if (!file.good()) {
    std::cout << "Could not read from file: " << std::endl;
    exit(1);
  }
}

template <typename T>
PointGenerator<T>::PointGenerator(BenchTag::DistributionGenerated) : benchmarkSize(T::size), offset(0) {
}

template <typename T>
PointGenerator<T>::PointGenerator(BenchTag::FileBackedReadAll) :
  benchmarkSize(T::size),
  backingFile(util::getenv_exc(T::filePathEnv)),
  offset(0) {
}

template <typename T>
PointGenerator<T>::PointGenerator(BenchTag::FileBackedReadChunksAtATime) :
  benchmarkSize(T::size),
  backingFile(util::getenv_exc(T::filePathEnv)),
  offset(0) {
}

template <typename T>
void PointGenerator<T>::reset(BenchTag::DistributionGenerated) {
  offset = 0;
}

template <typename T>
void PointGenerator<T>::reset(BenchTag::FileBackedReadAll) {
  offset = 0;
}

template <typename T>
void PointGenerator<T>::reset(BenchTag::FileBackedReadChunksAtATime) {
  offset = 0;
  backingFile.seekg(0);
}

template <typename T>
void PointGenerator<T>::reset() {
  reset(BenchDetail::getBenchTag<T>{});
}

template <>
inline std::optional<Point> PointGenerator<BenchTypeClasses::Uniform>::nextPoint(BenchTag::DistributionGenerated) {
  // We produce all of the points at once and shove them in the buffer.
  if (pointBuffer.empty()) {
    std::cout << "Produced generator using seed: " << BenchTypeClasses::Uniform::seed << std::endl;
    std::default_random_engine generator(BenchTypeClasses::Uniform::seed);
    std::uniform_real_distribution<double> pointDist(0.0, 1.0);
    double generated_point;

    pointBuffer.reserve(benchmarkSize);
    for (size_t i = 0; i < benchmarkSize; i++) {
      Point p;
      for (unsigned d = 0; d < BenchTypeClasses::Uniform::dimensions; d++) {
        generated_point = pointDist(generator);
        if (BenchTypeClasses::Uniform::precision != -1) {
          float power_of_10 = pow(10, BenchTypeClasses::Uniform::precision);
          generated_point = std::round(generated_point * power_of_10) / power_of_10;
          assert(std::isfinite(generated_point));
        }
        p[d] = generated_point;
      }
      pointBuffer.push_back(std::move(p));
    }

    std::sort(pointBuffer.begin(), pointBuffer.end(), [](const Point &l, const Point &r) {
      if (l[0] < r[0]) {
        return true;
      } else if (l[0] > r[0]) {
        return false;
      } else {
        return l[1] < r[1];
      }
    });
    pointBuffer.erase(std::unique(pointBuffer.begin(), pointBuffer.end()), pointBuffer.end());

  } // fall through
  if (offset < pointBuffer.size()) {
    return pointBuffer[offset++];
  }
  return std::nullopt;
}

template <>
inline std::optional<Point> PointGenerator<BenchTypeClasses::Zipf>::nextPoint(BenchTag::DistributionGenerated) {
  if (pointBuffer.empty()) {
    std::vector<double> cummulative;
    cummulative.reserve(BenchTypeClasses::Zipf::num_elements);
    double running_sum = 0.0;

    // For each i of num_elements, compute the inverse prob 1/(i+1)^alpha
    // Compute the cumulative probabilities up to i and store them in cummulative
    // Sum up the cumulative probabilities into running_sum
    for (unsigned pos = 0; pos < BenchTypeClasses::Zipf::num_elements; pos++) {
      double prob = BenchTypeClasses::Zipf::invert_prob(pos);
      running_sum += prob;
      cummulative.push_back(running_sum);
    }

    // For each cumulative probability per i, divide by running_sum to get the actual probability and store
    // it in cummulative.
    for (unsigned pos = 0; pos < BenchTypeClasses::Zipf::num_elements; pos++) {
      double cum_prob = cummulative.at(pos);
      double prob = cum_prob / running_sum;
      cummulative.at(pos) = prob;
    }

    // By now, cummulative is a vector containing increasing probabilities per i.
    // Generate a uniform_probability and use it to retrieve the smallest i such that
    // probability(i) < uniform_probability

    std::cout << "Produced generator using seed: " << BenchTypeClasses::Zipf::seed << " alpha: " << BenchTypeClasses::Zipf::alpha << std::endl;
    std::default_random_engine generator(BenchTypeClasses::Zipf::seed);
    std::uniform_real_distribution<double> pointDist(0.0, 1.0);

    pointBuffer.reserve(benchmarkSize);
    for (size_t i = 0; i < benchmarkSize; i++) {
      Point p;
      p[0] = BenchTypeClasses::Zipf::binary_search(pointDist(generator), cummulative);
      for (unsigned d = 1; d < BenchTypeClasses::Zipf::dimensions; d++) {
        p[d] = pointDist(generator);
      }
      pointBuffer.push_back(std::move(p));
    }
  }

  if (offset < pointBuffer.size()) {
    return pointBuffer[offset++];
  }
  return std::nullopt;
}

template <>
inline std::optional<Point> PointGenerator<BenchTypeClasses::Gauss>::nextPoint(BenchTag::DistributionGenerated) {
  // We produce all of the points at once and shove them in the buffer.
  if (pointBuffer.empty()) {
    std::cout << "Produced generator using seed: " << BenchTypeClasses::Gauss::seed << std::endl;
    std::default_random_engine generator(BenchTypeClasses::Gauss::seed);
    std::normal_distribution<double> pointDist(0.0, 1.0);

    pointBuffer.reserve(benchmarkSize);
    for (unsigned i = 0; i < benchmarkSize; i++) {
      Point p;
      for (unsigned d = 0; d < BenchTypeClasses::Gauss::dimensions; d++) {
        p[d] = pointDist(generator);
      }
      pointBuffer.push_back(std::move(p));
    }
  } // fall through
  if (offset < pointBuffer.size()) {
    return pointBuffer[offset++];
  }
  return std::nullopt;
}

template <typename T>
std::optional<Point> PointGenerator<T>::nextPoint() {
  return nextPoint(BenchDetail::getBenchTag<T>{});
}

// Materialize all points immediately
template <typename T>
void PointGenerator<T>::generate() {
  if (!pointBuffer.empty()) {
    std::cout << "Already generated! Exiting..." << std::endl;
    abort();
  }

  this->nextPoint();
  this->reset();
}

static std::vector<Rectangle> generateRectangles(
        size_t benchmarkSize, unsigned seed,
        unsigned numRectangles, double lengthMultiplier
) {
  std::default_random_engine generator(seed + benchmarkSize);
  std::uniform_real_distribution<double> pointDist(0.0, 1.0);
  unsigned lengthSeed = 2454;
  std::default_random_engine lengthGenerator(lengthSeed);

  double minLength = 1 * lengthMultiplier;
  double maxLength = 5 * lengthMultiplier;
  std::uniform_real_distribution<double> length(minLength, maxLength);

  // Initialize rectangles
  Point ll;
  Point ur;
  std::vector<Rectangle> rectangles;
  rectangles.reserve(numRectangles);

  for (unsigned i = 0; i < numRectangles; ++i) {
    // Generate a new point and then create a square from it that covers 5% of the total area
    for (unsigned d = 0; d < dimensions; ++d) {
      ll[d] = pointDist(generator);
      ur[d] = ll[d] + length(lengthGenerator);
    }

    rectangles.push_back(Rectangle(ll, ur));
  }
  std::cout << "Initialization OK." << std::endl;

  return rectangles;
}

// Idea:
// With uniform probability, pick a generated point. It is highly likely that this point
// belongs to the axis with the most number of points. Create a bounding box with a constant
// area using that point.
static std::vector<Rectangle> generateZipfRectangles(
  std::vector<Point> points, size_t benchmarkSize, unsigned seed,
  size_t numRectangles, size_t numElements
) {
  std::default_random_engine generator(seed + benchmarkSize);
  std::uniform_int_distribution<int> pointDist(0, points.size() - 1);
  std::vector<Rectangle> rectangles;
  // For every dimension other than the first one, zipf generates values
  // between 0 and 1. Therefore, length below has to be in that range.
  const double length = 0.3;

  for (unsigned i = 0; i < numRectangles; ++i) {
    int pointIndex = pointDist(generator);
    Point ll = points[pointIndex];
    Point ur = Point(ll);

    // Create a rectangle with the appropriate lengths
    ur[0] += 1000;
    for (unsigned d = 1; d < dimensions; d++) {
      ur[d] += length;
    }

    Rectangle rectangle(ll, ur);
    rectangles.emplace_back(rectangle);
  }

  return rectangles;
}

static std::vector<Rectangle> generateRectanglesFromFile(std::string fileName) {
  std::vector<Rectangle> rectangles;

  if (fileName.empty()) return rectangles;

  std::fstream file;
  file.open(fileName);
  fileGoodOrDie(file);

  Point lowerLeft, upperRight;

  std::string line;

  while(std::getline(file, line)) {
    std::istringstream ss(line);

    for (int i = 0; i < DIM; i++) {
      ss >> lowerLeft[i];
    }

    for (int i = 0; i < DIM; i++) {
      ss >> upperRight[i];
    }

    Rectangle rect(lowerLeft, upperRight);
    rectangles.push_back(rect);
  }

  return rectangles;
}

static bool is_already_loaded(std::map<std::string, uint64_t> &configU, Index *spatial_index) {
  if (configU["tree"] == NIR_TREE) {
    auto tree = (nirtreedisk::NIRTreeDisk<NIR_MIN_FANOUT, NIR_MAX_FANOUT> *) spatial_index;
    size_t existing_page_count = tree->node_allocator_->buffer_pool_.get_preexisting_page_count();

    if (existing_page_count > 0) {
      return true;
    }
  } else if (configU["tree"] == R_PLUS_TREE) {
    auto tree = (rplustreedisk::RPlusTreeDisk<R_PLUS_MIN_FANOUT, R_PLUS_MAX_FANOUT> *) spatial_index;
    size_t existing_page_count = tree->node_allocator_->buffer_pool_.get_preexisting_page_count();

    if (existing_page_count > 0) {
      return true;
    }
  } else if (configU["tree"] == R_STAR_TREE) {
    auto tree = (rstartreedisk::RStarTreeDisk<R_STAR_MIN_FANOUT, R_STAR_MAX_FANOUT> *) spatial_index;
    size_t existing_page_count = tree->node_allocator_->buffer_pool_.get_preexisting_page_count();

    if (existing_page_count > 0) {
      return true;
    }
  } else if (configU["tree"] == R_TREE) {
    auto tree = (rtreedisk::RTreeDisk<3, 6> *) spatial_index;
    size_t existing_page_count = tree->node_allocator_->buffer_pool_.get_preexisting_page_count();

    if (existing_page_count > 0) {
      return true;
    }
  }

  return false;
}

static void
runBench(
   std::vector<Point> points,
   std::map<std::string, uint64_t> &configU,
   std::map<std::string, double> &configD,
   std::map<std::string, std::string> &configS
) {
  std::cout << "Running benchmark." << std::endl;

  // Setup statistics
  double totalTimeInserts = 0.0;
  double totalTimeSearches = 0.0;
  double totalTimeRangeSearches = 0.0;
  double totalTimeDeletes = 0.0;
  unsigned totalInserts = 0;
  unsigned totalSearches = 0;
  double totalRangeSearches = 0.0;
  unsigned totalDeletes = 0.0;
  unsigned totalPageHits = 0;
  unsigned totalPageMisses = 0;

  // Grab a reference to the buffer pool to print out stats
  buffer_pool *bufferPool;

  // Initialize the index
  Index *spatialIndex;

  if (configU["tree"] == R_TREE) {
    auto tree = new rtreedisk::RTreeDisk<R_TREE_MIN_FANOUT, R_TREE_MAX_FANOUT>(
            configU["buffer_pool_memory"], configS["db_file_name"]
    );
    bufferPool = &(tree->node_allocator_->buffer_pool_);
    spatialIndex = tree;
  } else if (configU["tree"] == R_PLUS_TREE) {
    auto tree = new rplustreedisk::RPlusTreeDisk<R_PLUS_MIN_FANOUT, R_PLUS_MAX_FANOUT>(
            configU["buffer_pool_memory"], configS["db_file_name"]
    );
    bufferPool = &(tree->node_allocator_->buffer_pool_);
    spatialIndex = tree;
  } else if (configU["tree"] == R_STAR_TREE) {
    auto tree = new rstartreedisk::RStarTreeDisk<R_STAR_MIN_FANOUT, R_STAR_MAX_FANOUT>(
            configU["buffer_pool_memory"], configS["db_file_name"]
    );
    bufferPool = &(tree->node_allocator_->buffer_pool_);
    spatialIndex = tree;
  } else if (configU["tree"] == NIR_TREE) {
    auto tree = new nirtreedisk::NIRTreeDisk<NIR_MIN_FANOUT, NIR_MAX_FANOUT>(
            configU["buffer_pool_memory"], configS["db_file_name"],
            nirtreedisk::EXPERIMENTAL_STRATEGY
    );
    bufferPool = &(tree->node_allocator_->buffer_pool_);
    spatialIndex = tree;
  } else if (configU["tree"] == QUAD_TREE) {
    spatialIndex = new quadtree::QuadTree();
  } else if (configU["tree"] == REVISED_R_STAR_TREE) {
    spatialIndex = new revisedrstartree::RevisedRStarTree(
            configU["minfanout"], configU["maxfanout"]
    );
  } else {
    std::cout << "Unknown tree selected. Exiting." << std::endl;
    return;
  }

  // Initialize search rectangles
  std::vector<Rectangle> searchRectangles;
  if (configU["distribution"] == UNIFORM) {
    searchRectangles = generateRectangles(
            configU["size"], configU["seed"],
            configU["rectanglescount"], configD["length_multiplier"]
    );
  } else if (configU["distribution"] == ZIPF) {
    searchRectangles = generateZipfRectangles(
      points, configU["size"],
      configU["seed"], configU["rectanglescount"],
      configU["num_elements"]
    );
  } else if (configU["distribution"] == DATASET_FROM_FILE) {
    searchRectangles = generateRectanglesFromFile(configS["rects_file"]);
  } else {
    // Do nothing, rectangle searches are disabled for now...
  }

  std::optional<Point> nextPoint;

  if (not is_already_loaded(configU, spatialIndex)) {
    // If we read stuff from disk and don't need to reinsert, skip this.
    // Insert points and time their insertion
    std::cout << "Tree should already be generated! Exiting..." << std::endl;
    exit(1);
  } else {
    std::cout << "Tree is loaded. Running benchmarks..." << std::endl;
  }

  // Search for points and time their retrieval
  if (configU["num_points_to_search"] > 0) {
    std::cout << "Beginning search." << std::endl;

    std::mt19937 g;
    g.seed(9812);
    std::shuffle(points.begin(), points.end(), g);

    for (Point p : points) {
      // Search
      std::chrono::high_resolution_clock::time_point begin = std::chrono::high_resolution_clock::now();
      std::vector<Point> vals = spatialIndex->search(p);

      if (vals.empty() || vals[0] != p) {
        std::cout << "could not find " << p << std::endl;
        std::cout << "Total searches: " << totalSearches << std::endl;
        exit(1);
      }

      std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> delta = std::chrono::duration_cast<std::chrono::duration<double>>(end - begin);

      totalTimeSearches += delta.count();
      totalSearches += 1;
      totalPageHits += bufferPool->page_hits;
      totalPageMisses += bufferPool->page_misses;

      if (totalSearches % 10000 == 0) {
        std::cout << "Point[" << totalSearches << "] queried. " << delta.count() << " s" << std::endl;
        bufferPool->stat();
      }

      if (totalSearches >= 100000) {
        break;
      }

      bufferPool->resetStat();
    }
    std::cout << "Search OK." << std::endl;
  } else {
    std::cout << "Test for point search is disabled" << std::endl; 
  }

	// Search for rectangles
  if (!configS["rects_file"].empty()) {
    unsigned rangeSearchChecksum = 0;

    std::cout << "Beginning search for " << searchRectangles.size() << " rectangles..." << std::endl;
    for (unsigned i = 0; i < searchRectangles.size(); ++i)
    {
      // Search
      std::cout << "-------" << std::endl;
      std::cout << "Range Search #" << totalRangeSearches << std::endl;
      std::cout << "Searching for: " << searchRectangles.at(i) << std::endl;
      std::chrono::high_resolution_clock::time_point begin = std::chrono::high_resolution_clock::now();
      std::vector<Point> v = spatialIndex->search(searchRectangles[i]);
      std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
      std::cout << "Points: " << v.size() << std::endl;
      std::chrono::duration<double> delta = std::chrono::duration_cast<std::chrono::duration<double>>(end - begin);
      std::cout << "Latency: " << delta.count() << "s" << std::endl;
      totalPageHits += bufferPool->page_hits;
      totalPageMisses += bufferPool->page_misses;
      bufferPool->stat();
      bufferPool->resetStat();

      std::cout << "-------" << std::endl;

      totalTimeRangeSearches += delta.count();
      totalRangeSearches += 1;
      rangeSearchChecksum += v.size();
    }
    std::cout << "Range search OK. Checksum = " << rangeSearchChecksum << std::endl;
    std::cout << std::endl;
    std::cout << "Total Page Hit Per Level" << std::endl;
    for (unsigned i = 0; i < spatialIndex->stats.histogramHit.size(); i++){
      if (spatialIndex->stats.histogramHit.at(i) > 0){
        std::cout << "L-" << i << ": " << spatialIndex->stats.histogramHit.at(i) << std::endl;
      }
    }
    std::cout << "Average Page Hit Per Level Per Rectangle Search" << std::endl;
    for (unsigned i = 0; i < spatialIndex->stats.histogramHit.size(); i++){
      if (spatialIndex->stats.histogramHit.at(i) > 0){
        std::cout << "L-" << i << ": " << spatialIndex->stats.histogramHit.at(i) / searchRectangles.size() << std::endl;
      }
    }
  } else {
    std::cout << "Test for range search is disabled" << std::endl; 
  }

  // Delete for points and search after deletion
  if (configU["num_points_to_delete"] > 0) {
    std::cout << "Beginning delete." << std::endl;

    std::mt19937 g;
    g.seed(7167);
    std::shuffle(points.begin(), points.end(), g);

    for (Point p : points) {
      // Remove point
      std::chrono::high_resolution_clock::time_point begin = std::chrono::high_resolution_clock::now();
      spatialIndex->remove(p);
      std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> delta = std::chrono::duration_cast<std::chrono::duration<double>>(end - begin);

      totalTimeDeletes += delta.count();
      totalDeletes += 1;
      totalPageHits += bufferPool->page_hits;
      totalPageMisses += bufferPool->page_misses;
      
      // Search removed point
      std::vector<Point> vals = spatialIndex->search(p);

      if (not vals.empty()) {
        std::cout << "Removed Point " << p << " is found" << std::endl;
        std::cout << "Output size is " << vals.size() << std::endl;
        std::cout << "Total successful delete is: " << totalDeletes << std::endl;
        exit(1);
      }

      if (totalDeletes >= configU["num_points_to_delete"]) {
        break;
      }

      bufferPool->resetStat();
    }
    std::cout << "Delete OK." << std::endl;
  } else {
    std::cout << "Test for points deletion is disabled" << std::endl; 
  }
  // Gather statistics

  // Timing Statistics
  std::cout << "Total searches: " << totalSearches << std::endl;
  std::cout << "Total time to search: " << totalTimeSearches << "s" << std::endl;
  std::cout << "Avg time to search: " << totalTimeSearches / totalSearches << "s" << std::endl;
  std::cout << "Total time to range search: " << totalTimeRangeSearches << "s" << std::endl;
  std::cout << "Avg time to range search: " << totalTimeRangeSearches / totalRangeSearches << "s" << std::endl;
  std::cout << "Total time to delete: " << totalTimeDeletes << "s" << std::endl;
  std::cout << "Avg time to delete: " << totalTimeDeletes / (double)totalDeletes << "s" << std::endl;
  std::cout << "Total page misses: " << totalPageMisses << std::endl;
  std::cout << "Total page hits + misses: " << totalPageHits + totalPageMisses << std::endl;

  // Generate visualization
  if (configU["visualization"]) {
    spatialIndex->print();
  }
}
