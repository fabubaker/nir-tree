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

static const unsigned BitDataSize = 60000;
static const unsigned BitQuerySize = 3164;
static const unsigned HazeDataSize = 65365;
static const unsigned HazeQuerySize = 3244;
static const unsigned CaliforniaDataSize = 1875347; //1888012;
static const unsigned CaliforniaQuerySize = 5974;
static const unsigned BiologicalDataSize = 11958999;
static const unsigned BiologicalQuerySize = 37844;
static const unsigned ForestDataSize = 581012;
static const unsigned ForestQuerySize = 1838;
static const unsigned CanadaDataSize = 19371405;
static const unsigned CanadaQuerySize = 5000;
static const unsigned GaiaDataSize = 18084053;
static const unsigned GaiaQuerySize = 5000;
static const unsigned MicrosoftBuildingsDataSize = 752704741;
static const unsigned NYCTaxiDataSize = 81616580;

enum BenchType { UNIFORM,
                 SKEW,
                 CLUSTER,
                 CALIFORNIA,
                 BIOLOGICAL,
                 FOREST,
                 CANADA,
                 GAIA,
                 MICROSOFTBUILDINGS,
                 ZIPF,
                 GAUSS,
                 NYCTAXI };
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
  static constexpr char fileName[] = "";
  static int precision;
};

class Skew : public Benchmark {
public:
  static constexpr size_t size = BitDataSize;
  static constexpr unsigned querySize = BitQuerySize;
  static constexpr unsigned dimensions = 2;
  static constexpr char fileName[] = "/home/bjglasbe/Documents/code/nir-tree/data/bits02";
};

class Zipf : public Benchmark {
public:
  static size_t size;
  static double alpha;
  static unsigned dimensions;
  static unsigned seed;
  static constexpr char fileName[] = "";
  static unsigned num_elements;

  static unsigned binary_search(double needle, std::vector<double> &cummulative);
  static double invert_prob(unsigned position);
};

class Gauss : public Benchmark {
public:
  static size_t size;
  static unsigned dimensions;
  static unsigned seed;
  static constexpr char fileName[] = "";
};

class California : public Benchmark {
public:
  static constexpr size_t size = CaliforniaDataSize;
  static constexpr unsigned querySize = CaliforniaQuerySize;
  static constexpr unsigned dimensions = 2;
  static constexpr char fileName[] =
      "/home/bjglasbe/Documents/code/nir-tree/data/uniq_cali";
  //california";
};

class Biological : public Benchmark {
public:
  static constexpr size_t size = BiologicalDataSize;
  static constexpr unsigned querySize = BiologicalQuerySize;
  static constexpr unsigned dimensions = 3;
  static constexpr char fileName[] = "/home/bjglasbe/Documents/code/nir-tree/data/biological";
};

class Forest : public Benchmark {
public:
  static constexpr size_t size = ForestDataSize;
  static constexpr unsigned querySize = ForestQuerySize;
  static constexpr unsigned dimensions = 5;
  static constexpr char fileName[] = "/home/bjglasbe/Documents/code/nir-tree/data/forest";
};

class Canada : public Benchmark {
public:
  static constexpr size_t size = CanadaDataSize;
  static constexpr unsigned querySize = CanadaQuerySize;
  static constexpr unsigned dimensions = 2;
  static constexpr char fileName[] = "/home/bjglasbe/Documents/code/nir-tree/data/canada";
};

class Gaia : public Benchmark {
public:
  static constexpr size_t size = GaiaDataSize;
  static constexpr unsigned querySize = GaiaQuerySize;
  static constexpr unsigned dimensions = 3;
  static constexpr char fileName[] = "/home/bjglasbe/Documents/code/nir-tree/data/gaia";
};

class MicrosoftBuildings : public Benchmark {
public:
  static constexpr size_t size = MicrosoftBuildingsDataSize;
  static constexpr unsigned querySize = 0;
  static constexpr unsigned dimensions = 2;
  static constexpr char fileName[] = "/home/bjglasbe/Documents/code/nir-tree/data/microsoftbuildings";
};
class NYCTaxi : public Benchmark {
public:
  static constexpr size_t size = NYCTaxiDataSize;
  static constexpr unsigned querySize = 0;
  static constexpr unsigned dimensions = 2;
  static constexpr char fileName[] = "/hdd1/nir-tree++/data/NYCTaxi.csv";
};
}; // namespace BenchTypeClasses

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

template <>
struct getBenchTag<BenchTypeClasses::Skew> : BenchTag::FileBackedReadAll {};

template <>
struct getBenchTag<BenchTypeClasses::California> : BenchTag::FileBackedReadAll {};

template <>
struct getBenchTag<BenchTypeClasses::Biological> : BenchTag::FileBackedReadAll {};

template <>
struct getBenchTag<BenchTypeClasses::Forest> : BenchTag::FileBackedReadAll {};

template <>
struct getBenchTag<BenchTypeClasses::Canada> : BenchTag::FileBackedReadAll {};

template <>
struct getBenchTag<BenchTypeClasses::Gaia> : BenchTag::FileBackedReadAll {};

template <>
struct getBenchTag<BenchTypeClasses::MicrosoftBuildings> : BenchTag::FileBackedReadChunksAtATime {};

template <>
struct getBenchTag<BenchTypeClasses::NYCTaxi> : BenchTag::FileBackedReadAll {};
} // namespace BenchDetail

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
PointGenerator<T>::PointGenerator(BenchTag::FileBackedReadAll) : benchmarkSize(T::size), backingFile(T::fileName), offset(0) {
}

template <typename T>
PointGenerator<T>::PointGenerator(BenchTag::FileBackedReadChunksAtATime) : benchmarkSize(T::size), backingFile(T::fileName), offset(0) {
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
std::optional<Point> PointGenerator<T>::nextPoint(BenchTag::FileBackedReadAll) {
  if (pointBuffer.empty()) {
    // We produce all of the points at once, shove them into the buffer.
    fileGoodOrDie(backingFile);

    // Initialize points
    pointBuffer.reserve(benchmarkSize);
    for (size_t i = 0; i < benchmarkSize; ++i) {
      Point p;
      for (unsigned d = 0; d < T::dimensions; ++d) {
        fileGoodOrDie(backingFile);
        double dbl;
        backingFile >> dbl;
        p[d] = dbl;
      }
      pointBuffer.push_back(std::move(p));
    }
  }

  if (offset < pointBuffer.size()) {
    return pointBuffer[offset++];
  }
  return std::nullopt;
}

template <typename T>
std::optional<Point> PointGenerator<T>::nextPoint(BenchTag::FileBackedReadChunksAtATime) {
  if (offset >= benchmarkSize) {
    return std::nullopt;
  }
  if (pointBuffer.empty()) {
    // Fill it with 10k entries
    pointBuffer.resize(10000);
  }
  if (offset % 10000 == 0) {
    // Time to read 10k more things
    // Do the fstream song and dance

    // Initialize points
    for (size_t i = 0; i < 10000 and offset + i < benchmarkSize; ++i) {
      Point p;
      for (unsigned d = 0; d < T::dimensions; ++d) {
        fileGoodOrDie(backingFile);
        double dbl;
        backingFile >> dbl;
        p[d] = dbl;
      }
      pointBuffer[i] = p;
    }
  }

  return pointBuffer[offset++ % 10000];
}

template <typename T>
std::optional<Point> PointGenerator<T>::nextPoint() {
  return nextPoint(BenchDetail::getBenchTag<T>{});
}

static std::vector<Rectangle> generateRectangles(size_t benchmarkSize, unsigned seed, unsigned rectanglesSize, size_t pointsPerRectangle) {
  std::default_random_engine generator(seed + benchmarkSize);
  std::uniform_real_distribution<double> pointDist(0.0, 1.0);

  // Initialize rectangles
  Point ll;
  Point ur;
  std::vector<Rectangle> rectangles;
  rectangles.reserve(rectanglesSize);
  // Compute the dimensions-th root of a percentage that will give rectangles that in expectation return "pointsPerRectangle" points
  double requiredPercentage = pointsPerRectangle / (double)benchmarkSize;
  double root = std::pow(requiredPercentage, 1.0 / (double)dimensions);
  std::cout << "Beginning initialization of " << rectanglesSize << " rectangles with " << requiredPercentage << "% and " << root << "..." << std::endl;
  for (unsigned i = 0; i < rectanglesSize; ++i) {
    Rectangle rect;
    // Generate a new point and then create a square from it that covers 5% of the total area
    for (unsigned d = 0; d < dimensions; ++d) {
      ll[d] = pointDist(generator);
      ur[d] = ll[d] + root;
    }

    rectangles.push_back(Rectangle(ll, ur));
  }
  std::cout << "Initialization OK." << std::endl;

  return rectangles;
}

static std::vector<Rectangle> generateBitRectangles() {
  // Query set is pre-generated and requires 2 or 3 dimensions
  assert(dimensions == 2 || dimensions == 3);

  // Setup file reader and double buffer
  std::fstream file;
  std::string dataPath = dimensions == 2 ? "/home/bjglasbe/Documents/code/nir-tree/data/bit02.2" : "/home/bjglasbe/Documents/code/nir-tree/data/bit03.2";
  file.open(dataPath.c_str());
  fileGoodOrDie(file);
  char *buffer = new char[sizeof(double)];
  memset(buffer, 0, sizeof(double));
  double *doubleBuffer = (double *)buffer;

  // Initialize rectangles
  std::vector<Rectangle> rectangles;
  rectangles.reserve(BitQuerySize);
  std::cout << "Beginning initialization of " << BitQuerySize << " computer rectangles..." << std::endl;
  for (unsigned i = 0; i < BitQuerySize; ++i) {
    Rectangle rect;
    for (unsigned d = 0; d < dimensions; ++d) {
      fileGoodOrDie(file);
      file.read(buffer, sizeof(double));
      rect.lowerLeft[d] = *doubleBuffer;
      fileGoodOrDie(file);
      file.read(buffer, sizeof(double));
      rect.upperRight[d] = *doubleBuffer;
    }
    rectangles.push_back(rect);
  }
  std::cout << "Initialization OK." << std::endl;

  // Cleanup
  file.close();
  delete[] buffer;

  return rectangles;
}

static std::vector<Rectangle> generateHazeRectangles() {
  // Query set is pre-generated and requires 2 or 3 dimensions
  assert(dimensions == 2 || dimensions == 3);

  // Setup file reader and double buffer
  std::fstream file;
  std::string dataPath = dimensions == 2 ? "/home/bjglasbe/Documents/code/nir-tree/data/pha02.2" : "/home/bjglasbe/Documents/code/nir-tree/data/pha03.2";
  file.open(dataPath.c_str());
  fileGoodOrDie(file);
  char *buffer = new char[sizeof(double)];
  memset(buffer, 0, sizeof(double));
  double *doubleBuffer = (double *)buffer;

  // Initialize rectangles
  std::vector<Rectangle> rectangles;
  rectangles.reserve(HazeQuerySize);
  std::cout << "Beginning initialization of " << HazeQuerySize << " hazy rectangles..." << std::endl;
  for (unsigned i = 0; i < HazeQuerySize; ++i) {
    Rectangle rect;
    for (unsigned d = 0; d < dimensions; ++d) {
      fileGoodOrDie(file);
      file.read(buffer, sizeof(double));
      rect.lowerLeft[d] = *doubleBuffer;
      fileGoodOrDie(file);
      file.read(buffer, sizeof(double));
      rect.upperRight[d] = *doubleBuffer;
    }
    rectangles.push_back(rect);
  }
  std::cout << "Initialization OK." << std::endl;

  // Cleanup
  file.close();
  delete[] buffer;

  return rectangles;
}

static std::vector<Rectangle> generateCaliRectangles() {
  // Query set is pre-generated and requires 2 dimensions
  assert(dimensions == 2);

  // Setup file reader and double buffer
  std::fstream file;
  std::string dataPath = "/home/bjglasbe/Documents/code/nir-tree/data/rea02.2";
  file.open(dataPath);
  fileGoodOrDie(file);
  char *buffer = new char[sizeof(double)];
  memset(buffer, 0, sizeof(double));
  double *doubleBuffer = (double *)buffer;

  // Initialize rectangles
  std::vector<Rectangle> rectangles;
  rectangles.reserve(CaliforniaQuerySize);
  std::cout << "Beginning initialization of " << CaliforniaQuerySize << " california roll rectangles..." << std::endl;
  for (unsigned i = 0; i < CaliforniaQuerySize; ++i) {
    Rectangle loc_rect;

    for (unsigned d = 0; d < dimensions; ++d) {
      fileGoodOrDie(file);
      file.read(buffer, sizeof(double));
      loc_rect.lowerLeft[d] = *doubleBuffer;
      fileGoodOrDie(file);
      file.read(buffer, sizeof(double));
      loc_rect.upperRight[d] = *doubleBuffer;
    }
    rectangles.push_back(loc_rect);
  }
  std::cout << "Initialization OK." << std::endl;

  // Cleanup
  file.close();
  delete[] buffer;

  return rectangles;
}

static std::vector<Rectangle> generateBioRectangles() {
  // Query set is pre-generated and requires 3 dimensions
  assert(dimensions == 3);

  // Setup file reader and double buffer
  std::fstream file;
  std::string dataPath = "/home/bjglasbe/Documents/code/nir-tree/data/rea03.2";
  file.open(dataPath);
  fileGoodOrDie(file);
  char *buffer = new char[sizeof(double)];
  memset(buffer, 0, sizeof(double));
  double *doubleBuffer = (double *)buffer;

  // Initialize rectangles
  std::vector<Rectangle> rectangles;
  rectangles.reserve(BiologicalQuerySize);
  std::cout << "Beginning initialization of " << BiologicalQuerySize << " biological warfare rectangles..." << std::endl;
  for (unsigned i = 0; i < BiologicalQuerySize; ++i) {
    Rectangle rect;
    for (unsigned d = 0; d < dimensions; ++d) {
      fileGoodOrDie(file);
      file.read(buffer, sizeof(double));
      rect.lowerLeft[d] = *doubleBuffer;
      fileGoodOrDie(file);
      file.read(buffer, sizeof(double));
      rect.upperRight[d] = *doubleBuffer;
    }
    rectangles.push_back(rect);
  }
  std::cout << "Initialization OK." << std::endl;

  // Cleanup
  file.close();
  delete[] buffer;

  return rectangles;
}

static std::vector<Rectangle> generateForestRectangles() {
  // Query set is pre-generated and requires 5 dimensions
  assert(dimensions == 5);

  // Setup file reader and double buffer
  std::fstream file;
  std::string dataPath = "/home/bjglasbe/Documents/code/nir-tree/data/rea05.2";
  file.open(dataPath);
  fileGoodOrDie(file);
  char *buffer = new char[sizeof(double)];
  memset(buffer, 0, sizeof(double));
  double *doubleBuffer = (double *)buffer;

  // Initialize rectangles
  std::vector<Rectangle> rectangles;
  rectangles.reserve(ForestQuerySize);
  std::cout << "Beginning initialization of " << ForestQuerySize << " forest fire rectangles..." << std::endl;
  for (unsigned i = 0; i < ForestQuerySize; ++i) {
    Rectangle rect;
    for (unsigned d = 0; d < dimensions; ++d) {
      fileGoodOrDie(file);
      file.read(buffer, sizeof(double));
      rect.lowerLeft[d] = *doubleBuffer;
      fileGoodOrDie(file);
      file.read(buffer, sizeof(double));
      rect.upperRight[d] = *doubleBuffer;
    }

    rectangles.push_back(rect);
  }
  std::cout << "Initialization OK." << std::endl;

  // Cleanup
  file.close();
  delete[] buffer;

  return rectangles;
}

static std::vector<Rectangle> generateCanadaRectangles() {
  // Query set is pre-generated and requires 2 dimensions
  assert(dimensions == 2);

  // Setup file reader and double buffer
  std::fstream file;
  std::string dataPath = "/home/bjglasbe/Documents/code/nir-tree/data/canadaQ";
  file.open(dataPath);
  fileGoodOrDie(file);
  double bufferX, bufferY;

  // Initialize rectangles
  std::vector<Rectangle> rectangles;
  rectangles.reserve(CanadaQuerySize);
  std::cout << "Beginning initialization of " << CanadaQuerySize << " maple leaf rectangles..." << std::endl;
  for (unsigned i = 0; i < CanadaQuerySize; ++i) {
    Rectangle rect;
    // Read in lower left
    file >> bufferX;
    file >> bufferY;
    rect.lowerLeft[0] = bufferX;
    rect.lowerLeft[1] = bufferY;

    // Read in upper right
    file >> bufferX;
    file >> bufferY;
    rect.upperRight[0] = bufferX;
    rect.upperRight[1] = bufferY;

    rectangles.push_back(rect);
  }

  // Cleanup
  file.close();

  return rectangles;
}

static std::vector<Rectangle> generateGaiaRectangles() {
  // Query set is pre-generated and requires 3 dimensions
  assert(dimensions == 3);

  // Setup file reader and double buffer
  std::fstream file;
  std::string dataPath = "/home/bjglasbe/Documents/code/nir-tree/data/gaiaQ";
  file.open(dataPath);
  fileGoodOrDie(file);
  double buffer;

  // Initialize rectangles
  std::vector<Rectangle> rectangles;
  rectangles.reserve(GaiaQuerySize);
  std::cout << "Beginning initialization of " << GaiaQuerySize << " starry rectangles..." << std::endl;
  for (unsigned i = 0; i < GaiaQuerySize; ++i) {
    Rectangle rect;
    // Read in lower left
    file >> buffer;
    rect.lowerLeft[0] = buffer;
    file >> buffer;
    rect.lowerLeft[1] = buffer;
    file >> buffer;
    rect.lowerLeft[2] = buffer;

    // Read in upper right
    file >> buffer;
    rect.upperRight[0] = buffer;
    file >> buffer;
    rect.upperRight[1] = buffer;
    file >> buffer;
    rect.upperRight[2] = buffer;

    rectangles.push_back(rect);
  }

  // Cleanup
  file.close();

  return rectangles;
}

static bool is_already_loaded(std::map<std::string, uint64_t> &configU, Index *spatial_index) {
  if (configU["tree"] == NIR_TREE) {
    auto tree = (nirtreedisk::NIRTreeDisk<5, 9, nirtreedisk::ExperimentalStrategy> *) spatial_index;
    size_t existing_page_count = tree->node_allocator_->buffer_pool_.get_preexisting_page_count();

    if (existing_page_count > 0) {
      return true;
    }
  } else if (configU["tree"] == R_PLUS_TREE) {
    auto tree = (rplustreedisk::RPlusTreeDisk<5, 9> *) spatial_index;
    size_t existing_page_count = tree->node_allocator_.buffer_pool_.get_preexisting_page_count();

    if (existing_page_count > 0) {
      return true;
    }
  } else if (configU["tree"] == R_STAR_TREE) {
    auto tree = (rstartreedisk::RStarTreeDisk<5, 9> *) spatial_index;
    size_t existing_page_count = tree->node_allocator_->buffer_pool_.get_preexisting_page_count();

    if (existing_page_count > 0) {
      return true;
    }
  } else if (configU["tree"] == R_TREE) {
    auto tree = (rtreedisk::RTreeDisk<3, 6> *) spatial_index;
    size_t existing_page_count = tree->node_allocator_.buffer_pool_.get_preexisting_page_count();

    if (existing_page_count > 0) {
      return true;
    }
  }

  return false;
}

template <typename T>
static void
runBench(PointGenerator<T> &pointGen, std::map<std::string, uint64_t> &configU, std::map<std::string, double> &configD,
         std::map<std::string, std::string> &configS) {
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

  // Initialize the index
  Index *spatialIndex;
  if (configU["tree"] == R_TREE) {
    //spatialIndex = new rtree::RTree(configU["minfanout"], configU["maxfanout"]);
    spatialIndex = new rtreedisk::RTreeDisk<3, 6>(configU["buffer_pool_memory"], configS["db_file_name"]);
  } else if (configU["tree"] == R_PLUS_TREE) {
    spatialIndex = new rplustreedisk::RPlusTreeDisk<5, 9>(configU["buffer_pool_memory"], configS["db_file_name"]);
    //spatialIndex = new rplustree::RPlusTree(configU["minfanout"], configU["maxfanout"]);
  } else if (configU["tree"] == R_STAR_TREE) {
    //spatialIndex = new rstartree::RStarTree(configU["minfanout"], configU["maxfanout"]);
    spatialIndex = new rstartreedisk::RStarTreeDisk<5, 9>(configU["buffer_pool_memory"], configS["db_file_name"]);
  } else if (configU["tree"] == NIR_TREE) {
    //spatialIndex = new nirtree::NIRTree(configU["minfanout"], configU["maxfanout"]);
    //spatialIndex = new nirtree::NIRTree(5,9);
    spatialIndex = new nirtreedisk::NIRTreeDisk<5, 9, nirtreedisk::ExperimentalStrategy>(configU["buffer_pool_memory"], configS["db_file_name"]);
  } else if (configU["tree"] == QUAD_TREE) {
    spatialIndex = new quadtree::QuadTree();
  } else if (configU["tree"] == REVISED_R_STAR_TREE) {
    spatialIndex = new revisedrstartree::RevisedRStarTree(configU["minfanout"], configU["maxfanout"]);
  } else {
    std::cout << "Unknown tree selected. Exiting." << std::endl;
    return;
  }

  // Initialize search rectangles
  std::vector<Rectangle> searchRectangles;
  if (configU["distribution"] == UNIFORM) {
    searchRectangles = generateRectangles(configU["size"], configU["seed"], configU["rectanglescount"], configU["points_per_rectangle"]);
  } else if (configU["distribution"] == SKEW) {
    configU["rectanglescount"] = BitQuerySize;
    searchRectangles = generateBitRectangles();
  } else if (configU["distribution"] == CLUSTER) {
    configU["rectanglescount"] = HazeQuerySize;
    searchRectangles = generateHazeRectangles();
  } else if (configU["distribution"] == CALIFORNIA) {
    configU["rectanglescount"] = CaliforniaQuerySize;
    searchRectangles = generateCaliRectangles();
  } else if (configU["distribution"] == BIOLOGICAL) {
    configU["rectanglescount"] = BiologicalQuerySize;
    searchRectangles = generateBioRectangles();
  } else if (configU["distribution"] == FOREST) {
    configU["rectanglescount"] = ForestQuerySize;
    searchRectangles = generateForestRectangles();
  } else if (configU["distribution"] == CANADA) {
    configU["rectanglescount"] = CanadaQuerySize;
    searchRectangles = generateCanadaRectangles();
  } else if (configU["distribution"] == GAIA) {
    configU["rectanglescount"] = GaiaQuerySize;
    searchRectangles = generateGaiaRectangles();
  } else if (configU["distribution"] == ZIPF) {
    // This is not how we should generate rectangles for ZIPF.
    searchRectangles = generateRectangles(configU["size"], configU["seed"], configU["rectanglescount"], configU["points_per_rectangle"]);
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

  // Validate tree
  //spatialIndex->validate();
  //std::cout << "Validation OK." << std::endl;

  // Search for points and time their retrieval
  std::cout << "Beginning search." << std::endl;
  pointGen.reset();

#if 1
  while ((nextPoint = pointGen.nextPoint()) /* Intentional = not == */) {
    // Search
    Point p = nextPoint.value();
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

    if (totalSearches % 10000 == 0) {
      std::cout << "Point[" << totalSearches << "] queried. " << delta.count() << " s" << std::endl;
    }

    if (totalSearches >= configU["num_points_to_search"]) {
      break;
    }
  }

  std::cout << "Search OK." << std::endl;

#endif
  std::shuffle(searchRectangles.begin(), searchRectangles.end(), g);

#if 0
	// Search for rectangles
	unsigned rangeSearchChecksum = 0;
	std::cout << "Beginning search for " << searchRectangles.size() << " rectangles..." << std::endl;
	for (unsigned i = 0; i < searchRectangles.size(); ++i)
	{
		// Search
    std::cout << "Searching for: " << searchRectangles.at(i) << std::endl;
    std::chrono::high_resolution_clock::time_point begin = std::chrono::high_resolution_clock::now();
		std::vector<Point> v = spatialIndex->search(searchRectangles[i]);
		std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::cout << "Points: " << v.size() << std::endl;
		std::chrono::duration<double> delta = std::chrono::duration_cast<std::chrono::duration<double>>(end - begin);

    totalTimeRangeSearches += delta.count();
		totalRangeSearches += 1;
		rangeSearchChecksum += v.size();
		// std::cout << "searchRectangles[" << i << "] queried. " << delta.count() << " s" << std::endl;
		// std::cout << "searchRectangles[" << i << "] returned " << v.size() << " points" << std::endl;

	}
	std::cout << "Range search OK. Checksum = " << rangeSearchChecksum << std::endl;

#endif
  // Gather statistics

#ifdef STAT
  spatialIndex->stat();
  std::cout << "Statistics OK." << std::endl;
#endif

  // Timing Statistics
  std::cout << "Total time to insert: " << totalTimeInserts << "s" << std::endl;
  std::cout << "Avg time to insert: " << totalTimeInserts / (double)totalInserts << "s" << std::endl;
  std::cout << "Total searches: " << totalSearches << std::endl;
  std::cout << "Total time to search: " << totalTimeSearches << "s" << std::endl;
  std::cout << "Avg time to search: " << totalTimeSearches / totalSearches << "s" << std::endl;
  std::cout << "Total time to range search: " << totalTimeRangeSearches << "s" << std::endl;
  std::cout << "Avg time to range search: " << totalTimeRangeSearches / totalRangeSearches << "s" << std::endl;
  std::cout << "Total time to delete: " << totalTimeDeletes << "s" << std::endl;
  std::cout << "Avg time to delete: " << totalTimeDeletes / (double)totalDeletes << "s" << std::endl;

  // Generate visualization
  if (configU["visualization"]) {
    spatialIndex->print();
  }
}
