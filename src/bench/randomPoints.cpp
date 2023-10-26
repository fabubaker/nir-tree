#include <bench/randomPoints.h>

size_t BenchTypeClasses::Uniform::size = 10000;
unsigned BenchTypeClasses::Uniform::dimensions = dimensions;
unsigned BenchTypeClasses::Uniform::seed = 3141;
int BenchTypeClasses::Uniform::precision = -1;

size_t BenchTypeClasses::Zipf::size = 10000;
unsigned BenchTypeClasses::Zipf::dimensions = dimensions;
double BenchTypeClasses::Zipf::alpha = 1.0;
unsigned BenchTypeClasses::Zipf::seed = 3141;
unsigned BenchTypeClasses::Zipf::num_elements = 1000;

size_t BenchTypeClasses::Gauss::size = 10000;
unsigned BenchTypeClasses::Gauss::dimensions = dimensions;
unsigned BenchTypeClasses::Gauss::seed = 3141;

unsigned BenchTypeClasses::Zipf::binary_search(double needle, std::vector<double> &cummulative) {
  unsigned min_pos = 0;
  unsigned max_pos = num_elements;
  unsigned med_pos = 0;
  double elem = 0;

  while (min_pos <= max_pos) {
    med_pos = (min_pos + max_pos) / 2;
    if (med_pos == 0 || med_pos == UINT32_MAX) {
      break;
    }
    elem = cummulative.at(med_pos);
    if (elem == needle) {
      break;
    } else if (elem < needle) {
      min_pos = med_pos + 1;
    } else { // >
      if (min_pos == max_pos) {
        break;
      }
      max_pos = med_pos - 1;
    }
  }
  if (elem < needle) {
    unsigned old_pos = med_pos;
    med_pos = std::min(old_pos + 1, num_elements - 1);
  }
  return med_pos;
}

// 1/(position + 1)^alpha
double BenchTypeClasses::Zipf::invert_prob(unsigned position) {
  double p = pow((double)(position + 1), alpha);
  double inv_p = (1.0 / p);
  return inv_p;
}

void load_dataset(std::vector<Point> &points, std::string dataset_file_name) {
  std::ifstream file(dataset_file_name);

  if (!file.is_open()) {
    std::cerr << "Error: Couldn't open the file." << std::endl;
    abort();
  }

  std::string line;
  while (std::getline(file, line)) {
    std::ifstream line_stream(line.c_str());
    Point p;

    for (unsigned d = 0; d < DIM; ++d) {
      double dbl;
      line_stream >> dbl;
      p[d] = dbl;
    }

    points.push_back(p);
  }

  file.close();
}
