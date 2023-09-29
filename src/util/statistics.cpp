#include <ostream>
#include <vector>
#include <util/statistics.h>

void printFanoutHistogram(
        std::vector<std::vector<unsigned long>> histogramFanoutAtLevel,
        unsigned treeHeight
) {
  STATFANHIST();

  unsigned maxBucket = 0;
  // Get the max fanout for all levels
  for (unsigned lvl = 0; lvl < treeHeight; lvl++) {
    maxBucket = histogramFanoutAtLevel.at(lvl).size() > maxBucket ? histogramFanoutAtLevel.at(lvl).size() : maxBucket;
  }

  // Print histogram of fanout summed for all levels
  for (unsigned i = 0; i < maxBucket; i++) {
    unsigned totalCount = 0;
    // Sum up the count of branches with fanout=i at all levels
    for (unsigned lvl = 0; lvl < treeHeight; lvl++) {
      if (histogramFanoutAtLevel.at(lvl).size() > i && histogramFanoutAtLevel.at(lvl).at(i) > 0){
        totalCount += histogramFanoutAtLevel.at(lvl).at(i);
      }
    }
    if (totalCount > 0) {
      STATHIST(i, totalCount);
    }
  }

  STATFANHISTATLEVEL();

  // Print histogram of fanout at each level
  for (unsigned lvl = 0; lvl < treeHeight; lvl++) {
    STATEXEC(std::cout << "=== LEVEL: " << lvl << " ===" << std::endl);
    for (unsigned i = 0; i < histogramFanoutAtLevel.at(lvl).size(); i++) {
      if (histogramFanoutAtLevel.at(lvl).at(i) > 0) {
        STATHISTATLEVEL(i, histogramFanoutAtLevel.at(lvl).at(i), lvl);
      }
    }
  }
}

void printPolygonHistogram(
    std::vector<std::vector<unsigned long>> histogramPolygonAtLevel,
    unsigned treeHeight
) {
  STATPOLYHIST();

  unsigned maxBucket = 0;
  // Get the max polygon size for all levels
  for (unsigned lvl = 0; lvl < treeHeight; lvl++) {
    maxBucket = histogramPolygonAtLevel.at(lvl).size() > maxBucket ? histogramPolygonAtLevel.at(lvl).size() : maxBucket;
  }

  // Print histogram of polygon size summed for all levels
  for (unsigned i = 0; i < maxBucket; i++) {
    unsigned totalCount = 0;
    // Sum up the count of polygon with size i at all levels
    for (unsigned lvl = 0; lvl < treeHeight; lvl++) {
      if (histogramPolygonAtLevel.at(lvl).size() > i && histogramPolygonAtLevel.at(lvl).at(i) > 0){
        totalCount += histogramPolygonAtLevel.at(lvl).at(i);
      }
    }
    if (totalCount > 0){
      STATHIST(i, totalCount);
    }
  }

  STATPOLYHISTATLEVEL();

  // Print histogram of polygon size at each level
  for (unsigned lvl = 0; lvl < treeHeight; lvl++) {
    STATEXEC(std::cout << "=== LEVEL: " << lvl << " ===" << std::endl);
    for (unsigned i = 0; i < histogramPolygonAtLevel.at(lvl).size(); i++) {
      if (histogramPolygonAtLevel.at(lvl).at(i) > 0) {
        STATHISTATLEVEL(i, histogramPolygonAtLevel.at(lvl).at(i), lvl);
      }
    }
  }
}