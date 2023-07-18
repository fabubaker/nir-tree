#pragma once

// Copyright 2021 Kyle Langendoen

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//     http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstring>
#include <globals/globals.h>
#include <iostream>
#include <limits>
#include <list>
#include <omp.h>
#include <queue>
#include <stack>
#include <unordered_map>
#include <util/compression.h>
#include <util/debug.h>
#include <util/geometry.h>
#include <util/graph.h>
#include <util/repacking.h>
#include <util/statistics.h>
#include <utility>
#include <variant>
#include <vector>
#include <fstream>

namespace nirtreedisk {
struct BranchPartitionStrategy {};

struct LineMinimizeDownsplits : BranchPartitionStrategy {};

struct LineMinimizeDistanceFromMean : BranchPartitionStrategy {};

struct ExperimentalStrategy : BranchPartitionStrategy {};

template <int min_branch_factor, int max_branch_factor, class strategy>
class NIRTreeDisk;

template <int min_branch_factor, int max_branch_factor, class strategy>
requires(std::derived_from<strategy, BranchPartitionStrategy>)
tree_node_allocator *get_node_allocator(NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef) {
  return treeRef->node_allocator_.get();
}

struct Branch {
  Branch(Rectangle boundingBox, tree_node_handle child_handle):
    boundingBox(boundingBox), child(child_handle) {}

  Branch() = default;

  Branch(const Branch &other): boundingBox(other.boundingBox), child(other.child) {}

  bool operator==(const Branch &o) const = default;
  bool operator!=(const Branch &o) const = default;

  Rectangle boundingBox;
  tree_node_handle child;
};

struct SplitResult {
  Branch leftBranch;
  Branch rightBranch;
};

struct Partition {
  unsigned dimension;
  double location;
};

template <int min_branch_factor, int max_branch_factor, class strategy>
requires(std::derived_from<strategy, BranchPartitionStrategy>)
class LeafNode {
public:
  std::array<Point, max_branch_factor> entries;
  unsigned cur_offset_;

  LeafNode(): cur_offset_(0) {
    static_assert(sizeof(LeafNode<min_branch_factor, max_branch_factor, strategy>) <= PAGE_DATA_SIZE);
  }

  void addPoint(const Point &point) {
    entries.at(this->cur_offset_++) = point;
  }

  void deleteSubtrees();

  // Helper functions
  Rectangle boundingBox();

  void removePoint(const Point &point);
  void removeEntry(const tree_node_handle &handle);

  tree_node_handle chooseNode(Point givenPoint);
  tree_node_handle findLeaf(Point givenPoint);

  Partition partitionNode();
  Partition partitionLeafNode();

  void make_disjoint_from_children(IsotheticPolygon &polygon,
                                   tree_node_handle handle_to_skip);
  SplitResult splitNode(Partition p, bool is_downsplit);
  SplitResult splitNode();
  SplitResult adjustTree(std::vector<bool> &hasReinsertedOnLevel);
  void condenseTree();

  // Data structure interface functions
  tree_node_handle insert(Point givenPoint, std::vector<bool>
                                                &hasReinsertedOnLevel);
  void reInsert(std::vector<bool> &hasReinsertedOnLevel);

  tree_node_handle remove(Point givenPoint);

  // Miscellaneous
  unsigned checksum();
  bool validate(tree_node_handle expectedParent, unsigned index);
  std::vector<Point> bounding_box_validate();
  void print(unsigned n = 0);
  void printTree(unsigned n = 0);
  unsigned height();
};

template <int min_branch_factor, int max_branch_factor, class strategy>
requires(std::derived_from<strategy, BranchPartitionStrategy>)
class BranchNode {
public:
  std::array<Branch, max_branch_factor> entries;
  unsigned cur_offset_;

  // Constructors and destructors
  BranchNode(): cur_offset_(0) {}

  void deleteSubtrees();

  // Helper functions
  Rectangle boundingBox();
  Branch &locateBranch(tree_node_handle child) {
    for (unsigned i = 0; i < this->cur_offset_; i++) {
      Branch &b = entries.at(i);
      if (b.child == child) {
        return b;
      }
    }
    // If we are here, panic
    assert(false);
  };

  void updateBranch(tree_node_handle child, const InlineBoundedIsotheticPolygon &boundingPoly);
  void removeBranch(const tree_node_handle &handle);

  void addBranchToNode(const Branch &entry) {
    entries.at(this->cur_offset_++) = entry;
  }

  tree_node_handle chooseNode(std::variant<Branch, Point> &nodeEntry, uint8_t stopping_level);
  tree_node_handle findLeaf(Point givenPoint);
  Partition partitionNode();

  template <class S = strategy>
  Partition partitionBranchNode(typename std::enable_if<std::is_same<S, LineMinimizeDownsplits>::value, S>::type * = 0) {
#if 0
    Partition defaultPartition;

    tree_node_allocator *allocator = get_node_allocator(this->treeRef);
    std::vector<Rectangle> all_branch_polys;
    for (unsigned i = 0; i < this->cur_offset_; i++) {
      Branch &b_i = entries.at(i);
      all_branch_polys.push_back(b_i.get_summary_rectangle(allocator));
    }

    unsigned cut_off = min_branch_factor;
    double best_candidate = 0.0;
    double min_cost = std::numeric_limits<double>::max();
    unsigned best_dimension = 0;
    do {
      // D * ( M LOG M + M ) ~> O( D M LOG M )
      for (unsigned d = 0; d < dimensions; d++) {
        std::sort(all_branch_polys.begin(), all_branch_polys.end(),
                  [d](Rectangle &poly1, Rectangle &poly2) {
                    return poly1.upperRight[d] <
                           poly2.upperRight[d];
                  });

        std::vector<double> partition_candidates;

        for (unsigned i = 0; i < all_branch_polys.size(); i++) {
          partition_candidates.push_back(all_branch_polys.at(i).lowerLeft[d]);
          partition_candidates.push_back(all_branch_polys.at(i).upperRight[d]);
        }

        for (double partition_candidate : partition_candidates) {
          double cost = 0;
          // starts at 1 cause first goes left
          // Technically we should also walk the bottom bounds to
          // be sure, even in the non F, C case.
          unsigned left_count = 0;
          unsigned right_count = 0;
          double running_total = 0.0;
          // Existing metric wanted to avoid recursive splits
          // Let's try and do the same
          for (unsigned j = 0; j < all_branch_polys.size(); j++) {
            Rectangle &poly_ref = all_branch_polys.at(j);
            running_total += poly_ref.lowerLeft[d] +
                             poly_ref.upperRight[d];

            bool greater_than_left = poly_ref.lowerLeft[d] <
                                     partition_candidate;
            bool less_than_right = partition_candidate <
                                   poly_ref.upperRight[d];
            bool requires_split = greater_than_left and
                                  less_than_right;

            bool should_go_left = poly_ref.upperRight[d] <=
                                  partition_candidate;
            bool should_go_right = poly_ref.lowerLeft[d] >=
                                   partition_candidate;

            if (requires_split) {
              left_count++;
              right_count++;
              cost++;
            } else if (should_go_left) {
              // the zero area polys end up here too
              left_count++;
            } else if (should_go_right) {
              right_count++;
            }
          }
          if (cost < min_cost and left_count <= max_branch_factor and right_count <= max_branch_factor and
              left_count >= cut_off and right_count >= cut_off) {
            best_candidate = partition_candidate;
            best_dimension = d;
            min_cost = cost;
          }
        }
      }
      if (cut_off == 1 and min_cost ==
                               std::numeric_limits<double>::max()) {
        abort();
      }
      cut_off = 1;
    } while (min_cost == std::numeric_limits<double>::max());
    // Degenerate case
    assert(min_cost < std::numeric_limits<double>::max());

    defaultPartition.dimension = best_dimension;
    defaultPartition.location = best_candidate;

    return defaultPartition;
#endif

    // Unsupported
    abort();
  }

  template <class S = strategy>
  Partition partitionBranchNode(typename std::enable_if<std::is_same<S, LineMinimizeDistanceFromMean>::value,
                                                        S>::type * = 0) {
#if 0
    Partition defaultPartition;

    tree_node_allocator *allocator = get_node_allocator(
        this->treeRef);
    std::vector<Rectangle> all_branch_polys;
    for (unsigned i = 0; i < this->cur_offset_; i++) {
      Branch &b_i = entries.at(i);
      all_branch_polys.push_back(b_i.get_summary_rectangle(
          allocator));
    }

    double best_candidate = 0.0;
    double min_cost = std::numeric_limits<double>::max();
    unsigned best_dimension = 0;
    // D * ( M LOG M + M ) ~> O( D M LOG M )
    for (unsigned d = 0; d < dimensions; d++) {
      std::sort(all_branch_polys.begin(), all_branch_polys.end(),
                [d](Rectangle &poly1, Rectangle &poly2) {
                  if (poly1.upperRight[d] ==
                      poly2.upperRight[d]) {
                    for (unsigned i = 0; i < dimensions;
                         i++) {
                      if (poly1.upperRight[i] ==
                          poly2.upperRight[i]) {
                        continue;
                      }
                      return poly1.upperRight[i] <
                             poly2.upperRight[i];
                    }
                  }
                  return poly1.upperRight[d] <
                         poly2.upperRight[d];
                });
      for (unsigned i = 0; i < all_branch_polys.size(); i++) {
        double cost = 0;
        // starts at 1 cause {i} goes left
        // Technically we should also walk the bottom bounds to
        // be sure, even in the non F, C case.
        unsigned left_count = 0;
        unsigned right_count = 0;
        double partition_candidate =
            all_branch_polys.at(i).upperRight[d];
        double running_total = 0.0;
        // Existing metric wanted to avoid recursive splits
        // Let's try and do the same
        for (unsigned j = 0; j < all_branch_polys.size(); j++) {
          Rectangle &poly_ref = all_branch_polys.at(j);
          running_total += poly_ref.lowerLeft[d] +
                           poly_ref.upperRight[d];

          bool greater_than_left = poly_ref.lowerLeft[d] <
                                   partition_candidate;
          bool less_than_right = partition_candidate <
                                 poly_ref.upperRight[d];
          bool requires_split = greater_than_left and
                                less_than_right;

          bool should_go_left = poly_ref.upperRight[d] <=
                                partition_candidate;
          bool should_go_right = poly_ref.lowerLeft[d] >=
                                 partition_candidate;
          bool is_zero_area =
              poly_ref.lowerLeft[d] ==
              poly_ref.upperRight[d];

          if (requires_split) {
            left_count++;
            right_count++;
            cost++;
          } else if (is_zero_area and
                     poly_ref.upperRight[d] ==
                         partition_candidate) {
            // Partition on a zero-area thing, can
            // pick either side as convenient
            if (left_count <= right_count) {
              left_count++;
            } else {
              right_count++;
            }
          } else if (should_go_left) {
            left_count++;
          } else if (should_go_right) {
            right_count++;
          } else {
            assert(false);
          }
        }

        // Cost function 2
        // If the split will not overflow our children
        if (left_count <= max_branch_factor and right_count <= max_branch_factor and
            left_count > 0 and right_count > 0) {
          double mean_d_pt = running_total /
                             (2 * all_branch_polys.size());
          // Distance
          cost = (mean_d_pt - partition_candidate);
          cost = cost * cost;
          if (cost < min_cost) {
            best_candidate = partition_candidate;
            best_dimension = d;
            min_cost = cost;
          }
        }
      }
    }
    // Degenerate case
    assert(min_cost < std::numeric_limits<double>::max());

    defaultPartition.dimension = best_dimension;
    defaultPartition.location = best_candidate;

    // Sort per the dimension we need
    std::sort(entries.begin(), entries.begin() + this->cur_offset_,
              [this, allocator, best_dimension](Branch &b1, Branch &b2) {
                Rectangle poly1 = b1.get_summary_rectangle(
                    allocator);
                Rectangle poly2 = b2.get_summary_rectangle(
                    allocator);
                if (poly1.upperRight[best_dimension] ==
                    poly2.upperRight[best_dimension]) {
                  for (unsigned i = 0; i < dimensions; i++) {
                    if (poly1.upperRight[i] ==
                        poly2.upperRight[i]) {
                      continue;
                    }
                    return poly1.upperRight[i] <
                           poly2.upperRight[i];
                  }
                }
                return poly1.upperRight[best_dimension] <
                       poly2.upperRight[best_dimension];
              });

    return defaultPartition;
#endif

    // Unsupported
    abort();
  }

  std::pair<bool, Partition> try_cut_geo_mean(std::vector<Rectangle> &all_branch_polys) {
#if 0
    Partition defaultPartition;
    Point mean_point = Point::atOrigin;

    double mass = 0.0;
    for (auto &branch_bounding_box : all_branch_polys) {
      mean_point += branch_bounding_box.lowerLeft;
      mean_point += branch_bounding_box.upperRight;
      mass += 2.0;
    }

    mean_point /= mass;

    unsigned best_cost =
        std::numeric_limits<unsigned>::max();

    // Need to determine left and right count

    for (unsigned d = 0; d < dimensions; d++) {
      // Is this a valid split?
      double location = mean_point[d];
      unsigned cost = 0;
      unsigned left_count = 0;
      unsigned right_count = 0;
      for (auto &branch_bounding_box : all_branch_polys) {
        bool greater_than_left = branch_bounding_box.lowerLeft[d] <
                                 location;
        bool less_than_right = location <
                               branch_bounding_box.upperRight[d];
        bool requires_split = greater_than_left and
                              less_than_right;

        bool should_go_left = branch_bounding_box.upperRight[d] <= location;
        bool should_go_right = branch_bounding_box.lowerLeft[d] >= location;
        assert(not(should_go_left and should_go_right));

        if (requires_split) {
          left_count++;
          right_count++;
          cost++;
        } else if (should_go_left) {
          left_count++;
        } else if (should_go_right) {
          right_count++;
        } else {
          assert(false);
        }
      }

      if (left_count > 0 and right_count > 0 and
          left_count <= max_branch_factor and
          right_count <= max_branch_factor) {
        if (cost < best_cost) {
          best_cost = cost;
          defaultPartition.location = mean_point[d];
          defaultPartition.dimension = d;
        }
      }
    }

    if (best_cost < std::numeric_limits<unsigned>::max()) {
      return std::make_pair(true, defaultPartition);
    }
    return std::make_pair(false, defaultPartition);
#endif

    // Unsupported
    abort();
  }

  template <class S = strategy>
  Partition partitionBranchNode(typename std::enable_if<std::is_same<S, ExperimentalStrategy>::value,
                                                        S>::type * = 0) {
#if 0
    Partition defaultPartition;

    tree_node_allocator *allocator = get_node_allocator(
        this->treeRef);
    std::vector<Rectangle> all_branch_polys;
    for (unsigned i = 0; i < this->cur_offset_; i++) {
      Branch &b_i = entries.at(i);
      all_branch_polys.push_back(b_i.get_summary_rectangle(
          allocator));
    }

    auto geo_cut = try_cut_geo_mean(all_branch_polys);
    if (geo_cut.first) {
      return geo_cut.second;
    }
    // Can we cut along the geometric mean in any dimension
    // without overflowing our children?

    // If that didn't work, we gotta try something else.
    for (unsigned d = 0; d < dimensions; d++) {
      std::sort(all_branch_polys.begin(), all_branch_polys.end(),
                [d](Rectangle &poly1, Rectangle &poly2) {
                  if (poly1.upperRight[d] ==
                      poly2.upperRight[d]) {
                    for (unsigned i = 0; i < dimensions;
                         i++) {
                      if (poly1.upperRight[i] ==
                          poly2.upperRight[i]) {
                        continue;
                      }
                      return poly1.upperRight[i] <
                             poly2.upperRight[i];
                    }
                  }
                  return poly1.upperRight[d] <
                         poly2.upperRight[d];
                });
    }

    double best_candidate = 0.0;
    double min_cost = std::numeric_limits<double>::max();
    unsigned best_dimension = 0;
    // D * ( M LOG M + M ) ~> O( D M LOG M )
    for (unsigned d = 0; d < dimensions; d++) {
      std::sort(all_branch_polys.begin(), all_branch_polys.end(),
                [d](Rectangle &poly1, Rectangle &poly2) {
                  if (poly1.upperRight[d] ==
                      poly2.upperRight[d]) {
                    for (unsigned i = 0; i < dimensions;
                         i++) {
                      if (poly1.upperRight[i] ==
                          poly2.upperRight[i]) {
                        continue;
                      }
                      return poly1.upperRight[i] <
                             poly2.upperRight[i];
                    }
                  }
                  return poly1.upperRight[d] <
                         poly2.upperRight[d];
                });
      for (unsigned i = 0; i < all_branch_polys.size(); i++) {
        double cost = 0;
        // starts at 1 cause {i} goes left
        // Technically we should also walk the bottom bounds to
        // be sure, even in the non F, C case.
        unsigned left_count = 0;
        unsigned right_count = 0;
        double partition_candidate =
            all_branch_polys.at(i).upperRight[d];
        double running_total = 0.0;
        // Existing metric wanted to avoid recursive splits
        // Let's try and do the same
        for (unsigned j = 0; j < all_branch_polys.size(); j++) {
          Rectangle &poly_ref = all_branch_polys.at(j);
          running_total += poly_ref.lowerLeft[d] +
                           poly_ref.upperRight[d];

          bool greater_than_left = poly_ref.lowerLeft[d] <
                                   partition_candidate;
          bool less_than_right = partition_candidate <
                                 poly_ref.upperRight[d];
          bool requires_split = greater_than_left and
                                less_than_right;

          bool should_go_left = poly_ref.upperRight[d] <= partition_candidate;
          bool should_go_right = poly_ref.lowerLeft[d] >= partition_candidate;
          assert(not(should_go_left and
                     should_go_right));
          bool is_zero_area =
              poly_ref.lowerLeft[d] ==
              poly_ref.upperRight[d];

          if (requires_split) {
            //std::cout << "SIMUL: entry requires split." << std::endl;
            left_count++;
            right_count++;
            cost++;
          } else if (is_zero_area and
                     poly_ref.upperRight[d] ==
                         partition_candidate) {
            assert(false);
            // Partition on a zero-area thing, can
            // pick either side as convenient
            if (left_count <= right_count) {
              //std::cout << "SIMUL: entry contest, goes left." << std::endl;
              left_count++;
            } else {
              //std::cout << "SIMUL: entry contest, goes right." << std::endl;
              right_count++;
            }
          } else if (should_go_left) {
            //std::cout << "SIMUL: entry goes left." << std::endl;
            left_count++;
          } else if (should_go_right) {
            //std::cout << "SIMUL: entry goes right." << std::endl;
            right_count++;
          } else {
            assert(false);
          }
        }

        // Cost function 2
        // If the split will not overflow our children
        if (left_count <= max_branch_factor and right_count <= max_branch_factor and
            left_count > 0 and right_count > 0) {
          double mean_d_pt = running_total /
                             (2 * all_branch_polys.size());
          // Distance
          cost = (mean_d_pt - partition_candidate);
          cost = cost * cost;
          if (cost < min_cost) {
            best_candidate = partition_candidate;
            best_dimension = d;
            min_cost = cost;
            //std::cout << "Best Candidate LC: " <<
            //    left_count << " and RC: " <<
            //    right_count << std::endl;
          }
        }
      }
    }
    // Degenerate case
    assert(min_cost < std::numeric_limits<double>::max());

    defaultPartition.dimension = best_dimension;
    defaultPartition.location = best_candidate;

    // Sort per the dimension we need
    std::sort(entries.begin(), entries.begin() + this->cur_offset_,
              [this, allocator, best_dimension](Branch &b1, Branch &b2) {
                Rectangle poly1 = b1.get_summary_rectangle(
                    allocator);
                Rectangle poly2 = b2.get_summary_rectangle(
                    allocator);
                if (poly1.upperRight[best_dimension] ==
                    poly2.upperRight[best_dimension]) {
                  for (unsigned i = 0; i < dimensions; i++) {
                    if (poly1.upperRight[i] ==
                        poly2.upperRight[i]) {
                      continue;
                    }
                    return poly1.upperRight[i] <
                           poly2.upperRight[i];
                  }
                }
                return poly1.upperRight[best_dimension] <
                       poly2.upperRight[best_dimension];
              });

    return defaultPartition;
#endif

    // Unsupported
    abort();
  }

  void make_disjoint_from_children(IsotheticPolygon &polygon, tree_node_handle handle_to_skip);
  SplitResult splitNode(Partition p, bool is_downsplit);
  SplitResult splitNode();
  SplitResult adjustTree();
  void condenseTree();

  // Data structure interface functions
  void exhaustiveSearch(Point &requestedPoint, std::vector<Point> &accumulator);
  tree_node_handle insert(std::variant<Branch, Point> &nodeEntry, std::vector<bool> &hasReinsertedOnLevel);
  void reInsert(std::vector<bool> &hasReinsertedOnLevel);
  tree_node_handle remove(Point givenPoint);

  // Miscellaneous
  unsigned checksum();
  bool validate(tree_node_handle expectedParent, unsigned index);
  std::vector<Point> bounding_box_validate();
  void print(unsigned n = 0);
  void printTree(unsigned n = 0);
  unsigned height();

  std::pair<uint16_t, std::vector<std::optional<std::pair<char *, int>>>>
  compute_packed_size(
          tree_node_allocator *existing_allocator, tree_node_allocator *new_allocator,
          unsigned &maximum_repacked_rect_size
  );
  tree_node_handle repack(
          tree_node_allocator *existing_allocator,
          tree_node_allocator *new_allocator
  );
};

template <int min_branch_factor, int max_branch_factor, class strategy>
void LeafNode<min_branch_factor, max_branch_factor, strategy>::deleteSubtrees() {
  return;
}

template <int min_branch_factor, int max_branch_factor, class strategy>
Rectangle LeafNode<min_branch_factor, max_branch_factor, strategy>::boundingBox() {
  assert(this->cur_offset_ > 0);

  if (this->cur_offset_ <= 0) {
    std::cout << "Node has no entries!" << std::endl;
    abort();
  }

  Point &p = entries.at(0);
  Rectangle bb(p, Point::closest_larger_point(p));
  for (unsigned i = 1; i < this->cur_offset_; i++) {
    bb.expand(entries.at(i));
  }
  return bb;
}

// Point: removes the point from the node
template <int min_branch_factor, int max_branch_factor, class strategy>
void LeafNode<min_branch_factor, max_branch_factor, strategy>::removePoint(const Point &point) {
  // Locate the child
  unsigned childIndex;
  for (childIndex = 0; entries.at(childIndex) != point and
                       childIndex < this->cur_offset_;
       childIndex++) {
  }
  assert(entries.at(childIndex) == point);

  // Replace this index with whatever is in the last position
  entries.at(childIndex) = entries.at(this->cur_offset_ - 1);

  // Truncate array size
  this->cur_offset_--;
}

template <typename iter>
void shrink(IsotheticPolygon &polygon, iter begin, iter end, tree_node_allocator *allocator) {
  // Early exit
  if (polygon.basicRectangles.size() == 0 || begin == end) {
    return;
  }

  std::vector<Rectangle> rectangleSetShrunk;
  for (const Rectangle &basicRectangle : polygon.basicRectangles) {
    bool addRectangle = false;
    Rectangle shrunkRectangle = Rectangle(Point::atInfinity, Point::atNegInfinity);
    for (auto cur_iter = begin; cur_iter != end; cur_iter++) {
      Point &pinPoint = *cur_iter;
      if (basicRectangle.containsPoint(pinPoint)) {
        shrunkRectangle.expand(pinPoint);
        addRectangle = true;
        assert(shrunkRectangle.containsPoint(pinPoint));
      }
    }

    if (addRectangle) {
      rectangleSetShrunk.emplace_back(std::move(shrunkRectangle));
    }
  }

  assert(rectangleSetShrunk.size() > 0);

  polygon.basicRectangles.swap(rectangleSetShrunk);
  polygon.recomputeBoundingBox();
}

// Always called on root, this = root
// This top-to-bottom sweep is only for adjusting bounding boxes to contain the point and
// choosing a particular leaf
template <int min_branch_factor, int max_branch_factor, class strategy>
tree_node_handle LeafNode<min_branch_factor, max_branch_factor, strategy>::chooseNode(Point givenPoint) {
  return this->self_handle_;
}

template <int min_branch_factor, int max_branch_factor, class strategy>
tree_node_handle LeafNode<min_branch_factor, max_branch_factor, strategy>::findLeaf(Point givenPoint) {
  // FL2 [Search leaf node for record]
  // Check each entry to see if it matches E
  for (unsigned i = 0; i < this->cur_offset_; i++) {
    Point &p = std::get<Point>(entries.at(i));
    if (p == givenPoint) {
      return this->self_handle_;
    }
  }
  return tree_node_handle(nullptr);
}

template <int min_branch_factor, int max_branch_factor, class strategy>
Partition LeafNode<min_branch_factor, max_branch_factor, strategy>::partitionLeafNode() {
  Partition defaultPartition;
  double totalMass = 0.0;
  Point variance = Point::atOrigin;
  Point average = Point::atOrigin;
  Point sumOfSquares = Point::atOrigin;

  for (unsigned i = 0; i < this->cur_offset_; i++) {
    Point &dataPoint = entries.at(i);
    average += dataPoint;
    sumOfSquares += dataPoint * dataPoint;
    totalMass += 1.0;
  }

  // Compute final terms
  average /= totalMass;
  sumOfSquares /= totalMass;

  // Compute final variance
  variance = sumOfSquares - average * average;

  // Choose most variate dimension
  defaultPartition.dimension = 0;
  for (unsigned d = 0; d < dimensions; d++) {
    if (variance[d] > variance[defaultPartition.dimension]) {
      defaultPartition.dimension = d;
    }
  }
  defaultPartition.location = average[defaultPartition.dimension];

  return defaultPartition;
}

template <int min_branch_factor, int max_branch_factor, class strategy>
Partition LeafNode<min_branch_factor, max_branch_factor, strategy>::partitionNode() {
  return partitionLeafNode();
}

// We create two new nodes and free the old one.
// The old one is freed in adjustTree using removeEntry
// If we downsplit, then we won't call adjustTree for that split so we
// need to delete the node ourselves.
template <int min_branch_factor, int max_branch_factor, class strategy>
SplitResult LeafNode<min_branch_factor, max_branch_factor, strategy>::splitNode(Partition p, bool is_downsplit) {
#if 0
  tree_node_allocator *allocator = get_node_allocator(this->treeRef);

  auto alloc_data =
          allocator->create_new_tree_node<LeafNode<min_branch_factor, max_branch_factor, strategy>>(
                  NodeHandleType(LEAF_NODE));
  tree_node_handle left_handle = alloc_data.second;
  auto left_node = alloc_data.first; // take pin
  new (&(*left_node)) LeafNode<min_branch_factor, max_branch_factor, strategy>(this->treeRef,
                                            this->parent, left_handle, level_);
  assert(left_node->self_handle_ == left_handle);
  assert(left_handle.get_type() == LEAF_NODE);

  alloc_data = allocator->create_new_tree_node<LeafNode<min_branch_factor, max_branch_factor, strategy>>(
          NodeHandleType(LEAF_NODE));
  tree_node_handle right_handle = alloc_data.second;
  auto right_node = alloc_data.first; // take pin
  new (&(*right_node)) LeafNode<min_branch_factor, max_branch_factor, strategy>(this->treeRef,
                                             this->parent, right_handle, level_);
  assert(right_node->self_handle_ == right_handle);
  assert(right_handle.get_type() == LEAF_NODE);

  SplitResult split = {
          {tree_node_handle(nullptr), left_handle},
          {tree_node_handle(nullptr), right_handle}};

  bool containedLeft, containedRight;
  for (unsigned i = 0; i < this->cur_offset_; i++) {
    Point &dataPoint = entries.at(i);
    containedLeft = dataPoint[p.dimension] < p.location; // Not inclusive
    containedRight = dataPoint[p.dimension] >= p.location;
    assert(containedLeft or containedRight);

    if (containedLeft and not containedRight) {
      left_node->addPoint(dataPoint);
    } else if (not containedLeft and containedRight) {
      right_node->addPoint(dataPoint);
    }
  }

  // All points have been routed.

  IsotheticPolygon left_polygon(left_node->boundingBox());
  IsotheticPolygon right_polygon(right_node->boundingBox());
  assert(left_polygon.disjoint(right_polygon));

  assert(left_polygon.basicRectangles.size() <=
         MAX_RECTANGLE_COUNT);
  assert(right_polygon.basicRectangles.size() <=
         MAX_RECTANGLE_COUNT);

  // If we have a parent, we need to make these disjoint from our
  // siblings. If we don't, then we are automatically disjoint
  // from our siblings since these arethe only two polys and they
  // are disjoint from each other now.
  if (this->parent) {
    auto parent_node = this->treeRef->get_branch_node(this->parent);
    if (not is_downsplit) {
      parent_node->make_disjoint_from_children(left_polygon,
                                               this->self_handle_);
      assert(left_polygon.basicRectangles.size() > 0);
      left_polygon.refine();
      assert(left_polygon.basicRectangles.size() > 0);
      parent_node->make_disjoint_from_children(right_polygon,
                                               this->self_handle_);
      assert(right_polygon.basicRectangles.size() > 0);
      right_polygon.refine();
      assert(right_polygon.basicRectangles.size() > 0);
    } else {
      // Intersect with our existing poly to avoid intersect
      // with other children
      Branch &b = parent_node->locateBranch(this->self_handle_);
      IsotheticPolygon parent_poly = b.materialize_polygon(
              allocator);
      assert(parent_poly.basicRectangles.size() > 0);
      assert(left_polygon.basicRectangles.size() > 0);
      IsotheticPolygon poly_backup = left_polygon;
      left_polygon.intersection(parent_poly);
      if (left_polygon.basicRectangles.size() == 0) {
        std::cout << "Weird situation: " << poly_backup << " is disjoint from parent: " << parent_poly << std::endl;
      }
      assert(left_polygon.basicRectangles.size() > 0);
      left_polygon.refine();
      assert(left_polygon.basicRectangles.size() > 0);
      assert(right_polygon.basicRectangles.size() > 0);
      right_polygon.intersection(parent_poly);
      assert(right_polygon.basicRectangles.size() > 0);
      right_polygon.refine();
      assert(right_polygon.basicRectangles.size() > 0);
    }
  }
  this->cur_offset_ = 0;

  update_branch_polygon(split.leftBranch, left_polygon, treeRef, true);
  update_branch_polygon(split.rightBranch, right_polygon, treeRef, true);

  return split;
#endif

  // Unsupported
  abort();
}

// Splitting a node will remove it from its this->parent node and its memory will be freed
template <int min_branch_factor, int max_branch_factor, class strategy>
SplitResult LeafNode<min_branch_factor, max_branch_factor, strategy>::splitNode() {
      SplitResult returnSplit = splitNode(partitionNode(), false);
      return returnSplit;
}

template <int min_branch_factor, int max_branch_factor, class strategy>
void LeafNode<min_branch_factor, max_branch_factor, strategy>::reInsert(std::vector<bool> &hasReinsertedOnLevel) {
#if 0
  // Taken from R*
  hasReinsertedOnLevel.at(level_) = true;

  Point globalCenterPoint = boundingBox().centrePoint();

  std::sort(entries.begin(), entries.begin() + cur_offset_,
            [&globalCenterPoint](Point &a, Point &b) {
                Rectangle rectA(a, Point::closest_larger_point(a));
                Rectangle rectB(b, Point::closest_larger_point(b));
                return rectA.centrePoint().distance(globalCenterPoint) > rectB.centrePoint().distance(globalCenterPoint);
            });

  unsigned numNodesToReinsert = 0.3 * cur_offset_;
  unsigned remainder = cur_offset_ - numNodesToReinsert;

  std::vector<Point> entriesToReinsert;
  entriesToReinsert.reserve(numNodesToReinsert);

  std::copy(
          entries.begin() + remainder,
          entries.begin() + cur_offset_,
          std::back_inserter(entriesToReinsert));

  cur_offset_ = remainder;

  // These points will endup in exactly the same spot unless we
  // touch up the bounding boxes on the way up the tree to the root.
  // Each point path is unique, which means we would end up in the
  // same spot --- but it doesn't mean it is a good spot!

  // FIXME: Yolo for now, under the assumption that our siblings are
  // always disjoint so it should be fine.
  // Need to really think about whether this is correct
  // We want to adjust the paths on the way up to precisely reflect
  // what region we are in.

  fix_up_path_polys(self_handle_, treeRef);
  tree_node_handle root_handle = treeRef->root;

  for (const Point &entry : entriesToReinsert) {
    if (root_handle.get_type() == LEAF_NODE || root_handle.get_type() == REPACKED_LEAF_NODE) {

      auto root_node = treeRef->get_leaf_node(root_handle);
      root_handle = root_node->insert(entry, hasReinsertedOnLevel);
    } else {
      auto root_node = treeRef->get_branch_node(root_handle);
      std::variant<Branch, Point> ent = entry;
      root_handle = root_node->insert(ent, hasReinsertedOnLevel);
    }
  }
#endif

  // Unsupported
  abort();
}

template <int min_branch_factor, int max_branch_factor, class strategy>
IsotheticPolygon get_polygon_path_constraints(
        tree_node_handle start_handle,
        NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef)
{
#if 0
  tree_node_allocator *allocator = get_node_allocator(treeRef);
  tree_node_handle parent_handle;
  if (start_handle.get_type() == LEAF_NODE || start_handle.get_type() == REPACKED_LEAF_NODE) {
    auto leaf_node = treeRef->get_leaf_node(start_handle);
    parent_handle = leaf_node->parent;
  } else {
    auto branch_node = treeRef->get_branch_node(start_handle);
    parent_handle = branch_node->parent;
  }
  if (not parent_handle) {
    return IsotheticPolygon(Rectangle(Point::atNegInfinity,
                                      Point::atInfinity));
  }
  auto parent_node = treeRef->get_branch_node(parent_handle);
  Branch &b = parent_node->locateBranch(start_handle);
  IsotheticPolygon constraint_poly = b.materialize_polygon(allocator);
  tree_node_handle current_handle = parent_handle;

  while (current_handle) {
    auto current_node = treeRef->get_branch_node(current_handle);
    tree_node_handle parent_handle = current_node->parent;
    if (not parent_handle) {
      return constraint_poly;
    }
    auto parent_node = treeRef->get_branch_node(parent_handle);
    Branch &parent_branch = parent_node->locateBranch(
            current_handle);
    IsotheticPolygon parent_poly = parent_branch.materialize_polygon(
            allocator);

    constraint_poly.intersection(parent_poly);
    constraint_poly.recomputeBoundingBox();
    current_handle = parent_handle;
  }
  return constraint_poly;
#endif

  // Unsupported
  abort();
}

template <int min_branch_factor, int max_branch_factor, class strategy>
void update_branch_polygon(
        Branch &branch_to_update,
        IsotheticPolygon &polygon_to_write,
        NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
        bool force_create = false)
{
#if 0
  if (polygon_to_write.basicRectangles.size() <= MAX_RECTANGLE_COUNT) {
    // Could leak if we had an out of band rectangle before
    branch_to_update.boundingPoly = InlineBoundedIsotheticPolygon();
    std::get<InlineBoundedIsotheticPolygon>(branch_to_update.boundingPoly).push_polygon_to_disk(polygon_to_write);
  } else {
    tree_node_allocator *allocator = get_node_allocator(treeRef);
    if (std::holds_alternative<tree_node_handle>(
            branch_to_update.boundingPoly) and
        not force_create) {
      auto poly_pin =
              allocator->get_tree_node<InlineUnboundedIsotheticPolygon>(
                      std::get<tree_node_handle>(
                              branch_to_update.boundingPoly));
      if (poly_pin->get_max_rectangle_count_on_first_page() ==
          InlineUnboundedIsotheticPolygon::maximum_possible_rectangles_on_first_page()) {
        poly_pin->push_polygon_to_disk(polygon_to_write);
        return;
      }
    }

    uint8_t rect_count = force_create ? polygon_to_write.basicRectangles.size() : InlineUnboundedIsotheticPolygon::maximum_possible_rectangles_on_first_page();

    auto alloc_data =
            allocator->create_new_tree_node<InlineUnboundedIsotheticPolygon>(
                    compute_sizeof_inline_unbounded_polygon(
                            rect_count),
                    NodeHandleType(BIG_POLYGON));
    new (&(*alloc_data.first)) InlineUnboundedIsotheticPolygon(
            allocator,
            rect_count);
    alloc_data.first->push_polygon_to_disk(polygon_to_write);
    branch_to_update.boundingPoly = alloc_data.second;
  }
#endif

  // Unsupported
  abort();
}

template <int min_branch_factor, int max_branch_factor, class strategy>
void BranchNode<min_branch_factor, max_branch_factor, strategy>::reInsert(std::vector<bool> &hasReinsertedOnLevel) {
#if 0
  // Taken from R*
  hasReinsertedOnLevel.at(level_) = true;

  Point globalCenterPoint = boundingBox().centrePoint();

  std::sort(entries.begin(), entries.begin() + cur_offset_,
            [&globalCenterPoint, this](Branch &a, Branch &b) {
                Rectangle rectA = a.materialize_polygon(
                                treeRef->node_allocator_.get())
                        .boundingBox;
                Rectangle rectB = b.materialize_polygon(
                                treeRef->node_allocator_.get())
                        .boundingBox;

                return rectA.centrePoint().distance(
                        globalCenterPoint) > rectB.centrePoint().distance(globalCenterPoint);
            });

  unsigned numNodesToReinsert = 0.3 * cur_offset_;
  unsigned remainder = cur_offset_ - numNodesToReinsert;

  std::vector<Branch> entriesToReinsert;
  entriesToReinsert.reserve(numNodesToReinsert);

  std::copy(
          entries.begin() + remainder,
          entries.begin() + cur_offset_,
          std::back_inserter(entriesToReinsert));

  // We need to perfectly qualify what this branch holds so that other
  // people can fragment around it.
  // This needs to be here because otherwise the mask will be changed
  // by adjusting our bb.
  IsotheticPolygon true_region_mask = get_polygon_path_constraints(
          self_handle_, treeRef);

  cur_offset_ = remainder;

  auto tree_ref_backup = treeRef;
  tree_node_allocator *allocator = get_node_allocator(
          tree_ref_backup);

  fix_up_path_polys(self_handle_, treeRef);
  for (Branch &entry : entriesToReinsert) {
    // FIXME: consolidate with stuff below.
    IsotheticPolygon branch_poly = entry.materialize_polygon(
            allocator);
    branch_poly.intersection(true_region_mask);

    cut_out_branch_region_in_path(self_handle_, branch_poly, treeRef);
  }

  tree_node_handle root_handle = tree_ref_backup->root;

  for (Branch &entry : entriesToReinsert) {

    IsotheticPolygon branch_poly = entry.materialize_polygon(
            allocator);

    branch_poly.intersection(true_region_mask);
    update_branch_polygon(entry, branch_poly, tree_ref_backup);
    auto root_node = tree_ref_backup->get_branch_node(root_handle);

    std::variant<Branch, Point> ent = entry;
    root_handle = root_node->insert(ent, hasReinsertedOnLevel);
  }
#endif

  // Unsupported
  abort();
}

template <class NT, class TR>
std::pair<SplitResult, tree_node_handle> adjust_tree_bottom_half(
        NT current_node,
        TR *tree_ref_backup,
        int max_branch_factor,
        std::vector<bool> &hasReinsertedOnLevel)
{
#if 0
  SplitResult propagationSplit = {
          {InlineBoundedIsotheticPolygon(), tree_node_handle(nullptr)},
          {InlineBoundedIsotheticPolygon(), tree_node_handle(nullptr)}};

  if (current_node->cur_offset_ <= (unsigned)max_branch_factor) {
    return std::make_pair(propagationSplit,
                          tree_node_handle(nullptr));
  }

  // Otherwise, split node
  if (hasReinsertedOnLevel.at(current_node->level_) or true) {
    tree_node_handle parent = current_node->parent;
    auto propagationSplit = current_node->splitNode();

    // Cleanup before ascending
    if (parent != nullptr) {
      // This will probably destroy current_node, so if we need
      // current node for anything, need to do it before the
      // removeEntry call.

      auto parent_node = tree_ref_backup->get_branch_node(parent);
      parent_node->removeBranch(current_node->self_handle_);
    }

    // Ascend, propagating splits
    return std::make_pair(propagationSplit, parent);
  } else {
    // Nothing is real after you make this call
    // The reinsert might have come back around again and split this
    // node, or other nodes, or everyting
    // Signal to the caller that we shoudl stop
    current_node->reInsert(hasReinsertedOnLevel);
    return std::make_pair(propagationSplit,
                          tree_node_handle(nullptr));
  }
#endif

  // Unsupported
  abort();
}

// This bottom-to-top sweep is only for splitting bounding boxes as necessary
template <int min_branch_factor, int max_branch_factor, class strategy>
SplitResult LeafNode<min_branch_factor, max_branch_factor, strategy>::adjustTree(
        std::vector<bool> &hasReinsertedOnLevel)
{
#if 0
  // N.B., as we walk up the tree, we may perform a bunch of splits,
  // which is liable to destroy nodes that are downsplit. These
  // downsplit nodes' memory can then be re-used for other things,
  // like polygons. If we try to use treeRef (or other outside
  // pointers) from that node, it can be clobbered leading to amazing
  // segfaults. It is important that any variables we reference are
  // those we know are alive --- don't just rely on whatever the
  // leaf node class member to have reasonable things after split is
  // called!

  tree_node_handle current_handle = this->self_handle_;

  SplitResult propagationSplit = {
          {InlineBoundedIsotheticPolygon(), tree_node_handle(nullptr)},
          {InlineBoundedIsotheticPolygon(), tree_node_handle(nullptr)}};

  NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *tree_ref_backup =
          this->treeRef;

  // Loop from the bottom to the very top
  while (current_handle != nullptr) {

    // If there was a split we were supposed to propagate then propagate it
    if (propagationSplit.leftBranch.child != nullptr and propagationSplit.rightBranch.child != nullptr) {
      // We are at least one level up, so have to be a branch

      auto current_branch_node = tree_ref_backup->get_branch_node(
              current_handle);
      {
        if (propagationSplit.leftBranch.child.get_type() ==
            LEAF_NODE) {
          auto left_node =
                  tree_ref_backup->get_leaf_node(propagationSplit.leftBranch.child);
          if (left_node->cur_offset_ > 0) {
            current_branch_node->addBranchToNode(
                    propagationSplit.leftBranch);
          }
        } else {
          auto left_node =
                  tree_ref_backup->get_branch_node(propagationSplit.leftBranch.child);
          if (left_node->cur_offset_ > 0) {
            current_branch_node->addBranchToNode(
                    propagationSplit.leftBranch);
          }
        }
      }
      {

        if (propagationSplit.rightBranch.child.get_type() ==
            LEAF_NODE) {
          auto right_node = tree_ref_backup->get_leaf_node(propagationSplit.rightBranch.child);
          if (right_node->cur_offset_ > 0) {

            current_branch_node->addBranchToNode(propagationSplit.rightBranch);
          }
        } else {
          auto right_node = tree_ref_backup->get_branch_node(propagationSplit.rightBranch.child);
          if (right_node->cur_offset_ > 0) {

            current_branch_node->addBranchToNode(propagationSplit.rightBranch);
          }
        }
      }
    }

    std::pair<SplitResult, tree_node_handle>
            split_res_and_new_handle;
    if (current_handle.get_type() == LEAF_NODE || current_handle.get_type() == REPACKED_LEAF_NODE) {
      auto current_leaf_node = tree_ref_backup->get_leaf_node(
              current_handle);
      split_res_and_new_handle = adjust_tree_bottom_half(
              current_leaf_node, tree_ref_backup,
              max_branch_factor, hasReinsertedOnLevel);
    } else {
      auto current_branch_node = tree_ref_backup->get_branch_node(
              current_handle);
      split_res_and_new_handle = adjust_tree_bottom_half(
              current_branch_node, tree_ref_backup,
              max_branch_factor, hasReinsertedOnLevel);
    }

    propagationSplit = split_res_and_new_handle.first;
    current_handle = split_res_and_new_handle.second;
  }
  return propagationSplit;
#endif

    // Unsupported
    abort();
}

template <int min_branch_factor, int max_branch_factor, class strategy>
void fix_up_path_polys(
        tree_node_handle start_handle,
        NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef)
{
  tree_node_handle current_handle = start_handle;
  while (current_handle != nullptr) {
    // Regenerate this node's bounding polygon

    tree_node_handle parent_handle;
    IsotheticPolygon our_poly;
    if (current_handle.get_type() == LEAF_NODE || current_handle.get_type() == REPACKED_LEAF_NODE) {
      auto leaf_node = treeRef->get_leaf_node(current_handle);
      our_poly = IsotheticPolygon(leaf_node->boundingBox());
      parent_handle = leaf_node->parent;
    } else {
      auto branch_node = treeRef->get_branch_node(current_handle);
      our_poly = IsotheticPolygon(branch_node->boundingBox());
      parent_handle = branch_node->parent;
    }
    if (parent_handle) {
      // Q: Is it possible that this is bad?
      // Suppose we just transferred a branch from another region
      // to this spot. Then we might need to expand our polygon,
      // which intersects with other person's polygon who is
      // slightly more cavalier about what regions they think they
      // own. But we now *own* this space. So we need to make sure
      // they don't take it. How can we do that? Fragment their
      // rectangle on the way down.
      auto parent_node = treeRef->get_branch_node(parent_handle);

      // Make this polygon disjoint from its siblings
      parent_node->make_disjoint_from_children(our_poly, current_handle);

      // Now we need to store this poly
      Branch &parent_branch = parent_node->locateBranch(current_handle);
      update_branch_polygon(parent_branch, our_poly, treeRef);
    }
    current_handle = parent_handle;
  }
  // Hit the root, done!
}

template <int min_branch_factor, int max_branch_factor, class strategy>
void cut_out_branch_region_in_path(
        tree_node_handle start_handle,
        IsotheticPolygon &region_to_cut_out,
        NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef)
{
#if 0
  tree_node_handle current_handle = start_handle;
  while (current_handle != nullptr) {
    // Regenerate this node's bounding polygon

    tree_node_handle parent_handle;
    IsotheticPolygon our_poly;
    if (current_handle.get_type() == LEAF_NODE || current_handle.get_type() == REPACKED_LEAF_NODE) {
      auto leaf_node = treeRef->get_leaf_node(current_handle);
      our_poly = IsotheticPolygon(leaf_node->boundingBox());
      parent_handle = leaf_node->parent;
    } else {
      auto branch_node = treeRef->get_branch_node(current_handle);
      our_poly = IsotheticPolygon(branch_node->boundingBox());
      parent_handle = branch_node->parent;
    }
    if (parent_handle) {
      auto parent_node = treeRef->get_branch_node(parent_handle);

      Branch &my_branch = parent_node->locateBranch(current_handle);
      IsotheticPolygon our_poly = my_branch.materialize_polygon(
              treeRef->node_allocator_.get());
      our_poly.increaseResolution(Point::atInfinity, region_to_cut_out);
      our_poly.refine();
      our_poly.recomputeBoundingBox();

      Branch &parent_branch = parent_node->locateBranch(current_handle);
      update_branch_polygon(parent_branch, our_poly, treeRef);
    }
    current_handle = parent_handle;
  }
  // Hit the root, done!
#endif

  // Unsupported
  abort();
}

template <int min_branch_factor, int max_branch_factor, class strategy>
SplitResult adjustTreeSub(
        std::vector<bool> &hasReinsertedOnLevel,
        tree_node_handle start_handle,
        NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef)
{
#if 0
  // N.B., as we walk up the tree, we may perform a bunch of splits,
  // which is liable to destroy nodes that are downsplit. These
  // downsplit nodes' memory can then be re-used for other things,
  // like polygons. If we try to use variables stored in that
  // node, it can be clobbered leading to amazing
  // segfaults. It is important that any variables we reference are
  // those we know are alive.

  tree_node_handle current_handle = start_handle;

  SplitResult propagationSplit = {
          {InlineBoundedIsotheticPolygon(), tree_node_handle(nullptr)},
          {InlineBoundedIsotheticPolygon(), tree_node_handle(nullptr)}};

  // Loop from the bottom to the very top
  while (current_handle != nullptr) {

    // If there was a split we were supposed to propagate then propagate it
    if (propagationSplit.leftBranch.child != nullptr and propagationSplit.rightBranch.child != nullptr) {
      // We are at least one level up, so have to be a branch

      auto current_branch_node = treeRef->get_branch_node(current_handle, true);
      {
        if (propagationSplit.leftBranch.child.get_type() ==
            LEAF_NODE) {
          auto left_node =
                  treeRef->get_leaf_node(propagationSplit.leftBranch.child, false);
          if (left_node->cur_offset_ > 0) {
            current_branch_node->addBranchToNode(
                    propagationSplit.leftBranch);
          }
        } else {
          auto left_node =
                  treeRef->get_branch_node(propagationSplit.leftBranch.child, false);
          if (left_node->cur_offset_ > 0) {
            current_branch_node->addBranchToNode(
                    propagationSplit.leftBranch);
          }
        }
      }
      {

        if (propagationSplit.rightBranch.child.get_type() ==
            LEAF_NODE) {
          auto right_node = treeRef->get_leaf_node(propagationSplit.rightBranch.child, false);
          if (right_node->cur_offset_ > 0) {

            current_branch_node->addBranchToNode(propagationSplit.rightBranch);
          }
        } else {
          auto right_node = treeRef->get_branch_node(propagationSplit.rightBranch.child, false);
          if (right_node->cur_offset_ > 0) {

            current_branch_node->addBranchToNode(propagationSplit.rightBranch);
          }
        }
      }
    }

    std::pair<SplitResult, tree_node_handle>
            split_res_and_new_handle;
    if (current_handle.get_type() == LEAF_NODE || current_handle.get_type() == REPACKED_LEAF_NODE) {
      auto current_leaf_node = treeRef->get_leaf_node(
              current_handle, true);
      split_res_and_new_handle = adjust_tree_bottom_half(
              current_leaf_node, treeRef,
              max_branch_factor, hasReinsertedOnLevel);
    } else {
      auto current_branch_node = treeRef->get_branch_node(
              current_handle);
      split_res_and_new_handle = adjust_tree_bottom_half(
              current_branch_node, treeRef,
              max_branch_factor, hasReinsertedOnLevel);
    }

    propagationSplit = split_res_and_new_handle.first;
    current_handle = split_res_and_new_handle.second;
  }
  return propagationSplit;
#endif

  // Unsupported
  abort();
}

// This always get called on the root node. So if it got called on us,
// that's because we are the only node in the whole tree.
template <int min_branch_factor, int max_branch_factor, class strategy>
tree_node_handle LeafNode<min_branch_factor, max_branch_factor, strategy>::insert(
        Point givenPoint,
        std::vector<bool> &hasReinsertedOnLevel)
{
#if 0
  // This is a leaf, so we are the ONLY node.
  addPoint(givenPoint);

  auto tree_ref_backup = treeRef;

  SplitResult finalSplit = adjustTreeSub(hasReinsertedOnLevel,
                                         self_handle_, treeRef);

  // Grow the tree taller if we need to
  if (finalSplit.leftBranch.child != nullptr and finalSplit.rightBranch.child != nullptr) {
    tree_node_allocator *allocator = get_node_allocator(this->treeRef);
    auto alloc_data =
            allocator->create_new_tree_node<BranchNode<min_branch_factor, max_branch_factor, strategy>>(
                    NodeHandleType(BRANCH_NODE));
    new (&(*alloc_data.first)) BranchNode<min_branch_factor, max_branch_factor, strategy>(this->treeRef,
                                                       tree_node_handle(nullptr), alloc_data.second, level_ + 1);
    auto new_root_handle = alloc_data.second;
    assert(new_root_handle.get_type() == BRANCH_NODE);
    auto new_root_node = alloc_data.first;

    auto left_node = this->treeRef->get_leaf_node(
            finalSplit.leftBranch.child);
    left_node->parent = new_root_handle;
    new_root_node->addBranchToNode(finalSplit.leftBranch);

    auto right_node = this->treeRef->get_leaf_node(
            finalSplit.rightBranch.child);
    right_node->parent = new_root_handle;
    new_root_node->addBranchToNode(finalSplit.rightBranch);

    assert(new_root_handle.get_type() == BRANCH_NODE);
    treeRef->root = new_root_handle;

    assert(this->self_handle_.get_type() == LEAF_NODE);
    allocator->free(this->self_handle_, sizeof(LeafNode<min_branch_factor, max_branch_factor, strategy>));

    // Fix the reinserted length
    hasReinsertedOnLevel.push_back(false);

    return new_root_handle;
  }

  return tree_ref_backup->root;
#endif

    // Unsupported
    abort();
}

// To be called on a leaf
template <int min_branch_factor, int max_branch_factor, class strategy>
void LeafNode<min_branch_factor, max_branch_factor, strategy>::condenseTree() {
#if 0
  auto current_node_handle = this->self_handle_;
  auto previous_node_handle = tree_node_handle(nullptr);

  while (current_node_handle != nullptr) {
    if (current_node_handle.get_type() == LEAF_NODE || current_node_handle.get_type() == REPACKED_LEAF_NODE) {
      auto current_node = treeRef->get_leaf_node(current_node_handle);
      current_node_handle = current_node->parent;
    } else {
      auto current_node = treeRef->get_branch_node(current_node_handle);
      current_node_handle = current_node->parent;
    }
    if (previous_node_handle != nullptr) {
      auto current_parent_node = this->treeRef->get_branch_node(current_node_handle);
      unsigned loc_cur_offset = 0;
      if (previous_node_handle.get_type() == LEAF_NODE || previous_node_handle.get_type() == REPACKED_LEAF_NODE) {
        auto previous_node = treeRef->get_leaf_node(
                previous_node_handle);
        loc_cur_offset = previous_node->cur_offset_;
      } else {
        auto previous_node = treeRef->get_branch_node(
                previous_node_handle);
        loc_cur_offset = previous_node->cur_offset_;
      }
      if (loc_cur_offset == 0) {
        current_parent_node->removeBranch(previous_node_handle);
      }
    }
    previous_node_handle = current_node_handle;
  }
#endif

    // Unsupported
    abort();
}

// Always called on root, this = root
template <int min_branch_factor, int max_branch_factor, class strategy>
tree_node_handle LeafNode<min_branch_factor, max_branch_factor, strategy>::remove(Point givenPoint) {
#if 0
  removePoint(givenPoint);
  return this->self_handle_;
#endif

  // Unsupported
  abort();
}

template <int min_branch_factor, int max_branch_factor, class strategy>
unsigned LeafNode<min_branch_factor, max_branch_factor, strategy>::checksum() {
  unsigned sum = 0;

  for (unsigned i = 0; i < this->cur_offset_; i++) {
    Point &dataPoint = entries.at(i);
    for (unsigned d = 0; d < dimensions; d++) {
      sum += (unsigned)dataPoint[d];
    }
  }
  return sum;
}

template <int min_branch_factor, int max_branch_factor, class strategy>
std::vector<Point> LeafNode<min_branch_factor, max_branch_factor, strategy>::bounding_box_validate() {
  std::vector<Point> my_points;
  for (unsigned i = 0; i < this->cur_offset_; i++) {
    my_points.push_back(entries.at(i));
  }
  return my_points;
}

template <int min_branch_factor, int max_branch_factor, class strategy>
bool LeafNode<min_branch_factor, max_branch_factor, strategy>::validate(tree_node_handle expectedParent, unsigned index) {
#if 0
  if (expectedParent != nullptr and (this->parent != expectedParent ||
                                     this->cur_offset_ > max_branch_factor)) {
    std::cout << "node = " << (void *)this << std::endl;
    std::cout << "parent = " << this->parent << " expectedParent = " << expectedParent << std::endl;
    std::cout << "maxBranchFactor = " << max_branch_factor << std::endl;
    std::cout << "entries.size() = " << this->cur_offset_ << std::endl;
    assert(this->parent == expectedParent);
  }

  if (expectedParent != nullptr) {
    tree_node_allocator *allocator = get_node_allocator(this->treeRef);

    auto parent_node = treeRef->get_branch_node(parent, false);
    Branch &branch = parent_node->locateBranch(this->self_handle_);

    IsotheticPolygon poly;
    if (std::holds_alternative<InlineBoundedIsotheticPolygon>(
            branch.boundingPoly)) {
      poly = std::get<InlineBoundedIsotheticPolygon>(
              branch.boundingPoly)
              .materialize_polygon();
    } else {
      tree_node_handle poly_handle =
              std::get<tree_node_handle>(branch.boundingPoly);
      auto poly_pin =
              InlineUnboundedIsotheticPolygon::read_polygon_from_disk(
                      allocator, poly_handle);
      poly = poly_pin->materialize_polygon();
    }

    for (unsigned i = 0; i < this->cur_offset_; i++) {
      Point &dataPoint = entries.at(i);
      if (!poly.containsPoint(dataPoint)) {
        std::cout << poly << " fails to contain " << dataPoint << std::endl;
        std::cout << "Node: " << self_handle_ << std::endl;
        std::cout << "Parent: " << parent << std::endl;
        assert(false);
      }
    }
  }
  return true;
#endif

  // Unsupported
  abort();
}

template <int min_branch_factor, int max_branch_factor, class strategy>
void LeafNode<min_branch_factor, max_branch_factor, strategy>::print(unsigned n) {
  std::string indentation(n * 4, ' ');
  std::cout << indentation << "Node " << (void *)this << std::endl;
  std::cout << indentation << "    Parent: " << this->parent << std::endl;
  std::cout << indentation << "    Data: ";
  for (unsigned i = 0; i < this->cur_offset_; i++) {
    std::cout << entries.at(i);
  }
  std::cout << std::endl;
}

template <int min_branch_factor, int max_branch_factor, class strategy>
void LeafNode<min_branch_factor, max_branch_factor, strategy>::printTree(unsigned n) {
  // Print this node first
  print(n);
  std::cout << std::endl;
}

template <int min_branch_factor, int max_branch_factor, class strategy>
unsigned LeafNode<min_branch_factor, max_branch_factor, strategy>::height() {
  return 1;
}

template <int min_branch_factor, int max_branch_factor, class strategy>
void BranchNode<min_branch_factor, max_branch_factor, strategy>::deleteSubtrees() {
#if 0
  tree_node_allocator *allocator = get_node_allocator(this->treeRef);
  for (unsigned i = 0; i < this->cur_offset_; i++) {
    Branch &b = entries.at(i);
    tree_node_handle child_handle = b.child;

    if (child_handle.get_type() == REPACKED_LEAF_NODE) {
      child_handle = treeRef->get_leaf_node(child_handle, true);
    } else if (child_handle.get_type() == REPACKED_BRANCH_NODE) {
      child_handle = treeRef->get_branch_node(child_handle, true);
    }

    if (child_handle.get_type() == LEAF_NODE) {
      allocator->free(
              child_handle,
              sizeof(LeafNode<min_branch_factor, max_branch_factor, strategy>));
    } else if (child_handle.get_type() == BRANCH_NODE) {
      auto child = this->treeRef->get_branch_node(child_handle);
      child->deleteSubtrees();
      allocator->free(
              child_handle,
              sizeof(BranchNode<min_branch_factor, max_branch_factor, strategy>));
    }
  }
#endif

  // Unsupported
  abort();
}

template <int min_branch_factor, int max_branch_factor, class strategy>
Rectangle BranchNode<min_branch_factor, max_branch_factor, strategy>::boundingBox() {
  assert(cur_offset_ > 0);
  Rectangle boundingBox = entries.at(0).boundingBox;

  for (unsigned i = 1; i < cur_offset_; i++) {
    boundingBox.expand(entries.at(i).boundingBox);
  }

  return boundingBox;
}

template <int min_branch_factor, int max_branch_factor, class strategy>
void BranchNode<min_branch_factor, max_branch_factor, strategy>::updateBranch(
        tree_node_handle child_handle,
        const InlineBoundedIsotheticPolygon &boundingPoly)
{
  // Locate the child
  unsigned childIndex;
  for (childIndex = 0; entries.at(childIndex).child != child_handle &&
                       childIndex < this->cur_offset_;
       childIndex++) {
  }

  // Update the child
  entries.at(childIndex).boundingPoly = boundingPoly;
}

// Removes a child logical pointer from a this->parent node, freeing that
// child's memory and the memory of the associated polygon (if external
// to the node).
template <int min_branch_factor, int max_branch_factor, class strategy>
void BranchNode<min_branch_factor, max_branch_factor, strategy>::removeBranch(
        const tree_node_handle &entry)
{
#if 0
  unsigned found_index = 0;
  while (found_index < this->cur_offset_) {
    Branch &b = entries.at(found_index);
    if (b.child == entry) {
      break;
    }
    found_index++;
  }
  assert(entries.at(found_index).child == entry);
  Branch &b = entries.at(found_index);
  tree_node_allocator *allocator = get_node_allocator(this->treeRef);

  if (std::holds_alternative<tree_node_handle>(b.boundingPoly)) {
    tree_node_handle free_poly_handle = std::get<tree_node_handle>(
            b.boundingPoly);
    auto poly_pin =
            InlineUnboundedIsotheticPolygon::read_polygon_from_disk(
                    allocator, free_poly_handle);
    poly_pin->free_subpages(allocator);
    uint16_t alloc_size = compute_sizeof_inline_unbounded_polygon(
            poly_pin->get_max_rectangle_count_on_first_page());
    allocator->free(free_poly_handle, alloc_size);
  }

  if (b.child.get_type() == LEAF_NODE) {
    allocator->free(b.child, sizeof(LeafNode<min_branch_factor, max_branch_factor, strategy>));
  } else if (b.child.get_type() == BRANCH_NODE) {
    allocator->free(b.child, sizeof(BranchNode<min_branch_factor, max_branch_factor, strategy>));
  }
  b.child = tree_node_handle(nullptr);

  entries.at(found_index) = entries.at(this->cur_offset_ - 1);
  this->cur_offset_--;
#endif

  // Unsupported
  abort();
}

template <int min_branch_factor, int max_branch_factor, class strategy, typename functor>
void treeWalker(NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef, tree_node_handle root, functor &f) {
  std::stack<tree_node_handle> context;
  context.push(root);
  tree_node_handle currentContext;

  while (!context.empty()) {
    currentContext = context.top();
    context.pop();

    f(treeRef, currentContext);

    if (currentContext.get_type() == BRANCH_NODE) {
      auto node = treeRef->get_branch_node(currentContext);
      for (unsigned i = 0; i < node->cur_offset_; i++) {
        context.push(node->entries.at(i).child);
      }
    }
  }
}

template <int min_branch_factor, int max_branch_factor, class strategy>
void point_search_leaf_node(LeafNode<min_branch_factor, max_branch_factor, strategy> &node,
                            Point &requestedPoint,
                            std::vector<Point> &accumulator,
                            NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef)
{
  unsigned intersection_count = 0;

  for (unsigned i = 0; i < node.cur_offset_; i++) {
    const Point &p = node.entries.at(i);
    intersection_count++;
    if (p == requestedPoint) {
      accumulator.push_back(p);
      break;
    }
  }

  treeRef->stats.recordIntersectionCount(intersection_count);
}

template <int min_branch_factor, int max_branch_factor, class strategy>
void point_search_branch_node(BranchNode<min_branch_factor, max_branch_factor, strategy> &node,
                              Point &requestedPoint,
                              std::stack<tree_node_handle> &context,
                              NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef)
{
  unsigned matching_branch_counter = 0;
  unsigned intersection_count = 0;

  for (size_t i = 0; i < node.cur_offset_; i++) {
    Branch &b = node.entries.at(i);

    intersection_count++;
    if (b.boundingBox.containsPoint(requestedPoint)) {
      auto itr = treeRef->polygons.find(b.child);

      // This branch has no polygons, just use info from native MBR
      if (itr == treeRef->polygons.end()) {
        context.push(b.child);
        matching_branch_counter++;
        break;
      }

      IsotheticPolygon polygon = itr->second;

      if (polygon.containsPointWithMetrics(requestedPoint, intersection_count)) {
        context.push(b.child);
        matching_branch_counter++;
        break;
      }
    }
  }

  treeRef->stats.recordScatter(matching_branch_counter);
  treeRef->stats.recordIntersectionCount(intersection_count);
}


template <int min_branch_factor, int max_branch_factor, class strategy>
std::vector<Point> point_search(tree_node_handle start_point, Point &requestedPoint,
        NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef) {
  std::vector<Point> accumulator;
  std::stack<tree_node_handle> context;
  context.push(start_point);

  while (not context.empty()) {
    tree_node_handle current_handle = context.top();
    context.pop();

    if (current_handle.get_type() == LEAF_NODE) {
      auto current_node = treeRef->get_leaf_node(current_handle);
      point_search_leaf_node(*current_node, requestedPoint, accumulator, treeRef);
#ifdef STAT
      treeRef->stats.markLeafSearched();
#endif
    } else if (current_handle.get_type() == BRANCH_NODE) {
      auto current_node = treeRef->get_branch_node(current_handle);
      point_search_branch_node(*current_node, requestedPoint, context, treeRef);
#ifdef STAT
      treeRef->stats.markNonLeafNodeSearched();
#endif
    } else {
      assert(false);
    }
  }
#ifdef STAT
  treeRef->stats.resetSearchTracker(false);
#endif
  return accumulator;
}

template <int min_branch_factor, int max_branch_factor, class strategy>
tree_node_handle parent_handle_point_search(
        tree_node_handle start_point,
        Point &requestedPoint,
        NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
        tree_node_handle child_to_stop_at) {
#if 0
  std::stack<tree_node_handle> context;
  context.push(start_point);
  tree_node_allocator *allocator = treeRef->node_allocator_.get();

  while (not context.empty()) {
    tree_node_handle current_handle = context.top();
    context.pop();
    if (current_handle.get_type() == LEAF_NODE) {
      return tree_node_handle(nullptr);
    } else if (current_handle.get_type() == REPACKED_LEAF_NODE) {
      return tree_node_handle(nullptr);
    } else if (current_handle.get_type() == BRANCH_NODE) {
      parent_search_branch_handle(current_handle,
                                  requestedPoint, context, child_to_stop_at);
    } else if (current_handle.get_type() == REPACKED_BRANCH_NODE) {
      parent_search_packed_branch_handle(current_handle,
                                         requestedPoint, context, child_to_stop_at);
    } else {
      assert(false);
    }
  }

  return tree_node_handle(nullptr);
#endif

  // Unsupported
  abort();
}

template <int min_branch_factor, int max_branch_factor, class strategy>
void rectangle_search_leaf_node(LeafNode<min_branch_factor, max_branch_factor, strategy> &node,
                                Rectangle &requestedRectangle,
                                std::vector<Point> &accumulator,
                                NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef)
{
  unsigned intersection_count = 0;

  for (unsigned i = 0; i < node.cur_offset_; i++) {
    const Point &p = node.entries.at(i);
    intersection_count++;
    if (requestedRectangle.containsPoint(p)) {
      accumulator.push_back(p);
    }
  }

  treeRef->stats.recordIntersectionCount(intersection_count);
}

template <int min_branch_factor, int max_branch_factor, class strategy>
void rectangle_search_branch_node(BranchNode<min_branch_factor, max_branch_factor, strategy> &node,
                                  Rectangle &requestedRectangle,
                                  std::stack<tree_node_handle> &context,
                                  NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef)
{
  unsigned intersection_count = 0;

  for (size_t i = 0; i < node.cur_offset_; i++) {
    Branch &b = node.entries.at(i);

    intersection_count++;
    if (b.boundingBox.intersectsRectangle(requestedRectangle)) {
      auto itr = treeRef->polygons.find(b.child);

      // This branch has no polygons, just use info from native MBR
      if (itr == treeRef->polygons.end()) {
        context.push(b.child);
        continue;
      }

      IsotheticPolygon polygon = itr->second;

      if (polygon.intersectsRectangleWithMetrics(requestedRectangle, intersection_count)) {
        context.push(b.child);
      }
    }
  }

  treeRef->stats.recordIntersectionCount(intersection_count);
}

template <int min_branch_factor, int max_branch_factor, class strategy>
std::vector<Point> rectangle_search(
        tree_node_handle start_point,
        Rectangle &requestedRectangle,
        NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
        bool should_track_search = true)
{
  std::vector<Point> accumulator;

  std::stack<tree_node_handle> context;
  tree_node_allocator *allocator = treeRef->node_allocator_.get();
  context.push(start_point);

  while (not context.empty()) {
    tree_node_handle current_handle = context.top();
    context.pop();

    if (current_handle.get_type() == LEAF_NODE) {
      auto current_node = treeRef->get_leaf_node(current_handle);
      rectangle_search_leaf_node(*current_node, requestedRectangle, accumulator, treeRef);
#ifdef STAT
      if (should_track_search) {
        treeRef->stats.markLeafSearched();
      }
#endif
    } else if (current_handle.get_type() == BRANCH_NODE) {
      auto current_node = treeRef->get_branch_node(current_handle);
      rectangle_search_branch_node(*current_node, requestedRectangle, context, treeRef);
#ifdef STAT
      if (should_track_search) {
        treeRef->stats.markNonLeafNodeSearched();
      }
#endif
    } else {
      assert(false);
    }
  }
#ifdef STAT
  treeRef->stats.resetSearchTracker(true);
#endif
  return accumulator;
}

// Always called on root, this = root
// This top-to-bottom sweep is only for adjusting bounding boxes to contain the point and
// choosing a particular leaf
template <int min_branch_factor, int max_branch_factor, class strategy>
tree_node_handle
BranchNode<min_branch_factor, max_branch_factor, strategy>::chooseNode(
        std::variant<Branch, Point> &nodeEntry,
        uint8_t stopping_level) {
#if 0
  // FIXME: try and avoid all these materialize calls

  // CL1 [Initialize]
  tree_node_handle cur_node_handle = this->self_handle_;

  assert(cur_node_handle != nullptr);

  for (;;) {
    assert(cur_node_handle != nullptr);
    if (cur_node_handle.get_type() == REPACKED_LEAF_NODE) {
      cur_node_handle = treeRef->get_leaf_node(cur_node_handle, true)->self_handle_;
    } else if (cur_node_handle.get_type() == REPACKED_BRANCH_NODE) {
      cur_node_handle = treeRef->get_branch_node(cur_node_handle, true)->self_handle_;
    }

    if (cur_node_handle.get_type() == LEAF_NODE) {
      assert(std::holds_alternative<Point>(nodeEntry));
      return cur_node_handle;
    } else {
      assert(cur_node_handle.get_type() == BRANCH_NODE);

      auto cur_node = treeRef->get_branch_node(cur_node_handle);
      // Compute the smallest expansion
      assert(cur_node->cur_offset_ > 0);
      if (cur_node->level_ == stopping_level) {
        return cur_node_handle;
      }

      tree_node_allocator *allocator = get_node_allocator(this->treeRef);
      inline_poly node_poly = cur_node->entries.at(0).get_inline_polygon(
              allocator);
      bool node_poly_unbounded = std::holds_alternative<unbounded_poly_pin>(node_poly);
      if (node_poly_unbounded) {
        assert(std::get<unbounded_poly_pin>(node_poly)->get_total_rectangle_count() > 0);
      } else {
        assert(std::get<InlineBoundedIsotheticPolygon>(node_poly).get_rectangle_count() > 0);
      }

      // This is the minimum amount of additional area we need in
      // one of the branches to encode our expansion
      double minimal_area_expansion =
              std::numeric_limits<double>::max();
      double minimal_poly_area = std::numeric_limits<double>::max();

      // This is the branch that gives us that minimum area
      // expansion
      unsigned smallestExpansionBranchIndex = 0;

      // This is the list of optimal expansiosn we need to perform
      // to get the bounding box/bounding polygon
      std::vector<IsotheticPolygon::OptimalExpansion> expansions;

      if (std::holds_alternative<Branch>(nodeEntry)) {
        inline_poly branch_poly =
                std::get<Branch>(nodeEntry).get_inline_polygon(allocator);
        auto expansion_computation_results = computeExpansionArea(node_poly, branch_poly);
        minimal_area_expansion = expansion_computation_results.first;

        bool branch_poly_unbounded = std::holds_alternative<unbounded_poly_pin>(branch_poly);
        if (branch_poly_unbounded) {
          minimal_poly_area = std::get<unbounded_poly_pin>(branch_poly)->area();
        } else {
          minimal_poly_area = std::get<InlineBoundedIsotheticPolygon>(branch_poly).area();
        }

        expansions = expansion_computation_results.second;
        /*
                if( minimal_area_expansion == -1.0 ) {
                    abort();
                }
                */
      } else {
        IsotheticPolygon::OptimalExpansion exp;
        if (node_poly_unbounded) {
          auto unb_node_poly = std::get<unbounded_poly_pin>(node_poly);
          auto inline_exp = computeExpansionArea<InlineUnboundedIsotheticPolygon, InlineUnboundedIsotheticPolygon::Iterator>(
                  *unb_node_poly, unb_node_poly->begin(), unb_node_poly->end(), std::get<Point>(nodeEntry));
          exp = inline_exp;
        } else {
          auto b_node_poly = std::get<InlineBoundedIsotheticPolygon>(node_poly);
          exp = computeExpansionArea<InlineBoundedIsotheticPolygon, InlineBoundedIsotheticPolygon::iterator>(
                  b_node_poly, b_node_poly.begin(), b_node_poly.end(), std::get<Point>(nodeEntry));
        }
        minimal_area_expansion = exp.area;
        expansions.push_back(exp);
      }

      for (unsigned i = 1; i < cur_node->cur_offset_; i++) {
        Branch &b = cur_node->entries.at(i);
        node_poly = b.get_inline_polygon(allocator);
        bool node_poly_unbounded = std::holds_alternative<unbounded_poly_pin>(node_poly);
        if (node_poly_unbounded) {
          assert(std::get<unbounded_poly_pin>(node_poly)->get_total_rectangle_count() > 0);
        } else {
          assert(std::get<InlineBoundedIsotheticPolygon>(node_poly).get_rectangle_count() > 0);
        }

        if (std::holds_alternative<Branch>(nodeEntry)) {
          inline_poly branch_poly =
                  std::get<Branch>(nodeEntry).get_inline_polygon(
                          allocator);
          // Walk every rectangle in the branch's polygon
          // Find rectangle in our polygon that needs to be
          // expanded the least to fit the branch's rectangle
          // inside it.
          // N.B., this does not split the rectangle apart if
          // the expanded rectangle could be part of two
          // distinct polygons. So as a result of doing this
          // the polygon's constituent rectangles may now
          // overlap.
          auto expansion_computation_results = computeExpansionArea(node_poly, branch_poly);
          double poly_area;

          if (node_poly_unbounded) {
            poly_area = std::get<unbounded_poly_pin>(node_poly)->area();
          } else {
            poly_area = std::get<InlineBoundedIsotheticPolygon>(node_poly).area();
          }

          if (expansion_computation_results.first <
              minimal_area_expansion or
              (expansion_computation_results.first ==
               minimal_area_expansion and
               poly_area <
               minimal_poly_area)) {
            minimal_area_expansion =
                    expansion_computation_results.first;
            minimal_poly_area = poly_area;
            expansions = expansion_computation_results.second;
            smallestExpansionBranchIndex = i;
          }
          /*
                    if( minimal_area_expansion == -1.0 ) {
                        abort();
                    }
                    */
        } else {
          IsotheticPolygon::OptimalExpansion exp;
          if (node_poly_unbounded) {
            auto unb_node_poly = std::get<unbounded_poly_pin>(node_poly);
            auto inline_exp = computeExpansionArea<InlineUnboundedIsotheticPolygon, InlineUnboundedIsotheticPolygon::Iterator>(
                    *unb_node_poly, unb_node_poly->begin(), unb_node_poly->end(), std::get<Point>(nodeEntry));
            exp = inline_exp;
          } else {
            auto b_node_poly = std::get<InlineBoundedIsotheticPolygon>(node_poly);
            auto inline_exp = computeExpansionArea<InlineBoundedIsotheticPolygon, InlineBoundedIsotheticPolygon::iterator>(
                    b_node_poly, b_node_poly.begin(), b_node_poly.end(), std::get<Point>(nodeEntry));
            exp = inline_exp;
          }

          if (exp.area < minimal_area_expansion) {
            minimal_area_expansion = exp.area;
            expansions.clear();
            expansions.push_back(exp);
            smallestExpansionBranchIndex = i;
          }
        }
      }

      if (minimal_area_expansion != -1.0) {
#ifndef NDEBUG
        for (unsigned i = 0; i < cur_node->cur_offset_; i++) {
          for (unsigned j = 0; j < cur_node->cur_offset_; j++) {
            if (i == j) {
              continue;
            }
            Branch &b_i = cur_node->entries.at(i);
            Branch &b_j = cur_node->entries.at(j);
            IsotheticPolygon poly_i =
                    b_i.materialize_polygon(allocator);
            IsotheticPolygon poly_j =
                    b_j.materialize_polygon(allocator);
            assert(poly_i.disjoint(poly_j));
          }
        }
#endif

        Branch &chosen_branch =
                cur_node->entries.at(smallestExpansionBranchIndex);
        node_poly = chosen_branch.get_inline_polygon(allocator);

        IsotheticPolygon mat_node_poly = chosen_branch.materialize_polygon(allocator);

        if (std::holds_alternative<Branch>(nodeEntry)) {
          // Fragment them on the way down.
          Branch &inserting_branch = std::get<Branch>(nodeEntry);
          IsotheticPolygon insertion_poly = inserting_branch.materialize_polygon(allocator);
          for (unsigned i = 0; i < cur_node->cur_offset_; i++) {
            if (i == smallestExpansionBranchIndex) {
              continue;
            }
            Branch &other_branch = cur_node->entries.at(i);
            IsotheticPolygon other_poly =
                    other_branch.materialize_polygon(allocator);

            other_poly.increaseResolution(Point::atInfinity, insertion_poly);
            other_poly.refine();
            other_poly.recomputeBoundingBox();

            update_branch_polygon(other_branch, other_poly, treeRef);
          }
        }

        // We need to expand on way down so we know who is
        // responsible for the new point/branch
        // Everyone else needs to fragment around my nodeEntry,
        // then we expand and fragment around them.
        if (std::holds_alternative<Branch>(nodeEntry)) {
          Branch &inserting_branch = std::get<Branch>(nodeEntry);
          IsotheticPolygon insertion_poly = inserting_branch.materialize_polygon(allocator);
          assert(insertion_poly.basicRectangles.size() ==
                 expansions.size());
          for (unsigned i = 0; i <
                               insertion_poly.basicRectangles.size();
               i++) {
            auto &expansion = expansions.at(i);
            auto &insertion_rect = insertion_poly.basicRectangles.at(i);
            Rectangle &existing_rect =
                    mat_node_poly.basicRectangles.at(expansion.index);
            // Expand the existing rectangle. This rectangle
            // might now overlap with other rectangles in
            // the polygon. But if we make it not overlap,
            // then we alter the indices of the expansion
            // rectangles, which kind of sucks, So, leave it
            // for now.
            existing_rect.expand(insertion_rect);
            assert(existing_rect.containsRectangle(
                    insertion_rect));
          }
          mat_node_poly.recomputeBoundingBox();
        } else {
          assert(expansions.size() == 1);
          Point &p = std::get<Point>(nodeEntry);
          Rectangle &existing_rect =
                  mat_node_poly.basicRectangles.at(expansions.at(0).index);
          existing_rect.expand(p);
          mat_node_poly.recomputeBoundingBox();
          assert(mat_node_poly.containsPoint(p));
        }

        // Dodge all the other branches
        for (unsigned i = 0; i < cur_node->cur_offset_;
             i++) {
          if (i == smallestExpansionBranchIndex) {
            continue;
          }
          Branch &other_branch = cur_node->entries.at(i);
          mat_node_poly.increaseResolution(Point::atInfinity,
                                           other_branch.materialize_polygon(allocator));
        }

        /*
                if( std::holds_alternative<Branch>( nodeEntry ) ) {
                    Branch &inserting_branch = std::get<Branch>( nodeEntry );
                    IsotheticPolygon insertion_poly = inserting_branch.materialize_polygon( allocator );
                    for( unsigned i = 0; i < insertion_poly.basicRectangles.size(); i++ ) {
                        bool contained = false;
                        for( unsigned j = 0; j < node_poly.basicRectangles.size(); j++ ) {
                            if(
                                    node_poly.basicRectangles.at(j).containsRectangle(
                                        insertion_poly.basicRectangles.at(i)
                                        ) ) {
                                contained = true;
                                break;
                            }
                        }
                        assert( contained );
                    }

                }
                */

        mat_node_poly.refine();
        mat_node_poly.recomputeBoundingBox();

        update_branch_polygon(chosen_branch, mat_node_poly, treeRef);
      }

      // Descend
      Branch &b = cur_node->entries.at(smallestExpansionBranchIndex);
      cur_node_handle = b.child;
      assert(cur_node_handle != nullptr);
    }
  }
#endif

  // Unsupported
  abort();
}

template <int min_branch_factor, int max_branch_factor, class strategy>
tree_node_handle BranchNode<min_branch_factor, max_branch_factor, strategy>::findLeaf(Point givenPoint) {
#if 0
  // Initialize our context stack
  std::stack<tree_node_handle> context;
  context.push(this->self_handle_);
  tree_node_handle current_node_handle;

  tree_node_allocator *allocator = get_node_allocator(this->treeRef);

  while (!context.empty()) {
    current_node_handle = context.top();
    context.pop();

    if (current_node_handle.get_type() == REPACKED_LEAF_NODE) {
      current_node_handle = treeRef->get_leaf_node(current_node_handle, true)->self_handle_;
    } else if (current_node_handle.get_type() == REPACKED_BRANCH_NODE) {
      current_node_handle = treeRef->get_branch_node(current_node_handle, true)->self_handle_;
    }

    if (current_node_handle.get_type() == LEAF_NODE) {
      auto current_node = treeRef->get_leaf_node(
              current_node_handle);
      // FL2 [Search leaf node for record]
      // Check each entry to see if it matches E
      for (unsigned i = 0; i < current_node->cur_offset_; i++) {
        Point &p = current_node->entries.at(i);
        if (p == givenPoint) {
          return current_node_handle;
        }
      }
    } else {
      auto current_node = treeRef->get_branch_node(current_node_handle);
      // FL1 [Search subtrees]
      // Determine which branches we need to follow
      for (unsigned i = 0; i < current_node->cur_offset_; i++) {
        Branch &b = current_node->entries.at(i);
        IsotheticPolygon poly;
        if (std::holds_alternative<InlineBoundedIsotheticPolygon>(
                b.boundingPoly)) {
          InlineBoundedIsotheticPolygon &loc_poly =
                  std::get<InlineBoundedIsotheticPolygon>(
                          b.boundingPoly);
          // Quick check
          if (!loc_poly.get_summary_rectangle().containsPoint(
                  givenPoint)) {
            continue;
          }
          // Full containment check required
          poly = loc_poly.materialize_polygon();
        } else {
          tree_node_handle poly_handle =
                  std::get<tree_node_handle>(b.boundingPoly);
          auto poly_pin =
                  InlineUnboundedIsotheticPolygon::read_polygon_from_disk(
                          allocator, poly_handle);
          // Quick check
          if (!poly_pin->get_summary_rectangle().containsPoint(
                  givenPoint)) {
            continue;
          }
          // Full containment check required
          poly = poly_pin->materialize_polygon();
        }
        if (poly.containsPoint(givenPoint)) {
          // Add the child to the nodes we will consider
          context.push(b.child);
        }
      }
    }
  }

  return tree_node_handle(nullptr);
#endif

  // Unsupported
  abort();
}

template <int min_branch_factor, int max_branch_factor, class strategy>
Partition BranchNode<min_branch_factor, max_branch_factor, strategy>::partitionNode() {
  return partitionBranchNode();
}

struct summary_rectangle_sorter {

    enum class sort_point {
        LOWER_LEFT,
        UPPER_RIGHT
    };

    summary_rectangle_sorter(unsigned dimension, sort_point sort_on, tree_node_allocator *allocator) : dimension_(dimension),
                                                                                                       allocator_(allocator),
                                                                                                       sort_on_(sort_on) {
    }

    template <class NE>
    bool operator()(NE &n1, NE &n2) {
#if 0
      Branch &b1 = std::get<Branch>(n1);
      Branch &b2 = std::get<Branch>(n2);

      Rectangle bb1 = b1.get_summary_rectangle(
              allocator_);
      Rectangle bb2 = b2.get_summary_rectangle(
              allocator_);
      if (sort_on_ == sort_point::LOWER_LEFT) {
        return bb1.lowerLeft[dimension_] <= bb2.lowerLeft[dimension_];
      }
      return bb1.upperRight[dimension_] <= bb2.upperRight[dimension_];
#endif

      // Unsupported
      abort();
    }

    unsigned dimension_;
    tree_node_allocator *allocator_;
    sort_point sort_on_;
};

template <int min_branch_factor, int max_branch_factor, class strategy>
void BranchNode<min_branch_factor, max_branch_factor, strategy>::make_disjoint_from_children(
        IsotheticPolygon &polygon,
        tree_node_handle handle_to_skip) {
#if 0
  tree_node_allocator *allocator = get_node_allocator(this->treeRef);
  for (auto iter = entries.begin(); iter !=
                                    entries.begin() + this->cur_offset_;
       iter++) {
    Branch &b = *iter;
    if (b.child == handle_to_skip) {
      continue;
    }
    IsotheticPolygon child_poly = b.materialize_polygon(allocator);
    polygon.increaseResolution(Point::atInfinity, child_poly);
  }
  polygon.recomputeBoundingBox();
#endif

  // Unsupported
  abort();
}

// We create two new nodes and free the old one.
// The old one is freed in adjustTree using removeEntry
// If we downsplit, then we won't call adjustTree for that split so we
// need to delete the node ourselves.
template <int min_branch_factor, int max_branch_factor, class strategy>
SplitResult BranchNode<min_branch_factor, max_branch_factor, strategy>::splitNode(Partition p, bool is_downsplit) {
#if 0
  assert(this->self_handle_.get_type() == BRANCH_NODE);
  using NodeType = BranchNode<min_branch_factor, max_branch_factor, strategy>;
  tree_node_allocator *allocator = get_node_allocator(this->treeRef);

  auto alloc_data =
          allocator->create_new_tree_node<BranchNode<min_branch_factor, max_branch_factor, strategy>>(
                  NodeHandleType(BRANCH_NODE));
  tree_node_handle left_handle = alloc_data.second;
  auto left_node = alloc_data.first; // take pin
  new (&(*left_node)) NodeType(this->treeRef, this->parent,
                               left_handle, level_);
  assert(left_node->self_handle_ == left_handle);

  alloc_data =
          allocator->create_new_tree_node<BranchNode<min_branch_factor, max_branch_factor, strategy>>(
                  NodeHandleType(BRANCH_NODE));
  tree_node_handle right_handle = alloc_data.second;
  auto right_node = alloc_data.first; // take pin
  new (&(*right_node)) NodeType(this->treeRef, this->parent,
                                right_handle, level_);
  assert(right_node->self_handle_ == right_handle);

  SplitResult split = {
          {tree_node_handle(nullptr), left_handle},
          {tree_node_handle(nullptr), right_handle}};

  // So we are going to split this branch node.
  for (unsigned i = 0; i < this->cur_offset_; i++) {
    Branch &branch = entries.at(i);

    pinned_node_ptr<InlineUnboundedIsotheticPolygon>
            branch_poly_pin(allocator->buffer_pool_, nullptr,
                            nullptr);

    Rectangle summary_rectangle =
            branch.get_summary_rectangle(allocator);

    bool is_contained_left =
            summary_rectangle.upperRight[p.dimension] <= p.location;
    bool is_contained_right =
            summary_rectangle.lowerLeft[p.dimension] >= p.location;
    assert(not(is_contained_left and is_contained_right));

    IsotheticPolygon child_poly = branch.materialize_polygon(
            allocator);

    if (branch.child.get_type() == REPACKED_LEAF_NODE) {
      branch.child = treeRef->get_leaf_node(branch.child, true)->self_handle_;
    } else if (branch.child.get_type() == REPACKED_BRANCH_NODE) {
      branch.child = treeRef->get_branch_node(branch.child, true)->self_handle_;
    }

    // Entirely contained in the left polygon
    if (is_contained_left and not is_contained_right) {
      if (branch.child.get_type() == LEAF_NODE) {
        auto child = this->treeRef->get_leaf_node(branch.child);
        child->parent = split.leftBranch.child;
      } else {
        auto child = this->treeRef->get_branch_node(branch.child);
        child->parent = split.leftBranch.child;
      }

      assert(split.leftBranch.child == left_handle);
      left_node->addBranchToNode(branch);
      // Entirely contained in the right polygon
    } else if (is_contained_right and not is_contained_left) {
      if (branch.child.get_type() == LEAF_NODE) {
        auto child = this->treeRef->get_leaf_node(branch.child);
        child->parent = split.rightBranch.child;
      } else {
        auto child = this->treeRef->get_branch_node(branch.child);
        child->parent = split.rightBranch.child;
      }
      assert(split.rightBranch.child == right_handle);
      right_node->addBranchToNode(branch);
    } else if (summary_rectangle.upperRight[p.dimension] ==
               summary_rectangle.lowerLeft[p.dimension] and
               summary_rectangle.lowerLeft[p.dimension] ==
               p.location) {
      // These go left or right situationally
      if (left_node->cur_offset_ <= right_node->cur_offset_) {
        if (branch.child.get_type() == LEAF_NODE) {
          auto child = this->treeRef->get_leaf_node(branch.child);
          child->parent = split.leftBranch.child;
        } else {
          auto child = this->treeRef->get_branch_node(branch.child);
          child->parent = split.leftBranch.child;
        }
        assert(split.leftBranch.child == left_handle);
        left_node->addBranchToNode(branch);
      } else {
        if (branch.child.get_type() == LEAF_NODE) {
          auto child = this->treeRef->get_leaf_node(branch.child);
          child->parent = split.rightBranch.child;
        } else {
          auto child = this->treeRef->get_branch_node(branch.child);
          child->parent = split.rightBranch.child;
        }
        assert(split.rightBranch.child == right_handle);
        right_node->addBranchToNode(branch);
      }
      // Partially spanned by both nodes, need to downsplit
    } else {
      IsotheticPolygon branch_poly = branch.materialize_polygon(
              allocator);
      SplitResult downwardSplit;

      if (branch.child.get_type() == LEAF_NODE) {
        auto child = treeRef->get_leaf_node(branch.child);
        downwardSplit = child->splitNode(p, true);
        allocator->free(branch.child, sizeof(
                LeafNode<min_branch_factor, max_branch_factor, strategy>));
      } else {
        auto child = treeRef->get_branch_node(branch.child);
        downwardSplit = child->splitNode(p, true);
        allocator->free(branch.child, sizeof(
                BranchNode<min_branch_factor, max_branch_factor, strategy>));
      }

      if (std::holds_alternative<tree_node_handle>(
              branch.boundingPoly)) {
        // We can free this
        tree_node_handle free_poly_handle = std::get<tree_node_handle>(
                branch.boundingPoly);
        auto free_poly_pin = InlineUnboundedIsotheticPolygon::read_polygon_from_disk(
                allocator, free_poly_handle);
        free_poly_pin->free_subpages(allocator);
        uint16_t alloc_size = compute_sizeof_inline_unbounded_polygon(
                free_poly_pin->get_max_rectangle_count_on_first_page());

        allocator->free(free_poly_handle, alloc_size);
      }

      if (downwardSplit.leftBranch.child.get_type() == LEAF_NODE) {
        auto left_child =
                treeRef->get_leaf_node(downwardSplit.leftBranch.child);
        if (left_child->cur_offset_ > 0) {
          left_child->parent = split.leftBranch.child;
          left_node->addBranchToNode(downwardSplit.leftBranch);
        }
      } else {
        auto left_child =
                treeRef->get_branch_node(downwardSplit.leftBranch.child);
        if (left_child->cur_offset_ > 0) {
          left_child->parent = split.leftBranch.child;
          left_node->addBranchToNode(downwardSplit.leftBranch);
        }
      }

      if (downwardSplit.rightBranch.child.get_type() == LEAF_NODE) {
        auto right_child =
                treeRef->get_leaf_node(downwardSplit.rightBranch.child);
        if (right_child->cur_offset_ > 0) {
          right_child->parent = split.rightBranch.child;
          right_node->addBranchToNode(downwardSplit.rightBranch);
        }

      } else {
        auto right_child =
                treeRef->get_branch_node(downwardSplit.rightBranch.child);
        if (right_child->cur_offset_ > 0) {
          right_child->parent = split.rightBranch.child;
          right_node->addBranchToNode(downwardSplit.rightBranch);
        }
      }
    } //downsplit
  }   //split

  // It is possible that after splitting on the geometric median,
  // we still end up with an overfull node. This can happen
  // everything gets assigned to the left node except for one
  // branch that needs a downward split. That downward split
  // results in a node added to the left and to the right,
  // resulting in an overfull left node.
  // We have a heuristic that generally solves this problem, but
  // provably does not work in all cases. Fix in the future, alert
  // us if it happens
  // TODO: FIXME

  assert(left_node->cur_offset_ <= max_branch_factor and
         right_node->cur_offset_ <= max_branch_factor);

  IsotheticPolygon left_polygon(left_node->boundingBox());
  IsotheticPolygon right_polygon(right_node->boundingBox());

  assert(left_polygon.disjoint(right_polygon));

  // When downsplitting our node, one part of this node goes
  // to the "left parent", and one part of the node goes to
  // the "right parent". These node parts could revise their
  // bounding boxes to reflect only what is left in their nodes,
  // because part of their contents have been removed. These
  // revisions may include "simplifying" the bounding box, in
  // which case the space the polygon consumes becomes larger,
  // but it contains fewer rectangles.
  // E.g., space used be to be an "L" shape, but is now just
  // a rectangle.
  //
  // However, our siblings may also be revising their bounding
  // boxes, and they will do so based on what our bounding
  // box looks like. They need to ensure that the boxes do not
  // intersect. They could theoretically look at the
  // new parent node, figure out what all of their sibling's
  // boxes are, and then make themselves disjoint, but that's
  // a lot more work than just intersectingw ith our existing
  // polygon, which is already guaranteed to be disjoint from
  // our siblings. So, we do the latter.
  if (this->parent) {
    auto parent_node = treeRef->get_branch_node(parent);
    if (not is_downsplit) {
      parent_node->make_disjoint_from_children(left_polygon,
                                               this->self_handle_);
      assert(left_polygon.basicRectangles.size() > 0);
      left_polygon.refine();
      assert(left_polygon.basicRectangles.size() > 0);
      parent_node->make_disjoint_from_children(right_polygon,
                                               this->self_handle_);
      assert(right_polygon.basicRectangles.size() > 0);
      right_polygon.refine();
      assert(right_polygon.basicRectangles.size() > 0);
    } else {
      // Intersect with our existing poly to avoid intersect
      // with other children
      Branch &b = parent_node->locateBranch(this->self_handle_);
      IsotheticPolygon parent_poly = b.materialize_polygon(
              allocator);
      left_polygon.intersection(parent_poly);
      assert(left_polygon.basicRectangles.size() > 0);
      left_polygon.refine();
      assert(left_polygon.basicRectangles.size() > 0);
      right_polygon.intersection(parent_poly);
      assert(right_polygon.basicRectangles.size() > 0);
      right_polygon.refine();
      assert(right_polygon.basicRectangles.size() > 0);
    }
  }

  // Writeback our polygons
  if (left_polygon.basicRectangles.size() > MAX_RECTANGLE_COUNT) {
    unsigned overfull_rect_count = (unsigned)(2 *
                                              left_polygon.basicRectangles.size());

    overfull_rect_count =
            std::min(overfull_rect_count,
                     (unsigned)InlineUnboundedIsotheticPolygon::maximum_possible_rectangles_on_first_page());

    auto alloc_data =
            allocator->create_new_tree_node<InlineUnboundedIsotheticPolygon>(
                    compute_sizeof_inline_unbounded_polygon(
                            overfull_rect_count),
                    NodeHandleType(
                            BIG_POLYGON));
    new (&(*alloc_data.first))
            InlineUnboundedIsotheticPolygon(allocator,
                                            overfull_rect_count);
    split.leftBranch.boundingPoly = alloc_data.second;
    alloc_data.first->push_polygon_to_disk(left_polygon);
  } else {
    split.leftBranch.boundingPoly =
            InlineBoundedIsotheticPolygon();
    std::get<InlineBoundedIsotheticPolygon>(split.leftBranch.boundingPoly).push_polygon_to_disk(left_polygon);
    assert(std::get<InlineBoundedIsotheticPolygon>(
            split.leftBranch.boundingPoly)
                   .get_summary_rectangle() ==
           left_polygon.boundingBox);
  }
  if (right_polygon.basicRectangles.size() > MAX_RECTANGLE_COUNT) {
    unsigned overfull_rect_count = (unsigned)(2 *
                                              right_polygon.basicRectangles.size());

    overfull_rect_count =
            std::min(overfull_rect_count,
                     (unsigned)InlineUnboundedIsotheticPolygon::maximum_possible_rectangles_on_first_page());

    auto alloc_data =
            allocator->create_new_tree_node<InlineUnboundedIsotheticPolygon>(
                    compute_sizeof_inline_unbounded_polygon(
                            overfull_rect_count),
                    NodeHandleType(
                            BIG_POLYGON));
    new (&(*alloc_data.first))
            InlineUnboundedIsotheticPolygon(allocator,
                                            overfull_rect_count);
    split.rightBranch.boundingPoly = alloc_data.second;
    alloc_data.first->push_polygon_to_disk(right_polygon);
  } else {
    split.rightBranch.boundingPoly =
            InlineBoundedIsotheticPolygon();
    std::get<InlineBoundedIsotheticPolygon>(split.rightBranch.boundingPoly).push_polygon_to_disk(right_polygon);
  }
  this->cur_offset_ = 0;

  return split;
#endif

  // Unsupported
  abort();
}

// Splitting a node will remove it from its this->parent node and its memory will be freed
template <int min_branch_factor, int max_branch_factor, class strategy>
SplitResult BranchNode<min_branch_factor, max_branch_factor, strategy>::splitNode() {
  SplitResult returnSplit = splitNode(partitionNode(), false);
  return returnSplit;
}

template <int min_branch_factor, int max_branch_factor, class strategy>
tree_node_handle BranchNode<min_branch_factor, max_branch_factor, strategy>::insert(
        std::variant<Branch, Point> &nodeEntry,
        std::vector<bool> &hasReinsertedOnLevel)
{
#if 0
  bool givenIsPoint = std::holds_alternative<Point>(nodeEntry);
  uint8_t stopping_level = 0;
  if (not givenIsPoint) {
    Branch &b = std::get<Branch>(nodeEntry);
    auto child_handle = b.child;
    if (child_handle.get_type() == LEAF_NODE || child_handle.get_type() == REPACKED_LEAF_NODE) {
      stopping_level = 1;
    } else {
      auto child_node = treeRef->get_branch_node(child_handle, false);
      stopping_level = child_node->level_ + 1;
    }
  }

  assert(this->self_handle_.get_type() == BRANCH_NODE);

  // Find the appropriate position for the entry
  // Should stop at appropriate depth level
  tree_node_handle current_handle = chooseNode(nodeEntry, stopping_level);
  tree_node_allocator *allocator = get_node_allocator(this->treeRef);
  SplitResult finalSplit;
  auto tree_ref_backup = treeRef;

  if (current_handle.get_type() == REPACKED_LEAF_NODE) {
    current_handle = treeRef->get_leaf_node(current_handle, true)->self_handle_;
  } else if (current_handle.get_type() == REPACKED_BRANCH_NODE) {
    current_handle = treeRef->get_branch_node(current_handle, true)->self_handle_;
  }

  if (std::holds_alternative<Point>(nodeEntry)) {
    assert(current_handle.get_type() == LEAF_NODE);
    auto current_node = treeRef->get_leaf_node(current_handle);

    current_node->addPoint(std::get<Point>(nodeEntry));

    finalSplit = adjustTreeSub(hasReinsertedOnLevel,
                               current_handle, treeRef);

  } else {
    assert(current_handle.get_type() == BRANCH_NODE);
    auto current_node = treeRef->get_branch_node(current_handle);
    Branch &sub_branch = std::get<Branch>(nodeEntry);

    IsotheticPolygon insertion_polygon =
            sub_branch.materialize_polygon(allocator);

    // Before I add this node in, I need to fragment everyone else
    // around it
    for (unsigned int i = 0; i < current_node->cur_offset_; i++) {
      Branch &b = current_node->entries.at(i);
      IsotheticPolygon branch_polygon = b.materialize_polygon(
              allocator);
      branch_polygon.increaseResolution(Point::atInfinity, insertion_polygon);
      branch_polygon.recomputeBoundingBox();
      update_branch_polygon(b, branch_polygon, treeRef);
    }

    current_node->addBranchToNode(sub_branch);
    if (sub_branch.child.get_type() == LEAF_NODE || sub_branch.child.get_type() == REPACKED_LEAF_NODE) {
      auto child_node = treeRef->get_leaf_node(sub_branch.child, true);
      assert(child_node->level_ == current_node->level_ - 1);
      child_node->parent = current_handle;
    } else {
      auto child_node = treeRef->get_branch_node(sub_branch.child, true);
      assert(child_node->level_ == current_node->level_ - 1);
      child_node->parent = current_handle;
    }

    finalSplit = adjustTreeSub(hasReinsertedOnLevel,
                               current_handle, treeRef);
  }

  // Grow the tree taller if we need to
  if (finalSplit.leftBranch.child != nullptr and finalSplit.rightBranch.child != nullptr) {
    auto alloc_data =
            allocator->create_new_tree_node<BranchNode<min_branch_factor, max_branch_factor, strategy>>(
                    NodeHandleType(BRANCH_NODE));
    new (&(*alloc_data.first)) BranchNode<min_branch_factor, max_branch_factor, strategy>(
            tree_ref_backup, tree_node_handle(nullptr), alloc_data.second, level_ + 1);
    auto new_root_handle = alloc_data.second;
    auto new_root_node = alloc_data.first;

    if (finalSplit.leftBranch.child.get_type() == LEAF_NODE || finalSplit.leftBranch.child.get_type() == REPACKED_LEAF_NODE) {
      auto left_node = tree_ref_backup->get_leaf_node(finalSplit.leftBranch.child, true);
      left_node->parent = new_root_handle;
    } else {
      auto left_node = tree_ref_backup->get_branch_node(finalSplit.leftBranch.child, true);
      left_node->parent = new_root_handle;
    }
    new_root_node->addBranchToNode(finalSplit.leftBranch);

    if (finalSplit.rightBranch.child.get_type() == LEAF_NODE || finalSplit.rightBranch.child.get_type() == REPACKED_LEAF_NODE) {
      auto right_node = tree_ref_backup->get_leaf_node(finalSplit.rightBranch.child, true);
      right_node->parent = new_root_handle;
    } else {
      auto right_node = tree_ref_backup->get_branch_node(finalSplit.rightBranch.child, true);
      right_node->parent = new_root_handle;
    }
    new_root_node->addBranchToNode(finalSplit.rightBranch);

    tree_ref_backup->root = new_root_handle;

    assert(this->self_handle_.get_type() == BRANCH_NODE);
    allocator->free(this->self_handle_, sizeof(BranchNode<min_branch_factor, max_branch_factor, strategy>));
    this->self_handle_ = tree_node_handle(nullptr);

    // Fix the reinserted length
    hasReinsertedOnLevel.push_back(false);

    return new_root_handle;
  }

  return tree_ref_backup->root;
#endif
}

// Always called on root, this = root
template <int min_branch_factor, int max_branch_factor, class strategy>
tree_node_handle BranchNode<min_branch_factor, max_branch_factor, strategy>::remove(Point givenPoint) {
#if 0
  // D1 [Find node containing record]
  tree_node_handle leaf_handle = findLeaf(givenPoint);
  // Record not in the tree
  if (leaf_handle == nullptr) {
    return this->self_handle_;
  }

  // D2 [Delete record]
  auto leaf_node = treeRef->get_leaf_node(leaf_handle);
  leaf_node->removePoint(givenPoint);

  // D3 [Propagate changes]
  leaf_node->condenseTree();

  // D4 [Shorten tree]
  assert(this->self_handle_.get_type() == BRANCH_NODE);
  if (this->cur_offset_ == 1) {
    tree_node_handle new_root_handle = entries.at(0).child;
    tree_node_allocator *allocator = get_node_allocator(this->treeRef);
    allocator->free(this->self_handle_, sizeof(BranchNode<min_branch_factor, max_branch_factor, strategy>));
    if (new_root_handle.get_type() == LEAF_NODE || new_root_handle.get_type() == REPACKED_LEAF_NODE) {
      auto new_root = treeRef->get_leaf_node(new_root_handle, true);
      new_root->parent = tree_node_handle(nullptr);
    } else {
      auto new_root = this->treeRef->get_branch_node(new_root_handle, true);
      new_root->parent = tree_node_handle(nullptr);
    }

    return new_root_handle;
  }

  return this->self_handle_;
#endif

  // Unsupported
  abort();
}

template <int min_branch_factor, int max_branch_factor, class strategy>
unsigned BranchNode<min_branch_factor, max_branch_factor, strategy>::checksum() {
#if 0
  unsigned sum = 0;

  for (unsigned i = 0; i < this->cur_offset_; i++) {
    // Recurse
    Branch &branch = entries.at(i);
    if (branch.child.get_type() == LEAF_NODE || branch.child.get_type() == REPACKED_LEAF_NODE) {
      auto child = this->treeRef->get_leaf_node(branch.child);
      sum += child->checksum();
    } else {
      auto child = this->treeRef->get_branch_node(branch.child);
      sum += child->checksum();
    }
  }

  return sum;
#endif

  // Unsupported
  abort();
}

template <int min_branch_factor, int max_branch_factor, class strategy>
std::vector<Point> BranchNode<min_branch_factor, max_branch_factor, strategy>::bounding_box_validate() {
#if 0
  tree_node_allocator *allocator = get_node_allocator(this->treeRef);
  std::vector<Point> all_child_points;
  for (unsigned i = 0; i < this->cur_offset_; i++) {
    Branch &b_i = entries.at(i);
    if (b_i.child.get_type() == LEAF_NODE || b_i.child.get_type() == REPACKED_LEAF_NODE) {
      auto child_ptr = this->treeRef->get_leaf_node(b_i.child);
      std::vector<Point> child_points =
              child_ptr->bounding_box_validate();
      for (Point &p : child_points) {
        all_child_points.push_back(p);
      }
    } else {
      auto child_ptr = this->treeRef->get_branch_node(b_i.child);
      std::vector<Point> child_points =
              child_ptr->bounding_box_validate();
      for (Point &p : child_points) {
        all_child_points.push_back(p);
      }
    }
  }
  if (this->parent != nullptr) {
    auto parent_node = treeRef->get_branch_node(parent);
    Branch &parent_branch = parent_node->locateBranch(this->self_handle_);
    IsotheticPolygon parent_poly =
            parent_branch.materialize_polygon(allocator);
    Rectangle bounding_box =
            parent_branch.get_summary_rectangle(allocator);
    for (Point &p : all_child_points) {
      if (!parent_poly.containsPoint(p)) {
        std::cout << "Parent poly " << this->parent << "does not contain: " << p
                  << std::endl;
        std::cout << "Poly was: " << parent_poly << std::endl;
        std::cout << "BB was: " << bounding_box << std::endl;
        std::cout << "My node is: " << this->self_handle_ << std::endl;
        abort();
      }
      if (!bounding_box.containsPoint(p)) {
        std::cout << "Parent poly contains " << p << " but the box does not!" << std::endl;
        abort();
      }
    }
  }
  return all_child_points;
#endif

  // Unsupported
  abort();
}

template <int min_branch_factor, int max_branch_factor, class strategy>
bool BranchNode<min_branch_factor, max_branch_factor, strategy>::validate(tree_node_handle expectedParent, unsigned index) {
#if 0
  tree_node_allocator *allocator = get_node_allocator(this->treeRef);

  if (expectedParent != nullptr and (this->parent != expectedParent ||
                                     this->cur_offset_ > max_branch_factor)) {
    std::cout << "node = " << (void *)this << std::endl;
    std::cout << "parent = " << this->parent << " expectedParent = " << expectedParent << std::endl;
    std::cout << "maxBranchFactor = " << max_branch_factor << std::endl;
    std::cout << "entries.size() = " << this->cur_offset_ << std::endl;
    assert(this->parent == expectedParent);
  }

  if (expectedParent != nullptr) {
    for (unsigned i = 0; i < this->cur_offset_; i++) {
      for (unsigned j = 0; j < this->cur_offset_; j++) {
        if (i != j) {
          Branch &b_i = entries.at(i);
          Branch &b_j = entries.at(j);
          IsotheticPolygon poly;
          if (std::holds_alternative<InlineBoundedIsotheticPolygon>(
                  b_i.boundingPoly)) {
            poly = std::get<InlineBoundedIsotheticPolygon>(
                    b_i.boundingPoly)
                    .materialize_polygon();
          } else {
            tree_node_handle poly_handle =
                    std::get<tree_node_handle>(b_i.boundingPoly);
            auto poly_pin = InlineUnboundedIsotheticPolygon::read_polygon_from_disk(
                    allocator, poly_handle);
            poly = poly_pin->materialize_polygon();
          }
          IsotheticPolygon poly2;
          if (std::holds_alternative<InlineBoundedIsotheticPolygon>(
                  b_j.boundingPoly)) {
            poly2 = std::get<InlineBoundedIsotheticPolygon>(
                    b_j.boundingPoly)
                    .materialize_polygon();
          } else {
            tree_node_handle poly_handle =
                    std::get<tree_node_handle>(b_j.boundingPoly);
            auto poly_pin2 =
                    InlineUnboundedIsotheticPolygon::read_polygon_from_disk(
                            allocator, poly_handle);
            poly2 = poly_pin2->materialize_polygon();
          }

          if (!poly.disjoint(poly2)) {
            std::cout << "Branch " << i << " is not disjoint from sibling Branch " << j << std::endl;
            std::cout << "Branch " << i << " " << b_i.child << std::endl;
            std::cout << "Branch " << j << " " << b_j.child << std::endl;

            std::cout << "Poly1 is: " << poly << std::endl;
            std::cout << "Poly2 is: " << poly2 << std::endl;
            std::cout << "Parent is: " << this->self_handle_ << std::endl;

            /*
                        std::cout << "branches[" << i << "].boundingPoly= " << b_i_poly << std::endl;
                        std::cout << "branches[" << j << "].boundingPoly = " << b_j_poly << std::endl;
                        */
            assert(false);
          }
        }
      }
    }
  }

  bool valid = true;
  for (unsigned i = 0; i < this->cur_offset_; i++) {
    Branch &b = entries.at(i);
    if (b.child.get_type() == LEAF_NODE || b.child.get_type() == REPACKED_LEAF_NODE) {
      auto child = this->treeRef->get_leaf_node(b.child, false);
      valid = valid and child->validate(this->self_handle_, i);
    } else {
      auto child = this->treeRef->get_branch_node(b.child, false);
      valid = valid and child->validate(this->self_handle_, i);
    }
  }

  return valid;
#endif

  // Unsupported
  abort();
}

template <int min_branch_factor, int max_branch_factor, class strategy>
void BranchNode<min_branch_factor, max_branch_factor, strategy>::print(unsigned n) {
#if 0
  std::string indentation(n * 4, ' ');
  std::cout << indentation << "Node " << (void *)this << std::endl;
  std::cout << indentation << "    Parent: " << this->parent << std::endl;
  std::cout << indentation << "    Branches: " << std::endl;
  for (unsigned i = 0; i < this->cur_offset_; i++) {
    Branch &branch = entries.at(i);
    // FIXME: out of band poly
    auto poly = std::get<InlineBoundedIsotheticPolygon>(
            branch.boundingPoly)
            .materialize_polygon();
    std::cout << indentation << "		" << branch.child << std::endl;
    std::cout << indentation << "		" << poly << std::endl;
  }
  std::cout << std::endl;
#endif

  // Unsupported
  abort();
}

template <int min_branch_factor, int max_branch_factor, class strategy>
void BranchNode<min_branch_factor, max_branch_factor, strategy>::printTree(unsigned n) {
  // Print this node first
  print(n);

  // Print any of our children with one more level of indentation
  std::string indendtation(n * 4, ' ');
  for (unsigned i = 0; i < this->cur_offset_; i++) {
    Branch &branch = entries.at(i);
    if (branch.child.get_type() == LEAF_NODE || branch.child.get_type() == REPACKED_LEAF_NODE) {
      auto child = this->treeRef->get_leaf_node(branch.child, false);

      // Recurse
      child->printTree(n + 1);
    } else {
      auto child = this->treeRef->get_branch_node(branch.child, false);

      // Recurse
      child->printTree(n + 1);
    }
  }
  std::cout << std::endl;
}

template <int min_branch_factor, int max_branch_factor, class strategy>
unsigned BranchNode<min_branch_factor, max_branch_factor, strategy>::height() {
#if 0
  unsigned ret = 0;
  tree_node_handle current_handle = this->self_handle_;

  for (;;) {
    ret++;
    if (current_handle.get_type() == LEAF_NODE || current_handle.get_type() == REPACKED_LEAF_NODE) {
      return ret;
    }

    auto node = treeRef->get_branch_node(current_handle, false);
    current_handle = node->entries.at(0).child;
  }
#endif
}

template <int min_branch_factor, int max_branch_factor, class strategy>
void stat_node(tree_node_handle start_handle, NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef) {
  std::stack<tree_node_handle> context;

  // Initialize our context stack
  context.push(start_handle);
  unsigned long polygonSize;
  unsigned long totalPolygonSize = 0;
  unsigned long totalLines = 0;
  size_t memoryFootprint = 0;
  size_t memoryPolygons = 0;
  unsigned long totalNodes = 0;
  unsigned long totalLeaves = 0;
  size_t deadSpace = 0;

  std::vector<unsigned long> histogramPolygon;
  histogramPolygon.resize(10000, 0);
  std::vector<unsigned long> histogramFanout;
  histogramFanout.resize(max_branch_factor, 0);

  double coverage = 0.0;

  tree_node_allocator *allocator = get_node_allocator(treeRef);

  while (!context.empty()) {
    auto currentContext = context.top();

    context.pop();
    totalNodes++;

    if (currentContext.get_type() == LEAF_NODE) {
      auto current_node = treeRef->get_leaf_node(currentContext);
      unsigned fanout = current_node->cur_offset_;

      if (fanout >= histogramFanout.size()) {
        histogramFanout.resize(2 * fanout, 0);
      }

      histogramFanout[fanout]++;
      totalLeaves++;
      memoryFootprint += sizeof(LeafNode<min_branch_factor, max_branch_factor, strategy>);
    } else if (currentContext.get_type() == BRANCH_NODE) {
      auto current_branch_node = treeRef->get_branch_node(currentContext);

      unsigned fanout = current_branch_node->cur_offset_;
      if (fanout >= histogramFanout.size()) {
        histogramFanout.resize(2 * fanout, 0);
      }
      histogramFanout[fanout]++;

      // Compute the overlap and coverage of our children
      for (unsigned i = 0; i < current_branch_node->cur_offset_; i++) {
        Branch &b = current_branch_node->entries.at(i);
        IsotheticPolygon polygon;

        auto itr = treeRef->polygons.find(b.child);
        if (itr == treeRef->polygons.end()) {
          polygon = IsotheticPolygon(b.boundingBox);
        } else {
          polygon = itr->second;
        }

        coverage += polygon.area();

        polygonSize = polygon.basicRectangles.size();
        assert(polygonSize < histogramPolygon.size());
        histogramPolygon[polygonSize]++;
        totalPolygonSize += polygonSize;

        // Compute space occupied by polygons
        // Ignore polygons with just one rectangle, we don't actually store
        // them in the map, they are constructed at run-time.
        if (polygon.basicRectangles.size() > 1) {
          memoryPolygons += polygon.computeMemory();
        }

        // FIXME: these stats are all wrong now.
        for (Rectangle r : polygon.basicRectangles) {
          if (r.area() == 0.0) {
            totalLines++;
          }
        }

        context.push(b.child);
      }

      memoryFootprint += sizeof(BranchNode<min_branch_factor, max_branch_factor, strategy>); // +
      // other out of line polys

      deadSpace += (sizeof(Branch) * (max_branch_factor - current_branch_node->cur_offset_));
    }
  }

  // Print out what we have found
  STATEXEC(std::cout << "### Statistics ###" << std::endl);
  // Memory footprint is wrong!
  STATMEM(memoryFootprint);
  STATMEM(memoryPolygons);
  //STATHEIGHT(height());
  STATSIZE(totalNodes);
  //STATEXEC(std::cout << "DeadSpace: " << deadSpace << std::endl);
  //STATSINGULAR(singularBranches);
  STATLEAF(totalLeaves);
  STATBRANCH(totalNodes - totalLeaves);
  STATCOVER(coverage);
  STATFANHIST();
  for (unsigned i = 0; i < histogramFanout.size(); ++i) {
    if (histogramFanout[i] > 0) {
      STATHIST(i, histogramFanout[i]);
    }
  }
  STATLINES(totalLines);
  STATTOTALPOLYSIZE(totalPolygonSize);
  STATPOLYHIST();
  for (unsigned i = 0; i < histogramPolygon.size(); ++i) {
    if (histogramPolygon[i] > 0) {
      STATHIST(i, histogramPolygon[i]);
    }
  }
  std::cout << treeRef->stats;

  STATEXEC(std::cout << "### ### ### ###" << std::endl);
}

enum LookAheadMergeCmd {
    ADD = 0,
    STOP = 1,
    CREATE_VERTICAL = 2
};

inline bool is_vertical_stripe(
        std::vector<std::pair<Point, uint8_t>> &points_with_ownership,
        unsigned i,
        unsigned dimension) {
  static_assert(dimensions == 2);

  uint8_t existing_ownership = points_with_ownership.at(i).second;
  double existing_value =
          points_with_ownership.at(i).first[dimension];

  for (unsigned j = i + 1; j < points_with_ownership.size(); j++) {
    double next_value =
            points_with_ownership.at(j).first[dimension];
    uint8_t next_ownership = points_with_ownership.at(j).second;

    if (next_value != existing_value) {
      // Value changed before ownership did. It's not a vertical
      // stripe.
      return false;
    }

    if (next_ownership != existing_ownership) {
      // Ownership changed before value did. We have two owners of
      // values in the same coord, so it needs to split on another
      // dimension.
      return true;
    }

    // Both ownership and value are the same. Keep going.
  }
  return false;
}

inline LookAheadMergeCmd get_merge_cmd(
        std::vector<std::pair<Point, uint8_t>> &points_with_ownership,
        unsigned i,
        unsigned dimension,
        uint8_t cur_ownership) {
  // N.B.: valid ownerships are 0 or 1. A 2 implies we have no
  // operating point, so we can take either ownership.

  // Need to reconsider this with higher dimensions
  static_assert(dimensions == 2);

  // Can just add this to whatever ongoing rectangle we have
  // if it the last point and it is owned by us
  if (i == points_with_ownership.size() - 1 and
      (points_with_ownership.at(i).second == cur_ownership or
       cur_ownership == 2)) {
    return ADD;
  }

  if (is_vertical_stripe(points_with_ownership, i, dimension)) {
    return CREATE_VERTICAL;
  }

  if (points_with_ownership.at(i).second == cur_ownership or
      cur_ownership == 2) {
    return ADD;
  } else {
    return STOP;
  }
}

inline std::pair<std::vector<Rectangle>, std::vector<Rectangle>>
walk_and_dice_overlapping_region(
        std::vector<std::pair<Point, uint8_t>> &points_with_ownership,
        unsigned dimension_to_walk,
        Rectangle overlapping_dimensions) {
  std::vector<Rectangle> ret_a_rects;
  std::vector<Rectangle> ret_b_rects;

  Point lower_point = overlapping_dimensions.lowerLeft;
  Point upper_point = overlapping_dimensions.upperRight;

  // Need the largest point below the upper right corner
  // for the purposes of point-basd rectangle expansion
  Point max_included_point = Point::closest_smaller_point(upper_point);

  // Magic Values indicating that we are not currently operating on a
  // point
  Rectangle cur_rect = Rectangle::atInfinity;
  uint8_t cur_ownership = 2;

  // Containers for vertical splits
  std::vector<std::vector<std::pair<Point, uint8_t>>> vertical_splits;
  std::vector<Rectangle> vertical_rects;

  for (unsigned i = 0; i < points_with_ownership.size(); i++) {
    LookAheadMergeCmd cmd = get_merge_cmd(points_with_ownership, i,
                                          dimension_to_walk,
                                          cur_ownership);
    if (cmd == ADD) {
      if (cur_rect == Rectangle::atInfinity) {
        // If we aren't yet operating on a point, then start
        // with this point.
        cur_ownership = points_with_ownership.at(i).second;
        lower_point[dimension_to_walk] =
                points_with_ownership.at(i).first[dimension_to_walk];
        upper_point[dimension_to_walk] = nextafter(
                points_with_ownership.at(i).first[dimension_to_walk],
                DBL_MAX);
        cur_rect = Rectangle(lower_point, upper_point);
      } else {
        // If we are already operating on a point, then expand
        // the rectangle to include the next point
        max_included_point[dimension_to_walk] =
                points_with_ownership.at(i).first[dimension_to_walk];
        cur_rect.expand(max_included_point);
      }
    } else if (cmd == STOP) {
      // STOP means push the current rectangle on and start a new
      // one from the next point
      assert(cur_rect != Rectangle::atInfinity);
      // Adding existing rect
      if (cur_ownership == 0) {
        ret_a_rects.push_back(cur_rect);
      } else {
        ret_b_rects.push_back(cur_rect);
      }
      // Start new rect constrained to the starting point in this
      // dimension
      lower_point[dimension_to_walk] =
              points_with_ownership.at(i).first[dimension_to_walk];
      upper_point[dimension_to_walk] = nextafter(
              points_with_ownership.at(i).first[dimension_to_walk],
              DBL_MAX);
      cur_rect = Rectangle(lower_point, upper_point);
      cur_ownership = points_with_ownership.at(i).second;
    } else if (cmd == CREATE_VERTICAL) {
      // CREATE_VERTICAL means push the current rectangle on, if
      // we have one, and then create a vertical split record
      //Adding new rect
      if (cur_rect != Rectangle::atInfinity) {
        if (cur_ownership == 0) {
          ret_a_rects.push_back(cur_rect);
        } else {
          ret_b_rects.push_back(cur_rect);
        }
      } else {
        assert(cur_ownership == 2);
      }

      // Start new rect constrained to the starting point in this
      // dimension. This is a vertical stripe, which we will split
      std::vector<std::pair<Point, uint8_t>> vertical_data;
      vertical_data.push_back(points_with_ownership.at(i));
      double existing_val =
              points_with_ownership.at(i).first[dimension_to_walk];
      lower_point[dimension_to_walk] =
              points_with_ownership.at(i).first[dimension_to_walk];
      upper_point[dimension_to_walk] = nextafter(
              points_with_ownership.at(i).first[dimension_to_walk],
              DBL_MAX);
      Rectangle sub_overlap_rect(lower_point, upper_point);
      i++;
      while (i < points_with_ownership.size() and
             points_with_ownership.at(i).first[dimension_to_walk] == existing_val) {
        vertical_data.push_back(points_with_ownership.at(i));
        i++;
      }
      vertical_splits.push_back(vertical_data);
      vertical_rects.push_back(sub_overlap_rect);

      // Unset everything. Figure out what we will do next.
      cur_rect = Rectangle::atInfinity;
      cur_ownership = 2;
      i--;
    } else {
      assert(false);
    }
  }

  if (cur_ownership == 0) {
    ret_a_rects.push_back(cur_rect);
  } else if (cur_ownership == 1) {
    ret_b_rects.push_back(cur_rect);
  }

  // Do any subsplits if req'd
  for (unsigned i = 0; i < vertical_splits.size(); i++) {
    auto ret_data = walk_and_dice_overlapping_region(
            vertical_splits.at(i),
            dimension_to_walk + 1,
            vertical_rects.at(i));
    for (auto &rect : ret_data.first) {
      ret_a_rects.push_back(rect);
    }
    for (auto &rect : ret_data.second) {
      ret_b_rects.push_back(rect);
    }
  }
  return std::make_pair(ret_a_rects, ret_b_rects);
}

template <class NT>
std::pair<std::vector<Rectangle>, std::vector<Rectangle>> make_rectangles_disjoint_accounting_for_region_ownership(
        NT *treeRef,
        Rectangle &a,
        tree_node_handle &a_node,
        Rectangle &b,
        tree_node_handle &b_node) {

  std::vector<Rectangle> ret_a_rects;
  std::vector<Rectangle> ret_b_rects;
  if (not a.intersectsRectangle(b)) {
    ret_a_rects.push_back(a);
    ret_b_rects.push_back(b);
    return std::make_pair(ret_a_rects, ret_b_rects);
  }

  Rectangle overlapping_rect = a.intersection(b);
  assert(overlapping_rect != Rectangle::atInfinity);

  std::vector<Point> a_owned_points = rectangle_search(a_node,
                                                       overlapping_rect, treeRef, false /* don't track */);
  // If no owned points, then yield the region to b
  if (a_owned_points.empty()) {
    ret_a_rects = a.fragmentRectangle(overlapping_rect);
    ret_b_rects.push_back(b);
    return std::make_pair(ret_a_rects, ret_b_rects);
  }
  Rectangle a_owned_rect(a_owned_points[0],
                         Point::closest_larger_point(a_owned_points[0]));
  for (unsigned i = 1; i < a_owned_points.size(); i++) {
    a_owned_rect.expand(a_owned_points.at(i));
  }

  std::vector<Point> b_owned_points = rectangle_search(b_node,
                                                       overlapping_rect, treeRef, false /* don't track */);
  // If no owned points, yield the region to a
  if (b_owned_points.empty()) {
    ret_a_rects.push_back(a);
    ret_b_rects = b.fragmentRectangle(overlapping_rect);
    return std::make_pair(ret_a_rects, ret_b_rects);
  }

  Rectangle b_owned_rect(b_owned_points[0],
                         Point::closest_larger_point(b_owned_points[0]));
  for (unsigned i = 1; i < b_owned_points.size(); i++) {
    b_owned_rect.expand(b_owned_points.at(i));
  }

  Rectangle ownership_overlap = a_owned_rect.intersection(
          b_owned_rect);
  // No overlap in one dimension. We can just add the owned rectangles
  // to whatever we have after we fragment about the overalpping areas
  if (ownership_overlap == Rectangle::atInfinity) {
    // see above and return
    ret_a_rects = a.fragmentRectangle(overlapping_rect);
    ret_b_rects = b.fragmentRectangle(overlapping_rect);

    ret_a_rects.push_back(a_owned_rect);
    ret_b_rects.push_back(b_owned_rect);
    return std::make_pair(ret_a_rects, ret_b_rects);
  }

  // Now, its time to get fancy
  // Sort by each dimension, and create bounding boxes in that
  // dimension aligned with the overlapping region such that the boxes
  // start from the point of interest and stop just before the point
  // when sorted along that dimension. If two points have the same
  // value, note that and skip it for now, we will need to do a
  // 'horizontal split' after that.

  // Fragment out the overlapping region from each rectangle.
  ret_a_rects = a.fragmentRectangle(overlapping_rect);
  ret_b_rects = b.fragmentRectangle(overlapping_rect);

  std::vector<std::pair<Point, uint8_t>> points_with_ownership;
  for (const Point p : a_owned_points) {
    points_with_ownership.push_back(std::make_pair(p, 0));
  }
  for (const Point p : b_owned_points) {
    points_with_ownership.push_back(std::make_pair(p, 1));
  }

  // Sort on X
  std::sort(points_with_ownership.begin(),
            points_with_ownership.end(), [](std::pair<Point, uint8_t> &p1, std::pair<Point, uint8_t> &p2) {
              if (p1.first[0] != p2.first[0]) {
                return p1.first[0] < p2.first[0];
              }
              return p1.first[1] < p2.first[1];
          });

  auto ret_data = walk_and_dice_overlapping_region(
          points_with_ownership, 0, overlapping_rect);

  auto additional_a_rects = ret_data.first;
  for (auto &rect : additional_a_rects) {
    ret_a_rects.push_back(rect);
  }
  auto additional_b_rects = ret_data.second;
  for (auto &rect : additional_b_rects) {
    ret_b_rects.push_back(rect);
  }

  return std::make_pair(ret_a_rects, ret_b_rects);
}

static std::vector<Point> tree_validate_recursive(tree_node_handle current_handle, tree_node_allocator *allocator) {
  switch (current_handle.get_type()) {
    case LEAF_NODE:
    case BRANCH_NODE:
      assert(false);
    default: {
      assert(false);
    }
  }
}
} // namespace nirtreedisk
