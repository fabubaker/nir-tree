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

#define NDEBUG
// #define DEBUG_TEST
// #define DEBUG_TESTDISJOINT
// #define DEBUG_TESTCONTAINPOINTS
// #define DEBUG_TESTCOUNT


namespace nirtreedisk {

// Different strategies for partitioning BranchNode during splitNode
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


// Branch object contains child and boundingBox where child is the 
// tree_node_handler which points to the disk page, and boundingBox 
// is the MBB of this branch's polygon. child is also the key of this 
// branch's polygon within polygons map. 
struct Branch {
  Branch(Rectangle boundingBox, tree_node_handle child_handle):
    boundingBox(boundingBox), child(child_handle) {}
  Branch() = default;
  Branch(const Branch &other): boundingBox(other.boundingBox), child(other.child) {}
  bool operator==(const Branch &o) const = default;
  bool operator!=(const Branch &o) const = default;
  // members: 
  Rectangle boundingBox;
  tree_node_handle child;
};

// BranchAtLevel object is for inserting a branch at a specific
// tree level
// [LEGACY] for insert branch 
struct BranchAtLevel {
  BranchAtLevel(Branch branch, uint8_t level):
    branch(branch), level(level) {}
  BranchAtLevel() = default;
  BranchAtLevel(const BranchAtLevel &other): branch(other.branch), level(other.level) {}
  bool operator==(const BranchAtLevel &o) const = default;
  bool operator!=(const BranchAtLevel &o) const = default;
  // members: 
  Branch branch;
  uint8_t level;
};

// SplitResult object represents two newly created branches after 
// splitting an overflowed Branch
struct SplitResult {
  Branch leftBranch;
  Branch rightBranch;
};

// Partition object specifies on which dimension and which location
// the partition of Branch/Point should have be done. This is the 
// result of PartitionLeafNode() or PartitionBranchNode()
struct Partition {
  unsigned dimension;
  double location;
};

// Helper functions for polygon
template <int min_branch_factor, int max_branch_factor, class strategy>
IsotheticPolygon find_polygon(
        NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
        tree_node_handle node_handle,
        Rectangle rectangle); 

template <int min_branch_factor, int max_branch_factor, class strategy>
IsotheticPolygon find_polygon(
        NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
        Branch branch);

template <int min_branch_factor, int max_branch_factor, class strategy>
void remove_polygon(
        NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
        tree_node_handle node_handle);

template <int min_branch_factor, int max_branch_factor, class strategy>
void update_branch_polygon(
        NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
        Branch &branch_to_update,
        IsotheticPolygon &polygon_to_write);
  
// [LEGACY] for insert branch 
std::pair<double, std::vector<IsotheticPolygon::OptimalExpansion>>
computeExpansionArea(const IsotheticPolygon &this_poly, const IsotheticPolygon &other_poly);

// Helper functions for debugging:
template <int min_branch_factor, int max_branch_factor, class strategy>
void testDisjoint(
        NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
        tree_node_handle root,
        std::string msg="");
template <int min_branch_factor, int max_branch_factor, class strategy>
void testCount(
        NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
        tree_node_handle root,
        std::string msg="");
template <int min_branch_factor, int max_branch_factor, class strategy>
void testContainPoints(
        NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
        tree_node_handle root,
        std::string msg="");
  
// Helper function for finding a path from root to leaf node (both are inclusive)
template <int min_branch_factor, int max_branch_factor, class strategy>
bool find_path_to_leaf_helper(
        NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
        tree_node_handle start_handle,
        tree_node_handle leaf_handle,
        std::stack<tree_node_handle> &path_to_leaf);
template <int min_branch_factor, int max_branch_factor, class strategy>
std::stack<tree_node_handle> find_path_to_leaf(
        NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
        tree_node_handle start_handle,
        tree_node_handle leaf_handle);

// Helper function for adjust/split tree:
template <class NT, class TR>
std::pair<SplitResult, tree_node_handle> adjust_tree_bottom_half(
        NT current_node,
        TR *treeRef,
        tree_node_handle current_handle,
        std::stack<tree_node_handle> &parentHandles,
        std::vector<bool> &hasReinsertedOnLevel,
        int max_branch_factor);

template <int min_branch_factor, int max_branch_factor, class strategy>
SplitResult adjustTreeSub(
        NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
        tree_node_handle start_handle, 
        std::stack<tree_node_handle> &parentHandles,
        std::vector<bool> &hasReinsertedOnLevel);
// Functions for search: 
// [TODO]

// LeafNode object contains array of Points with name of entries   
// cur_offset_ specifis the current count of Points in array entries
  // [FIXME] [TODO]
  // reinsert()
  // validate()
template <int min_branch_factor, int max_branch_factor, class strategy>
requires(std::derived_from<strategy, BranchPartitionStrategy>)
class LeafNode {
public:
  // members: 
  // one extra space is saved for insertion which causes split 
  std::array<Point, max_branch_factor + 1> entries;
  unsigned cur_offset_;

  LeafNode(): cur_offset_(0) {
    static_assert(sizeof(LeafNode<min_branch_factor, max_branch_factor, strategy>) <= PAGE_DATA_SIZE);
  }
  void deleteSubtrees();

  // Data structure interface functions : 
  // Insert a Point on treeRef where selfHandle is the handle of *this node 
  tree_node_handle insert(
          NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
          tree_node_handle selfHandle,
          Point givenPoint, 
          std::vector<bool> &hasReinsertedOnLevel);
  // [FIXME]
  void reInsert(std::vector<bool> &hasReinsertedOnLevel);
  // Remove a Point on treeRef where selfHandle is the handle of *this node 
  tree_node_handle remove(
          NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
          tree_node_handle selfHandle,
          Point givenPoint);

  // Helper Functions: 
  // addPoints: add Point at the end of the array and increase cur_offset_
  void addPoint(const Point &point) {
    entries.at(this->cur_offset_++) = point;
  }
  // removePoint: remove the givenPoint from array 
  void removePoint(const Point &point);
  // chooseNode: chooses the Leaf Node to add given Point which has to be itself
  tree_node_handle chooseNode(tree_node_handle selfHandle, Point givenPoint);
  // findLeaf: returns itself if it contains givenPoint or nullptr if not 
  tree_node_handle findLeaf(
          NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
          tree_node_handle selfHandle,
          Point givenPoint);
  // Entry function for partitionLeafNode()
  Partition partitionNode();
  // partitionLeafNode: return the best dimension and location to partition 
  // on a LeafNode which is the dimension with the highest varaince and location
  // at average of all points 
  Partition partitionLeafNode();
  // Entry function for splitNode()
  SplitResult splitNode(
          NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
          tree_node_handle current_handle,
          tree_node_handle parent_handle);
  // splitNode: splits Leafnode with current_handle into two LeafNode objects 
  // according to partition p
  SplitResult splitNode(
          NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
          tree_node_handle current_handle,
          tree_node_handle parent_handle,
          Partition p, 
          bool is_downsplit);

  // this is called by BranchNode::remove() to shorten a tree
  void condenseTree(
          NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
          tree_node_handle selfHandle,
          std::stack<tree_node_handle> &parentHandles);

  // Miscellaneous
  // boundingBox() returns a MBB of all points on LeafNode 
  Rectangle boundingBox();
  // checksum: return sum of all points for all dimensions 
  unsigned checksum();
  // [FIXME] : called by branchnode validate [TODO]
  bool validate(tree_node_handle expectedParent, unsigned index);
  // bounding_box_validate: returns a copy of vector of Points
  std::vector<Point> bounding_box_validate();
  void print(tree_node_handle current_handle, tree_node_handle parent_handle, unsigned n = 0);
  void printTree(tree_node_handle current_handle, tree_node_handle parent_handle, unsigned n = 0);
  // height: returns 1 for Leaf Node
  unsigned height();
};


// BranchNode object contains array of Branches with name of entries   
// cur_offset_ specifies the current count of Branches in array entries
  // to be fixed: 
  // reinsert()  
  // validate()
  // bounding_box_validate()
template <int min_branch_factor, int max_branch_factor, class strategy>
requires(std::derived_from<strategy, BranchPartitionStrategy>)
class BranchNode {
public:
  // members: 
  // have a space for possible overflow
  std::array<Branch, max_branch_factor + 1> entries;
  unsigned cur_offset_;

  // Constructors and destructors
  BranchNode(): cur_offset_(0) {
  static_assert(sizeof(BranchNode<min_branch_factor, max_branch_factor, strategy>) <= PAGE_DATA_SIZE);
  }
  void deleteSubtrees(NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef);

  // Data structure interface functions: 
  // insert a Point/BranchAtLevel where selfHandle is root 
  //      BranchAtLevel is not considered for now 
  tree_node_handle insert(
        NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
        tree_node_handle selfHandle,
        std::variant<BranchAtLevel, Point> &nodeEntry, 
        std::vector<bool> &hasReinsertedOnLevel);
  // [FIXME]
  void reInsert(std::vector<bool> &hasReinsertedOnLevel);
  // remove a point from tree where selfHandle is root 
  tree_node_handle remove(
        NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
        tree_node_handle selfHandle,
        Point givenPoint);
  
  // Helper functions
  // locateBranch: locate Branch at BranchNode
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
  // addBranchToNode: add Branch at the end of array and increase cur_offset_
  void addBranchToNode(const Branch &entry) {
    entries.at(this->cur_offset_++) = entry;
  };
  // removeBranch: remove branch from BranchNode and free the memory associated with Branch
  // as well as removing polygon associated with this node from map 
  void removeBranch(
          NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
          const tree_node_handle entry);

  // choose a LeafNode for adding a point 
  // expansion and clipping of polygon are also done here 
  tree_node_handle chooseNodePoint(
          NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
          tree_node_handle selfHandle,
          std::stack<tree_node_handle> &parentHandles,
          Point &point);
  // this function is untested and not used for now 
  tree_node_handle chooseNodeBranch(
          NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
          tree_node_handle selfHandle,
          std::stack<tree_node_handle> &parentHandles,
          BranchAtLevel &branchLevel);
  // findLeaf: returns the LeafNode which contains the point or nullptr if none node contains it 
  tree_node_handle findLeaf(
          NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
          tree_node_handle selfHandle,
          Point givenPoint);
  // entry function for partitionPartitionNode()
  Partition partitionNode(NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef);

  // clip polygon if it overlaps with its siblings and ignore polygon 
  // assocaited with handle_to_skip
  void make_disjoint_from_children(
          NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
          tree_node_handle handle_to_skip,
          IsotheticPolygon &polygon);
  // entry function for splitNode()
  SplitResult splitNode(
          NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
          tree_node_handle current_handle,
          tree_node_handle parent_handle);
  // splitNode: splits BranchNode with current_handle into two BranchNode object according to p 
  SplitResult splitNode(
          NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
          tree_node_handle current_handle,
          tree_node_handle parent_handle,
          Partition p, 
          bool is_downsplit);

  //void condenseTree();

  // Miscellaneous
  // boundingBox() returns a MBB of all branches on BranchNode
  Rectangle boundingBox();
  // checksum: return sum of all points for each dimension of all subtrees 
  unsigned checksum(NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef);
  // [FIXME]
  bool validate(tree_node_handle expectedParent, unsigned index);
  // [FIXME]
  std::vector<Point> bounding_box_validate();
  void print(NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
            tree_node_handle current_handle, tree_node_handle parent_handle, unsigned n = 0);
  void printTree(NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
            tree_node_handle current_handle, tree_node_handle parent_handle, unsigned n = 0);
  // height: returns the height of subtree where LeafNode has height 1 
  unsigned height(NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
                  tree_node_handle selfHandle);


#if 0
template <class S = strategy>
  Partition partitionBranchNode(
            typename std::enable_if<std::is_same<S, LineMinimizeDownsplits>::value, S>::type * = 0) {
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
  }

template <class S = strategy>
  Partition partitionBranchNode(
            typename std::enable_if<std::is_same<S, LineMinimizeDistanceFromMean>::value,S>::type * = 0) {
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
  }
#endif 

std::pair<bool, Partition> try_cut_geo_mean(std::vector<Rectangle> &all_branch_polys) {
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
  }

template <class S = strategy>
Partition partitionBranchNode(
          NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
          typename std::enable_if<std::is_same<S, ExperimentalStrategy>::value,S>::type * = 0) {
    Partition defaultPartition;

    tree_node_allocator *allocator = get_node_allocator(treeRef);
    std::vector<Rectangle> all_branch_polys;
    for (unsigned i = 0; i < this->cur_offset_; i++) {
      Branch &b_i = entries.at(i);
      IsotheticPolygon b_poly = find_polygon(treeRef, b_i);
      all_branch_polys.push_back(b_poly.boundingBox);
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
                    for (unsigned i = 0; i < dimensions; i++) {
                      if (poly1.upperRight[i] == poly2.upperRight[i]) {
                        continue;
                      }
                      return poly1.upperRight[i] < poly2.upperRight[i];
                    }
                  }
                  return poly1.upperRight[d] < poly2.upperRight[d];
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
                Rectangle poly1 = b1.boundingBox;
                Rectangle poly2 = b2.boundingBox;
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
}

};


// root is freed separately, no work to do 
template <int min_branch_factor, int max_branch_factor, class strategy>
void LeafNode<min_branch_factor, max_branch_factor, strategy>::deleteSubtrees() {
  return;
}

template <int min_branch_factor, int max_branch_factor, class strategy>
Rectangle LeafNode<min_branch_factor, max_branch_factor, strategy>::boundingBox() {
  assert(this->cur_offset_ > 0);
  Point &p = entries.at(0);
  Rectangle bb(p, Point::closest_larger_point(p));
  for (unsigned i = 1; i < this->cur_offset_; i++) {
    bb.expand(entries.at(i));
  }
  return bb;
}

// Removes the point from the node
template <int min_branch_factor, int max_branch_factor, class strategy>
void LeafNode<min_branch_factor, max_branch_factor, strategy>::removePoint( const Point &point) {
  // Locate the child
  unsigned childIndex = 0;
  while( entries.at(childIndex) != point and childIndex < this->cur_offset_ ){
    childIndex++; 
  }
  // point is expected to exist on LeafNode 
  assert(entries.at(childIndex) == point);

  // Replace this index with whatever is in the last position
  entries.at(childIndex) = entries.at(this->cur_offset_ - 1);

  // Truncate array size
  this->cur_offset_--;
}

// ??? what is the function used for ???
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

// Always called on root, this = root and this is the only node in tree 
template <int min_branch_factor, int max_branch_factor, class strategy>
tree_node_handle LeafNode<min_branch_factor, max_branch_factor, strategy>::chooseNode(
        tree_node_handle selfHandle,
        Point givenPoint) {
  return selfHandle;
}

template <int min_branch_factor, int max_branch_factor, class strategy>
tree_node_handle LeafNode<min_branch_factor, max_branch_factor, strategy>::findLeaf(
        NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
        tree_node_handle selfHandle,
        Point givenPoint) {
  // Check each entry to see if it matches point
  for (unsigned i = 0; i < this->cur_offset_; i++) {
    Point &p = std::get<Point>(entries.at(i));
    if (p == givenPoint) {
      return selfHandle;
    }
  }
  return tree_node_handle(nullptr);
}

// Returns the best dimension and location to partition on a LeafNode 
//    best partition dimension is the dimension with the highest variance 
//    best partition location is the average of all points on the best dimension
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
// The old one is freed in adjust_tree_bottom_half using removeBranch
// If we downsplit, then we won't call adjustTree for that split so we
// need to delete the node ourselves.
template <int min_branch_factor, int max_branch_factor, class strategy>
SplitResult LeafNode<min_branch_factor, max_branch_factor, strategy>::splitNode(
              NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
              tree_node_handle current_handle,
              tree_node_handle parent_handle,
              Partition p, 
              bool is_downsplit) {

  assert(current_handle.get_type() == LEAF_NODE);
  tree_node_allocator *allocator = get_node_allocator(treeRef);

  // Allocate a leaf node for left node
  auto alloc_data = allocator->create_new_tree_node<LeafNode<min_branch_factor, max_branch_factor, strategy>>(
                    NodeHandleType(LEAF_NODE));
  tree_node_handle left_handle = alloc_data.second;
  auto left_node = alloc_data.first; // take pin
  new (&(*left_node)) LeafNode<min_branch_factor, max_branch_factor, strategy>();
  assert(left_handle.get_type() == LEAF_NODE);
  
  // Allocate a leaf node for right node
  alloc_data = allocator->create_new_tree_node<LeafNode<min_branch_factor, max_branch_factor, strategy>>(
                NodeHandleType(LEAF_NODE));
  tree_node_handle right_handle = alloc_data.second;
  auto right_node = alloc_data.first; // take pin
  new (&(*right_node)) LeafNode<min_branch_factor, max_branch_factor, strategy>();
  assert(right_handle.get_type() == LEAF_NODE);

  SplitResult split = {{Rectangle(), left_handle},
                       {Rectangle(), right_handle}};

  // Split all data points in current LeafNode into one of left or right node 
  bool containedLeft, containedRight;
  for (unsigned i = 0; i < this->cur_offset_; i++) {
    Point &dataPoint = entries.at(i);
    containedLeft = dataPoint[p.dimension] < p.location; // Not inclusive
    containedRight = dataPoint[p.dimension] >= p.location;
    if (containedLeft and not containedRight) {
      left_node->addPoint(dataPoint);
    } else if (not containedLeft and containedRight) {
      right_node->addPoint(dataPoint);
    } else {
      assert(false);
    }
  }

  // All points have been routed.
  IsotheticPolygon left_polygon(left_node->boundingBox());
  IsotheticPolygon right_polygon(right_node->boundingBox());
  // left and right polygon should be disjoint 
  assert(left_polygon.disjoint(right_polygon));

  // If we have a parent, we need to make these disjoint from our
  // siblings. If we don't, then we are automatically disjoint
  // from our siblings since these are the only two polys and they
  // are disjoint from each other now.
  if (parent_handle) {
    auto parent_node = treeRef->get_branch_node(parent_handle);
    if (not is_downsplit) {
      // make left_polygon disjoint from its siblings 
      parent_node->make_disjoint_from_children(treeRef, current_handle, left_polygon);
      assert(left_polygon.basicRectangles.size() > 0);
      left_polygon.refine(); //remove duplicated rectangles
      assert(left_polygon.basicRectangles.size() > 0);

      // make right_polygon disjoint from its siblings 
      parent_node->make_disjoint_from_children(treeRef, current_handle, right_polygon);
      assert(right_polygon.basicRectangles.size() > 0);
      right_polygon.refine(); //remove duplicated rectangles
      assert(right_polygon.basicRectangles.size() > 0);
    } else {
      // Intersect with our existing poly to avoid intersect
      // with other siblings, as the existing polygon is disjoint
      // with its siblings 
      IsotheticPolygon current_polygon = find_polygon(treeRef, current_handle, 
                                          this->boundingBox());
      assert(current_polygon.basicRectangles.size() > 0);
      
      // left side
      assert(left_polygon.basicRectangles.size() > 0);
      IsotheticPolygon poly_backup = left_polygon;
      left_polygon.intersection(current_polygon);
      if (left_polygon.basicRectangles.size() == 0) {
        std::cout << "Weird situation: " << poly_backup << " is disjoint from existing polygon: " << current_polygon << std::endl;
      }
      assert(left_polygon.basicRectangles.size() > 0);
      left_polygon.refine();
      assert(left_polygon.basicRectangles.size() > 0);

      // right side
      assert(right_polygon.basicRectangles.size() > 0);
      right_polygon.intersection(current_polygon);
      assert(right_polygon.basicRectangles.size() > 0);
      right_polygon.refine();
      assert(right_polygon.basicRectangles.size() > 0);
      assert(left_polygon.disjoint(right_polygon));
    }
  }
  update_branch_polygon(treeRef, split.leftBranch, left_polygon);
  update_branch_polygon(treeRef, split.rightBranch, right_polygon);
  // this will be deleted in caller adjust_tree_bottom_half()
  this->cur_offset_ = 0;
  return split;
}

template <int min_branch_factor, int max_branch_factor, class strategy>
SplitResult LeafNode<min_branch_factor, max_branch_factor, strategy>::splitNode(
              NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
              tree_node_handle current_handle,
              tree_node_handle parent_handle) {
  SplitResult returnSplit = splitNode(treeRef, current_handle, parent_handle, partitionNode(), false);
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

// ??? what is function used for ???
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
        NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
        Branch &branch_to_update,
        IsotheticPolygon &polygon_to_write)
{
  // Update boundingBox of Branch to be the same as bounding box of polygons
  branch_to_update.boundingBox = polygon_to_write.boundingBox;

  // If the MBR has not been split into a polygon, don't keep it in the map.
  if (polygon_to_write.basicRectangles.size() != 1) {
    auto insert_res = treeRef->polygons.insert({branch_to_update.child, polygon_to_write});
    if(not insert_res.second) {
      // Already exists in the map, update instead of insertion
      treeRef->polygons[branch_to_update.child] = polygon_to_write; 
    }
  } else {
    remove_polygon(treeRef, branch_to_update.child);
  }
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


// Split current_node if needed and return split result to caller 
template <class NT, class TR>
std::pair<SplitResult, tree_node_handle> adjust_tree_bottom_half(
        NT current_node,
        TR *treeRef,
        tree_node_handle current_handle,
        std::stack<tree_node_handle> &parentHandles,
        std::vector<bool> &hasReinsertedOnLevel,
        int max_branch_factor)
{
  SplitResult propagationSplit = {
          {Rectangle(), tree_node_handle(nullptr)},
          {Rectangle(), tree_node_handle(nullptr)}};

  if (current_node->cur_offset_ <= (unsigned)max_branch_factor) {
    // cur_offset_ is the index of the next inserted item which equals
    // to the number of items in entries
    // max_branch_factor is the max allowed number of elements in entries
    // which is inclusive, therefore, we allow equality relationship 
    // added branch fits into current branchNode, no split is needed 
    // tree_node_handle(nullptr) will break the while loop in caller 
    return std::make_pair(propagationSplit, tree_node_handle(nullptr));
  }

#ifdef DEBUG_TEST
  if (current_handle.get_type() == BRANCH_NODE ){
    testDisjoint(treeRef, current_handle, "before splitNode");
  }
#endif

  bool ignore_reinsert = true; 

  if ( ignore_reinsert ) {
    tree_node_handle parent_handle = tree_node_handle(nullptr);
    if ( !parentHandles.empty() ) {
      parent_handle = parentHandles.top();
      parentHandles.pop();
    } 
    // Split the Node
    propagationSplit = current_node->splitNode(treeRef, current_handle, parent_handle);

    // Cleanup before ascending
    if ( parent_handle != nullptr) {
      // This will probably destroy current_node, so if we need
      // current node for anything, need to do it before the
      // removeEntry call.
      auto parent_node = treeRef->get_branch_node(parent_handle);
      parent_node->removeBranch(treeRef, current_handle);
    }

    // Ascend, propagating splits
    return std::make_pair(propagationSplit, parent_handle);
  } else {
    abort(); 
    // Nothing is real after you make this call
    // The reinsert might have come back around again and split this
    // node, or other nodes, or everyting
    // Signal to the caller that we shoudl stop
    current_node->reInsert(hasReinsertedOnLevel);
    return std::make_pair(propagationSplit, tree_node_handle(nullptr));  
  }
}


template <int min_branch_factor, int max_branch_factor, class strategy>
void fix_up_path_polys(
        tree_node_handle start_handle,
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
#endif 
  abort();
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
        NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
        tree_node_handle start_handle, 
        std::stack<tree_node_handle> &parentHandles,
        std::vector<bool> &hasReinsertedOnLevel)
{
  // N.B., as we walk up the tree, we may perform a bunch of splits,
  // which is liable to destroy nodes that are downsplit. These
  // downsplit nodes' memory can then be re-used for other things,
  // like polygons. If we try to use variables stored in that
  // node, it can be clobbered leading to amazing
  // segfaults. It is important that any variables we reference are
  // those we know are alive.

  assert(start_handle.get_type() == LEAF_NODE);
  tree_node_handle current_handle = start_handle;
  
  SplitResult propagationSplit = {
          {Rectangle(), tree_node_handle(nullptr)},
          {Rectangle(), tree_node_handle(nullptr)}};
  
  // Loop from the bottom to the very top (root)
  while (current_handle != nullptr) {
    if (propagationSplit.leftBranch.child == nullptr){
      assert(propagationSplit.rightBranch.child == nullptr);
    }
    // If there was a split we were supposed to propagate
    if (propagationSplit.leftBranch.child != nullptr and propagationSplit.rightBranch.child != nullptr) {
      assert(current_handle.get_type() == BRANCH_NODE);
      auto current_branch_node = treeRef->get_branch_node(current_handle);
      
      if (propagationSplit.leftBranch.child.get_type() == LEAF_NODE){
        assert(propagationSplit.rightBranch.child.get_type() == LEAF_NODE);
        auto left_node = treeRef->get_leaf_node(propagationSplit.leftBranch.child);
        assert(left_node->cur_offset_ > 0);
        auto right_node = treeRef->get_leaf_node(propagationSplit.rightBranch.child);
        assert(right_node->cur_offset_ > 0);
      } else {
        assert(propagationSplit.rightBranch.child.get_type() == BRANCH_NODE);
        auto left_node = treeRef->get_branch_node(propagationSplit.leftBranch.child, false);
        assert(left_node->cur_offset_ > 0);
        auto right_node = treeRef->get_branch_node(propagationSplit.rightBranch.child, false);
        assert(right_node->cur_offset_ > 0);
      }
      // Add splitted branches to Parent node
      current_branch_node->addBranchToNode(propagationSplit.leftBranch);
      current_branch_node->addBranchToNode(propagationSplit.rightBranch);
    }

    std::pair<SplitResult, tree_node_handle> split_res_and_new_handle;
    if (current_handle.get_type() == LEAF_NODE ) {
      auto current_leaf_node = treeRef->get_leaf_node(current_handle, true);
      split_res_and_new_handle = adjust_tree_bottom_half(
                                  current_leaf_node, 
                                  treeRef,
                                  current_handle,
                                  parentHandles,
                                  hasReinsertedOnLevel,
                                  max_branch_factor);
    } else {
      auto current_branch_node = treeRef->get_branch_node(current_handle);
      split_res_and_new_handle = adjust_tree_bottom_half(
                                  current_branch_node,
                                  treeRef,
                                  current_handle,
                                  parentHandles,
                                  hasReinsertedOnLevel,
                                  max_branch_factor);
    }
    propagationSplit = split_res_and_new_handle.first;
    current_handle = split_res_and_new_handle.second;
  }
  return propagationSplit;
}

// This always get called on the root node. So if it got called on us,
// that's because we are the only node in the whole tree.
template <int min_branch_factor, int max_branch_factor, class strategy>
tree_node_handle LeafNode<min_branch_factor, max_branch_factor, strategy>::insert(
        NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
        tree_node_handle selfHandle,
        Point givenPoint,
        std::vector<bool> &hasReinsertedOnLevel)
{
  // This is a leaf, so we are the ONLY node.
  addPoint(givenPoint);
  // Empty parentHandles as we are the ONLY node. 
  std::stack<tree_node_handle> parentHandles; 
  // adjustTreeSub splits node if necessary 
  SplitResult finalSplit = adjustTreeSub(treeRef, 
                                         selfHandle, 
                                         parentHandles,
                                         hasReinsertedOnLevel);

  // Grow the tree taller if there are two branches at root level 
  if (finalSplit.leftBranch.child != nullptr and finalSplit.rightBranch.child != nullptr) {
    tree_node_allocator *allocator = get_node_allocator(treeRef);
    auto alloc_data = allocator->create_new_tree_node<BranchNode<min_branch_factor, max_branch_factor, strategy>>(
                      NodeHandleType(BRANCH_NODE));
    new (&(*alloc_data.first)) BranchNode<min_branch_factor, max_branch_factor, strategy>();
    auto new_root_handle = alloc_data.second;
    auto new_root_node = alloc_data.first;
    assert(new_root_handle.get_type() == BRANCH_NODE);

    // It is only possible for left and right tree node to be leaf 
    // as we are starting as one leaf node
    assert(finalSplit.leftBranch.child.get_type() == LEAF_NODE);
    assert(finalSplit.rightBranch.child.get_type() == LEAF_NODE);

    // Add to new root
    new_root_node->addBranchToNode(finalSplit.leftBranch);
    new_root_node->addBranchToNode(finalSplit.rightBranch);

    // update root in tree 
    treeRef->root = new_root_handle;
  
    // free the Leaf Node and remove the associated polygon from map 
    assert(selfHandle.get_type() == LEAF_NODE);
    allocator->free(selfHandle, sizeof(LeafNode<min_branch_factor, max_branch_factor, strategy>));
    remove_polygon(treeRef, selfHandle);

    // Fix the reinserted length
    hasReinsertedOnLevel.push_back(false);

    return new_root_handle;
  }

  // no need to grow the tree taller and the given Point is added to entries vector 
  return treeRef->root;
}

// This is called on a Leaf Node to shorten the height of tree
template <int min_branch_factor, int max_branch_factor, class strategy>
void LeafNode<min_branch_factor, max_branch_factor, strategy>::condenseTree(
        NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
        tree_node_handle selfHandle,
        std::stack<tree_node_handle> &parentHandles) {
  // Quick return as the current Leaf Node is not empty 
  if (this->cur_offset_ != 0) return; 
  // Working on root node which is Leaf Node 
  if (parentHandles.empty()) return; 

  tree_node_handle current_node_handle = selfHandle; 
  tree_node_handle parent_node_handle = parentHandles.top();
  parentHandles.pop();

  // remove current LeafNode Branch from parent 
  auto parent_node = treeRef->get_branch_node(parent_node_handle);
  assert(this->cur_offset_ == 0); 
  parent_node->removeBranch(treeRef, current_node_handle);

  // Loop from bottom to top 
  auto current_node = parent_node; 
  while (not parentHandles.empty()) {
    current_node = parent_node; 
    current_node_handle = parent_node_handle; 
    // no further condense work needed
    if (current_node->cur_offset_ != 0) break; 
    parent_node_handle = parentHandles.top();
    parentHandles.pop();
    // remove current LeafNode Branch from parent 
    parent_node = treeRef->get_branch_node(parent_node_handle);
    assert(current_node->cur_offset_ == 0); 
    parent_node->removeBranch(treeRef, current_node_handle);
  }
}

// Always called on root, this = root
template <int min_branch_factor, int max_branch_factor, class strategy>
tree_node_handle LeafNode<min_branch_factor, max_branch_factor, strategy>::remove(
        NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
        tree_node_handle selfHandle,
        Point givenPoint) 
{
  removePoint(givenPoint);
  return selfHandle;
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
void LeafNode<min_branch_factor, max_branch_factor, strategy>::print(tree_node_handle current_handle,
                                                                     tree_node_handle parent_handle,
                                                                     unsigned n) {
  std::string indentation(n * 4, ' ');
  std::cout << indentation << "Node " << (void *)this << std::endl;
  std::cout << indentation << "    Parent: " << parent_handle << std::endl;
  std::cout << indentation << "    Current: " << current_handle << std::endl;
  std::cout << indentation << "    Data: ";
  for (unsigned i = 0; i < this->cur_offset_; i++) {
    std::cout << entries.at(i);
  }
  std::cout << std::endl;
}

template <int min_branch_factor, int max_branch_factor, class strategy>
void LeafNode<min_branch_factor, max_branch_factor, strategy>::printTree(tree_node_handle current_handle,
                                                                         tree_node_handle parent_handle, 
                                                                         unsigned n) {
  // Print this node first
  print(current_handle, parent_handle, n);
  std::cout << std::endl;
}

template <int min_branch_factor, int max_branch_factor, class strategy>
unsigned LeafNode<min_branch_factor, max_branch_factor, strategy>::height() {
  return 1;
}

template <int min_branch_factor, int max_branch_factor, class strategy>
void BranchNode<min_branch_factor, max_branch_factor, strategy>::deleteSubtrees(
        NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef) {
  tree_node_allocator *allocator = get_node_allocator(treeRef);
  for (unsigned i = 0; i < this->cur_offset_; i++) {
    Branch &b = entries.at(i);
    tree_node_handle child_handle = b.child;

    if (child_handle.get_type() == LEAF_NODE) {
      allocator->free(
              child_handle,
              sizeof(LeafNode<min_branch_factor, max_branch_factor, strategy>));
    } else if (child_handle.get_type() == BRANCH_NODE) {
      auto child = treeRef->get_branch_node(child_handle);
      child->deleteSubtrees(treeRef);
      allocator->free(
              child_handle,
              sizeof(BranchNode<min_branch_factor, max_branch_factor, strategy>));
    }
    // remove polygon from map
    remove_polygon(treeRef, child_handle);
    b.child = tree_node_handle(nullptr);
  }
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


// Removes a child logical pointer from a this->parent node, freeing that
// child's memory and the memory of the associated polygon 
template <int min_branch_factor, int max_branch_factor, class strategy>
void BranchNode<min_branch_factor, max_branch_factor, strategy>::removeBranch(
          NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
          const tree_node_handle entry)
{
  // Locate the child
  unsigned childIndex = 0;
  while( entries.at(childIndex).child != entry and childIndex < this->cur_offset_ )
  { 
    childIndex++; 
  }
  assert(entries.at(childIndex).child == entry);

  // Free the branch/leaf Node 
  tree_node_allocator *allocator = get_node_allocator(treeRef);
  Branch &b = entries.at(childIndex);
  if (b.child.get_type() == LEAF_NODE) {
    allocator->free(b.child, sizeof(LeafNode<min_branch_factor, max_branch_factor, strategy>));
  } else if (b.child.get_type() == BRANCH_NODE) {
    allocator->free(b.child, sizeof(BranchNode<min_branch_factor, max_branch_factor, strategy>));
  }

  // remove the polygon from map 
  remove_polygon(treeRef, b.child);
  b.child = tree_node_handle(nullptr);

  // swap with last entry in array 
  entries.at(childIndex) = entries.at(this->cur_offset_ - 1);
  // Truncate array size
  this->cur_offset_--;
}

// ??? what is this function used for ???
template <int min_branch_factor, int max_branch_factor, class strategy, typename functor>
void is_vertical_stripe(NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef, tree_node_handle root, functor &f) {
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

// ??? what is this function used for ???
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

// ======================================================================================
// functions for searching methods:
// interface
//    point_search()
//    rectangle_search()
// search helper functions: 
//    point_search_leaf_node()
//    point_search_branch_node()
//    rectangle_search_leaf_node()
//    rectangle_search_branch_node()
// not used:
//    parent_handle_point_search 
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
// == SEARCH ====================================================================================


// Always called on root, this = root
// This top-to-bottom sweep is only for adjusting bounding boxes to contain the point and
// choosing a particular leaf
template <int min_branch_factor, int max_branch_factor, class strategy>
tree_node_handle
BranchNode<min_branch_factor, max_branch_factor, strategy>::chooseNodePoint(
          NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
          tree_node_handle selfHandle,
          std::stack<tree_node_handle> &parentHandles,
          Point &point) {
  
  // CN1 [Initialize]
  assert(selfHandle != nullptr);
  assert(treeRef->root == selfHandle); 
  tree_node_handle current_handle = selfHandle; 
  // Starting root, searching for a Leaf node for insertion
  for (;;) {
    assert(current_handle != nullptr);
    if (current_handle.get_type() == LEAF_NODE) {
      // Found Leaf Node, done
      return current_handle;
    } 
    assert(current_handle.get_type() == BRANCH_NODE);
    auto current_node = treeRef->get_branch_node(current_handle);
    assert(current_node->cur_offset_ > 0);

    // CN2 [initialize node search]
    // This is the minimum amount of additional area we need in
    // one of the branches to encode our expansion
    double minimal_area_expansion = std::numeric_limits<double>::max();
    // This is the list of optimal expansiosn we need to perform
    // to get the bounding box/bounding polygon
    std::vector<IsotheticPolygon::OptimalExpansion> expansions;

    // Find polygon of the first branch at current Branch Node
    IsotheticPolygon node_poly = find_polygon(treeRef, current_node->entries.at(0));
    IsotheticPolygon::OptimalExpansion exp = node_poly.computeExpansionArea(point);
    minimal_area_expansion = exp.area;
    expansions.push_back(exp);
    unsigned smallestExpansionBranchIndex = 0;

    // Checking other branch/child starting at index 1 for better option
    for (unsigned i = 1; i < current_node->cur_offset_; i++) {
      Branch &b = current_node->entries.at(i);
      node_poly = find_polygon(treeRef, b); 
      exp = node_poly.computeExpansionArea(point);
      if (exp.area < minimal_area_expansion) {
        minimal_area_expansion = exp.area;
        expansions.clear();
        expansions.push_back(exp);
        smallestExpansionBranchIndex = i;
      }
    }

    // The point has not been contained in a polygon 
    // expansion has to be made on some polygon 
    if (minimal_area_expansion != -1.0) {
#ifndef NDEBUG
      for (unsigned i = 0; i < current_node->cur_offset_; i++) {
        for (unsigned j = 0; j < current_node->cur_offset_; j++) {
          if (i == j) {
            continue;
          }
          Branch &b_i = current_node->entries.at(i);
          Branch &b_j = current_node->entries.at(j);

          IsotheticPolygon poly_i = find_polygon(treeRef, b_i); 
          IsotheticPolygon poly_j = find_polygon(treeRef, b_j);
          assert(poly_i.disjoint(poly_j));
        }
      }
#endif

      // expand chosen branch to include point
      Branch &chosen_branch = current_node->entries.at(smallestExpansionBranchIndex);
      IsotheticPolygon chosen_poly = find_polygon(treeRef, chosen_branch); 
      assert(expansions.size() == 1);
      Rectangle &existing_rect = chosen_poly.basicRectangles.at(expansions.at(0).index);
      existing_rect.expand(point);
      chosen_poly.recomputeBoundingBox();
      assert(chosen_poly.containsPoint(point));
        
      // Dodge all the other branches
      for (unsigned i = 0; i < current_node->cur_offset_; i++) {
        if (i == smallestExpansionBranchIndex) {
          continue;
        }
        Branch &other_branch = current_node->entries.at(i);
        IsotheticPolygon other_poly = find_polygon(treeRef, other_branch); 
        chosen_poly.increaseResolution(Point::atInfinity, other_poly);
      }
      chosen_poly.refine();
      chosen_poly.recomputeBoundingBox();
      assert(chosen_poly.containsPoint(point));
      update_branch_polygon(treeRef, chosen_branch, chosen_poly);
    }

    // Descend
    parentHandles.push(current_handle);
    Branch &b = current_node->entries.at(smallestExpansionBranchIndex);
    current_handle = b.child;
    assert(current_handle != nullptr);
  } // for 
} // chooseNodePoint

// Always called on root, this = root
// This top-to-bottom sweep is only for adjusting bounding boxes to contain the branch and
// at a particular level
// [TODO]
template <int min_branch_factor, int max_branch_factor, class strategy>
tree_node_handle
BranchNode<min_branch_factor, max_branch_factor, strategy>::chooseNodeBranch(
          NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
          tree_node_handle selfHandle,
          std::stack<tree_node_handle> &parentHandles,
          BranchAtLevel &branchLevel) 
{
  
  // CL1 [Initialize]
  assert(selfHandle != nullptr);
  // insert() should always be called at the root level
  assert(treeRef->root == selfHandle); 
  tree_node_handle current_handle = selfHandle; 
  // get tree height to calculate branch node level 
  uint8_t height = this->height(treeRef,selfHandle);
  // stop at parent level 
  uint8_t stopping_level = branchLevel.level + 1; 
  Branch &branch = branchLevel.branch; 
  // polygon of to-be-inserted branch 
  IsotheticPolygon branch_poly = find_polygon(treeRef, branch); 
  for (;;) {
    assert(current_handle != nullptr);
    assert(current_handle.get_type() != LEAF_NODE);
    // if we are at the stopping level, it means we have found the
    // Branch Node 
    uint8_t current_level = height - parentHandles.size(); 
    if (current_level == stopping_level) {
      return current_handle;
    }
    auto current_node = treeRef->get_branch_node(current_handle);
    assert(current_node->cur_offset_ > 0);

    // CN2 [initialize node search]
    // This is the minimum amount of additional area we need in
    // one of the branches to encode our expansion
    double minimal_area_expansion = std::numeric_limits<double>::max();
    double minimal_poly_area = std::numeric_limits<double>::max();
    // This is the list of optimal expansiosn we need to perform
    // to get the bounding box/bounding polygon
    std::vector<IsotheticPolygon::OptimalExpansion> expansions;
    // find polygon of the first branch at current Branch Node
    IsotheticPolygon node_poly = find_polygon(treeRef, current_node->entries.at(0));
    auto exp = computeExpansionArea(node_poly, branch_poly); 
    minimal_area_expansion = exp.first;
    expansions = exp.second;
    minimal_poly_area = branch_poly.area(); 
    // This is the branch that gives us that minimum area expansion
    unsigned smallestExpansionBranchIndex = 0;   
    // we are expecting expansion for branch insertion 
    assert(minimal_area_expansion == -1.0);

    // checking each branch/child of current Branch Node 
    for (unsigned i = 1; i < current_node->cur_offset_; i++) {
      Branch &b = current_node->entries.at(i);
      node_poly = find_polygon(treeRef, b); 
        // Walk every rectangle in the branch's polygon
        // Find rectangle in our polygon that needs to be
        // expanded the least to fit the branch's rectangle
        // inside it.
        // N.B., this does not split the rectangle apart if
        // the expanded rectangle could be part of two
        // distinct polygons. So as a result of doing this
        // the polygon's constituent rectangles may now
        // overlap.
      exp = computeExpansionArea(node_poly, branch_poly); 
      double poly_area = node_poly.area(); 
      if (exp.first < minimal_area_expansion or
          (exp.first == minimal_area_expansion and poly_area < minimal_poly_area)) {
        minimal_area_expansion = exp.first;
        minimal_poly_area = poly_area;
        expansions = exp.second;
        smallestExpansionBranchIndex = i;
      }
    }

    if (minimal_area_expansion != -1.0) {
#ifndef NDEBUG
      for (unsigned i = 0; i < current_node->cur_offset_; i++) {
        for (unsigned j = 0; j < current_node->cur_offset_; j++) {
          if (i == j) {
            continue;
          }
          Branch &b_i = current_node->entries.at(i);
          Branch &b_j = current_node->entries.at(j);
          IsotheticPolygon poly_i = find_polygon(treeRef, b_i); 
          IsotheticPolygon poly_j = find_polygon(treeRef, b_j); 
          assert(poly_i.disjoint(poly_j));
        }
      }
#endif

      Branch &chosen_branch = current_node->entries.at(smallestExpansionBranchIndex);
      IsotheticPolygon chosen_poly = find_polygon(treeRef, chosen_branch); 
      assert(branch_poly.basicRectangles.size() == expansions.size());
      // Fragment them on the way down.
      // clipping other poly according to the inserted branch 
      // why not the other way around ??? 
      for (unsigned i = 0; i < current_node->cur_offset_; i++) {
        if (i == smallestExpansionBranchIndex) {
          continue;
        }
        Branch &other_branch = current_node->entries.at(i);
        IsotheticPolygon other_poly = find_polygon(treeRef, other_branch); 
        other_poly.increaseResolution(Point::atInfinity, branch_poly);
        other_poly.refine();
        other_poly.recomputeBoundingBox();
        update_branch_polygon(treeRef, other_branch, other_poly);
      }
      // We need to expand on way down so we know who is
      // responsible for the new point/branch
      // Everyone else needs to fragment around my nodeEntry,
      // then we expand and fragment around them.
      for (unsigned i = 0; i < branch_poly.basicRectangles.size(); i++) {
        auto &expansion = expansions.at(i);
        auto &branch_rect = branch_poly.basicRectangles.at(i);
        Rectangle &existing_rect = chosen_poly.basicRectangles.at(expansion.index);
        // Expand the existing rectangle. This rectangle
        // might now overlap with other rectangles in
        // the polygon. But if we make it not overlap,
        // then we alter the indices of the expansion
        // rectangles, which kind of sucks, So, leave it
        // for now.
        existing_rect.expand(branch_rect);
        assert(existing_rect.containsRectangle(branch_rect));
      }
      chosen_poly.recomputeBoundingBox();

      // Dodge all the other branches
      for (unsigned i = 0; i < current_node->cur_offset_; i++) {
        if (i == smallestExpansionBranchIndex) {
          continue;
        }
        Branch &other_branch = current_node->entries.at(i);
        IsotheticPolygon other_poly = find_polygon(treeRef, other_branch); 
        chosen_poly.increaseResolution(Point::atInfinity, other_poly);
      }
      chosen_poly.refine();
      chosen_poly.recomputeBoundingBox();

      update_branch_polygon(treeRef, chosen_branch, chosen_poly);
    }

    // Descend
    parentHandles.push(current_handle);
    Branch &b = current_node->entries.at(smallestExpansionBranchIndex);
    current_handle = b.child;
    assert(current_handle != nullptr);
  
  } // for 
} // chooseNodeBranch


// Find which Leaf Node contains the Point or nullptr if none 
template <int min_branch_factor, int max_branch_factor, class strategy>
tree_node_handle BranchNode<min_branch_factor, max_branch_factor, strategy>::findLeaf(
     NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
     tree_node_handle selfHandle,
     Point givenPoint) {

  // Initialize our context stack
  std::stack<tree_node_handle> context;
  context.push(selfHandle);
  tree_node_handle current_node_handle;
  tree_node_allocator *allocator = get_node_allocator(treeRef);

  // Traverse through the tree with DFS
  while (!context.empty()) {
    current_node_handle = context.top();
    context.pop();

    if (current_node_handle.get_type() == LEAF_NODE) {
      auto current_node = treeRef->get_leaf_node(current_node_handle);
      // Check each entry to see if it matches p
      for (unsigned i = 0; i < current_node->cur_offset_; i++) {
        Point &p = current_node->entries.at(i);
        if (p == givenPoint) {
          return current_node_handle;
        }
      }
    } else {
      auto current_node = treeRef->get_branch_node(current_node_handle);
      // Determine which branches we need to follow
      for (unsigned i = 0; i < current_node->cur_offset_; i++) {
        Branch &b = current_node->entries.at(i);
        // Quick Check 
        if (!b.boundingBox.containsPoint(givenPoint)) {
          // if point is not contained in MBB, skip
          continue; 
        }
        IsotheticPolygon poly  = find_polygon(treeRef, b); 

        if (poly.containsPoint(givenPoint)) {
          // Add the child to the nodes we will consider
          context.push(b.child);
        }
      }
    }
  }

  return tree_node_handle(nullptr);
}

template <int min_branch_factor, int max_branch_factor, class strategy>
Partition BranchNode<min_branch_factor, max_branch_factor, strategy>::partitionNode(NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef) {
  return partitionBranchNode(treeRef);
}

// ??? what is this structure used for ???
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
          NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
          tree_node_handle handle_to_skip,
          IsotheticPolygon &polygon) {
  tree_node_allocator *allocator = get_node_allocator(treeRef);
  for (auto iter = entries.begin(); iter != entries.begin() + this->cur_offset_; iter++) {
    Branch &b = *iter;
    if (b.child == handle_to_skip) {
      continue;
    }
    // Get each siblings' polygon from map 
    IsotheticPolygon child_poly = find_polygon(treeRef, b); 
    // Polygon is clipped according to child_poly
    polygon.increaseResolution(Point::atInfinity, child_poly);
  }
  // refine and recomputeBoundingBox is usually called by caller 
  // polygon.refine();
  // polygon.recomputeBoundingBox();

}

// We create two new nodes and free the old one.
// The old one is freed in adjust_tree_bottom_half using removeBranch
// If we downsplit, then we won't call adjustTree for that split so we
// need to delete the node ourselves.
template <int min_branch_factor, int max_branch_factor, class strategy>
SplitResult BranchNode<min_branch_factor, max_branch_factor, strategy>::splitNode(
          NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
          tree_node_handle current_handle,
          tree_node_handle parent_handle,
          Partition p, 
          bool is_downsplit) {

  assert(current_handle.get_type() == BRANCH_NODE);
  tree_node_allocator *allocator = get_node_allocator(treeRef);
  
  // Allocate a branch node for left node
  auto alloc_data = allocator->create_new_tree_node<BranchNode<min_branch_factor, max_branch_factor, strategy>>(
                  NodeHandleType(BRANCH_NODE));
  tree_node_handle left_handle = alloc_data.second;
  auto left_node = alloc_data.first; // take pin
  new (&(*left_node)) BranchNode<min_branch_factor, max_branch_factor, strategy>();
  
  // Allocate a branch node for right node
  alloc_data = allocator->create_new_tree_node<BranchNode<min_branch_factor, max_branch_factor, strategy>>(
                  NodeHandleType(BRANCH_NODE));
  tree_node_handle right_handle = alloc_data.second;
  auto right_node = alloc_data.first; // take pin
  new (&(*right_node)) BranchNode<min_branch_factor, max_branch_factor, strategy>();

  SplitResult split = {{Rectangle(), left_handle},
                       {Rectangle(), right_handle}};
  
  // So we are going to split all branches at this branch node.
  for (unsigned i = 0; i < this->cur_offset_; i++) {
    Branch &branch = entries.at(i);
    Rectangle summary_rectangle = branch.boundingBox;

    bool is_contained_left = summary_rectangle.upperRight[p.dimension] <= p.location;
    bool is_contained_right = summary_rectangle.lowerLeft[p.dimension] >= p.location;
    assert( not(is_contained_left and is_contained_right));

    if (is_contained_left and not is_contained_right) {
      // Entirely contained in the left polygon
      assert(split.leftBranch.child == left_handle);
      left_node->addBranchToNode(branch);
    } else if (is_contained_right and not is_contained_left) {
      // Entirely contained in the right polygon
      assert(split.rightBranch.child == right_handle);
      right_node->addBranchToNode(branch);
    
    } else if (summary_rectangle.upperRight[p.dimension] ==
               summary_rectangle.lowerLeft[p.dimension] and
               summary_rectangle.lowerLeft[p.dimension] ==
               p.location) {
      // These go left or right situationally
      if (left_node->cur_offset_ <= right_node->cur_offset_) {
        assert(split.leftBranch.child == left_handle);
        left_node->addBranchToNode(branch);
      } else {
        assert(split.rightBranch.child == right_handle);
        right_node->addBranchToNode(branch);
      }
    } else {
      // Partially spanned by both nodes, need to downsplit
      IsotheticPolygon branch_poly = find_polygon(treeRef, branch); 
      SplitResult downwardSplit;

      // Downward Split
      if (branch.child.get_type() == LEAF_NODE) {
        auto child = treeRef->get_leaf_node(branch.child);
        downwardSplit = child->splitNode(treeRef, branch.child, current_handle, p, true);
        allocator->free(branch.child, sizeof(LeafNode<min_branch_factor, max_branch_factor, strategy>));
      } else {
        auto child = treeRef->get_branch_node(branch.child);
        downwardSplit = child->splitNode(treeRef, branch.child, current_handle, p, true);
        allocator->free(branch.child, sizeof(BranchNode<min_branch_factor, max_branch_factor, strategy>));
      }

      // DONT'T call removeBranch() which would mess up the index
      remove_polygon(treeRef, branch.child);
      branch.child = tree_node_handle(nullptr);

      // Add downward splitted branches to left and right nodes 
      if (downwardSplit.leftBranch.child.get_type() == LEAF_NODE) {
        auto left_child = treeRef->get_leaf_node(downwardSplit.leftBranch.child);
        if (left_child->cur_offset_ > 0) {
          left_node->addBranchToNode(downwardSplit.leftBranch);
        }
        assert(downwardSplit.rightBranch.child.get_type() == LEAF_NODE);
        auto right_child = treeRef->get_leaf_node(downwardSplit.rightBranch.child);
        if (right_child->cur_offset_ > 0) {
          right_node->addBranchToNode(downwardSplit.rightBranch);
        }
      } else {
        auto left_child = treeRef->get_branch_node(downwardSplit.leftBranch.child);
        if (left_child->cur_offset_ > 0) {
          left_node->addBranchToNode(downwardSplit.leftBranch);
        }
        assert(downwardSplit.rightBranch.child.get_type() == BRANCH_NODE);
        auto right_child = treeRef->get_branch_node(downwardSplit.rightBranch.child);
        if (right_child->cur_offset_ > 0) {
          right_node->addBranchToNode(downwardSplit.rightBranch);
        }
      }

#ifdef DEBUG_TEST
  IsotheticPolygon left_polygon = find_polygon(treeRef, downwardSplit.leftBranch);
  IsotheticPolygon right_polygon = find_polygon(treeRef, downwardSplit.rightBranch);
  assert(left_polygon.disjoint(right_polygon));
#endif

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
  // [TODO]: [FIXME]

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
  if (parent_handle) {
    auto parent_node = treeRef->get_branch_node(parent_handle);
    if (not is_downsplit) {
      // Left node
      assert(left_polygon.basicRectangles.size() > 0);
      parent_node->make_disjoint_from_children(treeRef,
                                               current_handle,
                                               left_polygon);
      assert(left_polygon.basicRectangles.size() > 0);
      left_polygon.refine();
      assert(left_polygon.basicRectangles.size() > 0);
      // do we need recomputeBoundingBox ???
      left_polygon.recomputeBoundingBox();

      // Right node
      assert(right_polygon.basicRectangles.size() > 0);
      parent_node->make_disjoint_from_children(treeRef,
                                               current_handle,
                                               right_polygon);
      assert(right_polygon.basicRectangles.size() > 0);
      right_polygon.refine();
      assert(right_polygon.basicRectangles.size() > 0);
      right_polygon.recomputeBoundingBox();
      
      assert(left_polygon.disjoint(right_polygon));
    } else { 
      // Intersect with our existing poly to avoid intersect
      // with other children
      IsotheticPolygon current_polygon = find_polygon(treeRef, current_handle, 
                                          this->boundingBox());
      assert(current_polygon.basicRectangles.size() > 0);
      
      // Left node: 
      assert(left_polygon.basicRectangles.size() > 0);
      left_polygon.intersection(current_polygon);
      assert(left_polygon.basicRectangles.size() > 0);
      left_polygon.refine();
      assert(left_polygon.basicRectangles.size() > 0);
      left_polygon.recomputeBoundingBox();
      
      // Right node: 
      right_polygon.intersection(current_polygon);
      assert(right_polygon.basicRectangles.size() > 0);
      right_polygon.refine();
      assert(right_polygon.basicRectangles.size() > 0);
      right_polygon.recomputeBoundingBox();
      assert(left_polygon.disjoint(right_polygon));
    }
  }
  // This BranchNode will be freed by caller 
  this->cur_offset_ = 0;
  update_branch_polygon(treeRef, split.leftBranch, left_polygon);
  update_branch_polygon(treeRef, split.rightBranch, right_polygon);

#ifdef DEBUG_TEST
  testDisjoint(treeRef, split.leftBranch.child, "left subtree after split");
  testDisjoint(treeRef, split.rightBranch.child, "right subtree after split");
#endif

  return split;

}


template <int min_branch_factor, int max_branch_factor, class strategy>
SplitResult BranchNode<min_branch_factor, max_branch_factor, strategy>::splitNode(
                        NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
                        tree_node_handle current_handle,
                        tree_node_handle parent_handle) {
  SplitResult returnSplit = splitNode(treeRef, current_handle, parent_handle, partitionNode(treeRef), false);
  return returnSplit;
}


// insert() is always called on root node 
template <int min_branch_factor, int max_branch_factor, class strategy>
tree_node_handle BranchNode<min_branch_factor, max_branch_factor, strategy>::insert(
        NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
        tree_node_handle selfHandle,
        std::variant<BranchAtLevel, Point> &nodeEntry, 
        std::vector<bool> &hasReinsertedOnLevel)
{
#ifdef DEBUG_TEST
  testDisjoint(treeRef, treeRef->root,"before insertion");
#endif

  bool givenIsPoint = std::holds_alternative<Point>(nodeEntry);
  assert(selfHandle.get_type() == BRANCH_NODE);
  
  // parentHandles are collected in chooseNodePoint 
  std::stack<tree_node_handle> parentHandles; 
  tree_node_handle current_handle;

  // choose the best node for insertion
  if (givenIsPoint){
    current_handle = chooseNodePoint(treeRef, 
                                     selfHandle, 
                                     parentHandles, 
                                     std::get<Point>(nodeEntry));
  } else {
    current_handle = chooseNodeBranch(treeRef, 
                                     selfHandle, 
                                     parentHandles, 
                                     std::get<BranchAtLevel>(nodeEntry));
  }

#ifdef DEBUG_TEST
  testDisjoint(treeRef, treeRef->root, "after chooseNode");
#endif

  tree_node_allocator *allocator = get_node_allocator(treeRef);
  SplitResult finalSplit;

  if ( givenIsPoint ) {
    assert(current_handle.get_type() == LEAF_NODE);
    auto current_node = treeRef->get_leaf_node(current_handle);
    
    // Add point to chosen node 
    current_node->addPoint(std::get<Point>(nodeEntry));

    // Split if needed 
    finalSplit = adjustTreeSub(treeRef,
                               current_handle,
                               parentHandles,
                               hasReinsertedOnLevel);
  } else {
    // working with inserting Branch to Branch Node 
    assert(current_handle.get_type() == BRANCH_NODE);
    auto current_node = treeRef->get_branch_node(current_handle);

    BranchAtLevel &sub_bl = std::get<BranchAtLevel>(nodeEntry);
    Branch &sub_branch = sub_bl.branch; 
    IsotheticPolygon insertion_polygon  = find_polygon(treeRef, sub_branch); 
    
    // Before I add this node in, I need to fragment everyone else
    // around it
    for (unsigned int i = 0; i < current_node->cur_offset_; i++) {
      Branch &b = current_node->entries.at(i);
      IsotheticPolygon branch_polygon = find_polygon(treeRef, b); 
      branch_polygon.increaseResolution(Point::atInfinity, insertion_polygon);
      branch_polygon.refine(); 
      branch_polygon.recomputeBoundingBox();
      update_branch_polygon(treeRef, b, branch_polygon);
    }

    // Add branch to chosen node 
    current_node->addBranchToNode(sub_branch);

    // Split if needed  
    finalSplit = adjustTreeSub(treeRef,
                               current_handle,
                               parentHandles,
                               hasReinsertedOnLevel);
  }

#ifdef DEBUG_TEST
  if (finalSplit.leftBranch.child == nullptr ) {
    testDisjoint(treeRef, treeRef->root, "after adjustTreeSub");
  }
#endif

  tree_node_handle ret_handle; 
  // Grow the tree taller if we need to
  if (finalSplit.leftBranch.child != nullptr and finalSplit.rightBranch.child != nullptr) {

#ifdef DEBUG_TEST
  testDisjoint(treeRef,finalSplit.rightBranch.child, "Right Subtree");
  testDisjoint(treeRef, finalSplit.leftBranch.child, "Left Subtree");
#endif
    // Allocate new root
    auto alloc_data = allocator->create_new_tree_node<BranchNode<min_branch_factor, max_branch_factor, strategy>>(
                    NodeHandleType(BRANCH_NODE));
    new (&(*alloc_data.first)) BranchNode<min_branch_factor, max_branch_factor, strategy>();
    auto new_root_handle = alloc_data.second;
    auto new_root_node = alloc_data.first;

    // Add to new root
    new_root_node->addBranchToNode(finalSplit.leftBranch);
    new_root_node->addBranchToNode(finalSplit.rightBranch);
    
    // update root in the tree
    treeRef->root = new_root_handle;

    assert(selfHandle.get_type() == BRANCH_NODE);
    allocator->free(selfHandle, sizeof(BranchNode<min_branch_factor, max_branch_factor, strategy>));
    remove_polygon(treeRef, selfHandle);
    selfHandle = tree_node_handle(nullptr);

    // Fix the reinserted length
    hasReinsertedOnLevel.push_back(false);
    ret_handle = new_root_handle;
#ifdef DEBUG_TEST
  std::cout << "created new root " << std::endl; 
  testDisjoint(treeRef, ret_handle, "after creating new root");
#endif 
  } else {
    // no need to grow the tree taller and the given Point is added to entries vector 
    ret_handle = treeRef->root;
  }

#ifdef DEBUG_TEST
  std::cout << "inserted " << std::get<Point>(nodeEntry) << std::endl; 
  //testDisjoint(treeRef, ret_handle);
  testCount(treeRef, ret_handle);
  testContainPoints(treeRef, ret_handle);
#endif 

  return ret_handle;
}

// Always called on root, this = root
template <int min_branch_factor, int max_branch_factor, class strategy>
tree_node_handle BranchNode<min_branch_factor, max_branch_factor, strategy>::remove(        
        NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
        tree_node_handle selfHandle,
        Point givenPoint) {

  // D1 [Locate record] 
  tree_node_handle leaf_handle = findLeaf(treeRef, selfHandle, givenPoint);
  // Record not in the tree ; done
  // ??? is this expected ??? 
  if (leaf_handle == nullptr) {
    return selfHandle;
  }

  // D2 [Delete record]
  auto leaf_node = treeRef->get_leaf_node(leaf_handle);
  leaf_node->removePoint(givenPoint);

  // D3 [Propagate changes]
  std::stack<tree_node_handle> path_handles = find_path_to_leaf(treeRef, selfHandle, leaf_handle); 
  assert(path_handles.top() == leaf_handle);
  // Pop the Leaf Node to make it just parent handles 
  path_handles.pop();
  leaf_node->condenseTree(treeRef, leaf_handle, path_handles);

  // D4 [Shorten tree]
  assert(selfHandle.get_type() == BRANCH_NODE);
  // There is an entry left in root, let the only entry be the new root
  if (this->cur_offset_ == 1) {
    tree_node_handle new_root_handle = entries.at(0).child;
    tree_node_allocator *allocator = get_node_allocator(treeRef);
    allocator->free(selfHandle, sizeof(BranchNode<min_branch_factor, max_branch_factor, strategy>));
    
    // remove the polygon associated with old root from map 
    remove_polygon(treeRef, selfHandle);
    return new_root_handle;
  }

  // no work to do
  return selfHandle;
}

template <int min_branch_factor, int max_branch_factor, class strategy>
unsigned BranchNode<min_branch_factor, max_branch_factor, strategy>::checksum(
        NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef) {
  unsigned sum = 0;

  for (unsigned i = 0; i < this->cur_offset_; i++) {
    // Recurse
    Branch &branch = entries.at(i);
    if (branch.child.get_type() == LEAF_NODE || branch.child.get_type() == REPACKED_LEAF_NODE) {
      auto child = treeRef->get_leaf_node(branch.child);
      sum += child->checksum();
    } else {
      auto child = treeRef->get_branch_node(branch.child);
      sum += child->checksum(treeRef);
    }
  }

  return sum;
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
void BranchNode<min_branch_factor, max_branch_factor, strategy>::print(NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
                                                                       tree_node_handle current_handle,
                                                                       tree_node_handle parent_handle, 
                                                                       unsigned n) {
  std::string indentation(n * 4, ' ');
  std::cout << indentation << "Node " << (void *)this << std::endl;
  std::cout << indentation << "    Parent: " << parent_handle << std::endl;  
  std::cout << indentation << "    Current: " << current_handle << std::endl;
  std::cout << indentation << "    Branches: " << std::endl;
  for (unsigned i = 0; i < this->cur_offset_; i++) {
    Branch &branch = entries.at(i);
    auto poly = find_polygon(treeRef, branch); 

    std::cout << indentation << "		" << branch.child << std::endl;
    std::cout << indentation << "		" << poly << std::endl;
  }
  std::cout << std::endl;

}

// if parent_handle == nullptr, then we know the current node is root 
template <int min_branch_factor, int max_branch_factor, class strategy>
void BranchNode<min_branch_factor, max_branch_factor, strategy>::printTree(NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
                                                                           tree_node_handle current_handle,
                                                                           tree_node_handle parent_handle, 
                                                                           unsigned n) {
  // Print this node first
  print(treeRef, current_handle, parent_handle, n);

  // Print any of our children with one more level of indentation
  std::string indendtation(n * 4, ' ');
  for (unsigned i = 0; i < this->cur_offset_; i++) {
    Branch &branch = entries.at(i);
    if (branch.child.get_type() == LEAF_NODE || branch.child.get_type() == REPACKED_LEAF_NODE) {
      auto child = treeRef->get_leaf_node(branch.child, false);

      // Recurse
      child->printTree(branch.child, current_handle, n + 1);
    } else {
      auto child = treeRef->get_branch_node(branch.child, false);

      // Recurse
      child->printTree(treeRef, branch.child, current_handle, n + 1);
    }
  }
  std::cout << std::endl;
}

template <int min_branch_factor, int max_branch_factor, class strategy>
unsigned BranchNode<min_branch_factor, max_branch_factor, strategy>::height(
                    NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
                    tree_node_handle selfHandle) {
  unsigned ret = 0;
  tree_node_handle current_handle = selfHandle;

  for (;;) {
    ret++;
    if (current_handle.get_type() == LEAF_NODE || current_handle.get_type() == REPACKED_LEAF_NODE) {
      return ret;
    }

    auto node = treeRef->get_branch_node(current_handle, false);
    current_handle = node->entries.at(0).child;
  }
}

// called by NIRTreeDisk<min_branch_factor, max_branch_factor, strategy>::stat()
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


template <int min_branch_factor, int max_branch_factor, class strategy>
IsotheticPolygon find_polygon(
                  NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
                  tree_node_handle node_handle,
                  Rectangle rectangle){
  assert(node_handle != nullptr);
  std::map<tree_node_handle, IsotheticPolygon>::iterator it;
  it = treeRef->polygons.find(node_handle);
  if(it != treeRef->polygons.end()){
    return it->second;
  } else {
    // If polygon is not found, polygon is the same as bounding box  
    return IsotheticPolygon(rectangle); 
  }
}

template <int min_branch_factor, int max_branch_factor, class strategy>
IsotheticPolygon find_polygon(
                  NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
                  Branch branch){
  assert(branch.child != nullptr);
  std::map<tree_node_handle, IsotheticPolygon>::iterator it;
  it = treeRef->polygons.find(branch.child);
  if(it != treeRef->polygons.end()){
    return it->second;
  } else {
    // If polygon is not found, polygon is the same as rectangle
    return IsotheticPolygon(branch.boundingBox); 
  }
}

// Remove the polygon associated with this node from map 
template <int min_branch_factor, int max_branch_factor, class strategy>
void remove_polygon(
        NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
        tree_node_handle node_handle){
  assert(node_handle != nullptr);
  std::map<tree_node_handle, IsotheticPolygon>::iterator it;
  it = treeRef->polygons.find(node_handle);
  if(it != treeRef->polygons.end()){
    // If this branch does have a polygon, remove it from map 
    treeRef->polygons.erase(it); 
  } 
}

// Helper function for find_parent_handles 
template <int min_branch_factor, int max_branch_factor, class strategy>
bool find_path_to_leaf_helper(
      NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
      tree_node_handle start_handle,
      tree_node_handle leaf_handle,
      std::stack<tree_node_handle> &path_to_leaf  
) {
  assert(leaf_handle != nullptr);
  if(start_handle.get_type() == LEAF_NODE){
    assert(start_handle == leaf_handle);
    return true; 
  }
  assert(start_handle.get_type() == BRANCH_NODE);
  std::stack<tree_node_handle> candidates; 
  auto current_node = treeRef->get_branch_node(start_handle); 
  auto leaf_node = treeRef->get_leaf_node(leaf_handle);
  Rectangle leaf_mbb = leaf_node->boundingBox(); 

  for (unsigned i = 0; i < current_node->cur_offset_; i++) {
    Branch &b = current_node->entries.at(i);
    // find all branches which contain the leaf node 
    if (b.boundingBox.containsRectangle(leaf_mbb)) {
      candidates.push(b.child); 
    }
  }
  while(not candidates.empty()){
    tree_node_handle candidate_handle = candidates.top();
    candidates.pop();
    path_to_leaf.push(candidate_handle);
    bool found_leaf = find_path_to_leaf_helper(treeRef, candidate_handle, leaf_handle, path_to_leaf);
    if (found_leaf) return true;
    // not found leaf in this subbranch, pop the parent_handle, try next candidate
    path_to_leaf.pop();
  }
  return false; 
}

// Find_parent_handles returns a stack of parent handles from start_handle to leave
template <int min_branch_factor, int max_branch_factor, class strategy>
std::stack<tree_node_handle> find_path_to_leaf(
      NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
      tree_node_handle start_handle,
      tree_node_handle leaf_handle) {
  assert(start_handle.get_type() == BRANCH_NODE);
  assert(leaf_handle.get_type() == LEAF_NODE);

  std::stack<tree_node_handle> path_to_leaf;
  auto current_node = treeRef->get_branch_node(start_handle); 
  auto leaf_node = treeRef->get_leaf_node(leaf_handle);
  Rectangle leaf_mbb = leaf_node->boundingBox(); 
  assert(current_node->boundingBox().containsRectangle(leaf_mbb)); 
  path_to_leaf.push(start_handle);
  // we expect to find the leaf 
  bool found_leaf = find_path_to_leaf_helper(treeRef, start_handle, leaf_handle, path_to_leaf);
  assert(found_leaf);
  return path_to_leaf; 
}

inline std::pair<double, std::vector<IsotheticPolygon::OptimalExpansion>>
computeExpansionArea( const IsotheticPolygon &this_poly, const IsotheticPolygon &other_poly )
{
    std::vector<IsotheticPolygon::OptimalExpansion> expansions;
    double totalAreas = 0.0;

    for (const Rectangle &rect : other_poly.basicRectangles) {
        IsotheticPolygon::OptimalExpansion exp = this_poly.computeExpansionArea(rect);
        if( exp.area != 0.0 and exp.area != -1.0 ) {
            totalAreas += exp.area;
        }
        expansions.push_back( exp );
    }
    return std::make_pair( totalAreas == 0.0 ? -1.0 : totalAreas, expansions );
}


template <int min_branch_factor, int max_branch_factor, class strategy>
void testDisjoint(NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef,
                  tree_node_handle root, 
                  std::string msg){
#ifdef DEBUG_TESTDISJOINT
  std::stack<tree_node_handle> context; 
  context.push(root);
  while(not context.empty()){
    tree_node_handle current_handle = context.top();
    context.pop();
    if(current_handle.get_type() == LEAF_NODE) continue;
    auto current_node = treeRef->get_branch_node(current_handle);
    for (int i = 0; i < current_node->cur_offset_; ++i){
      Branch bi = current_node->entries[i];
      context.push(bi.child);
      IsotheticPolygon pi = find_polygon(treeRef, bi);
      for (int j = i+1; j < current_node->cur_offset_; ++j){
        Branch bj = current_node->entries[j];
        IsotheticPolygon pj = find_polygon(treeRef, bj);
        //std::cout << "checking " << bi.child <<" and " << bj.child << std::endl;
        if (! pi.disjoint(pj)){
          int height = current_node->height(treeRef,current_handle);
          auto root_node = treeRef->get_branch_node(root);
          int root_height = current_node->height(treeRef, root);
          std::cout << bi.child <<" intersects " << bj.child << std::endl;
          std::cout << " parent is " << current_handle << std::endl;
          std::cout << " height of parent is " << height << std::endl;
          std::cout << " height of root is " << root_height << std::endl;
          std::cout << " root is " << root << std::endl;
          std::cout << msg << std::endl;
          root_node->printTree(treeRef, root, tree_node_handle(nullptr), 0);
          assert(pi.disjoint(pj));
        }
      }
    }
    //std::cout << "no conflict on node " << current_handle << std::endl;
    // for (int i = 0; i < current_node->cur_offset_; ++i){
    //   Branch bi = current_node->entries[i];
    //   IsotheticPolygon pi = find_polygon(treeRef, bi);
    //   std::cout << "polygon of " << bi.child << std::endl;
    //   std::cout << pi << std::endl;
    // }
  }
#endif
}

template <int min_branch_factor, int max_branch_factor, class strategy>
void testCount(NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef, 
               tree_node_handle root, 
               std::string msg){
#ifdef DEBUG_TESTCOUNT
  std::stack<tree_node_handle> context; 
  context.push(root);
  int total_count = 0; 
  while(not context.empty()){
    tree_node_handle current_handle = context.top();
    context.pop();
    if(current_handle.get_type() == LEAF_NODE) {
      auto current_node = treeRef->get_leaf_node(current_handle);
      total_count += current_node->cur_offset_;
    } else {
      auto current_node = treeRef->get_branch_node(current_handle);
      for (int i = 0; i < current_node->cur_offset_; ++i){
        Branch bi = current_node->entries[i];
        context.push(bi.child);
      }
    }
    std::cout << "total " << total_count << std::endl;
  }
#endif 
}

template <int min_branch_factor, int max_branch_factor, class strategy>
void testContainPoints(NIRTreeDisk<min_branch_factor, max_branch_factor, strategy> *treeRef, 
                       tree_node_handle root, 
                       std::string msg){
#ifdef DEBUG_TESTCONTAINPOINTS
  std::stack<tree_node_handle> context; 
  context.push(root);
  while(not context.empty()){
    tree_node_handle current_handle = context.top();
    context.pop();
    if(current_handle.get_type() == LEAF_NODE) {
      auto current_leaf = treeRef->get_leaf_node(current_handle);
      IsotheticPolygon current_poly = find_polygon(treeRef, current_handle, current_leaf->boundingBox()); 
      for (int i = 0; i < current_leaf->cur_offset_; ++i){
        Point p = current_leaf->entries[i];
        if(! current_poly.containsPoint(p)){
          std::cout << "**ERROR:" << current_handle << "doesn't contain the point "<< p << std::endl;
        }
      }
    } else {
      auto current_node = treeRef->get_branch_node(current_handle);
      for (int i = 0; i < current_node->cur_offset_; ++i){
        Branch bi = current_node->entries[i];
        context.push(bi.child);
      }
    }
  }
#endif 
}


} // namespace nirtreedisk

