#pragma once

#include <cassert>
#include <cmath>
#include <globals/globals.h>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <memory>
#include <numeric>
#include <stack>
#include <storage/tree_node_allocator.h>
#include <util/geometry.h>
#include <util/repacking.h>
#include <util/statistics.h>
#include <utility>
#include <variant>
#include <vector>
#include <fstream>

namespace rplustreedisk {
    template <int min_branch_factor, int max_branch_factor> class RPlusTreeDisk;

    template <int min_branch_factor, int max_branch_factor>
    tree_node_allocator *get_node_allocator(RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef) {
      return treeRef->node_allocator_.get();
    }

    template <int min_branch_factor, int max_branch_factor>
    tree_node_handle get_root_handle(RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef) {
      return treeRef->root;
    }

    class Branch {
    public:
        Branch(Rectangle boundingBox, tree_node_handle child_handle):
                boundingBox(boundingBox), child(child_handle) {}

        Branch() = default;

        Branch(const Branch &other): boundingBox(other.boundingBox), child(other.child) {}

        bool operator==(const Branch &o) const = default;

        Rectangle boundingBox;
        tree_node_handle child;
    };

    struct SplitResult
    {
        Branch leftBranch;
        Branch rightBranch;
    };

    struct Partition
    {
        unsigned dimension;
        double location;
    };

    template <int min_branch_factor, int max_branch_factor>
    class LeafNode {
    public:
        // Obnoxiously, this needs to have a +1 so we can overflow
        // by 1 entry and deal with it later.
        std::array<Point, max_branch_factor + 1> entries;
        unsigned cur_offset_;

        // Constructors and destructors
        LeafNode(): cur_offset_(0) {}

        void addPoint(const Point &p) {
          entries.at(cur_offset_++) = p;
        }

        void deleteSubtrees();

        // Helper functions
        Rectangle boundingBox() const;
        void removePoint(const Point &givenPoint);

        tree_node_handle findLeaf(
                RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
                tree_node_handle selfHandle,
                const Point &givenPoint
        );
        unsigned chooseSplitLeafAxis();
        unsigned chooseSplitNonLeafAxis();
        unsigned chooseSplitAxis();
        unsigned chooseSplitIndex(unsigned axis);
        SplitResult splitNode(
                RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
                tree_node_handle current_handle,
                Partition p
        );
        tree_node_handle adjustTree(
                RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
                tree_node_handle current_handle,
                tree_node_handle sibling_handle,
                std::stack<tree_node_handle> parentHandles,
                std::vector<bool> &hasReinsertedOnLevel
        );
        tree_node_handle overflowTreatment(
                RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
                tree_node_handle current_handle,
                std::stack<tree_node_handle> parentHandles,
                std::vector<bool> &hasReinsertedOnLevel
        );
        tree_node_handle condenseTree(std::vector<bool> &hasReinsertedOnLevel);
        Partition partitionNode();

        // Datastructure interface functions
        void exhaustiveSearch(const Point &requestedPoint, std::vector<Point> &accumulator) const;

        // These return the root of the tree.
        tree_node_handle insert(
                RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
                tree_node_handle current_handle,
                Point point
        );
        tree_node_handle remove(Point &givenPoint);

        // Miscellaneous
        unsigned checksum() const;
        void print() const;
        void printTree() const;
        unsigned height() const;

        // Operators
        bool operator<(const LeafNode &otherNode) const;
    };

    template <int min_branch_factor, int max_branch_factor>
    class BranchNode {
    private:
        void searchSub(const Point &requestedPoint, std::vector<Point> &accumulator);
        void searchSub(const Rectangle &rectangle, std::vector<Point> &accumulator);

    public:
        // Obnoxiously, this needs to have a +1 so we can overflow
        // by 1 entry and deal with it later.
        // Brad: This is needed for R-star overflow insertion.
        std::array<Branch, max_branch_factor + 1> entries;
        unsigned cur_offset_;

        // Constructors and destructors
        BranchNode() {
          cur_offset_ = 0;
        }

        void addBranchToNode(Rectangle boundingBox, tree_node_handle child) {
          Branch b(boundingBox, child);
          entries.at(cur_offset_++) = b;
        };

        void addBranchToNode(const Branch &b) {
          entries.at(cur_offset_++) = b;
        }

        void deleteSubtrees();

        // Helper functions
        Rectangle boundingBox() const;
        bool updateBoundingBox(tree_node_handle child, Rectangle updatedBoundingBox);
        void removeChild(tree_node_handle child);
        tree_node_handle chooseSubtree(
                RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
                tree_node_handle current_handle,
                std::stack<tree_node_handle> &parentHandles,
                const Point &givenPoint
        );
        tree_node_handle findLeaf(
                RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
                tree_node_handle selfHandle,
                const Point &givenPoint
        );
        unsigned chooseSplitLeafAxis();
        unsigned chooseSplitNonLeafAxis();
        unsigned chooseSplitAxis();
        unsigned chooseSplitIndex(unsigned axis);
        SplitResult splitNode(
                RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
                tree_node_handle current_handle,
                Partition p
        );
        tree_node_handle adjustTree(
                RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
                tree_node_handle current_handle,
                tree_node_handle sibling_handle,
                std::stack<tree_node_handle> parentHandles,
                std::vector<bool> &hasReinsertedOnLevel
        );
        tree_node_handle overflowTreatment(
                RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
                tree_node_handle current_handle,
                std::stack<tree_node_handle> parentHandles,
                std::vector<bool> &hasReinsertedOnLevel
        );
        tree_node_handle condenseTree(std::vector<bool> &hasReinsertedOnLevel);
        Partition partitionNode();

        // Datastructure interface functions
        void exhaustiveSearch(const Point &requestedPoint, std::vector<Point> &accumulator) const;

        // These return the root of the tree.
        tree_node_handle insert(
                RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
                tree_node_handle current_handle,
                Point point
        );
        tree_node_handle remove(Point &givenPoint);

        // Miscellaneous
        unsigned checksum() const;
        void print() const;
        void printTree() const;
        unsigned height(tree_node_handle self_handle) const;

        // Operators
        bool operator<(const BranchNode &otherNode) const;
    };

    template <class NE, class B, int N>
    double computeOverlapGrowth(unsigned index, const std::array<B, N + 1> &entries,
                                unsigned els_to_consider,
                                const Rectangle &givenBox)
    {
      // We cannot be a leaf
      assert(els_to_consider > 0);

      // 1. Make a test rectangle we will use to not modify the original
      const Rectangle &origRectangle = entries[index].boundingBox;
      Rectangle newRectangle = entries[index].boundingBox;

      // 2. Add the point to the copied Rectangle
      newRectangle.expand(givenBox);

      // 3. Compute the overlap expansion area
      double overlapDiff = 0;
      unsigned num_entries_els = els_to_consider;

      for (unsigned i = 0; i < num_entries_els; ++i) {
        const auto &entry = entries[i];

        if (i == index) {
          continue;
        }

        overlapDiff +=
                (newRectangle.computeIntersectionArea(entry.boundingBox) -
                 origRectangle.computeIntersectionArea(entry.boundingBox));
      }

      return overlapDiff;
    }

    template <int min_branch_factor, int max_branch_factor, typename functor>
    void treeWalker(RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef, tree_node_handle root, functor &f) {
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

    template <int min_branch_factor, int max_branch_factor>
    void LeafNode<min_branch_factor, max_branch_factor>::deleteSubtrees() {
      return;
    }

    template <int min_branch_factor, int max_branch_factor>
    Rectangle LeafNode<min_branch_factor, max_branch_factor>::boundingBox() const {
      assert(cur_offset_ > 0);
      Rectangle boundingBox(entries[0], Point::closest_larger_point(entries[0]));

      for (unsigned i = 0; i < cur_offset_; i++) {
        boundingBox.expand(entries.at(i));
      }

      return boundingBox;
    }

    template <int min_branch_factor, int max_branch_factor>
    void LeafNode<min_branch_factor, max_branch_factor>::removePoint(const Point &givenPoint) {
      unsigned i = 0;
      for (; i < cur_offset_; i++) {
        if (entries.at(i) == givenPoint) {
          break;
        }
      }
      assert(givenPoint == entries.at(i));

      entries.at(i) = entries.at(cur_offset_ - 1);
      cur_offset_--;
    }

    template <int min_branch_factor, int max_branch_factor>
    void LeafNode<min_branch_factor, max_branch_factor>::exhaustiveSearch(const Point &requestedPoint, std::vector<Point> &accumulator) const {
      for (unsigned i = 0; i < cur_offset_; i++) {
        const Point &p = entries.at(i);

        if (p == requestedPoint) {
          accumulator.push_back(p);
        }
      }
    }

    template <int min_branch_factor, int max_branch_factor>
    tree_node_handle LeafNode<min_branch_factor, max_branch_factor>::findLeaf(
            RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
            tree_node_handle selfHandle,
            const Point &givenPoint
    ) {
      for (unsigned i = 0; i < cur_offset_; i++) {
        if (entries.at(i) == givenPoint) {
          return selfHandle;
        }
      }

      return tree_node_handle(nullptr);
    }

    template <int min_branch_factor, int max_branch_factor>
    unsigned LeafNode<min_branch_factor, max_branch_factor>::chooseSplitLeafAxis() {
      unsigned optimalAxis = 0;
      double optimalMargin = std::numeric_limits<double>::infinity();

      // Make the entries easier to work with
      std::vector<Point *> entriesCopy;
      entriesCopy.reserve(cur_offset_);
      for (unsigned i = 0; i < cur_offset_; i++) {
        entriesCopy.push_back(&(entries.at(i)));
      }

      // Consider all M-2m+2 distributions in each dimension
      for (unsigned d = 0; d < dimensions; d++) {
        // First sort in the current dimension
        std::sort(entriesCopy.begin(), entriesCopy.end(), [d](Point *a, Point *b) {
            return (*a)[d] < (*b)[d];
        });

        // Setup groups
        std::vector<Point *> groupA(entriesCopy.begin(), entriesCopy.begin() + min_branch_factor);
        std::vector<Point *> groupB(entriesCopy.begin() + min_branch_factor, entriesCopy.end());

        // Cycle through all M-2m+2 distributions
        double totalMargin = 0.0;
        while (groupA.size() <= max_branch_factor and groupB.size() >= min_branch_factor) {
          // Compute the margin of groupA and groupB
          Rectangle boundingBoxA(*groupA[0], Point::closest_larger_point(*groupA[0]));
          for (unsigned i = 1; i < groupA.size(); i++) {
            boundingBoxA.expand(*groupA[i]);
          }

          Rectangle boundingBoxB(*groupB[0], Point::closest_larger_point(*groupB[0]));
          for (unsigned i = 1; i < groupB.size(); i++) {
            boundingBoxB.expand(*groupB[i]);
          }

          // Add to the total margin sum
          totalMargin += boundingBoxA.margin() + boundingBoxB.margin();

          // Add one new value to groupA and remove one from groupB to obtain next distribution
          Point *transferPoint = groupB.front();
          groupB.erase(groupB.begin());
          groupA.push_back(transferPoint);
        }

        if (totalMargin < optimalMargin) {
          optimalMargin = totalMargin;
          optimalAxis = d;
        }
      }

      // Sort along our best axis
      std::sort(entries.begin(), entries.begin() + cur_offset_,
                [optimalAxis](Point &a, Point &b) { return a[optimalAxis] < b[optimalAxis]; });

      return optimalAxis;
    }

// CSA1: Sort entries by lower and upper bound along each axis and compute S -> sum of all
//  margin values for the different distributions. This can be stored in a array of variable
//  that we keep in a loop -> and the just compare to the others?
// 	We can first call a helper function that returns an array of all possible distributions for it?
// CSA2: Return the Axis that has the minimum total sum of all the distributions
    template <int min_branch_factor, int max_branch_factor>
    unsigned LeafNode<min_branch_factor, max_branch_factor>::chooseSplitAxis() {
      return chooseSplitLeafAxis();
    }

// CSI1: Given the chosen split index
// 	group all the entries into multiple groups and choose the one that has the least
// 	overlap value; resolve ties with the minimum area
// 	returns tuple of best distribution group indices
    template <int min_branch_factor, int max_branch_factor>
    unsigned LeafNode<min_branch_factor, max_branch_factor>::chooseSplitIndex(unsigned axis) {
      // We assume this is called after we have sorted this->data according to axis.

      const auto groupABegin = entries.begin();
      const auto groupAEnd = entries.begin() + min_branch_factor;
      const auto groupBBegin = entries.begin() + min_branch_factor;
      const auto groupBEnd = entries.begin() + cur_offset_;

      std::vector<Point> groupA(groupABegin, groupAEnd);
      std::vector<Point> groupB(groupBBegin, groupBEnd);
      unsigned splitIndex = cur_offset_ / 2;

      // Find the best size out of all the distributions
      double minOverlap = std::numeric_limits<double>::infinity();
      double minArea = std::numeric_limits<double>::infinity();

      // Tracking what the current "cut" mark is
      unsigned currentSplitPoint = min_branch_factor;

      // Try each of the M-2m + 2 groups
      while (groupA.size() <= max_branch_factor and groupB.size() >= min_branch_factor) {
        // Compute the margin of groupA and groupB
        Rectangle boundingBoxA(groupA[0], Point::closest_larger_point(groupA[0]));
        for (unsigned i = 1; i < groupA.size(); i++) {
          boundingBoxA.expand(groupA[i]);
        }

        Rectangle boundingBoxB(groupB[0], Point::closest_larger_point(
                groupB[0]));
        for (unsigned i = 1; i < groupB.size(); i++) {
          boundingBoxB.expand(groupB[i]);
        }

        // Compute intersection area to determine best grouping of data points
        double evalDistOverlap = boundingBoxA.computeIntersectionArea(boundingBoxB);

        if (evalDistOverlap < minOverlap) {
          // We save this current distribution of indices to return
          minOverlap = evalDistOverlap;
          splitIndex = currentSplitPoint;

          // Set this if we haven't already
          if (minArea == std::numeric_limits<double>::infinity()) {
            minArea = boundingBoxA.area() + boundingBoxB.area();
          }
        } else if (evalDistOverlap == minOverlap) {
          // If overlap is equal, we use the distribution that creates the smallest areas
          double evalMinArea = boundingBoxA.area() + boundingBoxB.area();

          if (evalMinArea < minArea) {
            // Save this current distribution of indices to return
            minArea = evalMinArea;
            splitIndex = currentSplitPoint;
          }
        }

        // Add one new value to groupA and remove one from groupB to obtain next distribution
        Point transferPoint = groupB.front();
        groupB.erase(groupB.begin());
        groupA.push_back(transferPoint);

        // Push the split point forward.
        currentSplitPoint++;
      }

      return splitIndex;
    }

    template <int min_branch_factor, int max_branch_factor>
    Partition LeafNode<min_branch_factor, max_branch_factor>::partitionNode() {
      Partition defaultPartition;
      unsigned costMetric = std::numeric_limits<unsigned>::max();
      double location;

      for (unsigned d = 0; d < dimensions; d++) {
        // Sort along dimension d
        std::sort(
          entries.begin(), entries.begin() + cur_offset_,
          [d](Point &a, Point &b) { return a[d] < b[d]; }
        );

        // Pick at least half the data
        location = entries.at(cur_offset_ / 2 - 1)[d];

        // Compute cost, # of duplicates of this location
        unsigned duplicateCount = 0;
        for (unsigned i = 0; i < cur_offset_; i++) {
          Point &data_point = entries.at(i);

          if (location == data_point[d]) {
            duplicateCount++;
          }
        }

        // Compare cost
        if (duplicateCount < costMetric) {
          defaultPartition.dimension = d;
          // Set the partition location after the point.
          defaultPartition.location = nextafter(location, DBL_MAX);
          costMetric = duplicateCount;
        }
      }

      return defaultPartition;
    }

    template <int min_branch_factor, int max_branch_factor>
    SplitResult LeafNode<min_branch_factor, max_branch_factor>::splitNode(
      RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
      tree_node_handle current_handle,
      Partition p
    ) {
      using NodeType = LeafNode<min_branch_factor, max_branch_factor>;

      tree_node_allocator *allocator = get_node_allocator(treeRef);

      // Create a copy of entries
      std::vector<Point> entriesCopy(entries);

      // We are the left child
      auto left_handle = current_handle;
      auto left_node = this;
      left_node->cur_offset_ = 0; // Reset our entries and repopulate below

      // Create a new sibling
      auto alloc_data = allocator->create_new_tree_node<NodeType>(NodeHandleType(LEAF_NODE));
      new (&(*(alloc_data.first))) NodeType();
      auto right_handle = alloc_data.second;
      auto right_node = alloc_data.first;
      right_handle.set_level(left_handle.get_level());

      // Loop through entries and assign points to left or right nodes
      for (unsigned i = 0; i < entriesCopy.size(); i++) {
        Point data_point = entriesCopy.at(i);

        if (data_point[p.dimension] < p.location and left_node->cur_offset_ < max_branch_factor ) {
          left_node->entries.at(left_node->cur_offset_++) = data_point;
        } else {
          right_node->entries.at(right_node->cur_offset_++) = data_point;
        }

        assert( left_node->cur_offset_ <= max_branch_factor );
        assert( right_node->cur_offset_ <= max_branch_factor );
      }

      return {
        {left_handle, left_node->boundingBox()},
        {right_handle, right_node->boundingBox()}
      };
    }

// Note: This function modifies parentHandles.
    template <class NT, class TR>
    std::pair<tree_node_handle, tree_node_handle> adjustTreeBottomHalf(
            TR *treeRef,
            NT node,
            NT sibling,
            tree_node_handle current_handle,
            tree_node_handle sibling_handle,
            std::stack<tree_node_handle> &parentHandles,
            std::vector<bool> &hasReinsertedOnLevel,
            int max_branch_factor)
    {
      // The caller should ensure that parentHandles is not empty.
      assert(!parentHandles.empty());

      tree_node_handle node_handle = current_handle;
      tree_node_handle parent_handle = parentHandles.top();
      parentHandles.pop();

      /* After a split, we have two nodes at the same level: the current node
       * and the newly created sibling node. Current node has a new bounding box
       * after the split and needs to be updated in the parent. The sibling node
       * needs to be inserted into the parent. */
      auto parent_ptr = treeRef->get_branch_node(parent_handle);
      bool didUpdateBoundingBox = parent_ptr->updateBoundingBox(node_handle, node->boundingBox());

      // If we have a split then deal with it otherwise move up the tree
      if (sibling != nullptr) {
        Rectangle bb = sibling->boundingBox();
        // AT4 [Propogate the node split upwards]
        Branch b(bb, sibling_handle);
        parent_ptr->addBranchToNode(b);

        unsigned sz = parent_ptr->cur_offset_;

        if (sz >= (unsigned) max_branch_factor) {
          tree_node_handle parent_before_handle = parent_handle;
          tree_node_handle sibling_parent_handle = parent_ptr->overflowTreatment(
                  treeRef,
                  parent_handle,
                  parentHandles,
                  hasReinsertedOnLevel
          );

          if (sibling_parent_handle) {
            // We split our parent, so now we have two (possible) parents
            // Need to keep traversing up
            node_handle = parent_before_handle;
            sibling_handle = sibling_parent_handle;
            assert(node != sibling);

            return std::make_pair(node_handle, sibling_handle);
          }
        }

        node_handle = parent_handle;
        sibling_handle = tree_node_handle(nullptr);

        return std::make_pair(node_handle, sibling_handle);
      }

      // AT5 [Move up to next level]
      if (didUpdateBoundingBox) {
        node_handle = parent_handle;
      } else {
        // If we didn't update our bounding box and there was no split, no reason to keep
        // going.
        return std::make_pair(tree_node_handle(nullptr),
                              tree_node_handle(nullptr));
      }

      return std::make_pair(node_handle, tree_node_handle(nullptr));
    }

    template <class TR, class NT>
    SplitResult adjustTreeSubHelper(
      TR *treeRef,
      NT current_node,
      tree_node_handle current_handle,
      tree_node_handle parent_handle,
      SplitResult currentPropagationSplit,
      int max_branch_factor
    ) {
      // Returns a new split that needs to be propagated upwards
      SplitResult newPropagationSplit;

      if (currentPropagationSplit.leftBranch.child != nullptr and currentPropagationSplit.rightBranch.child != nullptr) {
        if (current_node->cur_offset_ >= max_branch_factor) {
          assert(false);
        }

        // If there was a split we were supposed to propagate then propagate it
        if (currentPropagationSplit.leftBranch.child.get_type() == LEAF_NODE) {
          auto left_node = treeRef->get_leaf_node(currentPropagationSplit.leftBranch.child);
          if (left_node->cur_offset_ > 0) {
            current_node->entries.at(current_node->cur_offset_++) = currentPropagationSplit.leftBranch;
          }
        } else {
          auto left_node = treeRef->get_branch_node(currentPropagationSplit.leftBranch.child);
          if (left_node->cur_offset_ > 0) {
            current_node->entries.at(current_node->cur_offset_++) = currentPropagationSplit.leftBranch;
          }
        }

        if (currentPropagationSplit.rightBranch.child.get_type() == LEAF_NODE) {
          auto right_node = treeRef->get_leaf_node(currentPropagationSplit.rightBranch.child);
          if (right_node->cur_offset_ > 0) {
            current_node->entries.at(current_node->cur_offset_++) = currentPropagationSplit.rightBranch;
          }
        } else {
          auto right_node = treeRef->get_branch_node(currentPropagationSplit.rightBranch.child);
          if (right_node->cur_offset_ > 0) {
            current_node->entries.at(current_node->cur_offset_++) = currentPropagationSplit.rightBranch;
          }
        }
      }

      // Early exit if this node does not overflow
      if (current_node->cur_offset_ <= max_branch_factor) {
        newPropagationSplit = {
          {Rectangle(), tree_node_handle( nullptr )},
          {Rectangle(), tree_node_handle( nullptr )}
        };
        return newPropagationSplit;
      }

      // Otherwise, split node
      newPropagationSplit = current_node->splitNode();

      // The current node has been split into two new nodes, remove it from the parent
      if (parent_handle != nullptr) {
        // parent is guaranteed to be a branch node
        auto parent_node = treeRef->get_branch_node(parent_handle);
        assert(parent_node->cur_offset_ <= max_branch_factor);

        parent_node->removeChild(current_handle);
        assert(parent_node->cur_offset_ <= max_branch_factor - 1);
      }

      return newPropagationSplit;
    }

    template <int min_branch_factor, int max_branch_factor>
    SplitResult adjustTreeSub(
      RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
      tree_node_handle current_handle,
      std::stack<tree_node_handle> parentHandles
    ) {
      SplitResult propagationSplit = {
        {Rectangle(), tree_node_handle(nullptr)},
        {Rectangle(), tree_node_handle(nullptr)}
      };

      return propagationSplit;
    }

    template <int min_branch_factor, int max_branch_factor>
    tree_node_handle LeafNode<min_branch_factor, max_branch_factor>::adjustTree(
            RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
            tree_node_handle current_handle,
            tree_node_handle sibling_handle,
            std::stack<tree_node_handle> parentHandles,
            std::vector<bool> &hasReinsertedOnLevel
    ) {
      return adjustTreeSub(treeRef, current_handle, sibling_handle, parentHandles, hasReinsertedOnLevel);
    }

    template <int min_branch_factor, int max_branch_factor>
    tree_node_handle BranchNode<min_branch_factor, max_branch_factor>::adjustTree(
            RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
            tree_node_handle current_handle,
            tree_node_handle sibling_handle,
            std::stack<tree_node_handle> parentHandles,
            std::vector<bool> &hasReinsertedOnLevel
    ) {
      return adjustTreeSub(treeRef, current_handle, sibling_handle, parentHandles, hasReinsertedOnLevel);
    }

// Overflow treatment for dealing with a node that is too big (overflow)
    template <int min_branch_factor, int max_branch_factor>
    tree_node_handle LeafNode<min_branch_factor, max_branch_factor>::overflowTreatment(
            RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
            tree_node_handle current_handle,
            std::stack<tree_node_handle> parentHandles,
            std::vector<bool> &hasReinsertedOnLevel
    ) {
      uint16_t current_level = current_handle.get_level();
      assert(hasReinsertedOnLevel.size() > current_level);

      if (hasReinsertedOnLevel.at(current_level)) {
        //std::cout << "Overflow treatment on leaf node, splitting." <<
        //    std::endl;
//        return splitNode(treeRef, current_handle, partitionNode());
      }
    }

// insert() is always called on the root node. If this gets called,
// that means the tree only has one leaf node.
    template <int min_branch_factor, int max_branch_factor>
    tree_node_handle LeafNode<min_branch_factor, max_branch_factor>::insert(
            RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
            tree_node_handle current_handle,
            Point point
    ) {
      std::vector<bool> hasReinsertedOnLevel;

      tree_node_allocator *allocator = get_node_allocator(treeRef);

      // I1 [Find position for new record]
      tree_node_handle sibling_handle = tree_node_handle(nullptr);

      // I2 [Add record to leaf node]
      addPoint(point);

      // Empty parentHandles since we are the only node in the tree.
      std::stack<tree_node_handle> parentHandles;

      uint16_t current_level = current_handle.get_level();
      assert(current_level == 0); // Leaf nodes have level = 0

      // If we exceed treeRef->maxBranchFactor we need to do something about it
      if (cur_offset_ > max_branch_factor) {
        // We call overflow treatment to determine how our sibling node is treated. If we do a
        // reInsert, sibling is nullptr. This is properly dealt with in adjustTree.
        sibling_handle = overflowTreatment(
                treeRef,
                current_handle,
                parentHandles,
                hasReinsertedOnLevel
        );
      }

      // I3 [Propagate overflow treatment changes upward]
      sibling_handle = adjustTree(
              treeRef,
              current_handle,
              sibling_handle,
              parentHandles,
              hasReinsertedOnLevel
      );

      // I4 [Grow tree taller]
      if (sibling_handle) {
        // Sanity check that we're still the root
        assert(treeRef->root == current_handle);

        auto alloc_data =
                allocator->create_new_tree_node<BranchNode<min_branch_factor, max_branch_factor>>(
                        NodeHandleType(BRANCH_NODE));
        auto newRoot = alloc_data.first;
        tree_node_handle new_root_handle = alloc_data.second;
        new_root_handle.set_level(current_level + 1);

        auto sibling = treeRef->get_leaf_node(sibling_handle);

        new (&(*(newRoot))) BranchNode<min_branch_factor, max_branch_factor>();

        // Make the existing root a child of newRoot
        Branch b1(boundingBox(), current_handle);
        newRoot->addBranchToNode(b1);

        // Make the new sibling node a child of newRoot
        Branch b2(sibling->boundingBox(), sibling_handle);
        newRoot->addBranchToNode(b2);

        // Ensure newRoot has both children
        assert(newRoot->cur_offset_ == 2);
        assert(sibling_handle.get_level() + 1 == new_root_handle.get_level());

        // Fix the reinserted length
        hasReinsertedOnLevel.push_back(false);

        return new_root_handle;
      }

      return treeRef->root;
    }

// To be called on a leaf
    template <int min_branch_factor, int max_branch_factor>
    tree_node_handle LeafNode<min_branch_factor, max_branch_factor>::condenseTree(std::vector<bool> &hasReinsertedOnLevel) {
#if 0
      // CT1 [Initialize]
  tree_node_handle node_handle = self_handle_;

  std::vector<NodeEntry> Q;

  // CT2 [Find parent entry]
  unsigned entriesSize;

  for (;;) {

    tree_node_handle parent_handle;
    Rectangle bb;
    if (node_handle.get_type() == LEAF_NODE) {
      auto node = treeRef->get_leaf_node(node_handle);
      parent_handle = node->parent;
      entriesSize = node->cur_offset_;
      bb = node->boundingBox();
    } else {
      auto node = treeRef->get_branch_node(node_handle);
      parent_handle = node->parent;
      entriesSize = node->cur_offset_;
      bb = node->boundingBox();
    }

    if (!parent_handle) {
      break;
    }

    // CT3 & CT4 [Eliminate under-full node. & Adjust covering rectangle.]
    if (entriesSize >= min_branch_factor) {
      auto parent = treeRef->get_branch_node(parent_handle);
      parent->updateBoundingBox(node_handle, bb);

      // CT5 [Move up one level in the tree]
      // Move up a level without deleting ourselves
      node_handle = parent->self_handle_;
    } else {
      auto parent = treeRef->get_branch_node(parent_handle);
      // Remove ourselves from our parent
      parent->removeChild(node_handle);

      if (node_handle.get_type() == LEAF_NODE) {
        auto node = treeRef->get_leaf_node(node_handle);
        // Push these entries into Q
        std::copy(node->entries.begin(),
                  node->entries.begin() + node->cur_offset_, std::back_inserter(Q));
      } else {
        auto node = treeRef->get_leaf_node(node_handle);
        // Push these entries into Q
        std::copy(node->entries.begin(),
                  node->entries.begin() + node->cur_offset_, std::back_inserter(Q));
      }

      // FIXME: Should garbage collect node_ptr, it is dead now
      //tree_node_handle garbage = node_handle;

      node_handle = parent->self_handle_;
      // Cleanup ourselves without deleting children b/c they will be reinserted
      // GarbageCollect( node_ptr );
    }
  }

  // CT6 [Re-insert oprhaned entries]
  for (const auto &entry : Q) {
    if (node_handle.get_type() == LEAF_NODE) {
      auto root = treeRef->get_leaf_node(node_handle);
      assert(!root->parent);
      node_handle = root->insert(std::get<Point>(entry), hasReinsertedOnLevel);
    } else {
      auto root = treeRef->get_branch_node(node_handle);
      assert(!root->parent);
      node_handle = root->insert(entry, hasReinsertedOnLevel);
    }
  }

  return node_handle;
#endif

      // Unsupported
      abort();
    }

// Always called on root, this = root
    template <int min_branch_factor, int max_branch_factor>
    tree_node_handle LeafNode<min_branch_factor, max_branch_factor>::remove(Point &givenPoint) {
#if 0
      removePoint(givenPoint);

  // D3 [Propagate changes]
  tree_node_handle root_handle = condenseTree(hasReinsertedOnLevel);
  if (root_handle.get_type() == LEAF_NODE) {
    return root_handle;
  }

  auto root = treeRef->get_branch_node(root_handle);

  // D4 [Shorten tree]
  if (root->cur_offset_ == 1) {
    // Slice the hasReinsertedOnLevel
    hasReinsertedOnLevel.pop_back();

    // We are removing the root to shorten the tree so we then decide to remove the root
    Branch &b = root->entries[0];

    // Get rid of the old root

    if (b.child.get_type() == LEAF_NODE) {
      auto child = treeRef->get_leaf_node(b.child);
      child->parent = tree_node_handle(nullptr);
      // Garbage Collect Root
      // FIXME GC(root);

      return b.child;
    } else {
      auto child = treeRef->get_branch_node(b.child);
      child->parent = tree_node_handle(nullptr);
      // Garbage Collect Root
      // FIXME GC(root);
      return b.child;
    }
  }
  return root_handle;
#endif

      // Unsupported
      abort();
    }

    template <int min_branch_factor, int max_branch_factor>
    void LeafNode<min_branch_factor, max_branch_factor>::print() const {

      std::string indentation(4, ' ');
      std::cout << indentation << "Node " << (void *)this << std::endl;
      std::cout << indentation << "{" << std::endl;
      std::cout << indentation << "    BoundingBox: " << boundingBox() << std::endl;
//  std::cout << indentation << "    Parent: " << parent << std::endl;
      std::cout << indentation << "    Entries: " << std::endl;

      for (unsigned i = 0; i < cur_offset_; i++) {
        std::cout << indentation << "		" << entries.at(i) << std::endl;
      }
      std::cout << std::endl
                << indentation << "}" << std::endl;
    }

    template <int min_branch_factor, int max_branch_factor>
    void LeafNode<min_branch_factor, max_branch_factor>::printTree() const {
      // No op
    }

    template <int min_branch_factor, int max_branch_factor>
    unsigned LeafNode<min_branch_factor, max_branch_factor>::checksum() const {
      unsigned checksum = 0;
      for (unsigned i = 0; i < cur_offset_; i++) {
        auto &p = entries.at(i);
        for (unsigned d = 0; d < dimensions; d++) {
          checksum += (unsigned)p[d];
        }
      }

      return checksum;
    }

    template <int min_branch_factor, int max_branch_factor>
    unsigned LeafNode<min_branch_factor, max_branch_factor>::height() const {
      return 1;
    }

//// BRANCH STARTS
    template <int min_branch_factor, int max_branch_factor>
    void BranchNode<min_branch_factor, max_branch_factor>::deleteSubtrees() {
#if 0
      for (unsigned i = 0; i < cur_offset_; i++) {
    const Branch &b = entries.at(i);
    tree_node_handle child_handle = b.child;
    if (child_handle.get_type() == LEAF_NODE) {
      auto child = treeRef->get_leaf_node(child_handle);
      child->deleteSubtrees();
    } else {
      auto child = treeRef->get_branch_node(child_handle);
      child->deleteSubtrees();
    }
  }
#endif
    }

    template <int min_branch_factor, int max_branch_factor>
    Rectangle BranchNode<min_branch_factor, max_branch_factor>::boundingBox() const {
      assert(cur_offset_ > 0);
      Rectangle boundingBox = entries.at(0).boundingBox;

      for (unsigned i = 1; i < cur_offset_; i++) {
        boundingBox.expand(entries.at(i).boundingBox);
      }

      return boundingBox;
    }

    template <int min_branch_factor, int max_branch_factor>
    bool BranchNode<min_branch_factor, max_branch_factor>::updateBoundingBox(tree_node_handle child, Rectangle updatedBoundingBox) {

      //std::cout << "Updating bounding box for to " <<
      //    updatedBoundingBox << std::endl;
      for (unsigned i = 0; i < cur_offset_; i++) {
        Branch &b = entries.at(i);
        if (b.child == child) {
          //std::cout << "Existing bounding box: " << b.boundingBox <<
          //    std::endl;
          if (b.boundingBox != updatedBoundingBox) {
            b.boundingBox = updatedBoundingBox;
            return true;
          }
          return false;
        }
      }

      /* If control-flow reaches here, that means the child wasn't found in
       * the parent and no update happened.
       *
       * One case where this can happen is when a reinsertion happens which causes
       * a split that changes the tree structure. In this case, the child may belong
       * to a new parent, but parentHandles has the old parent. Since the recursive
       * insert call from reinsertion already takes care of updating the tree correctly,
       * it's OK to return false and not throw an error. */
      return false;
    }

    template <int min_branch_factor, int max_branch_factor>
    void BranchNode<min_branch_factor, max_branch_factor>::removeChild(tree_node_handle child) {
      for (unsigned i = 0; i < cur_offset_; i++) {
        if (entries.at(i).child == child) {
          entries.at(i) = entries.at(cur_offset_ - 1);
          cur_offset_--;
          return;
        }
      }
    }

    template <int min_branch_factor, int max_branch_factor>
    void BranchNode<min_branch_factor, max_branch_factor>::exhaustiveSearch(const Point &requestedPoint, std::vector<Point> &accumulator) const {
#if 0
      for (unsigned i = 0; i < cur_offset_; i++) {
    tree_node_handle child_handle = entries.at(i).child;
    if (child_handle.get_type() == LEAF_NODE) {
      auto child = treeRef->get_leaf_node(child_handle);
      child->exhaustiveSearch(requestedPoint, accumulator);
    } else {
      auto child = treeRef->get_branch_node(child_handle);
      child->exhaustiveSearch(requestedPoint, accumulator);
    }
  }
#endif

      // Unsupported
      abort();
    }

    template <int min_branch_factor, int max_branch_factor>
    void point_search_leaf_node(LeafNode<min_branch_factor, max_branch_factor> &node,
                                Point &requestedPoint,
                                std::vector<Point> &accumulator,
                                RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef)
    {
      unsigned intersection_count = 0;

      for (unsigned i = 0; i < node.cur_offset_; i++) {
        const Point &p = node.entries.at(i);
        intersection_count++;
        if (p == requestedPoint) {
          accumulator.push_back(p);
        }
      }

      treeRef->stats.recordIntersectionCount(intersection_count);
    }

    template <int min_branch_factor, int max_branch_factor>
    void point_search_branch_node(BranchNode<min_branch_factor, max_branch_factor> &node,
                                  Point &requestedPoint,
                                  std::stack<tree_node_handle> &context,
                                  RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef)
    {
      unsigned matching_branch_counter = 0;
      unsigned intersection_count = 0;

      for (size_t i = 0; i < node.cur_offset_; i++) {
        Branch &b = node.entries.at(i);
        intersection_count++;
        if (b.boundingBox.containsPoint(requestedPoint)) {
          context.push(b.child);
          matching_branch_counter++;
        }
      }

      treeRef->stats.recordScatter(matching_branch_counter);
      treeRef->stats.recordIntersectionCount(intersection_count);
    }

    template <int min_branch_factor, int max_branch_factor>
    std::vector<Point> point_search(
            tree_node_handle start_point,
            Point &requestedPoint,
            RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef)
    {
      std::vector<Point> accumulator;
      std::stack<tree_node_handle> context;

      context.push(start_point);

      while (!context.empty()) {
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

    template <int min_branch_factor, int max_branch_factor>
    void rectangle_search_leaf_node(LeafNode<min_branch_factor, max_branch_factor> &node,
                                    Rectangle &requestedRectangle,
                                    std::vector<Point> &accumulator,
                                    RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef)
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

    template <int min_branch_factor, int max_branch_factor>
    void rectangle_search_branch_node(BranchNode<min_branch_factor, max_branch_factor> &node,
                                      Rectangle &requestedRectangle,
                                      std::stack<tree_node_handle> &context,
                                      RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef)
    {
      unsigned intersection_count = 0;

      for (size_t i = 0; i < node.cur_offset_; i++) {
        Branch &b = node.entries.at(i);
        intersection_count++;
        if (b.boundingBox.intersectsRectangle(requestedRectangle)) {
          context.push(b.child);
        }
      }

      treeRef->stats.recordIntersectionCount(intersection_count);
    }

    template <int min_branch_factor, int max_branch_factor>
    std::vector<Point> rectangle_search(
            tree_node_handle start_point,
            Rectangle &requestedRectangle,
            RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef)
    {
      std::vector<Point> accumulator;

      std::stack<tree_node_handle> context;
      context.push(start_point);

      while (not context.empty()) {
        tree_node_handle current_handle = context.top();
        context.pop();

        if (current_handle.get_type() == LEAF_NODE) {
          auto current_node = treeRef->get_leaf_node(current_handle);
          rectangle_search_leaf_node(*current_node, requestedRectangle, accumulator, treeRef);
#ifdef STAT
          treeRef->stats.markLeafSearched();
#endif
        } else if (current_handle.get_type() == BRANCH_NODE) {
          auto current_node = treeRef->get_branch_node(current_handle);
          rectangle_search_branch_node(*current_node, requestedRectangle, context, treeRef);
#ifdef STAT
          treeRef->stats.markNonLeafNodeSearched();
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

    // Populates parentHandles with the path taken to get to the leaf
    template <int min_branch_factor, int max_branch_factor>
    tree_node_handle BranchNode<min_branch_factor, max_branch_factor>::chooseSubtree(
            RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
            tree_node_handle current_handle,
            std::stack<tree_node_handle> &parentHandles,
            const Point &givenPoint
    ) {
      tree_node_handle node_handle = current_handle;

      for (;;) {
        if (node_handle.get_type() == LEAF_NODE) {
          return node_handle;
        }

        auto node = treeRef->get_branch_node(node_handle);

        // Compute the smallest expansion
        unsigned smallestExpansionIndex = 0;
        Branch &b_0 = node->entries.at(0);
        double smallestExpansionArea = b_0.boundingBox.computeExpansionArea(givenPoint);

        for (unsigned i = 1; i < node->cur_offset_ and smallestExpansionArea != -1.0; i++) {
          Branch &b_i = node->entries.at(i);
          double expansionArea = b_i.boundingBox.computeExpansionArea(givenPoint);

          if (expansionArea < smallestExpansionArea) {
            smallestExpansionIndex = i;
            smallestExpansionArea = expansionArea;
          }
        }

        if (smallestExpansionArea != -1.0) {
          Branch &b = node->entries.at(smallestExpansionIndex);
          b.boundingBox.expand(givenPoint);
        }

        // Descend
        // Keep track of the previous handle before descending
        parentHandles.push(node_handle);
        Branch &b = node->entries.at(smallestExpansionIndex);
        node_handle = b.child;
      }

      assert(false);
    }

    template <int min_branch_factor, int max_branch_factor>
    tree_node_handle BranchNode<min_branch_factor, max_branch_factor>::findLeaf(
            RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
            tree_node_handle selfHandle,
            const Point &givenPoint
    ) {
      assert(cur_offset_ > 0);

      for (unsigned i = 0; i < cur_offset_; i++) {
        const Branch &b = entries.at(i);

        if (b.boundingBox.containsPoint(givenPoint)) {
          tree_node_handle ret_handle(nullptr);
          tree_node_handle child_handle = b.child;

          if (child_handle.get_type() == LEAF_NODE) {
            auto child = treeRef->get_leaf_node(child_handle);
            ret_handle = child->findLeaf(treeRef, child_handle, givenPoint);
          } else {
            auto child = treeRef->get_branch_node(child_handle);
            ret_handle = child->findLeaf(treeRef, child_handle, givenPoint);
          }

          if (ret_handle) {
            return ret_handle;
          }

          // Nope, keep looking...
        }
      }

      return tree_node_handle(nullptr);
    }

    template <int min_branch_factor, int max_branch_factor>
    unsigned BranchNode<min_branch_factor, max_branch_factor>::chooseSplitNonLeafAxis() {
      unsigned optimalAxisLower = 0;
      unsigned optimalAxisUpper = 0;
      double optimalMarginLower = std::numeric_limits<double>::infinity();
      double optimalMarginUpper = std::numeric_limits<double>::infinity();

      // Make entries easier to work with
      std::vector<Branch *> lowerEntries;
      lowerEntries.reserve(cur_offset_);
      std::vector<Branch *> upperEntries;
      upperEntries.reserve(cur_offset_);
      for (unsigned i = 0; i < cur_offset_; i++) {
        lowerEntries.push_back(&(entries.at(i)));
        upperEntries.push_back(&(entries.at(i)));
      }

      // Consider all M-2m+2 distributions in each dimension
      for (unsigned d = 0; d < dimensions; d++) {
        // First sort in the current dimension sorting both the lower and upper arrays
        std::sort(lowerEntries.begin(), lowerEntries.end(), [d](Branch *a, Branch *b) {
            return a->boundingBox.lowerLeft[d] < b->boundingBox.lowerLeft[d];
        });
        std::sort(upperEntries.begin(), upperEntries.end(), [d](Branch *a, Branch *b) {
            return a->boundingBox.upperRight[d] < b->boundingBox.upperRight[d];
        });

        // Setup groups
        std::vector<Branch *> groupALower(lowerEntries.begin(),
                                          lowerEntries.begin() + min_branch_factor);
        std::vector<Branch *> groupAUpper(upperEntries.begin(),
                                          upperEntries.begin() + min_branch_factor);

        std::vector<Branch *> groupBLower(lowerEntries.begin() +
                                          min_branch_factor,
                                          lowerEntries.end());
        std::vector<Branch *> groupBUpper(upperEntries.begin() +
                                          min_branch_factor,
                                          upperEntries.end());

        // Cycle through all M-2m+2 distributions
        double totalMarginLower = 0.0;
        double totalMarginUpper = 0.0;
        while (groupALower.size() <= max_branch_factor and groupBLower.size() >= min_branch_factor) {
          // Compute the margin of groupA and groupB
          Rectangle boundingBoxALower = groupALower[0]->boundingBox;
          Rectangle boundingBoxAUpper = groupAUpper[0]->boundingBox;
          for (unsigned i = 1; i < groupALower.size(); i++) {
            boundingBoxALower.expand(groupALower[i]->boundingBox);
            boundingBoxAUpper.expand(groupAUpper[i]->boundingBox);
          }

          Rectangle boundingBoxBLower = groupBLower[0]->boundingBox;
          Rectangle boundingBoxBUpper = groupBUpper[0]->boundingBox;
          for (unsigned i = 1; i < groupBLower.size(); i++) {
            boundingBoxBLower.expand(groupBLower[i]->boundingBox);
            boundingBoxBUpper.expand(groupBUpper[i]->boundingBox);
          }

          // Add to the total margin sum
          totalMarginLower += boundingBoxALower.margin() + boundingBoxBLower.margin();
          totalMarginUpper += boundingBoxAUpper.margin() + boundingBoxBUpper.margin();

          // Add one new value to groupA and remove one from groupB to obtain next distribution
          Branch *transferPointLower = groupBLower.front();
          Branch *transferPointUpper = groupBUpper.front();
          groupBLower.erase(groupBLower.begin());
          groupBUpper.erase(groupBUpper.begin());
          groupALower.push_back(transferPointLower);
          groupAUpper.push_back(transferPointUpper);
        }

        if (totalMarginLower < optimalMarginLower) {
          optimalMarginLower = totalMarginLower;
          optimalAxisLower = d;
        }

        if (totalMarginUpper < optimalMarginUpper) {
          optimalMarginUpper = totalMarginUpper;
          optimalAxisUpper = d;
        }
      }

      bool sortLower = optimalMarginUpper > optimalMarginLower ? true : false;
      unsigned optimalAxis = sortLower ? optimalAxisLower : optimalAxisUpper;

      // Sort to match the optimal axis
      if (sortLower) {

        std::sort(entries.begin(), entries.begin() + cur_offset_,
                  [optimalAxis](Branch &a, Branch &b) {
                      return a.boundingBox.lowerLeft[optimalAxis] <
                             b.boundingBox.lowerLeft[optimalAxis];
                  });
      } else {
        std::sort(entries.begin(), entries.begin() + cur_offset_,
                  [optimalAxis](Branch &a, Branch &b) {
                      return a.boundingBox.upperRight[optimalAxis] <
                             b.boundingBox.upperRight[optimalAxis];
                  });
      }

      return optimalAxis;
    }

// CSA1: Sort entries by lower and upper bound along each axis and compute S -> sum of all
//  margin values for the different distributions. This can be stored in a array of variable
//  that we keep in a loop -> and the just compare to the others?
// 	We can first call a helper function that returns an array of all possible distributions for it?
// CSA2: Return the Axis that has the minimum total sum of all the distributions
    template <int min_branch_factor, int max_branch_factor>
    unsigned BranchNode<min_branch_factor, max_branch_factor>::chooseSplitAxis() {
      return chooseSplitNonLeafAxis();
    }

// CSI1: Given the chosen split index
// 	group all the entries into multiple groups and choose the one that has the least
// 	overlap value; resolve ties with the minimum area
// 	returns tuple of best distribution group indices
    template <int min_branch_factor, int max_branch_factor>
    unsigned BranchNode<min_branch_factor, max_branch_factor>::chooseSplitIndex(unsigned axis) {
      // We assume this is called after we have sorted this->data according to axis.

      const auto groupABegin = entries.begin();
      const auto groupAEnd = entries.begin() + min_branch_factor;
      const auto groupBBegin = entries.begin() + min_branch_factor;
      const auto groupBEnd = entries.begin() + cur_offset_;

      std::vector<Branch> groupA(groupABegin, groupAEnd);
      std::vector<Branch> groupB(groupBBegin, groupBEnd);
      unsigned splitIndex = cur_offset_ / 2;

      // Find the best size out of all the distributions
      double minOverlap = std::numeric_limits<double>::infinity();
      double minArea = std::numeric_limits<double>::infinity();

      // Tracking what the current "cut" mark is
      unsigned currentSplitPoint = min_branch_factor;

      // Try each of the M-2m + 2 groups
      while (groupA.size() <= max_branch_factor and groupB.size() >= min_branch_factor) {
        // Compute the margin of groupA and groupB
        Rectangle boundingBoxA = groupA[0].boundingBox;
        for (unsigned i = 1; i < groupA.size(); i++) {
          boundingBoxA.expand(groupA[i].boundingBox);
        }

        Rectangle boundingBoxB = groupB[0].boundingBox;
        for (unsigned i = 1; i < groupB.size(); i++) {
          boundingBoxB.expand(groupB[i].boundingBox);
        }

        // Compute intersection area to determine best grouping of data points
        double evalDistOverlap = boundingBoxA.computeIntersectionArea(boundingBoxB);

        if (evalDistOverlap < minOverlap) {
          // We save this current distribution of indices to return
          minOverlap = evalDistOverlap;
          splitIndex = currentSplitPoint;

          // Set this if we haven't already
          if (minArea == std::numeric_limits<double>::infinity()) {
            minArea = boundingBoxA.area() + boundingBoxB.area();
          }
        } else if (evalDistOverlap == minOverlap) {
          // If overlap is equal, we use the distribution that creates the smallest areas
          double evalMinArea = boundingBoxA.area() + boundingBoxB.area();

          if (evalMinArea < minArea) {
            // Save this current distribution of indices to return
            minArea = evalMinArea;
            splitIndex = currentSplitPoint;
          }
        }

        // Add one new value to groupA and remove one from groupB to obtain next distribution
        Branch transferPoint = groupB.front();
        groupB.erase(groupB.begin());
        groupA.push_back(transferPoint);

        // Push the split point forward.
        currentSplitPoint++;
      }

      return splitIndex;
    }

    template <int min_branch_factor, int max_branch_factor>
    Partition BranchNode<min_branch_factor, max_branch_factor>::partitionNode() {
      Partition defaultPartition;
      unsigned costMetric = std::numeric_limits<unsigned>::max();
      double location;
      std::vector<Rectangle> sortableBoundingBoxes;

      for (unsigned i = 0; i < cur_offset_; i++) {
        Branch &b = entries.at(i);
        sortableBoundingBoxes.push_back(b.boundingBox);
      }

      for (unsigned d = 0; d < dimensions; d++) {
        // Sort along d
        std::sort(sortableBoundingBoxes.begin(), sortableBoundingBoxes.end(),
                   [d](Rectangle a, Rectangle b){ return a.upperRight[d] < b.upperRight[d]; });


        // By picking a line based on the upper bounding point, we
        // guaranteed that at least some of the entries will go to the
        // left. But we can't guarantee that *all* entries won't go to the
        // left, because we might have to downsplit the remaining
        // entries in the array (R+ nodes are not guaranteed to be disjoint).
        // This would result in a split not actually reducing the number
        // of entries in our new split nodes, which defeats the whole
        // point.
        // I'm not sure we can even guarantee that a line partitions the
        // data in a way such that we DONT overflow --- presumably every
        // box isn't on top of each other but I'm not sure that's a
        // property you rely on.
        // But even if we could guarantee that, we still can't check
        // every cut point to figure it out because that would be N^2
        // and this is D N LOG N.

        location = sortableBoundingBoxes[sortableBoundingBoxes.size() / 2 - 1].upperRight[d];

        // Compute cost, # of splits if d is chosen
        unsigned currentInducedSplits = 0;
        unsigned left_count = 0;
        unsigned right_count = 0;
        for (unsigned i = 0; i < sortableBoundingBoxes.size(); i++) {
          bool is_contained_left = sortableBoundingBoxes[i].upperRight[d] <= location;
          bool is_contained_right = sortableBoundingBoxes[i].lowerLeft[d] >= location;

          if (is_contained_left) {
            left_count++;
          } else if(is_contained_right) {
            right_count++;
          } else if (sortableBoundingBoxes[i].lowerLeft[d] <= location and
                     location <= sortableBoundingBoxes[i].upperRight[d]) {
            currentInducedSplits++;
            left_count++;
            right_count++;
          }
        }

        // Compare cost
        if (left_count <= max_branch_factor and right_count <= max_branch_factor and
            currentInducedSplits < costMetric) {
          defaultPartition.dimension = d;
          defaultPartition.location = location;
          costMetric = currentInducedSplits;
        }
      }

      // If there was a default split point that didnt' overflow children,
      // use that
      if (costMetric < std::numeric_limits<unsigned>::max()) {
        return defaultPartition;
      }

      // It's time to get fancy
      for (unsigned d = 0; d < dimensions; d++) {
        for (unsigned i = 0; i < sortableBoundingBoxes.size(); i++) {
          unsigned left_count = 0;
          unsigned right_count = 0;
          double partition_candidate = sortableBoundingBoxes[i].upperRight[d];
          unsigned cost = 0;

          for (unsigned j = 0; j < sortableBoundingBoxes.size(); j++) {
            Rectangle &bounding_rect = sortableBoundingBoxes[j];
            bool greater_than_left = bounding_rect.lowerLeft[d] < partition_candidate;
            bool less_than_right = bounding_rect.upperRight[d] > partition_candidate;
            bool requires_split = greater_than_left and less_than_right;
            bool should_go_left = bounding_rect.upperRight[d] <= partition_candidate;
            bool should_go_right = bounding_rect.lowerLeft[d] >= partition_candidate;

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
          } //j

          if (left_count <= max_branch_factor and
              right_count <= max_branch_factor and
              left_count > 0 and right_count > 0 ) {
            if (cost < costMetric) {
              defaultPartition.dimension = d;
              defaultPartition.location = partition_candidate;
              costMetric = cost;
            }
          }
        } // i
      } // d

      assert(costMetric < std::numeric_limits<unsigned>::max());

      return defaultPartition;
    }

    template <int min_branch_factor, int max_branch_factor>
    SplitResult BranchNode<min_branch_factor, max_branch_factor>::splitNode(
            RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
            tree_node_handle current_handle,
            Partition p
    ) {
      using NodeType = BranchNode<min_branch_factor, max_branch_factor>;

      tree_node_allocator *allocator = get_node_allocator(treeRef);

      // Create a copy of entries
      std::vector<Branch> entriesCopy(entries);

      // We are the left child
      auto left_handle = current_handle;
      auto left_node = this;
      left_node->cur_offset_ = 0; // Reset our entries and repopulate below

      // Create a new sibling node
      auto alloc_data = allocator->create_new_tree_node<NodeType>(NodeHandleType(BRANCH_NODE));
      auto right_handle = alloc_data.second;
      auto right_node = alloc_data.first;
      new (&(*(right_node))) NodeType();
      right_handle.set_level(current_handle.get_level());

      assert(left_node->cur_offset_ == 0);
      assert(right_node->cur_offset_ == 0);

      // This is partitioning on a point that either shoves everything
      // left or cuts enough stuff that left gets everything.
      // Very confusing...
      for (unsigned i = 0; i < entriesCopy.size(); i++) {
        Branch b = entriesCopy.at(i);

        bool is_contained_left = b.boundingBox.upperRight[p.dimension] <= p.location;
        bool is_contained_right = b.boundingBox.lowerLeft[p.dimension] >= p.location;

        assert(not(is_contained_left and is_contained_right));

        if (is_contained_left) {
          left_node->entries.at(left_node->cur_offset_++) = b;
        } else if (is_contained_right) {
          right_node->entries.at(right_node->cur_offset_++) = b;
        } else { // Line cuts through child, so we downsplit
          SplitResult downwardSplit;

          if (b.child.get_type() == LEAF_NODE) {
            auto child_node = treeRef->get_leaf_node(b.child);
            downwardSplit = child_node->splitNode(treeRef, b.child, p);
          } else {
            auto child_node = treeRef->get_branch_node(b.child);
            downwardSplit = child_node->splitNode(treeRef, b.child, p);
          }

          if (downwardSplit.leftBranch.boundingBox != Rectangle::atInfinity) {
            left_node->entries.at(left_node->cur_offset_++) = downwardSplit.leftBranch;
          }

          if (downwardSplit.rightBranch.boundingBox != Rectangle::atInfinity) {
            right_node->entries.at(right_node->cur_offset_++) = downwardSplit.rightBranch;
          }
        }
      }

      assert(left_node->cur_offset_ <= max_branch_factor);
      assert(right_node->cur_offset_ <= max_branch_factor);
    }

// Overflow treatement for dealing with a node that is too big (overflow)
    template <int min_branch_factor, int max_branch_factor>
    tree_node_handle BranchNode<min_branch_factor, max_branch_factor>::overflowTreatment(
            RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
            tree_node_handle current_handle,
            std::stack<tree_node_handle> parentHandles,
            std::vector<bool> &hasReinsertedOnLevel
    ) {
      uint16_t current_level = current_handle.get_level();
      assert(hasReinsertedOnLevel.size() > current_level);

      if (hasReinsertedOnLevel.at(current_level)) {
        //std::cout << "Overflow treatment on branch node, splitting." <<
        //    std::endl;
//        return splitNode(treeRef, current_handle, partitionNode());
      }
    }

// insert() is always called on the root node.
    template <int min_branch_factor, int max_branch_factor>
    tree_node_handle BranchNode<min_branch_factor, max_branch_factor>::insert(
            RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
            tree_node_handle current_handle,
            Point point
    ) {
      std::vector<bool> hasReinsertedOnLevel;

      tree_node_allocator *allocator = get_node_allocator(treeRef);

      std::stack<tree_node_handle> parentHandles; // Populated by chooseSubtree
      uint16_t current_level = current_handle.get_level();
      assert(current_level > 0); // Branch nodes have level > 0

      // I1 [Find position for new record]
      tree_node_handle insertion_point_handle = chooseSubtree(
              treeRef,
              current_handle,
              parentHandles,
              point
      );
      uint16_t insertion_point_level = insertion_point_handle.get_level();

      tree_node_handle sibling_handle = tree_node_handle(nullptr);

      // I2 [Add record to leaf node]
      auto insertion_point = treeRef->get_leaf_node(insertion_point_handle);
      insertion_point->addPoint(point);

      /*std::cout << "Inserted Point: " << std::get<Point>( point ) <<
              std::endl;
          std::cout << "Insertion point now has points: { " << std::endl;
          for( size_t i = 0; i < insertion_point->cur_offset_; i++ ) {
              std::cout << insertion_point->entries.at(i) << std::endl;
          }
          std::cout << "}" << std::endl;
          */

      unsigned num_els = insertion_point->cur_offset_;

      // If we exceed treeRef->maxBranchFactor we need to do something about it
      if (num_els > max_branch_factor) {
        // We call overflow treatment to determine how our sibling node is treated. If we do a
        // reInsert, sibling is nullptr. This is properly dealt with in adjustTree.
        sibling_handle = insertion_point->overflowTreatment(
                treeRef,
                insertion_point_handle,
                parentHandles,
                hasReinsertedOnLevel
        );
      }

      // I3 [Propogate overflow treatment changes upward]
      sibling_handle = insertion_point->adjustTree(
              treeRef,
              insertion_point_handle,
              sibling_handle,
              parentHandles,
              hasReinsertedOnLevel
      );

      // I4 [Grow tree taller]
      if (sibling_handle) {
        // Sanity check that we're still the root
        assert(treeRef->root == current_handle);

        auto alloc_data =
                allocator->create_new_tree_node<BranchNode<min_branch_factor, max_branch_factor>>(
                        NodeHandleType(BRANCH_NODE));
        auto newRoot = alloc_data.first;
        tree_node_handle new_root_handle = alloc_data.second;
        new_root_handle.set_level(current_level + 1);

        auto sibling = treeRef->get_branch_node(sibling_handle);

        new (&(*(newRoot))) BranchNode<min_branch_factor, max_branch_factor>();

        // Make the existing root a child of newRoot
        Branch b1(boundingBox(), current_handle);
        newRoot->addBranchToNode(b1);

        // Make the new sibling node a child of newRoot
        Branch b2(sibling->boundingBox(), sibling_handle);
        newRoot->addBranchToNode(b2);

        // Ensure newRoot has both children
        assert(newRoot->cur_offset_ == 2);
        assert(sibling_handle.get_level() + 1 == new_root_handle.get_level());

        // Fix the reinserted length
        hasReinsertedOnLevel.push_back(false);

        return new_root_handle;
      }

      return treeRef->root;
    }

// Always called on root, this = root
    template <int min_branch_factor, int max_branch_factor>
    tree_node_handle BranchNode<min_branch_factor, max_branch_factor>::remove(Point &givenPoint) {
#if 0
      assert(!parent);

  // D1 [Find node containing record]
  tree_node_handle leaf_ptr = findLeaf(givenPoint);
  if (!leaf_ptr) {
    return leaf_ptr; /*nullptr*/
  }

  auto leaf = treeRef->get_leaf_node(leaf_ptr);

  // D2 [Delete record]
  leaf->removePoint(givenPoint);

  // D3 [Propagate changes]
  auto root_handle = leaf->condenseTree(hasReinsertedOnLevel);
  if (root_handle.get_type() == BRANCH_NODE) {
    auto root = treeRef->get_branch_node(root_handle);
    if (root->cur_offset_ == 1) {
      // Slice the hasReinsertedOnLevel
      hasReinsertedOnLevel.pop_back();

      // We are removing the root to shorten the tree so we then decide to remove the root
      Branch &b = root->entries[0];

      // Get rid of the old root
      if (b.child.get_type() == LEAF_NODE) {
        auto child = treeRef->get_leaf_node(b.child);
        child->parent = tree_node_handle(nullptr);
      } else {
        auto child = treeRef->get_branch_node(b.child);
        child->parent = tree_node_handle(nullptr);
      }

      // Garbage Collect Root
      // FIXME GC(root);

      return b.child;
    }
  }
  return root_handle;
#endif

      // Unsupported
      abort();
    }

    template <int min_branch_factor, int max_branch_factor>
    void BranchNode<min_branch_factor, max_branch_factor>::print() const {
      std::string indentation(4, ' ');
      std::cout << indentation << "Node " << (void *)this << std::endl;
      std::cout << indentation << "{" << std::endl;
      std::cout << indentation << "    BoundingBox: " << boundingBox() << std::endl;
//  std::cout << indentation << "    Parent: " << parent << std::endl;
      std::cout << indentation << "    Entries: " << std::endl;

      for (unsigned i = 0; i < cur_offset_; i++) {
        const Branch &b = entries.at(i);
        std::cout << indentation << "		" << b.boundingBox << ", ptr: " << b.child << std::endl;
      }
      std::cout << std::endl
                << indentation << "}" << std::endl;
    }

    template <int min_branch_factor, int max_branch_factor>
    void BranchNode<min_branch_factor, max_branch_factor>::printTree() const {
#if 0
      // Print this node first
  struct Printer {
    void operator()(RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef, tree_node_handle node_handle) {
      if (node_handle.get_type() == LEAF_NODE) {
        auto node = treeRef->get_leaf_node(node_handle);
        node->print();
      } else {
        auto node = treeRef->get_branch_node(node_handle);
        node->print();
      }
    }
  };

  Printer p;
  treeWalker<min_branch_factor, max_branch_factor>(treeRef, self_handle_, p);
#endif

      // Unsupported
      abort();
    }

    template <int min_branch_factor, int max_branch_factor>
    void printPackedNodes(
            RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
            tree_node_handle node_handle,
            std::ofstream &printFile
    ) {
#if 0
      tree_node_allocator *allocator = treeRef->node_allocator_.get();
  auto node = allocator->get_tree_node<packed_node>(node_handle);
  char *data = node->buffer_;
  decode_entry_count_and_offset_packed_node(data);

  if (node_handle.get_type() == REPACKED_BRANCH_NODE) {
    for (size_t i = 0; i < count; i++) {
      Branch *b = (Branch *)(data + offset);

      // Only print the last layer of branch nodes, this is where intersection
      // happens the most.
      if (b->child.get_type() == REPACKED_LEAF_NODE) {
        printFile << b->boundingBox << std::endl;
      }

      offset += sizeof(Branch);
    }
  } else if (node_handle.get_type() == REPACKED_LEAF_NODE) {
    for (size_t i = 0; i < count; i++) {
      Point *p = (Point *)(data + offset);
      printFile << *p << std::endl;
      offset += sizeof(Point);
    }
  }
#endif
    }

    template <int min_branch_factor, int max_branch_factor>
    unsigned BranchNode<min_branch_factor, max_branch_factor>::checksum() const {
#if 0
      struct ChecksumFunctor {
      unsigned checksum;

      ChecksumFunctor() {
        checksum = 0;
      }

      void operator()(
          RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
          tree_node_handle node_handle)
      {
        if (node_handle.get_type() == LEAF_NODE) {
          auto node = treeRef->get_leaf_node(node_handle);
          for (unsigned i = 0; i < node->cur_offset_; i++) {
            const Point &p = node->entries.at(i);
            for (unsigned d = 0; d < dimensions; d++) {
              checksum += (unsigned)p[d];
            }
          }
        }
      }
  };

  ChecksumFunctor cf;
  treeWalker<min_branch_factor, max_branch_factor>(treeRef, self_handle_, cf);
  return cf.checksum;
#endif

      // Unsupported
      abort();
    }

    template <int min_branch_factor, int max_branch_factor>
    unsigned BranchNode<min_branch_factor, max_branch_factor>::height(tree_node_handle self_handle) const {
      return self_handle.get_level() + 1;
    }

    template <int min_branch_factor, int max_branch_factor>
    void stat_node(tree_node_handle root_handle, RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef) {
      std::stack<tree_node_handle> context;

      context.push(root_handle);
      size_t memoryFootprint = 0;
      unsigned long totalNodes = 0;
      unsigned long totalLeaves = 0;
      unsigned treeHeight;
      size_t deadSpace = 0;
      // histogramFanoutAtLevel stores count of branch with different fanout at each level
      // Eg. histogramFanoutAtLevel.at(1).at(2) stores count of branch with fanout=2 at level 1
      std::vector<std::vector<unsigned long>> histogramFanoutAtLevel;

      tree_node_allocator *allocator = get_node_allocator(treeRef);

      if (root_handle.get_type() == LEAF_NODE) {
        auto root_node = treeRef->get_leaf_node(root_handle);
        treeHeight = root_node->height();
      } else {
        assert(root_handle.get_type() == BRANCH_NODE);
        auto root_node = treeRef->get_branch_node(root_handle);
        treeHeight = root_node->height(root_handle);
      }

      histogramFanoutAtLevel.resize(treeHeight);
      for (unsigned lvl = 0; lvl < treeHeight; lvl++) {
        histogramFanoutAtLevel.at(lvl).resize(10000,0);
      }

      while (!context.empty()) {
        auto currentContext = context.top();
        context.pop();
        totalNodes++;
        auto lvl = currentContext.get_level();

        if (currentContext.get_type() == LEAF_NODE) {
          auto current_node = treeRef->get_leaf_node(currentContext);
          unsigned fanout = current_node->cur_offset_;

          if (fanout >= histogramFanoutAtLevel.at(lvl).size()) {
            histogramFanoutAtLevel.at(lvl).resize(2 * fanout, 0);
          }

          histogramFanoutAtLevel.at(lvl).at(fanout)++;
          memoryFootprint += sizeof(LeafNode<min_branch_factor, max_branch_factor>);
          deadSpace += (sizeof(Point) * (max_branch_factor - current_node->cur_offset_));
          totalLeaves++;
        } else if (currentContext.get_type() == BRANCH_NODE) {
          auto current_node = treeRef->get_branch_node(currentContext);
          unsigned fanout = current_node->cur_offset_;

          if (fanout >= histogramFanoutAtLevel.at(lvl).size()) {
            histogramFanoutAtLevel.at(lvl).resize(2 * fanout, 0);
          }

          histogramFanoutAtLevel.at(lvl).at(fanout)++;
          memoryFootprint += sizeof(BranchNode<min_branch_factor, max_branch_factor>);
          deadSpace += (sizeof(Branch) * (max_branch_factor - current_node->cur_offset_));

          for (unsigned i = 0; i < current_node->cur_offset_; i++) {
            context.push(current_node->entries.at(i).child);
          }
        }
      }

      // Print out what we have found
      STATEXEC(std::cout << "### Statistics ###" << std::endl);
      // Memory footprint is wrong!
      STATMEM(memoryFootprint);
      //STATHEIGHT(height());
      STATSIZE(totalNodes);
      STATEXEC(std::cout << "DeadSpace: " << deadSpace << std::endl);
      //STATSINGULAR(singularBranches);
      STATLEAF(totalLeaves);
      printFanoutHistogram(histogramFanoutAtLevel, treeHeight);

      std::cout << treeRef->stats << std::endl;
    }
} // namespace rplustreedisk
