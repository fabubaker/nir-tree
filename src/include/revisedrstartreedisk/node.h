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

namespace revisedrstartreedisk {
    template <int min_branch_factor, int max_branch_factor> class RevisedRStarTreeDisk;

    template <int min_branch_factor, int max_branch_factor>
    tree_node_allocator *get_node_allocator(RevisedRStarTreeDisk<min_branch_factor, max_branch_factor> *treeRef) {
      return treeRef->node_allocator_.get();
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
                RevisedRStarTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
                tree_node_handle selfHandle,
                const Point &givenPoint
        );
        unsigned chooseSplitLeafAxis();
        unsigned chooseSplitNonLeafAxis();
        unsigned chooseSplitAxis();
        unsigned chooseSplitIndex(unsigned axis, double s);
        SplitResult splitNode(
                RevisedRStarTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
                tree_node_handle current_handle
        );
        tree_node_handle condenseTree();
        Partition partitionNode();

        // Datastructure interface functions
        void exhaustiveSearch(const Point &requestedPoint, std::vector<Point> &accumulator) const;

        // These return the root of the tree.
        tree_node_handle insert(
                RevisedRStarTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
                tree_node_handle current_handle,
                Point point
        );
        tree_node_handle remove(Point &givenPoint);

        // Miscellaneous
        unsigned checksum() const;
        void print() const;
        unsigned height() const;

        // Operators
        bool operator<(const LeafNode &otherNode) const;
    };

    template <int min_branch_factor, int max_branch_factor>
    class BranchNode {
    public:
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
        void removeChild(tree_node_handle child);
        void chooseNodeHelper(
                unsigned limitIndex, 
                Point &givenPoint, 
                unsigned &chosenIndex, 
                bool &success, 
                std::vector<bool> &candidates, 
                std::vector<double> &deltas, 
                unsigned startIndex, 
                bool useMarginDelta);
        tree_node_handle chooseSubtree(
                RevisedRStarTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
                tree_node_handle current_handle,
                std::stack<tree_node_handle> &parentHandles,
                const Point &givenPoint
        );
        tree_node_handle findLeaf(
                RevisedRStarTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
                tree_node_handle selfHandle,
                const Point &givenPoint
        );
        unsigned chooseSplitLeafAxis();
        unsigned chooseSplitNonLeafAxis();
        unsigned chooseSplitAxis();
        unsigned chooseSplitIndex(unsigned axis, double s);
        SplitResult splitNode(
                RevisedRStarTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
                tree_node_handle current_handle
        );
        tree_node_handle condenseTree();
        Partition partitionNode();

        // Datastructure interface functions
        void exhaustiveSearch(
          const Point &requestedPoint, std::vector<Point> &accumulator, RevisedRStarTreeDisk<min_branch_factor, max_branch_factor> *treeRef
        ) const;

        // These return the root of the tree.
        tree_node_handle insert(
                RevisedRStarTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
                tree_node_handle current_handle,
                Point point
        );
        tree_node_handle remove(Point &givenPoint);

        // Miscellaneous
        unsigned checksum() const;
        void print() const;
        unsigned height(tree_node_handle self_handle) const;
    };

    template <int min_branch_factor, int max_branch_factor, typename functor>
    void treeWalker(RevisedRStarTreeDisk<min_branch_factor, max_branch_factor> *treeRef, tree_node_handle root, functor &f) {
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
            RevisedRStarTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
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

          // get sum of margin on either side of the split
          double evalMargin = boundingBoxA.margin() + boundingBoxB.margin();
          if (evalMargin < optimalMargin) {
            optimalMargin = evalMargin;
            optimalAxis = d;
          }
          // Add one new value to groupA and remove one from groupB to obtain next distribution
          Point *transferPoint = groupB.front();
          groupB.erase(groupB.begin());
          groupA.push_back(transferPoint);
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
    unsigned LeafNode<min_branch_factor, max_branch_factor>::chooseSplitIndex(unsigned axis, double s) {
      // We assume this is called after we have sorted this->data according to axis.
      // Precompute the elements not dependant on candidateIndex
      Rectangle bb = this->boundingBox();
      double m = (double) min_branch_factor;
      double M = (double) max_branch_factor;
      double asym = 2 * (bb.centrePoint()[axis] - originalCentre[axis]) / (std::fabs(bb.upperRight[axis] - bb.lowerLeft[axis]));
      double u = (1 - (2 * m) / (M + 1)) * asym;
      double sigma = treeRef.s * (1 + std::fabs(u));
      double y1 = std::exp(-1 / std::pow(s, 2));
      double ys = 1 / (1 - y1);

      const auto groupABegin = entries.begin();
      const auto groupAEnd = entries.begin() + min_branch_factor;
      const auto groupBBegin = entries.begin() + min_branch_factor;
      const auto groupBEnd = entries.begin() + cur_offset_;

      std::vector<Point> groupA(groupABegin, groupAEnd);
      std::vector<Point> groupB(groupBBegin, groupBEnd);
      unsigned splitIndex = cur_offset_ / 2;

      // Find the best size out of all the distributions
      double minWeight = std::numeric_limits<double>::infinity();

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

        // Evaluate wf(i)
        double xi = ((2 * currentSplitPoint) / (M + 1)) - 1;
        double weightFunction = ys * (std::exp(-std::pow((xi - u) / sigma, 2)) - y1);

        // Compute intersection area to determine best grouping of data points
        double evalDistOverlap = boundingBoxA.computeIntersectionArea(boundingBoxB);
        
        // Evaluate wg(i)
        double evalWeight;
        if (evalDistOverlap == 0.0) {
          double weightGoal = boundingBoxA.margin() + boundingBoxB.margin() - bb.margin();
          evalWeight = weightGoal * weightFunction;
        } else {
          double weightGoal = evalDistOverlap;
          evalWeight = weightGoal / weightFunction;
        }

        if (evalWeight < minWeight) {
          splitIndex = currentSplitPoint;
          minWeight = evalWeight;
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
    SplitResult LeafNode<min_branch_factor, max_branch_factor>::splitNode(
      RevisedRStarTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
      tree_node_handle current_handle ) {
      using NodeType = LeafNode<min_branch_factor, max_branch_factor>;
      // S1: Call chooseSplitAxis to determine the axis perpendicular to which the split is performed
      // S2: Invoke chooseSplitIndex given the axis to determine the best distribution along this axis
      // S3: Distribute the entries among these two groups

      // Call chooseSplitAxis to determine the axis perpendicular to which the split is performed
      // For now we will save the axis as an int -> since this allows for room for growth in the future
      // Call ChooseSplitIndex to create optimal splitting of data array
      unsigned splitAxis = chooseSplitAxis();
      unsigned splitIndex = chooseSplitIndex(splitAxis, treeRef->s);

      // We are the left child
      auto left_handle = current_handle;
      auto left_node = this;   

      // Create a new sibling
      tree_node_allocator *allocator = get_node_allocator(treeRef);
      auto alloc_data = allocator->create_new_tree_node<NodeType>(NodeHandleType(LEAF_NODE));
      new (&(*(alloc_data.first))) NodeType();
      auto right_handle = alloc_data.second;
      auto right_node = alloc_data.first;
      right_handle.set_level(left_handle.get_level());
	
      // Copy everything to the right of the splitPoint (inclusive) to the new sibling
      std::copy(entries.begin() + splitIndex, entries.begin() + cur_offset_, right_node->entries.begin());
      right_node->cur_offset_ = cur_offset_ - splitIndex;

      // Chop our node's data down
      cur_offset_ = splitIndex;

      assert(left_node->cur_offset_ <= max_branch_factor);
      assert(right_node->cur_offset_ <= max_branch_factor);

      SplitResult split = {
        {left_node->boundingBox(), left_handle},
        {right_node->boundingBox(), right_handle}
      };

      return split;
    }

    template <int min_branch_factor, int max_branch_factor>
    SplitResult propagateSplit(
      RevisedRStarTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
      pinned_node_ptr<BranchNode<min_branch_factor, max_branch_factor>> current_node,
      tree_node_handle current_handle,
      tree_node_handle parent_handle,
      SplitResult currentPropagationSplit
    ) {
      // Returns a new split that needs to be propagated upwards
      SplitResult newPropagationSplit;

      // If there was a split we were supposed to propagate then propagate it
      if (currentPropagationSplit.leftBranch.child != nullptr and currentPropagationSplit.rightBranch.child != nullptr) {
        if (current_node->cur_offset_ >= max_branch_factor) {
          assert(false);
        }

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
      newPropagationSplit = current_node->splitNode(treeRef, current_handle);

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
    SplitResult propagateSplit(
      RevisedRStarTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
      pinned_node_ptr<LeafNode<min_branch_factor, max_branch_factor>> current_node,
      tree_node_handle current_handle,
      tree_node_handle parent_handle,
      SplitResult currentPropagationSplit
    ) {
      // Returns a new split that needs to be propagated upwards
      SplitResult newPropagationSplit;

      // Sanity check that there is no split to propagate on us, a leaf node
      assert(currentPropagationSplit.leftBranch.child == nullptr and
             currentPropagationSplit.rightBranch.child == nullptr);

      // Early exit if this node does not overflow
      if (current_node->cur_offset_ <= max_branch_factor) {
        newPropagationSplit = {
                {Rectangle(), tree_node_handle( nullptr )},
                {Rectangle(), tree_node_handle( nullptr )}
        };
        return newPropagationSplit;
      }

      // Otherwise, split node
      newPropagationSplit = current_node->splitNode(treeRef, current_handle);

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
    SplitResult adjustTree(
      RevisedRStarTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
      tree_node_handle current_handle,
      std::stack<tree_node_handle> parentHandles
    ) {
      SplitResult currentPropagationSplit = {
        {Rectangle(), tree_node_handle( nullptr )},
        {Rectangle(), tree_node_handle( nullptr )}
      };

      for (;;) {
        if (!current_handle) {
          break;
        }

        tree_node_handle parent_handle;

        if (parentHandles.empty()) {
          parent_handle = tree_node_handle(nullptr);
        } else {
          parent_handle = parentHandles.top();
          parentHandles.pop();
        }

        if (current_handle.get_type() == LEAF_NODE) {
          auto current_node = treeRef->get_leaf_node(current_handle);

          currentPropagationSplit = propagateSplit(
                  treeRef, current_node, current_handle, parent_handle, currentPropagationSplit
          );

          // Stop adjusting tree if there are no more splits to propagate
          if (currentPropagationSplit.leftBranch.child == nullptr and
              currentPropagationSplit.rightBranch.child == nullptr) {
            return currentPropagationSplit;
          }

          // Ascend
          current_handle = parent_handle;
        } else {
          auto current_node = treeRef->get_branch_node(current_handle);

          currentPropagationSplit = propagateSplit(
              treeRef, current_node, current_handle, parent_handle, currentPropagationSplit
          );

          // Stop adjusting tree if there are no more splits to propagate
          if (currentPropagationSplit.leftBranch.child == nullptr and
              currentPropagationSplit.rightBranch.child == nullptr) {
            return currentPropagationSplit;
          }

          // Ascend
          current_handle = parent_handle;
        }
      }

      // If the root has been split, this should contain the final split to grow the tree.
      return currentPropagationSplit;
    }

// insert() is always called on the root node. If this gets called,
// that means the tree only has one leaf node.
    template <int min_branch_factor, int max_branch_factor>
    tree_node_handle LeafNode<min_branch_factor, max_branch_factor>::insert(
            RevisedRStarTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
            tree_node_handle current_handle,
            Point point
    ) {
      tree_node_allocator *allocator = get_node_allocator(treeRef);

      // Empty parentHandles since we are the only node in the tree.
      std::stack<tree_node_handle> parentHandles;

      uint16_t current_level = current_handle.get_level();
      assert(current_level == 0); // Leaf nodes have level = 0

      addPoint(point);

      // If we exceed max_branch_factor we need to do something about it
      if (cur_offset_ > max_branch_factor) {
        SplitResult split;

        // Split ourselves where the new point was inserted
        split = splitNode(treeRef, current_handle);

        // Sanity check that we're still the root
        assert(treeRef->root == current_handle);
        // Sanity check that the split produced two nodes
        assert(split.leftBranch.child != nullptr and split.rightBranch.child != nullptr);

        // We need a new root
        auto node_type = NodeHandleType(BRANCH_NODE);
        auto alloc_data =
                allocator->create_new_tree_node<BranchNode<min_branch_factor, max_branch_factor>>(node_type);
        auto new_root_node = alloc_data.first;
        auto new_root_handle = alloc_data.second;
        new_root_handle.set_level(current_level + 1);

        new (&(*(new_root_node))) BranchNode<min_branch_factor, max_branch_factor>();

        new_root_node->entries.at(new_root_node->cur_offset_++) = split.leftBranch;
        new_root_node->entries.at(new_root_node->cur_offset_++) = split.rightBranch;

        treeRef->root = new_root_handle;
        return new_root_handle;
      }

      return treeRef->root;
    }

// To be called on a leaf
    template <int min_branch_factor, int max_branch_factor>
    tree_node_handle LeafNode<min_branch_factor, max_branch_factor>::condenseTree() {
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
    void BranchNode<min_branch_factor, max_branch_factor>::exhaustiveSearch(
        const Point &requestedPoint, std::vector<Point> &accumulator,
        RevisedRStarTreeDisk<min_branch_factor, max_branch_factor> *treeRef
    ) const {
      for (unsigned i = 0; i < cur_offset_; i++) {
        tree_node_handle child_handle = entries.at(i).child;

        if (child_handle.get_type() == LEAF_NODE) {
          auto child = treeRef->get_leaf_node(child_handle);
          child->exhaustiveSearch(requestedPoint, accumulator);
        } else {
          auto child = treeRef->get_branch_node(child_handle);
          child->exhaustiveSearch(requestedPoint, accumulator, treeRef);
        }
      }
    }

    template <int min_branch_factor, int max_branch_factor>
    void point_search_leaf_node(LeafNode<min_branch_factor, max_branch_factor> &node,
                                Point &requestedPoint,
                                std::vector<Point> &accumulator,
                                RevisedRStarTreeDisk<min_branch_factor, max_branch_factor> *treeRef)
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
                                  RevisedRStarTreeDisk<min_branch_factor, max_branch_factor> *treeRef)
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
            RevisedRStarTreeDisk<min_branch_factor, max_branch_factor> *treeRef)
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
                                    RevisedRStarTreeDisk<min_branch_factor, max_branch_factor> *treeRef)
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
                                      RevisedRStarTreeDisk<min_branch_factor, max_branch_factor> *treeRef)
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
            RevisedRStarTreeDisk<min_branch_factor, max_branch_factor> *treeRef)
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
    	
    template <int min_branch_factor, int max_branch_factor>
    void BranchNode<min_branch_factor, max_branch_factor>::chooseNodeHelper(unsigned limitIndex, 
                                                                            Point &givenPoint, 
                                                                            unsigned &chosenIndex, 
                                                                            bool &success, 
                                                                            std::vector<bool> &candidates, 
                                                                            std::vector<double> &deltas, 
                                                                            unsigned startIndex, 
                                                                            bool useMarginDelta)
    {
      candidates[startIndex] = true;
      deltas[startIndex] = 0.0;

      for (unsigned j = 0; j <= limitIndex; ++j)
      {
        if (j != startIndex)
        {
          unsigned additionalDelta = useMarginDelta ? entries.at(startIndex).boundingBox.marginDelta(givenPoint, entries.at(j).boundingBox) : entries.at(startIndex).boundingBox.areaDelta(givenPoint, entries.at(j0f32x).boundingBox);
          deltas[startIndex] += additionalDelta;

          if (additionalDelta != 0.0 && !candidates[j])
          {
            this->chooseNodeHelper(limitIndex, givenPoint, chosenIndex, success, candidates, deltas, j, useMarginDelta);
            if (success)
            {
              break;
            }
          }
        }
      }

      if (deltas[startIndex] == 0.0)
      {
        chosenIndex = startIndex;
        success = true;
      }
    }

    // Populates parentHandles with the path taken to get to the leaf
    template <int min_branch_factor, int max_branch_factor>
    tree_node_handle BranchNode<min_branch_factor, max_branch_factor>::chooseSubtree(
            RevisedRStarTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
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

				unsigned optimalBranchIndex = 0;

				// Find all rectangles that completely cover the given point
				std::vector<unsigned> covers(0);
				for (unsigned i = 0; i < node->cur_offset_; ++i) {
          Branch &b_i = node->entries.at(i);
					if (b_i.boundingBox.containsPoint(givenPoint)) {
						covers.push_back(i);
					}
				}

				if (!covers.empty()) {
					// If any rectangles cover the given point select the lowest volume, then lowest
					// margin rectangle among them
					double minVolume = std::numeric_limits<double>::infinity();
					double minMargin = std::numeric_limits<double>::infinity();
					unsigned minIndex = covers.front();

					for (unsigned i : covers) {
						double evalVolume = node->entries.at(i).boundingBox.computeExpansionArea(givenPoint);
            double evalMargin = node->entries.at(i).boundingBox.computeExpansionMargin(givenPoint);
						if (evalVolume < minVolume) {
							minVolume = evalVolume;
							minMargin = evalMargin;
							minIndex = i;
						} else if (evalVolume == minVolume && evalVolume == 0) {
							// Tie break using perimeter
							if (evalMargin < minMargin) {
								minVolume = evalVolume;
								minMargin = evalMargin;
								minIndex = i;
							}
						}
					}

					optimalBranchIndex = minIndex;
				} else {
          // Make the entries easier to work with
          // std::vector<Branch *> entriesCopy;
          // entriesCopy.reserve(node->cur_offset_);
          // for (unsigned i = 0; i < node->cur_offset_; i++) {
          //   entriesCopy.push_back(&(node->entries.at(i)));
          // }
					// Sort the entries in ascending order of their margin delta
					std::sort(node->entries.begin(), node->entries.begin() + node->cur_offset_, [givenPoint](Branch &a, Branch &b){
            return a.boundingBox.computeExpansionMargin(givenPoint) <= b.boundingBox.computeExpansionMargin(givenPoint);});

					// Look at the first entry's intersection margin with all the others
					double deltaWithAll = 0.0;
					for (Branch &branch : node->entries) {
						deltaWithAll += branch.boundingBox.marginDelta(givenPoint, node->entries[0].boundingBox);
					}

					if (deltaWithAll == 0.0) {
						optimalBranchIndex = 0;
					} else {
						// Set limitIndex based on margin deltas that are not 0
						unsigned limitIndex = 0;
						double maxMarginDelta = - std::numeric_limits<double>::infinity();

						for (unsigned i = 1; i < node->cur_offset_; ++i) {
							double evalMarginDelta = node->entries[0].boundingBox.marginDelta(givenPoint, node->entries[i].boundingBox);
							if (evalMarginDelta > maxMarginDelta) {
								maxMarginDelta = evalMarginDelta;
								limitIndex = i;
							}
						}

						// Consider branches only up to limitIndex
						std::vector<bool> candidate(limitIndex + 1, false);
						std::vector<double> deltas(limitIndex + 1, 0.0);
						bool success = false;
						unsigned chosenIndex;

						// Determine if there exists a rectangle with zero area containing given point
						bool zeroAreaContainer = false;
						for (unsigned i = 0; !zeroAreaContainer && i <= limitIndex; ++i) {
							zeroAreaContainer = 0.0 == node->entries[i].boundingBox.copyExpand(givenPoint).area();
						}

						if (zeroAreaContainer) {
							node->chooseNodeHelper(limitIndex, givenPoint, chosenIndex, success, candidate, deltas, 0, true);
						} else {
							node->chooseNodeHelper(limitIndex, givenPoint, chosenIndex, success, candidate, deltas, 0, false);
						}

						if (success) {
							optimalBranchIndex = chosenIndex;
						} else {
							double minDelta = std::numeric_limits<double>::infinity();

							for (unsigned i = 0; i <= limitIndex; ++i) {
								if (deltas[i] < minDelta && candidate[i]) {
									minDelta = deltas[i];
									optimalBranchIndex = i;
								}
							}
						}
					}
				}

				// Descend
        parentHandles.push(node_handle);
        Branch &b = node->entries.at(optimalBranchIndex);
        b.boundingBox.expand(givenPoint);
        node_handle = b.child;
      }        
      //   // Compute the smallest expansion
      //   unsigned smallestExpansionIndex = 0;
      //   Branch &b_0 = node->entries.at(0);
      //   double smallestExpansionArea = b_0.boundingBox.computeExpansionArea(givenPoint);

      //   for (unsigned i = 1; i < node->cur_offset_ and smallestExpansionArea != -1.0; i++) {
      //     Branch &b_i = node->entries.at(i);
      //     double expansionArea = b_i.boundingBox.computeExpansionArea(givenPoint);

      //     if (expansionArea < smallestExpansionArea) {
      //       smallestExpansionIndex = i;
      //       smallestExpansionArea = expansionArea;
      //     }
      //   }

      //   if (smallestExpansionArea != -1.0) {
      //     Branch &b = node->entries.at(smallestExpansionIndex);
      //     b.boundingBox.expand(givenPoint);
      //   }

      //   // Descend
      //   // Keep track of the previous handle before descending
      //   parentHandles.push(node_handle);
      //   Branch &b = node->entries.at(smallestExpansionIndex);
      //   node_handle = b.child;
      // }

      assert(false);
    }

    template <int min_branch_factor, int max_branch_factor>
    tree_node_handle BranchNode<min_branch_factor, max_branch_factor>::findLeaf(
            RevisedRStarTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
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
          double evalMarginLower = boundingBoxALower.margin() + boundingBoxBLower.margin();
          double evalMarginUpper = boundingBoxAUpper.margin() + boundingBoxBUpper.margin();
          if (evalMarginLower < optimalMarginLower) {
            optimalMarginLower = evalMarginLower;
            optimalAxisLower = d;
          }

          if (evalMarginUpper < optimalMarginUpper) {
            optimalMarginUpper = evalMarginUpper;
            optimalAxisUpper = d;
          }
          // Add one new value to groupA and remove one from groupB to obtain next distribution
          Branch *transferPointLower = groupBLower.front();
          Branch *transferPointUpper = groupBUpper.front();
          groupBLower.erase(groupBLower.begin());
          groupBUpper.erase(groupBUpper.begin());
          groupALower.push_back(transferPointLower);
          groupAUpper.push_back(transferPointUpper);
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
    unsigned BranchNode<min_branch_factor, max_branch_factor>::chooseSplitIndex(unsigned axis, double s) {
      // We assume this is called after we have sorted this->data according to axis.
      // Precompute the elements not dependant on candidateIndex
      Rectangle bb = this->boundingBox();
      double m = (double) min_branch_factor;
      double M = (double) max_branch_factor;
      double asym = 2 * (bb.centrePoint()[axis] - originalCentre[axis]) / (std::fabs(bb.upperRight[axis] - bb.lowerLeft[axis]));
      double u = (1 - (2 * m) / (M + 1)) * asym;
      double sigma = treeRef.s * (1 + std::fabs(u));
      double y1 = std::exp(-1 / std::pow(s, 2));
      double ys = 1 / (1 - y1);

      const auto groupABegin = entries.begin();
      const auto groupAEnd = entries.begin() + min_branch_factor;
      const auto groupBBegin = entries.begin() + min_branch_factor;
      const auto groupBEnd = entries.begin() + cur_offset_;

      std::vector<Branch> groupA(groupABegin, groupAEnd);
      std::vector<Branch> groupB(groupBBegin, groupBEnd);
      unsigned splitIndex = cur_offset_ / 2;

      // Find the best size out of all the distributions
      double minWeight = std::numeric_limits<double>::infinity();

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

        // Evaluate wf(i)
        double xi = ((2 * currentSplitPoint) / (M + 1)) - 1;
        double weightFunction = ys * (std::exp(-std::pow((xi - u) / sigma, 2)) - y1);

        // Compute intersection area to determine best grouping of data points
        double evalDistOverlap = boundingBoxA.computeIntersectionArea(boundingBoxB);
        
        // Evaluate wg(i)
        double evalWeight;
        if (evalDistOverlap == 0.0) {
          double weightGoal = boundingBoxA.margin() + boundingBoxB.margin() - bb.margin();
          evalWeight = weightGoal * weightFunction;
        } else {
          double weightGoal = evalDistOverlap;
          evalWeight = weightGoal / weightFunction;
        }

        if (evalWeight < minWeight) {
          splitIndex = currentSplitPoint;
          minWeight = evalWeight;
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
    SplitResult BranchNode<min_branch_factor, max_branch_factor>::splitNode(
            RevisedRStarTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
            tree_node_handle current_handle
    ) {
      using NodeType = BranchNode<min_branch_factor, max_branch_factor>;

      // S1: Call chooseSplitAxis to determine the axis perpendicular to which the split is performed
      // S2: Invoke chooseSplitIndex given the axis to determine the best distribution along this axis
      // S3: Distribute the entries among these two groups

      // Call chooseSplitAxis to determine the axis perpendicular to which the split is performed
      // For now we will save the axis as a int -> since this allows for room for growth in the future
      // Call ChooseSplitIndex to create optimal splitting of data array
      unsigned splitAxis = chooseSplitAxis();
      unsigned splitIndex = chooseSplitIndex(splitAxis, treeRef->s);

      // We are the left child
      auto left_handle = current_handle;
      auto left_node = this;

      // Create a new sibling node
      tree_node_allocator *allocator = get_node_allocator(treeRef);
      auto alloc_data = allocator->create_new_tree_node<NodeType>(NodeHandleType(BRANCH_NODE));
      auto right_handle = alloc_data.second;
      auto right_node = alloc_data.first;
      new (&(*(right_node))) NodeType();
      right_handle.set_level(current_handle.get_level());

      // Copy everything to the right of the splitPoint (inclusive) to the new sibling
      std::copy(entries.begin() + splitIndex, entries.begin() + cur_offset_, right_node->entries.begin());
      right_node->cur_offset_ = cur_offset_ - splitIndex;

      // Chop our node's data down
      cur_offset_ = splitIndex;

      assert(left_node->cur_offset_ <= max_branch_factor);
      assert(right_node->cur_offset_ <= max_branch_factor);

      SplitResult split = {
              {left_node->boundingBox(), left_handle},
              {right_node->boundingBox(), right_handle}
      };

      return split;
    }

// insert() is always called on the root node.
    template <int min_branch_factor, int max_branch_factor>
    tree_node_handle BranchNode<min_branch_factor, max_branch_factor>::insert(
            RevisedRStarTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
            tree_node_handle current_handle,
            Point point
    ) {
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

      // I2 [Add record to leaf node]
      auto insertion_point = treeRef->get_leaf_node(insertion_point_handle);
      insertion_point->addPoint(point);

      // If we exceed max_branch_factor we need to do something about it
      if (insertion_point->cur_offset_ > max_branch_factor) {
        SplitResult finalSplit;

        // Adjust the tree from bottom-to-top, splitting nodes if they are full
        finalSplit = adjustTree(treeRef, insertion_point_handle, parentHandles);

        // If the root was split in adjustTree, create a new root
        if (finalSplit.leftBranch.child != nullptr and finalSplit.rightBranch.child != nullptr) {
          // Sanity check that we're still the root
          assert(treeRef->root == current_handle);

          auto node_type = NodeHandleType(BRANCH_NODE);
          auto alloc_data =
                  allocator->create_new_tree_node<BranchNode<min_branch_factor, max_branch_factor>>(node_type);
          auto new_root_node = alloc_data.first;
          auto new_root_handle = alloc_data.second;
          new_root_handle.set_level(current_level + 1);

          new (&(*(new_root_node))) BranchNode<min_branch_factor, max_branch_factor>();

          new_root_node->entries.at(new_root_node->cur_offset_++) = finalSplit.leftBranch;
          new_root_node->entries.at(new_root_node->cur_offset_++) = finalSplit.rightBranch;

          treeRef->root = new_root_handle;
          return new_root_handle;
        }
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
      std::cout << indentation << "    Entries: " << std::endl;

      for (unsigned i = 0; i < cur_offset_; i++) {
        const Branch &b = entries.at(i);
        std::cout << indentation << "		" << b.boundingBox << ", ptr: " << b.child << std::endl;
      }
      std::cout << std::endl
                << indentation << "}" << std::endl;
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
          RevisedRStarTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
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
    void stat_node(tree_node_handle root_handle, RevisedRStarTreeDisk<min_branch_factor, max_branch_factor> *treeRef) {
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
        histogramFanoutAtLevel.at(lvl).resize(max_branch_factor + 1, 0);
      }

      while (!context.empty()) {
        auto currentContext = context.top();
        context.pop();
        totalNodes++;
        auto lvl = currentContext.get_level();

        if (currentContext.get_type() == LEAF_NODE) {
          auto current_node = treeRef->get_leaf_node(currentContext);
          unsigned fanout = current_node->cur_offset_;

          histogramFanoutAtLevel.at(lvl).at(fanout)++;
          memoryFootprint += sizeof(LeafNode<min_branch_factor, max_branch_factor>);
          deadSpace += (sizeof(Point) * (max_branch_factor - current_node->cur_offset_));
          totalLeaves++;
        } else if (currentContext.get_type() == BRANCH_NODE) {
          auto current_node = treeRef->get_branch_node(currentContext);
          unsigned fanout = current_node->cur_offset_;

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
      STATEXEC(std::cout << "Tree ");
      STATMEM(memoryFootprint);
      //STATHEIGHT(height());
      STATSIZE(totalNodes);
      STATEXEC(std::cout << "DeadSpace: " << deadSpace << std::endl);
      //STATSINGULAR(singularBranches);
      STATLEAF(totalLeaves);
      printFanoutHistogram(histogramFanoutAtLevel, treeHeight);

      std::cout << treeRef->stats << std::endl;
    }
} // namespace revisedrstartreedisk
