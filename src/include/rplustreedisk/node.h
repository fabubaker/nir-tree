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
        SplitResult splitNode(
                RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
                tree_node_handle current_handle,
                Partition p
        );
        tree_node_handle condenseTree();
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
        SplitResult splitNode(
                RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
                tree_node_handle current_handle,
                Partition p
        );
        tree_node_handle condenseTree();
        Partition partitionNode();

        // Datastructure interface functions
        void exhaustiveSearch(
          const Point &requestedPoint, std::vector<Point> &accumulator, RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef
        ) const;

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
        unsigned height(tree_node_handle self_handle) const;
    };

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

      std::vector<Point> entriesCopy;
      tree_node_allocator *allocator = get_node_allocator(treeRef);

      // Create a copy of entries
      for (unsigned i = 0; i < cur_offset_; i++) {
        entriesCopy.push_back(entries.at(i));
      }

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

      SplitResult split = {
        {left_node->boundingBox(), left_handle},
        {right_node->boundingBox(), right_handle}
      };

      return split;
    }

    template <int min_branch_factor, int max_branch_factor>
    SplitResult propagateSplit(
      RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
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
      newPropagationSplit = current_node->splitNode(treeRef, current_handle, current_node->partitionNode());

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
      RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
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
      newPropagationSplit = current_node->splitNode(treeRef, current_handle, current_node->partitionNode());

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
      RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
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
            RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
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
        split = splitNode(treeRef, current_handle, partitionNode());

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
        RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef
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

      std::vector<Branch> entriesCopy;
      tree_node_allocator *allocator = get_node_allocator(treeRef);

      // Create a copy of entries
      for (unsigned i = 0; i < cur_offset_; i++) {
        entriesCopy.push_back(entries.at(i));
      }

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

      SplitResult split = {
              {left_node->boundingBox(), left_handle},
              {right_node->boundingBox(), right_handle}
      };

      return split;
    }

// insert() is always called on the root node.
    template <int min_branch_factor, int max_branch_factor>
    tree_node_handle BranchNode<min_branch_factor, max_branch_factor>::insert(
            RPlusTreeDisk<min_branch_factor, max_branch_factor> *treeRef,
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
