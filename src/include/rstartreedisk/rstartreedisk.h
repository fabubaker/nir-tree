#pragma once
#include <cassert>
#include <index/index.h>
#include <iostream>
#include <memory>
#include <rstartreedisk/node.h>
#include <stack>
#include <storage/tree_node_allocator.h>
#include <util/bmpPrinter.h>
#include <util/geometry.h>
#include <vector>
#include <iostream>
#include <fstream>

#include <fcntl.h>
#include <unistd.h>

namespace rstartreedisk {
template <int min_branch_factor, int max_branch_factor>
class RStarTreeDisk : public Index {
public:
  static constexpr float p = 0.3; // For reinsertion entries. 0.3 by default

  tree_node_handle root;
  std::unique_ptr<tree_node_allocator> node_allocator_;
  std::string backing_file_;

  std::vector<bool> hasReinsertedOnLevel;

  // Constructors and destructors
  RStarTreeDisk(size_t memory_budget, std::string backing_file) : backing_file_(backing_file) {
    node_allocator_ = std::make_unique<tree_node_allocator>(memory_budget, backing_file);

    // Initialize buffer pool
    node_allocator_->initialize();

    hasReinsertedOnLevel = {false};

    /* We need to figure out if there was already data, and read
                 * that into memory if we have it. */
    size_t existing_page_count = node_allocator_->buffer_pool_.get_preexisting_page_count();

    // If this is a fresh tree, then make a fresh root
    if (existing_page_count == 0) {
      auto node_type = NodeHandleType(LEAF_NODE);
      auto alloc = node_allocator_->create_new_tree_node<LeafNode<min_branch_factor, max_branch_factor>>(node_type);
      root = alloc.second;
      root.set_level(0);
      new (&(*(alloc.first))) LeafNode<min_branch_factor, max_branch_factor>();
      return;
    }

    std::string meta_file = backing_file_ + ".meta";
    int fd = open(meta_file.c_str(), O_RDONLY);
    assert(fd >= 0);

    int rc = read(fd, (char *)&root, sizeof(root));
    assert(rc == sizeof(root));
  }

  ~RStarTreeDisk(){};

  // Datastructure interface
  std::vector<Point> exhaustiveSearch(Point requestedPoint);
  std::vector<Point> search(Point requestedPoint);
  std::vector<Point> search(Rectangle requestedRectangle);
  void insert(Point givenPoint);
  void remove(Point givenPoint);

  // Miscellaneous
  unsigned checksum();
  void print();
  bool validate();
  void stat();
  void visualize();

  inline pinned_node_ptr<LeafNode<min_branch_factor, max_branch_factor>> get_leaf_node(tree_node_handle node_handle) {
    assert(node_handle.get_type() == LEAF_NODE);
    auto ptr = node_allocator_->get_tree_node<LeafNode<min_branch_factor, max_branch_factor>>(node_handle);
    return ptr;
  }

  inline pinned_node_ptr<BranchNode<min_branch_factor, max_branch_factor>> get_branch_node(tree_node_handle node_handle) {
    assert(node_handle.get_type() == BRANCH_NODE);
    auto ptr = node_allocator_->get_tree_node<BranchNode<min_branch_factor, max_branch_factor>>(node_handle);
    return ptr;
  }

  void write_metadata() {
    // Step 1:
    // Writeback everything to disk
    node_allocator_->buffer_pool_.writeback_all_pages();

    // Step 2:
    // Write metadata file
    std::string meta_fname = node_allocator_->get_backing_file_name() + ".meta";

    int fd = open(meta_fname.c_str(), O_WRONLY | O_TRUNC | O_CREAT, S_IRUSR | S_IWUSR);
    assert(fd >= 0);
    // yes, yes, i should loop this whatever
    int rc = write(fd, (char *)&root, sizeof(root));
    assert(rc == sizeof(root));
    close(fd);
  }
};

template <int min_branch_factor, int max_branch_factor>
std::vector<Point> RStarTreeDisk<min_branch_factor, max_branch_factor>::exhaustiveSearch(Point requestedPoint) {
  std::vector<Point> v;
  if (root.get_type() == LEAF_NODE) {
    auto root_ptr = get_leaf_node(root);
    root_ptr->exhaustiveSearch(requestedPoint, v);
  } else {
    auto root_ptr = get_branch_node(root);
    root_ptr->exhaustiveSearch(requestedPoint, v);
  }

  return v;
}

template <int min_branch_factor, int max_branch_factor>
std::vector<Point> RStarTreeDisk<min_branch_factor, max_branch_factor>::search(Point requestedPoint) {
  return point_search(root, requestedPoint, this);
}

template <int min_branch_factor, int max_branch_factor>
std::vector<Point> RStarTreeDisk<min_branch_factor, max_branch_factor>::search(Rectangle
                                                                               requestedRectangle) {
  return rectangle_search(root, requestedRectangle, this);
}

template <int min_branch_factor, int max_branch_factor>
void RStarTreeDisk<min_branch_factor, max_branch_factor>::insert(Point givenPoint) {
  if (root.get_type() == LEAF_NODE) {
    auto root_ptr = get_leaf_node(root);
    std::fill(hasReinsertedOnLevel.begin(), hasReinsertedOnLevel.end(), false);
    root = root_ptr->insert(this, this->root, givenPoint, hasReinsertedOnLevel);
    return;
  }
  auto root_ptr = get_branch_node(root);
  std::fill(hasReinsertedOnLevel.begin(), hasReinsertedOnLevel.end(), false);
  root = root_ptr->insert(this, this->root, givenPoint, hasReinsertedOnLevel);
}

template <int min_branch_factor, int max_branch_factor>
void RStarTreeDisk<min_branch_factor, max_branch_factor>::remove(Point givenPoint) {
  std::fill(hasReinsertedOnLevel.begin(), hasReinsertedOnLevel.end(), false);
  if (root.get_type() == LEAF_NODE) {
    auto root_ptr = get_leaf_node(root);
    root = root_ptr->remove(givenPoint, hasReinsertedOnLevel);
  } else {
    auto root_ptr = get_branch_node(root);
    root = root_ptr->remove(givenPoint, hasReinsertedOnLevel);
  }
}

template <int min_branch_factor, int max_branch_factor>
unsigned RStarTreeDisk<min_branch_factor, max_branch_factor>::checksum() {
  if (root.get_type() == LEAF_NODE) {
    auto root_ptr = get_leaf_node(root);
    return root_ptr->checksum();
  }
  auto root_ptr = get_branch_node(root);
  return root_ptr->checksum();
}

template <int min_branch_factor, int max_branch_factor>
void RStarTreeDisk<min_branch_factor, max_branch_factor>::print() {
  std::ofstream outputFile("printed_rstar_tree.txt");

  struct Printer {
      Printer(std::ofstream &printFile): printFile(printFile) {}

      void operator()(RStarTreeDisk<min_branch_factor, max_branch_factor> *treeRef, tree_node_handle node_handle) {
        printPackedNodes<min_branch_factor, max_branch_factor>(treeRef, node_handle, printFile);
      }

      std::ofstream &printFile;
  };

  Printer printer(outputFile);
  treeWalker<min_branch_factor, max_branch_factor>(this, root, printer);

  outputFile.close();
}

template <int min_branch_factor, int max_branch_factor>
bool RStarTreeDisk<min_branch_factor, max_branch_factor>::validate() {
  return true;
}

template <int min_branch_factor, int max_branch_factor>
void RStarTreeDisk<min_branch_factor, max_branch_factor>::stat() {
  stat_node(root, this);
}

template <int min_branch_factor, int max_branch_factor>
void RStarTreeDisk<min_branch_factor, max_branch_factor>::visualize() {

}
} // namespace rstartreedisk
