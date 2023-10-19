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

#include <cassert>
#include <index/index.h>
#include <iostream>
#include <memory>
#include <nirtreedisk/node.h>
#include <queue>
#include <stack>
#include <storage/tree_node_allocator.h>
#include <string>
#include <util/bmpPrinter.h>
#include <util/geometry.h>
#include <util/statistics.h>
#include <utility>
#include <vector>

#include <fcntl.h>
#include <unistd.h>

#include <map>
#include <fstream>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <boost/serialization/map.hpp>

namespace nirtreedisk {

template <int min_branch_factor, int max_branch_factor>
class NIRTreeDisk : public Index {
public:
  tree_node_handle root;
  std::unique_ptr<tree_node_allocator> node_allocator_;
  std::vector<bool> hasReinsertedOnLevel;
  std::map<tree_node_handle, IsotheticPolygon> polygons;
  BranchPartitionStrategy strategy;

  // Constructors and destructors
  NIRTreeDisk(size_t memory_budget, std::string backing_file, BranchPartitionStrategy partition_strategy) : 
    node_allocator_(std::make_unique<tree_node_allocator>(memory_budget, backing_file)),
    strategy(partition_strategy) {
    node_allocator_->initialize();

    size_t existing_page_count = node_allocator_->buffer_pool_.get_preexisting_page_count();

    hasReinsertedOnLevel = {false};    
    
    // If this is a fresh tree, we need a root
    if (existing_page_count == 0) {
      auto node_type = NodeHandleType(LEAF_NODE);
      auto alloc = node_allocator_->create_new_tree_node<LeafNode<min_branch_factor, max_branch_factor>>(node_type);
      root = alloc.second;
      new (&(*(alloc.first))) LeafNode<min_branch_factor, max_branch_factor>();
      return;
    }

    std::string meta_file = backing_file + ".meta";
    int fd = open(meta_file.c_str(), O_RDONLY);
    assert(fd >= 0);

    int rc = read(fd, (char *)&root, sizeof(root));
    assert(rc == sizeof(root));

    std::string polygon_fname = backing_file + ".polygons";
    std::ifstream polygon_file(polygon_fname);
    boost::archive::text_iarchive polygon_archive(polygon_file);
    polygon_archive >> polygons;
    polygon_file.close();
  }

  ~NIRTreeDisk() {
    if (root.get_type() == LEAF_NODE){
      auto root_node = get_leaf_node(root);
      root_node->deleteSubtrees();
      node_allocator_->free(root, sizeof(LeafNode<min_branch_factor, max_branch_factor>));
    } else if (root.get_type() == BRANCH_NODE){
      auto root_node = get_branch_node(root);
      root_node->deleteSubtrees(this);
      node_allocator_->free(root, sizeof(BranchNode<min_branch_factor, max_branch_factor>));
    }
  }

  // Datastructure interface
  std::vector<Point> exhaustiveSearch(Point requestedPoint);
  std::vector<Point> search(Point requestedPoint);
  std::vector<Point> search(Rectangle requestedRectangle);
  void insert(Point givenPoint);
  void remove(Point givenPoint);

  // Miscellaneous
  unsigned checksum();
  bool validate();
  void stat();
  void print();
  void visualize();

  inline Point get_a_contained_point(tree_node_handle node_handle) {
    tree_node_allocator *allocator = node_allocator_.get();

    while (node_handle.get_type() != LEAF_NODE) {
      assert(node_handle.get_type() == BRANCH_NODE);

      auto branch_node =
              allocator->get_tree_node<BranchNode<min_branch_factor, max_branch_factor>>(node_handle);
      Branch &b = branch_node->entries.at(0);
      node_handle = b.child;
    }

    assert(node_handle.get_type() == LEAF_NODE);

    auto leaf_node =
        allocator->get_tree_node<LeafNode<min_branch_factor, max_branch_factor>>(node_handle);
    return leaf_node->entries.at(0);
  }

  /* If unpack_perm=true, the tree will be modified */
  inline pinned_node_ptr<LeafNode<min_branch_factor, max_branch_factor>>
  get_leaf_node(tree_node_handle node_handle, bool unpack_perm = true) {
    assert(node_handle.get_type() == LEAF_NODE);
    auto ptr =
        node_allocator_->get_tree_node<LeafNode<min_branch_factor, max_branch_factor>>(node_handle);
    return ptr;
  }

  /* If unpack_perm=true, the tree will be modified */
  inline pinned_node_ptr<BranchNode<min_branch_factor, max_branch_factor>>
  get_branch_node(tree_node_handle node_handle, bool unpack_perm = true) {
    assert(node_handle.get_type() == BRANCH_NODE);
    auto ptr =
        node_allocator_->get_tree_node<BranchNode<min_branch_factor, max_branch_factor>>(node_handle);
    return ptr;
  }

  void write_metadata() override {
    // Step 1:
    // Writeback everything to disk
    node_allocator_->buffer_pool_.writeback_all_pages();

    // Step 2:
    // Write metadata file
    std::string meta_fname = node_allocator_->get_backing_file_name() + ".meta";
    std::string polygon_fname = node_allocator_->get_backing_file_name() + ".polygons";


    int fd = open(meta_fname.c_str(), O_WRONLY | O_TRUNC | O_CREAT, S_IRUSR | S_IWUSR);
    assert(fd >= 0);
    // yes, yes, i should loop this whatever
    int rc = write(fd, (char *)&root, sizeof(root));
    assert(rc == sizeof(root));
    close(fd);

    std::ofstream polygon_file(polygon_fname);
    boost::archive::text_oarchive polygon_archive(polygon_file);
    polygon_archive << polygons;
    polygon_file.close();
  }
};

template <int min_branch_factor, int max_branch_factor>
std::vector<Point> NIRTreeDisk<min_branch_factor, max_branch_factor>::exhaustiveSearch(Point requestedPoint) {
  return point_search(root, requestedPoint, this);
}

template <int min_branch_factor, int max_branch_factor>
std::vector<Point>
NIRTreeDisk<min_branch_factor, max_branch_factor>::search(Point requestedPoint) {
  return point_search(root, requestedPoint, this);
}

template <int min_branch_factor, int max_branch_factor>
std::vector<Point>
NIRTreeDisk<min_branch_factor, max_branch_factor>::search(Rectangle requestedRectangle) {
  return rectangle_search(root, requestedRectangle, this, true /* track */);
}

template <int min_branch_factor, int max_branch_factor>
void NIRTreeDisk<min_branch_factor, max_branch_factor>::insert(Point givenPoint) {
  std::fill(hasReinsertedOnLevel.begin(), hasReinsertedOnLevel.end(), false);
  if (root.get_type() == LEAF_NODE) {
    auto root_node = get_leaf_node(root, true);
    root = root_node->insert(this, root, givenPoint, hasReinsertedOnLevel);
  } else {
    auto root_node = get_branch_node(root, true);
    std::variant<BranchAtLevel, Point> entry = givenPoint;
    root = root_node->insert(this, root, entry, hasReinsertedOnLevel);
  }
}

template <int min_branch_factor, int max_branch_factor>
void NIRTreeDisk<min_branch_factor, max_branch_factor>::remove(Point givenPoint) {
  if (root.get_type() == LEAF_NODE) {
    auto root_node = get_leaf_node(root);
    auto result = root_node->remove(this, root, givenPoint);
    if (result == nullptr) {
      std::cout << "Point " << givenPoint << " is not found in tree" << std::endl;
      // we expect remove() is only called on points which exist
      assert(result != nullptr);
    } else {
      root = result;
    }
  } else {
    auto root_node = get_branch_node(root);
    auto result = root_node->remove(this, root, givenPoint);
    if (result == nullptr) {
      std::cout << "Point " << givenPoint << " is not found in tree" << std::endl;
      assert(result != nullptr);
    } else {
      root = result;
    }
  }
}

template <int min_branch_factor, int max_branch_factor>
unsigned NIRTreeDisk<min_branch_factor, max_branch_factor>::checksum() {
  if (root.get_type() == LEAF_NODE) {
    auto root_node = get_leaf_node(root);
    return root_node->checksum();
  } else {
    auto root_node = get_branch_node(root);
    return root_node->checksum(this);
  }
}

template <int min_branch_factor, int max_branch_factor>
bool NIRTreeDisk<min_branch_factor, max_branch_factor>::validate() {
  if (root.get_type() == LEAF_NODE) {
    auto root_node = get_leaf_node(root);
    root_node->bounding_box_validate();
    return root_node->validate(tree_node_handle(nullptr), 0);
  } else {
    auto root_node = get_branch_node(root);
    root_node->bounding_box_validate();
    return root_node->validate(tree_node_handle(nullptr), 0);
  }
}

template <int min_branch_factor, int max_branch_factor>
void NIRTreeDisk<min_branch_factor, max_branch_factor>::stat() {
  stat_node(root, this);
}

template <int min_branch_factor, int max_branch_factor>
void NIRTreeDisk<min_branch_factor, max_branch_factor>::print() {
  std::ofstream outputFile("printed_nir_tree.txt");

  struct Printer {
      Printer(std::ofstream &printFile): printFile(printFile) {}

      void operator()(NIRTreeDisk<min_branch_factor, max_branch_factor> *treeRef, tree_node_handle node_handle) {
        if (node_handle.get_type() == BRANCH_NODE) {
          auto node = treeRef->get_branch_node(node_handle);
          printFile << node->boundingBox() << std::endl;
        } else {
          auto node = treeRef->get_leaf_node(node_handle);
          printFile << node->boundingBox() << std::endl;

          // Also print points
          for (size_t i = 0; i < node->cur_offset_; i++) {
            printFile << node->entries[i] << std::endl;
          }
        }
      }

      std::ofstream &printFile;
  };

  Printer printer(outputFile);
  treeWalker<min_branch_factor, max_branch_factor>(this, root, printer);

  outputFile.close();
}

template <int min_branch_factor, int max_branch_factor>
void NIRTreeDisk<min_branch_factor, max_branch_factor>::visualize() {
  //BMPPrinter p(1000, 1000);
  // p.printToBMP(root);
}
} // namespace nirtreedisk
