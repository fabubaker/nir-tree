#pragma once

#include <cassert>
#include <vector>
#include <stack>
#include <queue>
#include <string>
#include <iostream>
#include <utility>
#include <util/geometry.h>
#include <rplustreedisk/node.h>
#include <storage/tree_node_allocator.h>
#include <index/index.h>
#include <util/statistics.h>
#include <unistd.h>
#include <fcntl.h>

namespace rplustreedisk
{
  template <int min_branch_factor, int max_branch_factor> class RPlusTreeDisk: public Index
	{
		public:
			tree_node_handle root_;
      tree_node_allocator node_allocator_;
      std::string backing_file_;

      // Not required, but we mimic R* so file formats are the same.
      std::vector<bool> hasReinsertedOnLevel;

			// Constructors and destructors
			RPlusTreeDisk(size_t memory_budget, const std::string &backing_file):
        node_allocator_(memory_budget, backing_file), backing_file_(backing_file)
      {
        node_allocator_.initialize();

        size_t existing_page_count = node_allocator_.buffer_pool_.get_preexisting_page_count();

        // If existing page count is zero, then create a new
        // root node.
        if (existing_page_count == 0) {
            auto alloc_data = node_allocator_.create_new_tree_node<Node<min_branch_factor,max_branch_factor>>();
            root_ = alloc_data.second;

            new (&(*(alloc_data.first))) Node<min_branch_factor,max_branch_factor>(
              this, root_, tree_node_handle(nullptr)
            );

            return;
        }

        // Find existing root node.
        std::string meta_file = backing_file_ + ".meta";
        int fd = open( meta_file.c_str(), O_RDONLY );
        assert( fd >= 0 );

        int rc = read( fd, (char *) &root_, sizeof( root_ ) );
        assert( rc == sizeof( root_ ) );
      }

			~RPlusTreeDisk();

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

      inline pinned_node_ptr<Node<min_branch_factor,max_branch_factor>> get_node(tree_node_handle node_handle) {
          auto ptr =
              node_allocator_.get_tree_node<Node<min_branch_factor,max_branch_factor>>(node_handle);

          ptr->treeRef = this;
          return ptr;
      }

      void write_metadata() override {
        // Step 1:
        // Writeback everything to disk
        node_allocator_.buffer_pool_.writeback_all_pages();

        // Step 2:
        // Write metadata file

        auto root_node = get_node( root_ );
        assert( root_node->self_handle_ == root_ );
        std::string meta_fname = backing_file_ + ".meta";
        int fd = open( meta_fname.c_str(), O_WRONLY |
                O_TRUNC | O_CREAT, S_IRUSR | S_IWUSR );
        assert( fd >= 0 );
        // yes, yes, i should loop this whatever
        int rc = write( fd, (char *) &root_, sizeof(root_) );
        assert( rc == sizeof(root_) );
        close( fd );
      }
	};

  template <int min_branch_factor, int max_branch_factor>
  RPlusTreeDisk<min_branch_factor, max_branch_factor>::~RPlusTreeDisk() {
    auto root_node = get_node( root_ );
    root_node->deleteSubtrees();

    // FIXME GC root;
  }

  template <int min_branch_factor, int max_branch_factor>
  std::vector<Point> RPlusTreeDisk<min_branch_factor, max_branch_factor>::exhaustiveSearch(
          Point requestedPoint
  ) {
    std::vector<Point> v;
    auto root_node = get_node( root_ );
    root_node->exhaustiveSearch( requestedPoint, v );

    return v;
  }

  template <int min_branch_factor, int max_branch_factor>
  std::vector<Point> RPlusTreeDisk<min_branch_factor, max_branch_factor>::search(
          Point requestedPoint
  ) {
    auto root_node = get_node( root_ );
    return root_node->search( requestedPoint );
  }

  template <int min_branch_factor, int max_branch_factor>
  std::vector<Point> RPlusTreeDisk<min_branch_factor, max_branch_factor>::search(
          Rectangle requestedRectangle
  ) {
    auto root_node = get_node( root_ );
    return root_node->search( requestedRectangle );
  }

  template <int min_branch_factor, int max_branch_factor>
  void RPlusTreeDisk<min_branch_factor, max_branch_factor>::insert(
          Point givenPoint
  ) {
    auto root_node = get_node( root_ );
    root_ = root_node->insert( givenPoint );
  }

  template <int min_branch_factor, int max_branch_factor>
  void RPlusTreeDisk<min_branch_factor, max_branch_factor>::remove(
          Point givenPoint
  ) {
    auto root_node = get_node( root_ );
    root_ = root_node->remove( givenPoint );
  }

  template <int min_branch_factor, int max_branch_factor>
  unsigned RPlusTreeDisk<min_branch_factor, max_branch_factor>::checksum()
  {
    auto root_node = get_node( root_ );
    return root_node->checksum();
  }

  template <int min_branch_factor, int max_branch_factor>
  bool RPlusTreeDisk<min_branch_factor, max_branch_factor>::validate()
  {
    return true;
  }

  template <int min_branch_factor, int max_branch_factor>
  void RPlusTreeDisk<min_branch_factor, max_branch_factor>::stat()
  {
    auto root_node = get_node( root_ );
    root_node->stat();
  }

  template <int min_branch_factor, int max_branch_factor>
  void RPlusTreeDisk<min_branch_factor, max_branch_factor>::print()
  {
    auto root_node = get_node( root_ );
    root_node->printTree();
  }

  template <int min_branch_factor, int max_branch_factor>
  void RPlusTreeDisk<min_branch_factor, max_branch_factor>::visualize()
  {
    //BMPPrinter p( 1000, 1000 );
    //p.printToBMP( root_node );
  }
}
