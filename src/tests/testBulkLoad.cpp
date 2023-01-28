#include <bench/randomPoints.h>
#include <nirtreedisk/node.h>
#include <nirtreedisk/nirtreedisk.h>
#include <vector>
#include <cmath>
#include <unistd.h>
#include <random>
#include <catch2/catch.hpp>
#include <bulk_load.h>
#include <string>
#include <stdio.h>


std::vector<Point> generate_points_uniform( unsigned size ) {
    std::vector<Point> all_points;
    all_points.reserve( size );
    std::optional<Point> next;

    BenchTypeClasses::Uniform::size = size;
    BenchTypeClasses::Uniform::dimensions = dimensions;
    BenchTypeClasses::Uniform::seed = 0;
    PointGenerator<BenchTypeClasses::Uniform> points;
    while( (next = points.nextPoint() )) {
        all_points.push_back( next.value() );
    }

    return all_points;
}

std::vector<Point> generate_points_zipf( unsigned size ) {
    std::vector<Point> all_points;
    all_points.reserve( size );
    std::optional<Point> next;

    BenchTypeClasses::Zipf::size = size;
    BenchTypeClasses::Zipf::dimensions = dimensions;
    BenchTypeClasses::Zipf::seed = 0;
    BenchTypeClasses::Zipf::alpha = 1.0;
    BenchTypeClasses::Zipf::num_elements = 1000;
    PointGenerator<BenchTypeClasses::Zipf> points;
    while( (next = points.nextPoint() )) {
        all_points.push_back( next.value() );
    }

    return all_points;
}

std::vector<Point> generate_points_gauss( unsigned size ) {
    std::vector<Point> all_points;
    all_points.reserve( size );
    std::optional<Point> next;

    BenchTypeClasses::Gauss::size = size;
    BenchTypeClasses::Gauss::dimensions = dimensions;
    BenchTypeClasses::Gauss::seed = 0;
    PointGenerator<BenchTypeClasses::Gauss> points;
    while( (next = points.nextPoint() )) {
        all_points.push_back( next.value() );
    }

    return all_points;
}

void generate_tree(nirtreedisk::NIRTreeDisk<5,9,nirtreedisk::ExperimentalStrategy> *tree, std::vector<Point> &all_points, unsigned branch_factor) {
    double bulk_load_pct = 1.0;
    unsigned size = (unsigned) all_points.size();
    uint64_t cut_off_bulk_load = std::floor(bulk_load_pct*all_points.size());

    uint64_t max_depth =
        std::floor(log(size - 1)/log(branch_factor));
    auto ret = quad_tree_style_load( tree,
            all_points.begin(), all_points.begin() + cut_off_bulk_load,
            (unsigned) branch_factor, (unsigned) 0, (unsigned)
            max_depth, nullptr );
    tree->root = ret.first;
}

int count_underfill(tree_node_handle branch, nirtreedisk::NIRTreeDisk<5,9,nirtreedisk::ExperimentalStrategy> *tree) {
    if (branch.get_type() == LEAF_NODE || branch.get_type() == REPACKED_LEAF_NODE) {
        auto leaf_node = tree->get_leaf_node(branch);
        if (leaf_node->cur_offset_ != 9) {
            // printf("Leaf size: %d \n", leaf_node->cur_offset_);
            return 1;
        }
        return 0;
    } else if (branch.get_type() == BRANCH_NODE || branch.get_type() == REPACKED_BRANCH_NODE) {
        int count = 0;
        auto branch_node = tree->get_branch_node(branch);

        for( unsigned i = 0; i < branch_node->entries.size(); i++ ) {
            nirtreedisk::Branch &b = branch_node->entries.at(i);
            count += count_underfill(b.child, tree);
        }

        return count;
    }
    return 0;
}

void check_packing(nirtreedisk::NIRTreeDisk<5,9,nirtreedisk::ExperimentalStrategy> *tree) {
    tree_node_handle root = tree->root;
    int count = count_underfill(root, tree);
    assert(count <= 1);
}

TEST_CASE("gen_tree: test_uniform_generation") {
    unlink( "bulkloaded_tree_test.txt" );
    std::string file_name = "bulkloaded_tree_test.txt";
    nirtreedisk::NIRTreeDisk<5,9,nirtreedisk::ExperimentalStrategy> *tree =  new
            nirtreedisk::NIRTreeDisk<5,9,nirtreedisk::ExperimentalStrategy>(
                    40960UL*20000UL, file_name );
    unsigned size = 81;
    std::vector<Point> points = generate_points_uniform(size);
    generate_tree(tree, points, 9);

    for( Point p : points ) {
        std::vector<Point> out = tree->search(p);
        REQUIRE(out.size() == 1);
    }
    check_packing(tree);

    delete tree;
    unlink( "bulkloaded_tree_test.txt" );
}

TEST_CASE("gen_tree: test_zipf_generation") {
    unlink( "bulkloaded_tree_test.txt" );
    std::string file_name = "bulkloaded_tree_test.txt";
    nirtreedisk::NIRTreeDisk<5,9,nirtreedisk::ExperimentalStrategy> *tree =  new
            nirtreedisk::NIRTreeDisk<5,9,nirtreedisk::ExperimentalStrategy>(
                    40960UL*20000UL, file_name );
    unsigned size = 81;
    std::vector<Point> points = generate_points_zipf(size);
    generate_tree(tree, points, 9);

    for( Point p : points ) {
        std::vector<Point> out = tree->search(p);
        REQUIRE(out.size() == 1);
    }
    check_packing(tree);

    delete tree;
    unlink( "bulkloaded_tree_test.txt" );
}

TEST_CASE("gen_tree: test_gauss_generation") {
    unlink( "bulkloaded_tree_test.txt" );
    std::string file_name = "bulkloaded_tree_test.txt";
    nirtreedisk::NIRTreeDisk<5,9,nirtreedisk::ExperimentalStrategy> *tree =  new
            nirtreedisk::NIRTreeDisk<5,9,nirtreedisk::ExperimentalStrategy>(
                    40960UL*20000UL, file_name );
    unsigned size = 81;
    std::vector<Point> points = generate_points_gauss(size);
    generate_tree(tree, points, 9);

    for( Point p : points ) {
        std::vector<Point> out = tree->search(p);
        REQUIRE(out.size() == 1);
    }
    check_packing(tree);

    delete tree;
    unlink( "bulkloaded_tree_test.txt" );
}

TEST_CASE("gen_tree: test_non_perfect_depth") {
    unlink( "test_non_perfect_depth.txt" );
    std::string file_name = "test_non_perfect_depth.txt";
    nirtreedisk::NIRTreeDisk<5,9,nirtreedisk::ExperimentalStrategy> *tree =  new
            nirtreedisk::NIRTreeDisk<5,9,nirtreedisk::ExperimentalStrategy>(
                    4096UL*20, file_name );
    unsigned size = 1203;
    std::vector<Point> points = generate_points_uniform(size);
    generate_tree(tree, points, 9);

    for( Point p : points ) {
        std::vector<Point> out = tree->search(p);
        REQUIRE(out.size() == 1);
    }
    check_packing(tree);

    delete tree;
    unlink( "test_non_perfect_depth.txt" );

}

TEST_CASE("gen_tree: test_same_y") {
    unlink( "test_same_y.txt" );
    std::string file_name = "test_same_y.txt";
    nirtreedisk::NIRTreeDisk<5,9,nirtreedisk::ExperimentalStrategy> *tree =  new
            nirtreedisk::NIRTreeDisk<5,9,nirtreedisk::ExperimentalStrategy>(
                    4096UL*200, file_name);
    unsigned size = 81;
    std::vector<Point> points;
    for(unsigned i = 0; i < size; i++) {
        Point p;
        p[0] = (double) i;
        p[1] = 0;
        points.push_back(std::move(p));
    }
    generate_tree(tree, points, 9);

    for( Point p : points ) {
        std::vector<Point> out = tree->search(p);
        REQUIRE(out.size() == 1);
    }
    check_packing(tree);

    delete tree;
    unlink( "test_same_y.txt" );
}

TEST_CASE("gen_tree: test_same_x") {
    unlink( "test_same_x.txt" );
    std::string file_name = "test_same_x.txt";
    nirtreedisk::NIRTreeDisk<5,9,nirtreedisk::ExperimentalStrategy> *tree =  new
            nirtreedisk::NIRTreeDisk<5,9,nirtreedisk::ExperimentalStrategy>(
                    4096UL*200, file_name);
    unsigned size = 81;
    std::vector<Point> points;
    for(unsigned i = 0; i < size; i++) {
        Point p;
        p[1] = (double) i;
        p[0] = 0;
        points.push_back(std::move(p));
    }
    generate_tree(tree, points, 9);

    for( Point p : points ) {
        std::vector<Point> out = tree->search(p);
        REQUIRE(out.size() == 1);
    }
    check_packing(tree);

    delete tree;
    unlink( "test_same_x.txt" );
}

TEST_CASE("gen_tree: underfill_tree") {
    unlink( "underfill_tree.txt" );
    std::string file_name = "underfill_tree.txt";
    nirtreedisk::NIRTreeDisk<5,9,nirtreedisk::ExperimentalStrategy> *tree =  new
            nirtreedisk::NIRTreeDisk<5,9,nirtreedisk::ExperimentalStrategy>(
                    4096UL*200UL, file_name);
    unsigned size = 8;
    std::vector<Point> points = generate_points_uniform(size);
    generate_tree(tree, points, 9);

    for( Point p : points ) {
        std::vector<Point> out = tree->search(p);
        REQUIRE(out.size() == 1);
    }
    check_packing(tree);

    delete tree;
    unlink( "underfill_tree.txt" );
}



TEST_CASE("gen_tree: optimally_pack") {
    unlink( "optimally_pack.txt" );
    std::string file_name = "optimally_pack.txt";

    nirtreedisk::NIRTreeDisk<5,9,nirtreedisk::ExperimentalStrategy> *tree =  new
            nirtreedisk::NIRTreeDisk<5,9,nirtreedisk::ExperimentalStrategy>(
                    4096UL*20, file_name );
    unsigned size = 12039;
    std::vector<Point> points = generate_points_uniform(size);
    generate_tree(tree, points, 9);

    for( Point p : points ) {
        std::vector<Point> out = tree->search(p);
        REQUIRE(out.size() == 1);
    }
    check_packing(tree);

    unlink( "optimally_pack.txt" );
}


