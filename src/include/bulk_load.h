#pragma once

#include <bench/randomPoints.h>
#include <cmath>
#include <iostream>
#include <nirtreedisk/nirtreedisk.h>
#include <nirtreedisk/node.h>
#include <random>
#include <unistd.h>
#include <vector>

enum BulkLoadAlg { STR, QTS, TGS };

std::pair<uint64_t, double> compute_max_dist(const Point &point, std::vector<Point> &pts);

bool point_comparator(const Point &lhs, const Point &rhs);

template <class TreeType>
void make_all_rects_disjoint(
    TreeType *treeRef,
    std::vector<Rectangle> &rects_a,
    tree_node_handle a_node,
    std::vector<Rectangle> &rects_b,
    tree_node_handle b_node);

/* STR */

template <typename T, typename LN, typename BN>
void fill_branch(
    T *treeRef,
    pinned_node_ptr<BN> branch_node,
    tree_node_handle node_handle,
    std::vector<std::pair<Point, tree_node_handle>> &node_point_pairs,
    uint64_t &offset,
    unsigned branch_factor,
    LN *leaf_type);

template <typename T>
std::vector<tree_node_handle> str_packing_branch(
        T *tree,
        std::vector<tree_node_handle> &child_nodes,
        unsigned branch_factor,
        unsigned cur_depth);

template <typename T, typename LN, typename BN>
std::vector<tree_node_handle> str_packing_branch(
        T *tree,
        std::vector<tree_node_handle> &child_nodes,
        unsigned branch_factor,
        LN *leaf_node_type,
        BN *branch_node_type,
        unsigned cur_depth);

template <typename T>
std::vector<tree_node_handle> str_packing_leaf(
    T *tree,
    std::vector<Point>::iterator begin,
    std::vector<Point>::iterator end,
    unsigned branch_factor,
    unsigned cur_depth);

template <typename T, typename LN, typename BN>
std::vector<tree_node_handle> str_packing_leaf(
        T *tree,
        std::vector<Point>::iterator begin,
        std::vector<Point>::iterator end,
        unsigned branch_factor,
        LN *ln_type,
        BN *bn_type,
        unsigned cur_depth);

/* QTS */

std::vector<uint64_t> find_bounding_lines(
    std::vector<Point>::iterator start,
    std::vector<Point>::iterator stop,
    unsigned d,
    unsigned partitions);

std::pair<tree_node_handle, Rectangle> quad_tree_style_load(
    nirtreedisk::NIRTreeDisk<5, NIR_FANOUT> *tree,
    std::vector<Point>::iterator start,
    std::vector<Point>::iterator stop,
    unsigned branch_factor,
    unsigned cur_level);

std::pair<tree_node_handle, Rectangle> quad_tree_style_load(
    rstartreedisk::RStarTreeDisk<R_STAR_MIN_FANOUT, R_STAR_MAX_FANOUT> *tree,
    std::vector<Point>::iterator start,
    std::vector<Point>::iterator stop,
    unsigned branch_factor,
    unsigned cur_level);

/* TGS */

tree_node_handle tgs_load(
        rstartreedisk::RStarTreeDisk<R_STAR_MIN_FANOUT, R_STAR_MAX_FANOUT> *tree,
        std::vector<Point>::iterator begin,
        std::vector<Point>::iterator end,
        unsigned branch_factor
);

template <typename T>
void bulk_load_tree(
    T *tree,
    std::map<std::string, size_t> &configU,
    std::vector<Point>::iterator begin,
    std::vector<Point>::iterator end,
    unsigned max_branch_factor);

template <typename T>
void sequential_insert_tree(
    T *tree,
    std::map<std::string, size_t> &configU,
    std::vector<Point>::iterator begin,
    std::vector<Point>::iterator end,
    unsigned max_branch_factor);
