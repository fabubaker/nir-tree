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

/* STR */

template <typename T>
std::vector<tree_node_handle> str_packing_branch(
        T *tree,
        std::vector<tree_node_handle> &child_nodes,
        unsigned branch_factor,
        unsigned cur_depth);

template <typename T>
std::vector<tree_node_handle> str_packing_leaf(
    T *tree,
    std::vector<Point>::iterator begin,
    std::vector<Point>::iterator end,
    unsigned branch_factor,
    unsigned cur_depth);

/* QTS */

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

/* Bulk-load/insert entry points */

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
