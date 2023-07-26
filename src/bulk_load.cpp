#include <bulk_load.h>

unsigned intersection_count;

template <class TreeType>
void make_all_rects_disjoint(
    TreeType *treeRef,
    std::vector<Rectangle> &rects_a,
    tree_node_handle a_node,
    std::vector<Rectangle> &rects_b,
    tree_node_handle b_node) {
  std::vector<Rectangle> a_output;

  std::stack<Rectangle, std::vector<Rectangle>> remaining_a_rects(rects_a);

  while (not remaining_a_rects.empty()) {
    Rectangle a = remaining_a_rects.top();
    remaining_a_rects.pop();
    bool did_split = false;
    for (uint64_t i = 0; i < rects_b.size(); i++) {
      Rectangle &b = rects_b.at(i);
      // If there is no intersection with this rectangle, keep
      // going
      if (not a.intersectsRectangle(b)) {
        continue;
      }
      // If there is, we need to split it.
      auto ret = nirtreedisk::make_rectangles_disjoint_accounting_for_region_ownership(
          treeRef,
          a,
          a_node,
          b,
          b_node);

      IsotheticPolygon poly1;
      poly1.basicRectangles = ret.first;
      poly1.recomputeBoundingBox();

      IsotheticPolygon poly2;
      poly2.basicRectangles = ret.second;
      poly2.recomputeBoundingBox();

      if (not poly1.disjoint(poly2)) {
        std::cout << "A: " << a << std::endl;
        std::cout << "B: " << b << std::endl;
        std::cout << "Poly1: " << poly1 << std::endl;
        std::cout << "Poly2: " << poly2 << std::endl;

        std::cout << "Intersection." << std::endl;
        poly1.intersection(poly2);
        std::cout << poly1 << std::endl;

        assert(false);
        abort();
      }
      assert(poly1.disjoint(poly2));

      for (auto &ret_a_rect : ret.first) {
        remaining_a_rects.push(ret_a_rect);
      }
      rects_b.erase(rects_b.begin() + i);
      for (auto &ret_b_rect : ret.second) {
        rects_b.push_back(ret_b_rect);
      }
      // Need to loop around because we broke the iterator, and
      // both a and b have new sets of rectangles
      did_split = true;
      break;
    }
    if (not did_split) {
      a_output.push_back(a);
    }
  }
  rects_a = a_output;
}

template <>
void fill_branch(
    nirtreedisk::NIRTreeDisk<5, NIR_FANOUT, nirtreedisk::ExperimentalStrategy> *treeRef,
    pinned_node_ptr<nirtreedisk::BranchNode<5, NIR_FANOUT, nirtreedisk::ExperimentalStrategy>> branch_node,
    tree_node_handle node_handle,
    std::vector<std::pair<Point, tree_node_handle>> &node_point_pairs,
    uint64_t &offset,
    unsigned branch_factor,
    nirtreedisk::LeafNode<5, NIR_FANOUT, nirtreedisk::ExperimentalStrategy> *leaf_type) {
  using LN = nirtreedisk::LeafNode<5, NIR_FANOUT, nirtreedisk::ExperimentalStrategy>;
  using BN = nirtreedisk::BranchNode<5, NIR_FANOUT, nirtreedisk::ExperimentalStrategy>;

  std::vector<std::pair<IsotheticPolygon, tree_node_handle>> fixed_bb_and_handles;
  tree_node_allocator *allocator = treeRef->node_allocator_.get();

  // Add up to branch factor items to it
  for (uint64_t i = 0; i < branch_factor; i++) {
    tree_node_handle child_handle = node_point_pairs[offset++].second;
    Rectangle bbox;
    // Adjust parent
    if (child_handle.get_type() == LEAF_NODE) {
      auto node = allocator->get_tree_node<LN>(child_handle);
      bbox = node->boundingBox();
    } else {
      auto node = allocator->get_tree_node<BN>(child_handle);
      bbox = node->boundingBox();
    }

    fixed_bb_and_handles.push_back(std::make_pair(IsotheticPolygon(bbox), child_handle));

    if (offset == node_point_pairs.size()) {
      break;
    }
  }

  for (uint64_t i = 0; i < fixed_bb_and_handles.size(); i++) {
    for (uint64_t j = i + 1; j < fixed_bb_and_handles.size(); j++) {

      std::vector<Rectangle> &existing_rects_a = fixed_bb_and_handles.at(i).first.basicRectangles;
      std::vector<Rectangle> &existing_rects_b = fixed_bb_and_handles.at(j).first.basicRectangles;
      make_all_rects_disjoint(
          treeRef,
          existing_rects_a,
          fixed_bb_and_handles.at(i).second,
          existing_rects_b,
          fixed_bb_and_handles.at(j).second
      );
    }
  }

  // Now we have made all the BoundingRegions disjoint.
  // It is time to add our children
  for (uint64_t i = 0; i < fixed_bb_and_handles.size(); i++) {
    nirtreedisk::Branch b;
    IsotheticPolygon constructed_poly = fixed_bb_and_handles.at(i).first;

    b.child = fixed_bb_and_handles.at(i).second;
    b.boundingBox = constructed_poly.boundingBox;

    // If the MBR has not been split into a polygon, don't keep it in the map.
    if (constructed_poly.basicRectangles.size() != 1) {
      treeRef->polygons.insert({b.child, constructed_poly});
    }

    branch_node->addBranchToNode(b);
  }
}

template <>
void fill_branch(
    rstartreedisk::RStarTreeDisk<5, R_STAR_FANOUT> *treeRef,
    pinned_node_ptr<rstartreedisk::BranchNode<5, R_STAR_FANOUT>> branch_node,
    tree_node_handle node_handle,
    std::vector<std::pair<Point, tree_node_handle>> &node_point_pairs,
    uint64_t &offset,
    unsigned branch_factor,
    rstartreedisk::LeafNode<5, R_STAR_FANOUT> *leaf_type) {
  using LN = rstartreedisk::LeafNode<5, R_STAR_FANOUT>;
  using BN = rstartreedisk::BranchNode<5, R_STAR_FANOUT>;

  std::vector<std::pair<Rectangle, tree_node_handle>> bb_and_handles;
  tree_node_allocator *allocator = treeRef->node_allocator_.get();

  // Add up to branch factor items to it
  for (uint64_t i = 0; i < branch_factor; i++) {
    tree_node_handle child_handle = node_point_pairs[offset++].second;
    Rectangle bbox;

    // Adjust parent
    if (child_handle.get_type() == LEAF_NODE) {
      auto node = allocator->get_tree_node<LN>(child_handle);
      bbox = node->boundingBox();
    } else {
      auto node = allocator->get_tree_node<BN>(child_handle);
      bbox = node->boundingBox();
    }
    
    bb_and_handles.push_back(std::make_pair(bbox, child_handle));

    rstartreedisk::Branch b;
    b.child = child_handle;
    b.boundingBox = bbox;
    branch_node->addBranchToNode(b);

    if (offset == node_point_pairs.size()) {
      break;
    }
  }
}

template <typename T, typename LN, typename BN>
std::vector<tree_node_handle> str_packing_branch(
    T *tree,
    std::vector<tree_node_handle> &child_nodes,
    unsigned branch_factor,
    LN *leaf_node_type,
    BN *branch_node_type,
    unsigned cur_depth) {
  tree_node_allocator *allocator = tree->node_allocator_.get();

  // Get bbox once for everything so I'm not materializing it
  // constantly
  std::vector<std::pair<Point, tree_node_handle>> node_point_pairs;
  node_point_pairs.reserve(child_nodes.size());

  for (tree_node_handle &child_handle : child_nodes) {
    Rectangle bbox;

    if (child_handle.get_type() == LEAF_NODE) {
      auto child = allocator->get_tree_node<LN>(child_handle);
      bbox = child->boundingBox();
    } else {
      auto child = allocator->get_tree_node<BN>(child_handle);
      bbox = child->boundingBox();
    }

    node_point_pairs.push_back(std::make_pair(bbox.centrePoint(), child_handle));
  }

  uint64_t P = node_point_pairs.size() / branch_factor;

  if (node_point_pairs.size() % branch_factor != 0) P++;

  double S_dbl = std::ceil(sqrt(P));
  uint64_t S = (uint64_t)S_dbl;

  // Sort by X
  std::sort(node_point_pairs.begin(), node_point_pairs.end(), [](std::pair<Point, tree_node_handle> &l, std::pair<Point, tree_node_handle> &r) { return l.first[0] < r.first[0]; });

  // There are |S| vertical stripes.
  // Each has |S|*branch_factor points
  for (uint64_t i = 0; i < S; i++) {
    uint64_t start_offset = i * (S * branch_factor);
    uint64_t stop_offset = (i + 1) * (S * branch_factor);

    if (stop_offset > node_point_pairs.size()) {
      stop_offset = node_point_pairs.size();
    }

    auto start_point = node_point_pairs.begin() + start_offset;
    auto stop_point = node_point_pairs.begin() + stop_offset;

    std::sort(start_point, stop_point, [](std::pair<Point, tree_node_handle> &l, std::pair<Point, tree_node_handle> &r) { return l.first[1] < r.first[1]; });

    if (stop_offset == node_point_pairs.size()) break;
  }

  std::vector<tree_node_handle> branches;
  uint64_t offset = 0;

  while (offset < node_point_pairs.size()) {
    // Create the branch node
    auto alloc_data = allocator->create_new_tree_node<BN>(NodeHandleType(BRANCH_NODE));

    new (&(*alloc_data.first)) BN();

    auto branch_node = alloc_data.first;
    tree_node_handle branch_handle = alloc_data.second;
    branch_handle.set_level(cur_depth);

    fill_branch(
        tree,
        branch_node,
        branch_handle,
        node_point_pairs,
        offset,
        branch_factor,
        (LN *) nullptr
    );

    branches.push_back(branch_handle);
  }

  return branches;
}

std::pair<uint64_t, double> compute_max_dist(const Point &point, std::vector<Point> &pts) {
  uint64_t idx = 0;
  double max_dist = std::numeric_limits<double>::min();
  for (uint64_t i = 0; i < pts.size(); i++) {
    double dist = point.fast_distance(pts.at(i));
    if (dist > max_dist) {
      max_dist = dist;
      idx = i;
    }
  }
  return std::make_pair(idx, max_dist);
}

bool point_comparator(const Point &lhs, const Point &rhs) {
  for (unsigned d = 0; d < dimensions; d++) {
    if (lhs[d] != rhs[d]) {
      return lhs[d] < rhs[d];
    }
  }
  return false;
}

template <>
std::vector<tree_node_handle> str_packing_branch(
    nirtreedisk::NIRTreeDisk<5, NIR_FANOUT, nirtreedisk::ExperimentalStrategy> *tree,
    std::vector<tree_node_handle> &child_nodes,
    unsigned branch_factor,
    unsigned cur_depth) {
  nirtreedisk::LeafNode<5, NIR_FANOUT, nirtreedisk::ExperimentalStrategy> *targ = nullptr;
  nirtreedisk::BranchNode<5, NIR_FANOUT, nirtreedisk::ExperimentalStrategy> *targ2 = nullptr;

  return str_packing_branch(
    tree, child_nodes, branch_factor, targ,targ2, cur_depth
  );
}

template <>
std::vector<tree_node_handle> str_packing_branch(
    rstartreedisk::RStarTreeDisk<5, R_STAR_FANOUT> *tree,
    std::vector<tree_node_handle> &child_nodes,
    unsigned branch_factor,
    unsigned cur_depth
) {
  rstartreedisk::LeafNode<5, R_STAR_FANOUT> *targ = nullptr;
  rstartreedisk::BranchNode<5, R_STAR_FANOUT> *targ2 = nullptr;


  return str_packing_branch(
    tree, child_nodes, branch_factor, targ, targ2, cur_depth
  );
}

template <typename T, typename LN, typename BN>
std::vector<tree_node_handle> str_packing_leaf(
    T *tree,
    std::vector<Point>::iterator begin,
    std::vector<Point>::iterator end,
    unsigned branch_factor,
    LN *ln_type,
    BN *bn_type,
    unsigned cur_depth
) {
  uint64_t count = (end - begin);
  uint64_t P = count / branch_factor;
  if (count % branch_factor != 0) {
    P++;
  }

  double S_dbl = std::ceil(sqrt(P));
  uint64_t S = (uint64_t)S_dbl;

  // Sort on lower left
  std::sort(begin, end, [](Point &l, Point &r) { return l[0] < r[0]; });

  // There are |S| vertical stripes.
  // Each has |S|*branch_factor points
  for (uint64_t i = 0; i < S; i++) {
    uint64_t start_offset = i * (S * branch_factor);
    uint64_t stop_offset = (i + 1) * (S * branch_factor);
    auto start_point = begin + start_offset;
    auto stop_point = stop_offset >= count ? end : begin + stop_offset;

    std::sort(start_point, stop_point, [](Point &l, Point &r) { return l[1] < r[1]; });

    if (stop_point == end) {
      break;
    }
  }

  tree_node_allocator *allocator = tree->node_allocator_.get();
  std::vector<tree_node_handle> leaves;
  uint64_t offset = 0;

  while (offset < count) {
    auto alloc_data = allocator->create_new_tree_node<LN>(NodeHandleType(LEAF_NODE));

    new (&(*alloc_data.first)) LN();
    auto leaf_node = alloc_data.first;
    tree_node_handle leaf_handle = alloc_data.second;
    leaf_handle.set_level(cur_depth);

    for (uint64_t i = 0; i < branch_factor; i++) {
      leaf_node->addPoint(*(begin + offset));
      offset++;
      if (offset == count) {
        break;
      }
    }

    leaves.push_back(leaf_handle);
  }

  return leaves;
}

template <>
std::vector<tree_node_handle> str_packing_leaf(
    nirtreedisk::NIRTreeDisk<5, NIR_FANOUT, nirtreedisk::ExperimentalStrategy> *tree,
    std::vector<Point>::iterator begin,
    std::vector<Point>::iterator end,
    unsigned branch_factor,
    unsigned cur_depth) {
  nirtreedisk::LeafNode<5, NIR_FANOUT, nirtreedisk::ExperimentalStrategy> *targ = nullptr;
  nirtreedisk::BranchNode<5, NIR_FANOUT, nirtreedisk::ExperimentalStrategy>
      *targ2 = nullptr;

  return str_packing_leaf(
    tree, begin, end, branch_factor,targ, targ2, cur_depth
  );
}

template <>
std::vector<tree_node_handle> str_packing_leaf(
    rstartreedisk::RStarTreeDisk<5, R_STAR_FANOUT> *tree,
    std::vector<Point>::iterator begin,
    std::vector<Point>::iterator end,
    unsigned branch_factor,
    unsigned cur_depth) {
  rstartreedisk::LeafNode<5, R_STAR_FANOUT> *targ = nullptr;
  rstartreedisk::BranchNode<5, R_STAR_FANOUT> *targ2 = nullptr;

  return str_packing_leaf(
    tree, begin, end, branch_factor,targ, targ2, cur_depth
  );
}

std::vector<uint64_t> find_bounding_lines(
    std::vector<Point>::iterator start,
    std::vector<Point>::iterator stop,
    unsigned d, // dimension
    unsigned branch_factor,
    unsigned partitions,
    unsigned sub_partitions,
    unsigned length) {
  std::sort(start, stop, [d](const Point &p1, const Point &p2) {
                         if (d == 0 && p1[0] == p2[0]) {return p1[1] < p2[1]; } 
                         else if (d == 1 && p1[1] == p2[1]) {return p1[0] < p2[0]; }
                         return p1[d] < p2[d]; });

  std::vector<uint64_t> lines;
  lines.reserve(partitions + 2); // We include x = 0 and x = count
  lines.push_back(0);

  uint64_t count = stop - start;
  unsigned optimally_filled_branch = pow(branch_factor, length); // Total number of children (points) starting from this level and below
  unsigned min_branches = std::ceil((double) count / (double) optimally_filled_branch); // Minimum branches at this level
  unsigned branches_per_partition = min_branches / partitions; // Divide branches at this level to partitions
  unsigned remaining_branches = min_branches % partitions; // Any leftover branches?

  unsigned partition_index = 0;
  for (unsigned partition = 0; partition < partitions - 1; partition++) {
    uint64_t partition_size = optimally_filled_branch * branches_per_partition;

    if (remaining_branches != 0 && partition_size + optimally_filled_branch <= optimally_filled_branch * sub_partitions) {
      partition_size += optimally_filled_branch;
      remaining_branches -= 1;
    }

    partition_index += partition_size;
    if (partition_index > count) {
      partition_index = count;
    }
    lines.push_back(partition_index);
  }
  lines.push_back(count);

  // If there are remaining branches, there are not enough partitions to host
  // all branches.
  assert(remaining_branches == 0);

  return lines;
}

std::pair<tree_node_handle, Rectangle> quad_tree_style_load(
  nirtreedisk::NIRTreeDisk<5, NIR_FANOUT, nirtreedisk::ExperimentalStrategy> *tree,
  std::vector<Point>::iterator start,
  std::vector<Point>::iterator stop,
  unsigned branch_factor,
  unsigned cur_depth,
  unsigned max_depth,
  tree_node_handle parent_handle
) {
  uint64_t num_els = (stop - start);
  tree_node_allocator *allocator = tree->node_allocator_.get();
  if (cur_depth == max_depth) {
    if (num_els > branch_factor) {
      std::cout << "NUM ELS: " << num_els << std::endl;
    }

    assert(num_els <= branch_factor);

    num_els = (stop - start);
    auto alloc_data =
      allocator->create_new_tree_node<nirtreedisk::LeafNode<5, NIR_FANOUT, nirtreedisk::ExperimentalStrategy>>(
                    NodeHandleType(LEAF_NODE));
    new (&(*(alloc_data.first))) nirtreedisk::LeafNode<5, NIR_FANOUT, nirtreedisk::ExperimentalStrategy>();

    auto leaf_node = alloc_data.first;
    auto leaf_handle = alloc_data.second;
    leaf_handle.set_level(cur_depth);

    for (auto iter = start; iter != stop; iter++) {
      leaf_node->addPoint(*iter);
    }

    Rectangle bbox = leaf_node->boundingBox();
    return std::make_pair(leaf_handle, bbox);
  }

  // Return a tree node handle with pointers to all of its necessary
  // children.
  auto alloc_data =
          allocator->create_new_tree_node<nirtreedisk::BranchNode<5, NIR_FANOUT, nirtreedisk::ExperimentalStrategy>>(
                  NodeHandleType(BRANCH_NODE));
  new (&(*(alloc_data.first))) nirtreedisk::BranchNode<5, NIR_FANOUT, nirtreedisk::ExperimentalStrategy>();

  auto branch_node = alloc_data.first;
  tree_node_handle branch_handle = alloc_data.second;
  branch_handle.set_level(cur_depth);

  uint64_t partitions = std::ceil(sqrt(branch_factor));
  uint64_t sub_partitions = std::ceil(branch_factor / (float) partitions);

  std::vector<uint64_t> x_lines = find_bounding_lines(start, stop, 0, branch_factor, partitions, sub_partitions, max_depth - cur_depth);

  std::vector<std::pair<IsotheticPolygon, tree_node_handle>> branch_handles;
  branch_handles.reserve(NIR_FANOUT);

  for (uint64_t i = 0; i < x_lines.size() - 1; i++) {
    uint64_t x_start = x_lines.at(i);
    uint64_t x_end = x_lines.at(i + 1); /* not inclusive */

    std::vector<uint64_t> y_lines = find_bounding_lines(
        start + x_start, start + x_end, 1, branch_factor, sub_partitions, 1, max_depth - cur_depth);
    for (uint64_t j = 0; j < y_lines.size() - 1; j++) {
      uint64_t y_start = y_lines.at(j);
      uint64_t y_end = y_lines.at(j + 1); /* not inclusive */

      std::vector<Point>::iterator sub_start = start + x_start + y_start;
      std::vector<Point>::iterator sub_stop = start + x_start + y_end;

      if (sub_start == sub_stop) {
        // I think this can happen when we run out of points in
        // the lowest layer to split among children.
        continue;
      }

      // This needs to return the box and handle.
      auto ret = quad_tree_style_load(
          tree,
          sub_start,
          sub_stop,
          branch_factor,
          cur_depth + 1,
          max_depth,
          branch_handle
      );

      tree_node_handle child_handle = ret.first;
      Rectangle bbox = ret.second;

      // Produce a polygon using bbox.
      IsotheticPolygon ip = IsotheticPolygon(bbox);

      // Look at each of the existing children's polygons, and figure out if
      // any of their polygons overlap with our polygon. If so, we need to
      // figure out who owns the overlapping region.
      for (uint64_t i = 0; i < branch_handles.size(); i++) {
        std::vector<Rectangle> &existing_rects = branch_handles.at(i).first.basicRectangles;

        if (branch_handles.at(i).first.intersectsPolygon(ip)) {
          // It is imperative that the references here are set carefully.
          // make_all_rects_disjoint updates the second and fourth arguments,
          // so these should point to the basicRectangles vector of the polygons.
          intersection_count += 1;
          make_all_rects_disjoint(
              tree,
              existing_rects,              // Existing rects of sibling poly
              branch_handles.at(i).second, // Sibling handle
              ip.basicRectangles,          // Existing rects of our poly
              child_handle                 // our handle
          );

          // These may have been updated. Update the metadata.
          ip.recomputeBoundingBox();
          branch_handles.at(i).first.recomputeBoundingBox();
        }
      }

      branch_handles.push_back(std::make_pair(ip, child_handle));

#ifndef NDEBUG
      // Double check non-intersection - Really inefficient, but I don't see a better way
      // of doing this.
      for (uint64_t i = 0; i < branch_handles.size(); i++) {
        for (uint64_t j = i + 1; j < branch_handles.size(); j++) {
          std::vector<Rectangle> &existing_rects_a =
              branch_handles.at(i).first.basicRectangles;
          std::vector<Rectangle> &existing_rects_b =
              branch_handles.at(j).first.basicRectangles;
          for (Rectangle &rect_a : existing_rects_a) {
            for (Rectangle &rect_b : existing_rects_b) {
              if (rect_b.intersectsRectangle(rect_a)) {
                std::cout << rect_b << " intersects: " << rect_a << std::endl;
                std::cout << "Cur y_lo: " << (*(start + x_start + y_start))[1] << std::endl;
                std::cout << "Prev y: " << (*(start + x_start + y_start - 1))[1] << std::endl;
                abort();
              }
            }
          }
        }
      }
#endif
    }
  }

  for (uint64_t i = 0; i < branch_handles.size(); i++) {
    nirtreedisk::Branch b;
    IsotheticPolygon constructed_poly = branch_handles.at(i).first;

    b.child = branch_handles.at(i).second;
    b.boundingBox = constructed_poly.boundingBox;

    // If the MBR has not been split into a polygon, don't keep it in the map.
    if (constructed_poly.basicRectangles.size() != 1) {
      tree->polygons.insert({b.child, constructed_poly});
    }

    branch_node->addBranchToNode(b);
  }

  Rectangle bbox = branch_node->boundingBox();
  return std::make_pair(branch_handle, bbox);
}


// loop from point from begin to end once 
// and create a MBB for points in this range
Rectangle min_bounding_box(
    std::vector<Point>::iterator begin, 
    std::vector<Point>::iterator end )
{
  if (end == begin) {
    // box contains only 1 point 
    assert(end != begin );
  }
  // assumes 2 dimensions 
  double smallest_x = (*begin)[0]; 
  double biggest_x = (*begin)[0]; 
  double smallest_y = (*begin)[1]; 
  double biggest_y = (*begin)[1]; 
  for (auto iter = begin; iter != end; iter++) {
    if ((*iter)[0] < smallest_x){ smallest_x = (*iter)[0]; }
    if ((*iter)[0] > biggest_x){ biggest_x = (*iter)[0]; }
    if ((*iter)[1] < smallest_y){ smallest_y = (*iter)[0]; }
    if ((*iter)[1] > biggest_y){ biggest_y = (*iter)[0]; }
  }
  Rectangle mbb = Rectangle(smallest_x, smallest_y, biggest_x, biggest_y);
  return mbb;
}

// cost_function calcualtes the weighted cost for cut 
double cost_function(Rectangle& B0, Rectangle& B1, double area_weight)
{
  // option 1: just area (area_weight = 1)
  // option 2: weighted area + perimeter  (area_weight < 1)
  assert(area_weight >= 0.0 && area_weight <= 1.0);
  double sum_area =  B0.area() + B1.area(); 
  double sum_margin = B0.margin() + B1.margin();
  return area_weight * sum_area + (1 - area_weight) * sum_margin; 
}

// basic_split returns a cut based on cost function 
std::pair<uint64_t, double> find_best_cut(
        std::vector<Point>::iterator begin,
        std::vector<Point>::iterator end,
        unsigned dimension,
        uint64_t M,
        double lowest_cost
) {
  uint64_t best_cut = 0; 
  double curr_lowest_cost = lowest_cost; 
  uint64_t range = std::ceil((end-begin)/M);
  for ( uint64_t i = 1; i < range; i++ ) {
    // cut is at [1 - M*i] | [M*i + 1, n]
    Rectangle B0 = min_bounding_box(begin, begin+i*M);
    Rectangle B1 = min_bounding_box(begin+i*M+1, end);
    // 1.0 is the weight of area cost 
    double cut_cost = cost_function(B0, B1, 1.0); 
    if (lowest_cost == 0 || cut_cost < lowest_cost){
      best_cut = i * M; 
      curr_lowest_cost = cut_cost;
    }
  }
  return std::make_pair(curr_lowest_cost, best_cut); 
}

// Basic Split: 
// 1. if n <= M (fill a leaf node or create a leaf node or create a branch node )
// 2. for each dimension and for each ordering
// 3.   for i from 1 to ceil(n/M) - 1
//        B0 = MBB(1 to i*M)
//        B1 = MBB(i*M + 1 to n)
// 4.     compute f(B0, B1)
// 5. split the input set based on i which has the highest f()

void basic_split_leaf(
        rstartreedisk::RStarTreeDisk<5, R_STAR_FANOUT> *tree,
        std::vector<std::vector<Point>::iterator> begins,
        std::vector<std::vector<Point>::iterator> ends,
        unsigned branch_factor,
        unsigned height, 
        uint64_t M, 
        pinned_node_ptr<rstartreedisk::LeafNode<5, R_STAR_FANOUT>> leaf_node)
{
  using LN = rstartreedisk::LeafNode<5, R_STAR_FANOUT>;
  using BN = rstartreedisk::BranchNode<5, R_STAR_FANOUT>;
  // count is number of points in this range 
  uint64_t count = ends[0] - begins[0]; 
  if (count <= M) {
    // there is no more split for this subrange 
    // if leaf node:
    // height = 0
      // if we are height 0, then add all points to the current node 
    for (auto iter = begins[0]; begins[0] != ends[0]; iter++) {
      leaf_node->addPoint(*iter);

    }
    return; 
  }

  // find the best one cut from all dimensions: 
  // starts with dimension x 
  unsigned dimension = 0; 
  uint64_t best_cut; 
  double lowest_cost;
  // checking x dimension 
  std::tie(lowest_cost, best_cut) = find_best_cut(begins[0], ends[0], 0, M, 0); 
  // checking y dimension 
  auto ret = find_best_cut(begins[1], ends[1], 1, M, lowest_cost);
  // if cut on y dimension is better: 
  if (ret.first < lowest_cost) {
    lowest_cost = ret.first; 
    best_cut = ret.second; 
    dimension = 1; 
  }

  // split the points into left and right based on the cut above 
  std::vector<std::vector<Point>::iterator> new_begins_left; 
  std::vector<std::vector<Point>::iterator> new_ends_left; 
  std::vector<std::vector<Point>::iterator> new_begins_right; 
  std::vector<std::vector<Point>::iterator> new_ends_right;
  std::vector<Point> points_copy_left;
  std::vector<Point> points_copy_right; 
  std::vector<Point>::iterator begin;
  std::vector<Point>::iterator end;
  if (dimension == 0){
    // cut is in dimension x 
    begin = begins[0];
    end = ends[0];
    std::copy(begin, begin+best_cut, back_inserter(points_copy_left)); 
    std::copy(begin+best_cut+1, end, back_inserter(points_copy_right)); 

    // the original list doesn't need to be sorted 
    std::vector<Point>::iterator begin_left = points_copy_left.begin();
    std::vector<Point>::iterator end_left = begin_left + points_copy_left.size();
    std::sort(begin_left, end_left, [](Point &l, Point &r) { return l[1] < r[1]; });

    std::vector<Point>::iterator begin_right = points_copy_right.begin();
    std::vector<Point>::iterator end_right = begin_right + points_copy_right.size();
    std::sort(begin_right, end_right, [](Point &l, Point &r) { return l[1] < r[1]; });

    new_begins_left.push_back(begin);
    new_begins_left.push_back(begin_left);
    new_ends_left.push_back(begin+best_cut);
    new_ends_left.push_back(end_left);
    
    new_begins_right.push_back(begin+best_cut +1);
    new_begins_right.push_back(begin_right);
    new_ends_right.push_back(end);
    new_ends_right.push_back(end_right);

    basic_split_leaf(tree, new_begins_left, new_ends_left, branch_factor, height, M, leaf_node);
    basic_split_leaf(tree, new_begins_right, new_ends_right, branch_factor, height, M, leaf_node);
    //basic_split(tree, new_begins_left, new_ends_left, branch_factor, height, M, branch_node, child_node);
    //basic_split(tree, new_begins_right, new_ends_right, branch_factor, height, M, branch_node, child_node);
  } else if(dimension == 1){
    // cut is in dimension y 
    begin = begins[1];
    end = ends[1];
    std::copy(begin, begin+best_cut, back_inserter(points_copy_left)); 
    std::copy(begin+best_cut+1, end, back_inserter(points_copy_right)); 

    // the original list doesn't need to be sorted 
    std::vector<Point>::iterator begin_left = points_copy_left.begin();
    std::vector<Point>::iterator end_left = begin_left + points_copy_left.size();
    std::sort(begin_left, end_left, [](Point &l, Point &r) { return l[0] < r[0]; });

    std::vector<Point>::iterator begin_right = points_copy_right.begin();
    std::vector<Point>::iterator end_right = begin_right + points_copy_right.size();
    std::sort(begin_right, end_right, [](Point &l, Point &r) { return l[0] < r[0]; });

    new_begins_left.push_back(begin_left);
    new_begins_left.push_back(begin);
    new_ends_left.push_back(end_left);
    new_ends_left.push_back(begin+best_cut);
    
    new_begins_right.push_back(begin_right);
    new_begins_right.push_back(begin+best_cut +1);
    new_ends_right.push_back(end_right);
    new_ends_right.push_back(end);

    basic_split_leaf(tree, new_begins_left, new_ends_left, branch_factor, height, M, leaf_node);
    basic_split_leaf(tree, new_begins_right, new_ends_right, branch_factor, height, M, leaf_node);
    //basic_split(tree, new_begins_left, new_ends_left, branch_factor, height, M, branch_node, child_node);
    //basic_split(tree, new_begins_right, new_ends_right, branch_factor, height, M, branch_node, child_node);
  } else{
    assert (dimension != 0 && dimension != 1);
  }

}

void basic_split_branch(
        rstartreedisk::RStarTreeDisk<5, R_STAR_FANOUT> *tree,
        std::vector<std::vector<Point>::iterator> begins,
        std::vector<std::vector<Point>::iterator> ends,
        unsigned branch_factor,
        unsigned height, 
        uint64_t M, 
        pinned_node_ptr<rstartreedisk::BranchNode<5, R_STAR_FANOUT>> branch_node)
{
        //NT tree_node){
  using LN = rstartreedisk::LeafNode<5, R_STAR_FANOUT>;
  using BN = rstartreedisk::BranchNode<5, R_STAR_FANOUT>;
  // count is number of points in this range 
  uint64_t count = ends[0] - begins[0]; 
  if (count <= M) {
    // there is no more split for this subrange 
    //if branch node: 
    tree_node_allocator *allocator = tree->node_allocator_.get();
    if (height == 1) {
      // next level should be leaf node 
      // allocate a new leaf node at next level: 
      auto alloc_data = allocator->create_new_tree_node<LN>(NodeHandleType(LEAF_NODE));
      new (&(*(alloc_data.first))) LN();
      auto leaf_node = alloc_data.first;
      tree_node_handle leaf_handle = alloc_data.second;

      // recursion 
      uint64_t new_M = pow(branch_factor, (std::ceil(log(count) / log(branch_factor)) - 1));
      basic_split_leaf(tree, begins, ends, branch_factor, height - 1, new_M, leaf_node);
      //basic_split(tree, begins, ends, branch_factor, height - 1, new_M, NULL, child_node);
      
      // ready to create the current branch 
      rstartreedisk::Branch b;
      b.child = leaf_handle; 
      b.boundingBox = leaf_node->boundingBox();
      branch_node->addBranchToNode(b);
      
    } else {
      // next level is still branch node:  
      // allocate a new branchnode at next level: 
      auto alloc_data = allocator->create_new_tree_node<BN>(NodeHandleType(BRANCH_NODE));
      new (&(*(alloc_data.first))) BN();
      auto child_node = alloc_data.first;
      tree_node_handle child_handle = alloc_data.second; 

      // recursion 
      uint64_t new_M = pow(branch_factor, (std::ceil(log(count) / log(branch_factor)) - 1));
      basic_split_branch(tree, begins, ends, branch_factor, height - 1, new_M, child_node);
      //basic_split(tree, begins, ends, branch_factor, height - 1, new_M, child_node, NULL);
      
      // ready to create the current branch 
      rstartreedisk::Branch b;
      b.child = child_handle; 
      b.boundingBox = child_node->boundingBox();
      branch_node->addBranchToNode(b);
    }
    return; 
  }

  // find the best one cut from all dimensions: 
  // starts with dimension x 
  unsigned dimension = 0; 
  uint64_t best_cut; 
  double lowest_cost;
  // checking x dimension 
  std::tie(lowest_cost, best_cut) = find_best_cut(begins[0], ends[0], 0, M, 0); 
  // checking y dimension 
  auto ret = find_best_cut(begins[1], ends[1], 1, M, lowest_cost);
  // if cut on y dimension is better: 
  if (ret.first < lowest_cost) {
    lowest_cost = ret.first; 
    best_cut = ret.second; 
    dimension = 1; 
  }

  // split the points into left and right based on the cut above 
  std::vector<std::vector<Point>::iterator> new_begins_left; 
  std::vector<std::vector<Point>::iterator> new_ends_left; 
  std::vector<std::vector<Point>::iterator> new_begins_right; 
  std::vector<std::vector<Point>::iterator> new_ends_right;
  std::vector<Point> points_copy_left;
  std::vector<Point> points_copy_right; 
  std::vector<Point>::iterator begin;
  std::vector<Point>::iterator end;
  if (dimension == 0){
    // cut is in dimension x 
    begin = begins[0];
    end = ends[0];
    std::copy(begin, begin+best_cut, back_inserter(points_copy_left)); 
    std::copy(begin+best_cut+1, end, back_inserter(points_copy_right)); 

    // the original list doesn't need to be sorted 
    std::vector<Point>::iterator begin_left = points_copy_left.begin();
    std::vector<Point>::iterator end_left = begin_left + points_copy_left.size();
    std::sort(begin_left, end_left, [](Point &l, Point &r) { return l[1] < r[1]; });

    std::vector<Point>::iterator begin_right = points_copy_right.begin();
    std::vector<Point>::iterator end_right = begin_right + points_copy_right.size();
    std::sort(begin_right, end_right, [](Point &l, Point &r) { return l[1] < r[1]; });

    new_begins_left.push_back(begin);
    new_begins_left.push_back(begin_left);
    new_ends_left.push_back(begin+best_cut);
    new_ends_left.push_back(end_left);
    
    new_begins_right.push_back(begin+best_cut +1);
    new_begins_right.push_back(begin_right);
    new_ends_right.push_back(end);
    new_ends_right.push_back(end_right);

    basic_split_branch(tree, new_begins_left, new_ends_left, branch_factor, height, M, branch_node);
    basic_split_branch(tree, new_begins_right, new_ends_right, branch_factor, height, M, branch_node);
    //basic_split(tree, new_begins_left, new_ends_left, branch_factor, height, M, branch_node, child_node);
    //basic_split(tree, new_begins_right, new_ends_right, branch_factor, height, M, branch_node, child_node);
  } else if(dimension == 1){
    // cut is in dimension y 
    begin = begins[1];
    end = ends[1];
    std::copy(begin, begin+best_cut, back_inserter(points_copy_left)); 
    std::copy(begin+best_cut+1, end, back_inserter(points_copy_right)); 

    // the original list doesn't need to be sorted 
    std::vector<Point>::iterator begin_left = points_copy_left.begin();
    std::vector<Point>::iterator end_left = begin_left + points_copy_left.size();
    std::sort(begin_left, end_left, [](Point &l, Point &r) { return l[0] < r[0]; });

    std::vector<Point>::iterator begin_right = points_copy_right.begin();
    std::vector<Point>::iterator end_right = begin_right + points_copy_right.size();
    std::sort(begin_right, end_right, [](Point &l, Point &r) { return l[0] < r[0]; });

    new_begins_left.push_back(begin_left);
    new_begins_left.push_back(begin);
    new_ends_left.push_back(end_left);
    new_ends_left.push_back(begin+best_cut);
    
    new_begins_right.push_back(begin_right);
    new_begins_right.push_back(begin+best_cut +1);
    new_ends_right.push_back(end_right);
    new_ends_right.push_back(end);

    basic_split_branch(tree, new_begins_left, new_ends_left, branch_factor, height, M, branch_node);
    basic_split_branch(tree, new_begins_right, new_ends_right, branch_factor, height, M, branch_node);
    //basic_split(tree, new_begins_left, new_ends_left, branch_factor, height, M, branch_node, child_node);
    //basic_split(tree, new_begins_right, new_ends_right, branch_factor, height, M, branch_node, child_node);
  } else{
    assert (dimension != 0 && dimension != 1);
  }

}

// Bulk load method of a tree with Top Down Greedy Split
// algorithm:
// 1. sort data on all dimensions individually with c orderings 
//    points: just 1 ordering available 
//    rectangles: (1) min coord (2) max coord (3) center 
// 2. run basic_split to generate branch nodes and leaf nodes
// Shirley's Notes: 
// this is an entry point 
// for now, this is done for 2 dimensions individually, should be genralized 
// to multi-dimensions later 

tree_node_handle tgs_load(
        rstartreedisk::RStarTreeDisk<5, R_STAR_FANOUT> *tree,
        std::vector<Point>::iterator begin,
        std::vector<Point>::iterator end,
        unsigned branch_factor
) {
  using LN = rstartreedisk::LeafNode<5, R_STAR_FANOUT>;
  using BN = rstartreedisk::BranchNode<5, R_STAR_FANOUT>;

  // allocate root branch node on buffer pool 
  // todo: check for small tree size which can fit into a leaf node? 
  tree_node_allocator *allocator = tree->node_allocator_.get();
  auto alloc_data = allocator->create_new_tree_node<BN>(NodeHandleType(BRANCH_NODE));
  new (&(*(alloc_data.first))) BN();
  auto root_node = alloc_data.first;
  tree_node_handle root_handle = alloc_data.second;
  
  // make a copy of the Points vector in order to keep one for each dimension
  std::vector<Point> points_copy;
  std::copy(begin, end, back_inserter(points_copy));  
  std::vector<Point>::iterator begin_y = points_copy.begin();
  std::vector<Point>::iterator end_y = begin_y + points_copy.size();
  // sort on x dimension 
  std::sort(begin, end, [](Point &l, Point &r) { return l[0] < r[0]; });
  // sort on y dimension 
  std::sort(begin_y, end_y, [](Point &l, Point &r) { return l[1] < r[1]; });
  std::vector<std::vector<Point>::iterator> begins; 
  std::vector<std::vector<Point>::iterator> ends; 
  // begins has [sorted_x.begin, sorted_y.begin]
  // ends has [sorted_x.end, sorted_y.end]
  begins.push_back(begin);
  begins.push_back(begin_y);
  ends.push_back(end);
  ends.push_back(end_y);

  // height of the Rtree root node (leaf has height 0)
  unsigned height = std::ceil(log(end - begin) / log(branch_factor)) - 1;
  uint64_t M = pow(branch_factor, (std::ceil(log(end-begin) / log(branch_factor)) - 1));
  // basic_split does the most work 
  basic_split_branch(tree, begins, ends, branch_factor, height, M, root_node); 
  return root_handle;
}



template <>
void bulk_load_tree(
    nirtreedisk::NIRTreeDisk<5, NIR_FANOUT, nirtreedisk::ExperimentalStrategy> *tree,
    std::map<std::string, size_t> &configU,
    std::vector<Point>::iterator begin,
    std::vector<Point>::iterator end,
    unsigned max_branch_factor
) {
  intersection_count = 0;
  auto tree_ptr = tree;
  uint64_t num_els = (end - begin);
  // Keep in mind there is a 0th level, so floor is correct
  uint64_t max_depth = std::floor(log(num_els) / log(max_branch_factor));

  std::cout << "Num els: " << num_els << std::endl;
  std::cout << "Max depth required: " << max_depth << std::endl;
  std::cout << "Size of NIR branch node: " <<
    sizeof(nirtreedisk::BranchNode<5, NIR_FANOUT, nirtreedisk::ExperimentalStrategy>) << std::endl;

  std::chrono::high_resolution_clock::time_point begin_time = std::chrono::high_resolution_clock::now();
  auto ret = quad_tree_style_load(
  tree_ptr, begin, end,
      max_branch_factor, 0, max_depth, nullptr
  );
  tree->root = ret.first;
  std::cout << "Out of line size: " << tree->node_allocator_.get()->out_of_line_nodes_size << std::endl;
  std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> delta = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - begin_time);

  std::cout << "Bulk loading NIRTree took: " << delta.count() << std::endl;
  std::cout << "Completed with " << intersection_count << " intersections" << std::endl;
  std::cout << "Total pages occupied: " << tree->node_allocator_->cur_page_ << std::endl;

  tree->write_metadata();
}

template <>
void bulk_load_tree(
    rstartreedisk::RStarTreeDisk<5, R_STAR_FANOUT> *tree,
    std::map<std::string, size_t> &configU,
    std::vector<Point>::iterator begin,
    std::vector<Point>::iterator end,
    unsigned max_branch_factor
) {
  uint64_t num_els = (end - begin);
  // Keep in mind there is a 0th level, so floor is correct
  uint64_t max_depth = std::floor(log(num_els) / log(max_branch_factor));
  uint64_t cur_depth = max_depth;

  std::cout << "Size of R* branch node: " << sizeof(rstartreedisk::BranchNode<5, R_STAR_FANOUT>) << std::endl;

  /* Start measuring bulk load time */
  std::chrono::high_resolution_clock::time_point begin_time = std::chrono::high_resolution_clock::now();

  if (configU["bulk_load_alg"] == STR) {
    std::cout << "Bulk-loading R* using Sort-Tile-Recursive..." << std::endl;
    std::vector<tree_node_handle> leaves = str_packing_leaf(tree, begin, end,max_branch_factor, cur_depth);
    cur_depth--;
    std::vector<tree_node_handle> branches = str_packing_branch(tree, leaves, max_branch_factor, cur_depth);
    cur_depth--;

    while (branches.size() > 1) {
      branches = str_packing_branch(tree, branches, max_branch_factor, cur_depth);
      cur_depth--;
    }

    tree->root = branches.at(0);
  } else if (configU["bulk_load_alg"] == TGS) {
    std::cout << "Bulk-loading R* using Top-Down Greedy Splitting..." << std::endl;
    tree->root = tgs_load(tree, begin, end, max_branch_factor);
  }

  std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> delta = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - begin_time);
  /* End measuring bulk load time */

  std::cout << "Bulk loading tree took: " << delta.count() << std::endl;
  std::cout << "Total pages occupied: " << tree->node_allocator_->cur_page_ << std::endl;
  tree->write_metadata();
}
