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

template <typename T, typename LN, typename BN>
void fill_branch(
        T *treeRef,
        pinned_node_ptr<BN> branch_node,
        std::vector<std::pair<Point, tree_node_handle>> &node_point_pairs,
        uint64_t &offset,
        unsigned branch_factor,
        LN *leaf_type
) {
  std::vector<std::pair<Rectangle, tree_node_handle>> bb_and_handles;
  tree_node_allocator *allocator = treeRef->node_allocator_.get();

  // Add up to branch factor items to it
  for (uint64_t i = 0; i < branch_factor; i++) {
    Rectangle bbox;
    tree_node_handle child_handle = node_point_pairs[offset++].second;

    if (child_handle.get_type() == LEAF_NODE) {
      auto node = allocator->get_tree_node<LN>(child_handle);
      bbox = node->boundingBox();
    } else {
      auto node = allocator->get_tree_node<BN>(child_handle);
      bbox = node->boundingBox();
    }

    bb_and_handles.push_back(std::make_pair(bbox, child_handle));

    branch_node->addBranchToNode(bbox, child_handle);

    if (offset == node_point_pairs.size()) {
      break;
    }
  }
}

template <>
void fill_branch(
    nirtreedisk::NIRTreeDisk<NIR_MIN_FANOUT, NIR_MAX_FANOUT> *treeRef,
    pinned_node_ptr<nirtreedisk::BranchNode<NIR_MIN_FANOUT, NIR_MAX_FANOUT>> branch_node,
    std::vector<std::pair<Point, tree_node_handle>> &node_point_pairs,
    uint64_t &offset,
    unsigned branch_factor,
    nirtreedisk::LeafNode<NIR_MIN_FANOUT, NIR_MAX_FANOUT> *leaf_type) {
  using LN = nirtreedisk::LeafNode<NIR_MIN_FANOUT, NIR_MAX_FANOUT>;
  using BN = nirtreedisk::BranchNode<NIR_MIN_FANOUT, NIR_MAX_FANOUT>;

  std::vector<std::pair<IsotheticPolygon, tree_node_handle>> fixed_bb_and_handles;
  tree_node_allocator *allocator = treeRef->node_allocator_.get();

  // Add up to branch factor items to it
  for (uint64_t i = 0; i < branch_factor; i++) {
    tree_node_handle child_handle = node_point_pairs[offset++].second;
    Rectangle bbox;

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
    IsotheticPolygon constructed_poly = fixed_bb_and_handles.at(i).first;

    Rectangle boundingBox = constructed_poly.boundingBox;
    tree_node_handle child = fixed_bb_and_handles.at(i).second;

    // If the MBR has not been split into a polygon, don't keep it in the map.
    if (constructed_poly.basicRectangles.size() != 1) {
      treeRef->polygons.insert({child, constructed_poly});
    }

    branch_node->addBranchToNode(boundingBox, child);
  }
}

template <typename T, typename LN, typename BN>
std::vector<tree_node_handle> str_packing_branch(
    T *tree,
    std::vector<tree_node_handle> &child_nodes,
    unsigned branch_factor,
    unsigned cur_depth,
    LN *leaf_node_type,
    BN *branch_node_type
) {
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

template <typename T, typename LN, typename BN>
std::vector<tree_node_handle> str_packing_leaf(
    T *tree,
    std::vector<Point>::iterator begin,
    std::vector<Point>::iterator end,
    unsigned branch_factor,
    unsigned cur_depth,
    LN *leaf_node_type,
    BN *branch_node_type
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

std::vector<uint64_t> find_bounding_lines(
    std::vector<Point>::iterator start,
    std::vector<Point>::iterator stop,
    unsigned d, // dimension
    unsigned branch_factor,
    unsigned partitions,
    unsigned sub_partitions,
    unsigned length) {
  std::sort(start, stop, [d](const Point &p1, const Point &p2) {
                         return p1.orderedCompare(p2, d); });

  std::vector<uint64_t> lines;
  lines.reserve(partitions + 2); // We include x = 0 and x = count
  lines.push_back(0);

  uint64_t count = stop - start;
  unsigned optimally_filled_branch = pow(branch_factor, length); // Total number of children (points) starting from this level and below
  unsigned min_branches = std::ceil((double) count / (double) optimally_filled_branch); // Minimum branches at this level
  unsigned branches_per_partition = min_branches / partitions; // Divide branches at this level to partitions
  unsigned remaining_branches = min_branches % partitions; // Any leftover branches?

  unsigned partition_index = 0;
  // is partitions - 1 correct ? 
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

template <typename T, typename LN, typename BN>
void quad_tree_style_fill_branch(
  T *tree,
  std::vector<Point>::iterator start,
  std::vector<Point>::iterator stop,
  unsigned branch_factor,
  unsigned cur_level,
  unsigned cur_dimension,
  unsigned partitions,
  unsigned sub_partitions,
  LN *leaf_node_type,
  pinned_node_ptr<BN> branch_node
) {
  std::vector<uint64_t> line_candidates = find_bounding_lines(start, 
                                                              stop, 
                                                              cur_dimension, 
                                                              branch_factor, 
                                                              partitions, 
                                                              sub_partitions, 
                                                              cur_level);
  for (unsigned i = 0; i < line_candidates.size() - 1; i++) {
    
    // find next sub range
    unsigned start_ind = line_candidates.at(i);
    unsigned end_ind = line_candidates.at(i + 1); /* not inclusive */
    std::vector<Point>::iterator sub_start = start + start_ind;
    std::vector<Point>::iterator sub_stop = start + end_ind;

    // most inner dimension 
    if (cur_dimension == dimensions - 1) {
      if (sub_start == sub_stop) {
        // I think this can happen when we run out of points in
        // the lowest layer to split among children.
        return;
      }

      // This needs to return the box and handle.
      auto ret = quad_tree_style_load(
        tree, sub_start, sub_stop,
        branch_factor, cur_level - 1,
        leaf_node_type, branch_node_type
      );

      Rectangle boundingBox = ret.second;
      tree_node_handle child = ret.first;
      branch_node->addBranchToNode(boundingBox, child);

    } else {
      // recursive on next dimension
      // how to get sub_partitions for next dimension ??? 
      unsigned sub_sub_partions = 1;
      quad_tree_style_fill_branch(tree,
                                  sub_start,
                                  sub_stop,
                                  branch_factor,
                                  cur_level,
                                  cur_dimension + 1,
                                  sub_partitions,
                                  sub_sub_partions,
                                  leaf_node_type,
                                  branch_node);

    }
  }
}

template <typename T, typename LN, typename BN>
std::pair<tree_node_handle, Rectangle> quad_tree_style_load(
  T *tree,
  std::vector<Point>::iterator start,
  std::vector<Point>::iterator stop,
  unsigned branch_factor,
  unsigned cur_level,
  LN *leaf_node_type,
  BN *branch_node_type
) {
  uint64_t num_els = (stop - start);
  tree_node_allocator *allocator = tree->node_allocator_.get();
  if (cur_level == 0) {
    // Leaf Node
    if (num_els > branch_factor) {
      std::cout << "NUM ELS: " << num_els << std::endl;
    }

    assert(num_els <= branch_factor);

    num_els = (stop - start);
    auto alloc_data = allocator->create_new_tree_node<LN>(NodeHandleType(LEAF_NODE));
    new (&(*(alloc_data.first))) LN();

    auto leaf_node = alloc_data.first;
    auto leaf_handle = alloc_data.second;
    leaf_handle.set_level(cur_level);
    assert(leaf_handle.get_level() == 0);

    for (auto iter = start; iter != stop; iter++) {
      leaf_node->addPoint(*iter);
    }

    Rectangle bbox = leaf_node->boundingBox();
    return std::make_pair(leaf_handle, bbox);
  }
  // branch Node case
  // Return a tree node handle with pointers to all of its necessary
  // children.
  auto alloc_data = allocator->create_new_tree_node<BN>(NodeHandleType(BRANCH_NODE));
  new (&(*(alloc_data.first))) BN();

  auto branch_node = alloc_data.first;
  tree_node_handle branch_handle = alloc_data.second;
  branch_handle.set_level(cur_level);

  uint64_t partitions = std::ceil(sqrt(branch_factor));
  uint64_t sub_partitions = std::ceil(branch_factor / (float) partitions);

  unsigned start_dimension = 0;
  quad_tree_style_fill_branch(tree,
                              start,
                              stop,
                              branch_factor,
                              cur_level,
                              start_dimension,
                              partitions,
                              sub_partitions,
                              leaf_node_type,
                              branch_node);

  Rectangle bbox = branch_node->boundingBox();
  return std::make_pair(branch_handle, bbox);
}

template <>
std::pair<tree_node_handle, Rectangle> quad_tree_style_load(
  nirtreedisk::NIRTreeDisk<NIR_MIN_FANOUT, NIR_MAX_FANOUT> *tree,
  std::vector<Point>::iterator start,
  std::vector<Point>::iterator stop,
  unsigned branch_factor,
  unsigned cur_level,
  nirtreedisk::LeafNode<NIR_MIN_FANOUT, NIR_MAX_FANOUT> *leaf_node_type,
  nirtreedisk::BranchNode<NIR_MIN_FANOUT, NIR_MAX_FANOUT> *branch_node_type
) {
  using LN = nirtreedisk::LeafNode<NIR_MIN_FANOUT, NIR_MAX_FANOUT>;
  using BN = nirtreedisk::BranchNode<NIR_MIN_FANOUT, NIR_MAX_FANOUT>;

  uint64_t num_els = (stop - start);
  tree_node_allocator *allocator = tree->node_allocator_.get();
  if (cur_level == 0) {
    if (num_els > branch_factor) {
      std::cout << "NUM ELS: " << num_els << std::endl;
    }

    assert(num_els <= branch_factor);

    num_els = (stop - start);
    auto alloc_data = allocator->create_new_tree_node<LN>(NodeHandleType(LEAF_NODE));
    new (&(*(alloc_data.first))) LN();

    auto leaf_node = alloc_data.first;
    auto leaf_handle = alloc_data.second;
    leaf_handle.set_level(cur_level);
    assert(leaf_handle.get_level() == 0);

    for (auto iter = start; iter != stop; iter++) {
      leaf_node->addPoint(*iter);
    }

    Rectangle bbox = leaf_node->boundingBox();
    return std::make_pair(leaf_handle, bbox);
  }

  // Return a tree node handle with pointers to all of its necessary
  // children.
  auto alloc_data = allocator->create_new_tree_node<BN>(NodeHandleType(BRANCH_NODE));
  new (&(*(alloc_data.first))) BN();

  auto branch_node = alloc_data.first;
  tree_node_handle branch_handle = alloc_data.second;
  branch_handle.set_level(cur_level);

  uint64_t partitions = std::ceil(sqrt(branch_factor));
  uint64_t sub_partitions = std::ceil(branch_factor / (float) partitions);

  unsigned start_dimension = 0;
  quad_tree_style_fill_branch(tree,
                              start,
                              stop,
                              branch_factor,
                              cur_level,
                              start_dimension,
                              partitions,
                              sub_partitions,
                              leaf_node_type,
                              branch_node);

  if (child_handle.get_level() != 0) {
    std::vector<std::pair<IsotheticPolygon, tree_node_handle>> fixed_bb_and_handles;
    //branch_handles.reserve(NIR_MAX_FANOUT);

    for (unsigned i = 0; i < branch_node->cur_offset_; i++) {
      tree_node_handle child_handle = branch_node->entries.at(i).child;
      Rectangle bbox = branch_node->entries.at(i).boundingBox;
      // Produce a polygon using bbox.
      IsotheticPolygon ip = IsotheticPolygon(bbox);
      fixed_bb_and_handles.push_back(std::make_pair(ip, child_handle));
    }

    // make all branches polygon disjoint
    for (unsigned i = 0; i < fixed_bb_and_handles.size(); i++) {
      for (unsigned j = i + 1; j < fixed_bb_and_handles.size(); j++) {
        IsotheticPolygon &poly_a = fixed_bb_and_handles.at(i);
        IsotheticPolygon &poly_b = fixed_bb_and_handles.at(j);
        if (poly_a.intersectsPolygon(poly_b)) {
          std::vector<Rectangle> &existing_rects_a = poly_a.basicRectangles;
          std::vector<Rectangle> &existing_rects_b = poly_b.basicRectangles; 
          // It is imperative that the references here are set carefully.
          // make_all_rects_disjoint updates the second and fourth arguments,
          // so these should point to the basicRectangles vector of the polygons.
          make_all_rects_disjoint(
              treeRef,
              existing_rects_a,
              poly_a.second,
              existing_rects_b,
              poly_b.second
          );
          intersection_count += 1;
          // These may have been updated. Update the metadata.
          poly_a.recomputeBoundingBox();
          poly_b.recomputeBoundingBox();
        }
      }
    }

    // save polygon in map and update boundingBox
    for(unsigned i = 0; i < fixed_bb_and_handles.size(); i++) {
      IsotheticPolygon constructed_poly = fixed_bb_and_handles.at(i).first;
      Rectangle boundingBox = constructed_poly.boundingBox;
      tree_node_handle child = fixed_bb_and_handles.at(i).second;
      
      // If the MBR has not been split into a polygon, don't keep it in the map.
      if (constructed_poly.basicRectangles.size() != 1) {
        tree->polygons.insert({child, constructed_poly});
      }
      branch_node->entries.at(i).boundingBox = constructed_poly.boundingBox;
    }

  #ifndef NDEBUG
    // Double check non-intersection - Really inefficient, but I don't see a better way
    // of doing this.
    for (uint64_t i = 0; i < fixed_bb_and_handles.size(); i++) {
      for (uint64_t j = i + 1; j < fixed_bb_and_handles.size(); j++) {
        std::vector<Rectangle> &existing_rects_a =
                fixed_bb_and_handles.at(i).first.basicRectangles;
        std::vector<Rectangle> &existing_rects_b =
                fixed_bb_and_handles.at(j).first.basicRectangles;
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
  } // if (child_handle.get_level() != 0)

  Rectangle bbox = branch_node->boundingBox();
  return std::make_pair(branch_handle, bbox);
}

// cost_function calculates the weighted cost for a cut
// for the TGS bulk-loading algorithm.
double cost_function(Rectangle& B0, Rectangle& B1, double area_weight)
{
  // option 1: just area (area_weight = 1)
  // option 2: weighted area + perimeter  (area_weight < 1)
  assert(area_weight >= 0.0 && area_weight <= 1.0);
  double sum_area =  B0.area() + B1.area(); 
  double sum_margin = B0.margin() + B1.margin();
  return area_weight * sum_area + (1 - area_weight) * sum_margin; 
}

// find_best_cut returns a cut based on an arbitrary cost function
std::pair<double, uint64_t> find_best_cut(
        std::vector<Point>::iterator begin,
        std::vector<Point>::iterator end,
        uint64_t M,
        double lowest_cost
) {
  uint64_t best_cut = 0; 
  double curr_lowest_cost = lowest_cost; 
  uint64_t num_cuts = std::ceil((end - begin) / (double) M) - 1;

  for (uint64_t i = 1; i <= num_cuts; i++) {
    // Cut is applied at [1, M * i] | [M * i + 1, n]
    Rectangle B0 = Rectangle(begin, begin + i * M);
    Rectangle B1 = Rectangle(begin + i * M, end);

    // For now, we set area_weight to 1.0
    double area_weight = 1.0;
    double cut_cost = cost_function(B0, B1, 1.0); 

    if (cut_cost < curr_lowest_cost){
      best_cut = i * M; 
      curr_lowest_cost = cut_cost;
    }
  }

  return std::make_pair(curr_lowest_cost, best_cut);
}

// Basic TGS Split:
// 1. if n <= M (fill a leaf node or create a leaf node or create a branch node)
// 2. for each dimension and for each ordering
// 3.   for i from 1 to ceil(n/M) - 1
//        B0 = MBB of [1, M * i]
//        B1 = MBB of [M * i + 1, n]
// 4.     compute cost_function(B0, B1)
// 5. split the input set based on i which has the lowest cost
template <typename T, typename LN, typename BN>
void basic_split_leaf(
        T *tree,
        std::vector<std::vector<Point>::iterator> begins,
        std::vector<std::vector<Point>::iterator> ends,
        unsigned branch_factor,
        unsigned height, 
        uint64_t M, 
        pinned_node_ptr<LN> leaf_node,
        BN *branch_node_type
) {
  // count is number of points in this range
  uint64_t count = ends[0] - begins[0];

  if (count <= M) {
    // There are no more splits for this subrange, and we're a leaf node,
    // so add all points to the current node
    for (auto iter = begins[0]; iter != ends[0]; iter++) {
      leaf_node->addPoint(*iter);
    }

    return; 
  }

  // find the best one cut from all dimensions: 
  // starts with dimension x 
  unsigned dimension = 0; 
  uint64_t best_cut; 
  double lowest_cost = std::numeric_limits<double>::max();

  // Checking x dimension...
  auto begin_x = begins[0];
  auto end_x = ends[0];
  std::tie(lowest_cost, best_cut) = find_best_cut(begin_x, end_x, M, lowest_cost);

  // Checking y dimension...
  auto begin_y = begins[1];
  auto end_y = ends[1];
  auto ret = find_best_cut(begin_y, end_y, M, lowest_cost);

  // if cut on y dimension is better:
  if (ret.first < lowest_cost) {
    lowest_cost = ret.first; 
    best_cut = ret.second; 
    dimension = 1; 
  }

  // Split the points into left and right based on the cut above
  std::vector<std::vector<Point>::iterator> new_begins_left; 
  std::vector<std::vector<Point>::iterator> new_ends_left; 
  std::vector<std::vector<Point>::iterator> new_begins_right; 
  std::vector<std::vector<Point>::iterator> new_ends_right;
  std::vector<Point> points_left;
  std::vector<Point> points_right;

  if (dimension == 0) {
    // The best cut is in dimension x
    auto begin_x_left = begins[0];
    auto end_x_left = begin_x_left + best_cut;
    auto begin_x_right = end_x_left;
    auto end_x_right = ends[0];

    // Create copies of the points vector to sort in the y-dimension
    std::copy(begin_x_left, end_x_left, back_inserter(points_left));
    std::copy(begin_x_right, end_x_right, back_inserter(points_right));

    // The original vector is already sorted in the x-dimension, now we sort in the y-dimension
    auto begin_y_left = points_left.begin();
    auto end_y_left = points_left.end();
    std::sort(begin_y_left, end_y_left, [](Point &l, Point &r) { return l[1] < r[1]; });

    auto begin_y_right = points_right.begin();
    auto end_y_right = points_right.end();
    std::sort(begin_y_right, end_y_right, [](Point &l, Point &r) { return l[1] < r[1]; });

    new_begins_left.push_back(begin_x_left);
    new_begins_left.push_back(begin_y_left);
    new_ends_left.push_back(end_x_left);
    new_ends_left.push_back(end_y_left);

    new_begins_right.push_back(begin_x_right);
    new_begins_right.push_back(begin_y_right);
    new_ends_right.push_back(end_x_right);
    new_ends_right.push_back(end_y_right);

    basic_split_leaf(tree, new_begins_left, new_ends_left, branch_factor, height, M, leaf_node, branch_node_type);
    basic_split_leaf(tree, new_begins_right, new_ends_right, branch_factor, height, M, leaf_node, branch_node_type);
  } else if (dimension == 1){
    // The best cut is in dimension y
    auto begin_y_left = begins[1];
    auto end_y_left = begin_y_left + best_cut;
    auto begin_y_right = end_y_left;
    auto end_y_right = ends[1];

    // Create copies of the points vector to sort in the x-dimension
    std::copy(begin_y_left, end_y_left, back_inserter(points_left));
    std::copy(begin_y_right, end_y_right, back_inserter(points_right));

    // The original vector is already sorted in the y-dimension, we sort in the x-dimension
    auto begin_x_left = points_left.begin();
    auto end_x_left = points_left.end();
    std::sort(begin_x_left, end_x_left, [](Point &l, Point &r) { return l[0] < r[0]; });

    auto begin_x_right = points_right.begin();
    auto end_x_right = begin_x_right + points_right.size();
    std::sort(begin_x_right, end_x_right, [](Point &l, Point &r) { return l[0] < r[0]; });

    new_begins_left.push_back(begin_x_left);
    new_begins_left.push_back(begin_y_left);
    new_ends_left.push_back(end_x_left);
    new_ends_left.push_back(end_y_left);

    new_begins_right.push_back(begin_x_right);
    new_begins_right.push_back(begin_y_right);
    new_ends_right.push_back(end_x_right);
    new_ends_right.push_back(end_y_right);

    basic_split_leaf(tree, new_begins_left, new_ends_left, branch_factor, height, M, leaf_node, branch_node_type);
    basic_split_leaf(tree, new_begins_right, new_ends_right, branch_factor, height, M, leaf_node, branch_node_type);
  } else {
    assert (dimension != 0 && dimension != 1);
  }
}

template <typename T, typename LN, typename BN>
void basic_split_branch(
        T *tree,
        std::vector<std::vector<Point>::iterator> begins,
        std::vector<std::vector<Point>::iterator> ends,
        unsigned branch_factor,
        unsigned height, 
        uint64_t M, 
        pinned_node_ptr<BN> branch_node,
        LN *leaf_node_type
) {
  // count is number of points in this range
  uint64_t count = ends[0] - begins[0]; 

  if (count <= M) {
    // There are no more splits for this subrange
    tree_node_allocator *allocator = tree->node_allocator_.get();

    if (height == 1) {
      // The next level is a leaf node, so we allocate a new leaf node
      auto alloc_data = allocator->create_new_tree_node<LN>(NodeHandleType(LEAF_NODE));
      new (&(*(alloc_data.first))) LN();
      auto leaf_node = alloc_data.first;
      tree_node_handle leaf_handle = alloc_data.second;
      unsigned next_level = height - 1;
      // set level for leaf node handle
      leaf_handle.set_level(next_level);

      // Recurse onto next level with M = branch factor
      uint64_t new_M = branch_factor;
      basic_split_leaf(
          tree, begins, ends, branch_factor, height - 1, new_M, leaf_node, (BN *) nullptr
      );
      
      // Create the current branch after recursing
      assert(leaf_handle.get_level() == 0);
      branch_node->addBranchToNode(leaf_node->boundingBox(), leaf_handle);
    } else {
      // The next level is a branch node, so we allocate a branch node
      auto alloc_data = allocator->create_new_tree_node<BN>(NodeHandleType(BRANCH_NODE));
      new (&(*(alloc_data.first))) BN();
      auto child_node = alloc_data.first;
      tree_node_handle child_handle = alloc_data.second;
      unsigned next_level = height - 1;
      // set level for branch node handle
      child_handle.set_level(next_level);

      // Recurse onto next level
      uint64_t new_M = pow(branch_factor, (std::ceil(log(count) / log(branch_factor)) - 1));
      basic_split_branch(
          tree, begins, ends, branch_factor, height - 1, new_M, child_node, leaf_node_type
      );

      // Create the current branch after recursing
      assert(child_handle.get_level() == (height - 1));
      branch_node->addBranchToNode(child_node->boundingBox(), child_handle);
    }

    return; 
  }

  // Find the best one split across all dimensions:
  // Start with dimension x
  unsigned dimension = 0; 
  uint64_t best_cut; 
  double lowest_cost = std::numeric_limits<double>::max();

  // Checking x dimension...
  auto begin_x = begins[0];
  auto end_x = ends[0];
  std::tie(lowest_cost, best_cut) = find_best_cut(begin_x, end_x, M, lowest_cost);

  // Checking y dimension...
  auto begin_y = begins[1];
  auto end_y = ends[1];
  auto ret = find_best_cut(begin_y, end_y, M, lowest_cost);

  // if cut on y dimension is better:
  if (ret.first < lowest_cost) {
    lowest_cost = ret.first; 
    best_cut = ret.second; 
    dimension = 1; 
  }

  // Split the points into left and right based on the cut above
  std::vector<std::vector<Point>::iterator> new_begins_left; 
  std::vector<std::vector<Point>::iterator> new_ends_left; 
  std::vector<std::vector<Point>::iterator> new_begins_right; 
  std::vector<std::vector<Point>::iterator> new_ends_right;
  std::vector<Point> points_left;
  std::vector<Point> points_right;

  if (dimension == 0) {
    // The best cut is in dimension x
    auto begin_x_left = begins[0];
    auto end_x_left = begin_x_left + best_cut;
    auto begin_x_right = end_x_left;
    auto end_x_right = ends[0];

    // Create copies of the points vector to sort in the y-dimension
    std::copy(begin_x_left, end_x_left, back_inserter(points_left));
    std::copy(begin_x_right, end_x_right, back_inserter(points_right));

    // The original vector is already sorted in the x-dimension, now we sort in the y-dimension
    auto begin_y_left = points_left.begin();
    auto end_y_left = points_left.end();
    std::sort(begin_y_left, end_y_left, [](Point &l, Point &r) { return l[1] < r[1]; });

    auto begin_y_right = points_right.begin();
    auto end_y_right = points_right.end();
    std::sort(begin_y_right, end_y_right, [](Point &l, Point &r) { return l[1] < r[1]; });

    new_begins_left.push_back(begin_x_left);
    new_begins_left.push_back(begin_y_left);
    new_ends_left.push_back(end_x_left);
    new_ends_left.push_back(end_y_left);

    new_begins_right.push_back(begin_x_right);
    new_begins_right.push_back(begin_y_right);
    new_ends_right.push_back(end_x_right);
    new_ends_right.push_back(end_y_right);

    basic_split_branch(tree, new_begins_left, new_ends_left, branch_factor, height, M, branch_node, leaf_node_type);
    basic_split_branch(tree, new_begins_right, new_ends_right, branch_factor, height, M, branch_node, leaf_node_type);
  } else if (dimension == 1) {
    // The best cut is in dimension y
    auto begin_y_left = begins[1];
    auto end_y_left = begin_y_left + best_cut;
    auto begin_y_right = end_y_left;
    auto end_y_right = ends[1];

    // Create copies of the points vector to sort in the x-dimension
    std::copy(begin_y_left, end_y_left, back_inserter(points_left));
    std::copy(begin_y_right, end_y_right, back_inserter(points_right));

    // The original vector is already sorted in the y-dimension, we sort in the x-dimension
    auto begin_x_left = points_left.begin();
    auto end_x_left = points_left.end();
    std::sort(begin_x_left, end_x_left, [](Point &l, Point &r) { return l[0] < r[0]; });

    auto begin_x_right = points_right.begin();
    auto end_x_right = begin_x_right + points_right.size();
    std::sort(begin_x_right, end_x_right, [](Point &l, Point &r) { return l[0] < r[0]; });

    new_begins_left.push_back(begin_x_left);
    new_begins_left.push_back(begin_y_left);
    new_ends_left.push_back(end_x_left);
    new_ends_left.push_back(end_y_left);

    new_begins_right.push_back(begin_x_right);
    new_begins_right.push_back(begin_y_right);
    new_ends_right.push_back(end_x_right);
    new_ends_right.push_back(end_y_right);

    basic_split_branch(tree, new_begins_left, new_ends_left, branch_factor, height, M, branch_node, leaf_node_type);
    basic_split_branch(tree, new_begins_right, new_ends_right, branch_factor, height, M, branch_node, leaf_node_type);
  } else {
    assert (dimension != 0 && dimension != 1);
  }
}

// Bulk load method of an R-tree with Top Down Greedy Split
// algorithm:
// 1. sort data on all dimensions individually with c orderings 
//    points: just 1 ordering available 
//    rectangles: (1) min coord (2) max coord (3) center 
// 2. run basic_split to generate branch nodes and leaf nodes
template <typename T, typename LN, typename BN>
tree_node_handle tgs_load(
        T *tree,
        std::vector<Point>::iterator begin,
        std::vector<Point>::iterator end,
        unsigned branch_factor,
        LN *leaf_node_type,
        BN *branch_node_type
) {
  // TODO: The number of points may be small enough to fit into a single leaf node.
  //       We do not handle this case.
  tree_node_allocator *allocator = tree->node_allocator_.get();
  auto alloc_data = allocator->create_new_tree_node<BN>(NodeHandleType(BRANCH_NODE));
  new (&(*(alloc_data.first))) BN();
  auto root_node = alloc_data.first;
  tree_node_handle root_handle = alloc_data.second;
  
  // Copy the Points vector in order to sort it separately for each dimension
  std::vector<Point> points;
  std::copy(begin, end, back_inserter(points));

  auto begin_x = begin;
  auto end_x = end;
  auto begin_y = points.begin();
  auto end_y = points.end();

  // Sort on x dimension
  std::sort(begin_x, end_x, [](Point &l, Point &r) { return l[0] < r[0]; });
  // Sort on y dimension
  std::sort(begin_y, end_y, [](Point &l, Point &r) { return l[1] < r[1]; });

  // begins has [sorted_x.begin, sorted_y.begin]
  // ends has [sorted_x.end, sorted_y.end]
  std::vector<std::vector<Point>::iterator> begins;
  std::vector<std::vector<Point>::iterator> ends;

  begins.push_back(begin_x);
  begins.push_back(begin_y);
  ends.push_back(end_x);
  ends.push_back(end_y);

  // Height of the R-tree root node (leaf has height 0)
  unsigned height = std::ceil(log(end - begin) / log(branch_factor)) - 1;
  unsigned root_level = height; 
  uint64_t M = pow(branch_factor, height);

  basic_split_branch(tree, begins, ends, branch_factor, height, M, root_node, leaf_node_type);
  // Set level for root node
  root_handle.set_level(root_level);
  return root_handle;
}

template <class TR>
void testLevels(TR *tree, tree_node_handle root, unsigned root_level){
  std::stack<std::pair<tree_node_handle, unsigned>> context; 
  context.push(std::make_pair(root, root_level));
  while(not context.empty()){
    auto handle_level = context.top();
    context.pop();
    tree_node_handle current_handle = handle_level.first;
    uint8_t current_level = handle_level.second;
    if(current_handle.get_level() != current_level){
      std::cout << "current_handle.get_level() is " << current_handle.get_level() << std::endl;
      std::cout << "current_level is expected to be " << current_level << std::endl;
    }
    if(current_handle.get_type() == LEAF_NODE) {
      if(current_handle.get_level() != 0){
        std::cout << "current_handle.get_level() is " << current_handle.get_level() << std::endl;
        std::cout << "Leaf Node is expected to be 0" << std::endl;
      }
    } else {
      auto current_node = tree->get_branch_node(current_handle);
      for (int i = 0; i < current_node->cur_offset_; ++i){
        auto bi = current_node->entries[i];
        context.push(std::make_pair(bi.child, current_level - 1));
      }
    }
  }
  std::cout << "Pass Level Test" << std::endl;
}

template <>
void bulk_load_tree(
    nirtreedisk::NIRTreeDisk<NIR_MIN_FANOUT, NIR_MAX_FANOUT> *tree,
    std::map<std::string, size_t> &configU,
    std::vector<Point>::iterator begin,
    std::vector<Point>::iterator end,
    unsigned max_branch_factor,
    nirtreedisk::LeafNode<NIR_MIN_FANOUT, NIR_MAX_FANOUT> *leaf_node_type,
    nirtreedisk::BranchNode<NIR_MIN_FANOUT, NIR_MAX_FANOUT> *branch_node_type
) {
  intersection_count = 0;
  auto tree_ptr = tree;
  uint64_t num_els = (end - begin);
  // Leaf is at 0th level
  uint64_t max_depth = std::ceil(log(num_els) / log(max_branch_factor)) - 1;
  // QTS is top down
  uint64_t root_level = max_depth;

  std::cout << "Num els: " << num_els << std::endl;
  std::cout << "Max depth required: " << max_depth << std::endl;
  std::cout << "Size of NIR branch node: " <<
    sizeof(nirtreedisk::BranchNode<NIR_MIN_FANOUT, NIR_MAX_FANOUT>) << std::endl;
  
  std::cout << "Bulk-loading NIRTree using Quad Tree Style Load..." << std::endl;
  std::chrono::high_resolution_clock::time_point begin_time = std::chrono::high_resolution_clock::now();
  auto ret = quad_tree_style_load(
    tree_ptr, begin, end, max_branch_factor, root_level,
    leaf_node_type, branch_node_type
  );
  tree->root = ret.first;

  std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> delta = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - begin_time);

#ifndef NDEBUG
  // run level test
  testLevels(tree, tree->root, root_level);
#endif

  std::cout << "Bulk loading NIRTree took: " << delta.count() << std::endl;
  std::cout << "Completed with " << intersection_count << " intersections" << std::endl;
  std::cout << "Total pages occupied: " << tree->node_allocator_->cur_page_ << std::endl;

  tree->write_metadata();
}

template <typename T, typename LN, typename BN>
void bulk_load_tree(
    T *tree,
    std::map<std::string, size_t> &configU,
    std::vector<Point>::iterator begin,
    std::vector<Point>::iterator end,
    unsigned max_branch_factor,
    LN *leaf_node_type,
    BN *branch_node_type
) {
  uint64_t num_els = (end - begin);
  // Leaf is at 0th level
  uint64_t max_depth = std::ceil(log(num_els) / log(max_branch_factor)) - 1;

  std::cout << "Num els: " << num_els << std::endl;
  std::cout << "Max depth required: " << max_depth << std::endl;

  /* Start measuring bulk load time */
  std::chrono::high_resolution_clock::time_point begin_time = std::chrono::high_resolution_clock::now();

  switch(configU["bulk_load_alg"]) {
    case STR: {
      // STR is bottom up
      uint64_t cur_level = 0;
      std::cout << "Bulk-loading using Sort-Tile-Recursive..." << std::endl;

      // Bulk load the leaf level first
      std::vector<tree_node_handle> leaves = str_packing_leaf(
        tree, begin, end, max_branch_factor, cur_level,
        (LN *) nullptr, (BN *) nullptr
      );
      cur_level++;

      // Bulk load all branch levels next
      std::vector<tree_node_handle> branches = str_packing_branch(
        tree, leaves, max_branch_factor, cur_level,
        (LN *) nullptr, (BN *) nullptr
      );
      cur_level++;

      while (branches.size() > 1) {
        assert(cur_level <= max_depth);
        branches = str_packing_branch(
          tree, branches, max_branch_factor, cur_level,
          (LN *) nullptr, (BN *) nullptr
        );
        cur_level++;
      }

      tree->root = branches.at(0);

      break;
    }
    case TGS: {
      std::cout << "Bulk-loading using Top-Down Greedy Splitting..." << std::endl;
      tree->root = tgs_load(
        tree, begin, end, max_branch_factor,
        (LN *) nullptr, (BN *) nullptr
      );

      break;
    }
    case QTS: {
      // QTS is top down
      uint64_t root_level = max_depth;
      std::cout << "Bulk-loading using Quad Tree Style Load..." << std::endl;
      auto ret = quad_tree_style_load(
        tree, begin, end, max_branch_factor, root_level,
        leaf_node_type, branch_node_type
      );
      tree->root = ret.first;

      break;
    }
    default: {
      std::cout << "Bulk-loading algorithm is not recognized..." << std::endl;
      exit(1);
    }
  }

  std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> delta = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - begin_time);
  /* End measuring bulk load time */

#ifndef NDEBUG
  // run level test
  // Level of the R-tree root node (leaf has height 0)
  unsigned root_level = max_depth;
  testLevels(tree, tree->root, root_level);
#endif

  std::cout << "Bulk loading tree took: " << delta.count() << std::endl;
  std::cout << "Total pages occupied: " << tree->node_allocator_->cur_page_ << std::endl;

  tree->hasReinsertedOnLevel.resize(max_depth + 1, false);
  tree->write_metadata();
}

template <typename T>
void sequential_insert_tree(
    T *tree,
    std::map<std::string, size_t> &configU,
    std::vector<Point>::iterator begin,
    std::vector<Point>::iterator end,
    unsigned max_branch_factor
) {
  // begin is inclusive, end is exclusive
  uint64_t num_els = (end - begin);
  std::cout << "Num els to insert: " << num_els << std::endl;
  uint64_t total_insert = 0;
  uint64_t print_count = pow(10, int(log10(num_els)) - 1);
  std::chrono::high_resolution_clock::time_point begin_time = std::chrono::high_resolution_clock::now();
  std::chrono::high_resolution_clock::time_point section_begin_time = begin_time;
  for(auto iter = begin ; iter < end; iter++){
      tree->insert(*iter);
      total_insert ++;
    if (total_insert % print_count == 0) {
      std::chrono::high_resolution_clock::time_point section_end_time = std::chrono::high_resolution_clock::now();
      auto delta =  std::chrono::duration_cast<std::chrono::duration<double>>(section_end_time - section_begin_time);
      std::cout << "Finished insertion for " << total_insert << " points with " << delta.count() << "s..."<< std::endl;
      section_begin_time = section_end_time;
    }
  }
  std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> delta = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - begin_time);

  std::cout << "Sequentially Inserting "<< total_insert << " points took: " << delta.count() << std::endl;
  std::cout << "Total pages occupied now: " << tree->node_allocator_->get_total_pages_occupied() << std::endl;
  tree->write_metadata();
}

/* Generalized template definitions are supposed to be exposed in the header file in order
 * for the compiler to generate code for specific instances. Alternatively, we can explicitly
 * instantiate functions below for only those types that we care about and still keep everything
 * in the .cpp file. */

/* sequential_insert_tree */
template void sequential_insert_tree(
        nirtreedisk::NIRTreeDisk<NIR_MIN_FANOUT, NIR_MAX_FANOUT> *tree,
        std::map<std::string, size_t> &configU,
        std::vector<Point>::iterator begin,
        std::vector<Point>::iterator end,
        unsigned max_branch_factor
);

template void sequential_insert_tree(
        rstartreedisk::RStarTreeDisk<R_STAR_MIN_FANOUT, R_STAR_MAX_FANOUT> *tree,
        std::map<std::string, size_t> &configU,
        std::vector<Point>::iterator begin,
        std::vector<Point>::iterator end,
        unsigned max_branch_factor
);

/* bulk_load_tree */
template void bulk_load_tree(
        rstartreedisk::RStarTreeDisk<R_STAR_MIN_FANOUT, R_STAR_MAX_FANOUT> *tree,
        std::map<std::string, size_t> &configU,
        std::vector<Point>::iterator begin,
        std::vector<Point>::iterator end,
        unsigned max_branch_factor,
        rstartreedisk::LeafNode<R_STAR_MIN_FANOUT, R_STAR_MAX_FANOUT> *leaf_node_type,
        rstartreedisk::BranchNode<R_STAR_MIN_FANOUT, R_STAR_MAX_FANOUT> *branch_node_type
);
