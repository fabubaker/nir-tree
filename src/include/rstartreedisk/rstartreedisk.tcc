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
    root = root_ptr->insert(givenPoint, hasReinsertedOnLevel);
    return;
  }
  auto root_ptr = get_branch_node(root);
  std::fill(hasReinsertedOnLevel.begin(), hasReinsertedOnLevel.end(), false);
  root = root_ptr->insert(givenPoint, hasReinsertedOnLevel);
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
  if (root.get_type() == LEAF_NODE) {
    auto root_ptr = get_leaf_node(root);
    root_ptr->printTree();
  }
  auto root_ptr = get_branch_node(root);
  root_ptr->printTree();
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
