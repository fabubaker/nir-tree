#include <catch2/catch.hpp>
#include <nirtreedisk/nirtreedisk.h>
#include <storage/page.h>
#include <util/geometry.h>
#include <iostream>
#include <unistd.h>

#define DefaultLeafNodeType nirtreedisk::LeafNode<3,7,nirtreedisk::LineMinimizeDownsplits>
#define DefaultBranchNodeType nirtreedisk::BranchNode<3,7,nirtreedisk::LineMinimizeDownsplits>
#define DefaulTreeType nirtreedisk::NIRTreeDisk<3,7, nirtreedisk::LineMinimizeDownsplits>
#define DefaultTemplateParams 3,7, nirtreedisk::LineMinimizeDownsplits

#define MeanBalancedNNType nirtreedisk::Node<3,7,nirtreedisk::LineMinimizeDistanceFromMean>
#define MeanBalancedLeafNodeType nirtreedisk::LeafNode<3,7,nirtreedisk::LineMinimizeDistanceFromMean>
#define MeanBalancedBranchNodeType nirtreedisk::BranchNode<3,7,nirtreedisk::LineMinimizeDistanceFromMean>
#define MeanBalancedTreeType nirtreedisk::NIRTreeDisk<3,7, nirtreedisk::LineMinimizeDistanceFromMean>

static nirtreedisk::Branch createBranchEntry(
    const InlineBoundedIsotheticPolygon &boundingBox,
    tree_node_handle child
) {
    nirtreedisk::Branch b(boundingBox, child);
	return b;
}

static nirtreedisk::Branch createBranchEntry(
    tree_node_handle poly_handle,
    tree_node_handle child
) {
    nirtreedisk::Branch b(poly_handle, child);
	return b;
}

static tree_node_handle
createFullLeafNode(DefaulTreeType &tree, tree_node_handle parent, Point p=Point::atOrigin)
{
    // Allocate new node
    auto alloc_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                NodeHandleType( LEAF_NODE ) );
    tree_node_handle node_handle = alloc_data.second;
    auto node = alloc_data.first;
    new (&(*node)) DefaultLeafNodeType( &tree, tree_node_handle(nullptr),
            node_handle, 0 );

    std::vector<bool> hasReinsertedOnLevel = {false};
	for (unsigned i = 0; i < 7; i++) {
		node->insert(p, hasReinsertedOnLevel );
	}

    node->parent = parent;
	return node_handle;
}


TEST_CASE("NIRTreeDisk: testBoundingBox")
{
	// Test set one

    unlink( "nirdiskbacked.txt" );
    {
        DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );

        // The root node is a LeafNode by default. To avoid the
        // annoyance of making it into a Branch node, just make a fresh
        // node usingthe tree's allocator.

        auto root_alloc_data =
            tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>(
                    NodeHandleType( BRANCH_NODE ) );
        new (&(*root_alloc_data.first)) DefaultBranchNodeType( &tree,
                 tree_node_handle( nullptr ), root_alloc_data.second, 1 );

        pinned_node_ptr<DefaultBranchNodeType> rootNode =
            root_alloc_data.first;
        REQUIRE( rootNode->parent == tree_node_handle( nullptr ) );
        tree_node_handle root = root_alloc_data.second;
        
        // Make a bunch of Leaf Nodes
        std::pair<pinned_node_ptr<DefaultLeafNodeType>, tree_node_handle> alloc_data =
            tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                    NodeHandleType( LEAF_NODE ) );
        tree_node_handle child0 = alloc_data.second;
        new (&(*alloc_data.first)) DefaultLeafNodeType( &tree, root,
                child0, 0);
        nirtreedisk::Branch entry =
            createBranchEntry( InlineBoundedIsotheticPolygon(), child0);

        std::get<InlineBoundedIsotheticPolygon>( entry.boundingPoly ).push_polygon_to_disk(
                    IsotheticPolygon( Rectangle(8.0, 1.0, 12.0, 5.0) ) );

        rootNode->addBranchToNode( entry );

        alloc_data =
            tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                    NodeHandleType( LEAF_NODE ) );
        tree_node_handle child1 = alloc_data.second;
        new (&(*alloc_data.first)) DefaultLeafNodeType( &tree, root,
                child1, 0 );
        entry = createBranchEntry( InlineBoundedIsotheticPolygon(), child1);

        std::get<InlineBoundedIsotheticPolygon>( entry.boundingPoly ).push_polygon_to_disk(
                    IsotheticPolygon( Rectangle(12.0, -4.0, 16.0, -2.0) ) );

        rootNode->addBranchToNode( entry );

        alloc_data =
            tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                    NodeHandleType( LEAF_NODE ) );
        tree_node_handle child2 = alloc_data.second;
        new (&(*alloc_data.first)) DefaultLeafNodeType( &tree, root,
                child2, 0 );
        entry = createBranchEntry( InlineBoundedIsotheticPolygon( ), child2 );

        std::get<InlineBoundedIsotheticPolygon>( entry.boundingPoly ).push_polygon_to_disk(
                    IsotheticPolygon( Rectangle(8.0, -6.0, 10.0, -4.0) ) );
        rootNode->addBranchToNode( entry );

        REQUIRE( rootNode->cur_offset_ == 3 );
        REQUIRE( rootNode->boundingBox() == Rectangle(8.0, -6.0, 16.0, 5.0) );

    }

    unlink( "nirdiskbacked.txt" );

    {
        // Test set two
        DefaulTreeType tree(4096 * 5, "nirdiskbacked.txt" );

        auto alloc_root_data =
            tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>(
                    NodeHandleType( BRANCH_NODE ) );
        new (&(*alloc_root_data.first)) DefaultBranchNodeType( &tree,
                 tree_node_handle( nullptr ), alloc_root_data.second, 1);
        tree_node_handle root = alloc_root_data.second;
        auto rootNode = alloc_root_data.first;

        auto alloc_data =
            tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                    NodeHandleType( LEAF_NODE ) );
        tree_node_handle child0 = alloc_data.second;
        new (&(*alloc_data.first)) DefaultLeafNodeType( &tree, child0,
                root, 0 );

        nirtreedisk::Branch entry = createBranchEntry( InlineBoundedIsotheticPolygon(), child0); 

        std::get<InlineBoundedIsotheticPolygon>( entry.boundingPoly ).push_polygon_to_disk(
                    IsotheticPolygon( Rectangle(8.0, 12.0, 10.0, 14.0) ) );
        rootNode->addBranchToNode( entry );

        alloc_data =
            tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                    NodeHandleType( LEAF_NODE ) );
        tree_node_handle child1 = alloc_data.second;
        new (&(*alloc_data.first)) DefaultLeafNodeType( &tree, root,
                child1, 0 );

        entry = createBranchEntry(
                InlineBoundedIsotheticPolygon(), child1);
        std::get<InlineBoundedIsotheticPolygon>( entry.boundingPoly ).push_polygon_to_disk(
                    IsotheticPolygon( Rectangle(10.0, 12.0, 12.0, 14.0) ) );
        rootNode->addBranchToNode( entry );

        alloc_data =
            tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                    NodeHandleType( LEAF_NODE ) );
        tree_node_handle child2 = alloc_data.second;
        new (&(*alloc_data.first)) DefaultLeafNodeType( &tree, root,
                child2, 0 );
        entry = createBranchEntry( InlineBoundedIsotheticPolygon(), child2);
        std::get<InlineBoundedIsotheticPolygon>( entry.boundingPoly ).push_polygon_to_disk(
                    IsotheticPolygon( Rectangle(12.0, 12.0, 14.0, 14.0) ) );

        rootNode->addBranchToNode( entry );
            
        REQUIRE( rootNode->cur_offset_ ==  3 );

        REQUIRE(rootNode->boundingBox() == Rectangle(8.0, 12.0, 14.0, 14.0));
    }

    unlink( "nirdiskbacked.txt" );
}
TEST_CASE("NIRTreeDisk: testUpdateBoundingBox") {

    unlink( "nirdiskbacked.txt" );
	DefaulTreeType tree(4096*5, "nirdiskbacked.txt");

    auto alloc_root_data =
        tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>(
                NodeHandleType( BRANCH_NODE ) );
    new (&(*alloc_root_data.first)) DefaultBranchNodeType( &tree,
             tree_node_handle( nullptr ), alloc_root_data.second, 1 );
    tree_node_handle root = alloc_root_data.second;
    auto parentNode = alloc_root_data.first;

    auto alloc_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                NodeHandleType( LEAF_NODE ) );
    tree_node_handle child0 = alloc_data.second;
    new (&(*alloc_data.first)) DefaultLeafNodeType( &tree, root, child0,
            0 );
    auto child0Node = alloc_data.first;
	REQUIRE( child0Node->parent == root );

    auto entry = createBranchEntry(InlineBoundedIsotheticPolygon(), child0);

    std::get<InlineBoundedIsotheticPolygon>( entry.boundingPoly ).push_polygon_to_disk(
                    IsotheticPolygon( Rectangle(8.0, -6.0, 10.0, -4.0) ) );

    parentNode->addBranchToNode( entry );

    alloc_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                NodeHandleType( LEAF_NODE ) );
    tree_node_handle child1 = alloc_data.second;
    new (&(*alloc_data.first)) DefaultLeafNodeType( &tree, root,
            child1,0 );
    auto child1Node = alloc_data.first;
    REQUIRE( child1Node->parent == root );

    entry = createBranchEntry(
            InlineBoundedIsotheticPolygon(), child1);
    std::get<InlineBoundedIsotheticPolygon>( entry.boundingPoly ).push_polygon_to_disk(
                    IsotheticPolygon( Rectangle(12.0, -4.0, 16.0, -2.0) ) );
    parentNode->addBranchToNode(entry);

    alloc_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                NodeHandleType( LEAF_NODE ) );
    tree_node_handle child2 = alloc_data.second;
    new (&(*alloc_data.first)) DefaultLeafNodeType( &tree, root, child2,
            0);
    auto child2Node = alloc_data.first;
    REQUIRE( child2Node->parent == root );
    
    entry = createBranchEntry(
            InlineBoundedIsotheticPolygon(), child2);
    std::get<InlineBoundedIsotheticPolygon>( entry.boundingPoly ).push_polygon_to_disk(
                    IsotheticPolygon( Rectangle(10.0, 12.0, 12.0, 14.0) ) );
    parentNode->addBranchToNode( entry );

    alloc_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                NodeHandleType( LEAF_NODE ) );
    tree_node_handle child3 = alloc_data.second;
    new (&(*alloc_data.first)) DefaultLeafNodeType( &tree, root, child3,
            0);
    auto child3Node = alloc_data.first;
    REQUIRE( child3Node->parent == root );
    entry = createBranchEntry( InlineBoundedIsotheticPolygon(), child3);
    std::get<InlineBoundedIsotheticPolygon>( entry.boundingPoly ).push_polygon_to_disk(
                    IsotheticPolygon( Rectangle(12.0, 12.0, 14.0, 14.0) ) );

    parentNode->addBranchToNode( entry );
    REQUIRE( parentNode->cur_offset_ == 4 );

    InlineBoundedIsotheticPolygon stack_poly;
    IsotheticPolygon loc_poly(Rectangle(3.0, 3.0, 5.0,5.0) );
    stack_poly.push_polygon_to_disk( loc_poly );
	parentNode->updateBranch(child3, stack_poly);

	auto &b = parentNode->entries[3];
    auto &poly = std::get<InlineBoundedIsotheticPolygon>( b.boundingPoly );
	REQUIRE(poly.materialize_polygon().boundingBox == Rectangle(3.0, 3.0, 5.0, 5.0));
    unlink( "nirdiskbacked.txt" );
}

TEST_CASE( "NIRTreeDisk: testRemoveChild" ) {
    unlink( "nirdiskbacked.txt" );
	DefaulTreeType tree(4096*5, "nirdiskbacked.txt");

    auto alloc_root_data =
        tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>(
                NodeHandleType( BRANCH_NODE ) );
    new (&(*alloc_root_data.first)) DefaultBranchNodeType( &tree,
             tree_node_handle( nullptr ), alloc_root_data.second, 1 );
    tree_node_handle root = alloc_root_data.second;
    auto parentNode = alloc_root_data.first;

    auto alloc_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                NodeHandleType( LEAF_NODE ) );
    tree_node_handle child0 = alloc_data.second;
    new (&(*alloc_data.first)) DefaultLeafNodeType( &tree, root, child0,
            0
            );
    auto child0Node = alloc_data.first;
	REQUIRE( child0Node->parent == root );
	parentNode->addBranchToNode( createBranchEntry(
            InlineBoundedIsotheticPolygon(
                Rectangle(8.0, -6.0, 10.0, -4.0)), child0) );

    alloc_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                NodeHandleType( LEAF_NODE ) );
    tree_node_handle child1 = alloc_data.second;
    new (&(*alloc_data.first)) DefaultLeafNodeType( &tree, root, child1,
            0 );
    auto child1Node = alloc_data.first;
    REQUIRE( child1Node->parent == root );
    parentNode->addBranchToNode( createBranchEntry(
            InlineBoundedIsotheticPolygon(Rectangle(12.0, -4.0, 16.0,
                    -2.0)), child1) );

    alloc_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                NodeHandleType( LEAF_NODE ) );
    tree_node_handle child2 = alloc_data.second;
    new (&(*alloc_data.first)) DefaultLeafNodeType( &tree, root, child2,
            0
            );
    auto child2Node = alloc_data.first;
    REQUIRE( child2Node->parent == root );
    parentNode->addBranchToNode( createBranchEntry(
            InlineBoundedIsotheticPolygon(Rectangle(10.0, 12.0, 12.0,
                    14.0)), child2) );

    alloc_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                NodeHandleType( LEAF_NODE ) );
    tree_node_handle child3 = alloc_data.second;
    new (&(*alloc_data.first)) DefaultLeafNodeType( &tree, root, child3,
            0);
    auto child3Node = alloc_data.first;
    REQUIRE( child3Node->parent == root );
    parentNode->addBranchToNode( createBranchEntry(
            InlineBoundedIsotheticPolygon(Rectangle(12.0, 12.0, 14.0, 14.0)), child3) );

    REQUIRE( parentNode->cur_offset_ ==  4 );

	// Remove one of the children
	parentNode->removeBranch(child3);
    REQUIRE( parentNode->cur_offset_ == 3 );

    unlink( "nirdiskbacked.txt" );

}

TEST_CASE("NIRTreeDisk: testRemoveData")
{

    unlink( "nirdiskbacked.txt" );

	// Setup a rtree::Node with some data
	DefaulTreeType tree( 4096, "nirdiskbacked.txt" );
	tree_node_handle root = tree.root;
    auto parentNode = tree.get_leaf_node( root );

    parentNode->addPoint( Point(9.0, -5.0) );
	parentNode->addPoint( Point(14.0, -3.0) );
	parentNode->addPoint( Point(11.0, 13.0) );
	parentNode->addPoint( Point(13.0, 13.0) );

	REQUIRE(parentNode->cur_offset_ == 4);
	// Remove some of the data
	parentNode->removePoint( Point(13.0, 13.0) );

	// Test the removal
	REQUIRE(parentNode->cur_offset_ == 3);

    unlink( "nirdiskbacked.txt" );
}

TEST_CASE("NIRTreeDisk: testFindLeaf")
{
	// Create rtree::Nodes
    unlink( "nirdiskbacked.txt" );


    // Need a bunch of pages so we don't run out of memory while
    // everything is pinned
	DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );
    auto root_alloc_data =
        tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>(
                NodeHandleType( BRANCH_NODE ) );
    new (&(*root_alloc_data.first)) DefaultBranchNodeType( &tree,
             tree_node_handle( nullptr ), root_alloc_data.second, 2 );

    pinned_node_ptr<DefaultBranchNodeType> rootNode =
        root_alloc_data.first;
    REQUIRE( rootNode->parent == tree_node_handle( nullptr ) );
    tree_node_handle root = root_alloc_data.second;

    auto alloc_branch_data =
        tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>(
                NodeHandleType( BRANCH_NODE ) );
    auto leftNode = alloc_branch_data.first;
    tree_node_handle left = alloc_branch_data.second;
    new (&(*leftNode)) DefaultBranchNodeType( &tree, root, left, 1 );

    alloc_branch_data =
        tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>(
                NodeHandleType( BRANCH_NODE ) );
    auto rightNode = alloc_branch_data.first;
    tree_node_handle right = alloc_branch_data.second;
    new (&(*rightNode)) DefaultBranchNodeType( &tree, root, right, 1 );

    tree_node_handle leftChild0 = createFullLeafNode(
            tree, left, Point(8.5, 12.5));
    tree_node_handle leftChild1 = createFullLeafNode(tree, left,
            Point(11.0,15.0));
    tree_node_handle leftChild2 = createFullLeafNode(tree, left,
            Point(13.5,13.5));
    leftNode->addBranchToNode(
        createBranchEntry( InlineBoundedIsotheticPolygon(Rectangle(8.0, 12.0,
                    nextafter(10.0, DBL_MAX),
                    nextafter(14.0, DBL_MAX))), leftChild0 ) );
    leftNode->addBranchToNode( createBranchEntry(
            InlineBoundedIsotheticPolygon(Rectangle(7.0, 12.0,
                    nextafter(12.0, DBL_MAX),
                    nextafter(15.0, DBL_MAX))), leftChild1 ) );
    leftNode->addBranchToNode( createBranchEntry(
            InlineBoundedIsotheticPolygon(Rectangle(12.0, 12.0,
                    nextafter(14.0, DBL_MAX),
                    nextafter(14.0, DBL_MAX))), leftChild2 ) );
    rootNode->addBranchToNode( createBranchEntry(
            InlineBoundedIsotheticPolygon(Rectangle(7.0, 12.0,
                    nextafter(14.0, DBL_MAX),
                    nextafter(15.0, DBL_MAX))), left ) );

    tree_node_handle rightChild0 = createFullLeafNode(tree,right,
            Point(7.0, 3.0) );
    tree_node_handle rightChild1 = createFullLeafNode(tree,right,
            Point(13.0,-3.0));
    tree_node_handle rightChild2 = createFullLeafNode(tree,right);
    REQUIRE( rightNode->parent == root );
    rightNode->addBranchToNode( createBranchEntry(
            InlineBoundedIsotheticPolygon(Rectangle(7.0, 1.0,
                    nextafter(12.0, DBL_MAX),
                    nextafter(5.0, DBL_MAX))), rightChild0 ));
    rightNode->addBranchToNode( createBranchEntry(
            InlineBoundedIsotheticPolygon(Rectangle(12.0, -4.0,
                    nextafter(16.0, DBL_MAX),
                    nextafter(-2.0, DBL_MAX))), rightChild1 ));
    rightNode->addBranchToNode( createBranchEntry(
            InlineBoundedIsotheticPolygon(Rectangle(8.0, -6.0,
                    nextafter(10.0, DBL_MAX),
                    nextafter(-4.0, DBL_MAX))), rightChild2 ));
    rootNode->addBranchToNode( createBranchEntry(
            InlineBoundedIsotheticPolygon(Rectangle(7.0, -6.0,
                    nextafter(16.0, DBL_MAX),
                    nextafter(5.0, DBL_MAX))), right ));

	// Test that we get the correct child for the given point
	REQUIRE(rightChild1 == rootNode->findLeaf(Point(13.0, -3.0)));
	REQUIRE(leftChild0 == rootNode->findLeaf(Point(8.5, 12.5)));
	REQUIRE(leftChild2 == rootNode->findLeaf(Point(13.5, 13.5)));
	REQUIRE(rightChild0 == rootNode->findLeaf(Point(7.0, 3.0)));
	REQUIRE(leftChild1 == rootNode->findLeaf(Point(11.0, 15.0)));

    unlink( "nirdiskbacked.txt" );
}

TEST_CASE("NIRTreeDisk: testFindLeaf2")
{
	// Setup the tree

	// Cluster 4, n = 7
	// (-10, -2), (-12, -3), (-11, -3), (-10, -3), (-9, -3), (-7, -3), (-10, -5)
	// Organized into two rtree::Nodes
    unlink( "nirdiskbacked.txt" );

	DefaulTreeType tree( 4096 * 10, "nirdiskbacked.txt" );
    auto alloc_root_data =
        tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>(
                NodeHandleType( BRANCH_NODE ) );
    new (&(*alloc_root_data.first)) DefaultBranchNodeType( &tree,
             tree_node_handle( nullptr ), alloc_root_data.second, 2 );
    tree_node_handle root = alloc_root_data.second;
    auto rootNode= alloc_root_data.first;


    auto alloc_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                NodeHandleType( LEAF_NODE ) );
    auto cluster4aNode = alloc_data.first;
    tree_node_handle cluster4a = alloc_data.second;
    new (&(*cluster4aNode)) DefaultLeafNodeType( &tree,
            tree_node_handle(nullptr), cluster4a, 0 );
    cluster4aNode->addPoint( Point(-10.0, -2.0) );
    cluster4aNode->addPoint( Point(-12.0, -3.0) );
    cluster4aNode->addPoint( Point(-11.0, -3.0) );
    cluster4aNode->addPoint( Point(-10.0, -3.0) );
	
    alloc_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                NodeHandleType( LEAF_NODE ) );
    auto cluster4bNode = alloc_data.first;
    tree_node_handle cluster4b = alloc_data.second;
    new (&(*cluster4bNode)) DefaultLeafNodeType( &tree,
            tree_node_handle(nullptr), cluster4b, 0 );

	cluster4bNode->addPoint( Point(-9.0, -3.0) );
	cluster4bNode->addPoint( Point(-7.0, -3.0) );
	cluster4bNode->addPoint( Point(-10.0, -5.0) );

    auto alloc_branch_data =
        tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>(
                NodeHandleType( BRANCH_NODE ) );
    auto cluster4Node = alloc_branch_data.first;
    tree_node_handle cluster4 = alloc_branch_data.second;
    new (&(*cluster4Node)) DefaultBranchNodeType( &tree, root, cluster4,
            1);

	cluster4aNode->parent = cluster4;
	cluster4Node->addBranchToNode(
        createBranchEntry( InlineBoundedIsotheticPolygon(cluster4aNode->boundingBox()), cluster4a) );
	cluster4bNode->parent = cluster4;
	cluster4Node->addBranchToNode(
        createBranchEntry( InlineBoundedIsotheticPolygon(cluster4bNode->boundingBox()), cluster4b));

	// Cluster 5, n = 16
	// (-14.5, -13), (-14, -13), (-13.5, -13.5), (-15, -14), (-14, -14), (-13, -14), (-12, -14),
	// (-13.5, -16), (-15, -14.5), (-14, -14.5), (-12.5, -14.5), (-13.5, -15.5), (-15, -15),
	// (-14, -15), (-13, -15), (-12, -15)

    alloc_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                NodeHandleType( LEAF_NODE ) );
    auto cluster5aNode = alloc_data.first;
    tree_node_handle cluster5a = alloc_data.second;
    new (&(*cluster5aNode)) DefaultLeafNodeType( &tree,
            tree_node_handle(nullptr), cluster5a, 0 );
	cluster5aNode->addPoint( Point(-14.5, -13.0) );
	cluster5aNode->addPoint( Point(-14.0, -13.0) );
	cluster5aNode->addPoint( Point(-13.5, -13.5) );
	cluster5aNode->addPoint( Point(-15.0, -14.0) );

    alloc_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                NodeHandleType( LEAF_NODE ) );
    auto cluster5bNode = alloc_data.first;
    tree_node_handle cluster5b = alloc_data.second;
    new (&(*cluster5bNode)) DefaultLeafNodeType( &tree,
            tree_node_handle(nullptr), cluster5b, 0 );
	cluster5bNode->addPoint( Point(-14.0, -14.0) );
	cluster5bNode->addPoint( Point(-13.0, -14.0) );
	cluster5bNode->addPoint( Point(-12.0, -14.0) );
	cluster5bNode->addPoint( Point(-13.5, -16.0) );

    alloc_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                NodeHandleType( LEAF_NODE ) );
    auto cluster5cNode = alloc_data.first;
    tree_node_handle cluster5c = alloc_data.second;
    new (&(*cluster5cNode)) DefaultLeafNodeType( &tree,
                tree_node_handle(nullptr), cluster5c, 0 );

	cluster5cNode->addPoint( Point(-15.0, -14.5) );
	cluster5cNode->addPoint( Point(-14.0, -14.5) );
	cluster5cNode->addPoint( Point(-12.5, -14.5) );
	cluster5cNode->addPoint( Point(-13.5, -15.5) );

    alloc_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                NodeHandleType( LEAF_NODE ) );
    auto cluster5dNode = alloc_data.first;
    tree_node_handle cluster5d = alloc_data.second;
    new (&(*cluster5dNode)) DefaultLeafNodeType( &tree, tree_node_handle( nullptr
                ), cluster5d, 0 );
	cluster5dNode->addPoint( Point(-15.0, -15.0));
	cluster5dNode->addPoint( Point(-14.0, -15.0));
	cluster5dNode->addPoint( Point(-13.0, -15.0));
	cluster5dNode->addPoint( Point(-12.0, -15.0));
	cluster5dNode->addPoint( Point(-15.0, -15.0));

    alloc_branch_data =
        tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>(
                NodeHandleType( BRANCH_NODE ) );
    auto cluster5Node = alloc_branch_data.first;
    tree_node_handle cluster5 = alloc_branch_data.second;
    new (&(*cluster5Node)) DefaultBranchNodeType( &tree, root, cluster5,
            1 );

	cluster5aNode->parent = cluster5;
	cluster5Node->addBranchToNode(
        createBranchEntry(
            InlineBoundedIsotheticPolygon(cluster5aNode->boundingBox()),cluster5a));
	cluster5bNode->parent = cluster5;
	cluster5Node->addBranchToNode(
        createBranchEntry(
            InlineBoundedIsotheticPolygon(cluster5bNode->boundingBox()), cluster5b));
	cluster5cNode->parent = cluster5;
	cluster5Node->addBranchToNode(
        createBranchEntry(
            InlineBoundedIsotheticPolygon(cluster5cNode->boundingBox()), cluster5c));
	cluster5dNode->parent = cluster5;
	cluster5Node->addBranchToNode(
        createBranchEntry(
            InlineBoundedIsotheticPolygon(cluster5dNode->boundingBox()), cluster5d));

	// Root
	rootNode->addBranchToNode(
        createBranchEntry(
            InlineBoundedIsotheticPolygon(cluster4Node->boundingBox()), cluster4));
	rootNode->addBranchToNode(
        createBranchEntry(
            InlineBoundedIsotheticPolygon(cluster5Node->boundingBox()), cluster5));

	// Test finding leaves
	REQUIRE(rootNode->findLeaf(Point(-11.0, -3.0)) == cluster4a);
	REQUIRE(rootNode->findLeaf(Point(-9.0, -3.0)) == cluster4b);
	REQUIRE(rootNode->findLeaf(Point(-13.5, -13.5)) == cluster5a);
	REQUIRE(rootNode->findLeaf(Point(-12.0, -14.0)) == cluster5b);
	REQUIRE(rootNode->findLeaf(Point(-12.5, -14.5)) == cluster5c);
	REQUIRE(rootNode->findLeaf(Point(-13.0, -15.0)) == cluster5d);

    unlink( "nirdiskbacked.txt" );
}

TEST_CASE("NIRTreeDisk: testFindLeaf2 ON DISK")
{
	// Setup the tree

	// Cluster 4, n = 7
	// (-10, -2), (-12, -3), (-11, -3), (-10, -3), (-9, -3), (-7, -3), (-10, -5)
	// Organized into two rtree::Nodes
    unlink( "nirdiskbacked.txt" );
    tree_node_handle root;
    tree_node_handle cluster4a;
    tree_node_handle cluster4b;
    tree_node_handle cluster4;
    tree_node_handle cluster5a;
    tree_node_handle cluster5b;
    tree_node_handle cluster5c;
    tree_node_handle cluster5d;
    tree_node_handle cluster5;

    {
        DefaulTreeType tree( 4096 * 10, "nirdiskbacked.txt" );
        auto alloc_root_data =
            tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>(
                    NodeHandleType( BRANCH_NODE ) );
        new (&(*alloc_root_data.first)) DefaultBranchNodeType( &tree,
                 tree_node_handle( nullptr ), alloc_root_data.second, 2 );
        root = alloc_root_data.second;
        auto rootNode= alloc_root_data.first;


        auto alloc_data =
            tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                    NodeHandleType( LEAF_NODE ) );
        auto cluster4aNode = alloc_data.first;
        cluster4a = alloc_data.second;
        new (&(*cluster4aNode)) DefaultLeafNodeType( &tree,
                tree_node_handle(nullptr), cluster4a, 0 );
        cluster4aNode->addPoint( Point(-10.0, -2.0) );
        cluster4aNode->addPoint( Point(-12.0, -3.0) );
        cluster4aNode->addPoint( Point(-11.0, -3.0) );
        cluster4aNode->addPoint( Point(-10.0, -3.0) );
        
        alloc_data =
            tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                    NodeHandleType( LEAF_NODE ) );
        auto cluster4bNode = alloc_data.first;
        cluster4b = alloc_data.second;
        new (&(*cluster4bNode)) DefaultLeafNodeType( &tree,
                tree_node_handle(nullptr), cluster4b, 0 );

        cluster4bNode->addPoint( Point(-9.0, -3.0) );
        cluster4bNode->addPoint( Point(-7.0, -3.0) );
        cluster4bNode->addPoint( Point(-10.0, -5.0) );

        auto alloc_branch_data =
            tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>(
                    NodeHandleType( BRANCH_NODE ) );
        auto cluster4Node = alloc_branch_data.first;
        cluster4 = alloc_branch_data.second;
        new (&(*cluster4Node)) DefaultBranchNodeType( &tree, root,
                cluster4, 1);

        cluster4aNode->parent = cluster4;
        cluster4Node->addBranchToNode(
            createBranchEntry( InlineBoundedIsotheticPolygon(cluster4aNode->boundingBox()), cluster4a) );
        cluster4bNode->parent = cluster4;
        cluster4Node->addBranchToNode(
            createBranchEntry( InlineBoundedIsotheticPolygon(cluster4bNode->boundingBox()), cluster4b));

        // Cluster 5, n = 16
        // (-14.5, -13), (-14, -13), (-13.5, -13.5), (-15, -14), (-14, -14), (-13, -14), (-12, -14),
        // (-13.5, -16), (-15, -14.5), (-14, -14.5), (-12.5, -14.5), (-13.5, -15.5), (-15, -15),
        // (-14, -15), (-13, -15), (-12, -15)

        alloc_data =
            tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                    NodeHandleType( LEAF_NODE ) );
        auto cluster5aNode = alloc_data.first;
        cluster5a = alloc_data.second;
        new (&(*cluster5aNode)) DefaultLeafNodeType( &tree,
                tree_node_handle(nullptr), cluster5a, 0 );
        cluster5aNode->addPoint( Point(-14.5, -13.0) );
        cluster5aNode->addPoint( Point(-14.0, -13.0) );
        cluster5aNode->addPoint( Point(-13.5, -13.5) );
        cluster5aNode->addPoint( Point(-15.0, -14.0) );

        alloc_data =
            tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                    NodeHandleType( LEAF_NODE ) );
        auto cluster5bNode = alloc_data.first;
        cluster5b = alloc_data.second;
        new (&(*cluster5bNode)) DefaultLeafNodeType( &tree,
                tree_node_handle(nullptr), cluster5b, 0 );
        cluster5bNode->addPoint( Point(-14.0, -14.0) );
        cluster5bNode->addPoint( Point(-13.0, -14.0) );
        cluster5bNode->addPoint( Point(-12.0, -14.0) );
        cluster5bNode->addPoint( Point(-13.5, -16.0) );

        alloc_data =
            tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                    NodeHandleType( LEAF_NODE ) );
        auto cluster5cNode = alloc_data.first;
        cluster5c = alloc_data.second;
        new (&(*cluster5cNode)) DefaultLeafNodeType( &tree,
                    tree_node_handle(nullptr), cluster5c, 0 );

        cluster5cNode->addPoint( Point(-15.0, -14.5) );
        cluster5cNode->addPoint( Point(-14.0, -14.5) );
        cluster5cNode->addPoint( Point(-12.5, -14.5) );
        cluster5cNode->addPoint( Point(-13.5, -15.5) );

        alloc_data =
            tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                    NodeHandleType( LEAF_NODE ) );
        auto cluster5dNode = alloc_data.first;
        cluster5d = alloc_data.second;
        new (&(*cluster5dNode)) DefaultLeafNodeType( &tree, tree_node_handle( nullptr
                    ), cluster5d, 0 );
        cluster5dNode->addPoint( Point(-15.0, -15.0));
        cluster5dNode->addPoint( Point(-14.0, -15.0));
        cluster5dNode->addPoint( Point(-13.0, -15.0));
        cluster5dNode->addPoint( Point(-12.0, -15.0));
        cluster5dNode->addPoint( Point(-15.0, -15.0));

        alloc_branch_data =
            tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>(
                    NodeHandleType( BRANCH_NODE ) );
        auto cluster5Node = alloc_branch_data.first;
        cluster5 = alloc_branch_data.second;
        new (&(*cluster5Node)) DefaultBranchNodeType( &tree, root,
                cluster5, 0 );

        cluster5aNode->parent = cluster5;
        cluster5Node->addBranchToNode(
            createBranchEntry(
                InlineBoundedIsotheticPolygon(cluster5aNode->boundingBox()),cluster5a));
        cluster5bNode->parent = cluster5;
        cluster5Node->addBranchToNode(
            createBranchEntry(
                InlineBoundedIsotheticPolygon(cluster5bNode->boundingBox()), cluster5b));
        cluster5cNode->parent = cluster5;
        cluster5Node->addBranchToNode(
            createBranchEntry(
                InlineBoundedIsotheticPolygon(cluster5cNode->boundingBox()), cluster5c));
        cluster5dNode->parent = cluster5;
        cluster5Node->addBranchToNode(
            createBranchEntry(
                InlineBoundedIsotheticPolygon(cluster5dNode->boundingBox()), cluster5d));

        // Root
        rootNode->addBranchToNode(
            createBranchEntry(
                InlineBoundedIsotheticPolygon(cluster4Node->boundingBox()), cluster4));
        rootNode->addBranchToNode(
            createBranchEntry(
                InlineBoundedIsotheticPolygon(cluster5Node->boundingBox()), cluster5));

        // Test finding leaves
        REQUIRE(rootNode->findLeaf(Point(-11.0, -3.0)) == cluster4a);
        REQUIRE(rootNode->findLeaf(Point(-9.0, -3.0)) == cluster4b);
        REQUIRE(rootNode->findLeaf(Point(-13.5, -13.5)) == cluster5a);
        REQUIRE(rootNode->findLeaf(Point(-12.0, -14.0)) == cluster5b);
        REQUIRE(rootNode->findLeaf(Point(-12.5, -14.5)) == cluster5c);
        REQUIRE(rootNode->findLeaf(Point(-13.0, -15.0)) == cluster5d);

        tree.root = root;

        // Destroy tree
        tree.write_metadata();
    }

    // Read existing tree from disk
    DefaulTreeType tree( 4096 * 5, "nirdiskbacked.txt" );
    auto rootNode = tree.get_branch_node( tree.root );

	// Test finding leaves
	REQUIRE(rootNode->findLeaf(Point(-11.0, -3.0)) == cluster4a);
	REQUIRE(rootNode->findLeaf(Point(-9.0, -3.0)) == cluster4b);
	REQUIRE(rootNode->findLeaf(Point(-13.5, -13.5)) == cluster5a);
	REQUIRE(rootNode->findLeaf(Point(-12.0, -14.0)) == cluster5b);
	REQUIRE(rootNode->findLeaf(Point(-12.5, -14.5)) == cluster5c);
	REQUIRE(rootNode->findLeaf(Point(-13.0, -15.0)) == cluster5d);


    unlink( "nirdiskbacked.txt" );
    unlink( "nirdiskbacked.txt.meta" );
}

TEST_CASE("NIRTreeDisk: testInsertGrowTreeHeight")
{
    unlink( "nirdiskbacked.txt" );
    {
        unsigned maxBranchFactor = 7;
        DefaulTreeType tree(4096*5, "nirdiskbacked.txt");

        for( unsigned i = 0; i < maxBranchFactor + 1; i++) {
            tree.insert( Point(i,i) );
        }


        auto rootNode =
            tree.node_allocator_->get_tree_node<DefaultBranchNodeType>(
                    tree.root );
        REQUIRE(rootNode->cur_offset_ == 2);
        REQUIRE( rootNode->level_ == 1 );

        nirtreedisk::Branch bLeft = rootNode->entries[0];
        nirtreedisk::Branch bRight = rootNode->entries[1];
        REQUIRE( bLeft.child.get_type() ==  LEAF_NODE );
        REQUIRE( bRight.child.get_type() ==  LEAF_NODE );

        auto left = tree.get_leaf_node( bLeft.child );
        auto right = tree.get_leaf_node( bRight.child );


        REQUIRE(left->cur_offset_ == 4);
        REQUIRE(right->cur_offset_ == 4);
    }
    unlink( "nirdiskbacked.txt" );
}

TEST_CASE("NIRTreeDisk: doubleGrowTreeHeight")
{
    unlink( "nirdiskbacked.txt" );
    {
        unsigned max_branch_factor = 7;
        unsigned insertion_count = max_branch_factor*max_branch_factor + 1;

        MeanBalancedTreeType tree(4096*20, "nirdiskbacked.txt");
        for( unsigned i = 0; i < insertion_count; i++) {
            tree.insert(Point(i,i));
        }

        tree_node_handle root = tree.root;
        auto root_node = tree.get_branch_node( root );

        for( unsigned i = 0; i < insertion_count; i++) {
            REQUIRE( tree.search( Point(i,i) ).size() == 1 );
        }

        REQUIRE( root_node->cur_offset_ == 3 );
        REQUIRE( root_node->level_ == 2 );

        nirtreedisk::Branch bLeft = root_node->entries[0];
        nirtreedisk::Branch bRight = root_node->entries[1];

        auto left = tree.get_branch_node( bLeft.child );
        auto right = tree.get_branch_node( bRight.child );

        REQUIRE(left->cur_offset_ == 4);
        REQUIRE(right->cur_offset_ == 4);
        REQUIRE( left->level_ == 1 );
        REQUIRE( right->level_ == 1 );
    }
    unlink( "nirdiskbacked.txt" );

}

TEST_CASE( "NIRTreeDisk: grow well-beyond memory provisions" )
{
    unlink( "nirdiskbacked.txt" );
    {
        unsigned max_branch_factor = 7;
        // Guestimate way more nodes.
        size_t nodes_per_page = PAGE_DATA_SIZE / sizeof(
                DefaultLeafNodeType);

        // We need a decent number of pages in memory because during
        // searches the whole path down to the leaf is pinned.
        size_t page_count = 20;

        size_t insertion_count = max_branch_factor * (nodes_per_page *
                page_count) * 4;

        MeanBalancedTreeType tree(4096*page_count, "nirdiskbacked.txt");
        for( unsigned i = 0; i < insertion_count; i++) {
            tree.insert(Point(i,i));
        }

        for( unsigned i = 0; i < insertion_count; i++) {
            REQUIRE( tree.search( Point(i,i) ).size() == 1);
        }
    }

    unlink( "nirdiskbacked.txt" );
}

TEST_CASE("NIRTreeDisk: pack simple leaf node") {
    unlink( "nirdiskbacked.txt" );
    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );
    for( unsigned i = 0; i < 5; i++ ) {
        tree.insert( Point(i,i) );
    }

    auto root_leaf_node = tree.get_leaf_node( tree.root );
    REQUIRE( root_leaf_node->cur_offset_ == 5 );

    REQUIRE( root_leaf_node->compute_packed_size() <
            sizeof(DefaultLeafNodeType) );
    REQUIRE( root_leaf_node->compute_packed_size() ==
            sizeof(unsigned) + sizeof(Point) * 5 );

    tree_node_handle repacked_handle = root_leaf_node->repack(
            tree.node_allocator_.get() );

    REQUIRE( repacked_handle != nullptr );
    auto packed_leaf =
        tree.node_allocator_->get_tree_node<packed_node>(
                repacked_handle );

    REQUIRE( * (unsigned *) (packed_leaf->buffer_) ==
        root_leaf_node->cur_offset_ );

    Point *p = (Point *) (packed_leaf->buffer_ + sizeof(unsigned));
    REQUIRE( *(p++) == root_leaf_node->entries.at(0) );
    REQUIRE( *(p++) == root_leaf_node->entries.at(1) );
    REQUIRE( *(p++) == root_leaf_node->entries.at(2) );
    REQUIRE( *(p++) == root_leaf_node->entries.at(3) );
    REQUIRE( *(p++) == root_leaf_node->entries.at(4) );
    unlink( "nirdiskbacked.txt" );
}

TEST_CASE("NIRTreeDisk: pack branch node all inline") {
    unlink( "nirdiskbacked.txt" );
    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );
    for( unsigned i = 0; i < 5; i++ ) {
        tree.insert( Point(i,i) );
    }

    auto alloc_branch_data =
        tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>();
    new (&(*alloc_branch_data.first)) DefaultBranchNodeType( &tree,
            tree_node_handle(nullptr), alloc_branch_data.second, 1 );

    auto alloc_leaf_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>();
    new (&(*alloc_leaf_data.first)) DefaultLeafNodeType( &tree,
            alloc_branch_data.second, alloc_leaf_data.second, 0 );

    auto branch_node = alloc_branch_data.first;
    auto leaf_node = alloc_leaf_data.first;
    auto leaf_handle = alloc_leaf_data.second;

    // Add five different branches for the same leaf.
    // This isn't valid in practice but for testing purposes it is fine.
    nirtreedisk::Branch b = createBranchEntry(
            InlineBoundedIsotheticPolygon(
                Rectangle(0.0, 0.0, 1.0, 1.0)), leaf_handle );
    branch_node->addBranchToNode( b );
    b = createBranchEntry( InlineBoundedIsotheticPolygon(
                Rectangle(0.0, 0.0, 2.0, 2.0)), leaf_handle );
    branch_node->addBranchToNode( b );
    b = createBranchEntry( InlineBoundedIsotheticPolygon(
                Rectangle(0.0, 0.0, 3.0, 3.0)), leaf_handle );
    branch_node->addBranchToNode( b );
    b = createBranchEntry( InlineBoundedIsotheticPolygon(
                Rectangle(0.0, 0.0, 4.0, 4.0)), leaf_handle );
    branch_node->addBranchToNode( b );
    b = createBranchEntry( InlineBoundedIsotheticPolygon(
                Rectangle(0.0, 0.0, 5.0, 5.0)), leaf_handle );
    branch_node->addBranchToNode( b );

    REQUIRE( branch_node->cur_offset_ == 5 );

    tree_node_handle packed_handle = branch_node->repack(
            tree.node_allocator_.get(), tree.node_allocator_.get() );
    REQUIRE( packed_handle != nullptr );

    auto packed_branch= tree.node_allocator_->get_tree_node<packed_node>(
            packed_handle );
    
    size_t offset = 0;
    REQUIRE( * (unsigned *) (packed_branch->buffer_ + offset) ==
            branch_node->cur_offset_ );
    offset += sizeof(unsigned);

    // There are 5 entries

    //Entry 1
    REQUIRE( * (tree_node_handle *) (packed_branch->buffer_ + offset) ==
            leaf_handle );
    offset += sizeof(tree_node_handle);
    REQUIRE( * (unsigned *) (packed_branch->buffer_ + offset) ==
            1U );
    offset += sizeof(unsigned);
    REQUIRE( * (Rectangle *) (packed_branch->buffer_ + offset) ==
            Rectangle(0.0,0.0,1.0,1.0) );
    offset += sizeof(Rectangle);

    //Entry 2
    REQUIRE( * (tree_node_handle *) (packed_branch->buffer_ + offset) ==
            leaf_handle );
    offset += sizeof(tree_node_handle);
    REQUIRE( * (unsigned *) (packed_branch->buffer_ + offset) ==
            1U );
    offset += sizeof(unsigned);
    REQUIRE( * (Rectangle *) (packed_branch->buffer_ + offset) ==
            Rectangle(0.0,0.0,2.0,2.0) );
    offset += sizeof(Rectangle);

    //Entry 3
    REQUIRE( * (tree_node_handle *) (packed_branch->buffer_ + offset) ==
            leaf_handle );
    offset += sizeof(tree_node_handle);
    REQUIRE( * (unsigned *) (packed_branch->buffer_ + offset) ==
            1U );
    offset += sizeof(unsigned);
    REQUIRE( * (Rectangle *) (packed_branch->buffer_ + offset) ==
            Rectangle(0.0,0.0,3.0,3.0) );
    offset += sizeof(Rectangle);

    //Entry 4
    REQUIRE( * (tree_node_handle *) (packed_branch->buffer_ + offset) ==
            leaf_handle );
    offset += sizeof(tree_node_handle);
    REQUIRE( * (unsigned *) (packed_branch->buffer_ + offset) ==
            1U );
    offset += sizeof(unsigned);
    REQUIRE( * (Rectangle *) (packed_branch->buffer_ + offset) ==
            Rectangle(0.0,0.0,4.0,4.0) );
    offset += sizeof(Rectangle);

    //Entry 5
    REQUIRE( * (tree_node_handle *) (packed_branch->buffer_ + offset) ==
            leaf_handle );
    offset += sizeof(tree_node_handle);
    REQUIRE( * (unsigned *) (packed_branch->buffer_ + offset) ==
            1U );
    offset += sizeof(unsigned);
    REQUIRE( * (Rectangle *) (packed_branch->buffer_ + offset) ==
            Rectangle(0.0,0.0,5.0,5.0) );
    offset += sizeof(Rectangle);

    unlink( "nirdiskbacked.txt" );
}

TEST_CASE("NIRTreeDisk: pack complex inline polygon") {
    unlink( "nirdiskbacked.txt" );
    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );
    for( unsigned i = 0; i < 5; i++ ) {
        tree.insert( Point(i,i) );
    }

    auto alloc_branch_data =
        tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>();
    new (&(*alloc_branch_data.first)) DefaultBranchNodeType( &tree,
            tree_node_handle(nullptr), alloc_branch_data.second, 1 );

    auto alloc_leaf_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>();
    new (&(*alloc_leaf_data.first)) DefaultLeafNodeType( &tree,
            alloc_branch_data.second, alloc_leaf_data.second, 0 );

    auto branch_node = alloc_branch_data.first;
    auto leaf_node = alloc_leaf_data.first;
    auto leaf_handle = alloc_leaf_data.second;

    IsotheticPolygon polygon( Rectangle(0.0,0.0,1.0,1.0) );
    polygon.basicRectangles.push_back( Rectangle( 1.0, 2.0, 3.0, 4.0 ) );
    polygon.basicRectangles.push_back( Rectangle( -1.0, -2.0, -3.0, -4.0 ) );
    polygon.basicRectangles.push_back( Rectangle( 10.0, 10.0, 30.0, 40.0 ) );
    polygon.basicRectangles.push_back( Rectangle( 1.0, 3.0, 3.0, 7.0 ) );
    polygon.recomputeBoundingBox();

    nirtreedisk::Branch b = createBranchEntry(
            InlineBoundedIsotheticPolygon(), leaf_handle );
    std::get<InlineBoundedIsotheticPolygon>(b.boundingPoly).push_polygon_to_disk( polygon );
    branch_node->addBranchToNode( b );

    tree_node_handle packed_handle = branch_node->repack(
            tree.node_allocator_.get(), tree.node_allocator_.get() ); 

    auto packed_branch= tree.node_allocator_->get_tree_node<packed_node>(
            packed_handle );
    
    size_t offset = 0;
    REQUIRE( * (unsigned*) (packed_branch->buffer_ + offset) ==
            branch_node->cur_offset_ );
    offset += sizeof(unsigned);

    //Entry 1
    REQUIRE( * (tree_node_handle *) (packed_branch->buffer_ + offset) ==
            leaf_handle );
    offset += sizeof(tree_node_handle);
    REQUIRE( * (unsigned *) (packed_branch->buffer_ + offset) ==
            5U );
    offset += sizeof(unsigned);
    REQUIRE( * (Rectangle *) (packed_branch->buffer_ + offset) ==
            Rectangle(0.0,0.0,1.0,1.0) );
    offset += sizeof(Rectangle);
    REQUIRE( * (Rectangle *) (packed_branch->buffer_ + offset) ==
            Rectangle(1.0,2.0,3.0,4.0) );
    offset += sizeof(Rectangle);
    REQUIRE( * (Rectangle *) (packed_branch->buffer_ + offset) ==
            Rectangle(-1.0,-2.0,-3.0,-4.0) );
    offset += sizeof(Rectangle);
    REQUIRE( * (Rectangle *) (packed_branch->buffer_ + offset) ==
            Rectangle(10.0,10.0,30.0,40.0) );
    offset += sizeof(Rectangle);
    REQUIRE( * (Rectangle *) (packed_branch->buffer_ + offset) ==
            Rectangle(1.0,3.0,3.0,7.0) );
    offset += sizeof(Rectangle);
}

TEST_CASE("NIRTreeDisk: pack out of line polygon") {
    unlink( "nirdiskbacked.txt" );
    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );
    for( unsigned i = 0; i < 5; i++ ) {
        tree.insert( Point(i,i) );
    }

    auto alloc_branch_data =
        tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>();
    new (&(*alloc_branch_data.first)) DefaultBranchNodeType( &tree,
            tree_node_handle(nullptr), alloc_branch_data.second, 1 );

    auto alloc_leaf_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>();
    new (&(*alloc_leaf_data.first)) DefaultLeafNodeType( &tree,
            alloc_branch_data.second, alloc_leaf_data.second, 0 );

    auto branch_node = alloc_branch_data.first;
    auto leaf_node = alloc_leaf_data.first;
    auto leaf_handle = alloc_leaf_data.second;

    IsotheticPolygon polygon( Rectangle(0.0, 0.0, 1.0, 1.0) );
    for( double i = 1.0; i < 30.0; i += 1.0 ) {
        polygon.basicRectangles.push_back( Rectangle( i, i, i+1.0,
                    i+1.0) );
    }
    polygon.recomputeBoundingBox();

    auto alloc_poly_data = tree.node_allocator_->create_new_tree_node<InlineUnboundedIsotheticPolygon>(
            compute_sizeof_inline_unbounded_polygon(
                polygon.basicRectangles.size()
                 ), NodeHandleType(BIG_POLYGON));
    new (&(*alloc_poly_data.first)) InlineUnboundedIsotheticPolygon(
            tree.node_allocator_.get(), polygon.basicRectangles.size() );
    alloc_poly_data.first->push_polygon_to_disk( polygon );

    nirtreedisk::Branch b = createBranchEntry( alloc_poly_data.second, leaf_handle );
    branch_node->addBranchToNode( b );

    tree_node_handle packed_handle = branch_node->repack(
            tree.node_allocator_.get(), tree.node_allocator_.get() ); 

    auto packed_branch = tree.node_allocator_->get_tree_node<packed_node>(
            packed_handle );

    size_t offset = 0;
    REQUIRE( * (unsigned*) (packed_branch->buffer_ + offset) ==
            branch_node->cur_offset_ );
    offset += sizeof(unsigned);

    // Entry 1
    tree_node_handle child_handle = * (tree_node_handle *)
        (packed_branch->buffer_ + offset);
    REQUIRE( child_handle == leaf_handle ); 
    offset += sizeof( tree_node_handle );
    bool is_compressed =
        child_handle.get_associated_poly_is_compressed();
    REQUIRE( is_compressed == true );
    IsotheticPolygon decoded_poly =
        decompress_polygon( (packed_branch->buffer_ + offset) );
    REQUIRE( decoded_poly.basicRectangles.size() == 30 );
    for( unsigned i = 0; i < decoded_poly.basicRectangles.size(); i++ ) {
        REQUIRE( decoded_poly.basicRectangles.at(i) == Rectangle( i, i,
                    i+1, i+1) );
    }
    unlink("nirdiskbacked.txt");
}

TEST_CASE("NIRTreeDisk: pack in a small out of band polygon") {
    unlink( "nirdiskbacked.txt" );
    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );
    for( unsigned i = 0; i < 5; i++ ) {
        tree.insert( Point(i,i) );
    }

    auto alloc_branch_data =
        tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>();
    new (&(*alloc_branch_data.first)) DefaultBranchNodeType( &tree,
            tree_node_handle(nullptr), alloc_branch_data.second, 1 );

    auto alloc_leaf_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>();
    new (&(*alloc_leaf_data.first)) DefaultLeafNodeType( &tree,
            alloc_branch_data.second, alloc_leaf_data.second, 0 );

    auto branch_node = alloc_branch_data.first;
    auto leaf_node = alloc_leaf_data.first;
    auto leaf_handle = alloc_leaf_data.second;

    // This is small enough to fit
    IsotheticPolygon polygon( Rectangle(0.0, 0.0, 1.0, 1.0) );
    for( double i = 1.0; i < 20.0; i += 1.0 ) {
        polygon.basicRectangles.push_back( Rectangle( i, i, i+1.0,
                    i+1.0) );
    }
    polygon.recomputeBoundingBox();

    auto alloc_poly_data = tree.node_allocator_->create_new_tree_node<InlineUnboundedIsotheticPolygon>(
            compute_sizeof_inline_unbounded_polygon(
                polygon.basicRectangles.size()
                 ), NodeHandleType(BIG_POLYGON));
    new (&(*alloc_poly_data.first)) InlineUnboundedIsotheticPolygon(
            tree.node_allocator_.get(), polygon.basicRectangles.size() );
    alloc_poly_data.first->push_polygon_to_disk( polygon );

    nirtreedisk::Branch b = createBranchEntry( alloc_poly_data.second, leaf_handle );
    branch_node->addBranchToNode( b );

    // too big, should be out of line
    IsotheticPolygon polygon2( Rectangle(0.0, 0.0, 1.0, 1.0) );
    for( double i = 1.0; i < 30.0; i += 1.0 ) {
        polygon2.basicRectangles.push_back( Rectangle( i, i, i+1.0,
                    i+1.0) );
    }
    polygon2.recomputeBoundingBox();

    alloc_poly_data = tree.node_allocator_->create_new_tree_node<InlineUnboundedIsotheticPolygon>(
            compute_sizeof_inline_unbounded_polygon(
                polygon2.basicRectangles.size()
                 ), NodeHandleType(BIG_POLYGON));
    new (&(*alloc_poly_data.first)) InlineUnboundedIsotheticPolygon(
            tree.node_allocator_.get(), polygon2.basicRectangles.size() );
    alloc_poly_data.first->push_polygon_to_disk( polygon2 );

    b = createBranchEntry( alloc_poly_data.second, leaf_handle );
    branch_node->addBranchToNode( b );

    tree_node_handle packed_handle = branch_node->repack(
            tree.node_allocator_.get(), tree.node_allocator_.get() ); 

    auto packed_branch= tree.node_allocator_->get_tree_node<packed_node>(
            packed_handle );

    size_t offset = 0;
    REQUIRE( * (unsigned*) (packed_branch->buffer_ + offset) ==
            branch_node->cur_offset_ );
    offset += sizeof(unsigned);

    // Entry 1
    int new_offset = 0;
    tree_node_handle child_handle = * (tree_node_handle *)
        (packed_branch->buffer_ + offset);
    REQUIRE( child_handle == leaf_handle ); 
    REQUIRE( child_handle.get_associated_poly_is_compressed() == true );
    offset += sizeof( tree_node_handle );

    IsotheticPolygon child_poly = decompress_polygon(
            packed_branch->buffer_ + offset, &new_offset );
    REQUIRE( child_poly.basicRectangles.size() == 20 );
    offset += new_offset;

    child_handle = * (tree_node_handle *)
        (packed_branch->buffer_ + offset);
    REQUIRE( child_handle == leaf_handle ); 
    REQUIRE( child_handle.get_associated_poly_is_compressed() == true );
    offset += sizeof( tree_node_handle );

    child_poly = decompress_polygon( packed_branch->buffer_ + offset );
    REQUIRE( child_poly.basicRectangles.size() == 30 );
    unlink("nirdiskbacked.txt");
}



TEST_CASE("NIRTreeDisk: Search packed leaf from branch." ) {

    unlink( "nirdiskbacked.txt" );
    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );
    for( unsigned i = 0; i < 5; i++ ) {
        tree.insert( Point(i,i) );
    }

    auto root_leaf_node = tree.get_leaf_node( tree.root );
    REQUIRE( root_leaf_node->cur_offset_ == 5 );

    REQUIRE( root_leaf_node->compute_packed_size() <
            sizeof(DefaultLeafNodeType) );
    REQUIRE( root_leaf_node->compute_packed_size() ==
            sizeof(unsigned) + sizeof(Point) * 5 );

    tree_node_handle repacked_handle = root_leaf_node->repack(
            tree.node_allocator_.get() );
    repacked_handle.set_type(
            NodeHandleType(REPACKED_LEAF_NODE) );

    REQUIRE( repacked_handle != nullptr );
    auto packed_leaf =
        tree.node_allocator_->get_tree_node<packed_node>(
                repacked_handle );

    auto alloc_branch_data =
        tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>(
                NodeHandleType(BRANCH_NODE));
    new (&(*alloc_branch_data.first)) DefaultBranchNodeType( &tree,
            tree_node_handle(nullptr), alloc_branch_data.second, 1 );

    auto branch_node = alloc_branch_data.first;

    nirtreedisk::Branch b = createBranchEntry(
            InlineBoundedIsotheticPolygon(root_leaf_node->boundingBox()),
            tree.root );
    branch_node->addBranchToNode( b );

    for( int i = 0; i < 7; i++ ) {
        Point p( i, i );
        auto vec = point_search( branch_node->self_handle_, p, &tree );
        if( i < 5 ) {
            REQUIRE( vec.size() == 1 );
        } else {
            REQUIRE( vec.size() == 0 );
        }
    }
    unlink( "nirdiskbacked.txt");
}

TEST_CASE("NIRTreeDisk: Search packed leaf direct." ) {

    unlink( "nirdiskbacked.txt" );
    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );
    for( unsigned i = 0; i < 5; i++ ) {
        tree.insert( Point(i,i) );
    }

    auto root_leaf_node = tree.get_leaf_node( tree.root );
    REQUIRE( root_leaf_node->cur_offset_ == 5 );

    REQUIRE( root_leaf_node->compute_packed_size() <
            sizeof(DefaultLeafNodeType) );
    REQUIRE( root_leaf_node->compute_packed_size() ==
            sizeof(unsigned) + sizeof(Point) * 5 );

    tree_node_handle repacked_handle = root_leaf_node->repack(
            tree.node_allocator_.get() );

    REQUIRE( repacked_handle != nullptr );
    auto packed_leaf =
        tree.node_allocator_->get_tree_node<packed_node>(
                repacked_handle );

    for( int i = 0; i < 7; i++ ) {
        Point p( i, i );
        auto vec = point_search(repacked_handle, p, &tree );
        if( i < 5 ) {
            REQUIRE( vec.size() == 1 );
        } else {
            REQUIRE( vec.size() == 0 );
        }
    }
    unlink( "nirdiskbacked.txt");
}

TEST_CASE( "NIRTreeDisk: Point search packed branch" ) {
    unlink( "nirdiskbacked.txt" );

    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );
    auto alloc_branch_data = tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>(
            NodeHandleType(BRANCH_NODE));
    auto branch_handle = alloc_branch_data.second;
    new (&(*alloc_branch_data.first)) DefaultBranchNodeType( &tree,
            tree_node_handle(nullptr), alloc_branch_data.second, 1 );
    auto branch_node = alloc_branch_data.first;

    auto alloc_leaf_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
            NodeHandleType(LEAF_NODE));
    new (&(*alloc_leaf_data.first)) DefaultLeafNodeType( &tree,
            branch_handle, alloc_leaf_data.second, 0 );
    auto leaf_node1 = alloc_leaf_data.first;

    alloc_leaf_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
            NodeHandleType(LEAF_NODE));
    new (&(*alloc_leaf_data.first)) DefaultLeafNodeType( &tree,
            branch_handle, alloc_leaf_data.second, 0 );
    auto leaf_node2 = alloc_leaf_data.first;

    alloc_leaf_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
            NodeHandleType(LEAF_NODE));
    new (&(*alloc_leaf_data.first)) DefaultLeafNodeType( &tree,
            branch_handle, alloc_leaf_data.second, 0 );
    auto leaf_node3 = alloc_leaf_data.first;

    leaf_node1->addPoint( Point(1,1) );
    leaf_node1->addPoint( Point(2,2) );
    leaf_node1->addPoint( Point(3,3) );
    leaf_node1->addPoint( Point(4,4) );
    leaf_node1->addPoint( Point(5,5) );
    nirtreedisk::Branch b = createBranchEntry(
            InlineBoundedIsotheticPolygon(leaf_node1->boundingBox()),
            leaf_node1->self_handle_ );
    branch_node->addBranchToNode( b );

    leaf_node2->addPoint( Point(6,6) );
    leaf_node2->addPoint( Point(7,7) );
    leaf_node2->addPoint( Point(8,8) );
    leaf_node2->addPoint( Point(9,9) );
    leaf_node2->addPoint( Point(10,10) );
    b = createBranchEntry(
            InlineBoundedIsotheticPolygon(leaf_node2->boundingBox()),
            leaf_node2->self_handle_ );
    branch_node->addBranchToNode( b );

    leaf_node3->addPoint( Point(11,11) );
    leaf_node3->addPoint( Point(12,12) );
    leaf_node3->addPoint( Point(13,13) );
    leaf_node3->addPoint( Point(14,14) );
    leaf_node3->addPoint( Point(15,15) );
    b = createBranchEntry(
            InlineBoundedIsotheticPolygon(leaf_node3->boundingBox()),
            leaf_node3->self_handle_ );
    branch_node->addBranchToNode( b );

    auto packed_branch_handle = branch_node->repack(
            tree.node_allocator_.get(), tree.node_allocator_.get() );

    for( unsigned i = 1; i < 20; i++ ) {
        Point p(i,i);
        if( i < 16 ) {
            REQUIRE( point_search( packed_branch_handle, p, &tree ).size() == 1 );
        } else {
            REQUIRE( point_search( packed_branch_handle, p, &tree ).size() == 0 );
        }
    }

    unlink( "nirdiskbacked.txt" );
}

TEST_CASE( "NIRTreeDisk: Point search packed branch multi-rect poly" ) {

    unlink( "nirdiskbacked.txt" );

    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );
    auto alloc_branch_data = tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>(
            NodeHandleType(BRANCH_NODE));
    auto branch_handle = alloc_branch_data.second;
    new (&(*alloc_branch_data.first)) DefaultBranchNodeType( &tree,
            tree_node_handle(nullptr), alloc_branch_data.second, 1 );
    auto branch_node = alloc_branch_data.first;

    auto alloc_leaf_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
            NodeHandleType(LEAF_NODE));
    new (&(*alloc_leaf_data.first)) DefaultLeafNodeType( &tree,
            branch_handle, alloc_leaf_data.second, 0 );
    auto leaf_node1 = alloc_leaf_data.first;

    alloc_leaf_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
            NodeHandleType(LEAF_NODE));
    new (&(*alloc_leaf_data.first)) DefaultLeafNodeType( &tree,
            branch_handle, alloc_leaf_data.second, 0 );
    auto leaf_node2 = alloc_leaf_data.first;

    alloc_leaf_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
            NodeHandleType(LEAF_NODE));
    new (&(*alloc_leaf_data.first)) DefaultLeafNodeType( &tree,
            branch_handle, alloc_leaf_data.second, 0 );
    auto leaf_node3 = alloc_leaf_data.first;

    leaf_node1->addPoint( Point(1,1) );
    leaf_node1->addPoint( Point(2,2) );
    leaf_node1->addPoint( Point(3,3) );
    leaf_node1->addPoint( Point(4,4) );
    leaf_node1->addPoint( Point(5,5) );

    Rectangle rect1(1,1,nextafter(1,DBL_MAX),nextafter(5, DBL_MAX));
    Rectangle rect2(7,1,nextafter(20, DBL_MAX),nextafter(1,DBL_MAX));
    Rectangle rect3(-100,-100, -20, -20);
    Rectangle rect4(1,1,nextafter(5,DBL_MAX),nextafter(5,DBL_MAX));
    IsotheticPolygon polygon( rect1 );
    polygon.basicRectangles.push_back( rect2 );
    polygon.basicRectangles.push_back( rect3 );
    polygon.basicRectangles.push_back( rect4 );
    polygon.recomputeBoundingBox();

    nirtreedisk::Branch b = createBranchEntry(
            InlineBoundedIsotheticPolygon(),
            leaf_node1->self_handle_ );
    std::get<InlineBoundedIsotheticPolygon>( b.boundingPoly
            ).push_polygon_to_disk( polygon );
    branch_node->addBranchToNode( b );

    leaf_node2->addPoint( Point(6,6) );
    leaf_node2->addPoint( Point(7,7) );
    leaf_node2->addPoint( Point(8,8) );
    leaf_node2->addPoint( Point(9,9) );
    leaf_node2->addPoint( Point(10,10) );
    b = createBranchEntry(
            InlineBoundedIsotheticPolygon(leaf_node2->boundingBox()),
            leaf_node2->self_handle_ );
    branch_node->addBranchToNode( b );

    leaf_node3->addPoint( Point(11,11) );
    leaf_node3->addPoint( Point(12,12) );
    leaf_node3->addPoint( Point(13,13) );
    leaf_node3->addPoint( Point(14,14) );
    leaf_node3->addPoint( Point(15,15) );
    b = createBranchEntry(
            InlineBoundedIsotheticPolygon(leaf_node3->boundingBox()),
            leaf_node3->self_handle_ );
    branch_node->addBranchToNode( b );

    auto packed_branch_handle = branch_node->repack(
            tree.node_allocator_.get(), tree.node_allocator_.get() );

    for( unsigned i = 1; i < 20; i++ ) {
        Point p(i,i);
        if( i < 16 ) {
            REQUIRE( point_search( packed_branch_handle, p, &tree ).size() == 1 );
        } else {
            REQUIRE( point_search( packed_branch_handle, p, &tree ).size() == 0 );
        }
    }

    unlink( "nirdiskbacked.txt" );
}


TEST_CASE( "NIRTreeDisk: Point search packed branch out of band poly" ) {

    unlink( "nirdiskbacked.txt" );

    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );
    auto alloc_branch_data = tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>(
            NodeHandleType(BRANCH_NODE));
    auto branch_handle = alloc_branch_data.second;
    new (&(*alloc_branch_data.first)) DefaultBranchNodeType( &tree,
            tree_node_handle(nullptr), alloc_branch_data.second, 1 );
    auto branch_node = alloc_branch_data.first;

    auto alloc_leaf_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
            NodeHandleType(LEAF_NODE));
    new (&(*alloc_leaf_data.first)) DefaultLeafNodeType( &tree,
            branch_handle, alloc_leaf_data.second, 0 );
    auto leaf_node1 = alloc_leaf_data.first;

    alloc_leaf_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
            NodeHandleType(LEAF_NODE));
    new (&(*alloc_leaf_data.first)) DefaultLeafNodeType( &tree,
            branch_handle, alloc_leaf_data.second, 0 );
    auto leaf_node2 = alloc_leaf_data.first;

    alloc_leaf_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
            NodeHandleType(LEAF_NODE));
    new (&(*alloc_leaf_data.first)) DefaultLeafNodeType( &tree,
            branch_handle, alloc_leaf_data.second, 0 );
    auto leaf_node3 = alloc_leaf_data.first;

    leaf_node1->addPoint( Point(1,1) );
    leaf_node1->addPoint( Point(2,2) );
    leaf_node1->addPoint( Point(3,3) );
    leaf_node1->addPoint( Point(4,4) );
    leaf_node1->addPoint( Point(5,5) );


    IsotheticPolygon polygon( Rectangle(-100,-100,-50,-50) );
    for( int i = -50; i < 1; i++ ) {
        polygon.basicRectangles.push_back( Rectangle(i, i, i+6, i+6) );
    }
    polygon.recomputeBoundingBox();


    auto alloc_poly_data =
        tree.node_allocator_->create_new_tree_node<InlineUnboundedIsotheticPolygon>(
                compute_sizeof_inline_unbounded_polygon(
                    polygon.basicRectangles.size() ),
                NodeHandleType(BIG_POLYGON) );
    new (&(*alloc_poly_data.first)) InlineUnboundedIsotheticPolygon(
            tree.node_allocator_.get(), polygon.basicRectangles.size() );
    alloc_poly_data.first->push_polygon_to_disk( polygon );
    

    nirtreedisk::Branch b = createBranchEntry(
            alloc_poly_data.second,
            leaf_node1->self_handle_ );
    branch_node->addBranchToNode( b );

    leaf_node2->addPoint( Point(6,6) );
    leaf_node2->addPoint( Point(7,7) );
    leaf_node2->addPoint( Point(8,8) );
    leaf_node2->addPoint( Point(9,9) );
    leaf_node2->addPoint( Point(10,10) );
    b = createBranchEntry(
            InlineBoundedIsotheticPolygon(leaf_node2->boundingBox()),
            leaf_node2->self_handle_ );
    branch_node->addBranchToNode( b );

    leaf_node3->addPoint( Point(11,11) );
    leaf_node3->addPoint( Point(12,12) );
    leaf_node3->addPoint( Point(13,13) );
    leaf_node3->addPoint( Point(14,14) );
    leaf_node3->addPoint( Point(15,15) );
    b = createBranchEntry(
            InlineBoundedIsotheticPolygon(leaf_node3->boundingBox()),
            leaf_node3->self_handle_ );
    branch_node->addBranchToNode( b );

    auto packed_branch_handle = branch_node->repack(
            tree.node_allocator_.get(), tree.node_allocator_.get() );

    for( unsigned i = 1; i < 20; i++ ) {
        Point p(i,i);
        if( i < 16 ) {
            REQUIRE( point_search( packed_branch_handle, p, &tree ).size() == 1 );
        } else {
            REQUIRE( point_search( packed_branch_handle, p, &tree ).size() == 0 );
        }
    }

    unlink( "nirdiskbacked.txt" );
}

TEST_CASE( "NIRTreeDisk: Repack subtree" ) {

    unlink( "nirdiskbacked.txt" );

    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );
    auto alloc_branch_data = tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>(
            NodeHandleType(BRANCH_NODE));
    auto branch_handle = alloc_branch_data.second;
    new (&(*alloc_branch_data.first)) DefaultBranchNodeType( &tree,
            tree_node_handle(nullptr), alloc_branch_data.second, 1 );
    auto branch_node = alloc_branch_data.first;

    auto alloc_leaf_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
            NodeHandleType(LEAF_NODE));
    new (&(*alloc_leaf_data.first)) DefaultLeafNodeType( &tree,
            branch_handle, alloc_leaf_data.second, 0 );
    auto leaf_node1 = alloc_leaf_data.first;

    alloc_leaf_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
            NodeHandleType(LEAF_NODE));
    new (&(*alloc_leaf_data.first)) DefaultLeafNodeType( &tree,
            branch_handle, alloc_leaf_data.second, 0 );
    auto leaf_node2 = alloc_leaf_data.first;

    alloc_leaf_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
            NodeHandleType(LEAF_NODE));
    new (&(*alloc_leaf_data.first)) DefaultLeafNodeType( &tree,
            branch_handle, alloc_leaf_data.second, 0 );
    auto leaf_node3 = alloc_leaf_data.first;

    leaf_node1->addPoint( Point(1,1) );
    leaf_node1->addPoint( Point(2,2) );
    leaf_node1->addPoint( Point(3,3) );
    leaf_node1->addPoint( Point(4,4) );
    leaf_node1->addPoint( Point(5,5) );

    Rectangle rect1(1,1,nextafter(1,DBL_MAX),nextafter(5, DBL_MAX));
    Rectangle rect2(7,1,nextafter(20, DBL_MAX),nextafter(1,DBL_MAX));
    Rectangle rect3(-100,-100, -20, -20);
    Rectangle rect4(1,1,nextafter(5,DBL_MAX),nextafter(5,DBL_MAX));
    IsotheticPolygon polygon( rect1 );
    polygon.basicRectangles.push_back( rect2 );
    polygon.basicRectangles.push_back( rect3 );
    polygon.basicRectangles.push_back( rect4 );
    polygon.recomputeBoundingBox();

    nirtreedisk::Branch b = createBranchEntry(
            InlineBoundedIsotheticPolygon(),
            leaf_node1->self_handle_ );
    std::get<InlineBoundedIsotheticPolygon>( b.boundingPoly
            ).push_polygon_to_disk( polygon );
    branch_node->addBranchToNode( b );

    leaf_node2->addPoint( Point(6,6) );
    leaf_node2->addPoint( Point(7,7) );
    leaf_node2->addPoint( Point(8,8) );
    leaf_node2->addPoint( Point(9,9) );
    leaf_node2->addPoint( Point(10,10) );
    b = createBranchEntry(
            InlineBoundedIsotheticPolygon(leaf_node2->boundingBox()),
            leaf_node2->self_handle_ );
    branch_node->addBranchToNode( b );

    leaf_node3->addPoint( Point(11,11) );
    leaf_node3->addPoint( Point(12,12) );
    leaf_node3->addPoint( Point(13,13) );
    leaf_node3->addPoint( Point(14,14) );
    leaf_node3->addPoint( Point(15,15) );
    b = createBranchEntry(
            InlineBoundedIsotheticPolygon(leaf_node3->boundingBox()),
            leaf_node3->self_handle_ );
    branch_node->addBranchToNode( b );

    auto packed_branch_handle =
        nirtreedisk::repack_subtree<3,7,nirtreedisk::LineMinimizeDownsplits>( branch_handle,
            tree.node_allocator_.get(), tree.node_allocator_.get() );

    for( unsigned i = 1; i < 20; i++ ) {
        Point p(i,i);
        if( i < 16 ) {
            REQUIRE( point_search( packed_branch_handle, p, &tree ).size() == 1 );
        } else {
            REQUIRE( point_search( packed_branch_handle, p, &tree ).size() == 0 );
        }
    }

    unlink( "nirdiskbacked.txt" );
}

TEST_CASE( "NIRTreeDisk: Compress Simple Branch" ) {
    unlink( "nirdiskbacked.txt" );

    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );
    auto alloc_branch_data = tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>(
            NodeHandleType(BRANCH_NODE));
    auto branch_handle = alloc_branch_data.second;
    new (&(*alloc_branch_data.first)) DefaultBranchNodeType( &tree,
            tree_node_handle(nullptr), alloc_branch_data.second, 1 );
    auto branch_node = alloc_branch_data.first;

    auto alloc_leaf_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
            NodeHandleType(LEAF_NODE));
    new (&(*alloc_leaf_data.first)) DefaultLeafNodeType( &tree,
            branch_handle, alloc_leaf_data.second, 0 );
    auto leaf_node = alloc_leaf_data.first;

    // Create highly compressible polygon
    IsotheticPolygon branch_polygon;
    branch_polygon.basicRectangles.push_back( Rectangle( 1, 1, 1, 1 ) );
    branch_polygon.basicRectangles.push_back( Rectangle( 1, 1, 1, 1 ) );
    branch_polygon.basicRectangles.push_back( Rectangle( 1, 1, 1, 1 ) );
    branch_polygon.basicRectangles.push_back( Rectangle( 1, 1, 1, 1 ) );

    nirtreedisk::Branch b = createBranchEntry(
            InlineBoundedIsotheticPolygon(),
            leaf_node->self_handle_ );
    std::get<InlineBoundedIsotheticPolygon>( b.boundingPoly
            ).push_polygon_to_disk( branch_polygon );
    branch_node->addBranchToNode( b );

    tree_node_allocator *allocator = tree.node_allocator_.get();

    auto compression_data = b.compute_compression_data( allocator );
    REQUIRE( compression_data.has_value() );
    REQUIRE( compression_data.value().second == 20 );

    IsotheticPolygon decomp_poly = decompress_polygon(
            compression_data.value().first );
    REQUIRE( decomp_poly.basicRectangles.size() ==
            branch_polygon.basicRectangles.size() );
    for( unsigned i = 0; i < decomp_poly.basicRectangles.size(); i++ ) {
        REQUIRE( decomp_poly.basicRectangles.at(i) ==
                branch_polygon.basicRectangles.at(i) );
    }

    unlink( "nirdiskbacked.txt" );

}

TEST_CASE( "NIRTreeDisk: Compress/Repack Single Branch" ) {

    unlink( "nirdiskbacked.txt" );

    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );
    auto alloc_branch_data = tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>(
            NodeHandleType(BRANCH_NODE));
    auto branch_handle = alloc_branch_data.second;
    new (&(*alloc_branch_data.first)) DefaultBranchNodeType( &tree,
            tree_node_handle(nullptr), alloc_branch_data.second, 1 );
    auto branch_node = alloc_branch_data.first;

    auto alloc_leaf_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
            NodeHandleType(LEAF_NODE));
    new (&(*alloc_leaf_data.first)) DefaultLeafNodeType( &tree,
            branch_handle, alloc_leaf_data.second, 0 );
    auto leaf_node = alloc_leaf_data.first;

    // Create highly compressible polygon
    IsotheticPolygon branch_polygon;
    branch_polygon.basicRectangles.push_back( Rectangle( 1, 1, 1, 1 ) );
    branch_polygon.basicRectangles.push_back( Rectangle( 1, 1, 1, 1 ) );
    branch_polygon.basicRectangles.push_back( Rectangle( 1, 1, 1, 1 ) );
    branch_polygon.basicRectangles.push_back( Rectangle( 1, 1, 1, 1 ) );

    nirtreedisk::Branch b = createBranchEntry(
            InlineBoundedIsotheticPolygon(),
            leaf_node->self_handle_ );
    std::get<InlineBoundedIsotheticPolygon>( b.boundingPoly
            ).push_polygon_to_disk( branch_polygon );
    branch_node->addBranchToNode( b );

    tree_node_allocator *allocator = tree.node_allocator_.get();
    tree_node_handle compressed_handle = branch_node->repack( allocator, allocator );

    auto compressed_branch = allocator->get_tree_node<packed_node>(
            compressed_handle );
    char *buffer = compressed_branch->buffer_;
    unsigned offset = 0;
    unsigned count = * (unsigned *) (buffer +offset); \
    offset += sizeof( unsigned );

    REQUIRE( count == 1 );
    tree_node_handle *child = (tree_node_handle *) (buffer + offset );
    REQUIRE( child->get_associated_poly_is_compressed() == true );
    offset += sizeof( tree_node_handle );

    IsotheticPolygon decomp_poly = decompress_polygon( buffer + offset );
    REQUIRE( decomp_poly.basicRectangles.size() ==
            branch_polygon.basicRectangles.size() );
    for( unsigned i = 0; i < decomp_poly.basicRectangles.size(); i++ ) {
        REQUIRE( decomp_poly.basicRectangles.at(i) ==
                branch_polygon.basicRectangles.at(i) );
    }

    unlink( "nirdiskbacked.txt" );
}

TEST_CASE( "NIRTreeDisk: Compress/Repack Multiple-Branch" ) {

    unlink( "nirdiskbacked.txt" );

    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );
    auto alloc_branch_data = tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>(
            NodeHandleType(BRANCH_NODE));
    auto branch_handle = alloc_branch_data.second;
    new (&(*alloc_branch_data.first)) DefaultBranchNodeType( &tree,
            tree_node_handle(nullptr), alloc_branch_data.second, 1 );
    auto branch_node = alloc_branch_data.first;

    auto alloc_leaf_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
            NodeHandleType(LEAF_NODE));
    new (&(*alloc_leaf_data.first)) DefaultLeafNodeType( &tree,
            branch_handle, alloc_leaf_data.second, 0 );
    auto leaf_node = alloc_leaf_data.first;

    // Create highly compressible polygon
    IsotheticPolygon branch_polygon;
    branch_polygon.basicRectangles.push_back( Rectangle( 1, 1, 2, 1 ) );
    branch_polygon.basicRectangles.push_back( Rectangle( 2, 1, 3, 1 ) );
    branch_polygon.basicRectangles.push_back( Rectangle( 3, 1, 4, 1 ) );
    branch_polygon.basicRectangles.push_back( Rectangle( 4, 1, 1, 1 ) );

    nirtreedisk::Branch b = createBranchEntry(
            InlineBoundedIsotheticPolygon(),
            leaf_node->self_handle_ );
    std::get<InlineBoundedIsotheticPolygon>( b.boundingPoly
            ).push_polygon_to_disk( branch_polygon );
    branch_node->addBranchToNode( b );

    IsotheticPolygon branch_polygon2;
    branch_polygon2.basicRectangles.push_back( Rectangle( 5, 1, 1, 2 ) );
    branch_polygon2.basicRectangles.push_back( Rectangle( 5, 2, 1, 3 ) );
    branch_polygon2.basicRectangles.push_back( Rectangle( 5, 3, 1, 4 ) );
    branch_polygon2.basicRectangles.push_back( Rectangle( 5, 4, 1, 1 ) );

    b = createBranchEntry(
            InlineBoundedIsotheticPolygon(),
            leaf_node->self_handle_ );
    std::get<InlineBoundedIsotheticPolygon>( b.boundingPoly
            ).push_polygon_to_disk( branch_polygon2 );
    branch_node->addBranchToNode( b );


    tree_node_allocator *allocator = tree.node_allocator_.get();
    tree_node_handle compressed_handle = branch_node->repack( allocator, allocator );

    auto compressed_branch = allocator->get_tree_node<packed_node>(
            compressed_handle );
    char *buffer = compressed_branch->buffer_;
    unsigned offset = 0;
    unsigned count = * (unsigned *) (buffer +offset); \
    offset += sizeof( unsigned );

    REQUIRE( count == 2 );
    tree_node_handle *child = (tree_node_handle *) (buffer + offset );
    REQUIRE( child->get_associated_poly_is_compressed() == true );
    offset += sizeof( tree_node_handle );

    int new_offset;
    IsotheticPolygon decomp_poly = decompress_polygon( buffer + offset,
            &new_offset );
    REQUIRE( decomp_poly.basicRectangles.size() ==
            branch_polygon.basicRectangles.size() );
    for( unsigned i = 0; i < decomp_poly.basicRectangles.size(); i++ ) {
        REQUIRE( decomp_poly.basicRectangles.at(i) ==
                branch_polygon.basicRectangles.at(i) );
    }
    offset += new_offset;
    child = (tree_node_handle *) (buffer + offset );
    REQUIRE( child->get_associated_poly_is_compressed() == true );
    offset += sizeof( tree_node_handle );

    decomp_poly = decompress_polygon( buffer + offset,
            &new_offset );
    REQUIRE( decomp_poly.basicRectangles.size() ==
            branch_polygon2.basicRectangles.size() );
    for( unsigned i = 0; i < decomp_poly.basicRectangles.size(); i++ ) {
        REQUIRE( decomp_poly.basicRectangles.at(i) ==
                branch_polygon2.basicRectangles.at(i) );
    }

    unlink( "nirdiskbacked.txt" );
}

TEST_CASE( "NIRTreeDisk: Some branches compressed" ) {

    unlink( "nirdiskbacked.txt" );

    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );
    auto alloc_branch_data = tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>(
            NodeHandleType(BRANCH_NODE));
    auto branch_handle = alloc_branch_data.second;
    new (&(*alloc_branch_data.first)) DefaultBranchNodeType( &tree,
            tree_node_handle(nullptr), alloc_branch_data.second, 1 );
    auto branch_node = alloc_branch_data.first;

    auto alloc_leaf_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
            NodeHandleType(LEAF_NODE));
    new (&(*alloc_leaf_data.first)) DefaultLeafNodeType( &tree,
            branch_handle, alloc_leaf_data.second, 0 );
    auto leaf_node = alloc_leaf_data.first;

    // Create highly compressible polygon
    IsotheticPolygon branch_polygon;
    branch_polygon.basicRectangles.push_back( Rectangle( 5, 1, 1, 2 ) );
    branch_polygon.basicRectangles.push_back( Rectangle( 5, 2, 1, 3 ) );
    branch_polygon.basicRectangles.push_back( Rectangle( 5, 3, 1, 4 ) );
    branch_polygon.basicRectangles.push_back( Rectangle( 5, 4, 1, 1 ) );

    nirtreedisk::Branch b = createBranchEntry(
            InlineBoundedIsotheticPolygon(),
            leaf_node->self_handle_ );
    std::get<InlineBoundedIsotheticPolygon>( b.boundingPoly
            ).push_polygon_to_disk( branch_polygon );

    branch_node->addBranchToNode( b );

    IsotheticPolygon branch_polygon2;
    branch_polygon2.basicRectangles.push_back( Rectangle( 1.324234525342,
                -12.34925289, -354.95892761, 12.592089053 ) );
    branch_polygon2.basicRectangles.push_back( Rectangle( 23.8954980323,
                2093.729, 98.6442, 43.942222 ) );
    branch_polygon2.basicRectangles.push_back( Rectangle(
                17.293432908571, -31.2890808942, 1.980809808,
                4.2345982582 ) );

    b = createBranchEntry(
            InlineBoundedIsotheticPolygon(),
            leaf_node->self_handle_ );
    std::get<InlineBoundedIsotheticPolygon>( b.boundingPoly
            ).push_polygon_to_disk( branch_polygon2 );
    branch_node->addBranchToNode( b );

    tree_node_allocator *allocator = tree.node_allocator_.get();
    tree_node_handle compressed_handle = branch_node->repack( allocator, allocator );

    auto compressed_branch = allocator->get_tree_node<packed_node>(
            compressed_handle );
    char *buffer = compressed_branch->buffer_;
    unsigned offset = 0;
    unsigned count = * (unsigned *) (buffer +offset); \
    offset += sizeof( unsigned );

    REQUIRE( count == 2 );
    tree_node_handle *child = (tree_node_handle *) (buffer + offset );
    REQUIRE( child->get_associated_poly_is_compressed() == true );
    offset += sizeof( tree_node_handle );

    int new_offset;
    IsotheticPolygon decomp_poly = decompress_polygon( buffer + offset,
            &new_offset );
    REQUIRE( decomp_poly.basicRectangles.size() ==
            branch_polygon.basicRectangles.size() );
    for( unsigned i = 0; i < decomp_poly.basicRectangles.size(); i++ ) {
        REQUIRE( decomp_poly.basicRectangles.at(i) ==
                branch_polygon.basicRectangles.at(i) );
    }

    offset += new_offset;
    child = (tree_node_handle *) (buffer + offset );
    REQUIRE( child->get_associated_poly_is_compressed() == false );
    offset += sizeof( tree_node_handle );

    uint8_t rect_count = * (uint8_t *) (buffer+offset);
    REQUIRE( rect_count == branch_polygon2.basicRectangles.size() );
    offset += sizeof(uint8_t);

    for (uint8_t i = 0; i < rect_count; i++) {
        Rectangle *rect = (Rectangle *) (buffer + offset);
        REQUIRE( *rect == branch_polygon2.basicRectangles.at(i));
        offset += sizeof( Rectangle );
    }

    unlink( "nirdiskbacked.txt" );
}

TEST_CASE( "NIRTreeDisk: Compress/Repack weird alignment" ) {

    unlink( "nirdiskbacked.txt" );

    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );
    auto alloc_branch_data = tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>(
            NodeHandleType(BRANCH_NODE));
    auto branch_handle = alloc_branch_data.second;
    new (&(*alloc_branch_data.first)) DefaultBranchNodeType( &tree,
            tree_node_handle(nullptr), alloc_branch_data.second, 1 );
    auto branch_node = alloc_branch_data.first;

    auto alloc_leaf_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
            NodeHandleType(LEAF_NODE));
    new (&(*alloc_leaf_data.first)) DefaultLeafNodeType( &tree,
            branch_handle, alloc_leaf_data.second, 0 );
    auto leaf_node = alloc_leaf_data.first;

    // Create highly compressible polygon
    IsotheticPolygon branch_polygon;
    branch_polygon.basicRectangles.push_back( Rectangle( 1, 1, 1, 1 ) );
    branch_polygon.basicRectangles.push_back( Rectangle( 1, 1, 1, 1 ) );
    branch_polygon.basicRectangles.push_back( Rectangle( 1, 1, 1, 1 ) );

    nirtreedisk::Branch b = createBranchEntry(
            InlineBoundedIsotheticPolygon(),
            leaf_node->self_handle_ );
    std::get<InlineBoundedIsotheticPolygon>( b.boundingPoly
            ).push_polygon_to_disk( branch_polygon );
    branch_node->addBranchToNode( b );

    IsotheticPolygon branch_polygon2;
    branch_polygon2.basicRectangles.push_back( Rectangle( 2, 2, 2, 2 ) );
    branch_polygon2.basicRectangles.push_back( Rectangle( 2, 2, 2, 2 ) );
    branch_polygon2.basicRectangles.push_back( Rectangle( 2, 2, 2, 2 ) );

    b = createBranchEntry(
            InlineBoundedIsotheticPolygon(),
            leaf_node->self_handle_ );
    std::get<InlineBoundedIsotheticPolygon>( b.boundingPoly
            ).push_polygon_to_disk( branch_polygon2 );
    branch_node->addBranchToNode( b );

    tree_node_allocator *allocator = tree.node_allocator_.get();
    tree_node_handle compressed_handle = branch_node->repack( allocator, allocator );

    auto compressed_branch = allocator->get_tree_node<packed_node>(
            compressed_handle );
    char *buffer = compressed_branch->buffer_;
    unsigned offset = 0;
    unsigned count = * (unsigned *) (buffer +offset); \
    offset += sizeof( unsigned );

    REQUIRE( count == 2 );
    tree_node_handle *child = (tree_node_handle *) (buffer + offset );
    REQUIRE( child->get_associated_poly_is_compressed() == true );
    offset += sizeof( tree_node_handle );

    int new_offset;
    IsotheticPolygon decomp_poly = decompress_polygon( buffer + offset,
            &new_offset );
    REQUIRE( decomp_poly.basicRectangles.size() ==
            branch_polygon.basicRectangles.size() );
    for( unsigned i = 0; i < decomp_poly.basicRectangles.size(); i++ ) {
        REQUIRE( decomp_poly.basicRectangles.at(i) ==
                branch_polygon.basicRectangles.at(i) );
    }

    offset += new_offset;
    child = (tree_node_handle *) (buffer + offset );
    REQUIRE( child->get_associated_poly_is_compressed() == true );
    offset += sizeof(tree_node_handle);
    decomp_poly = decompress_polygon( buffer + offset,
            &new_offset );
    REQUIRE( decomp_poly.basicRectangles.size() ==
            branch_polygon2.basicRectangles.size() );
    for( unsigned i = 0; i < decomp_poly.basicRectangles.size(); i++ ) {
        REQUIRE( decomp_poly.basicRectangles.at(i) ==
                branch_polygon2.basicRectangles.at(i) );
    }

    unlink( "nirdiskbacked.txt" );
}

TEST_CASE( "NIRTreeDisk: DodgeRectangle NoIntersect" ) {
    unlink( "nirdiskbacked.txt" );

    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );

    tree_node_handle a_handle = nullptr;
    tree_node_handle b_handle = nullptr;
    Rectangle a_rect( 1.0, 1.0, 3.0, 3.0 );
    Rectangle b_rect( 3.1, 3.1, 5.0, 5.0 );
    auto res = nirtreedisk::make_rectangles_disjoint_accounting_for_region_ownership(
        &tree, a_rect, a_handle, b_rect, b_handle );
    REQUIRE( res.first.size() == 1 );
    REQUIRE( res.second.size() == 1 );
    REQUIRE( res.first.at(0) == a_rect );
    REQUIRE( res.second.at(0) == b_rect );

    unlink( "nirdiskbacked.txt" );
}

TEST_CASE( "NIRTreeDisk: DodgeRectangle NoOwnerIntersect A Yields" ) {
    unlink( "nirdiskbacked.txt" );

    // Two rectanglesthat overlap, but there are no points in the
    // overlapping region
    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );

    auto alloc_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                NodeHandleType( LEAF_NODE ) );
    tree_node_handle a_handle = alloc_data.second;
    auto a_leaf_node = alloc_data.first;
    new (&(*a_leaf_node)) DefaultLeafNodeType( &tree, tree_node_handle(nullptr),
            a_handle, 0 );
    a_leaf_node->addPoint( Point(0,1) );
    a_leaf_node->addPoint( Point(2,3) );
    alloc_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                NodeHandleType( LEAF_NODE ) );
    tree_node_handle b_handle = alloc_data.second;
    auto b_leaf_node = alloc_data.first;
    new (&(*b_leaf_node)) DefaultLeafNodeType( &tree, tree_node_handle(nullptr),
            b_handle, 0 );
    b_leaf_node->addPoint( Point(1,0) );
    b_leaf_node->addPoint( Point(3,2) );

    Rectangle a_rect = a_leaf_node->boundingBox();
    Rectangle b_rect = b_leaf_node->boundingBox();
    auto res = nirtreedisk::make_rectangles_disjoint_accounting_for_region_ownership(
        &tree, a_rect, a_handle, b_rect, b_handle );
    REQUIRE( res.first.size() == 2 );
    REQUIRE( res.second.size() == 1 );
    Rectangle decomp_a_first =
        Rectangle(0,1,1,nextafter(3.0, DBL_MAX));
    Rectangle decomp_a_second =
        Rectangle(1,nextafter(2,DBL_MAX),
                nextafter(2, DBL_MAX), nextafter(3.0, DBL_MAX));
    REQUIRE( std::find( res.first.begin(), res.first.end(),
                decomp_a_first ) != res.first.end() );
    REQUIRE( std::find( res.first.begin(), res.first.end(),
                decomp_a_second ) != res.first.end() );
    REQUIRE( res.second.at(0) == b_rect );

    IsotheticPolygon poly;
    poly.basicRectangles = res.first;
    poly.recomputeBoundingBox();
    REQUIRE( not poly.intersectsRectangle( res.second.at(0) ) );

    unlink( "nirdiskbacked.txt" );
}


TEST_CASE("NIRTreeDisk: Unpack simple leaf node") {
    unlink( "nirdiskbacked.txt" );
    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );
    for( unsigned i = 0; i < 5; i++ ) {
        tree.insert( Point(i,i) );
    }

    auto root_leaf_node = tree.get_leaf_node( tree.root );
    REQUIRE( root_leaf_node->cur_offset_ == 5 );

    REQUIRE( root_leaf_node->compute_packed_size() <
            sizeof(DefaultLeafNodeType) );
    REQUIRE( root_leaf_node->compute_packed_size() ==
            sizeof(unsigned) + sizeof(Point) * 5 );

    tree_node_handle repacked_handle = root_leaf_node->repack(
            tree.node_allocator_.get() );

    REQUIRE( repacked_handle != nullptr );
    auto packed_leaf =
        tree.node_allocator_->get_tree_node<packed_node>(
                repacked_handle );
    REQUIRE( * (unsigned *) (packed_leaf->buffer_) ==
        root_leaf_node->cur_offset_ );

    Point *p = (Point *) (packed_leaf->buffer_ + sizeof(unsigned));
    REQUIRE( *(p++) == root_leaf_node->entries.at(0) );
    REQUIRE( *(p++) == root_leaf_node->entries.at(1) );
    REQUIRE( *(p++) == root_leaf_node->entries.at(2) );
    REQUIRE( *(p++) == root_leaf_node->entries.at(3) );
    REQUIRE( *(p++) == root_leaf_node->entries.at(4) );

    tree_node_handle parent_handle( nullptr );
    tree_node_handle unpacked_node_handle =
        nirtreedisk::unpack<DefaultTemplateParams>( repacked_handle,
                tree.node_allocator_.get(), &tree, parent_handle );
    auto unpacked_node = tree.node_allocator_.get()->get_tree_node<DefaultLeafNodeType>( unpacked_node_handle );
    REQUIRE( root_leaf_node->treeRef == &tree );
    REQUIRE( unpacked_node->treeRef == root_leaf_node->treeRef );
    REQUIRE( unpacked_node->parent == root_leaf_node->parent );
    REQUIRE( unpacked_node->cur_offset_ == root_leaf_node->cur_offset_ );
    REQUIRE( unpacked_node->level_ == root_leaf_node->level_ );

    for (unsigned i = 0; i < unpacked_node->cur_offset_; i++) {
        REQUIRE( unpacked_node->entries[i] == root_leaf_node->entries[i] );
    }
    unlink( "nirdiskbacked.txt" );
}

TEST_CASE("NIRTreeDisk: unpack complex inline polygon") {
    unlink( "nirdiskbacked.txt" );
    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );
    for( unsigned i = 0; i < 5; i++ ) {
        tree.insert( Point(i,i) );
    }

    auto alloc_branch_data =
        tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>( NodeHandleType( BRANCH_NODE ) );
    new (&(*alloc_branch_data.first)) DefaultBranchNodeType( &tree,
            tree_node_handle(nullptr), alloc_branch_data.second, 1 );

    auto alloc_leaf_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>( NodeHandleType( LEAF_NODE ) );
    new (&(*alloc_leaf_data.first)) DefaultLeafNodeType( &tree,
            alloc_branch_data.second, alloc_leaf_data.second, 0 );

    auto branch_node = alloc_branch_data.first;
    auto leaf_node = alloc_leaf_data.first;
    auto leaf_handle = alloc_leaf_data.second;

    IsotheticPolygon polygon( Rectangle(0.0,0.0,1.0,1.0) );
    polygon.basicRectangles.push_back( Rectangle( 1.0, 2.0, 3.0, 4.0 ) );
    polygon.basicRectangles.push_back( Rectangle( -1.0, -2.0, -3.0, -4.0 ) );
    polygon.basicRectangles.push_back( Rectangle( 10.0, 10.0, 30.0, 40.0 ) );
    polygon.basicRectangles.push_back( Rectangle( 1.0, 3.0, 3.0, 7.0 ) );
    polygon.recomputeBoundingBox();

    nirtreedisk::Branch b = createBranchEntry(
            InlineBoundedIsotheticPolygon(), leaf_handle );
    std::get<InlineBoundedIsotheticPolygon>(b.boundingPoly).push_polygon_to_disk( polygon );
    branch_node->addBranchToNode( b );

    tree_node_handle packed_handle = branch_node->repack(
            tree.node_allocator_.get(), tree.node_allocator_.get() ); 

    auto packed_branch= tree.node_allocator_->get_tree_node<packed_node>(
            packed_handle );
    
    size_t offset = 0;
    REQUIRE( * (unsigned*) (packed_branch->buffer_ + offset) ==
            branch_node->cur_offset_ );
    offset += sizeof(unsigned);

    //Entry 1
    REQUIRE( * (tree_node_handle *) (packed_branch->buffer_ + offset) ==
            leaf_handle );
    offset += sizeof(tree_node_handle);
    REQUIRE( * (unsigned *) (packed_branch->buffer_ + offset) ==
            5U );
    offset += sizeof(unsigned);
    REQUIRE( * (Rectangle *) (packed_branch->buffer_ + offset) ==
            Rectangle(0.0,0.0,1.0,1.0) );
    offset += sizeof(Rectangle);
    REQUIRE( * (Rectangle *) (packed_branch->buffer_ + offset) ==
            Rectangle(1.0,2.0,3.0,4.0) );
    offset += sizeof(Rectangle);
    REQUIRE( * (Rectangle *) (packed_branch->buffer_ + offset) ==
            Rectangle(-1.0,-2.0,-3.0,-4.0) );
    offset += sizeof(Rectangle);
    REQUIRE( * (Rectangle *) (packed_branch->buffer_ + offset) ==
            Rectangle(10.0,10.0,30.0,40.0) );
    offset += sizeof(Rectangle);
    REQUIRE( * (Rectangle *) (packed_branch->buffer_ + offset) ==
            Rectangle(1.0,3.0,3.0,7.0) );
    offset += sizeof(Rectangle);

    tree_node_handle parent_handle( nullptr );
    tree_node_handle unpacked_node_handle =
        nirtreedisk::unpack<DefaultTemplateParams>( packed_handle,
                tree.node_allocator_.get(), &tree, parent_handle );
    auto unpacked_node = tree.node_allocator_.get()->get_tree_node<DefaultBranchNodeType>( unpacked_node_handle );
    REQUIRE( unpacked_node->treeRef == branch_node->treeRef );
    REQUIRE( unpacked_node->parent == branch_node->parent );
    REQUIRE( unpacked_node->cur_offset_ == branch_node->cur_offset_ );
    REQUIRE( unpacked_node->level_ == branch_node->level_ );

    for (unsigned i = 0; i < unpacked_node->cur_offset_; i++) {
        nirtreedisk::Branch &b1 = unpacked_node->entries[i];
        nirtreedisk::Branch &b2 = branch_node->entries[i];

        REQUIRE( b1.child  ==  b2.child );
        REQUIRE( std::holds_alternative<tree_node_handle>(b1.boundingPoly)
                ==  std::holds_alternative<tree_node_handle>(b2.boundingPoly) );

        if ( std::holds_alternative<tree_node_handle>(b1.boundingPoly) ) {
            REQUIRE( std::get<tree_node_handle>(b1.boundingPoly)
                ==  std::get<tree_node_handle>(b2.boundingPoly) );
        } else {
            REQUIRE( std::get<InlineBoundedIsotheticPolygon>(b1.boundingPoly)
                ==  std::get<InlineBoundedIsotheticPolygon>(b2.boundingPoly) );
        }
    }
    unlink( "nirdiskbacked.txt" );
}

TEST_CASE( "NIRTreeDisk: DodgeRectangle NoOwnerIntersect B Yields" ) {
    unlink( "nirdiskbacked.txt" );

    // Two rectanglesthat overlap, but there are no points in the
    // overlapping region
    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );

    auto alloc_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                NodeHandleType( LEAF_NODE ) );
    tree_node_handle a_handle = alloc_data.second;
    auto a_leaf_node = alloc_data.first;
    new (&(*a_leaf_node)) DefaultLeafNodeType( &tree, tree_node_handle(nullptr),
            a_handle, 0 );
    a_leaf_node->addPoint( Point(0,1) );
    a_leaf_node->addPoint( Point(2,3) );
    a_leaf_node->addPoint( Point(1.5,1.5) );
    alloc_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                NodeHandleType( LEAF_NODE ) );
    tree_node_handle b_handle = alloc_data.second;
    auto b_leaf_node = alloc_data.first;
    new (&(*b_leaf_node)) DefaultLeafNodeType( &tree, tree_node_handle(nullptr),
            b_handle, 0 );
    b_leaf_node->addPoint( Point(1,0) );
    b_leaf_node->addPoint( Point(3,2) );

    Rectangle a_rect = a_leaf_node->boundingBox();
    Rectangle b_rect = b_leaf_node->boundingBox();
    auto res = nirtreedisk::make_rectangles_disjoint_accounting_for_region_ownership(
        &tree, a_rect, a_handle, b_rect, b_handle );
    REQUIRE( res.first.size() == 1 );
    REQUIRE( res.second.size() == 2 );
    REQUIRE( res.first.at(0) == a_rect );

    Rectangle decomp_b_first =
        Rectangle(1,0,nextafter(2,DBL_MAX),1);
    Rectangle decomp_b_second =
        Rectangle(nextafter(2,DBL_MAX),0,nextafter(3,DBL_MAX),nextafter(2,DBL_MAX));

    REQUIRE( std::find( res.second.begin(), res.second.end(),
                decomp_b_first ) != res.second.end() );
    REQUIRE( std::find( res.second.begin(), res.second.end(),
                decomp_b_second ) != res.second.end() );

    IsotheticPolygon poly;
    poly.basicRectangles = res.second;
    poly.recomputeBoundingBox();
    REQUIRE( not poly.intersectsRectangle( res.first.at(0) ) );

    unlink( "nirdiskbacked.txt" );
}

TEST_CASE( "NIRTreeDisk: DodgeRectangle OwnershipIntersect but boxable" ) {
    unlink( "nirdiskbacked.txt" );

    // Two rectanglesthat overlap, but there are no points in the
    // overlapping region
    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );

    auto alloc_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                NodeHandleType( LEAF_NODE ) );
    tree_node_handle a_handle = alloc_data.second;
    auto a_leaf_node = alloc_data.first;
    new (&(*a_leaf_node)) DefaultLeafNodeType( &tree, tree_node_handle(nullptr),
            a_handle, 0 );
    a_leaf_node->addPoint( Point(0,1) );
    a_leaf_node->addPoint( Point(2,3) );
    a_leaf_node->addPoint( Point(1.5,1.5) );
    alloc_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                NodeHandleType( LEAF_NODE ) );
    tree_node_handle b_handle = alloc_data.second;
    auto b_leaf_node = alloc_data.first;
    new (&(*b_leaf_node)) DefaultLeafNodeType( &tree, tree_node_handle(nullptr),
            b_handle, 0 );
    b_leaf_node->addPoint( Point(1,0) );
    b_leaf_node->addPoint( Point(1.7,1.7) );
    b_leaf_node->addPoint( Point(3,2) );

    Rectangle a_rect = a_leaf_node->boundingBox();
    Rectangle b_rect = b_leaf_node->boundingBox();
    auto res = nirtreedisk::make_rectangles_disjoint_accounting_for_region_ownership(
        &tree, a_rect, a_handle, b_rect, b_handle );

    Rectangle decomp_a_first =
        Rectangle(0,1,1,nextafter(3.0, DBL_MAX));
    Rectangle decomp_a_second =
        Rectangle(1,nextafter(2,DBL_MAX),
                nextafter(2, DBL_MAX), nextafter(3.0, DBL_MAX));
    Rectangle decomp_a_third =
        Rectangle(1.5, 1.5, nextafter(1.5,DBL_MAX), nextafter(1.5,
                    DBL_MAX) );
    REQUIRE( res.first.size() == 3 );
    REQUIRE( std::find( res.first.begin(), res.first.end(),
                decomp_a_first ) != res.first.end() );
    REQUIRE( std::find( res.first.begin(), res.first.end(),
                decomp_a_second ) != res.first.end() );
    REQUIRE( std::find( res.first.begin(), res.first.end(),
                decomp_a_third ) != res.first.end() );

    REQUIRE( res.second.size() == 3 );
    Rectangle decomp_b_first =
        Rectangle(1,0,nextafter(2,DBL_MAX),1);
    Rectangle decomp_b_second =
        Rectangle(nextafter(2,DBL_MAX),0,nextafter(3,DBL_MAX),nextafter(2,DBL_MAX));
    Rectangle decomp_b_third=
        Rectangle(1.7, 1.7, nextafter(1.7, DBL_MAX), nextafter(1.7,
                    DBL_MAX));

    REQUIRE( std::find( res.second.begin(), res.second.end(),
                decomp_b_first ) != res.second.end() );
    REQUIRE( std::find( res.second.begin(), res.second.end(),
                decomp_b_second ) != res.second.end() );
    REQUIRE( std::find( res.second.begin(), res.second.end(),
                decomp_b_third ) != res.second.end() );
    IsotheticPolygon poly;
    poly.basicRectangles = res.first;
    poly.recomputeBoundingBox();

    IsotheticPolygon poly2;
    poly2.basicRectangles = res.second;
    poly2.recomputeBoundingBox();

    REQUIRE( poly.disjoint( poly2 ) );


    unlink( "nirdiskbacked.txt" );
}

TEST_CASE( "NIRTreeDisk: DodgeRectangle OwnershipIntersect space slice required" ) {
    unlink( "nirdiskbacked.txt" );

    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );

    auto alloc_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                NodeHandleType( LEAF_NODE ) );
    tree_node_handle a_handle = alloc_data.second;
    auto a_leaf_node = alloc_data.first;
    new (&(*a_leaf_node)) DefaultLeafNodeType( &tree, tree_node_handle(nullptr),
            a_handle, 0 );
    a_leaf_node->addPoint( Point(0,1) );
    a_leaf_node->addPoint( Point(2,3) );
    a_leaf_node->addPoint( Point(1.5,1.5) );
    a_leaf_node->addPoint( Point(1.7,1.7) );
    alloc_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                NodeHandleType( LEAF_NODE ) );
    tree_node_handle b_handle = alloc_data.second;
    auto b_leaf_node = alloc_data.first;
    new (&(*b_leaf_node)) DefaultLeafNodeType( &tree, tree_node_handle(nullptr),
            b_handle, 0 );
    b_leaf_node->addPoint( Point(1,0) );
    b_leaf_node->addPoint( Point(1.6,1.6) );
    b_leaf_node->addPoint( Point(3,2) );

    Rectangle a_rect = a_leaf_node->boundingBox();
    Rectangle b_rect = b_leaf_node->boundingBox();
    auto res = nirtreedisk::make_rectangles_disjoint_accounting_for_region_ownership(
        &tree, a_rect, a_handle, b_rect, b_handle );

    Rectangle decomp_a_first =
        Rectangle(0,1,1,nextafter(3.0, DBL_MAX));
    Rectangle decomp_a_second =
        Rectangle(1,nextafter(2,DBL_MAX),
                nextafter(2, DBL_MAX), nextafter(3.0, DBL_MAX));

    Rectangle decomp_a_third =
        Rectangle(1.5, 1, nextafter(1.5,DBL_MAX), nextafter(2,DBL_MAX) );

    Rectangle decomp_a_fourth =
        Rectangle(1.7, 1, nextafter(1.7,DBL_MAX), nextafter(2,DBL_MAX) );

    REQUIRE( res.first.size() == 4 );
    REQUIRE( std::find( res.first.begin(), res.first.end(),
                decomp_a_first ) != res.first.end() );
    REQUIRE( std::find( res.first.begin(), res.first.end(),
                decomp_a_second ) != res.first.end() );
    REQUIRE( std::find( res.first.begin(), res.first.end(),
                decomp_a_third ) != res.first.end() );
    REQUIRE( std::find( res.first.begin(), res.first.end(),
                decomp_a_fourth ) != res.first.end() );


    REQUIRE( res.second.size() == 3 );
    Rectangle decomp_b_first =
        Rectangle(1,0,nextafter(2,DBL_MAX),1);
    Rectangle decomp_b_second =
        Rectangle(nextafter(2,DBL_MAX),0,nextafter(3,DBL_MAX),nextafter(2,DBL_MAX));
    Rectangle decomp_b_third=
        Rectangle(1.6, 1, nextafter(1.6, DBL_MAX), nextafter(2,
                    DBL_MAX));

    REQUIRE( std::find( res.second.begin(), res.second.end(),
                decomp_b_first ) != res.second.end() );
    REQUIRE( std::find( res.second.begin(), res.second.end(),
                decomp_b_second ) != res.second.end() );
    REQUIRE( std::find( res.second.begin(), res.second.end(),
                decomp_b_third ) != res.second.end() );
    IsotheticPolygon poly;
    poly.basicRectangles = res.first;
    poly.recomputeBoundingBox();

    IsotheticPolygon poly2;
    poly2.basicRectangles = res.second;
    poly2.recomputeBoundingBox();

    REQUIRE( poly.disjoint( poly2 ) );

    unlink( "nirdiskbacked.txt" );
}

TEST_CASE("NIRTreeDisk: unpack out of line polygon") {
    unlink( "nirdiskbacked.txt" );
    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );
    for( unsigned i = 0; i < 5; i++ ) {
        tree.insert( Point(i,i) );
    }

    auto alloc_branch_data =
        tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>( NodeHandleType( BRANCH_NODE ) );
    new (&(*alloc_branch_data.first)) DefaultBranchNodeType( &tree,
            tree_node_handle(nullptr), alloc_branch_data.second, 1 );

    auto alloc_leaf_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>( NodeHandleType( LEAF_NODE ) );
    new (&(*alloc_leaf_data.first)) DefaultLeafNodeType( &tree,
            alloc_branch_data.second, alloc_leaf_data.second, 0 );

    auto branch_node = alloc_branch_data.first;
    auto leaf_node = alloc_leaf_data.first;
    auto leaf_handle = alloc_leaf_data.second;

    IsotheticPolygon polygon( Rectangle(0.0, 0.0, 1.0, 1.0) );
    for( double i = 1.0; i < 30.0; i += 1.0 ) {
        polygon.basicRectangles.push_back( Rectangle( i, i, i+1.0, i+1.0) );
    }
    polygon.recomputeBoundingBox();

    auto alloc_poly_data = tree.node_allocator_->create_new_tree_node<InlineUnboundedIsotheticPolygon>(
            compute_sizeof_inline_unbounded_polygon(
                polygon.basicRectangles.size()
                 ), NodeHandleType(BIG_POLYGON));
    new (&(*alloc_poly_data.first)) InlineUnboundedIsotheticPolygon(
            tree.node_allocator_.get(), polygon.basicRectangles.size() );
    alloc_poly_data.first->push_polygon_to_disk( polygon );

    nirtreedisk::Branch b = createBranchEntry( alloc_poly_data.second, leaf_handle );
    branch_node->addBranchToNode( b );

    tree_node_handle packed_handle = branch_node->repack(
            tree.node_allocator_.get(), tree.node_allocator_.get() ); 

    auto packed_branch= tree.node_allocator_->get_tree_node<packed_node>(
            packed_handle );

    size_t offset = 0;
    REQUIRE( * (unsigned*) (packed_branch->buffer_ + offset) ==
            branch_node->cur_offset_ );
    offset += sizeof(unsigned);

    // Entry 1
    tree_node_handle child_handle = * (tree_node_handle *)
        (packed_branch->buffer_ + offset);
    REQUIRE( child_handle == leaf_handle ); 
    offset += sizeof( tree_node_handle );
    bool is_compressed =
        child_handle.get_associated_poly_is_compressed();
    REQUIRE( is_compressed == true );
    IsotheticPolygon decoded_poly =
        decompress_polygon( (packed_branch->buffer_ + offset) );
    REQUIRE( decoded_poly.basicRectangles.size() == 30 );
    for( unsigned i = 0; i < decoded_poly.basicRectangles.size(); i++ ) {
        REQUIRE( decoded_poly.basicRectangles.at(i) == Rectangle( i, i,
                    i+1, i+1) );
    }

    tree_node_handle parent_handle( nullptr );
    tree_node_handle unpacked_node_handle =
        nirtreedisk::unpack<DefaultTemplateParams>( packed_handle,
                tree.node_allocator_.get(), &tree, parent_handle );
    auto unpacked_node = tree.node_allocator_.get()->get_tree_node<DefaultBranchNodeType>( unpacked_node_handle );
    REQUIRE( unpacked_node->treeRef == branch_node->treeRef );
    REQUIRE( unpacked_node->parent == branch_node->parent );
    REQUIRE( unpacked_node->cur_offset_ == branch_node->cur_offset_ );
    REQUIRE( unpacked_node->level_ == branch_node->level_ );

    for (unsigned i = 0; i < unpacked_node->cur_offset_; i++) {
        nirtreedisk::Branch &b1 = unpacked_node->entries[i];
        nirtreedisk::Branch &b2 = branch_node->entries[i];

        REQUIRE( b1.child  ==  b2.child );
        REQUIRE( std::holds_alternative<tree_node_handle>(b1.boundingPoly)
                ==  std::holds_alternative<tree_node_handle>(b2.boundingPoly) );

        if ( std::holds_alternative<tree_node_handle>(b1.boundingPoly) ) {
            auto poly1 = tree.node_allocator_->get_tree_node<InlineUnboundedIsotheticPolygon>( std::get<tree_node_handle>(b1.boundingPoly) );
            auto poly2 = tree.node_allocator_->get_tree_node<InlineUnboundedIsotheticPolygon>( std::get<tree_node_handle>(b2.boundingPoly) );

            auto it1 = poly1->begin();
            auto it2 = poly2->begin();

            while (it1 != poly1->end() && it2 != poly2->end()) {
                REQUIRE( *it1 == *it2 );
                it1++;
                it2++;
            }

            REQUIRE( it1 == poly1->end() );
            REQUIRE( it2 == poly2->end() );
        } else {
            REQUIRE( std::get<InlineBoundedIsotheticPolygon>(b1.boundingPoly)
                ==  std::get<InlineBoundedIsotheticPolygon>(b2.boundingPoly) );
        }
    }

    unlink("nirdiskbacked.txt");
}

TEST_CASE("NIRTreeDisk: unpack branch node all inline") {
    unlink( "nirdiskbacked.txt" );
    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );
    for( unsigned i = 0; i < 5; i++ ) {
        tree.insert( Point(i,i) );
    }

    auto alloc_branch_data =
        tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>( NodeHandleType( BRANCH_NODE ) );
    new (&(*alloc_branch_data.first)) DefaultBranchNodeType( &tree,
            tree_node_handle(nullptr), alloc_branch_data.second, 1 );

    auto alloc_leaf_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>( NodeHandleType( LEAF_NODE ) );
    new (&(*alloc_leaf_data.first)) DefaultLeafNodeType( &tree,
            alloc_branch_data.second, alloc_leaf_data.second, 0 );

    auto branch_node = alloc_branch_data.first;
    auto leaf_node = alloc_leaf_data.first;
    auto leaf_handle = alloc_leaf_data.second;

    // Add five different branches for the same leaf.
    // This isn't valid in practice but for testing purposes it is fine.
    nirtreedisk::Branch b = createBranchEntry(
            InlineBoundedIsotheticPolygon(
                Rectangle(0.0, 0.0, 1.0, 1.0)), leaf_handle );
    branch_node->addBranchToNode( b );
    b = createBranchEntry( InlineBoundedIsotheticPolygon(
                Rectangle(0.0, 0.0, 2.0, 2.0)), leaf_handle );
    branch_node->addBranchToNode( b );
    b = createBranchEntry( InlineBoundedIsotheticPolygon(
                Rectangle(0.0, 0.0, 3.0, 3.0)), leaf_handle );
    branch_node->addBranchToNode( b );
    b = createBranchEntry( InlineBoundedIsotheticPolygon(
                Rectangle(0.0, 0.0, 4.0, 4.0)), leaf_handle );
    branch_node->addBranchToNode( b );
    b = createBranchEntry( InlineBoundedIsotheticPolygon(
                Rectangle(0.0, 0.0, 5.0, 5.0)), leaf_handle );
    branch_node->addBranchToNode( b );

    REQUIRE( branch_node->cur_offset_ == 5 );

    tree_node_handle packed_handle = branch_node->repack(
            tree.node_allocator_.get(), tree.node_allocator_.get() );
    REQUIRE( packed_handle != nullptr );

    auto packed_branch= tree.node_allocator_->get_tree_node<packed_node>(
            packed_handle );
    
    size_t offset = 0;
    REQUIRE( * (unsigned *) (packed_branch->buffer_ + offset) ==
            branch_node->cur_offset_ );
    offset += sizeof(unsigned);

    // There are 5 entries

    //Entry 1
    REQUIRE( * (tree_node_handle *) (packed_branch->buffer_ + offset) ==
            leaf_handle );
    offset += sizeof(tree_node_handle);
    REQUIRE( * (unsigned *) (packed_branch->buffer_ + offset) ==
            1U );
    offset += sizeof(unsigned);
    REQUIRE( * (Rectangle *) (packed_branch->buffer_ + offset) ==
            Rectangle(0.0,0.0,1.0,1.0) );
    offset += sizeof(Rectangle);

    //Entry 2
    REQUIRE( * (tree_node_handle *) (packed_branch->buffer_ + offset) ==
            leaf_handle );
    offset += sizeof(tree_node_handle);
    REQUIRE( * (unsigned *) (packed_branch->buffer_ + offset) ==
            1U );
    offset += sizeof(unsigned);
    REQUIRE( * (Rectangle *) (packed_branch->buffer_ + offset) ==
            Rectangle(0.0,0.0,2.0,2.0) );
    offset += sizeof(Rectangle);

    //Entry 3
    REQUIRE( * (tree_node_handle *) (packed_branch->buffer_ + offset) ==
            leaf_handle );
    offset += sizeof(tree_node_handle);
    REQUIRE( * (unsigned *) (packed_branch->buffer_ + offset) ==
            1U );
    offset += sizeof(unsigned);
    REQUIRE( * (Rectangle *) (packed_branch->buffer_ + offset) ==
            Rectangle(0.0,0.0,3.0,3.0) );
    offset += sizeof(Rectangle);

    //Entry 4
    REQUIRE( * (tree_node_handle *) (packed_branch->buffer_ + offset) ==
            leaf_handle );
    offset += sizeof(tree_node_handle);
    REQUIRE( * (unsigned *) (packed_branch->buffer_ + offset) ==
            1U );
    offset += sizeof(unsigned);
    REQUIRE( * (Rectangle *) (packed_branch->buffer_ + offset) ==
            Rectangle(0.0,0.0,4.0,4.0) );
    offset += sizeof(Rectangle);

    //Entry 5
    REQUIRE( * (tree_node_handle *) (packed_branch->buffer_ + offset) ==
            leaf_handle );
    offset += sizeof(tree_node_handle);
    REQUIRE( * (unsigned *) (packed_branch->buffer_ + offset) ==
            1U );
    offset += sizeof(unsigned);
    REQUIRE( * (Rectangle *) (packed_branch->buffer_ + offset) ==
            Rectangle(0.0,0.0,5.0,5.0) );
    offset += sizeof(Rectangle);

    tree_node_handle parent_handle( nullptr );
    tree_node_handle unpacked_node_handle =
        nirtreedisk::unpack<DefaultTemplateParams>( packed_handle,
                tree.node_allocator_.get(), &tree, parent_handle );
    auto unpacked_node = tree.node_allocator_.get()->get_tree_node<DefaultBranchNodeType>( unpacked_node_handle );
    REQUIRE( unpacked_node->treeRef == branch_node->treeRef );
    REQUIRE( unpacked_node->parent == branch_node->parent );
    REQUIRE( unpacked_node->cur_offset_ == branch_node->cur_offset_ );
    REQUIRE( unpacked_node->level_ == branch_node->level_ );

    for (unsigned i = 0; i < unpacked_node->cur_offset_; i++) {
        nirtreedisk::Branch &b1 = unpacked_node->entries[i];
        nirtreedisk::Branch &b2 = branch_node->entries[i];

        REQUIRE( b1.child  ==  b2.child );
        REQUIRE( std::holds_alternative<tree_node_handle>(b1.boundingPoly)
                ==  std::holds_alternative<tree_node_handle>(b2.boundingPoly) );

        if ( std::holds_alternative<tree_node_handle>(b1.boundingPoly) ) {
            auto poly1 = tree.node_allocator_->get_tree_node<InlineUnboundedIsotheticPolygon>( std::get<tree_node_handle>(b1.boundingPoly) );
            auto poly2 = tree.node_allocator_->get_tree_node<InlineUnboundedIsotheticPolygon>( std::get<tree_node_handle>(b2.boundingPoly) );

            auto it1 = poly1->begin();
            auto it2 = poly2->begin();

            while (it1 != poly1->end() && it2 != poly2->end()) {
                REQUIRE( *it1 == *it2 );
                it1++;
                it2++;
            }

            REQUIRE( it1 == poly1->end() );
            REQUIRE( it2 == poly2->end() );
        } else {
            REQUIRE( std::get<InlineBoundedIsotheticPolygon>(b1.boundingPoly)
                ==  std::get<InlineBoundedIsotheticPolygon>(b2.boundingPoly) );
        }
    }

    unlink( "nirdiskbacked.txt" );
}

TEST_CASE( "NIRTreeDisk: MergeCmd Simple Tests" ) {
    std::vector<std::pair<Point,uint8_t>> points_with_ownership =
    { {Point(0,0), 0}, {Point(1,0), 0}, {Point(2,0), 1}, {Point(3,0),0}
        };

    REQUIRE( nirtreedisk::get_merge_cmd( points_with_ownership, 1, 0, 0 ) ==
            nirtreedisk::ADD );
    REQUIRE( nirtreedisk::get_merge_cmd( points_with_ownership, 2, 0, 0 ) ==
            nirtreedisk::STOP );
    REQUIRE( nirtreedisk::get_merge_cmd( points_with_ownership, 3, 0, 1 ) ==
            nirtreedisk::STOP );

    points_with_ownership =
    { {Point(0,0), 0}, {Point(1,0), 1}, {Point(2,0), 1}, {Point(3,0),1} };
    REQUIRE( nirtreedisk::get_merge_cmd( points_with_ownership, 1, 0, 0 ) ==
            nirtreedisk::STOP );
    REQUIRE( nirtreedisk::get_merge_cmd( points_with_ownership, 2, 0, 1 ) ==
            nirtreedisk::ADD );
    REQUIRE( nirtreedisk::get_merge_cmd( points_with_ownership, 3, 0, 1 ) ==
            nirtreedisk::ADD );
}

TEST_CASE( "NIRTreeDisk: MergeCmd Vertical Tests" ) {
    std::vector<std::pair<Point,uint8_t>> points_with_ownership =
    { {Point(0,0), 0}, {Point(1,0), 0}, {Point(2,0), 1}, {Point(2,1),1},
        {Point(2,2),1}, {Point(3,0),0}, {Point(4,0),0} };

    REQUIRE( nirtreedisk::get_merge_cmd( points_with_ownership, 1, 0, 0 ) ==
            nirtreedisk::ADD );
    REQUIRE( nirtreedisk::get_merge_cmd( points_with_ownership, 2, 0, 0 ) ==
            nirtreedisk::STOP );
    REQUIRE( nirtreedisk::get_merge_cmd( points_with_ownership, 3, 0, 1 ) ==
            nirtreedisk::ADD );
    REQUIRE( nirtreedisk::get_merge_cmd( points_with_ownership, 4, 0, 1 ) ==
            nirtreedisk::ADD );
    REQUIRE( nirtreedisk::get_merge_cmd( points_with_ownership, 5, 0, 1 ) ==
            nirtreedisk::STOP );
    REQUIRE( nirtreedisk::get_merge_cmd( points_with_ownership, 6, 0, 0 ) ==
            nirtreedisk::ADD );

    points_with_ownership =
    { {Point(0,0), 0}, {Point(1,0), 0}, {Point(2,0), 1}, {Point(2,1),0},
        {Point(2,2),1}, {Point(3,0),0}, {Point(4,0),0} };

    REQUIRE( nirtreedisk::get_merge_cmd( points_with_ownership, 1, 0, 0 ) ==
            nirtreedisk::ADD );
    REQUIRE( nirtreedisk::get_merge_cmd( points_with_ownership, 2, 0, 0 ) ==
            nirtreedisk::CREATE_VERTICAL);
    REQUIRE( nirtreedisk::get_merge_cmd( points_with_ownership, 5, 0, 0 ) ==
            nirtreedisk::ADD );
    REQUIRE( nirtreedisk::get_merge_cmd( points_with_ownership, 5, 0, 0 ) ==
            nirtreedisk::ADD );

}

TEST_CASE( "NIRTreeDisk: DodgeRectangle OwnershipIntersect vertical slice required" ) {
    unlink( "nirdiskbacked.txt" );

    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );

    auto alloc_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                NodeHandleType( LEAF_NODE ) );
    tree_node_handle a_handle = alloc_data.second;
    auto a_leaf_node = alloc_data.first;
    new (&(*a_leaf_node)) DefaultLeafNodeType( &tree, tree_node_handle(nullptr),
            a_handle, 0 );
    a_leaf_node->addPoint( Point(0,1) );
    a_leaf_node->addPoint( Point(2,3) );
    a_leaf_node->addPoint( Point(1.5,1.5) );
    a_leaf_node->addPoint( Point(1.6,1.65) );
    a_leaf_node->addPoint( Point(1.7,1.7) );
    alloc_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                NodeHandleType( LEAF_NODE ) );
    tree_node_handle b_handle = alloc_data.second;
    auto b_leaf_node = alloc_data.first;
    new (&(*b_leaf_node)) DefaultLeafNodeType( &tree, tree_node_handle(nullptr),
            b_handle, 0 );
    b_leaf_node->addPoint( Point(1,0) );
    b_leaf_node->addPoint( Point(1.6,1.6) );
    b_leaf_node->addPoint( Point(1.6,1.7) );
    b_leaf_node->addPoint( Point(3,2) );

    Rectangle a_rect = a_leaf_node->boundingBox();
    Rectangle b_rect = b_leaf_node->boundingBox();
    auto res = nirtreedisk::make_rectangles_disjoint_accounting_for_region_ownership(
        &tree, a_rect, a_handle, b_rect, b_handle );

    Rectangle decomp_a_first =
        Rectangle(0,1,1,nextafter(3.0, DBL_MAX));
    Rectangle decomp_a_second =
        Rectangle(1,nextafter(2,DBL_MAX),
                nextafter(2, DBL_MAX), nextafter(3.0, DBL_MAX));

    Rectangle decomp_a_third =
        Rectangle(1.5, 1, nextafter(1.5,DBL_MAX), nextafter(2,DBL_MAX) );

    Rectangle decomp_a_fourth =
        Rectangle(1.7, 1, nextafter(1.7,DBL_MAX), nextafter(2,DBL_MAX) );

    Rectangle decomp_a_fifth =
        Rectangle(1.6, 1.65, nextafter(1.6,DBL_MAX), nextafter(1.65,DBL_MAX) );

    REQUIRE( res.first.size() == 5 );
    REQUIRE( std::find( res.first.begin(), res.first.end(),
                decomp_a_first ) != res.first.end() );
    REQUIRE( std::find( res.first.begin(), res.first.end(),
                decomp_a_second ) != res.first.end() );
    REQUIRE( std::find( res.first.begin(), res.first.end(),
                decomp_a_third ) != res.first.end() );
    REQUIRE( std::find( res.first.begin(), res.first.end(),
                decomp_a_fourth ) != res.first.end() );
    REQUIRE( std::find( res.first.begin(), res.first.end(),
                decomp_a_fifth ) != res.first.end() );

    REQUIRE( res.second.size() == 4 );
    Rectangle decomp_b_first =
        Rectangle(1,0,nextafter(2,DBL_MAX),1);

    Rectangle decomp_b_second =
        Rectangle(nextafter(2,DBL_MAX),0,nextafter(3,DBL_MAX),nextafter(2,DBL_MAX));

    Rectangle decomp_b_third =
        Rectangle(1.6, 1.6, nextafter(1.6, DBL_MAX), nextafter(1.6,
                    DBL_MAX) );

    Rectangle decomp_b_fourth =
        Rectangle(1.6, 1.7, nextafter(1.6, DBL_MAX), nextafter(1.7,
                    DBL_MAX) );


    REQUIRE( std::find( res.second.begin(), res.second.end(),
                decomp_b_first ) != res.second.end() );
    REQUIRE( std::find( res.second.begin(), res.second.end(),
                decomp_b_second ) != res.second.end() );
    REQUIRE( std::find( res.second.begin(), res.second.end(),
                decomp_b_third ) != res.second.end() );
    REQUIRE( std::find( res.second.begin(), res.second.end(),
                decomp_b_fourth ) != res.second.end() );

    IsotheticPolygon poly;
    poly.basicRectangles = res.first;
    poly.recomputeBoundingBox();

    IsotheticPolygon poly2;
    poly2.basicRectangles = res.second;
    poly2.recomputeBoundingBox();

    REQUIRE( poly.disjoint( poly2 ) );


    unlink( "nirdiskbacked.txt" );
}

TEST_CASE( "NIRTreeDisk: DodgeRectangle OwnershipIntersect vertical slice last el" ) {
    unlink( "nirdiskbacked.txt" );

    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );

    auto alloc_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                NodeHandleType( LEAF_NODE ) );
    tree_node_handle a_handle = alloc_data.second;
    auto a_leaf_node = alloc_data.first;
    new (&(*a_leaf_node)) DefaultLeafNodeType( &tree, tree_node_handle(nullptr),
            a_handle, 0 );
    a_leaf_node->addPoint( Point(0,0) );
    a_leaf_node->addPoint( Point(1.5,0.5) );
    a_leaf_node->addPoint( Point(2,1.8) );
    a_leaf_node->addPoint( Point(2,2) );
    alloc_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                NodeHandleType( LEAF_NODE ) );
    tree_node_handle b_handle = alloc_data.second;
    auto b_leaf_node = alloc_data.first;
    new (&(*b_leaf_node)) DefaultLeafNodeType( &tree, tree_node_handle(nullptr),
            b_handle, 0 );
    b_leaf_node->addPoint( Point(1,1) );
    b_leaf_node->addPoint( Point(1.6,1.6) );
    b_leaf_node->addPoint( Point(2,1.9) );
    b_leaf_node->addPoint( Point(3,1) );
    b_leaf_node->addPoint( Point(3,3) );

    Rectangle a_rect = a_leaf_node->boundingBox();
    Rectangle b_rect = b_leaf_node->boundingBox();
    auto res = nirtreedisk::make_rectangles_disjoint_accounting_for_region_ownership(
        &tree, a_rect, a_handle, b_rect, b_handle );

    REQUIRE( res.first.size() == 4 );
    Rectangle decomp_a_first( 0.0, 0.0, 1.0, nextafter(2.0,DBL_MAX) );
    Rectangle decomp_a_second( 1.0, 0.0, nextafter(2.0,DBL_MAX), 1.0 );
    Rectangle decomp_a_third( 2.0, 2.0, nextafter(2.0,DBL_MAX),
            nextafter(2.0,DBL_MAX) );
    Rectangle decomp_a_fourth( 2.0, 1.8, nextafter(2.0,DBL_MAX),
            nextafter(1.8,DBL_MAX) );


    REQUIRE( std::find( res.first.begin(), res.first.end(),
                decomp_a_first ) != res.first.end() );
    REQUIRE( std::find( res.first.begin(), res.first.end(),
                decomp_a_second ) != res.first.end() );
    REQUIRE( std::find( res.first.begin(), res.first.end(),
                decomp_a_third ) != res.first.end() );
    REQUIRE( std::find( res.first.begin(), res.first.end(),
                decomp_a_fourth ) != res.first.end() );

    REQUIRE( res.second.size() == 4 );
    Rectangle decomp_b_first( 1.0, nextafter(2.0,DBL_MAX), nextafter(2.0, DBL_MAX), nextafter(3.0,DBL_MAX) );
    Rectangle decomp_b_second( nextafter(2.0, DBL_MAX), 1, nextafter(3.0, DBL_MAX), nextafter(3.0,DBL_MAX) );

    Rectangle decomp_b_third( 1.0, 1.0, nextafter(1.6, DBL_MAX),
            nextafter(2.0,DBL_MAX) );
    Rectangle decomp_b_fourth( 2.0, 1.9, nextafter(2.0, DBL_MAX), nextafter(1.9,DBL_MAX) );

    REQUIRE( std::find( res.second.begin(), res.second.end(),
                decomp_b_second) != res.second.end() );
    REQUIRE( std::find( res.second.begin(), res.second.end(),
                decomp_b_second ) != res.second.end() );
    REQUIRE( std::find( res.second.begin(), res.second.end(),
                decomp_b_third ) != res.second.end() );
    REQUIRE( std::find( res.second.begin(), res.second.end(),
                decomp_b_fourth ) != res.second.end() );

    IsotheticPolygon poly;
    poly.basicRectangles = res.first;
    poly.recomputeBoundingBox();

    IsotheticPolygon poly2;
    poly2.basicRectangles = res.second;
    poly2.recomputeBoundingBox();

    REQUIRE( poly.disjoint( poly2 ) );


    unlink( "nirdiskbacked.txt" );
}

TEST_CASE( "NIRTreeDisk: DodgeRectangle Bug1" ) {
    unlink( "nirdiskbacked.txt" );

    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );

    auto alloc_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                NodeHandleType( LEAF_NODE ) );
    tree_node_handle a_handle = alloc_data.second;
    auto a_leaf_node = alloc_data.first;
    new (&(*a_leaf_node)) DefaultLeafNodeType( &tree, tree_node_handle(nullptr),
            a_handle, 0 );
    a_leaf_node->addPoint( Point(0,0) );
    a_leaf_node->addPoint( Point(1.5,0.5) );
    a_leaf_node->addPoint( Point(2,1.8) );
    a_leaf_node->addPoint( Point(2,2) );
    alloc_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                NodeHandleType( LEAF_NODE ) );
    tree_node_handle b_handle = alloc_data.second;
    auto b_leaf_node = alloc_data.first;
    new (&(*b_leaf_node)) DefaultLeafNodeType( &tree, tree_node_handle(nullptr),
            b_handle, 0 );
    b_leaf_node->addPoint( Point(1,1) );
    b_leaf_node->addPoint( Point(1.6,1.6) );
    b_leaf_node->addPoint( Point(2,1.9) );
    b_leaf_node->addPoint( Point(3,1) );
    b_leaf_node->addPoint( Point(3,3) );

    Rectangle a_rect = a_leaf_node->boundingBox();
    Rectangle b_rect = b_leaf_node->boundingBox();
    auto res = nirtreedisk::make_rectangles_disjoint_accounting_for_region_ownership(
        &tree, a_rect, a_handle, b_rect, b_handle );

    REQUIRE( res.first.size() == 4 );
    Rectangle decomp_a_first( 0.0, 0.0, 1.0, nextafter(2.0,DBL_MAX) );
    Rectangle decomp_a_second( 1.0, 0.0, nextafter(2.0,DBL_MAX), 1.0 );
    Rectangle decomp_a_third( 2.0, 2.0, nextafter(2.0,DBL_MAX),
            nextafter(2.0,DBL_MAX) );
    Rectangle decomp_a_fourth( 2.0, 1.8, nextafter(2.0,DBL_MAX),
            nextafter(1.8,DBL_MAX) );


    REQUIRE( std::find( res.first.begin(), res.first.end(),
                decomp_a_first ) != res.first.end() );
    REQUIRE( std::find( res.first.begin(), res.first.end(),
                decomp_a_second ) != res.first.end() );
    REQUIRE( std::find( res.first.begin(), res.first.end(),
                decomp_a_third ) != res.first.end() );
    REQUIRE( std::find( res.first.begin(), res.first.end(),
                decomp_a_fourth ) != res.first.end() );

    REQUIRE( res.second.size() == 4 );
    Rectangle decomp_b_first( 1.0, nextafter(2.0,DBL_MAX), nextafter(2.0, DBL_MAX), nextafter(3.0,DBL_MAX) );
    Rectangle decomp_b_second( nextafter(2.0, DBL_MAX), 1, nextafter(3.0, DBL_MAX), nextafter(3.0,DBL_MAX) );

    Rectangle decomp_b_third( 1.0, 1.0, nextafter(1.6, DBL_MAX),
            nextafter(2.0,DBL_MAX) );
    Rectangle decomp_b_fourth( 2.0, 1.9, nextafter(2.0, DBL_MAX), nextafter(1.9,DBL_MAX) );

    REQUIRE( std::find( res.second.begin(), res.second.end(),
                decomp_b_second) != res.second.end() );
    REQUIRE( std::find( res.second.begin(), res.second.end(),
                decomp_b_second ) != res.second.end() );
    REQUIRE( std::find( res.second.begin(), res.second.end(),
                decomp_b_third ) != res.second.end() );
    REQUIRE( std::find( res.second.begin(), res.second.end(),
                decomp_b_fourth ) != res.second.end() );

    IsotheticPolygon poly;
    poly.basicRectangles = res.first;
    poly.recomputeBoundingBox();

    IsotheticPolygon poly2;
    poly2.basicRectangles = res.second;
    poly2.recomputeBoundingBox();

    REQUIRE( poly.disjoint( poly2 ) );


    unlink( "nirdiskbacked.txt" );
}

TEST_CASE("NIRTreeDisk: unpack in a small out of band polygon") {
    unlink( "nirdiskbacked.txt" );
    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );
    for( unsigned i = 0; i < 5; i++ ) {
        tree.insert( Point(i,i) );
    }

    auto alloc_branch_data =
        tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>( NodeHandleType( BRANCH_NODE ) );
    new (&(*alloc_branch_data.first)) DefaultBranchNodeType( &tree,
            tree_node_handle(nullptr), alloc_branch_data.second, 1 );

    auto alloc_leaf_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>( NodeHandleType( LEAF_NODE ) );
    new (&(*alloc_leaf_data.first)) DefaultLeafNodeType( &tree,
            alloc_branch_data.second, alloc_leaf_data.second, 0 );

    auto branch_node = alloc_branch_data.first;
    auto leaf_node = alloc_leaf_data.first;
    auto leaf_handle = alloc_leaf_data.second;

    // This is small enough to fit
    IsotheticPolygon polygon( Rectangle(0.0, 0.0, 1.0, 1.0) );
    for( double i = 1.0; i < 20.0; i += 1.0 ) {
        polygon.basicRectangles.push_back( Rectangle( i, i, i+1.0,
                    i+1.0) );
    }
    polygon.recomputeBoundingBox();

    auto alloc_poly_data = tree.node_allocator_->create_new_tree_node<InlineUnboundedIsotheticPolygon>(
            compute_sizeof_inline_unbounded_polygon(
                polygon.basicRectangles.size()
                 ), NodeHandleType(BIG_POLYGON));
    new (&(*alloc_poly_data.first)) InlineUnboundedIsotheticPolygon(
            tree.node_allocator_.get(), polygon.basicRectangles.size() );
    alloc_poly_data.first->push_polygon_to_disk( polygon );

    nirtreedisk::Branch b = createBranchEntry( alloc_poly_data.second, leaf_handle );
    branch_node->addBranchToNode( b );

    // too big, should be out of line
    IsotheticPolygon polygon2( Rectangle(0.0, 0.0, 1.0, 1.0) );
    for( double i = 1.0; i < 30.0; i += 1.0 ) {
        polygon2.basicRectangles.push_back( Rectangle( i, i, i+1.0,
                    i+1.0) );
    }
    polygon2.recomputeBoundingBox();

    alloc_poly_data = tree.node_allocator_->create_new_tree_node<InlineUnboundedIsotheticPolygon>(
            compute_sizeof_inline_unbounded_polygon(
                polygon2.basicRectangles.size()
                 ), NodeHandleType(BIG_POLYGON));
    new (&(*alloc_poly_data.first)) InlineUnboundedIsotheticPolygon(
            tree.node_allocator_.get(), polygon2.basicRectangles.size() );
    alloc_poly_data.first->push_polygon_to_disk( polygon2 );

    b = createBranchEntry( alloc_poly_data.second, leaf_handle );
    branch_node->addBranchToNode( b );

    tree_node_handle packed_handle = branch_node->repack(
            tree.node_allocator_.get(), tree.node_allocator_.get() ); 

    auto packed_branch= tree.node_allocator_->get_tree_node<packed_node>(
            packed_handle );

    size_t offset = 0;
    REQUIRE( * (unsigned*) (packed_branch->buffer_ + offset) ==
            branch_node->cur_offset_ );
    offset += sizeof(unsigned);

    // Entry 1
    int new_offset = 0;
    tree_node_handle child_handle = * (tree_node_handle *)
        (packed_branch->buffer_ + offset);
    REQUIRE( child_handle == leaf_handle ); 
    REQUIRE( child_handle.get_associated_poly_is_compressed() == true );
    offset += sizeof( tree_node_handle );

    IsotheticPolygon child_poly = decompress_polygon(
            packed_branch->buffer_ + offset, &new_offset );
    REQUIRE( child_poly.basicRectangles.size() == 20 );
    offset += new_offset;

    child_handle = * (tree_node_handle *)
        (packed_branch->buffer_ + offset);
    REQUIRE( child_handle == leaf_handle ); 
    REQUIRE( child_handle.get_associated_poly_is_compressed() == true );
    offset += sizeof( tree_node_handle );

    child_poly = decompress_polygon( packed_branch->buffer_ + offset );
    REQUIRE( child_poly.basicRectangles.size() == 30 );


    tree_node_handle parent_handle(nullptr);
    tree_node_handle unpacked_node_handle =
        nirtreedisk::unpack<DefaultTemplateParams>( packed_handle,
                tree.node_allocator_.get(), &tree, parent_handle  );
    auto unpacked_node = tree.node_allocator_.get()->get_tree_node<DefaultBranchNodeType>( unpacked_node_handle );
    REQUIRE( unpacked_node->treeRef == branch_node->treeRef );
    REQUIRE( unpacked_node->parent == branch_node->parent );
    REQUIRE( unpacked_node->cur_offset_ == branch_node->cur_offset_ );
    REQUIRE( unpacked_node->level_ == branch_node->level_ );

    for (unsigned i = 0; i < unpacked_node->cur_offset_; i++) {
        nirtreedisk::Branch &b1 = unpacked_node->entries[i];
        nirtreedisk::Branch &b2 = branch_node->entries[i];

        REQUIRE( b1.child  ==  b2.child );
        REQUIRE( std::holds_alternative<tree_node_handle>(b1.boundingPoly)
                ==  std::holds_alternative<tree_node_handle>(b2.boundingPoly) );

        if ( std::holds_alternative<tree_node_handle>(b1.boundingPoly) ) {
            auto poly1 = tree.node_allocator_->get_tree_node<InlineUnboundedIsotheticPolygon>( std::get<tree_node_handle>(b1.boundingPoly) );
            auto poly2 = tree.node_allocator_->get_tree_node<InlineUnboundedIsotheticPolygon>( std::get<tree_node_handle>(b2.boundingPoly) );

            auto it1 = poly1->begin();
            auto it2 = poly2->begin();

            while (it1 != poly1->end() && it2 != poly2->end()) {
                REQUIRE( *it1 == *it2 );
                it1++;
                it2++;
            }

            REQUIRE( it1 == poly1->end() );
            REQUIRE( it2 == poly2->end() );
        } else {
            REQUIRE( std::get<InlineBoundedIsotheticPolygon>(b1.boundingPoly)
                ==  std::get<InlineBoundedIsotheticPolygon>(b2.boundingPoly) );
        }
    }

    unlink("nirdiskbacked.txt");
}

TEST_CASE( "NIRTreeDisk: Repack subtree and then insert" ) {

    unlink( "nirdiskbacked.txt" );

    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );
    auto alloc_branch_data = tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>(
            NodeHandleType(BRANCH_NODE));
    auto branch_handle = alloc_branch_data.second;
    new (&(*alloc_branch_data.first)) DefaultBranchNodeType( &tree,
            tree_node_handle(nullptr), alloc_branch_data.second, 1 );
    auto branch_node = alloc_branch_data.first;

    auto alloc_leaf_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
            NodeHandleType(LEAF_NODE));
    new (&(*alloc_leaf_data.first)) DefaultLeafNodeType( &tree,
            branch_handle, alloc_leaf_data.second, 0 );
    auto leaf_node1 = alloc_leaf_data.first;

    alloc_leaf_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
            NodeHandleType(LEAF_NODE));
    new (&(*alloc_leaf_data.first)) DefaultLeafNodeType( &tree,
            branch_handle, alloc_leaf_data.second, 0 );
    auto leaf_node2 = alloc_leaf_data.first;

    alloc_leaf_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
            NodeHandleType(LEAF_NODE));
    new (&(*alloc_leaf_data.first)) DefaultLeafNodeType( &tree,
            branch_handle, alloc_leaf_data.second, 0 );
    auto leaf_node3 = alloc_leaf_data.first;

    leaf_node1->addPoint( Point(1,1) );
    leaf_node1->addPoint( Point(2,2) );
    leaf_node1->addPoint( Point(3,3) );
    leaf_node1->addPoint( Point(4,4) );
    leaf_node1->addPoint( Point(5,5) );

    Rectangle rect1(1,1,nextafter(1,DBL_MAX),nextafter(5, DBL_MAX));
    Rectangle rect2(7,1,nextafter(20, DBL_MAX),nextafter(1,DBL_MAX));
    Rectangle rect3(-100,-100, -20, -20);
    Rectangle rect4(1,1,nextafter(5,DBL_MAX),nextafter(5,DBL_MAX));
    IsotheticPolygon polygon( rect1 );
    polygon.basicRectangles.push_back( rect2 );
    polygon.basicRectangles.push_back( rect3 );
    polygon.basicRectangles.push_back( rect4 );
    polygon.recomputeBoundingBox();

    nirtreedisk::Branch b = createBranchEntry(
            InlineBoundedIsotheticPolygon(),
            leaf_node1->self_handle_ );
    std::get<InlineBoundedIsotheticPolygon>( b.boundingPoly
            ).push_polygon_to_disk( polygon );
    branch_node->addBranchToNode( b );

    leaf_node2->addPoint( Point(6,6) );
    leaf_node2->addPoint( Point(7,7) );
    leaf_node2->addPoint( Point(8,8) );
    leaf_node2->addPoint( Point(9,9) );
    leaf_node2->addPoint( Point(10,10) );
    b = createBranchEntry(
            InlineBoundedIsotheticPolygon(leaf_node2->boundingBox()),
            leaf_node2->self_handle_ );
    branch_node->addBranchToNode( b );

    leaf_node3->addPoint( Point(11,11) );
    leaf_node3->addPoint( Point(12,12) );
    leaf_node3->addPoint( Point(13,13) );
    leaf_node3->addPoint( Point(14,14) );
    leaf_node3->addPoint( Point(15,15) );
    b = createBranchEntry(
            InlineBoundedIsotheticPolygon(leaf_node3->boundingBox()),
            leaf_node3->self_handle_ );
    branch_node->addBranchToNode( b );

    auto packed_branch_handle =
        nirtreedisk::repack_subtree<3,7,nirtreedisk::LineMinimizeDownsplits>( branch_handle,
            tree.node_allocator_.get(), tree.node_allocator_.get() );

    for( unsigned i = 1; i < 20; i++ ) {
        Point p(i,i);
        if( i < 16 ) {
            REQUIRE( point_search( packed_branch_handle, p, &tree ).size() == 1 );
        } else {
            REQUIRE( point_search( packed_branch_handle, p, &tree ).size() == 0 );
        }
    }

    // Note: Branch handle is the root
    tree.root = packed_branch_handle;
    auto b_node = tree.get_branch_node(packed_branch_handle);
    Point p(68, 68);
    tree.insert(p);

    REQUIRE( point_search( tree.root, p, &tree ).size() == 1 );
    unlink( "nirdiskbacked.txt" );
}

TEST_CASE("NIRTreeDisk: insert after pack simple leaf node") {
    unlink( "nirdiskbacked.txt" );
    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );
    for( unsigned i = 0; i < 5; i++ ) {
        tree.insert( Point(i,i) );
    }

    auto root_leaf_node = tree.get_leaf_node( tree.root );
    REQUIRE( root_leaf_node->cur_offset_ == 5 );

    REQUIRE( root_leaf_node->compute_packed_size() <
            sizeof(DefaultLeafNodeType) );
    REQUIRE( root_leaf_node->compute_packed_size() ==
            sizeof(unsigned) + sizeof(Point) * 5 );

    tree_node_handle repacked_handle = root_leaf_node->repack(
            tree.node_allocator_.get() );

    REQUIRE( repacked_handle != nullptr );
    auto packed_leaf =
        tree.node_allocator_->get_tree_node<packed_node>(
                repacked_handle );
    REQUIRE( * (unsigned *) (packed_leaf->buffer_ ) ==
        root_leaf_node->cur_offset_ );

    Point *p = (Point *) (packed_leaf->buffer_ + sizeof(unsigned) );
    REQUIRE( *(p++) == root_leaf_node->entries.at(0) );
    REQUIRE( *(p++) == root_leaf_node->entries.at(1) );
    REQUIRE( *(p++) == root_leaf_node->entries.at(2) );
    REQUIRE( *(p++) == root_leaf_node->entries.at(3) );
    REQUIRE( *(p++) == root_leaf_node->entries.at(4) );


    Point p1(42, 42);
    auto in = std::variant<nirtreedisk::Branch, Point>(p1);

    std::vector<bool> hasReinsertedOnLevel = {false};
    tree.root = repacked_handle;
    repacked_handle = tree.get_leaf_node(repacked_handle, true)->insert(p1, hasReinsertedOnLevel);
    REQUIRE( tree.search(p1).size() == 1 );
    REQUIRE( point_search( repacked_handle, p1, &tree ).size() == 1 );
    unlink( "nirdiskbacked.txt" );
}

TEST_CASE("NIRTreeDisk: Insert after searching packed leaf from branch." ) {

    unlink( "nirdiskbacked.txt" );
    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );
    for( unsigned i = 0; i < 5; i++ ) {
        tree.insert( Point(i,i) );
    }

    auto root_leaf_node = tree.get_leaf_node( tree.root );
    REQUIRE( root_leaf_node->cur_offset_ == 5 );

    REQUIRE( root_leaf_node->compute_packed_size() <
            sizeof(DefaultLeafNodeType) );
    REQUIRE( root_leaf_node->compute_packed_size() ==
            sizeof(unsigned) + sizeof(Point) * 5 );

    tree_node_handle repacked_handle = root_leaf_node->repack(
            tree.node_allocator_.get() );
    repacked_handle.set_type(
            NodeHandleType(REPACKED_LEAF_NODE) );

    REQUIRE( repacked_handle != nullptr );
    auto packed_leaf =
        tree.node_allocator_->get_tree_node<packed_node>(
                repacked_handle );

    auto alloc_branch_data =
        tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>(
                NodeHandleType(BRANCH_NODE));
    new (&(*alloc_branch_data.first)) DefaultBranchNodeType( &tree,
            tree_node_handle(nullptr), alloc_branch_data.second, 1 );

    auto branch_node = alloc_branch_data.first;

    nirtreedisk::Branch b = createBranchEntry(
            InlineBoundedIsotheticPolygon(root_leaf_node->boundingBox()),
            tree.root );
    branch_node->addBranchToNode( b );

    for( int i = 0; i < 7; i++ ) {
        Point p( i, i );
        auto vec = point_search( branch_node->self_handle_, p, &tree );
        if( i < 5 ) {
            REQUIRE( vec.size() == 1 );
        } else {
            REQUIRE( vec.size() == 0 );
        }
    }

    // So the branch is unpacked but the leaf is repacked.
    Point p(68, 68);
    auto in = std::variant<nirtreedisk::Branch, Point>(p);
    std::vector<bool> hasReinsertedOnLevel = {false, false};

    branch_node->insert(in, hasReinsertedOnLevel);
    REQUIRE( point_search( branch_node->self_handle_, p, &tree ).size() == 1 );
    unlink( "nirdiskbacked.txt");
}

TEST_CASE("NIRTreeDisk: Find parent tests repack." ) {
    unlink( "nirdiskbacked.txt" );
    DefaulTreeType tree( 4096*5, "nirdiskbacked.txt" );

    auto alloc_leaf_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                NodeHandleType(LEAF_NODE) );
    new (&(*(alloc_leaf_data.first))) DefaultLeafNodeType( &tree,
            tree_node_handle(nullptr), alloc_leaf_data.second, 0 );

    auto leaf_node1 = alloc_leaf_data.first;
    leaf_node1->addPoint( Point(0.0,0.0) );
    leaf_node1->addPoint( Point(1.0,1.0) );
    leaf_node1->addPoint( Point(2.0,2.0) );
    leaf_node1->addPoint( Point(3.0,3.0) );
    leaf_node1->addPoint( Point(4.0,4.0) );

    alloc_leaf_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                NodeHandleType(LEAF_NODE) );
    new (&(*(alloc_leaf_data.first))) DefaultLeafNodeType( &tree,
            tree_node_handle(nullptr), alloc_leaf_data.second, 0 );

    auto leaf_node2 = alloc_leaf_data.first;
    leaf_node2->addPoint( Point(10.0,10.0) );
    leaf_node2->addPoint( Point(11.0,11.0) );
    leaf_node2->addPoint( Point(12.0,12.0) );
    leaf_node2->addPoint( Point(13.0,13.0) );
    leaf_node2->addPoint( Point(14.0,14.0) );



    auto alloc_branch_data =
        tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>(
                NodeHandleType(BRANCH_NODE));
    new (&(*alloc_branch_data.first)) DefaultBranchNodeType( &tree,
            tree_node_handle(nullptr), alloc_branch_data.second, 1 );
    auto branch_node1 = alloc_branch_data.first;

    leaf_node1->parent = branch_node1->self_handle_;
    leaf_node2->parent = branch_node1->self_handle_;

    {
        nirtreedisk::Branch b;
        b.child = leaf_node1->self_handle_;
        std::get<InlineBoundedIsotheticPolygon>(b.boundingPoly).push_polygon_to_disk(IsotheticPolygon(leaf_node1->boundingBox()));

        branch_node1->addBranchToNode( b );
    }
    {
        nirtreedisk::Branch b;
        b.child = leaf_node2->self_handle_;
        std::get<InlineBoundedIsotheticPolygon>(b.boundingPoly).push_polygon_to_disk(IsotheticPolygon(leaf_node2->boundingBox()));

        branch_node1->addBranchToNode( b );
    }

    alloc_branch_data =
        tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>(
                NodeHandleType(BRANCH_NODE));
    new (&(*alloc_branch_data.first)) DefaultBranchNodeType( &tree,
            tree_node_handle(nullptr), alloc_branch_data.second, 1 );
    auto wrong_branch = alloc_branch_data.first;

    alloc_leaf_data =
        tree.node_allocator_->create_new_tree_node<DefaultLeafNodeType>(
                NodeHandleType(LEAF_NODE) );
    new (&(*(alloc_leaf_data.first))) DefaultLeafNodeType( &tree,
            tree_node_handle(nullptr), alloc_leaf_data.second, 0 );

    // This will match to force multi-node search
    // when trying to find the parent
    auto leaf_node3 = alloc_leaf_data.first;
    leaf_node3->addPoint( Point(-1.0,-1.0) );
    leaf_node3->addPoint( Point(20.0,20.0) );

    leaf_node3->parent = wrong_branch->self_handle_;
    {
        nirtreedisk::Branch b;
        b.child = leaf_node3->self_handle_;
        std::get<InlineBoundedIsotheticPolygon>(b.boundingPoly).push_polygon_to_disk(IsotheticPolygon(leaf_node3->boundingBox()));

        wrong_branch->addBranchToNode( b );
    }

    alloc_branch_data =
        tree.node_allocator_->create_new_tree_node<DefaultBranchNodeType>(
                NodeHandleType(BRANCH_NODE));
    new (&(*alloc_branch_data.first)) DefaultBranchNodeType( &tree,
            tree_node_handle(nullptr), alloc_branch_data.second, 1 );
    auto branch_node2 = alloc_branch_data.first;

    branch_node1->parent = branch_node2->self_handle_;

    {
        nirtreedisk::Branch b;
        b.child = wrong_branch->self_handle_;
        std::get<InlineBoundedIsotheticPolygon>(b.boundingPoly).push_polygon_to_disk(IsotheticPolygon(wrong_branch->boundingBox()));

        branch_node2->addBranchToNode( b );
    }
    {
        nirtreedisk::Branch b;
        b.child = branch_node1->self_handle_;
        std::get<InlineBoundedIsotheticPolygon>(b.boundingPoly).push_polygon_to_disk(IsotheticPolygon(branch_node1->boundingBox()));
        branch_node2->addBranchToNode( b );
    }

    auto packed_branch_handle =
        nirtreedisk::repack_subtree<3,7,nirtreedisk::LineMinimizeDownsplits>(
                branch_node2->self_handle_,
            tree.node_allocator_.get(), tree.node_allocator_.get() );

    tree.root = packed_branch_handle; 

    // Get child handles
    auto packed_branch =
        tree.node_allocator_->get_tree_node<packed_node>(
                packed_branch_handle );
    char *buffer = packed_branch->buffer_;
    decode_entry_count_and_offset_packed_node( buffer );
    auto packed_wrong_branch_handle = * (tree_node_handle *) (buffer + offset);
    offset += sizeof( tree_node_handle );
    REQUIRE( packed_wrong_branch_handle.get_associated_poly_is_compressed() == false );
    REQUIRE( * (unsigned *) (buffer+offset) == 1U ); // 1 rect
    offset += sizeof( unsigned );
    offset += sizeof( Rectangle ); // skip rectangle
    auto packed_correct_branch_handle = * (tree_node_handle *) (buffer + offset);

    packed_branch =
        tree.node_allocator_->get_tree_node<packed_node>(
                packed_correct_branch_handle );
    buffer = packed_branch->buffer_;
    count = 2;
    offset = sizeof(unsigned);
    auto packed_leaf_handle1 = * (tree_node_handle *) (buffer + offset);
    offset += sizeof( tree_node_handle );
    REQUIRE( packed_leaf_handle1.get_associated_poly_is_compressed() == false );
    REQUIRE( * (unsigned *) (buffer+offset) == 1U ); // 1 rect
    offset += sizeof( unsigned );
    offset += sizeof( Rectangle ); // skip rectangle
    auto packed_leaf_handle2= * (tree_node_handle *) (buffer + offset);

    REQUIRE( tree.get_parent_handle_for_repacked_node( packed_leaf_handle1 )
            == packed_correct_branch_handle );
    REQUIRE( tree.get_parent_handle_for_repacked_node(
                packed_leaf_handle2 ) == packed_correct_branch_handle );
    REQUIRE( tree.get_parent_handle_for_repacked_node(
                packed_correct_branch_handle ) == packed_branch_handle );
    REQUIRE( tree.get_parent_handle_for_repacked_node(
                packed_wrong_branch_handle ) == packed_branch_handle );
    REQUIRE( tree.get_parent_handle_for_repacked_node(
                packed_branch_handle ) == tree_node_handle(nullptr) );

    unlink( "nirdiskbacked.txt" );
}
