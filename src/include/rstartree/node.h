#ifndef __RSTARTREE_NODE__
#define __RSTARTREE_NODE__

#include <cassert>
#include <vector>
#include <stack>
#include <map>
#include <list>
#include <utility>
#include <cmath>
#include <numeric>
#include <iostream>
#include <limits>
#include <memory>
#include <variant>
#include <globals/globals.h>
#include <util/geometry.h>
#include <util/statistics.h>

namespace rstartree
{
	class RStarTree;

	class Node
	{
		private:

			RStarTree &treeRef;

			void searchSub(const Point &requestedPoint, std::vector<Point> &accumulator) CONST_IF_NOT_STAT;
			void searchSub(const Rectangle &rectangle, std::vector<Point> &accumulator) CONST_IF_NOT_STAT;

		public:
			class Branch
			{
				public:
					Rectangle boundingBox;
					Node *child;

					Branch(Rectangle boundingBox, Node *child) : boundingBox(boundingBox), child(child) {}
					Branch(const Branch &other) : boundingBox(other.boundingBox), child(other.child) {}

					bool operator==(const Branch &o) const;
			};
			typedef std::variant<Point, Branch> NodeEntry;

			Node *parent;
			std::vector<NodeEntry> entries;
			unsigned level;

			// Constructors and destructors
			Node(RStarTree &treeRef, Node *p=nullptr, unsigned level=0);
			void deleteSubtrees();

			// Helper functions
			Rectangle boundingBox() const;
			bool updateBoundingBox(Node *child, Rectangle updatedBoundingBox);
			void removeChild(Node *child);
			void removeData(const Point &givenPoint);
			Node *chooseSubtree(const NodeEntry &nodeEntry);
			Node *findLeaf(const Point &givenPoint);
			inline bool isLeafNode() const { return level == 0; }
			unsigned chooseSplitLeafAxis();
			unsigned chooseSplitNonLeafAxis();
			unsigned chooseSplitAxis();
			unsigned chooseSplitIndex(unsigned axis);
			Node *splitNode();
			Node *adjustTree(Node *siblingLeaf, std::vector<bool> &hasReinsertedOnLevel);
			Node *reInsert(std::vector<bool> &hasReinsertedOnLevel);
			Node *overflowTreatment(std::vector<bool> &hasReinsertedOnLevel);
			Node *condenseTree(std::vector<bool> &hasReinsertedOnLevel);

			// Datastructure interface functions
			void exhaustiveSearch(const Point &requestedPoint, std::vector<Point> &accumulator) const;

			std::vector<Point> search(const Point &requestedPoint) CONST_IF_NOT_STAT;
			std::vector<Point> search(const Rectangle &requestedRectangle) CONST_IF_NOT_STAT;

			// These return the root of the tree.
			Node *insert(NodeEntry nodeEntry, std::vector<bool> &hasReinsertedOnLevel);
			Node *remove(Point &givenPoint, std::vector<bool> hasReinsertedOnLevel);

			// Miscellaneous
			unsigned checksum() const;
			void print() const;
			void printTree() const;
			unsigned height() const;
			void stat() const;

			// Operators
			bool operator<(const Node &otherNode) const;
	};

	Rectangle boxFromNodeEntry(const Node::NodeEntry &entry);
	double computeOverlapGrowth(unsigned index, const std::vector<Node::NodeEntry> &entries, const Rectangle &rect);
}

#endif
