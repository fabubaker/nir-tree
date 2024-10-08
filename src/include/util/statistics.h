#ifndef __STATISTICS__
#define __STATISTICS__

#define unlikely(x) __builtin_expect((x),0)

#ifndef STAT
	#define CONST_IF_NOT_STAT const
#else
	#define CONST_IF_NOT_STAT 
#endif

class Statistics {
	public:
		Statistics() :
			histogramSearch(1000000, 0),
			histogramLeaves(1000000, 0),
			histogramRangeSearch(1000000, 0),
			histogramRangeLeaves(100000, 0),
			histogramScatter(100, 0),
    			histogramOutOfLine(100000, 0),
	  		histogramHit(100, 0)
		{
			nodesSearched = 0;
			leavesSearched = 0;
      outOfLineSearched = 0;
		}

		inline void resetSearchTracker( bool isRange )
		{
			// Store stuff in the array
			if (isRange)
			{
				if (unlikely(leavesSearched >= histogramRangeLeaves.size()))
				{
					// 2x so we don't have to resize these often
					histogramRangeLeaves.resize(2*leavesSearched);
				}
				histogramRangeLeaves[leavesSearched]++;

				if (unlikely(nodesSearched >= histogramRangeSearch.size()))
				{
					// 2x so we don't have to resize these often
					histogramRangeSearch.resize(2*nodesSearched);
				}
				histogramRangeSearch.at(nodesSearched)++;

        if (unlikely(outOfLineSearched >= histogramOutOfLine.size()))
        {
          // 2x so we don't have to resize these often
          histogramOutOfLine.resize(2*outOfLineSearched);
        }
        histogramOutOfLine[outOfLineSearched]++;
			}
			else
			{
				if (unlikely(leavesSearched >= histogramLeaves.size()))
				{
					// 2x so we don't have to resize these often
					histogramLeaves.resize(2*leavesSearched);
				}
				histogramLeaves[leavesSearched]++;
				if (unlikely(nodesSearched >= histogramSearch.size()))
				{
					// 2x so we don't have to resize these often
					histogramSearch.resize(2*nodesSearched);
				}

				histogramSearch[nodesSearched]++;
			}

			nodesSearched = 0;
			leavesSearched = 0;
      outOfLineSearched = 0;
		}

		inline void markLeafSearched()
		{
			nodesSearched++;
			leavesSearched++;
		}

		inline void markNonLeafNodeSearched()
		{
			nodesSearched++;
		}

		inline void recordScatter( unsigned scatter ) {
			histogramScatter.at( scatter )++;
		}

    inline void recordIntersectionCount(unsigned intersection_count) {
      intersectionCount += intersection_count;
    }

    inline void recordOutOfLineSearched()
    {
      outOfLineSearched++;
    }

		friend std::ostream& operator<<(std::ostream &os, const Statistics &stats)
		{
			os << "Histogram of Searched Nodes Follows:" << std::endl;
			for (unsigned i = 0; i < stats.histogramSearch.size(); i++)
			{
				if (stats.histogramSearch[i] > 0)
				{
					os << "  " << i << " : " << stats.histogramSearch[i] << std::endl;
				}
			}
			os << "Histogram of Searched Leaves Follows:" << std::endl;
			for (unsigned i = 0; i < stats.histogramLeaves.size(); i++)
			{
				if (stats.histogramLeaves[i] > 0)
				{
					os << "  " << i << " : " << stats.histogramLeaves[i] << std::endl;
				}
			}
			os << "Histogram of Range Searched Nodes Follows:" << std::endl;
			for (unsigned i = 0; i < stats.histogramRangeSearch.size(); i++)
			{
				if (stats.histogramRangeSearch[i] > 0)
				{
					os << "  " << i << " : " << stats.histogramRangeSearch[i] << std::endl;
				}
			}
			os << "Histogram of Range Searched Leaves Follows:" << std::endl;
			for (unsigned i = 0; i < stats.histogramRangeLeaves.size(); i++)
			{
				if (stats.histogramRangeLeaves[i] > 0)
				{
					os << "  " << i << " : " << stats.histogramRangeLeaves[i] << std::endl;
				}
			}

			os << "Histogram of Scatter Follows:" << std::endl;
			for (unsigned i = 0; i < stats.histogramScatter.size(); i++)
			{
				if (stats.histogramScatter[i] > 0)
				{
					os << "  " << i << " : " << stats.histogramScatter[i] << std::endl;
				}
			}

      os << "Histogram of Out of Line Follows:" << std::endl;
      for (unsigned i = 0; i < stats.histogramOutOfLine.size(); i++)
      {
        if (stats.histogramOutOfLine[i] > 0)
        {
          os << "  " << i << " : " << stats.histogramOutOfLine[i] << std::endl;
        }
      }

			return os;
		}

  size_t intersectionCount = 0;
  std::vector<unsigned> histogramHit;
	private:
		std::vector<unsigned> histogramSearch;
		std::vector<unsigned> histogramLeaves;
		std::vector<unsigned> histogramRangeSearch;
		std::vector<unsigned> histogramRangeLeaves;
		std::vector<unsigned> histogramScatter;
    std::vector<unsigned> histogramOutOfLine;
		unsigned nodesSearched;
		unsigned leavesSearched;
    unsigned outOfLineSearched;
};

void printFanoutHistogram(
        std::vector<std::vector<unsigned long>> histogramFanoutAtLevel,
        unsigned treeHeight
);

void printPolygonHistogram(
        std::vector<std::vector<unsigned long>> histogramPolygonAtLevel,
        unsigned treeHeight
);

#ifdef STAT
	#include <iostream>

	#define STATEXEC(e) e
	#define STATMEM(mem) std::cout << "Memory Usage: " << (mem / 1024) << "KB, " << (mem / (1024 * 1024)) << "MB, " << (mem / (1024 * 1024 * 1024)) << "GB" << std::endl
	#define STATHEIGHT(height) std::cout << "Tree Height: " << height << std::endl
	#define STATSIZE(n) std::cout << "Tree Nodes: " << n << std::endl
	#define STATSINGULAR(n) std::cout << "Tree Nodes w/fanout=1: " << n << std::endl
	#define STATLEAF(n) std::cout << "Tree Leaves: " << n << std::endl
	#define STATBRANCH(branches) std::cout << "Tree Branches: " << branches << std::endl
	#define STATCOVER(c) std::cout << "Total Coverage: " << c << std::endl;
	#define STATOVERLAP(o) std::cout << "Total Overlap: " << o << std::endl;
	#define STATAVGCOVER(c) std::cout << "Avg Coverage Per Node: " << c << std::endl;
	#define STATAVGOVERLAP(o) std::cout << "Avg Overlap Per Node: " << o << std::endl;
	#define STATFANHIST() std::cout << "Histogram of Fanout Follows: " << std::endl
	#define STATLINES(n) std::cout << "Bounding Lines: " << n << std::endl
	#define STATTOTALPOLYSIZE(n) std::cout << "Total Polygon Size: " << n << std::endl
	#define STATPOLYHIST() std::cout << "Histogram of Polygon Sizes Follows:" << std::endl
	#define STATSEARCHHIST() std::cout << "Histogram of Searched Nodes Follows:" << std::endl
	#define STATLEAVESHIST() std::cout << "Histogram of Searched Leaves Follows:" << std::endl
	#define STATRANGESEARCHHIST() std::cout << "Histogram of Range Searched Nodes Follows:" << std::endl
	#define STATRANGELEAVESHIST() std::cout << "Histogram of Range Searched Leaves Follows:" << std::endl
	#define STATHIST(bucket, count) std::cout << "  " << bucket << " : " << count << std::endl
	#define STATHISTATLEVEL(bucket, count, level) std::cout << " L-" << level << "  " << bucket << " : " << count << std::endl
	#define STATFANHISTATLEVEL() std::cout << "Histogram of Fanout for each level Follows: " << std::endl
	#define STATPOLYHISTATLEVEL() std::cout << "Histogram of Polygon Sizes for each level Follows: " << std::endl
#else 
	#define STATEXEC(e)
	#define STATMEM(mem)
	#define STATHEIGHT(height)
	#define STATSIZE(n)
	#define STATSINGULAR(n)
	#define STATLEAF(n)
	#define STATBRANCH(branches)
	#define STATFANHIST()
	#define STATLINES(n)
	#define STATTOTALPOLYSIZE(n)
	#define STATPOLYHIST()
	#define STATCOVER(c)
	#define STATOVERLAP(o)
	#define STATAVGCOVER(c)
	#define STATAVGOVERLAP(o)
	#define STATSEARCHHIST()
	#define STATLEAVESHIST()
	#define STATRANGESEARCHHIST()
	#define STATRANGELEAVESHIST()
	#define STATHIST(bucket, count)
	#define STATHISTATLEVEL(bucket, count, level)
	#define STATFANHISTATLEVEL()
	#define STATPOLYHISTATLEVEL()

#endif

#endif
