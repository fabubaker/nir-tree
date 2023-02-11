#include <bulk_load.h>


void generate_tree( std::map<std::string, size_t> &configU ) {

    std::string backing_file = "bulkloaded_tree.txt";
    unlink( backing_file.c_str() );

    std::vector<Point> all_points;
    std::optional<Point> next;

    if( configU["distribution"] == CALIFORNIA ) {
        PointGenerator<BenchTypeClasses::California> points;
        while( (next = points.nextPoint() )) {
            all_points.push_back( next.value() );
        }
    } else if( configU["distribution"] == UNIFORM ) {
        BenchTypeClasses::Uniform::size = configU["size"];
        BenchTypeClasses::Uniform::dimensions = dimensions;
        BenchTypeClasses::Uniform::seed = configU["seed"];
        if (configU.count("precision")) {
            BenchTypeClasses::Uniform::precision = configU["precision"];
            std::cout << "Using precision " << configU["precision"] << std::endl;
        }
        PointGenerator<BenchTypeClasses::Uniform> points;
        while( (next = points.nextPoint() )) {
            all_points.push_back( next.value() );
        }
    } else if (configU["distribution"] == ZIPF) {
        BenchTypeClasses::Zipf::size = configU["size"];
        BenchTypeClasses::Zipf::dimensions = dimensions;
        BenchTypeClasses::Zipf::seed = configU["seed"];
        BenchTypeClasses::Zipf::num_elements = configU["num_elements"];
        BenchTypeClasses::Zipf::alpha = configU["alpha"];

        PointGenerator<BenchTypeClasses::Zipf> points;
        while( (next = points.nextPoint() )) {
            std::cout << next.value() << std::endl;
            all_points.push_back( next.value() );
        }
    } else if( configU["distribution"] == GAUSS ) {
        BenchTypeClasses::Gauss::size = configU["size"];
        BenchTypeClasses::Gauss::dimensions = dimensions;
        BenchTypeClasses::Gauss::seed = configU["seed"];

        PointGenerator<BenchTypeClasses::Gauss> points;
        while( (next = points.nextPoint() )) {
            all_points.push_back( next.value() );
        }
    }

    double bulk_load_pct = 1.0;
    uint64_t cut_off_bulk_load = std::floor(bulk_load_pct*all_points.size());
    std::cout << "Bulk loading " << cut_off_bulk_load << " points." << std::endl;
    std::cout << "Sequential Inserting " << all_points.size() - cut_off_bulk_load << " points." << std::endl;

    Index *spatialIndex;
    if( configU["tree"] == NIR_TREE ) {
        // FIXME: one of the main slow downs of bulk loading in the NIR-Tree is going to be that
        // it looks for compression opportunities during the bulk load. This makes no sense because
        // we are guaranteed that each generated rectangle is disjoint.
        nirtreedisk::NIRTreeDisk<5,9,nirtreedisk::ExperimentalStrategy> *tree =  new
            nirtreedisk::NIRTreeDisk<5,9,nirtreedisk::ExperimentalStrategy>(
                    40960UL*130000UL, backing_file );
        std::cout << "Bulk Loading..." << std::endl;
        std::cout << "Creating tree with " << 40960UL *130000UL << "bytes" << std::endl;
        bulk_load_tree( tree, configU, all_points.begin(), all_points.begin() + cut_off_bulk_load, 9 );
        std::cout << "Created NIRTree." << std::endl;
        
        std::cout << "Creating consolidator..." << std::endl;
        auto consolidated_allocator = std::make_unique<tree_node_allocator>( 40960UL * 1300UL, "consolidated_nirtree.txt" );
        consolidated_allocator->initialize();
        std::cout << "Repacking into consolidator..." << std::endl;
        tree_node_handle new_root = nirtreedisk::repack_subtree<5,9,nirtreedisk::ExperimentalStrategy>( tree->root, get_node_allocator( tree ), consolidated_allocator.get() );
        std::cout << "Done." << std::endl;

        std::cout << "Swapping out allocator..." << std::endl;
        tree->node_allocator_ = std::move( consolidated_allocator );
        tree->root = new_root;
        tree->write_metadata();

        tree->stat();

        nirtreedisk::tree_validate_recursive( tree->root, tree->node_allocator_.get() );

        /*
        if( !tree->validate() ) {
            std::cout << "Tree Validation Failed" << std::endl;
        }
        */
        spatialIndex = tree;
    } else if( configU["tree"] == R_STAR_TREE ) {
        rstartreedisk::RStarTreeDisk<5,9> *tree = new rstartreedisk::RStarTreeDisk<5,9>(
                    40960UL*130000UL, backing_file );
        std::cout << "Bulk Loading..." << std::endl;
        bulk_load_tree( tree, configU, all_points.begin(), all_points.begin() + cut_off_bulk_load, 9 );
        std::cout << "Created R*Tree" << std::endl;
        spatialIndex = tree;

        tree->stat();
        exit(0);

    } else {
        abort();
    }

    std::mt19937 g;
    g.seed(0);

    std::shuffle( all_points.begin(), all_points.end(), g );

    unsigned totalSearches  = 0;
	double totalTimeSearches = 0.0;
    for( Point p : all_points ) {
        std::chrono::high_resolution_clock::time_point begin = std::chrono::high_resolution_clock::now();
        std::cout << "Searching: for " << p << std::endl;
        std::vector<Point> out = spatialIndex->search(p);
        if( out.size() != 1 ) {
            std::cout << "Could not find " << p << std::endl;
            std::cout << out.size() << std::endl;
            std::cout << "Total successful searches: " << totalSearches << std::endl;
            abort();
        }
        std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> delta = std::chrono::duration_cast<std::chrono::duration<double>>(end - begin);
        totalTimeSearches += delta.count();
        totalSearches += 1;
        if( totalSearches >= 10000 ) {
            break;
        }
    }

    spatialIndex->stat();

	std::cout << "Total time to search: " << totalTimeSearches << "s" << std::endl;
	std::cout << "Avg time to search: " << totalTimeSearches / totalSearches << "s" << std::endl;

    return;
}



int main( int argc, char **argv ) {

    std::map<std::string, size_t> configU;
    configU.emplace( "tree", NIR_TREE );
    configU.emplace( "distribution", CALIFORNIA);
    configU.emplace( "seed", 0 );

    int option;

    while( (option = getopt(argc,argv, "t:m:n:s:p:g:z:")) != -1 ) {
        switch( option ) {
            case 't': {
                configU["tree"] = (TreeType)atoi(optarg);
                break;
            }
            case 'm': {
                configU["distribution"] = (BenchType)atoi(optarg);
                break;
            }
            case 'n': {
                configU["size"] = atoi(optarg);
                break;
            }
            case 'g': {
                configU["num_elements"] = atoi(optarg);
                break;
            }
            case 'z': {
                configU["alpha"] = std::stod(optarg);
                break;
            }
            case 's': {
                configU["seed"] = atoi(optarg);
                break;
            }
            case 'p': {
                configU["precision"] = atoi(optarg);
                break;
            }
        }
    }

    generate_tree( configU );

    return 0;
}
