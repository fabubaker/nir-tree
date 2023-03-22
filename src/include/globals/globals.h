#ifndef __GLOBALS__
#define __GLOBALS__

#ifndef DIM
#define DIM 2
#endif

const unsigned dimensions = DIM;
// Buffer pool memory for bulk-loading.
#define GEN_TREE_BUFFER_POOL_MEMORY 40960 * 130000

// Buffer pool memory for running main benchmarks.
#define MAIN_BUFFER_POOL_MEMORY 4096 * 13000

// The two above are separate since you might wanna benchmark
// both with different buffer pool sizes.

#endif
