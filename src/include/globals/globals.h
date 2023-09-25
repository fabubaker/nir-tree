#ifndef __GLOBALS__
#define __GLOBALS__

#ifndef DIM
#define DIM 2
#endif

const unsigned dimensions = DIM;

// one extra space is saved for insertion which causes split
const unsigned R_STAR_FANOUT = 83;
const unsigned NIR_FANOUT = 83;
#endif
