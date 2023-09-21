#ifndef __GLOBALS__
#define __GLOBALS__

#ifndef DIM
#define DIM 2
#endif

const unsigned dimensions = DIM;

const unsigned R_STAR_FANOUT = 84;
// one extra space is saved for insertion which causes split 
const unsigned NIR_FANOUT = 83;
#endif
