#ifndef __GLOBALS__
#define __GLOBALS__

#ifndef DIM
#define DIM 2
#endif

const unsigned dimensions = DIM;

// one extra space is saved for insertion which causes split
const unsigned R_STAR_MIN_FANOUT = 5;
const unsigned R_STAR_MAX_FANOUT = 83;
const unsigned R_PLUS_MIN_FANOUT = 5;
const unsigned R_PLUS_MAX_FANOUT = 83;
const unsigned NIR_MIN_FANOUT = 5;
const unsigned NIR_MAX_FANOUT = 83;
#endif
