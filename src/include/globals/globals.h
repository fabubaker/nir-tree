#ifndef __GLOBALS__
#define __GLOBALS__

#ifndef DIM
#define DIM 2
#endif

const unsigned dimensions = DIM;

const unsigned R_STAR_FANOUT = 84;
//const unsigned NIR_FANOUT = 84;
// shirley: leave a slot for possible overflow caused by insertion 
const unsigned NIR_FANOUT = 83;
#endif
