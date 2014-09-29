
#pragma once


typedef enum{
    TC_ER = 0,
    TC_RMAT,
    TC_REAL_DATA,       // DOES NOT USE ADJACENCY MATRIX
    TC_CC_ER,
    TC_SAME_GRAPH_ER,
    TC_SAME_GRAPH_RMAT,
} testCase;

typedef enum{
    UP_INSERT = 0,
    UP_DELETE,
} updateType;


extern int64_t noComputeDiffComp;

extern int64_t NV;
extern int64_t NE;
extern int64_t NT;
extern int64_t NK;
extern int64_t randomSeed;
extern int64_t iterationCount;
char initial_graph_name[1024];
extern updateType opType;



extern testCase graphTestCase; // GraphType ER=0,RMAT=1,real graph=2, connec

typedef struct
{
    uint64_t movement;
    uint64_t adjacent;
    uint64_t sameLevel;
    uint64_t connectedComponents;
}StreamingExtraInfo;




