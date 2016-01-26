
#pragma once

#include <stdint.h> 


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



typedef struct
{
    uint64_t movement;
    uint64_t adjacent;
    uint64_t sameLevel;
    uint64_t connectedComponents;
}StreamingExtraInfo;




