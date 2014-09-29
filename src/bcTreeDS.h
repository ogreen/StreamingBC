#pragma once

#include "streaming_utils.h"

#if PARENT_INFO_ON
#include "list.h"
#endif

typedef struct{
    int64_t level;
    int64_t pathsToRoot;
    bc_t delta;


}bcV;

typedef struct {
    int64_t NV;
    bcV* vArr;
	double dummy[4];
} bcTree;

typedef bcTree* bcTreePtr;
typedef bcTreePtr* bcTreePtrPtr;

typedef struct {
	bcTreePtr* forest;
    bc_t* totalBC;
    int64_t NV;
} bcForest;

typedef bcForest* bcForestPtr;
typedef bcForestPtr* bcForestPtrPtr;


bcForest* CreateForestForApproxCase(bcForest** newForest, int64_t numVertices, uint64_t* rootArray, uint64_t rootArraySize);
void DestroyForestForApproxCase(bcForest** deadForest, uint64_t* rootArray, uint64_t rootArraySize);

