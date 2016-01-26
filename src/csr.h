
#pragma once

#include "stinger.h"

typedef struct{
    int64_t NV;
    int64_t NE;
    int64_t* edgeArray;
    int64_t* vertexPointerArray;
} csrGraph;


// NE number of undirected edges.
csrGraph* CreateCSRFromStinger(struct stinger* stingerGraph,int64_t NV,int64_t NE);
void CreateStingerFromCSR(csrGraph* csr,struct stinger** stingerGraph);
//csrGraph* CreateEmptyCSR(int64_t NV,int64_t NE);
csrGraph* CreateEmptyUndirectedCSR(int64_t NV,int64_t NE);
csrGraph* CreateEmptyDirectedCSR(int64_t NV,int64_t NE);

void FreeCSR(csrGraph* graph);
