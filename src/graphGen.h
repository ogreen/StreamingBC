#pragma once


#include "stinger.h"
#include "csr.h"

typedef enum{
    GT_ER=0,
    GT_RMAT,
}graphType;

void uniformRandGraph(uint64_t NV, uint64_t NE, struct stinger ** stingerGraph, uint64_t randomSeed);
void FillAdjacencyMat(uint64_t NV, uint64_t NE, uint64_t * graph[],struct stinger * stingerGraph, uint64_t randomSeed);
void uniformRandGraphMin2Components(uint64_t NV, uint64_t NE, uint64_t * graph[], struct stinger ** stingerGraph, uint64_t randomSeed, float percentageOfBigComponent);
void CreateUniformEdge(uint64_t* out_u, uint64_t* out_v, uint64_t NV,uint64_t* graph[]);


void rMatRandGraph(uint64_t NV, uint64_t NE, uint64_t * graph[],struct stinger ** stingerGraph, uint64_t randomSeed,int64_t SCALE);

void CreateRMatEdge(uint64_t* out_u, uint64_t* out_v, uint64_t NV,int64_t SCALE);
void CreateCSRGraph(csrGraph* csrG,uint64_t NV, uint64_t NE,uint64_t randomSeed,graphType gt, uint64_t SCALE);
