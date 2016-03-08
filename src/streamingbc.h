#pragma once


#include "stinger.h"
#include "streamingbc_aux.h"


// Note, if you are executing the vertex and edge UPDATES you are REQUIRED to pass on the array of the roots.
// If you are doing an exact computation, then you need to create an array of size NV such that arrayRoot[i]=i;

// Forest creation and destruction. The forest is used for both the static and dynamic computations
bcForest* streamingBCCreateForestExact(uint64_t NV);
bcForest* streamingBCCreateForestApproximate(uint64_t NV, uint64_t NK, uint64_t *rootArrayForApproximation);

void streamingBCDeleteForestExact(bcForestPtr* deadForest);
void streamingBCDeleteForestApproximate(bcForestPtr* deadForest, uint64_t rootArraySize,uint64_t* rootArray);

// Auxillary data structure - creation and destruction
extraArraysPerThread** streamingBCCreateAuxilary(int64_t threadCount,int64_t NV);	
void streamingBCDeleteAuxilary(extraArraysPerThread** parallelExtra, int64_t threadCount, int64_t NV);

// Static BC computation used for the initializations
void streamingBCInitStaticExact(bcForest* forest, struct stinger*stingerGraph, uint64_t NT,extraArraysPerThread** auxilary);
void streamingBCInitStaticApproximate(bcForest* forest, struct stinger*stingerGraph, uint64_t NT,extraArraysPerThread** auxilary,uint64_t NK, uint64_t *rootArrayForApproximation);

// Edge updates
StreamingExtraInfo insertEdgeStreamingBC(bcForest* forest, struct stinger* sStinger,
                       uint64_t newU, uint64_t newV, uint64_t * rootArrayForApproximation,int64_t NK, int64_t NV, int64_t NT, extraArraysPerThread** eAPT, uint32_t loadBalancing);
StreamingExtraInfo deleteEdgeStreamingBC(bcForest *forest, struct stinger *sStinger,
                       uint64_t oldU, uint64_t oldV, uint64_t *rootArrayForApproximation, int64_t NK, int64_t Nv, int64_t NT, extraArraysPerThread **eAPT, uint32_t loadBalancing);

// Vertex updates
StreamingExtraInfo insertVertexStreamingBC(bcForest* forest, struct stinger* sStinger, uint64_t src,
                       uint64_t* adjacencyArray,uint64_t adjacencySize, uint64_t * rootArrayForApproximation,int64_t NK, int64_t NV, int64_t NT, extraArraysPerThread** eAPT);
StreamingExtraInfo deleteVertexStreamingBC(bcForest* forest, struct stinger* sStinger, uint64_t src,
                       uint64_t* adjacencyArray,uint64_t* adjacencySize, uint64_t * rootArrayForApproximation,int64_t NK, int64_t NV, int64_t NT, extraArraysPerThread** eAPT);

