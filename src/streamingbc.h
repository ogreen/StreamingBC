#pragma once


#include "stinger.h"
#include "streamingbc_aux.h"


typedef enum {
    COARSE,
    FINE,
} granularity_t;

typedef enum {
    DELETE,
    INSERT,
} operation_t;

typedef enum {
    NO_BALANCE,
    BALANCE,
} lbalance_t;

typedef enum {
    UP_INSERT = 0,
    UP_DELETE,
} updateType;

typedef struct {
    int64_t NV;
    int64_t NE;
    int64_t * edgeArray;
    int64_t * vertexPointerArray;
} csrGraph;

// Note, if you are executing the vertex and edge UPDATES you are REQUIRED to pass on the array of the roots.
// If you are doing an exact computation, then you need to create an array of size NV such that arrayRoot[i]=i;

// Forest creation and destruction. The forest is used for both the static and dynamic computations
bcForest * streamingBCCreateForestExact(uint64_t NV);
bcForest * streamingBCCreateForestApproximate(uint64_t NV, uint64_t NK, uint64_t * rootArrayForApproximation);

void streamingBCDeleteForestExact(bcForestPtr * deadForest);
void streamingBCDeleteForestApproximate(bcForestPtr * deadForest, uint64_t rootArraySize, uint64_t * rootArray);

// Auxillary data structure - creation and destruction
extraArraysPerThread ** streamingBCCreateAuxilary(int64_t threadCount, int64_t NV);
void streamingBCDeleteAuxilary(extraArraysPerThread ** parallelExtra, int64_t threadCount, int64_t NV);

// Static BC computation used for the initializations
void streamingBCInitStaticExact(bcForest * forest, struct stinger * stingerGraph, uint64_t NT, extraArraysPerThread ** auxilary);
void streamingBCInitStaticApproximate(bcForest * forest, struct stinger * stingerGraph, uint64_t NT, extraArraysPerThread ** auxilary, uint64_t NK, uint64_t * rootArrayForApproximation);

// Edge updates
StreamingExtraInfo insertEdgeStreamingBC(bcForest * forest, struct stinger * sStinger,
        uint64_t newU, uint64_t newV, uint64_t * rootArrayForApproximation, int64_t NK, int64_t NV, int64_t NT, extraArraysPerThread ** eAPT, uint32_t loadBalancing, uint32_t granularity);
StreamingExtraInfo deleteEdgeStreamingBC(bcForest * forest, struct stinger * sStinger,
        uint64_t oldU, uint64_t oldV, uint64_t * rootArrayForApproximation, int64_t NK, int64_t Nv, int64_t NT, extraArraysPerThread ** eAPT, uint32_t loadBalancing, uint32_t granularity);

// Vertex updates
StreamingExtraInfo insertVertexStreamingBC(bcForest * forest, struct stinger * sStinger, uint64_t src,
        uint64_t * adjacencyArray, uint64_t adjacencySize, uint64_t * rootArrayForApproximation, int64_t NK, int64_t NV, int64_t NT, extraArraysPerThread ** eAPT);
StreamingExtraInfo deleteVertexStreamingBC(bcForest * forest, struct stinger * sStinger, uint64_t src,
        uint64_t * adjacencyArray, uint64_t * adjacencySize, uint64_t * rootArrayForApproximation, int64_t NK, int64_t NV, int64_t NT, extraArraysPerThread ** eAPT);

// sbcMain methods.
csrGraph * CreateCSRFromStinger(struct stinger * stingerGraph, int64_t NV, int64_t NE);
void CreateStingerFromCSR(csrGraph * csr, struct stinger ** stingerGraph);
csrGraph * CreateEmptyUndirectedCSR(int64_t NV, int64_t NE);
csrGraph * CreateEmptyDirectedCSR(int64_t NV, int64_t NE);

void FreeCSR(csrGraph * graph);

void hostParseArgsVitalUpdate(int argc, char ** argv, int64_t * NV, int64_t * NE, int64_t * NK, int64_t * NT,
                              int64_t * randomSeed, int64_t * iterationCount, char * initial_graph_name[1024],
                              operation_t * operation, lbalance_t * loadBalancing, granularity_t * granularity, int64_t * edgeCount);

void CreateRandomEdgeListFromGraph(struct stinger * stingerGraph, int64_t NV, int64_t * insertionArraySrc,
                                   int64_t * insertionArrayDest, int64_t insertionCount);

void ParseDimacsGraph(char initial_graph_name[1024], struct stinger ** stingerGraph, int64_t * NV, int64_t * NE);

void PickRoots(uint64_t * rootArrayForApproximation, int64_t NK, int64_t NV);

float InsertEdge(int64_t count, int64_t * insertionArraySrc, int64_t * insertionArrayDest, struct stinger ** stingerGraph, 
                StreamingExtraInfo * oneSEI, extraArraysPerThread ** eAPT_perThread, int64_t * rootArrayForApproximation, int64_t NK, int64_t NV, int64_t NT, 
                bcForest * beforeBCForest, lbalance_t loadBalancing, granularity_t granularity, int64_t * iterationCount);

float DeleteEdge(int64_t count, int64_t * deletionArraySrc, int64_t * deletionArrayDest, struct stinger ** stingerGraph, 
                StreamingExtraInfo * oneSEI, extraArraysPerThread ** eAPT_perThread, int64_t * rootArrayForApproximation, int64_t NK, int64_t NV, int64_t NT, 
                bcForest * beforeBCForest, lbalance_t loadBalancing, granularity_t granularity, int64_t * iterationCount);

void CompareDynamicWithExactResults(bcForest * beforeBCForest, struct stinger * stingerGraph, int64_t * rootArrayForApproximation, int64_t NK, int64_t NT, int64_t NV);

double updateEdgeNEW(struct stinger * stingerGraph, StreamingExtraInfo * oneSEI,
                     extraArraysPerThread ** eAPT_perThread, uint64_t * rootArrayForApproximation, int64_t NK,
                     int64_t NV, int64_t NT, bcForest * beforeBCForest, int64_t u_new, int64_t v_new, int64_t * iterationCount,
                     lbalance_t loadBalancing, granularity_t granularity);


double deleteEdgeNEW(struct stinger * stingerGraph, StreamingExtraInfo * oneSEI, extraArraysPerThread ** eAPT_perThread,
                     uint64_t * rootArrayForApproximation, int64_t NK, int64_t NV, int64_t NT,
                     bcForest * beforeBCForest, int64_t u_old, int64_t v_old, int64_t * iterationCount, lbalance_t loadBalancing, granularity_t granularity);
