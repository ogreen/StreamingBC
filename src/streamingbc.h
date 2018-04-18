#pragma once


#include "stinger.h"
#include "streamingbc_internal.h"
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
bcForest * streamingBCCreateForestExact(int64_t NV);
bcForest * streamingBCCreateForestApproximate(int64_t NV, int64_t NK, int64_t * rootArrayForApproximation);

void streamingBCDeleteForestExact(bcForestPtr * deadForest);
void streamingBCDeleteForestApproximate(bcForestPtr * deadForest, int64_t rootArraySize, int64_t * rootArray);

// Auxillary data structure - creation and destruction
extraArraysPerThread ** streamingBCCreateAuxilary(int64_t threadCount, int64_t NV);
void streamingBCDeleteAuxilary(extraArraysPerThread ** parallelExtra, int64_t threadCount, int64_t NV);

// Static BC computation used for the initializations
void streamingBCInitStaticExact(bcForest * forest, struct stinger * stingerGraph, int64_t NT, extraArraysPerThread ** auxilary);
void streamingBCInitStaticApproximate(bcForest * forest, struct stinger * stingerGraph, int64_t NT, extraArraysPerThread ** auxilary, int64_t NK, int64_t * rootArrayForApproximation);

// Edge updates
StreamingExtraInfo insertEdgeStreamingBC(bcForest * forest, struct stinger * sStinger,
        int64_t newU, int64_t newV, int64_t * rootArrayForApproximation, int64_t NK, int64_t NV, int64_t NT, extraArraysPerThread ** eAPT, uint32_t loadBalancing, uint32_t granularity);
StreamingExtraInfo deleteEdgeStreamingBC(bcForest * forest, struct stinger * sStinger,
        int64_t oldU, int64_t oldV, int64_t * rootArrayForApproximation, int64_t NK, int64_t Nv, int64_t NT, extraArraysPerThread ** eAPT, uint32_t loadBalancing, uint32_t granularity);

// Vertex updates
StreamingExtraInfo insertVertexStreamingBC(bcForest * forest, struct stinger * sStinger, int64_t src,
        int64_t * adjacencyArray, int64_t adjacencySize, int64_t * rootArrayForApproximation, int64_t NK, int64_t NV, int64_t NT, extraArraysPerThread ** eAPT);
StreamingExtraInfo deleteVertexStreamingBC(bcForest * forest, struct stinger * sStinger, int64_t src,
        int64_t * adjacencyArray, int64_t * adjacencySize, int64_t * rootArrayForApproximation, int64_t NK, int64_t NV, int64_t NT, extraArraysPerThread ** eAPT);

// sbcMain methods.
csrGraph * CreateCSRFromStinger(struct stinger * stingerGraph, int64_t NV, int64_t NE);
void CreateStingerFromCSR(csrGraph * csr, struct stinger ** stingerGraph);
csrGraph * CreateEmptyUndirectedCSR(int64_t NV, int64_t NE);
csrGraph * CreateEmptyDirectedCSR(int64_t NV, int64_t NE);

void FreeCSR(csrGraph * graph);

void hostParseArgsVitalUpdate(int argc, char ** argv, int64_t * NV, int64_t * NE, int64_t * NK, int64_t * NT,
                              int64_t * randomSeed, int64_t * iterationCount, char * initial_graph_name,
                              operation_t * operation, lbalance_t * loadBalancing, granularity_t * granularity, int64_t * edgeCount);

void CreateRandomEdgeListFromGraph(struct stinger * stingerGraph, int64_t NV, int64_t * insertionArraySrc,
                                   int64_t * insertionArrayDest, int64_t insertionCount);

void CreateRandomEdgeListFromGraphDeleting(struct stinger * stingerGraph, int64_t NV, int64_t * deletionArraySrc,
        int64_t * deletionArrayDest, int64_t deletionCount);

void ParseDimacsGraph(char initial_graph_name[1024], struct stinger ** stingerGraph, int64_t * NV, int64_t * NE);

void PickRoots(int64_t * rootArrayForApproximation, int64_t NK, int64_t NV);

float InsertEdge(int64_t count, int64_t * insertionArraySrc, int64_t * insertionArrayDest, struct stinger ** stingerGraph, 
                StreamingExtraInfo * oneSEI, extraArraysPerThread ** eAPT_perThread, int64_t * rootArrayForApproximation, int64_t NK, int64_t NV, int64_t NT, 
                bcForest * beforeBCForest, lbalance_t loadBalancing, granularity_t granularity, int64_t * iterationCount);

float DeleteEdge(int64_t count, int64_t * deletionArraySrc, int64_t * deletionArrayDest, struct stinger ** stingerGraph, 
                StreamingExtraInfo * oneSEI, extraArraysPerThread ** eAPT_perThread, int64_t * rootArrayForApproximation, int64_t NK, int64_t NV, int64_t NT, 
                bcForest * beforeBCForest, lbalance_t loadBalancing, granularity_t granularity, int64_t * iterationCount);

void CompareDynamicWithExactResults(bcForest * beforeBCForest, struct stinger * stingerGraph, int64_t * rootArrayForApproximation, int64_t NK, int64_t NT, int64_t NV);

double updateEdgeNEW(struct stinger * stingerGraph, StreamingExtraInfo * oneSEI,
                     extraArraysPerThread ** eAPT_perThread, int64_t * rootArrayForApproximation, int64_t NK,
                     int64_t NV, int64_t NT, bcForest * beforeBCForest, int64_t u_new, int64_t v_new, int64_t * iterationCount,
                     lbalance_t loadBalancing, granularity_t granularity);


double deleteEdgeNEW(struct stinger * stingerGraph, StreamingExtraInfo * oneSEI, extraArraysPerThread ** eAPT_perThread,
                     int64_t * rootArrayForApproximation, int64_t NK, int64_t NV, int64_t NT,
                     bcForest * beforeBCForest, int64_t u_old, int64_t v_old, int64_t * iterationCount, lbalance_t loadBalancing, granularity_t granularity);
