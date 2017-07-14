#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <assert.h>

#include <time.h>
#include <math.h>

#include "streamingbc.h"

#include "timer.h"
#include "stinger-utils.h"
#include "stinger-traversal.h"
#include "xmalloc.h"

#include <omp.h>

int main(int argc, char * argv[])
{

    int64_t NV;
    int64_t NE;
    int64_t NK;
    int64_t NT;
    int64_t edgeCount = -1;
    int64_t randomSeed;
    char initial_graph_name[1024];
    int64_t iterationCount;
    operation_t operation = INSERT; 
    lbalance_t loadBalancing = BALANCE; 
    granularity_t granularity = COARSE; 

    hostParseArgsVitalUpdate(argc, argv, &NV, &NE, &NK, &NT, &randomSeed,
                             &iterationCount, &initial_graph_name, &operation, &loadBalancing, &granularity, &edgeCount);

    int64_t insertionArraySrc[edgeCount];
    int64_t insertionArrayDest[edgeCount];
    int64_t deletionArraySrc[edgeCount];
    int64_t deletionArrayDest[edgeCount];


    printf("initial_graph_name: %s\n", initial_graph_name);
    double timingDynamic[edgeCount];
    double timingStatic;
    double timingDynamicTotal;

    StreamingExtraInfo seiDynamic[edgeCount];
    StreamingExtraInfo seiDynamicTotal;

    int64_t staticTraverseVerticeCounter;
    int64_t staticTraverseEdgeCounter;

    int64_t dynamicTraverseVerticeCounter[edgeCount];
    int64_t dynamicTraverseEdgeCounter[edgeCount];
    int64_t dynamicTraverseVerticeCounterTotal;
    int64_t dynamicTraverseEdgeCounterTotal;
    int64_t dynamicTraverseEdgeCounterMax[edgeCount];

    struct stinger * stingerGraph;
    ParseDimacsGraph(initial_graph_name, &stingerGraph, &NV, &NE);

    if (randomSeed == 0) {
        srand(time(NULL));
    } else {
        srand(randomSeed);
    }

    uint64_t * rootArrayForApproximation = NULL;
    rootArrayForApproximation = (uint64_t *) xmalloc(sizeof(uint64_t) * NK);

    PickRoots(rootArrayForApproximation, NK, NV);

    bcForest * beforeBCForest = NULL;

    // Create batch of edges to delete or insert from the graph depending on the operation.
    if (operation == INSERT) {
        // Remove edges from the graph to re-insert and update BC values.
        CreateRandomEdgeListFromGraph(stingerGraph, NV, insertionArraySrc, insertionArrayDest, edgeCount);
    } else {
        // Insert edges into the graph to re-delete and update BC values.
        CreateRandomEdgeListFromGraphDeleting(stingerGraph, NV, deletionArraySrc, deletionArrayDest, edgeCount);
    }


    //-------Compute static BC - Parallel
    //-------START

    beforeBCForest = CreateForestApproximate(NV, rootArrayForApproximation, NK);
    extraArraysPerThread ** eAPT_perThread2 = createExtraArraysForThreads(NT, NV);

    tic();
    BrandesApproxCaseParallel(beforeBCForest, stingerGraph, rootArrayForApproximation, NK, eAPT_perThread2, NT);

    // Timing analysis.
    timingStatic = toc();
    staticTraverseVerticeCounter = eAPT_perThread2[0]->staticTraverseVerticeCounter;
    staticTraverseEdgeCounter = eAPT_perThread2[0]->staticTraverseEdgeCounter;
    destroyExtraArraysForThreads(eAPT_perThread2, NT, NV);

    //------- static computation END


    //-------Compute streaming BC
    //-------START
    StreamingExtraInfo oneSEI;
    StreamingExtraInfo globalSEI = {0, 0, 0, 0};

    extraArraysPerThread ** eAPT_perThread = createExtraArraysForThreads(NT, NV);
    timingDynamicTotal = 0.0;
    dynamicTraverseVerticeCounterTotal = 0;
    dynamicTraverseEdgeCounterTotal = 0;

    for (int64_t i = 0; i < edgeCount; i++) {
        if (operation == INSERT) {
            float time_taken = InsertEdge(i, insertionArraySrc, insertionArrayDest, &stingerGraph, &oneSEI, eAPT_perThread, 
                                    rootArrayForApproximation, NK, NV, NT, beforeBCForest, loadBalancing, granularity, &iterationCount);
            
            timingDynamic[i] = time_taken;
        } else {
            float time_taken = DeleteEdge(i, deletionArraySrc, deletionArrayDest, &stingerGraph, &oneSEI, eAPT_perThread, 
                                    rootArrayForApproximation, NK, NV, NT, beforeBCForest, loadBalancing, granularity, &iterationCount);
            
            timingDynamic[i] = time_taken;
        }

        // Timing analysis begin.
        timingDynamicTotal += timingDynamic[i];

        globalSEI.adjacent += oneSEI.adjacent;
        globalSEI.movement += oneSEI.movement;
        globalSEI.sameLevel += oneSEI.sameLevel;
        seiDynamic[i].movement = oneSEI.movement;
        seiDynamic[i].adjacent = oneSEI.adjacent;
        seiDynamic[i].sameLevel = oneSEI.sameLevel;
        oneSEI.adjacent = 0;
        oneSEI.movement = 0;
        oneSEI.sameLevel = 0;

        dynamicTraverseVerticeCounter[i] = eAPT_perThread[0]->dynamicTraverseVerticeCounter;
        dynamicTraverseEdgeCounter[i] = eAPT_perThread[0]->dynamicTraverseEdgeCounter;
        dynamicTraverseVerticeCounterTotal += dynamicTraverseVerticeCounter[i];
        dynamicTraverseEdgeCounterTotal += dynamicTraverseEdgeCounter[i];

        dynamicTraverseEdgeCounterMax[i] = eAPT_perThread[0]->dynamicTraverseEdgeCounter;

        for (int t = 0; t < NT; t++) {
            if (dynamicTraverseEdgeCounterMax[i] < eAPT_perThread[t]->dynamicTraverseEdgeCounter) {
                dynamicTraverseEdgeCounterMax[i] = eAPT_perThread[t]->dynamicTraverseEdgeCounter;
            }

            ClearCounters(eAPT_perThread[t]);
        }
        // Timing analysis end.
    }

    seiDynamicTotal.movement = globalSEI.movement;
    seiDynamicTotal.adjacent = globalSEI.adjacent;
    seiDynamicTotal.sameLevel = globalSEI.sameLevel;


    // Not necessary unless you want to check what the static algorithm would output on the final graph.
    //CompareDynamicWithExactResults(beforeBCForest, stingerGraph, rootArrayForApproximation, NK, NT, NV);
    //------- streaming computation END

    free(rootArrayForApproximation);
    stinger_free(stingerGraph);

    // Print time taken for each edge to be updated.
    for (int64_t i = 0; i < edgeCount; i++) {
        printf("%ld, %9lf\n", i, (double)(timingDynamic[i])); // Min speedup
    }
}

