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

#define LINE_SIZE 100000

void CreateRandomEdgeListFromGraph(struct stinger * stingerGraph, int64_t NV, int64_t * insertionArraySrc,
                                   int64_t * insertionArrayDest, int64_t insertionCount)
{

    int64_t ins = 0, src, dest, srcAdj, destInAdj, destCounter;

    while (ins < insertionCount) {
        src = rand() % NV;

        srcAdj = stinger_typed_outdegree(stingerGraph, src, 0);

        // No adjacenies
        if (srcAdj == 0) {
            continue;
        }

        destInAdj = rand() % srcAdj;
        destCounter = 0;
        dest = 0;

        STINGER_FORALL_EDGES_OF_VTX_BEGIN(stingerGraph, src) {
            dest = STINGER_EDGE_DEST;

            if (destInAdj == destCounter) {
                break;
            }

            destCounter++;
        }
        STINGER_FORALL_EDGES_OF_VTX_END();

        if (src == dest) {
            continue;
        }

        stinger_remove_edge(stingerGraph, 0, src, dest);
        stinger_remove_edge(stingerGraph, 0, dest, src);

        insertionArraySrc[ins] = src;
        insertionArrayDest[ins] = dest;
        ins++;
    }
}

void CreateRandomEdgeListFromGraphDeleting(struct stinger * stingerGraph, int64_t NV, int64_t * deletionArraySrc,
        int64_t * deletionArrayDest, int64_t deletionCount)
{

    int64_t del = 0, src, dest;

    while (del < deletionCount) {
        src = rand() % NV;
        dest = rand() % NV;

        if (src == dest) {
            continue;
        }

        if (src == 0 || dest == 0) {
            continue;
        }

        int result = stinger_insert_edge(stingerGraph, 0, src, dest, 0, 0);

        if (result < 1) {
            continue;
        }

        stinger_insert_edge(stingerGraph, 0, dest, src, 0, 0);
        deletionArraySrc[del] = src;
        deletionArrayDest[del] = dest;
        del++;
    }
}


void hostParseArgsVitalUpdate(int argc, char ** argv, int64_t * NV, int64_t * NE, int64_t * NK, int64_t * NT,
                              int64_t * randomSeed, int64_t * iterationCount, char * initial_graph_name[1024],
                              operation_t * operation, lbalance_t * loadBalancing, granularity_t * granularity, int64_t * edgeCount)
{
    updateType opType;
    static struct option long_options[] = {
        {"VCount", optional_argument, 0, 'V'},
        {"Edge Count", optional_argument, 0, 'E'},
        {"RandomSeed", required_argument, 0, 'R'},
        {"Granularity", required_argument, 0, 'G'},
        {"GraphType", optional_argument, 0, 'I'},
        {"GraphName", optional_argument, 0, 'N'},
        {"OperationType", optional_argument, 0, 'O'},
        {"Threads", optional_argument, 0, 'T'},
        {"ApproximateVCount", optional_argument, 0, 'K'},
        {"LoadBalancing", optional_argument, 0, 'L'},
        {0, 0, 0, 0}
    };

    while (1) {
        int32_t option_index = 0;
        int32_t c = getopt_long(argc, argv, "V:E:R:G:I:N:O:T:K:L:?h", long_options, &option_index);
        extern char * optarg;
        extern int32_t    optind, opterr, optopt;
        int64_t intout = 0;

        if (-1 == c) {
            if (*edgeCount == -1) {
                printf("Error - Must specify number of edges to operate on.\n");
                exit(-1);
            }

            break;
        }

        switch (c) {
        default:
            printf("Unrecognized option: %c\n\n", c);

        case '?':
        case 'h':
            exit(0);
            break;

        case 'V':
            errno = 0;
            intout = strtol(optarg, NULL, 10);

            if (errno || intout < 0) {
                printf("Error - Number of vertices needs to be positive %s\n", optarg);
                exit(-1);
            }

            *NV = intout;
            break;

        case 'E':
            errno = 0;
            intout = strtol(optarg, NULL, 10);

            if (errno || intout < 0) {
                printf("Error - Number of edges needs to be positive %s\n", optarg);
                exit(-1);
            }

            *edgeCount = intout;
            break;

        case 'R':
            errno = 0;
            intout = strtol(optarg, NULL, 10);

            if (errno || intout < 0) {
                printf("Error - Seed needs to be positive. R=0 uses a random seed based on the time%s\n", optarg);
                exit(-1);
            }

            *randomSeed = intout;
            break;

        case 'I':
            errno = 0;
            intout = strtol(optarg, NULL, 10);

            if (errno || intout < 0) {
                printf("Error - Number of iterations needs to be positive %s\n", optarg);
                exit(-1);
            }

            iterationCount = intout;
            break;

        case 'G':
            errno = 0;

            if (errno || (strcmp(optarg, "COARSE") && strcmp(optarg, "FINE"))) {
                printf("Error - Granularity must be COARSE or FINE %s Invalid\n", optarg);
                exit(-1);
            }

            if (!strcmp(optarg, "COARSE")) {
                *granularity = COARSE;
            } else {
                *granularity = FINE;
            }
            break;

        case 'N':
            strcpy(initial_graph_name, optarg);
            break;

        case 'O':
            errno = 0;
            //intout = strtol(optarg, NULL, 10);

            if (errno || (strcmp(optarg, "INSERT") && strcmp(optarg, "DELETE"))) {
                printf("Error - Operation can be either INSERT or DELETE. %s Invalid\n", optarg);
                exit(-1);
            }

            if (!strcmp(optarg, "INSERT")) {
                *operation = INSERT;
            } else {
                *operation = DELETE;
            }
            //*operation = intout;
            break;

        case 'L':
            errno = 0;
            intout = strtol(optarg, NULL, 10);

            if (errno || (strcmp(optarg, "BALANCE") && strcmp(optarg, "NO_BALANCE"))) {
                printf("Error - Load balancing can either be BALANCE or NO_BALANCE. %s Invalid\n", optarg);
                exit(-1);
            }

            if (!strcmp(optarg, "BALANCE")) {
                *loadBalancing = BALANCE;
            } else {
                *loadBalancing = NO_BALANCE;
            }
            //*loadBalancing = intout;
            break;

        case 'T':
            errno = 0;
            intout = strtol(optarg, NULL, 10);

            if (errno || intout < 0) {
                printf("Error - Number of threads needs to be positive %s\n", optarg);
                exit(-1);
            }

            *NT = intout;

            if (*NT < 1) {
                *NT = 1;
            }

            break;

        case 'K':
            errno = 0;
            intout = strtol(optarg, NULL, 10);

            if (errno || intout < 0) {
                printf("Error - Number of approximate vertices needs to be positive %s\n", optarg);
                exit(-1);
            }

            *NK = intout;
            break;
        }
    }
}
//------------------------------------------
//------------------------------------------
// CSR Graph representation.
// Used only for graph ingestion.

// NE number of undirected edges.
csrGraph * CreateCSRFromStinger(struct stinger * stingerGraph, int64_t NV, int64_t NE)
{
    csrGraph * newGraph = (csrGraph *)malloc(sizeof(csrGraph));
    newGraph->NV = NV;
    newGraph->NE = NE * 2;

    newGraph->vertexPointerArray = (int64_t *)malloc(sizeof(int64_t) * (NV + 1));
    newGraph->edgeArray = (int64_t *)malloc(sizeof(int64_t) * NE);

    int64_t edgeCounter = 0;
    newGraph->vertexPointerArray[0] = 0;

    for (int64_t v = 0; v < NV; v++) {
        newGraph->vertexPointerArray[v + 1] = newGraph->vertexPointerArray[v] + stinger_outdegree(stingerGraph, v);

        STINGER_FORALL_EDGES_OF_VTX_BEGIN(stingerGraph, v) {
            newGraph->edgeArray[edgeCounter] = STINGER_EDGE_DEST;
            edgeCounter++;
        }
        STINGER_FORALL_EDGES_OF_VTX_END();
    }

    return newGraph;
}

// NE number of undirected edges.
void CreateStingerFromCSR(csrGraph * csr, struct stinger ** stingerGraph)
{
    *stingerGraph = stinger_new();

    int64_t * offset = (int64_t *)malloc(sizeof(int64_t) * (csr->NV + 1));
    int64_t * edges = (int64_t *)malloc(sizeof(int64_t) * (csr->NE));
    int64_t * weight = (int64_t *)malloc(sizeof(int64_t) * (csr->NE));

    for (int64_t v = 0; v <= csr->NV; v++) {
        offset[v] = csr->vertexPointerArray[v];
    }

    for (int64_t e = 0; e < csr->NE; e++) {
        edges[e] = csr->edgeArray[e];
        weight[e] = 0;
    }


    stinger_set_initial_edges (*stingerGraph /* G */,
                               csr->NV,
                               0,
                               offset,
                               edges,
                               weight,
                               NULL,
                               NULL,
                               0);

    free(offset);
    free(edges);
    free(weight);
    return;
}



// NE number of undirected edges.
csrGraph * CreateEmptyUndirectedCSR(int64_t NV, int64_t NE)
{

    csrGraph * newGraph = (csrGraph *)malloc(sizeof(csrGraph));
    newGraph->NV = NV;
    newGraph->NE = NE * 2;

    newGraph->vertexPointerArray = (int64_t *)malloc(sizeof(int64_t) * (newGraph->NV + 1));
    newGraph->edgeArray = (int64_t *)malloc(sizeof(int64_t) * (newGraph->NE));

    newGraph->vertexPointerArray[0] = 0;

    for (int64_t v = 0; v < NV; v++) {
        newGraph->vertexPointerArray[v + 1] = 0;
    }

    return newGraph;
}

csrGraph * CreateEmptyDirectedCSR(int64_t NV, int64_t NE)
{

    csrGraph * newGraph = (csrGraph *)malloc(sizeof(csrGraph));
    newGraph->NV = NV;
    newGraph->NE = NE;

    newGraph->vertexPointerArray = (int64_t *)malloc(sizeof(int64_t) * (newGraph->NV + 1));
    newGraph->edgeArray = (int64_t *)malloc(sizeof(int64_t) * (newGraph->NE));

    newGraph->vertexPointerArray[0] = 0;

    for (int64_t v = 0; v < NV; v++) {
        newGraph->vertexPointerArray[v + 1] = 0;
    }

    return newGraph;
}

/* Takes a path to a dimacs .graph file and parses it into a STINGER */
void ParseDimacsGraph(char initial_graph_name[1024], struct stinger ** stingerGraph, int64_t * NV_main, int64_t * NE_main)
{
    FILE * fp = fopen(initial_graph_name, "r");

    char line[LINE_SIZE];

    int64_t NV, NE;
    fgets(line, LINE_SIZE, fp);
    sscanf(line, "%ld %ld", &NV, &NE);

    NV = NV + 1;
    NE = NE * 2;

    csrGraph * CSRGraph = CreateEmptyDirectedCSR(NV, NE);

    CSRGraph->vertexPointerArray[0] = 0;
    CSRGraph->vertexPointerArray[1] = 0;
    int64_t counter = 0;

    for (int64_t u = 1; fgets(line, LINE_SIZE, fp); u++) {
        uint64_t neigh = 0;
        uint64_t v = 0;
        char * ptr = line;
        int read = 0;

        while (sscanf(ptr, "%ld%n", &v, &read) > 0) {
            ptr += read;
            neigh++;
            CSRGraph->edgeArray[counter++] = v;
        }

        CSRGraph->vertexPointerArray[u + 1] = CSRGraph->vertexPointerArray[u] + neigh;
    }

    CreateStingerFromCSR(CSRGraph, stingerGraph);
    fclose(fp);
    FreeCSR(CSRGraph);

    *NV_main = NV;
    *NE_main= NE;
}

/* Randomly picks NK roots to use for the algorithm (static and dynamic) */
void PickRoots(uint64_t * rootArrayForApproximation, int64_t NK, int64_t NV) 
{
    if (NK == NV) {
        for (int64_t vr = 0; vr < NK; vr++) {
            rootArrayForApproximation[vr] = vr;
        }
    }
    else if (NK > NV) {
        printf("Error - Number of roots is greater than number of vertices\n");
    } 
    else {
        int64_t * flag = (uint64_t *)xmalloc(sizeof(uint64_t) * NV);

        for (int64_t vv = 0; vv < NV; vv++) {
            flag[vv] = 0;
        }

        for (int64_t vr = 0; vr < NK; vr++) {
            int tempV = rand() % NV;

            if (flag[tempV] == 1) {
                vr--;
                continue;
            }

            rootArrayForApproximation[vr] = tempV;
            flag[tempV] = 1;
        }

        free(flag);
    }
}

float InsertEdge(int64_t count, int64_t * insertionArraySrc, int64_t * insertionArrayDest, struct stinger ** stingerGraph, 
                StreamingExtraInfo * oneSEI, extraArraysPerThread ** eAPT_perThread, int64_t * rootArrayForApproximation, int64_t NK, int64_t NV, int64_t NT, 
                bcForest * beforeBCForest, lbalance_t loadBalancing, granularity_t granularity, int64_t * iterationCount)
{
    int64_t src = insertionArraySrc[count];
    int64_t dest = insertionArrayDest[count];

    stinger_insert_edge(*stingerGraph, 0, src, dest, 0, 0);
    stinger_insert_edge(*stingerGraph, 0, dest, src, 0, 0);

    return  updateEdgeNEW(*stingerGraph, oneSEI, eAPT_perThread,
                                        rootArrayForApproximation, NK, NV, NT, beforeBCForest, src,
                                        dest, iterationCount, loadBalancing, granularity);
}

float DeleteEdge(int64_t count, int64_t * deletionArraySrc, int64_t * deletionArrayDest, struct stinger ** stingerGraph, 
                StreamingExtraInfo * oneSEI, extraArraysPerThread ** eAPT_perThread, int64_t * rootArrayForApproximation, int64_t NK, int64_t NV, int64_t NT, 
                bcForest * beforeBCForest, lbalance_t loadBalancing, granularity_t granularity, int64_t * iterationCount)
{
    int64_t src = deletionArraySrc[count];
    int64_t dest = deletionArrayDest[count];

    stinger_remove_edge(*stingerGraph, 0, src, dest);
    stinger_remove_edge(*stingerGraph, 0, dest, src);

    return  deleteEdgeNEW(*stingerGraph, oneSEI, eAPT_perThread,
                                        rootArrayForApproximation, NK, NV, NT, beforeBCForest, src,
                                        dest, iterationCount, loadBalancing, granularity);
}

void CompareDynamicWithExactResults(bcForest * beforeBCForest, struct stinger * stingerGraph, int64_t * rootArrayForApproximation, int64_t NK, int64_t NT, int64_t NV)
{
    // Remember to comment out the summing phase in streamingbc.c
    bcForest * afterBCForest = NULL;

    afterBCForest = CreateForestApproximate(NV, rootArrayForApproximation, NK);
    extraArraysPerThread ** eAPT_perThreadAfter = createExtraArraysForThreads(NT, NV);
    BrandesApproxCaseParallel(afterBCForest, stingerGraph, rootArrayForApproximation, NK, eAPT_perThreadAfter, NT);

    for (int a = 0; a < beforeBCForest->NV; a++) {
        if (beforeBCForest->totalBC[a] - afterBCForest->totalBC[a] > 0.001 ||
                afterBCForest->totalBC[a] - beforeBCForest->totalBC[a] > 0.001) {
            printf("Error in computation %d, before: %lf  after: %lf\n", a, beforeBCForest->totalBC[a], afterBCForest->totalBC[a]);
        }
    }


    for (int i = 0; i < NK; i++) {
        int root = rootArrayForApproximation[i];
        bcTree * beforeTree = beforeBCForest->forest[root];
        bcTree * afterTree  = afterBCForest->forest[root];

        for (int j = 0; j < beforeBCForest->NV; j++) {
            if (beforeTree->vArr[j].edgesAbove != afterTree->vArr[j].edgesAbove) {
                printf("Error in edges above value - root: %ld, vertex: %ld, before: %ld, after: %ld\n",
                       root, j, beforeTree->vArr[j].edgesAbove, afterTree->vArr[j].edgesAbove);
            }
        }

        for (int j = 0; j < beforeBCForest->NV; j++) {
            if (beforeTree->vArr[j].edgesBelow != afterTree->vArr[j].edgesBelow) {
                printf("Error in edges below value - root: %ld, vertex: %ld, before: %ld, after %ld\n",
                       root, j, beforeTree->vArr[j].edgesBelow, afterTree->vArr[j].edgesBelow);
            }
        }

        for (int j = 0; j < beforeBCForest->NV; j++) {
            if (beforeTree->vArr[j].sigma != afterTree->vArr[j].sigma) {
                printf("Error in sigma value - root: %ld, vertex: %ld, before: %ld, after: %ld\n",
                       root, j, beforeTree->vArr[j].sigma, afterTree->vArr[j].sigma);
            }
        }

        for (int j = 0; j < beforeBCForest->NV; j++) {
            if (beforeTree->vArr[j].level != afterTree->vArr[j].level) {
                printf("Error in level value - root: %ld, vertex: %ld, before: %ld, after: %ld\n",
                       root, j, beforeTree->vArr[j].level, afterTree->vArr[j].level);
            }
        }

        for (int j = 0; j < beforeBCForest->NV; j++) {
            if ((beforeTree->vArr[j].delta - afterTree->vArr[j].delta) > 0.001 ||
                    (afterTree->vArr[j].delta - beforeTree->vArr[j].delta) > 0.001) {
                printf("Error in delta value - root: %ld, vertex: %ld, before: %lf, after: %lf\n",
                       root, j, beforeTree->vArr[j].delta, afterTree->vArr[j].delta);
            }
        }
    }

    destroyExtraArraysForThreads(eAPT_perThreadAfter, NT, NV);
}
void FreeCSR(csrGraph * graph)
{
    free(graph->vertexPointerArray);
    free(graph->edgeArray);
    free(graph);
}
// NE number of undirected edges.
double updateEdgeNEW(struct stinger * stingerGraph, StreamingExtraInfo * oneSEI,
                     extraArraysPerThread ** eAPT_perThread, uint64_t * rootArrayForApproximation, int64_t NK,
                     int64_t NV, int64_t NT, bcForest * beforeBCForest, int64_t u_new, int64_t v_new, int64_t * iterationCount,
                     lbalance_t loadBalancing, granularity_t granularity)
{

    uint64_t iterator;
    float timeFullBeforeMulti = 0, timeFullAfterMulti = 0, timeStreamMulti = 0;
    float_t avgComponents = 0;
    iterationCount = 1;
    double  timeFullBefore = 0, timeFullAfter = 0, timeStream = 0;

    tic();
    *oneSEI = insertEdgeStreamingBC(beforeBCForest, stingerGraph, u_new, v_new,
                                    rootArrayForApproximation, NK, NV, NT, eAPT_perThread, loadBalancing, granularity);
    timeStreamMulti = toc();

    return timeStreamMulti;
}

double deleteEdgeNEW(struct stinger * stingerGraph, StreamingExtraInfo * oneSEI, extraArraysPerThread ** eAPT_perThread,
                     uint64_t * rootArrayForApproximation, int64_t NK, int64_t NV, int64_t NT,
                     bcForest * beforeBCForest, int64_t u_old, int64_t v_old, int64_t * iterationCount, lbalance_t loadBalancing, granularity_t granularity)
{
    uint64_t iterator;
    float timeFullBeforeMulti = 0, timeFullAfterMulti = 0, timeStreamMulti = 0;
    float_t avgComponents = 0;

    iterationCount = 1;

    double timeFullBefore = 0, timeFullAfter = 0, timeStream = 0;

    tic();
    *oneSEI = deleteEdgeStreamingBC(beforeBCForest, stingerGraph, u_old, v_old,
                                    rootArrayForApproximation, NK, NV, NT,  eAPT_perThread, loadBalancing, granularity);
    timeStreamMulti = toc();

    return timeStreamMulti;
}

int doubleCompareforMinSort (const void * a, const void * b)
{
    if ( *(double *)a <  * (double *)b ) return 1;

    if ( *(double *)a == *(double *)b ) return 0;

    if ( *(double *)a >  *(double *)b ) return -1;
}
