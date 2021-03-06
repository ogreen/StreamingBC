#include <omp.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#include "streamingbc_aux.h"

#define PARENT_ANCHORED 3
#define SIBLING_ANCHORED 4

void addEdgeWithoutMovementBrandesFG(bcForest * forest, struct stinger * sStinger,
                                     int64_t currRoot, int64_t startVertex, int64_t parentVertex,
                                     int64_t addedPathsToRoot, extraArraysPerThread * eAPT, int64_t cores)
{

    bcTree * tree = forest->forest[currRoot];
    int64_t * QueueDown = eAPT->QueueDown;
    int64_t * QueueUp = eAPT->QueueUp;
    int64_t * QueueDownBorders = eAPT->QueueSame;

    int64_t NV = forest->NV;


    eAPT->sV[startVertex].newSigma = tree->vArr[startVertex].sigma;
    eAPT->sV[parentVertex].newEdgesBelow = tree->vArr[parentVertex].edgesBelow;
    eAPT->sV[startVertex].newEdgesBelow = tree->vArr[startVertex].edgesBelow;
    eAPT->sV[parentVertex].newEdgesAbove = tree->vArr[parentVertex].edgesAbove;
    eAPT->sV[startVertex].newEdgesAbove = tree->vArr[startVertex].edgesAbove;
    eAPT->sV[startVertex].touched = 1;
    eAPT->sV[startVertex].newSigma += addedPathsToRoot;
    eAPT->sV[startVertex].diffPath = addedPathsToRoot;
    eAPT->sV[startVertex].newEdgesAbove += eAPT->sV[parentVertex].newEdgesAbove + 1;
    eAPT->sV[parentVertex].newEdgesBelow += eAPT->sV[startVertex].newEdgesBelow + 1;

    QueueDown[0] = startVertex;

    int64_t * qStart = &(eAPT->qStart);
    int64_t * qEnd = &(eAPT->qEnd);
    int64_t * qStart_nxt = &(eAPT->qStart_nxt);
    int64_t * qEnd_nxt = &(eAPT->qEnd_nxt);

    int64_t qDownBIndex = 0;
    *qEnd = 1;
    *qStart_nxt = 1;
    *qEnd_nxt = 1;

    int64_t deepestLevel = tree->vArr[startVertex].level;
    int64_t intialLevel = tree->vArr[startVertex].level;

    int64_t qDownEndMarker = -1;
    int64_t qDownEnd = -1;

    // Starting BFS decent from "startVertex", down to all the vertices that have shortest paths through "startVertex"
    // All elements that will be touched will receive a positive value in their touched field.
    // In this implementation, "STACKS" are not used for the "moving up" stage. Rather, a multi-level queue is used.
    // Each level in the tree(max depth NV) has a queue and a counter specifiying how deep a specific deepth-queue is.
    // For simplicity, all elements are pushed both into the multi-level queue and into the regular queue which is used
    // for the BFS traversal.
    #pragma omp parallel num_threads(cores)
    {
        while (*qStart < *qEnd) {
            #pragma omp master
            {
                QueueDownBorders[qDownBIndex++] = *qStart;
                QueueDownBorders[qDownBIndex++] = *qEnd;
            }
            #pragma omp barrier

            #pragma omp for

            for (int64_t i = *qStart; i < *qEnd; i++) {
                int64_t currElement = QueueDown[i];
                int64_t levelCurrPlusOne = tree->vArr[currElement].level + 1;
                int64_t touchedCurrPlusOne = eAPT->sV[currElement].touched + 1;

                STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement) {
                    int64_t k = STINGER_EDGE_DEST;

                    // if this vertex has not been added yet
                    if (levelCurrPlusOne == (tree->vArr[k].level)) {

                        if (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 0, currElement)) {
                            // Checking if a "deeper level" has been reached.
                            if (deepestLevel < tree->vArr[k].level) {
                                deepestLevel = tree->vArr[k].level;
                            }

                            __atomic_fetch_add(&(eAPT->sV[k].newEdgesAbove), tree->vArr[k].edgesAbove - tree->vArr[currElement].edgesAbove +
                                               eAPT->sV[currElement].newEdgesAbove, __ATOMIC_RELAXED);

                            __atomic_fetch_add(&(eAPT->sV[k].newSigma), tree->vArr[k].sigma, __ATOMIC_RELAXED);

                            // insert this vertex into the BFS queue
                            QueueDown[__atomic_fetch_add(qEnd_nxt, 1, __ATOMIC_RELAXED)] = k;
                            // indicate that it is in the next level of the BFS
                            // add new paths to root that go through current BFS Vertex
                            __atomic_fetch_add(&(eAPT->sV[k].newSigma), eAPT->sV[currElement].diffPath, __ATOMIC_RELAXED);
                            // pass on my new paths to root for its search
                            __atomic_fetch_add(&(eAPT->sV[k].diffPath), eAPT->sV[currElement].diffPath, __ATOMIC_RELAXED);
                        }
                        // otherwise if it has been touched, but is specifically in the next level
                        // of the search (meaning it has more than one edge to the current level)
                        else if (eAPT->sV[k].touched != currElement) {
                            __atomic_fetch_add(&(eAPT->sV[k].newEdgesAbove), -tree->vArr[currElement].edgesAbove, __ATOMIC_RELAXED);
                            __atomic_fetch_add(&(eAPT->sV[k].newEdgesAbove), eAPT->sV[currElement].newEdgesAbove, __ATOMIC_RELAXED);

                            // add new paths to root that go through current BFS Vertex
                            __atomic_fetch_add(&(eAPT->sV[k].newSigma), eAPT->sV[currElement].diffPath, __ATOMIC_RELAXED);
                            // pass on my new paths to root for its search
                            __atomic_fetch_add(&(eAPT->sV[k].diffPath), eAPT->sV[currElement].diffPath, __ATOMIC_RELAXED);
                        }
                    }
                }
                STINGER_FORALL_EDGES_OF_VTX_END();

#if COUNT_TRAVERSALS==1
                eAPT->dynamicTraverseEdgeCounter += stinger_typed_outdegree(sStinger, currElement, 0);
                eAPT->dynamicTraverseVerticeCounter++;
#endif
            }

            #pragma omp master
            {
                *qStart = *qStart_nxt;
                *qEnd = *qEnd_nxt;
                *qStart_nxt = *qEnd;
                *qEnd_nxt = *qStart_nxt;
            }
            #pragma omp barrier
        }

        #pragma omp master
        {
            qDownEndMarker = *qEnd - 1;
        }
        #pragma omp barrier

        #pragma omp master
        {
            *qStart = 0;
            *qEnd = 0;
            *qStart_nxt = 0;
            *qEnd_nxt = 0;
        }
        #pragma omp barrier

        // Starting Multi-Level "BFS" ascent.
        // The ascent continues going up as long as the root has not been reached and that there
        // are elements in the current level of the ascent. The ascent starts in the deepest level
        // of the graph.
        // It was worth noting that in the ascent stage:
        // 1) All previously untouched elements that are touched are marked with "-1".
        // 2) On the way up, it is possible that elements that were touched in the BFS decent, will
        // touch elements that were not touchded in the decsent and that are below "vertex". These
        // are elements that do not have shortest paths going through "vertex" ,however, there BC
        // values have changed due to the changes occuring below them. Because of this, they are
        // placed in the Multi-level queue.
        //#pragma omp parallel num_threads(cores)
        //{
        while (!(qDownBIndex <= 0 && *qStart >= *qEnd && *qStart_nxt >= *qEnd_nxt)) {
            if (qDownBIndex >= 2) {

                #pragma omp for
                for (int64_t i = QueueDownBorders[qDownBIndex - 2]; i < QueueDownBorders[qDownBIndex - 1]; i++) {

                    int64_t currElement = QueueDown[i];
                    int64_t levelCurrMinusOne = tree->vArr[currElement].level - 1;
                    STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement) {
                        int64_t k = STINGER_EDGE_DEST;

                        if (tree->vArr[k].level == levelCurrMinusOne) {
                            // Checking to see if "k" has been touched before.
                            if (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 0, -1)) {
                                eAPT->sV[k].newDelta += tree->vArr[k].delta;
                                // Marking element as touched in the ascent stage.
                                QueueUp[__atomic_fetch_add(qEnd_nxt, 1, __ATOMIC_RELAXED)] = k;
                                __atomic_fetch_add(&(eAPT->sV[k].newSigma), tree->vArr[k].sigma, __ATOMIC_RELAXED);

                                if (k != parentVertex) {
                                    __atomic_fetch_add(&(eAPT->sV[k].newEdgesBelow), tree->vArr[k].edgesBelow, __ATOMIC_RELAXED);
                                }
                            }

                            if (k != parentVertex && tree->vArr[k].level <= tree->vArr[parentVertex].level) {
                                __atomic_fetch_add(&(eAPT->sV[k].newEdgesBelow), -tree->vArr[currElement].edgesBelow, __ATOMIC_RELAXED);
                                __atomic_fetch_add(&(eAPT->sV[k].newEdgesBelow), eAPT->sV[currElement].newEdgesBelow, __ATOMIC_RELAXED);
                            }
                        }


                        if (tree->vArr[k].level == tree->vArr[currElement].level + 1 && eAPT->sV[k].touched != 0) {
                            eAPT->sV[currElement].newDelta +=
                                ((bc_t)eAPT->sV[currElement].newSigma / (bc_t)eAPT->sV[k].newSigma) *
                                (bc_t)(eAPT->sV[k].newDelta + 1);


                            // For the elements that are touched in the ascent stage it is necessary to
                            // to reduce the values that they previously had.
                            // In addition to this, the "parentVertex" that is connected to "vertex", i.e.
                            // the vertices of the new edge, needs to increase its betweenness centrality
                            // following the new connection, without removing the old delta value.
                            if (eAPT->sV[currElement].touched < 0 && (currElement != parentVertex || k != startVertex)) {
                                eAPT->sV[currElement].newDelta -=
                                    ((bc_t)tree->vArr[currElement].sigma / (bc_t)tree->vArr[k].sigma) *
                                    (bc_t)(tree->vArr[k].delta + 1);
                            }
                        }
                    }
                    STINGER_FORALL_EDGES_OF_VTX_END();

#if COUNT_TRAVERSALS==1
                    eAPT->dynamicTraverseEdgeCounter += stinger_typed_outdegree(sStinger, currElement, 0);
                    eAPT->dynamicTraverseVerticeCounter++;
#endif

                    if (currElement != currRoot) {
                        eAPT->sV[currElement].totalBC += eAPT->sV[currElement].newDelta - tree->vArr[currElement].delta;
                    }
                }
            }

            #pragma omp master
            {
                qDownBIndex -= 2;
            }
            #pragma omp barrier

            #pragma omp master
            {
                *qStart = *qStart_nxt;
                *qEnd = *qEnd_nxt;
                *qStart_nxt = *qEnd;
                *qEnd_nxt = *qStart_nxt;
            }
            #pragma omp barrier

            #pragma omp for

            for (int64_t i = *qStart; i < *qEnd; i++) {
                int64_t currElement = QueueUp[i];
                int64_t levelCurrMinusOne = tree->vArr[currElement].level - 1;
                STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement) {
                    int64_t k = STINGER_EDGE_DEST;

                    if (tree->vArr[k].level == levelCurrMinusOne) {
                        // Checking to see if "k" has been touched before.
                        if (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 0, -1)) {
                            eAPT->sV[k].newDelta += tree->vArr[k].delta;
                            // Marking element as touched in the ascent stage.
                            QueueUp[__atomic_fetch_add(qEnd_nxt, 1, __ATOMIC_RELAXED)] = k;
                            __atomic_fetch_add(&(eAPT->sV[k].newSigma), tree->vArr[k].sigma, __ATOMIC_RELAXED);

                            if (k != parentVertex) {
                                __atomic_fetch_add(&(eAPT->sV[k].newEdgesBelow), tree->vArr[k].edgesBelow, __ATOMIC_RELAXED);
                            }
                        }

                        if (k != parentVertex && tree->vArr[k].level <= tree->vArr[parentVertex].level) {
                            __atomic_fetch_add(&(eAPT->sV[k].newEdgesBelow), -tree->vArr[currElement].edgesBelow, __ATOMIC_RELAXED);
                            __atomic_fetch_add(&(eAPT->sV[k].newEdgesBelow), eAPT->sV[currElement].newEdgesBelow, __ATOMIC_RELAXED);
                        }

                    }

                    if (tree->vArr[k].level == tree->vArr[currElement].level + 1 && eAPT->sV[k].touched != 0) {
                        eAPT->sV[currElement].newDelta +=
                            ((bc_t)eAPT->sV[currElement].newSigma / (bc_t)eAPT->sV[k].newSigma) *
                            (bc_t)(eAPT->sV[k].newDelta + 1);

                        // For the elements that are touched in the ascent stage it is necessary to
                        // to reduce the values that they previously had.
                        // In addition to this, the "parentVertex" that is connected to "vertex", i.e.
                        // the vertices of the new edge, needs to increase its betweenness centrality
                        // following the new connection, without removing the old delta value.
                        if (eAPT->sV[currElement].touched < 0 && (currElement != parentVertex || k != startVertex)) {
                            eAPT->sV[currElement].newDelta -=
                                ((bc_t)tree->vArr[currElement].sigma / (bc_t)tree->vArr[k].sigma) *
                                (bc_t)(tree->vArr[k].delta + 1);
                        }
                    }
                }
                STINGER_FORALL_EDGES_OF_VTX_END();

#if COUNT_TRAVERSALS==1
                eAPT->dynamicTraverseEdgeCounter += stinger_typed_outdegree(sStinger, currElement, 0);
                eAPT->dynamicTraverseVerticeCounter++;
#endif

                if (currElement != currRoot) {
                    eAPT->sV[currElement].totalBC += eAPT->sV[currElement].newDelta - tree->vArr[currElement].delta;
                }
            }

            #pragma omp master
            {
                *qStart = *qEnd;
            }
            #pragma omp barrier
        }

        #pragma omp for

        for (int64_t c = 0; c <= qDownEndMarker; c++) {
            int64_t k = QueueDown[c];
            tree->vArr[k].delta = eAPT->sV[k].newDelta;
            tree->vArr[k].sigma = eAPT->sV[k].newSigma;
            tree->vArr[k].edgesAbove = eAPT->sV[k].newEdgesAbove;
            eAPT->sV[k].diffPath = 0;
            eAPT->sV[k].touched = 0;
            eAPT->sV[k].newDelta = 0.0;
            eAPT->sV[k].newSigma = 0;
            eAPT->sV[k].newEdgesAbove = 0;
        }

        eAPT->sV[startVertex].newEdgesAbove = 0;
        eAPT->sV[parentVertex].newEdgesAbove = 0;

        #pragma omp for

        for (int64_t c = 0; c < *qEnd; c++) {
            int64_t k = QueueUp[c];
            tree->vArr[k].delta = eAPT->sV[k].newDelta;
            tree->vArr[k].edgesBelow = eAPT->sV[k].newEdgesBelow;
            eAPT->sV[k].diffPath = 0;
            eAPT->sV[k].touched = 0;
            eAPT->sV[k].newDelta = 0.0;
            eAPT->sV[k].newEdgesBelow = 0;
            eAPT->sV[k].newSigma = 0;
            eAPT->sV[k].newEdgesBelow = 0;
        }

    } 
    eAPT->sV[parentVertex].newEdgesBelow = 0;
    eAPT->sV[startVertex].newEdgesBelow = 0;

    eAPT->qStart = 0;
    eAPT->qEnd = 0;
    eAPT->qStart_nxt = 0;
    eAPT->qEnd_nxt = 0;
}


void moveUpTreeBrandesFG(bcForest * forest, struct stinger * sStinger,
                         int64_t currRoot, int64_t startVertex, int64_t parentVertex,
                         int64_t prevDist, extraArraysPerThread * eAPT, int64_t cores)
{

    bcTree * tree = forest->forest[currRoot];
    int64_t NV = forest->NV;
    int64_t * QueueDown = eAPT->QueueDown;
    int64_t * QueueUp = eAPT->QueueUp;
    int64_t * QueueSame = eAPT->QueueSame;
    int64_t * QueueDownBorders = eAPT->Stack;

    list_ptr * multiLevelQueues = eAPT->multiLevelQueues;
    queue_t * queue = eAPT->queue;
    level_node_t * levelIndices = eAPT->levelIndices;

    eAPT->sV[parentVertex].newSigma = tree->vArr[parentVertex].sigma;
    eAPT->sV[startVertex].newSigma = tree->vArr[startVertex].sigma;

    int64_t * qStart = &(eAPT->qStart);
    int64_t * qEnd = &(eAPT->qEnd);
    int64_t * qStart_nxt = &(eAPT->qStart_nxt);
    int64_t * qEnd_nxt = &(eAPT->qEnd_nxt);
    int64_t qDownBIndex = 0;

    int64_t * qStartSame = &(eAPT->qStartSame);
    int64_t * qEndSame = &(eAPT->qEndSame);
    int64_t * qStartSame_nxt = &(eAPT->qStartSame_nxt);
    int64_t * qEndSame_nxt = &(eAPT->qEndSame_nxt);

    *qEnd = 1;
    *qStart_nxt = 1;
    *qEnd_nxt = 1;


    int64_t qDownEndMarker = -1;
    int64_t depthDown = -1, depthUp = -1, depthSame = -1;
    int64_t upCounter = -1;
    int64_t currElement = 0; //dummy initilization - variable will be initialized in function.
    int operation = -1; // 0 - down, 1 - up, 2 - same for dependency accumulation.

    QueueDown[0] = startVertex;
    eAPT->sV[startVertex].touched = 1;
    eAPT->sV[startVertex].newSigma = eAPT->sV[parentVertex].newSigma;
    eAPT->sV[startVertex].diffPath = eAPT->sV[parentVertex].newSigma;

    eAPT->sV[startVertex].movementDelta = prevDist;
    eAPT->sV[startVertex].IMoved = 1;

    eAPT->sV[parentVertex].newEdgesAbove = tree->vArr[parentVertex].edgesAbove;
    eAPT->sV[startVertex].newEdgesAbove = eAPT->sV[parentVertex].newEdgesAbove + 1;
    int64_t deepestLevel = tree->vArr[parentVertex].level + 1;

    // Starting BFS decent from "startVertex", down to all the vertices that have shortest paths through "startVertex"
    // All elements that will be touched will receive a positive value in their touched field.
    // In this implementation, "STACKS" are not used for the "moving up" stage. Rather, a multi-level queue is used.
    // Each level in the tree(max depth NV) has a queue and a counter specifiying how deep a specific deepth-queue is.
    // For simplicity, all elements are pushed both into the multi-level queue and into the regular queue which is used
    // for the BFS traversal.
    #pragma omp parallel num_threads(cores)
    {
        while (*qStart < *qEnd) {
            #pragma omp master
            {
                QueueDownBorders[qDownBIndex++] = *qStart;
                QueueDownBorders[qDownBIndex++] = *qEnd;
            }
            #pragma omp barrier


            #pragma omp for

            for (int64_t i = *qStart; i < *qEnd; i++) {

                int64_t currElement = QueueDown[i];
                int64_t touchedCurrPlusOne = eAPT->sV[currElement].touched + 1;
                __atomic_fetch_add(&(eAPT->sV[currElement].newEdgesAbove), -eAPT->sV[currElement].newEdgesAbove, __ATOMIC_RELAXED);

                STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement) {
                    int64_t k = STINGER_EDGE_DEST;

                    int64_t computedDelta = eAPT->sV[currElement].movementDelta -
                                            (tree->vArr[currElement].level - tree->vArr[k].level + 1);

                    int64_t newCurrLevel = 0;
                    __atomic_fetch_add(&newCurrLevel, tree->vArr[currElement].level, __ATOMIC_RELAXED);
                    __atomic_fetch_add(&newCurrLevel, -eAPT->sV[currElement].movementDelta, __ATOMIC_RELAXED);

                    int64_t newKLevel = 0;
                    __atomic_fetch_add(&newKLevel, tree->vArr[k].level, __ATOMIC_RELAXED);
                    __atomic_fetch_add(&newKLevel, -computedDelta, __ATOMIC_RELAXED);

                    if (computedDelta < 0 && eAPT->sV[k].touched == 0) {
                        if (computedDelta >= 0 && newKLevel < newCurrLevel) {
                            __atomic_fetch_add(&(eAPT->sV[k].newEdgesAbove), tree->vArr[k].edgesAbove, __ATOMIC_RELAXED);
                            __atomic_fetch_add(&(eAPT->sV[currElement].newEdgesAbove), eAPT->sV[k].newEdgesAbove + 1, __ATOMIC_RELAXED);
                        }

                        if (computedDelta < 0 && tree->vArr[k].level < newCurrLevel) {
                            __atomic_fetch_add(&(eAPT->sV[currElement].newEdgesAbove), tree->vArr[k].edgesAbove + 1, __ATOMIC_RELAXED);
                        }
                    }
                    else if (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 0, touchedCurrPlusOne)) {
                        if (computedDelta >= 0 && newKLevel < newCurrLevel) {
                            __atomic_fetch_add(&(eAPT->sV[k].newEdgesAbove), tree->vArr[k].edgesAbove, __ATOMIC_RELAXED);
                            __atomic_fetch_add(&(eAPT->sV[currElement].newEdgesAbove), eAPT->sV[k].newEdgesAbove + 1, __ATOMIC_RELAXED);
                        }

                        if (computedDelta < 0 && tree->vArr[k].level < newCurrLevel) {
                            __atomic_fetch_add(&(eAPT->sV[currElement].newEdgesAbove), tree->vArr[k].edgesAbove + 1, __ATOMIC_RELAXED);
                        }

                        // if the adjacent vertex should be moved, put it in the queue
                        if (computedDelta > 0) {
                            __atomic_fetch_add(&(eAPT->sV[k].newSigma), eAPT->sV[currElement].diffPath, __ATOMIC_RELAXED);
                            __atomic_fetch_add(&(eAPT->sV[k].diffPath), eAPT->sV[currElement].diffPath, __ATOMIC_RELAXED);
                            __atomic_fetch_add(&(eAPT->sV[k].movementDelta), computedDelta, __ATOMIC_RELAXED);
                            __atomic_fetch_add(&(eAPT->sV[k].IMoved), 2, __ATOMIC_RELAXED);
                            QueueDown[__atomic_fetch_add(qEnd_nxt, 1, __ATOMIC_RELAXED)] = k;
                        }
                        // Vertex that will not be moved has been found.
                        else if (computedDelta == 0) {
                            __atomic_fetch_add(&(eAPT->sV[k].newSigma), tree->vArr[k].sigma, __ATOMIC_RELAXED);
                            __atomic_fetch_add(&(eAPT->sV[k].newSigma), eAPT->sV[currElement].diffPath, __ATOMIC_RELAXED);
                            __atomic_fetch_add(&(eAPT->sV[k].diffPath), eAPT->sV[currElement].diffPath, __ATOMIC_RELAXED);
                            __atomic_fetch_add(&(eAPT->sV[k].movementDelta), computedDelta, __ATOMIC_RELAXED);
                            __atomic_fetch_add(&(eAPT->sV[k].IMoved), -eAPT->sV[k].IMoved, __ATOMIC_RELAXED);
                            QueueDown[__atomic_fetch_add(qEnd_nxt, 1, __ATOMIC_RELAXED)] = k;
                        }

                        // Vertex that the number of shortest path to the root does not change has been found.
                        // This vertex is not marked as it might be touched on the way up.

                        // if adjacent and in the next level
                    }
                    else if (eAPT->sV[k].touched == touchedCurrPlusOne) {
                        if (computedDelta >= 0 && newKLevel < newCurrLevel) {
                            __atomic_fetch_add(&(eAPT->sV[currElement].newEdgesAbove), eAPT->sV[k].newEdgesAbove + 1, __ATOMIC_RELAXED);
                        }

                        if (computedDelta < 0 && tree->vArr[k].level < newCurrLevel) {
                            __atomic_fetch_add(&(eAPT->sV[currElement].newEdgesAbove), eAPT->sV[k].newEdgesAbove + 1, __ATOMIC_RELAXED);
                        }

                        if (computedDelta >= 0) {
                            __atomic_fetch_add(&(eAPT->sV[k].newSigma), eAPT->sV[currElement].diffPath, __ATOMIC_RELAXED);
                            __atomic_fetch_add(&(eAPT->sV[k].diffPath), eAPT->sV[currElement].diffPath, __ATOMIC_RELAXED);
                        }


                    } else if (computedDelta >= 0 && newKLevel < newCurrLevel) {
                        __atomic_fetch_add(&(eAPT->sV[currElement].newEdgesAbove), eAPT->sV[k].newEdgesAbove + 1, __ATOMIC_RELAXED);
                    } else if (computedDelta < 0 && tree->vArr[k].level < newCurrLevel) {
                        __atomic_fetch_add(&(eAPT->sV[currElement].newEdgesAbove), eAPT->sV[k].newEdgesAbove + 1, __ATOMIC_RELAXED);
                    }
                }
                STINGER_FORALL_EDGES_OF_VTX_END();


#if COUNT_TRAVERSALS==1
                eAPT->dynamicTraverseEdgeCounter += stinger_typed_outdegree(sStinger, currElement, 0);
                eAPT->dynamicTraverseVerticeCounter++;
#endif
                // move ourself and retire
                __atomic_fetch_add(&(tree->vArr[currElement].level), -eAPT->sV[currElement].movementDelta, __ATOMIC_RELAXED);
                appendDS2(queue, levelIndices, tree->vArr[currElement].level, currElement, omp_get_thread_num());

                // Checking if a "deeper level" has been reached.
                if (deepestLevel < tree->vArr[currElement].level) {
                    deepestLevel = tree->vArr[currElement].level;
                }
            }

            #pragma omp master
            {
                *qStart = *qStart_nxt;
                *qEnd = *qEnd_nxt;
                *qStart_nxt = *qEnd;
                *qEnd_nxt = *qStart_nxt;
            }
            #pragma omp barrier
        }

        // Starting Multi-Level "BFS" ascent.
        #pragma omp master
        {
            *qEnd = 0;
            queue_node_t * temp_node;

            for (int lev = tree->vArr[startVertex].level; lev < NV; lev++) {
                temp_node = getFirstDS(queue, levelIndices, lev);

                while (temp_node != NULL) {
                    QueueDown[(*qEnd)++] = temp_node->data;
                    deleteFirstDS(queue, levelIndices, lev);
                    temp_node = getFirstDS(queue, levelIndices, lev);
                }
            }
        }
        #pragma omp barrier

        #pragma omp master
        {
            (*qEnd)--;
            qDownEndMarker = *qEnd;
        }
        #pragma omp barrier

        int64_t QUpStart = 0, QUpEnd = 0;
        int64_t QSameStart = 0, QSameEnd = 0;
        currElement = 0; //dummy initilization - variable will be initialized in function.
        upCounter = 0;
        depthDown = tree->vArr[QueueDown[*qEnd]].level, depthUp = -1, depthSame = -1;

        *qStart = 0;
        *qEnd = 0;
        *qStart_nxt = 0;
        *qEnd_nxt = 0;

        *qStartSame = 0;
        *qEndSame = 0;
        *qStartSame_nxt = 0;
        *qEndSame_nxt = 0;

        // The ascent continues going up as long as the root has not been reached and that there
        // are elements in the current level of the ascent. The ascent starts in the deepest level
        // of the graph.
        // It was worth noting that in the ascent stage:
        // 1) All previously untouched elements that are touched are marked with "-1".
        // 2) On the way up, it is possible that elements that were touched in the BFS decent, will
        // touch elements that were not touchded in the decsent and that are below "startVertex". These
        // are elements that do not have shortest paths going through "startVertex" ,however, there BC
        // values have changed due to the changes occuring below them. Because of this, they are
        // placed in the Multi-level queue.
        // 3) There are vertices that did not move and that one of their neighbors move up(such that
        // the vertices are now in the same level). Consequently, the number of shortest path going
        // through the vertex that did not move was reduced. These vertices will be touched as -2
        // and added to the queue and the "BFS ascent" will continue from these vertices as well.

        while (!(qDownBIndex <= 0 && *qStart >= *qEnd && *qStart_nxt >= *qEnd_nxt && *qStartSame >= *qEndSame && *qStartSame_nxt >= *qEndSame_nxt)
                && !(depthDown == -1 && depthSame == -1 && depthUp == -1)) {

            #pragma omp master
            {
                operation = -1; // 0 - down, 1 - up, 2 - same
            }
            #pragma omp barrier

            #pragma omp master
            {
                if (depthUp >= depthSame && depthUp >= depthDown)
                {
                    operation = 1;

                    if (*qEnd_nxt > *qStart_nxt)
                        depthUp = -1;
                    else {
                        depthUp = tree->vArr[QueueUp[*qStart_nxt]].level;
                    }
                } else if (depthDown >= depthSame && depthDown >= depthUp)
                {
                    operation = 0;

                    if (qDownBIndex < 2 || QueueDownBorders[qDownBIndex - 2] > QueueDownBorders[qDownBIndex - 1])
                        depthDown = -1;
                    else if (qDownBIndex > 2) {
                        depthDown = tree->vArr[QueueDown[QueueDownBorders[qDownBIndex - 2] - 1]].level;
                    }
                }

                else if (depthDown <= depthSame && depthUp <= depthSame)
                {
                    operation = 2;

                    if (*qEndSame_nxt > *qStartSame_nxt)
                        depthSame = -1;
                    else
                        depthSame = tree->vArr[QueueSame[*qStartSame_nxt]].level;
                }
            }
            #pragma omp barrier

            if (operation == 0 && qDownBIndex >= 2) {

                #pragma omp for
                for (int64_t i = QueueDownBorders[qDownBIndex - 2]; i < QueueDownBorders[qDownBIndex - 1]; i++) {

                    int64_t currElement = QueueDown[i];
                    int64_t levelCurrMinusOne = tree->vArr[currElement].level - 1;
                    eAPT->sV[currElement].newEdgesBelow = 0;
                    STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement) {
                        int64_t k = STINGER_EDGE_DEST;
                        // Checking that the vertices are in different levels.

                        if (tree->vArr[k].level == tree->vArr[currElement].level + 1) {

                            if (eAPT->sV[k].touched == 0) {
                                eAPT->sV[currElement].newEdgesBelow += tree->vArr[k].edgesBelow + 1;
                            } else {
                                eAPT->sV[currElement].newEdgesBelow += eAPT->sV[k].newEdgesBelow + 1;
                            }
                        }

                        if (tree->vArr[k].level == levelCurrMinusOne) {
                            // Checking to see if "k" has been touched before.
                            if (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 0, -1)) {
                                eAPT->sV[k].newDelta += tree->vArr[k].delta;
                                upCounter++;
                                // Marking element as touched in the ascent stage.
                                eAPT->sV[k].touched = -1;

                                __sync_bool_compare_and_swap(&depthUp, -1, tree->vArr[k].level);
                                QueueUp[__atomic_fetch_add(qEnd_nxt, 1, __ATOMIC_RELAXED)] = k;

                                if (k != parentVertex)
                                    eAPT->sV[k].newSigma += tree->vArr[k].sigma;
                            }
                        }

                        if (tree->vArr[k].level == tree->vArr[currElement].level + 1 && eAPT->sV[k].touched != 0) {
                            eAPT->sV[currElement].newDelta +=
                                ((bc_t)eAPT->sV[currElement].newSigma / (bc_t)eAPT->sV[k].newSigma) *
                                (bc_t)(eAPT->sV[k].newDelta + 1);

                            // For the elements that are touched in the ascent stage it is necessary to
                            // to reduce the values that they previously had.
                            // In addition to this, the "parentVertex" that is connected to "startVertex", i.e.
                            // the vertices of the new edge, needs to increase its betweenness centrality
                            // following the new connection, without removing the old delta value.
                            if (eAPT->sV[currElement].touched < 0 && ( currElement != parentVertex || k != startVertex)) {
                                eAPT->sV[currElement].newDelta -=
                                    ((bc_t)tree->vArr[currElement].sigma / (bc_t)tree->vArr[k].sigma) *
                                    (bc_t)(tree->vArr[k].delta + 1);

                            }
                        }


                        // Vertices that did not move and that one of their neighbors move up(such that
                        // the vertices are now in the same level).
                        if (tree->vArr[k].level == tree->vArr[currElement].level && ((eAPT->sV[currElement].IMoved == 1 && eAPT->sV[k].IMoved < 0) )) {
                            // Checking to see if "k" has been touched before.
                            if (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 0, -1)) {
                                eAPT->sV[k].newDelta += tree->vArr[k].delta;

                                upCounter++;
                                // Marking element as touched in the ascent stage.
                                eAPT->sV[k].touched = -2;
                                __sync_bool_compare_and_swap(&depthSame, -1, tree->vArr[k].level);
                                QueueSame[__atomic_fetch_add(qEndSame_nxt, 1, __ATOMIC_RELAXED)] = k;
                                eAPT->sV[k].newSigma += tree->vArr[k].sigma;
                            }

                            // Paths that previosul went through this vertex no longer go through them, thus the
                            // shortest path count(BC) is reduced.

                        }

                        if (tree->vArr[k].level == tree->vArr[currElement].level && ((eAPT->sV[k].IMoved == 1 && eAPT->sV[currElement].IMoved < 0) )) {
                            eAPT->sV[currElement].newDelta -=
                                ((bc_t)tree->vArr[currElement].sigma / (bc_t)tree->vArr[k].sigma) *
                                (bc_t)(tree->vArr[k].delta + 1);

                        }
                    }
                    STINGER_FORALL_EDGES_OF_VTX_END();

#if COUNT_TRAVERSALS==1
                    eAPT->dynamicTraverseEdgeCounter += stinger_typed_outdegree(sStinger, currElement, 0);
                    eAPT->dynamicTraverseVerticeCounter++;
#endif

                    if (currElement != currRoot) {
                        eAPT->sV[currElement].totalBC += eAPT->sV[currElement].newDelta - tree->vArr[currElement].delta;
                    }
                }

                #pragma omp master
                {
                    qDownBIndex -= 2;
                }
                #pragma omp barrier

            }

            if (operation == 1) {
                #pragma omp master
                {
                    *qStart = *qStart_nxt;
                    * qEnd = *qEnd_nxt;
                    * qStart_nxt = *qEnd;
                    * qEnd_nxt = *qStart_nxt;
                }
                #pragma omp barrier

                #pragma omp for

                for (int64_t i = *qStart; i < *qEnd; i++) {
                    int64_t currElement = QueueUp[i];
                    int64_t levelCurrMinusOne = tree->vArr[currElement].level - 1;
                    eAPT->sV[currElement].newEdgesBelow = 0;
                    STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement) {
                        int64_t k = STINGER_EDGE_DEST;
                        // Checking that the vertices are in different levels.

                        if (tree->vArr[k].level == tree->vArr[currElement].level + 1) {

                            if (eAPT->sV[k].touched == 0) {
                                eAPT->sV[currElement].newEdgesBelow += tree->vArr[k].edgesBelow + 1;
                            } else {
                                eAPT->sV[currElement].newEdgesBelow += eAPT->sV[k].newEdgesBelow + 1;
                            }
                        }

                        if (tree->vArr[k].level == levelCurrMinusOne) {
                            // Checking to see if "k" has been touched before.
                            if (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 0, -1)) {
                                eAPT->sV[k].newDelta += tree->vArr[k].delta;
                                upCounter++;
                                // Marking element as touched in the ascent stage.
                                eAPT->sV[k].touched = -1;

                                __sync_bool_compare_and_swap(&depthUp, -1, tree->vArr[k].level);
                                QueueUp[__atomic_fetch_add(qEnd_nxt, 1, __ATOMIC_RELAXED)] = k;

                                if (k != parentVertex)
                                    eAPT->sV[k].newSigma += tree->vArr[k].sigma;
                            }
                        }

                        if (tree->vArr[k].level == tree->vArr[currElement].level + 1 && eAPT->sV[k].touched != 0) {
                            eAPT->sV[currElement].newDelta +=
                                ((bc_t)eAPT->sV[currElement].newSigma / (bc_t)eAPT->sV[k].newSigma) *
                                (bc_t)(eAPT->sV[k].newDelta + 1);

                            // For the elements that are touched in the ascent stage it is necessary to
                            // to reduce the values that they previously had.
                            // In addition to this, the "parentVertex" that is connected to "startVertex", i.e.
                            // the vertices of the new edge, needs to increase its betweenness centrality
                            // following the new connection, without removing the old delta value.
                            if (eAPT->sV[currElement].touched < 0 && ( currElement != parentVertex || k != startVertex)) {
                                eAPT->sV[currElement].newDelta -=
                                    ((bc_t)tree->vArr[currElement].sigma / (bc_t)tree->vArr[k].sigma) *
                                    (bc_t)(tree->vArr[k].delta + 1);

                            }
                        }


                        // Vertices that did not move and that one of their neighbors move up(such that
                        // the vertices are now in the same level).
                        if (tree->vArr[k].level == tree->vArr[currElement].level && ((eAPT->sV[currElement].IMoved == 1 && eAPT->sV[k].IMoved < 0) )) {
                            // Checking to see if "k" has been touched before.
                            if (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 0, -1)) {
                                eAPT->sV[k].newDelta += tree->vArr[k].delta;

                                upCounter++;
                                // Marking element as touched in the ascent stage.
                                eAPT->sV[k].touched = -2;

                                __sync_bool_compare_and_swap(&depthSame, -1, tree->vArr[k].level);
                                QueueSame[__atomic_fetch_add(qEndSame_nxt, 1, __ATOMIC_RELAXED)] = k;
                                eAPT->sV[k].newSigma += tree->vArr[k].sigma;
                            }

                            // Paths that previosul went through this vertex no longer go through them, thus the
                            // shortest path count(BC) is reduced.

                        }

                        if (tree->vArr[k].level == tree->vArr[currElement].level && ((eAPT->sV[k].IMoved == 1 && eAPT->sV[currElement].IMoved < 0) )) {
                            eAPT->sV[currElement].newDelta -=
                                ((bc_t)tree->vArr[currElement].sigma / (bc_t)tree->vArr[k].sigma) *
                                (bc_t)(tree->vArr[k].delta + 1);

                        }
                    }
                    STINGER_FORALL_EDGES_OF_VTX_END();

#if COUNT_TRAVERSALS==1
                    eAPT->dynamicTraverseEdgeCounter += stinger_typed_outdegree(sStinger, currElement, 0);
                    eAPT->dynamicTraverseVerticeCounter++;
#endif

                    if (currElement != currRoot) {
                        eAPT->sV[currElement].totalBC += eAPT->sV[currElement].newDelta - tree->vArr[currElement].delta;
                    }
                }

                #pragma omp master
                {
                    *qStart = *qEnd;
                }
                #pragma omp barrier
            }

            if (operation == 2) {
                #pragma omp master
                {
                    *qStartSame = *qStartSame_nxt;
                    * qEndSame = *qEndSame_nxt;
                    * qStartSame_nxt = *qEndSame;
                    * qEndSame_nxt = *qStartSame_nxt;
                }
                #pragma omp barrier

                #pragma omp for

                for (int64_t i = *qStartSame; i < *qEndSame; i++) {
                    int64_t currElement = QueueSame[i];
                    int64_t levelCurrMinusOne = tree->vArr[currElement].level - 1;
                    eAPT->sV[currElement].newEdgesBelow = 0;
                    STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement) {
                        int64_t k = STINGER_EDGE_DEST;
                        // Checking that the vertices are in different levels.

                        if (tree->vArr[k].level == tree->vArr[currElement].level + 1) {

                            if (eAPT->sV[k].touched == 0) {
                                eAPT->sV[currElement].newEdgesBelow += tree->vArr[k].edgesBelow + 1;
                            } else {
                                eAPT->sV[currElement].newEdgesBelow += eAPT->sV[k].newEdgesBelow + 1;
                            }
                        }

                        if (tree->vArr[k].level == levelCurrMinusOne) {
                            // Checking to see if "k" has been touched before.
                            if (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 0, -1)) {
                                eAPT->sV[k].newDelta += tree->vArr[k].delta;
                                upCounter++;
                                // Marking element as touched in the ascent stage.
                                eAPT->sV[k].touched = -1;

                                __sync_bool_compare_and_swap(&depthUp, -1, tree->vArr[k].level);
                                QueueUp[__atomic_fetch_add(qEnd_nxt, 1, __ATOMIC_RELAXED)] = k;

                                if (k != parentVertex)
                                    eAPT->sV[k].newSigma += tree->vArr[k].sigma;
                            }
                        }

                        if (tree->vArr[k].level == tree->vArr[currElement].level + 1 && eAPT->sV[k].touched != 0) {
                            eAPT->sV[currElement].newDelta +=
                                ((bc_t)eAPT->sV[currElement].newSigma / (bc_t)eAPT->sV[k].newSigma) *
                                (bc_t)(eAPT->sV[k].newDelta + 1);

                            // For the elements that are touched in the ascent stage it is necessary to
                            // to reduce the values that they previously had.
                            // In addition to this, the "parentVertex" that is connected to "startVertex", i.e.
                            // the vertices of the new edge, needs to increase its betweenness centrality
                            // following the new connection, without removing the old delta value.
                            if (eAPT->sV[currElement].touched < 0 && ( currElement != parentVertex || k != startVertex)) {
                                eAPT->sV[currElement].newDelta -=
                                    ((bc_t)tree->vArr[currElement].sigma / (bc_t)tree->vArr[k].sigma) *
                                    (bc_t)(tree->vArr[k].delta + 1);

                            }
                        }


                        // Vertices that did not move and that one of their neighbors move up(such that
                        // the vertices are now in the same level).
                        if (tree->vArr[k].level == tree->vArr[currElement].level && ((eAPT->sV[currElement].IMoved == 1 && eAPT->sV[k].IMoved < 0) )) {
                            // Checking to see if "k" has been touched before.
                            if (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 0, -1)) {
                                eAPT->sV[k].newDelta += tree->vArr[k].delta;

                                upCounter++;
                                // Marking element as touched in the ascent stage.
                                eAPT->sV[k].touched = -2;
                                __sync_bool_compare_and_swap(&depthSame, -1, tree->vArr[k].level);
                                QueueSame[__atomic_fetch_add(qEndSame_nxt, 1, __ATOMIC_RELAXED)] = k;
                                eAPT->sV[k].newSigma += tree->vArr[k].sigma;
                            }

                            // Paths that previosul went through this vertex no longer go through them, thus the
                            // shortest path count(BC) is reduced.

                        }

                        if (tree->vArr[k].level == tree->vArr[currElement].level && ((eAPT->sV[k].IMoved == 1 && eAPT->sV[currElement].IMoved < 0) )) {
                            eAPT->sV[currElement].newDelta -=
                                ((bc_t)tree->vArr[currElement].sigma / (bc_t)tree->vArr[k].sigma) *
                                (bc_t)(tree->vArr[k].delta + 1);

                        }
                    }
                    STINGER_FORALL_EDGES_OF_VTX_END();

#if COUNT_TRAVERSALS==1
                    eAPT->dynamicTraverseEdgeCounter += stinger_typed_outdegree(sStinger, currElement, 0);
                    eAPT->dynamicTraverseVerticeCounter++;
#endif

                    if (currElement != currRoot) {
                        eAPT->sV[currElement].totalBC += eAPT->sV[currElement].newDelta - tree->vArr[currElement].delta;
                    }
                }

                #pragma omp master
                {
                    *qStartSame = *qEndSame;
                }
                #pragma omp barrier
            }
        }

        #pragma omp for

        for (int64_t c = 0; c <= qDownEndMarker; c++) {
            int64_t k = QueueDown[c];
            tree->vArr[k].delta = eAPT->sV[k].newDelta;
            tree->vArr[k].sigma = eAPT->sV[k].newSigma;
            tree->vArr[k].edgesBelow = eAPT->sV[k].newEdgesBelow;
            tree->vArr[k].edgesAbove = eAPT->sV[k].newEdgesAbove;
            eAPT->sV[k].diffPath = 0;
            eAPT->sV[k].touched = 0;
            eAPT->sV[k].newDelta = 0.0;
            eAPT->sV[k].movementDelta = 0;
            eAPT->sV[k].IMoved = -1;
            eAPT->sV[k].newSigma = 0;
            eAPT->sV[k].newEdgesAbove = 0;
            eAPT->sV[k].newEdgesBelow = 0;
        }

        eAPT->sV[startVertex].newEdgesAbove = 0;
        eAPT->sV[parentVertex].newEdgesAbove = 0;

        #pragma omp for

        for (int64_t c = 0; c < *qEndSame; c++) {
            int64_t k = QueueSame[c];
            tree->vArr[k].delta = eAPT->sV[k].newDelta;
            tree->vArr[k].edgesBelow = eAPT->sV[k].newEdgesBelow;
            eAPT->sV[k].diffPath = 0;
            eAPT->sV[k].touched = 0;
            eAPT->sV[k].newDelta = 0.0;
            eAPT->sV[k].movementDelta = 0;
            eAPT->sV[k].IMoved = -1;
            eAPT->sV[k].newSigma = 0;
            eAPT->sV[k].newEdgesBelow = 0;
        }

        #pragma omp for

        for (int64_t c = 0; c < *qEnd; c++) {
            int64_t k = QueueUp[c];
            tree->vArr[k].delta = eAPT->sV[k].newDelta;
            tree->vArr[k].edgesBelow = eAPT->sV[k].newEdgesBelow;
            eAPT->sV[k].diffPath = 0;
            eAPT->sV[k].touched = 0;
            eAPT->sV[k].newDelta = 0.0;
            eAPT->sV[k].movementDelta = 0;
            eAPT->sV[k].IMoved = -1;
            eAPT->sV[k].newSigma = 0;
            eAPT->sV[k].newEdgesBelow = 0;
        }

    } 
    eAPT->sV[startVertex].newEdgesBelow = 0;
    eAPT->sV[parentVertex].newEdgesBelow = 0;

    queue->size = 0;

    eAPT->qStart = 0;
    eAPT->qEnd = 0;
    eAPT->qStart_nxt = 0;
    eAPT->qEnd_nxt = 0;

    eAPT->qStartSame = 0;
    eAPT->qEndSame = 0;
    eAPT->qStartSame_nxt = 0;
    eAPT->qEndSame_nxt = 0;
}

// Case 2
void removeEdgeWithoutMovementBrandesFG(bcForest * forest, struct stinger * sStinger, int64_t currRoot,
                                        int64_t startVertex, int64_t parentVertex, int64_t deletedPathsFromRoot,
                                        extraArraysPerThread * eAPT, int64_t cores)
{
    bcTree * tree = forest->forest[currRoot];

    int64_t NV = forest->NV;
    int64_t * Queue = eAPT->QueueSame;
    int64_t * QueueDown = eAPT->QueueDown;
    int64_t * QueueUp = eAPT->QueueUp;
    int64_t * QueueDownBorders = eAPT->Stack;

    eAPT->sV[startVertex].newEdgesBelow = tree->vArr[startVertex].edgesBelow;
    eAPT->sV[parentVertex].newEdgesBelow = tree->vArr[parentVertex].edgesBelow;
    eAPT->sV[startVertex].newEdgesAbove = tree->vArr[startVertex].edgesAbove;
    eAPT->sV[parentVertex].newEdgesAbove = tree->vArr[parentVertex].edgesAbove;
    eAPT->sV[startVertex].newSigma = tree->vArr[startVertex].sigma;
    eAPT->sV[startVertex].touched = 1;
    eAPT->sV[startVertex].newSigma -= deletedPathsFromRoot;
    eAPT->sV[startVertex].diffPath = deletedPathsFromRoot;
    eAPT->sV[startVertex].newEdgesAbove -= eAPT->sV[parentVertex].newEdgesAbove + 1;
    eAPT->sV[parentVertex].newEdgesBelow -= eAPT->sV[startVertex].newEdgesBelow + 1;

    QueueDown[0] = startVertex;

    int64_t * qDownStart = &(eAPT->qStart);
    int64_t * qDownEnd = &(eAPT->qEnd);
    int64_t * qDownStart_nxt = &(eAPT->qStart_nxt);
    int64_t * qDownEnd_nxt = &(eAPT->qEnd_nxt);

    int64_t qDownBIndex = 0;
    *qDownEnd = 1;
    *qDownStart_nxt = 1;
    *qDownEnd_nxt = 1;

    int64_t deepestLevel = tree->vArr[startVertex].level;

    queue_t * queue = eAPT->queue;
    level_node_t * levelIndices = eAPT->levelIndices;

    // Starting BFS decent from "startVertex", down to all the vertices that have shortest paths through "startVertex"
    // All elements that will be touched will receive a positive value in their touched field.
    // In this implementation, "STACKS" are not used for the "moving up" stage. Rather, a multi-level queue is used.
    // Each level in the tree(max depth NV) has a queue and a counter specifiying how deep a specific deepth-queue is.
    // For simplicity, all elements are pushed both into the multi-level queue and into the regular queue which is used
    // for the BFS traversal.
    while (*qDownStart != *qDownEnd) {
        QueueDownBorders[qDownBIndex++] = *qDownStart;
        QueueDownBorders[qDownBIndex++] = *qDownEnd;
        int64_t thread_nums = cores;

        if ((*qDownEnd - *qDownStart) < cores) {
            thread_nums = *qDownEnd - *qDownStart;
        }

        #pragma omp parallel num_threads(thread_nums)
        {
            #pragma omp for

            for (int64_t i = *qDownStart; i < *qDownEnd; i++) {
                int64_t currElement = QueueDown[i];

                if (currElement != startVertex) {
                    __atomic_fetch_add(&(eAPT->sV[currElement].newEdgesAbove), tree->vArr[currElement].edgesAbove, __ATOMIC_RELAXED);
                }

                STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement) {
                    int64_t k = STINGER_EDGE_DEST;

                    if (currElement != startVertex
                            && tree->vArr[currElement].level - 1 == tree->vArr[k].level
                            && tree->vArr[currElement].level >= tree->vArr[startVertex].level) {

                        if (eAPT->sV[k].touched != 0) {
                            __atomic_fetch_add(&(eAPT->sV[currElement].newEdgesAbove), -tree->vArr[k].edgesAbove, __ATOMIC_RELAXED);
                            __atomic_fetch_add(&(eAPT->sV[currElement].newEdgesAbove), eAPT->sV[k].newEdgesAbove, __ATOMIC_RELAXED);
                        }
                    }

                    // if this vertex has not been added yet
                    if ((tree->vArr[currElement].level + 1) == (tree->vArr[k].level)) {
                        if (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 0, currElement)) {
                            // Checking if a "deeper level" has been reached.
                            if (deepestLevel < tree->vArr[k].level)
                                deepestLevel = tree->vArr[k].level;

                            // insert this vertex into the BFS queue
                            QueueDown[__atomic_fetch_add(qDownEnd_nxt, 1, __ATOMIC_RELAXED)] = k;

                            __atomic_fetch_add(&(eAPT->sV[k].newSigma), tree->vArr[k].sigma, __ATOMIC_RELAXED);

                            // indicate that it is in the next level of the BFS

                            // add new paths to root that go through current BFS Vertex
                            __atomic_fetch_add(&(eAPT->sV[k].newSigma), -eAPT->sV[currElement].diffPath, __ATOMIC_RELAXED);
                            // pass on my new paths to root for its search
                            __atomic_fetch_add(&(eAPT->sV[k].diffPath), eAPT->sV[currElement].diffPath, __ATOMIC_RELAXED);
                        }
                        // otherwise if it has been touched, but is specifically in the next level
                        // of the search (meaning it has more than one edge to the current level)
                        else if (eAPT->sV[k].touched != currElement) {
                            // add new paths to root that go through current BFS Vertex
                            __atomic_fetch_add(&(eAPT->sV[k].newSigma), -eAPT->sV[currElement].diffPath, __ATOMIC_RELAXED);

                            // pass on my new paths to root for its search
                            __atomic_fetch_add(&(eAPT->sV[k].diffPath), eAPT->sV[currElement].diffPath, __ATOMIC_RELAXED);
                        }
                    }
                }
                STINGER_FORALL_EDGES_OF_VTX_END();

#if COUNT_TRAVERSALS==1
                eAPT->dynamicTraverseEdgeCounter += stinger_typed_outdegree(sStinger, currElement, 0);
                eAPT->dynamicTraverseVerticeCounter++;
#endif
            }
        }
        *qDownStart = *qDownStart_nxt;
        *qDownEnd = *qDownEnd_nxt;
        *qDownStart_nxt = *qDownEnd;
        *qDownEnd_nxt = *qDownStart_nxt;
    }

    // The parent vertex needs to be placed in the queue for the dependency accumulation stage.
    // Also, it no longer has a child and so the delta from the child needs to be removed.

    int64_t qUpStart = 0, qUpEnd = 0;
    (*qDownEnd)--;
    int64_t qDownEndMarker = *qDownEnd;
    // Starting Multi-Level "BFS" ascent.

    // The ascent continues going up as long as the root has not been reached and that there
    // are elements in the current level of the ascent. The ascent starts in the deepest level
    // of the graph.
    // It was worth noting that in the ascent stage:
    // 1) All previously untouched elements that are touched are marked with "-1".
    // 2) On the way up, it is possible that elements that were touched in the BFS decent, will
    // touch elements that were not touchded in the decsent and that are below "vertex". These
    // are elements that do not have shortest paths going through "vertex" ,however, there BC
    // values have changed due to the changes occuring below them. Because of this, they are
    // placed in the Multi-level queue.

    *qDownStart = 0;
    *qDownEnd = 0;
    *qDownStart_nxt = 0;
    *qDownEnd_nxt = 0;

    while (!(qDownBIndex <= 0 && *qDownStart >= *qDownEnd && *qDownStart_nxt >= *qDownEnd_nxt)) {

        if (qDownBIndex >= 2) {

            int64_t thread_nums = cores;

            if ((QueueDownBorders[qDownBIndex - 1] - QueueDownBorders[qDownBIndex - 2]) < cores) {
                thread_nums = QueueDownBorders[qDownBIndex - 1] - QueueDownBorders[qDownBIndex - 2];
            }

            #pragma omp parallel num_threads(thread_nums)
            {
                #pragma omp for

                for (int64_t i = QueueDownBorders[qDownBIndex - 2]; i < QueueDownBorders[qDownBIndex - 1]; i++) {
                    int64_t currElement = QueueDown[i];

                    if (currElement != parentVertex && currElement != startVertex) {
                        __atomic_fetch_add(&(eAPT->sV[currElement].newEdgesBelow), tree->vArr[currElement].edgesBelow, __ATOMIC_RELAXED);
                    }

                    STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement) {
                        int64_t k = STINGER_EDGE_DEST;


                        if (currElement != parentVertex
                                && tree->vArr[currElement].level <= tree->vArr[parentVertex].level
                                && tree->vArr[k].level > tree->vArr[currElement].level) {


                            if (eAPT->sV[k].touched != 0) {
                                __atomic_fetch_add(&(eAPT->sV[k].newEdgesBelow), -tree->vArr[k].edgesBelow, __ATOMIC_RELAXED);
                                __atomic_fetch_add(&(eAPT->sV[k].newEdgesBelow), eAPT->sV[k].newEdgesBelow, __ATOMIC_RELAXED);
                            }
                        }


                        if (tree->vArr[k].level == tree->vArr[parentVertex].level && __sync_bool_compare_and_swap(&(eAPT->sV[parentVertex].touched), 0, -1)) {
                            QueueUp[__atomic_fetch_add(qDownEnd_nxt, 1, __ATOMIC_RELAXED)] = parentVertex;
                            __atomic_fetch_add(&(eAPT->sV[parentVertex].newSigma), tree->vArr[parentVertex].sigma, __ATOMIC_RELAXED);
                            eAPT->sV[parentVertex].newDelta = tree->vArr[parentVertex].delta -
                                                              ((bc_t)tree->vArr[parentVertex].sigma / (bc_t)tree->vArr[startVertex].sigma) *
                                                              (bc_t)(tree->vArr[startVertex].delta + 1);


                        }

                        // Checking that the vertices are in different levels.
                        if (tree->vArr[k].level == (tree->vArr[currElement].level - 1)) {
                            // Checking to see if "k" has been touched before.
                            if (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 0, -1)) {
                                eAPT->sV[k].newDelta = tree->vArr[k].delta;
                                // Marking element as touched in the ascent stage.
                                __atomic_fetch_add(&(eAPT->sV[k].newSigma), tree->vArr[k].sigma, __ATOMIC_RELAXED);
                                QueueUp[__atomic_fetch_add(qDownEnd_nxt, 1, __ATOMIC_RELAXED)] = k;
                            }
                        }

                        if (tree->vArr[k].level == tree->vArr[currElement].level + 1 && eAPT->sV[k].touched != 0) {

                            eAPT->sV[currElement].newDelta +=
                                ((bc_t)eAPT->sV[currElement].newSigma / (bc_t)eAPT->sV[k].newSigma) *
                                (bc_t)(eAPT->sV[k].newDelta + 1);

                            // For the elements that are touched in the ascent stage it is necessary to
                            // to reduce the values that they previously had.
                            // In addition to this, the "parentVertex" that is connected to "vertex", i.e.
                            // the vertices of the new edge, needs to increase its betweenness centrality
                            // following the new connection, without removing the old delta value.
                            if (eAPT->sV[currElement].touched < 0 && ( currElement != parentVertex || k != startVertex)) {
                                eAPT->sV[currElement].newDelta -=
                                    ((bc_t)tree->vArr[currElement].sigma / (bc_t)tree->vArr[k].sigma) *
                                    (bc_t)(tree->vArr[k].delta + 1);

                            }
                        }
                    }
                    STINGER_FORALL_EDGES_OF_VTX_END();

#if COUNT_TRAVERSALS==1
                    eAPT->dynamicTraverseEdgeCounter += stinger_typed_outdegree(sStinger, currElement, 0);
                    eAPT->dynamicTraverseVerticeCounter++;
#endif

                    if (currElement != currRoot) {
                        eAPT->sV[currElement].totalBC += eAPT->sV[currElement].newDelta - tree->vArr[currElement].delta;
                    }
                }
            }
        }

        qDownBIndex -= 2;

        *qDownStart = *qDownStart_nxt;
        *qDownEnd = *qDownEnd_nxt;
        *qDownStart_nxt = *qDownEnd;
        *qDownEnd_nxt = *qDownStart_nxt;

        int64_t thread_nums = cores;

        if (*qDownEnd - *qDownStart < cores) {
            thread_nums = *qDownEnd - *qDownStart;
        }

        #pragma omp parallel num_threads(thread_nums)
        {
            #pragma omp for

            for (int64_t i = *qDownStart; i < *qDownEnd; i++) {
                int64_t currElement = QueueUp[i];

                if (currElement != parentVertex && currElement != startVertex) {
                    __atomic_fetch_add(&(eAPT->sV[currElement].newEdgesBelow), tree->vArr[currElement].edgesBelow, __ATOMIC_RELAXED);
                }

                STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement) {
                    int64_t k = STINGER_EDGE_DEST;


                    if (currElement != parentVertex
                            && tree->vArr[currElement].level <= tree->vArr[parentVertex].level
                            && tree->vArr[k].level > tree->vArr[currElement].level) {

                        if (eAPT->sV[k].touched != 0) {
                            __atomic_fetch_add(&(eAPT->sV[currElement].newEdgesBelow), -tree->vArr[k].edgesBelow, __ATOMIC_RELAXED);
                            __atomic_fetch_add(&(eAPT->sV[currElement].newEdgesBelow), eAPT->sV[k].newEdgesBelow, __ATOMIC_RELAXED);
                        }
                    }


                    if (tree->vArr[k].level == tree->vArr[parentVertex].level && __sync_bool_compare_and_swap(&(eAPT->sV[parentVertex].touched), 0, -1)) {
                        QueueUp[__atomic_fetch_add(qDownEnd_nxt, 1, __ATOMIC_RELAXED)] = parentVertex;
                        __atomic_fetch_add(&(eAPT->sV[parentVertex].newSigma), tree->vArr[parentVertex].sigma, __ATOMIC_RELAXED);
                        eAPT->sV[parentVertex].newDelta = tree->vArr[parentVertex].delta -
                                                          ((bc_t)tree->vArr[parentVertex].sigma / (bc_t)tree->vArr[startVertex].sigma) *
                                                          (bc_t)(tree->vArr[startVertex].delta + 1);


                    }

                    // Checking that the vertices are in different levels.
                    if (tree->vArr[k].level == (tree->vArr[currElement].level - 1)) {
                        // Checking to see if "k" has been touched before.
                        if (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 0, -1)) {
                            eAPT->sV[k].newDelta = tree->vArr[k].delta;
                            // Marking element as touched in the ascent stage.
                            __atomic_fetch_add(&(eAPT->sV[k].newSigma), tree->vArr[k].sigma, __ATOMIC_RELAXED);
                            QueueUp[__atomic_fetch_add(qDownEnd_nxt, 1, __ATOMIC_RELAXED)] = k;
                        }
                    }

                    if (tree->vArr[k].level == tree->vArr[currElement].level + 1 && eAPT->sV[k].touched != 0) {

                        eAPT->sV[currElement].newDelta +=
                            ((bc_t)eAPT->sV[currElement].newSigma / (bc_t)eAPT->sV[k].newSigma) *
                            (bc_t)(eAPT->sV[k].newDelta + 1);

                        // For the elements that are touched in the ascent stage it is necessary to
                        // to reduce the values that they previously had.
                        // In addition to this, the "parentVertex" that is connected to "vertex", i.e.
                        // the vertices of the new edge, needs to increase its betweenness centrality
                        // following the new connection, without removing the old delta value.
                        if (eAPT->sV[currElement].touched < 0 && ( currElement != parentVertex || k != startVertex)) {
                            eAPT->sV[currElement].newDelta -=
                                ((bc_t)tree->vArr[currElement].sigma / (bc_t)tree->vArr[k].sigma) *
                                (bc_t)(tree->vArr[k].delta + 1);

                        }
                    }
                }
                STINGER_FORALL_EDGES_OF_VTX_END();

#if COUNT_TRAVERSALS==1
                eAPT->dynamicTraverseEdgeCounter += stinger_typed_outdegree(sStinger, currElement, 0);
                eAPT->dynamicTraverseVerticeCounter++;
#endif

                if (currElement != currRoot) {
                    eAPT->sV[currElement].totalBC += eAPT->sV[currElement].newDelta - tree->vArr[currElement].delta;
                }
            }
        }
        *qDownStart = *qDownEnd;
    }


    for (int64_t q = 0; q <= qDownEndMarker; q++) {
        int64_t k = QueueDown[q];

        if (eAPT->sV[k].touched != 0) {
            tree->vArr[k].delta = eAPT->sV[k].newDelta;
            tree->vArr[k].sigma = eAPT->sV[k].newSigma;
        }

        tree->vArr[k].edgesAbove = eAPT->sV[k].newEdgesAbove;
        eAPT->sV[k].newEdgesAbove = 0;
        eAPT->sV[k].newEdgesBelow = 0;
        eAPT->sV[k].diffPath = 0;
        eAPT->sV[k].touched = 0;
        eAPT->sV[k].newDelta = 0.0;
        eAPT->sV[k].newSigma = 0;
    }

    eAPT->sV[startVertex].newEdgesAbove = 0;
    eAPT->sV[parentVertex].newEdgesAbove = 0;

    for (int64_t q = 0; q < *qDownEnd; q++) {
        int64_t k = QueueUp[q];

        if (eAPT->sV[k].touched != 0) {
            tree->vArr[k].delta = eAPT->sV[k].newDelta;
            tree->vArr[k].sigma = eAPT->sV[k].newSigma;
        }

        tree->vArr[k].edgesBelow = eAPT->sV[k].newEdgesBelow;
        eAPT->sV[k].newEdgesBelow = 0;
        eAPT->sV[k].newEdgesAbove = 0;
        eAPT->sV[k].diffPath = 0;
        eAPT->sV[k].touched = 0;
        eAPT->sV[k].newDelta = 0.0;
        eAPT->sV[k].newSigma = 0;
    }

    eAPT->sV[startVertex].newEdgesBelow = 0;
    eAPT->sV[parentVertex].newEdgesBelow = 0;

    queue->size = 0;

    eAPT->qStart = 0;
    eAPT->qEnd = 0;
    eAPT->qStart_nxt = 0;
    eAPT->qEnd_nxt = 0;

}

void moveDownTreeBrandesFG(bcForest * forest, struct stinger * sStinger, int64_t currRoot,
                           int64_t startVertex, int64_t parentVertex, extraArraysPerThread * eAPT, int64_t cores)
{
    bcTree * tree = forest->forest[currRoot];

    int64_t NV = forest->NV;
    int64_t * Queue = eAPT->QueueDown;
    int64_t * QueueUp = eAPT->QueueUp;
    int64_t * topQueue = eAPT->QueueSame;
    int64_t * QueueDownBorders = eAPT->Stack;
    int64_t * tqBorders = eAPT->tqBorders;

    int64_t * touchedVerticesDown = eAPT->touchedVerticesDown;
    int64_t * touchedVerticesUp = eAPT->touchedVerticesUp;
    queue_t * queue = eAPT->queue;
    level_node_t * levelIndices = eAPT->levelIndices;

    Queue[0] = startVertex;
    int64_t tvDownEnd = 0, tvUpEnd = 0;
    int64_t stopLevel = tree->vArr[startVertex].level;

    int64_t * qStart = &(eAPT->qStart);
    int64_t * qEnd = &(eAPT->qEnd);
    int64_t * qStart_nxt = &(eAPT->qStart_nxt);
    int64_t * qEnd_nxt = &(eAPT->qEnd_nxt);

    int64_t * tqStart = &(eAPT->tqStart);
    int64_t * tqEnd = &(eAPT->tqEnd);
    int64_t * tqStart_nxt = &(eAPT->tqStart_nxt);
    int64_t * tqEnd_nxt = &(eAPT->tqEnd_nxt);

    int64_t qDownBIndex = 0, tqBIndex = 0;
    *qEnd = 1;
    *qStart_nxt = 1;
    *qEnd_nxt = 1;

    *tqStart = 0;
    *tqEnd = 0;
    *tqStart_nxt = 0;
    *tqEnd_nxt = 0;

    eAPT->sV[startVertex].newLevel = INFINITY_MY;
    eAPT->sV[startVertex].newSigma = INFINITY_MY;
    eAPT->sV[startVertex].newDelta = 0.0;
    eAPT->sV[startVertex].newEdgesAbove = INFINITY_MY;
    eAPT->sV[startVertex].touched = 1;
    touchedVerticesDown[tvDownEnd++] = startVertex;

    int64_t deepestLevel = stopLevel;

    *qStart = 0;

    while (*qStart != *qEnd) {
        int64_t thread_nums = cores;

        if (*qEnd - *qStart < cores) {
            thread_nums = *qEnd - *qStart;
        }

        #pragma omp parallel num_threads(thread_nums)
        {
            #pragma omp for

            for (int64_t i = *qStart; i < *qEnd; i++) {
                int64_t currElement = Queue[i];

                STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement) {
                    int64_t k = STINGER_EDGE_DEST;

                    if (tree->vArr[k].level == tree->vArr[currElement].level + 1 && __sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 0, -1)) {
                        __atomic_fetch_add(&(eAPT->sV[k].newEdgesAbove), INFINITY_MY, __ATOMIC_RELAXED);
                        __atomic_fetch_add(&(eAPT->sV[k].newLevel), INFINITY_MY, __ATOMIC_RELAXED);
                        __atomic_fetch_add(&(eAPT->sV[k].newSigma), INFINITY_MY, __ATOMIC_RELAXED);
                        touchedVerticesDown[__atomic_fetch_add(&tvDownEnd, 1, __ATOMIC_RELAXED)] = k;
                        Queue[__atomic_fetch_add(qEnd_nxt, 1, __ATOMIC_RELAXED)] = k;
                        eAPT->sV[k].newDelta = 0.0;
                    }
                }
                STINGER_FORALL_EDGES_OF_VTX_END();

                int64_t parentOutsideSubtree = 0;
                int64_t siblingOutsideSubtree = 0;
                int64_t parentPathsToRoot = 0;
                int64_t siblingPathsToRoot = 0;
                int64_t parentEdgesAbove = 0;
                int64_t siblingEdgesAbove = 0;

                STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement) {
                    int64_t l = STINGER_EDGE_DEST;

                    if (tree->vArr[l].level == tree->vArr[currElement].level - 1) {
                        if (eAPT->sV[l].touched == 0) {
                            parentOutsideSubtree = l;
                            parentPathsToRoot += tree->vArr[l].sigma;
                        } else {
                            parentPathsToRoot += tree->vArr[l].sigma - 1;
                        }

                        parentEdgesAbove += tree->vArr[l].edgesAbove + 1;
                    } else if (tree->vArr[l].level == tree->vArr[currElement].level) {
                        if (eAPT->sV[l].touched == 0) {
                            siblingOutsideSubtree = l;
                            siblingPathsToRoot += tree->vArr[l].sigma;
                        } else {
                            siblingPathsToRoot += tree->vArr[l].sigma - 1;
                        }

                        siblingEdgesAbove += tree->vArr[l].edgesAbove + 1;
                    }
                }
                STINGER_FORALL_EDGES_OF_VTX_END();

                if (parentOutsideSubtree) {

                    if (eAPT->sV[currElement].touched == 1 || eAPT->sV[currElement].touched == SIBLING_ANCHORED) {
                        topQueue[__atomic_fetch_add(tqEnd_nxt, 1, __ATOMIC_RELAXED)] = currElement;
                        eAPT->sV[currElement].newLevel = tree->vArr[parentOutsideSubtree].level + 1;
                    }

                    eAPT->sV[currElement].touched = PARENT_ANCHORED;
                    eAPT->sV[currElement].newDelta = 0.0;
                } else if (siblingOutsideSubtree) {

                    if (eAPT->sV[currElement].touched == 1) {
                        topQueue[__atomic_fetch_add(tqEnd_nxt, 1, __ATOMIC_RELAXED)] = currElement;
                        eAPT->sV[currElement].newLevel = tree->vArr[siblingOutsideSubtree].level + 1;
                    }

                    if (eAPT->sV[currElement].touched != PARENT_ANCHORED) {
                        eAPT->sV[currElement].touched = SIBLING_ANCHORED;
                    }

                    eAPT->sV[currElement].newDelta = 0.0;
                }

                STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement) {
                    int64_t k = STINGER_EDGE_DEST;

                    if (tree->vArr[k].level == tree->vArr[currElement].level + 1 && (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), -1, -2)
                            || __sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 1, -2))) {

                        if (eAPT->sV[currElement].touched == PARENT_ANCHORED && (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), -2, PARENT_ANCHORED)
                                || __sync_bool_compare_and_swap(&(eAPT->sV[k].touched), SIBLING_ANCHORED, PARENT_ANCHORED))) {}

                        else if (eAPT->sV[currElement].touched == SIBLING_ANCHORED && __sync_bool_compare_and_swap(&(eAPT->sV[k].touched), -2, SIBLING_ANCHORED)) {}

                        else {
                            __atomic_fetch_add(&(eAPT->sV[k].touched), 3, __ATOMIC_RELAXED);
                        }
                    }
                }
                STINGER_FORALL_EDGES_OF_VTX_END();

            }
        }
        *qStart = *qStart_nxt;
        *qEnd = *qEnd_nxt;
        *qStart_nxt = *qEnd;
        *qEnd_nxt = *qStart_nxt;

        *tqStart = *tqStart_nxt;
        *tqEnd = *tqEnd_nxt;
        *tqStart_nxt = *tqEnd;
        *tqEnd_nxt = *tqStart_nxt;
    }

    *qEnd = 1;
    *qStart_nxt = 1;
    *qEnd_nxt = 1;
    *qStart = 0;

    *tqStart = 0;

    int64_t key, j;

    for (int64_t i = 1; i < *tqEnd; i++) {
        key = topQueue[i];
        j = i - 1;

        while (j >= 0 && eAPT->sV[topQueue[j]].newLevel > eAPT->sV[key].newLevel) {
            topQueue[j + 1] = topQueue[j];
            j = j - 1;
        }

        topQueue[j + 1] = key;
    }

    int64_t lo = 0;
    int64_t hi = 0;

    while (lo < *tqEnd && hi < *tqEnd) {
        while (lo < *tqEnd && hi < *tqEnd && eAPT->sV[topQueue[lo]].newLevel == eAPT->sV[topQueue[hi]].newLevel) {
            hi++;
        }

        tqBorders[tqBIndex++] = lo;
        tqBorders[tqBIndex++] = hi;
        lo = hi;
    }

    // While queue is not empty

    int64_t tqLevel = 0;

    if (*tqEnd != 0) {
        appendDS(queue, levelIndices, eAPT->sV[topQueue[*tqStart]].newLevel, topQueue[*tqStart]);
        Queue[0] = topQueue[(*tqStart)++];
        tqBorders[0]++;

        if (tqBorders[0] == tqBorders[1])
            tqLevel += 2;

        eAPT->sV[Queue[0]].touched = 5;
    } else {
        Queue[0] = startVertex;
        eAPT->sV[Queue[0]].touched = 5;
        *qStart = 0;
        *qEnd = 1;
    }

    while (*qStart != *qEnd) {

        if (tqLevel < tqBIndex && tqBorders[tqLevel] < *tqEnd && eAPT->sV[topQueue[tqBorders[tqLevel]]].newLevel <= eAPT->sV[Queue[*qStart]].newLevel) {
            int64_t thread_nums = cores;

            if (tqBorders[tqLevel + 1] - tqBorders[tqLevel] < cores) {
                thread_nums = tqBorders[tqLevel + 1] - tqBorders[tqLevel];
            }

            #pragma omp parallel num_threads(thread_nums)
            {
                #pragma omp for

                for (int64_t i = tqBorders[tqLevel]; i < tqBorders[tqLevel + 1]; i++) {

                    if (__sync_bool_compare_and_swap(&(eAPT->sV[topQueue[i]].touched), 1, 5) ||
                            __sync_bool_compare_and_swap(&(eAPT->sV[topQueue[i]].touched), PARENT_ANCHORED, 5) ||
                            __sync_bool_compare_and_swap(&(eAPT->sV[topQueue[i]].touched), SIBLING_ANCHORED, 5)) {

                        Queue[__atomic_fetch_add(qEnd, 1, __ATOMIC_RELAXED)] = topQueue[i];
                        appendDS(queue, levelIndices, eAPT->sV[topQueue[i]].newLevel, topQueue[i]);
                    }
                }
            }
            tqLevel += 2;
        }

        *qStart_nxt = *qEnd;
        *qEnd_nxt = *qStart_nxt;

        QueueDownBorders[qDownBIndex++] = *qStart;
        QueueDownBorders[qDownBIndex++] = *qEnd;

        int64_t thread_nums = cores;

        if (*qEnd - *qStart < cores) {
            thread_nums = *qEnd - *qStart;
        }

        #pragma omp parallel num_threads(thread_nums)
        {
            #pragma omp for

            for (int64_t i = *qStart; i < *qEnd; i++) {
                int64_t currElement = Queue[i];

                eAPT->sV[currElement].newEdgesAbove = 0;

                if (deepestLevel < eAPT->sV[currElement].newLevel) {
                    deepestLevel = eAPT->sV[currElement].newLevel;
                }

                STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement) {
                    int64_t k = STINGER_EDGE_DEST;

                    if (eAPT->sV[k].newLevel > eAPT->sV[currElement].newLevel) {
                        // Checking if "k" has been found.

                        __sync_bool_compare_and_swap(&(eAPT->sV[k].newLevel), INFINITY_MY, eAPT->sV[currElement].newLevel + 1);

                        if (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 1, 5) ||
                                __sync_bool_compare_and_swap(&(eAPT->sV[k].touched), PARENT_ANCHORED, 5) ||
                                __sync_bool_compare_and_swap(&(eAPT->sV[k].touched), SIBLING_ANCHORED, 5)) {

                            Queue[__atomic_fetch_add(qEnd_nxt, 1, __ATOMIC_RELAXED)] = k;
                            eAPT->sV[k].newDelta = 0.0;

                            if (deepestLevel < eAPT->sV[k].newLevel)
                                deepestLevel = eAPT->sV[k].newLevel;

                            appendDS(queue, levelIndices, eAPT->sV[k].newLevel, k);
                        }
                    }

                    if (eAPT->sV[currElement].newLevel == tree->vArr[k].level + 1 && eAPT->sV[k].touched == 0) {
                        if (__sync_bool_compare_and_swap(&(eAPT->sV[currElement].newSigma), INFINITY_MY, tree->vArr[k].sigma)) {}
                        else {
                            __atomic_fetch_add(&(eAPT->sV[currElement].newSigma), tree->vArr[k].sigma, __ATOMIC_RELAXED);
                        }

                        if (__sync_bool_compare_and_swap(&(eAPT->sV[currElement].newEdgesAbove), INFINITY_MY, tree->vArr[k].edgesAbove + 1)) {}
                        else {
                            __atomic_fetch_add(&(eAPT->sV[currElement].newEdgesAbove), tree->vArr[k].edgesAbove + 1, __ATOMIC_RELAXED);
                        }
                    } else if (eAPT->sV[currElement].newLevel == eAPT->sV[k].newLevel + 1 && eAPT->sV[k].touched != 0) {

                        if (__sync_bool_compare_and_swap(&(eAPT->sV[currElement].newSigma), INFINITY_MY, eAPT->sV[k].newSigma)) {
                        } else {
                            __atomic_fetch_add(&(eAPT->sV[currElement].newSigma), eAPT->sV[k].newSigma, __ATOMIC_RELAXED);
                        }

                        if (__sync_bool_compare_and_swap(&(eAPT->sV[currElement].newEdgesAbove), INFINITY_MY, eAPT->sV[k].newEdgesAbove)) {
                        } else {
                            __atomic_fetch_add(&(eAPT->sV[currElement].newEdgesAbove), eAPT->sV[k].newEdgesAbove + 1, __ATOMIC_RELAXED);
                        }
                    }
                }


                STINGER_FORALL_EDGES_OF_VTX_END();
#if COUNT_TRAVERSALS==1
                eAPT->dynamicTraverseEdgeCounter += stinger_typed_outdegree(sStinger, currElement, 0);
                eAPT->dynamicTraverseVerticeCounter++;
#endif
            }
        }
        *qStart = *qStart_nxt;
        *qEnd = *qEnd_nxt;
        *qStart_nxt = *qEnd;
        *qEnd_nxt = *qStart_nxt;
    }

    *qEnd = 0;

    // If it is not a case 4
    if (deepestLevel != INFINITY_MY) {
        for (int64_t lev = stopLevel; lev <= deepestLevel; lev++) {
            int64_t index = levelIndices[lev].front;
            int64_t levelEmpty = 1;

            while (index != -1) {
                levelEmpty = 0;
                queue_node_t * temp_node = queue->nodes + index;
                Queue[(*qEnd)++] = temp_node->data;
                index = temp_node->next;
            }

            levelIndices[lev].front = -1;
            levelIndices[lev].back = -1;
        }

        queue->size = 0;
    }

    int64_t * qUpStart = &(eAPT->qStartSame);
    int64_t * qUpEnd = &(eAPT->qEndSame);
    int64_t * qUpStart_nxt = &(eAPT->qStartSame_nxt);
    int64_t * qUpEnd_nxt = &(eAPT->qEndSame_nxt);
    *qUpStart = 0;
    *qUpEnd = 0;
    *qUpStart_nxt = 0;
    *qUpEnd_nxt = 0;

    (*qEnd)--;

    int case4 = 0;

    if (eAPT->sV[startVertex].newLevel == INFINITY_MY) {
        if (eAPT->sV[parentVertex].touched == 0) {
            eAPT->sV[parentVertex].newDelta = tree->vArr[parentVertex].delta -
                                              ((bc_t)tree->vArr[parentVertex].sigma / (bc_t)tree->vArr[startVertex].sigma) *
                                              (bc_t)(tree->vArr[startVertex].delta + 1);

            eAPT->sV[parentVertex].newSigma = tree->vArr[parentVertex].sigma;
            eAPT->sV[parentVertex].touched = -1;
            QueueUp[(*qUpEnd_nxt)++] = parentVertex;
            case4 = 1;
        }
    }

    // Starting Multi-Level "BFS" ascent.

    // The ascent continues going up as long as the root has not been reached and that there
    // are elements in the current level of the ascent. The ascent starts in the deepest level
    // of the graph.
    // It was worth noting that in the ascent stage:
    // 1) All previously untouched elements that are touched are marked with "-1".
    // 2) On the way up, it is possible that elements that were touched in the BFS decent, will
    // touch elements that were not touchded in the decsent and that are below "vertex". These
    // are elements that do not have shortest paths going through "vertex" ,however, there BC
    // values have changed due to the changes occuring below them. Because of this, they are
    // placed in the Multi-level queue.

    while (!(qDownBIndex <= 0 && *qUpStart >= *qUpEnd && *qUpStart_nxt >= *qUpEnd_nxt)) {

        if (qDownBIndex >= 2 && !case4) {

            int64_t thread_nums = cores;

            if (QueueDownBorders[qDownBIndex - 1] - QueueDownBorders[qDownBIndex - 2] < cores) {
                thread_nums = QueueDownBorders[qDownBIndex - 1] - QueueDownBorders[qDownBIndex - 2];
            }

            #pragma omp parallel num_threads(thread_nums)
            {
                #pragma omp for

                for (int64_t i = QueueDownBorders[qDownBIndex - 2]; i < QueueDownBorders[qDownBIndex - 1]; i++) {
                    int64_t currElement = Queue[i];
                    eAPT->sV[currElement].newEdgesBelow = 0;
                    touchedVerticesUp[__atomic_fetch_add(&tvUpEnd, 1, __ATOMIC_RELAXED)] = currElement;

                    int64_t currElementLevel = eAPT->sV[currElement].newLevel;
                    __sync_bool_compare_and_swap(&currElementLevel, 0, tree->vArr[currElement].level);

                    int64_t parentVertexLevel = eAPT->sV[parentVertex].newLevel;
                    __sync_bool_compare_and_swap(&parentVertexLevel, 0, tree->vArr[parentVertex].level);

                    STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement) {
                        int64_t k = STINGER_EDGE_DEST;

                        int64_t kLevel = eAPT->sV[k].newLevel;
                        __sync_bool_compare_and_swap(&kLevel, 0, tree->vArr[k].level);

                        if (kLevel == currElementLevel + 1) {
                            if (eAPT->sV[k].touched != 0) {
                                __atomic_fetch_add(&(eAPT->sV[currElement].newEdgesBelow), eAPT->sV[k].newEdgesBelow + 1, __ATOMIC_RELAXED);
                            } else if (eAPT->sV[k].touched == 0) {
                                __atomic_fetch_add(&(eAPT->sV[currElement].newEdgesBelow), tree->vArr[k].edgesBelow + 1, __ATOMIC_RELAXED);
                            }
                        }

                        if (kLevel == parentVertexLevel) {
                            if (__sync_bool_compare_and_swap(&(eAPT->sV[parentVertex].touched), 0, -1)) {
                                eAPT->sV[parentVertex].newDelta = tree->vArr[parentVertex].delta -
                                                                  ((bc_t)tree->vArr[parentVertex].sigma / (bc_t)tree->vArr[startVertex].sigma) *
                                                                  (bc_t)(tree->vArr[startVertex].delta + 1);

                                __atomic_fetch_add(&(eAPT->sV[parentVertex].newSigma), tree->vArr[parentVertex].sigma, __ATOMIC_RELAXED);
                                QueueUp[__atomic_fetch_add(qUpEnd_nxt, 1, __ATOMIC_RELAXED)] = parentVertex;
                            }
                        }

                        // Checking that the vertices are in different levels.
                        if (kLevel == currElementLevel - 1) {

                            // Checking to see if "k" has been touched before.
                            if (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 0, -1)) {
                                eAPT->sV[k].newDelta = tree->vArr[k].delta;

                                // Marking element as touched in the ascent stage.
                                __atomic_fetch_add(&(eAPT->sV[k].newSigma), tree->vArr[k].sigma, __ATOMIC_RELAXED);
                                QueueUp[__atomic_fetch_add(qUpEnd_nxt, 1, __ATOMIC_RELAXED)] = k;
                            }

                        }

                        if (kLevel == currElementLevel + 1 && eAPT->sV[k].touched != 0) {

                            eAPT->sV[currElement].newDelta +=
                                ((bc_t)eAPT->sV[currElement].newSigma / (bc_t)eAPT->sV[k].newSigma) *
                                (bc_t)(eAPT->sV[k].newDelta + 1);

                            // For the elements that are touched in the ascent stage it is necessary to
                            // to reduce the values that they previously had.
                            // In addition to this, the "parentVertex" that is connected to "vertex", i.e.
                            // the vertices of the new edge, needs to increase its betweenness centrality
                            // following the new connection, without removing the old delta value.

                            if (eAPT->sV[currElement].touched < 0 && tree->vArr[currElement].level < tree->vArr[k].level) { 
                                eAPT->sV[currElement].newDelta -=
                                    ((bc_t)tree->vArr[currElement].sigma / (bc_t)tree->vArr[k].sigma) *
                                    (bc_t)(tree->vArr[k].delta + 1);
                            }
                        }
                    }
                    STINGER_FORALL_EDGES_OF_VTX_END();
#if COUNT_TRAVERSALS==1
                    eAPT->dynamicTraverseEdgeCounter += stinger_typed_outdegree(sStinger, currElement, 0);
                    eAPT->dynamicTraverseVerticeCounter++;
#endif

                    if (currElement != currRoot) {
                        eAPT->sV[currElement].totalBC += eAPT->sV[currElement].newDelta - tree->vArr[currElement].delta;
                    }
                }
            }
        }

        qDownBIndex -= 2;
        *qUpStart = *qUpStart_nxt;
        *qUpEnd = *qUpEnd_nxt;
        *qUpStart_nxt = *qUpEnd;
        *qUpEnd_nxt = *qUpStart_nxt;

        int64_t thread_nums = cores;

        if (*qUpEnd - *qUpStart < cores) {
            thread_nums = *qUpEnd - *qUpStart;
        }

        #pragma omp parallel num_threads(thread_nums)
        {
            #pragma omp for

            for (int64_t i = *qUpStart; i < *qUpEnd; i++) {
                int64_t currElement = QueueUp[i];
                eAPT->sV[currElement].newEdgesBelow = 0;
                touchedVerticesUp[__atomic_fetch_add(&tvUpEnd, 1, __ATOMIC_RELAXED)] = currElement;

                int64_t currElementLevel = eAPT->sV[currElement].newLevel;
                __sync_bool_compare_and_swap(&currElementLevel, 0, tree->vArr[currElement].level);

                int64_t parentVertexLevel = eAPT->sV[parentVertex].newLevel;
                __sync_bool_compare_and_swap(&parentVertexLevel, 0, tree->vArr[parentVertex].level);

                STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement) {
                    int64_t k = STINGER_EDGE_DEST;

                    int64_t kLevel = eAPT->sV[k].newLevel;
                    __sync_bool_compare_and_swap(&kLevel, 0, tree->vArr[k].level);

                    if (kLevel == currElementLevel + 1) {
                        if (eAPT->sV[k].touched != 0) {
                            __atomic_fetch_add(&(eAPT->sV[currElement].newEdgesBelow), eAPT->sV[k].newEdgesBelow + 1, __ATOMIC_RELAXED);
                        } else if (eAPT->sV[k].touched == 0) {
                            __atomic_fetch_add(&(eAPT->sV[currElement].newEdgesBelow), tree->vArr[k].edgesBelow + 1, __ATOMIC_RELAXED);
                        }
                    }

                    if (kLevel == parentVertexLevel) {
                        if (__sync_bool_compare_and_swap(&(eAPT->sV[parentVertex].touched), 0, -1)) {
                            eAPT->sV[parentVertex].newDelta = tree->vArr[parentVertex].delta -
                                                              ((bc_t)tree->vArr[parentVertex].sigma / (bc_t)tree->vArr[startVertex].sigma) *
                                                              (bc_t)(tree->vArr[startVertex].delta + 1);

                            __atomic_fetch_add(&(eAPT->sV[parentVertex].newSigma), tree->vArr[parentVertex].sigma, __ATOMIC_RELAXED);
                            QueueUp[__atomic_fetch_add(qUpEnd_nxt, 1, __ATOMIC_RELAXED)] = parentVertex;
                        }
                    }

                    // Checking that the vertices are in different levels.
                    if (kLevel == currElementLevel - 1) {

                        // Checking to see if "k" has been touched before.
                        if (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 0, -1)) {
                            eAPT->sV[k].newDelta = tree->vArr[k].delta;

                            // Marking element as touched in the ascent stage.
                            __atomic_fetch_add(&(eAPT->sV[k].newSigma), tree->vArr[k].sigma, __ATOMIC_RELAXED);
                            QueueUp[__atomic_fetch_add(qUpEnd_nxt, 1, __ATOMIC_RELAXED)] = k;
                        }

                    }

                    if (kLevel == currElementLevel + 1 && eAPT->sV[k].touched != 0) {

                        eAPT->sV[currElement].newDelta +=
                            ((bc_t)eAPT->sV[currElement].newSigma / (bc_t)eAPT->sV[k].newSigma) *
                            (bc_t)(eAPT->sV[k].newDelta + 1);

                        // For the elements that are touched in the ascent stage it is necessary to
                        // to reduce the values that they previously had.
                        // In addition to this, the "parentVertex" that is connected to "vertex", i.e.
                        // the vertices of the new edge, needs to increase its betweenness centrality
                        // following the new connection, without removing the old delta value.

                        if (eAPT->sV[currElement].touched < 0 && tree->vArr[currElement].level < tree->vArr[k].level) { 
                            eAPT->sV[currElement].newDelta -=
                                ((bc_t)tree->vArr[currElement].sigma / (bc_t)tree->vArr[k].sigma) *
                                (bc_t)(tree->vArr[k].delta + 1);
                        }
                    }
                }
                STINGER_FORALL_EDGES_OF_VTX_END();
#if COUNT_TRAVERSALS==1
                eAPT->dynamicTraverseEdgeCounter += stinger_typed_outdegree(sStinger, currElement, 0);
                eAPT->dynamicTraverseVerticeCounter++;
#endif

                if (currElement != currRoot) {
                    eAPT->sV[currElement].totalBC += eAPT->sV[currElement].newDelta - tree->vArr[currElement].delta;
                }
            }
        }
        *qUpStart = *qUpEnd;
    }

    for (int64_t k = 0; k < tvDownEnd; k++) {
        int64_t vertex = touchedVerticesDown[k];
        tree->vArr[vertex].level = eAPT->sV[vertex].newLevel;
    }

    // Handles case where edge deletion creates new connected component.
    if (tree->vArr[startVertex].level == INFINITY_MY) {
        *qStart = 0;
        *qEnd = 1;
        *qStart_nxt = 1;
        *qEnd_nxt = 1;
        Queue[0] = startVertex;

        eAPT->sV[startVertex].touched = -2;

        while (*qStart != *qEnd) {
            int64_t thread_nums = cores;

            if (*qEnd - *qStart < cores) {
                thread_nums = *qEnd - *qStart;
            }

            #pragma omp parallel num_threads(thread_nums)
            {
                #pragma omp for

                for (int64_t i = *qStart; i < *qEnd; i++) {
                    int64_t currElement = Queue[i];

                    eAPT->sV[currElement].totalBC -= tree->vArr[currElement].delta;

                    tree->vArr[currElement].edgesBelow = 0;
                    tree->vArr[currElement].edgesAbove = 0;
                    eAPT->sV[currElement].newEdgesAbove = 0;
                    eAPT->sV[currElement].newEdgesBelow = 0;
                    tree->vArr[currElement].sigma = INFINITY_MY;
                    eAPT->sV[currElement].newSigma = 0;
                    STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement) {
                        int64_t k = STINGER_EDGE_DEST;

                        if (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 1, -2)) {
                            touchedVerticesUp[__atomic_fetch_add(&tvUpEnd, 1, __ATOMIC_RELAXED)] =  k;
                            Queue[__atomic_fetch_add(qEnd_nxt, 1, __ATOMIC_RELAXED)] = k;
                        }
                    }
                    STINGER_FORALL_EDGES_OF_VTX_END();
                }
            }
            *qStart = *qStart_nxt;
            *qEnd = *qEnd_nxt;
            *qStart_nxt = *qEnd;
            *qEnd_nxt = *qStart_nxt;
        }
    }

    for (int64_t q = 0; q < tvDownEnd; q++) {
        int64_t k = touchedVerticesDown[q];

        if (eAPT->sV[k].touched > 0) {
            tree->vArr[k].sigma = eAPT->sV[k].newSigma;
        }

        if (eAPT->sV[k].touched != 0) {
            tree->vArr[k].delta = eAPT->sV[k].newDelta;
        }

        tree->vArr[k].edgesAbove = eAPT->sV[k].newEdgesAbove;
        eAPT->sV[k].diffPath = 0;
        eAPT->sV[k].touched = 0;
        eAPT->sV[k].newDelta = 0.0;
        eAPT->sV[k].newSigma = 0;
        eAPT->sV[k].newLevel = 0;
        eAPT->sV[k].newEdgesAbove = 0;
    }

    for (int64_t q = 0; q < tvUpEnd; q++) {
        int64_t k = touchedVerticesUp[q];

        if (eAPT->sV[k].touched > 0) {
            tree->vArr[k].sigma = eAPT->sV[k].newSigma;
        }

        if (eAPT->sV[k].touched != 0) {
            tree->vArr[k].delta = eAPT->sV[k].newDelta;
        }

        tree->vArr[k].edgesBelow = eAPT->sV[k].newEdgesBelow;
        eAPT->sV[k].diffPath = 0;
        eAPT->sV[k].touched = 0;
        eAPT->sV[k].newDelta = 0.0;
        eAPT->sV[k].newSigma = 0;
        eAPT->sV[k].newLevel = 0;

        eAPT->sV[k].newEdgesAbove = 0;
        eAPT->sV[k].newEdgesBelow = 0;
    }

    eAPT->qStart = 0;
    eAPT->qEnd = 0;
    eAPT->qStart_nxt = 0;
    eAPT->qEnd_nxt = 0;

    eAPT->tqStart = 0;
    eAPT->tqEnd = 0;
    eAPT->tqStart_nxt = 0;
    eAPT->tqEnd_nxt = 0;

    eAPT->qStartSame = 0;
    eAPT->qEndSame = 0;
    eAPT->qStartSame_nxt = 0;
    eAPT->qEndSame_nxt = 0;
}

