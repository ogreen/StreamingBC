#include <omp.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#include "streamingbc_aux.h"

#define PARENT_ANCHORED 3
#define SIBLING_ANCHORED 4

void addEdgeWithoutMovementBrandes(bcForest * forest, struct stinger * sStinger,
                                   int64_t currRoot, int64_t startVertex, int64_t parentVertex,
                                   int64_t addedPathsToRoot,  extraArraysPerThread * eAPT)
{

    bcTree * tree = forest->forest[currRoot];
    int64_t * QueueDown = eAPT->QueueDown;
    int64_t * QueueUp = eAPT->QueueUp;
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

    *qEnd = 1;
    *qStart_nxt = 1;
    *qEnd_nxt = 1;

    int64_t deepestLevel = tree->vArr[startVertex].level;
    int64_t intialLevel = tree->vArr[startVertex].level;

    // Starting BFS decent from "startVertex", down to all the vertices that have shortest paths through "startVertex"
    // All elements that will be touched will receive a positive value in their touched field.
    // In this implementation, "STACKS" are not used for the "moving up" stage. Rather, a multi-level queue is used.
    // Each level in the tree(max depth NV) has a queue and a counter specifiying how deep a specific deepth-queue is.
    // For simplicity, all elements are pushed both into the multi-level queue and into the regular queue which is used
    // for the BFS traversal.
    while (*qStart < *qEnd) {
        int64_t currElement = QueueDown[*qStart];
        int64_t levelCurrPlusOne = tree->vArr[currElement].level + 1;
        int64_t touchedCurrPlusOne = eAPT->sV[currElement].touched + 1;

        if (currElement != startVertex) {
            eAPT->sV[currElement].newEdgesAbove = tree->vArr[currElement].edgesAbove;
        }

        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement) {
            int64_t k = STINGER_EDGE_DEST;

            // if this vertex has not been added yet
            if (currElement != startVertex
                    && tree->vArr[currElement].level >= tree->vArr[startVertex].level
                    && tree->vArr[k].level < tree->vArr[currElement].level) {

                if (eAPT->sV[k].touched != 0) {
                    eAPT->sV[currElement].newEdgesAbove -= tree->vArr[k].edgesAbove;
                    eAPT->sV[currElement].newEdgesAbove += eAPT->sV[k].newEdgesAbove;
                }
            }

            if (levelCurrPlusOne == (tree->vArr[k].level)) {

                if (eAPT->sV[k].touched == 0) {
                    // Checking if a "deeper level" has been reached.
                    if (deepestLevel < tree->vArr[k].level) {
                        deepestLevel = tree->vArr[k].level;
                    }

                    eAPT->sV[k].newSigma += tree->vArr[k].sigma;
                    // insert this vertex into the BFS queue
                    QueueDown[(*qEnd_nxt)++] = k;
                    // indicate that it is in the next level of the BFS
                    eAPT->sV[k].touched += touchedCurrPlusOne;
                    // add new paths to root that go through current BFS Vertex
                    eAPT->sV[k].newSigma += eAPT->sV[currElement].diffPath;
                    // pass on my new paths to root for its search
                    eAPT->sV[k].diffPath += eAPT->sV[currElement].diffPath;

                } else {
                    // otherwise if it has been touched, but is specifically in the next level
                    // of the search (meaning it has more than one edge to the current level)

                    // add new paths to root that go through current BFS Vertex
                    eAPT->sV[k].newSigma += eAPT->sV[currElement].diffPath;
                    // pass on my new paths to root for its search
                    eAPT->sV[k].diffPath += eAPT->sV[currElement].diffPath;
                }
            }
        }
        STINGER_FORALL_EDGES_OF_VTX_END();

#if COUNT_TRAVERSALS==1
        eAPT->dynamicTraverseEdgeCounter += stinger_typed_outdegree(sStinger, currElement, 0);
        eAPT->dynamicTraverseVerticeCounter++;
#endif
        (*qStart)++;

        if (*qStart == *qEnd) {
            *qStart = *qStart_nxt;
            *qEnd = *qEnd_nxt;
            *qStart_nxt = *qEnd;
            *qEnd_nxt = *qStart_nxt;
        }
    }

    int64_t QUpStart = 0, QUpEnd = 0;
    int64_t currElement;
    (*qEnd)--;

    int64_t qDownEndMarker = *qEnd;

    // Starting Multi-Level "BFS" ascent.
    // The ascent continues going up as long as the root has not been reached and that there
    // are elements in the current level of the ascent. The ascent starts in the deepest level
    // of the graph.
    // It was worth noting that in the ascent stage:
    // 1) All previously untouched elements that are touched are marked with "-1".
    // 2) On the way up, it is possible that elements that were touched in the BFS decent, will
    //    touch elements that were not touchded in the decsent and that are below "vertex". These
    //    are elements that do not have shortest paths going through "vertex" ,however, there BC
    //    values have changed due to the changes occuring below them. Because of this, they are
    //    placed in the Multi-level queue.
    while (1) {
        if (*qEnd < 0 && QUpStart >= QUpEnd) {
            break;
        } else if (*qEnd < 0) {
            currElement = QueueUp[QUpStart++];
        } else if (*qEnd >= 0) {
            if (QUpEnd > 0) {
                if (tree->vArr[QueueUp[QUpStart]].level < tree->vArr[QueueDown[*qEnd]].level) {
                    currElement = QueueDown[(*qEnd)--];
                } else {
                    currElement = QueueUp[QUpStart++];
                }
            } else {
                currElement = QueueDown[(*qEnd)--];
            }
        }

        int64_t levelCurrMinusOne = tree->vArr[currElement].level - 1;

        if (currElement != parentVertex) {
            eAPT->sV[currElement].newEdgesBelow = tree->vArr[currElement].edgesBelow;
        }

        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement) {
            int64_t k = STINGER_EDGE_DEST;

            if (currElement != parentVertex
                    && tree->vArr[currElement].level <= tree->vArr[parentVertex].level
                    && tree->vArr[k].level > tree->vArr[currElement].level) {

                if (eAPT->sV[k].touched != 0) {
                    eAPT->sV[currElement].newEdgesBelow -= tree->vArr[k].edgesBelow;
                    eAPT->sV[currElement].newEdgesBelow += eAPT->sV[k].newEdgesBelow;
                }
            }

            if (tree->vArr[k].level == levelCurrMinusOne) {
                // Checking to see if "k" has been touched before.
                if (eAPT->sV[k].touched == 0) {
                    eAPT->sV[k].newDelta += tree->vArr[k].delta;
                    // Marking element as touched in the ascent stage.
                    eAPT->sV[k].touched = -1;
                    QueueUp[QUpEnd] = k;
                    QUpEnd++;
                    eAPT->sV[k].newSigma += tree->vArr[k].sigma;
                }

                eAPT->sV[k].newDelta +=
                    ((bc_t)eAPT->sV[k].newSigma / (bc_t)eAPT->sV[currElement].newSigma) *
                    (bc_t)(eAPT->sV[currElement].newDelta + 1);

                // For the elements that are touched in the ascent stage it is necessary to
                // to reduce the values that they previously had.
                // In addition to this, the "parentVertex" that is connected to "vertex", i.e.
                // the vertices of the new edge, needs to increase its betweenness centrality
                // following the new connection, without removing the old delta value.
                if (eAPT->sV[k].touched < 0 && (k != parentVertex  || currElement != startVertex)) {
                    eAPT->sV[k].newDelta -=
                        ((bc_t)tree->vArr[k].sigma / (bc_t)tree->vArr[currElement].sigma) *
                        (bc_t)(tree->vArr[currElement].delta + 1);
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
        eAPT->sV[k].newEdgesBelow = 0;
    }

    eAPT->sV[startVertex].newEdgesAbove = 0;
    eAPT->sV[parentVertex].newEdgesAbove = 0;

    for (int64_t c = 0; c < QUpEnd; c++) {
        int64_t k = QueueUp[c];
        tree->vArr[k].delta = eAPT->sV[k].newDelta;
        tree->vArr[k].edgesBelow = eAPT->sV[k].newEdgesBelow;
        eAPT->sV[k].diffPath = 0;
        eAPT->sV[k].touched = 0;
        eAPT->sV[k].newDelta = 0.0;
        eAPT->sV[k].newEdgesAbove = 0;
        eAPT->sV[k].newEdgesBelow = 0;
        eAPT->sV[k].newSigma = 0;
    }

    eAPT->sV[startVertex].newEdgesBelow = 0;
    eAPT->sV[parentVertex].newEdgesBelow = 0;

    eAPT->qStart = 0;
    eAPT->qEnd = 0;
    eAPT->qStart_nxt = 0;
    eAPT->qEnd_nxt = 0;
}


void moveUpTreeBrandes(bcForest * forest, struct stinger * sStinger,
                       int64_t currRoot, int64_t startVertex, int64_t parentVertex,
                       int64_t prevDist,  extraArraysPerThread * eAPT)
{

    bcTree * tree = forest->forest[currRoot];
    int64_t NV = forest->NV;
    int64_t * QueueDown = eAPT->QueueDown;
    int64_t * QueueUp = eAPT->QueueUp;
    int64_t * QueueSame = eAPT->QueueSame;

    list_ptr * multiLevelQueues = eAPT->multiLevelQueues;
    queue_t * queue = eAPT->queue;
    level_node_t * levelIndices = eAPT->levelIndices;

    eAPT->sV[parentVertex].newSigma = tree->vArr[parentVertex].sigma;
    eAPT->sV[startVertex].newSigma = tree->vArr[startVertex].sigma;

    int64_t * qStart = &(eAPT->qStart);
    int64_t * qEnd = &(eAPT->qEnd);
    int64_t * qStart_nxt = &(eAPT->qStart_nxt);
    int64_t * qEnd_nxt = &(eAPT->qEnd_nxt);

    *qEnd = 1;
    *qStart_nxt = 1;
    *qEnd_nxt = 1;

    QueueDown[0] = startVertex;
    eAPT->sV[startVertex].touched = 1;
    eAPT->sV[startVertex].newSigma = eAPT->sV[parentVertex].newSigma;
    eAPT->sV[startVertex].diffPath = eAPT->sV[parentVertex].newSigma;
    eAPT->sV[startVertex].movementDelta = prevDist;
    eAPT->sV[startVertex].IMoved = 1;

    eAPT->sV[startVertex].newEdgesAbove = eAPT->sV[parentVertex].newEdgesAbove + 1;
    int64_t deepestLevel = tree->vArr[parentVertex].level + 1;

    // Starting BFS decent from "startVertex", down to all the vertices that have shortest paths through "startVertex"
    // All elements that will be touched will receive a positive value in their touched field.
    // In this implementation, "STACKS" are not used for the "moving up" stage. Rather, a multi-level queue is used.
    // Each level in the tree(max depth NV) has a queue and a counter specifiying how deep a specific deepth-queue is.
    // For simplicity, all elements are pushed both into the multi-level queue and into the regular queue which is used
    // for the BFS traversal.
    while (*qStart != *qEnd) {
        int64_t currElement = QueueDown[*qStart];

        int64_t touchedCurrPlusOne = eAPT->sV[currElement].touched + 1;
        eAPT->sV[currElement].newEdgesAbove = 0;
        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement) {
            int64_t k = STINGER_EDGE_DEST;

            int64_t computedDelta = eAPT->sV[currElement].movementDelta -
                                    (tree->vArr[currElement].level - tree->vArr[k].level + 1);

            int64_t newCurrLevel = tree->vArr[currElement].level - eAPT->sV[currElement].movementDelta;
            int64_t newKLevel = tree->vArr[k].level - computedDelta;


            if (computedDelta >= 0 && newKLevel < newCurrLevel) {
                if (eAPT->sV[k].touched == 0) {
                    eAPT->sV[currElement].newEdgesAbove += tree->vArr[k].edgesAbove + 1;
                } else {
                    eAPT->sV[currElement].newEdgesAbove += eAPT->sV[k].newEdgesAbove + 1;
                }
            } else if (computedDelta < 0 && tree->vArr[k].level < newCurrLevel) {
                if (eAPT->sV[k].touched == 0) {
                    eAPT->sV[currElement].newEdgesAbove += tree->vArr[k].edgesAbove + 1;
                } else {
                    eAPT->sV[currElement].newEdgesAbove += eAPT->sV[k].newEdgesAbove + 1;
                }
            }

            if (eAPT->sV[k].touched == 0) {
                // compute distance the adjacent vertex should be moved
                int64_t computedDelta = eAPT->sV[currElement].movementDelta -
                                        (tree->vArr[currElement].level - tree->vArr[k].level + 1);

                eAPT->sV[k].touched = touchedCurrPlusOne;

                // if the adjacent vertex should be moved, put it in the queue
                if (computedDelta > 0) {
                    eAPT->sV[k].newSigma = eAPT->sV[currElement].diffPath;
                    eAPT->sV[k].diffPath = eAPT->sV[currElement].diffPath;
                    eAPT->sV[k].movementDelta = computedDelta;
                    eAPT->sV[k].IMoved = 1;
                    QueueDown[(*qEnd_nxt)++] = k;
                } else if (computedDelta == 0) {
                    // Vertex that will not be moved has been found.
                    eAPT->sV[k].newSigma = tree->vArr[k].sigma;
                    eAPT->sV[k].newSigma += eAPT->sV[currElement].diffPath;
                    eAPT->sV[k].diffPath += eAPT->sV[currElement].diffPath;
                    eAPT->sV[k].movementDelta = computedDelta;
                    eAPT->sV[k].IMoved = 0;
                    QueueDown[(*qEnd_nxt)++] = k;
                } else {
                    // Vertex that the number of shortest path to the root does not change has been found.
                    // This vertex is not marked as it might be touched on the way up.
                    eAPT->sV[k].touched = 0;
                }

                // if adjacent and in the next level
            } else if (eAPT->sV[k].touched == touchedCurrPlusOne) {
                int64_t computedDelta = eAPT->sV[currElement].movementDelta - (tree->vArr[currElement].level - tree->vArr[k].level + 1);

                if (computedDelta >= 0) {
                    eAPT->sV[k].newSigma += eAPT->sV[currElement].diffPath;
                    eAPT->sV[k].diffPath += eAPT->sV[currElement].diffPath;
                }
            }

        }
        STINGER_FORALL_EDGES_OF_VTX_END();


#if COUNT_TRAVERSALS==1
        eAPT->dynamicTraverseEdgeCounter += stinger_typed_outdegree(sStinger, currElement, 0);
        eAPT->dynamicTraverseVerticeCounter++;
#endif
        // move ourself and retire
        tree->vArr[currElement].level -= eAPT->sV[currElement].movementDelta;
        appendDS(queue, levelIndices, tree->vArr[currElement].level, currElement);

        // Checking if a "deeper level" has been reached.
        if (deepestLevel < tree->vArr[currElement].level) {
            deepestLevel = tree->vArr[currElement].level;
        }

        (*qStart)++;

        if (*qStart == *qEnd) {
            *qStart = *qStart_nxt;
            *qEnd = *qEnd_nxt;
            *qStart_nxt = *qEnd;
            *qEnd_nxt = *qStart_nxt;
        }
    }

    // Starting Multi-Level "BFS" ascent.
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

    (*qEnd)--;

    int64_t qDownEndMarker = *qEnd;

    int64_t QUpStart = 0, QUpEnd = 0;
    int64_t QSameStart = 0, QSameEnd = 0;
    int64_t currElement = 0;  //dummy initilization - variable will be initialized in function.
    int64_t upCounter = 0;
    int64_t depthDown = tree->vArr[QueueDown[*qEnd]].level, depthUp = -1, depthSame = - 1;

    // The ascent continues going up as long as the root has not been reached and that there
    // are elements in the current level of the ascent. The ascent starts in the deepest level
    // of the graph.
    // It was worth noting that in the ascent stage:
    // 1) All previously untouched elements that are touched are marked with "-1".
    // 2) On the way up, it is possible that elements that were touched in the BFS decent, will
    //    touch elements that were not touchded in the decsent and that are below "startVertex". These
    //    are elements that do not have shortest paths going through "startVertex" ,however, there BC
    //    values have changed due to the changes occuring below them. Because of this, they are
    //    placed in the Multi-level queue.
    // 3) There are vertices that did not move and that one of their neighbors move up(such that
    //    the vertices are now in the same level). Consequently, the number of shortest path going
    //    through the vertex that did not move was reduced. These vertices will be touched  as -2
    //    and added to the queue and the "BFS ascent" will continue from these vertices as well.

    while (1) {
        if (depthDown == -1 && depthSame == -1 && depthUp == -1) {
            break;
        }

        if (depthDown >= depthSame && depthDown >= depthUp) {
            currElement = QueueDown[*qEnd];
            (*qEnd)--;

            if (*qEnd < 0) {
                depthDown = -1;
            } else {
                depthDown = tree->vArr[QueueDown[*qEnd]].level;
            }
        } else if (depthUp >= depthSame && depthUp >= depthDown) {
            currElement = QueueUp[QUpStart];
            QUpStart++;

            if (QUpStart >= QUpEnd) {
                depthUp = -1;
            } else {
                depthUp = tree->vArr[QueueUp[QUpStart]].level;
            }
        } else if (depthDown <= depthSame && depthUp <= depthSame) {
            currElement = QueueSame[QSameStart];
            QSameStart++;

            if (QSameStart >= QSameEnd) {
                depthSame = -1;
            } else {
                depthSame = tree->vArr[QueueSame[QSameStart]].level;
            }
        }


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
                if (eAPT->sV[k].touched == 0) {
                    eAPT->sV[k].newDelta = tree->vArr[k].delta;
                    upCounter++;
                    // Marking element as touched in the ascent stage.
                    eAPT->sV[k].touched = -1;

                    QueueUp[QUpEnd] = k;

                    if (depthUp == -1) {
                        depthUp = tree->vArr[QueueUp[QUpStart]].level;
                    }

                    QUpEnd++;
                    eAPT->sV[k].newSigma = tree->vArr[k].sigma;
                }

                eAPT->sV[k].newDelta +=
                    ((bc_t)eAPT->sV[k].newSigma / (bc_t)eAPT->sV[currElement].newSigma) *
                    (bc_t)(eAPT->sV[currElement].newDelta + 1);

                // For the elements that are touched in the ascent stage it is necessary to
                // to reduce the values that they previously had.
                // In addition to this, the "parentVertex" that is connected to "startVertex", i.e.
                // the vertices of the new edge, needs to increase its betweenness centrality
                // following the new connection, without removing the old delta value.
                if (eAPT->sV[k].touched < 0 && (k != parentVertex  || currElement != startVertex)) {
                    eAPT->sV[k].newDelta -=
                        ((bc_t)tree->vArr[k].sigma / (bc_t)tree->vArr[currElement].sigma) *
                        (bc_t)(tree->vArr[currElement].delta + 1);
                }
            } else if (tree->vArr[k].level == tree->vArr[currElement].level && ((eAPT->sV[currElement].IMoved == 1 && eAPT->sV[k].IMoved < 0) )) {
                // Vertices that did not move and that one of their neighbors move up(such that
                // the vertices are now in the same level).
                // Checking to see if "k" has been touched before.
                if (eAPT->sV[k].touched == 0) {
                    eAPT->sV[k].newDelta = tree->vArr[k].delta;
                    upCounter++;
                    // Marking element as touched in the ascent stage.
                    eAPT->sV[k].touched = -2;
                    QueueSame[QSameEnd] = k;

                    if (depthSame == -1) {
                        depthSame = tree->vArr[QueueSame[QSameStart]].level;
                    }

                    QSameEnd++;
                    eAPT->sV[k].newSigma = tree->vArr[k].sigma;
                }

                // Paths that previosul went through this vertex no longer go through them, thus the
                // shortest path count(BC) is reduced.
                eAPT->sV[k].newDelta -=
                    ((bc_t)tree->vArr[k].sigma / (bc_t)tree->vArr[currElement].sigma) *
                    (bc_t)(tree->vArr[currElement].delta + 1);
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


    for (int64_t c = 0; c < QSameEnd; c++) {
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

    for (int64_t c = 0; c < QUpEnd; c++) {
        int64_t k = QueueUp[c];
        tree->vArr[k].delta = eAPT->sV[k].newDelta;
        tree->vArr[k].edgesBelow = eAPT->sV[k].newEdgesBelow;
        eAPT->sV[k].diffPath = 0;
        eAPT->sV[k].touched = 0;
        eAPT->sV[k].newDelta = 0.0;
        eAPT->sV[k].movementDelta = 0;
        eAPT->sV[k].IMoved = -1;
        eAPT->sV[k].newEdgesBelow = 0;
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


// Case 2
void removeEdgeWithoutMovementBrandes(bcForest * forest, struct stinger * sStinger, int64_t currRoot,
                                      int64_t startVertex, int64_t parentVertex, int64_t deletedPathsFromRoot,
                                      extraArraysPerThread * eAPT)
{

    bcTree * tree = forest->forest[currRoot];

    int64_t NV = forest->NV;
    int64_t * Queue = eAPT->QueueSame;
    int64_t * QueueDown = eAPT->QueueDown;
    int64_t * QueueUp = eAPT->QueueUp;

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
        int64_t currElement = QueueDown[*qDownStart];

        if (currElement != startVertex)
            eAPT->sV[currElement].newEdgesAbove = tree->vArr[currElement].edgesAbove;

        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement) {
            int64_t k = STINGER_EDGE_DEST;

            if (currElement != startVertex
                    && tree->vArr[currElement].level - 1 == tree->vArr[k].level
                    && tree->vArr[currElement].level >= tree->vArr[startVertex].level) {

                if (eAPT->sV[k].touched == 0)
                    eAPT->sV[k].newEdgesAbove = tree->vArr[k].edgesAbove;

                eAPT->sV[currElement].newEdgesAbove -= tree->vArr[k].edgesAbove;
                eAPT->sV[currElement].newEdgesAbove += eAPT->sV[k].newEdgesAbove;
            }

            // if this vertex has not been added yet
            if ((tree->vArr[currElement].level + 1) == (tree->vArr[k].level)) {
                if (eAPT->sV[k].touched == 0) {
                    // Checking if a "deeper level" has been reached.
                    if (deepestLevel < tree->vArr[k].level)
                        deepestLevel = tree->vArr[k].level;

                    // insert this vertex into the BFS queue
                    QueueDown[(*qDownEnd_nxt)++] = k;

                    eAPT->sV[k].newSigma = tree->vArr[k].sigma;
                    // indicate that it is in the next level of the BFS
                    eAPT->sV[k].touched = eAPT->sV[currElement].touched + 1;
                    // add new paths to root that go through current BFS Vertex
                    eAPT->sV[k].newSigma -= eAPT->sV[currElement].diffPath;
                    // pass on my new paths to root for its search
                    eAPT->sV[k].diffPath = eAPT->sV[currElement].diffPath;
                }
                // otherwise if it has been touched, but is specifically in the next level
                // of the search (meaning it has more than one edge to the current level)
                else if (eAPT->sV[k].touched == eAPT->sV[currElement].touched + 1) {
                    // add new paths to root that go through current BFS Vertex
                    eAPT->sV[k].newSigma -= eAPT->sV[currElement].diffPath;
                    // pass on my new paths to root for its search
                    eAPT->sV[k].diffPath += eAPT->sV[currElement].diffPath;
                }
            }
        }
        STINGER_FORALL_EDGES_OF_VTX_END();

#if COUNT_TRAVERSALS==1
        eAPT->dynamicTraverseEdgeCounter += stinger_typed_outdegree(sStinger, currElement, 0);
        eAPT->dynamicTraverseVerticeCounter++;
#endif
        (*qDownStart)++;

        if (*qDownStart == *qDownEnd) {
            *qDownStart = *qDownStart_nxt;
            *qDownEnd = *qDownEnd_nxt;
            *qDownStart_nxt = *qDownEnd;
            *qDownEnd_nxt = *qDownStart_nxt;
        }
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
    //    touch elements that were not touchded in the decsent and that are below "vertex". These
    //    are elements that do not have shortest paths going through "vertex" ,however, there BC
    //    values have changed due to the changes occuring below them. Because of this, they are
    //    placed in the Multi-level queue.

    while (1) {

        int64_t currElement = -1;

        if (*qDownEnd < 0 && qUpStart >= qUpEnd) {
            break;
        } else if (*qDownEnd < 0) {
            currElement = QueueUp[qUpStart++];
        } else if (*qDownEnd >= 0) {
            if (qUpEnd > 0) {
                if (tree->vArr[QueueUp[qUpStart]].level < tree->vArr[QueueDown[*qDownEnd]].level) {
                    currElement = QueueDown[(*qDownEnd)--];
                } else {
                    currElement = QueueUp[qUpStart++];
                }
            } else {
                currElement = QueueDown[(*qDownEnd)--];
            }
        } 

        if (currElement != parentVertex) {
            eAPT->sV[currElement].newEdgesBelow = tree->vArr[currElement].edgesBelow;
        }

        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement) {
            int64_t k = STINGER_EDGE_DEST;


            if (currElement != parentVertex
                    && tree->vArr[currElement].level <= tree->vArr[parentVertex].level
                    && tree->vArr[k].level > tree->vArr[currElement].level) {

                if (eAPT->sV[k].touched == 0)
                    eAPT->sV[k].newEdgesBelow = tree->vArr[k].edgesBelow;

                eAPT->sV[currElement].newEdgesBelow -= tree->vArr[k].edgesBelow;
                eAPT->sV[currElement].newEdgesBelow += eAPT->sV[k].newEdgesBelow;
            }

            if (eAPT->sV[parentVertex].touched != -1 && tree->vArr[k].level == tree->vArr[parentVertex].level) {
                eAPT->sV[parentVertex].touched = -1;
                QueueUp[qUpEnd++] = parentVertex;
                eAPT->sV[parentVertex].newSigma = tree->vArr[parentVertex].sigma;
                eAPT->sV[parentVertex].newDelta = tree->vArr[parentVertex].delta -
                                                  ((bc_t)tree->vArr[parentVertex].sigma / (bc_t)tree->vArr[startVertex].sigma) *
                                                  (bc_t)(tree->vArr[startVertex].delta + 1);
            }

            // Checking that the vertices are in different levels.
            if (tree->vArr[k].level == (tree->vArr[currElement].level - 1)) {
                // Checking to see if "k" has been touched before.
                if (eAPT->sV[k].touched == 0) {
                    eAPT->sV[k].newDelta = tree->vArr[k].delta;
                    // Marking element as touched in the ascent stage.
                    eAPT->sV[k].touched = -1;
                    eAPT->sV[k].newSigma = tree->vArr[k].sigma;
                    QueueUp[qUpEnd++] = k;
                }

                eAPT->sV[k].newDelta +=
                    ((bc_t)eAPT->sV[k].newSigma / (bc_t)eAPT->sV[currElement].newSigma) *
                    (bc_t)(eAPT->sV[currElement].newDelta + 1);

                // For the elements that are touched in the ascent stage it is necessary to
                // to reduce the values that they previously had.
                // In addition to this, the "parentVertex" that is connected to "vertex", i.e.
                // the vertices of the new edge, needs to increase its betweenness centrality
                // following the new connection, without removing the old delta value.
                if (eAPT->sV[k].touched < 0 && (k != parentVertex  || currElement != startVertex)) {
                    eAPT->sV[k].newDelta -=
                        ((bc_t)tree->vArr[k].sigma / (bc_t)tree->vArr[currElement].sigma) *
                        (bc_t)(tree->vArr[currElement].delta + 1);
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

    for (int64_t q = 0; q < qUpEnd; q++) {
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

int compare_levels(const void * a, const void * b, void * arg)
{
    int64_t * u = (int64_t *) a;
    int64_t * v = (int64_t *) b;
    extraArraysPerThread * eAPT = (extraArraysPerThread *) arg;

    return eAPT->sV[*u].newLevel - eAPT->sV[*v].newLevel;
}

// Case 3
void moveDownTreeBrandes(bcForest * forest, struct stinger * sStinger, int64_t currRoot,
                         int64_t startVertex, int64_t parentVertex, extraArraysPerThread * eAPT)
{
    bcTree * tree = forest->forest[currRoot];

    int64_t NV = forest->NV;
    int64_t * Queue = eAPT->QueueDown;
    int64_t * QueueUp = eAPT->QueueUp;
    int64_t * topQueue = eAPT->QueueSame;

    int64_t * touchedVerticesDown = eAPT->touchedVerticesDown;
    int64_t * touchedVerticesUp = eAPT->touchedVerticesUp;
    queue_t * queue = eAPT->queue;
    level_node_t * levelIndices = eAPT->levelIndices;

    Queue[0] = startVertex;
    int64_t tqStart = 0, tqEnd = 0, tvDownEnd = 0, tvUpEnd = 0;
    int64_t stopLevel = tree->vArr[startVertex].level;

    int64_t * qStart = &(eAPT->qStart);
    int64_t * qEnd = &(eAPT->qEnd);
    int64_t * qStart_nxt = &(eAPT->qStart_nxt);
    int64_t * qEnd_nxt = &(eAPT->qEnd_nxt);

    *qEnd = 1;
    *qStart_nxt = 1;
    *qEnd_nxt = 1;

    eAPT->sV[startVertex].newLevel = INFINITY_MY;
    eAPT->sV[startVertex].newSigma = INFINITY_MY;
    eAPT->sV[startVertex].newDelta = 0.0;
    eAPT->sV[startVertex].newEdgesAbove = INFINITY_MY;
    eAPT->sV[startVertex].touched = 1;
    touchedVerticesDown[tvDownEnd++] = startVertex;

    int64_t deepestLevel = stopLevel;

    *qStart = 0;

    while (*qStart != *qEnd) {
        int64_t currElement = Queue[(*qStart)++];

        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement) {
            int64_t k = STINGER_EDGE_DEST;

            if (tree->vArr[k].level == tree->vArr[currElement].level + 1 && eAPT->sV[k].touched == 0) {
                eAPT->sV[k].newEdgesAbove = INFINITY_MY;
                eAPT->sV[k].newLevel = INFINITY_MY;
                eAPT->sV[k].newSigma = INFINITY_MY;
                touchedVerticesDown[tvDownEnd++] = k;
                Queue[(*qEnd)++] = k;
                eAPT->sV[k].newDelta = 0.0;

            }
        }
        STINGER_FORALL_EDGES_OF_VTX_END();

        int64_t parentOutsideSubtree = 0;
        int64_t siblingOutsideSubtree = 0;
        int64_t parentPathsToRoot = 0;
        int64_t siblingPathsToRoot = 0;

        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement) {
            int64_t l = STINGER_EDGE_DEST;

            if (tree->vArr[l].level == tree->vArr[currElement].level - 1 && eAPT->sV[l].touched == 0) {
                parentOutsideSubtree = l;
                parentPathsToRoot += tree->vArr[l].sigma;
            } else if (tree->vArr[l].level == tree->vArr[currElement].level && eAPT->sV[l].touched == 0) {
                siblingOutsideSubtree = l;
                siblingPathsToRoot += tree->vArr[l].sigma;
            }
        }
        STINGER_FORALL_EDGES_OF_VTX_END();

        if (parentOutsideSubtree) {

            if (eAPT->sV[currElement].touched == 1 || eAPT->sV[currElement].touched == SIBLING_ANCHORED) {
                topQueue[tqEnd++] = currElement;
                eAPT->sV[currElement].newLevel = tree->vArr[parentOutsideSubtree].level + 1;
            }

            eAPT->sV[currElement].touched = PARENT_ANCHORED;
            eAPT->sV[currElement].newDelta = 0.0;
        } else if (siblingOutsideSubtree) {

            if (eAPT->sV[currElement].touched == 1) {
                topQueue[tqEnd++] = currElement;
                eAPT->sV[currElement].newLevel = tree->vArr[siblingOutsideSubtree].level + 1;
            }

            if (eAPT->sV[currElement].touched != PARENT_ANCHORED) {
                eAPT->sV[currElement].touched = SIBLING_ANCHORED;
            }

            eAPT->sV[currElement].newDelta = 0.0;
        }

        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement) {
            int64_t k = STINGER_EDGE_DEST;

            if (tree->vArr[k].level == tree->vArr[currElement].level + 1 && eAPT->sV[k].touched == 0) {
                if (eAPT->sV[currElement].touched == PARENT_ANCHORED) {
                    eAPT->sV[k].touched = PARENT_ANCHORED;
                    eAPT->sV[k].newLevel = eAPT->sV[currElement].newLevel + 1;
                } else if (eAPT->sV[currElement].touched == SIBLING_ANCHORED && eAPT->sV[k].touched != PARENT_ANCHORED) {
                    eAPT->sV[k].touched = SIBLING_ANCHORED;
                    eAPT->sV[k].newLevel = eAPT->sV[currElement].newLevel + 1;
                } else {
                    eAPT->sV[k].touched = 1;
                }
            }
        }
        STINGER_FORALL_EDGES_OF_VTX_END();


    }

    *qEnd = 1;
    *qStart = 0;
    tqStart = 0;

    int64_t key, j;

    // Insertion sort on topQueue by levels of the vertices in the queue.
    for (int64_t i = 1; i < tqEnd; i++) {
        key = topQueue[i];
        j = i - 1;

        while (j >= 0 && eAPT->sV[topQueue[j]].newLevel > eAPT->sV[key].newLevel) {
            topQueue[j + 1] = topQueue[j];
            j = j - 1;
        }

        topQueue[j + 1] = key;
    }

    if (tqEnd != 0) {
        appendDS(queue, levelIndices, eAPT->sV[topQueue[tqStart]].newLevel, topQueue[tqStart]);
        Queue[0] = topQueue[tqStart++];
        eAPT->sV[Queue[0]].touched = 5;
    } else {
        Queue[0] = startVertex;
        eAPT->sV[Queue[0]].touched = 5;
        *qStart = 0;
        *qEnd = 1;
    }

    // While queue is not empty
    while (*qStart != *qEnd) {
        int64_t currElement = Queue[(*qStart)++];

        while (tqStart != tqEnd && (eAPT->sV[currElement].newLevel >= eAPT->sV[topQueue[tqStart]].newLevel)) {
            if (eAPT->sV[topQueue[tqStart]].touched != 5) {
                Queue[(*qEnd_nxt)++] = topQueue[tqStart];
                eAPT->sV[topQueue[tqStart]].touched = 5;
                appendDS(queue, levelIndices, eAPT->sV[topQueue[tqStart]].newLevel, topQueue[tqStart]);
            }

            tqStart++;
        }

        eAPT->sV[currElement].newEdgesAbove = 0;

        if (deepestLevel < eAPT->sV[currElement].newLevel) {
            deepestLevel = eAPT->sV[currElement].newLevel;
        }

        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement) {
            int64_t k = STINGER_EDGE_DEST;

            if (eAPT->sV[currElement].newLevel == tree->vArr[k].level + 1 && eAPT->sV[k].touched == 0) {
                if (eAPT->sV[currElement].newSigma == INFINITY_MY) {
                    // k has not been found and therefore its paths to the roots are through its parent.
                    eAPT->sV[currElement].newSigma = tree->vArr[k].sigma;
                } else {
                    // k has been found and has multiple paths to the root as it has multiple parents.
                    eAPT->sV[currElement].newSigma += tree->vArr[k].sigma;
                }

                if (eAPT->sV[currElement].newEdgesAbove == INFINITY_MY) {
                    eAPT->sV[currElement].newEdgesAbove = tree->vArr[k].edgesAbove + 1;
                } else {
                    eAPT->sV[currElement].newEdgesAbove += tree->vArr[k].edgesAbove + 1;
                }
            } else if (eAPT->sV[currElement].newLevel == eAPT->sV[k].newLevel + 1 && eAPT->sV[k].touched != 0) {
                if (eAPT->sV[currElement].newSigma == INFINITY_MY) {
                    eAPT->sV[currElement].newSigma = eAPT->sV[k].newSigma;
                } else {
                    eAPT->sV[currElement].newSigma += eAPT->sV[k].newSigma;
                }

                if (eAPT->sV[currElement].newEdgesAbove == INFINITY_MY) {
                    eAPT->sV[currElement].newEdgesAbove = eAPT->sV[k].newEdgesAbove + 1;
                } else {
                    eAPT->sV[currElement].newEdgesAbove += eAPT->sV[k].newEdgesAbove + 1;
                }
            }

        }
        STINGER_FORALL_EDGES_OF_VTX_END();

        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement) {
            int64_t k = STINGER_EDGE_DEST;

            if (eAPT->sV[k].newEdgesAbove == 0) {
                eAPT->sV[k].newEdgesAbove = tree->vArr[k].edgesAbove;
            }

            if (eAPT->sV[k].newLevel > eAPT->sV[currElement].newLevel) {
                // Checking if "k" has been found.

                if (eAPT->sV[k].newLevel == INFINITY_MY) {
                    eAPT->sV[k].newLevel = eAPT->sV[currElement].newLevel + 1;
                } else if (eAPT->sV[currElement].newLevel + 1 < eAPT->sV[k].newLevel) {
                    eAPT->sV[k].newLevel = eAPT->sV[currElement].newLevel + 1;
                }

                if (eAPT->sV[k].touched != 5) {
                    eAPT->sV[k].touched = 5;
                    Queue[(*qEnd_nxt)++] = k;
                    eAPT->sV[k].newDelta = 0.0;

                    if (deepestLevel < eAPT->sV[k].newLevel) {
                        deepestLevel = eAPT->sV[k].newLevel;
                    }

                    appendDS(queue, levelIndices, eAPT->sV[k].newLevel, k);
                }
            }
        }


        STINGER_FORALL_EDGES_OF_VTX_END();
#if COUNT_TRAVERSALS==1
        eAPT->dynamicTraverseEdgeCounter += stinger_typed_outdegree(sStinger, currElement, 0);
        eAPT->dynamicTraverseVerticeCounter++;
#endif
        int dequed = 0;

        while (!dequed && tqStart != tqEnd && *qStart == *qEnd && *qStart_nxt == *qEnd_nxt) {
            if (eAPT->sV[topQueue[tqStart]].touched != 5) {
                dequed = 1;
                Queue[(*qEnd_nxt)++] = topQueue[tqStart];
                eAPT->sV[topQueue[tqStart]].touched = 5;
                appendDS(queue, levelIndices, eAPT->sV[topQueue[tqStart]].newLevel, topQueue[tqStart]);
            }

            tqStart++;
        }


        if (*qStart == *qEnd) {
            *qStart = *qStart_nxt;
            *qEnd = *qEnd_nxt;
            *qStart_nxt = *qEnd;
            *qEnd_nxt = *qStart_nxt;
        }
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

    int64_t qUpStart = 0, qUpEnd = 0;
    (*qEnd)--;

    if (eAPT->sV[startVertex].newLevel == INFINITY_MY) {
        if (eAPT->sV[parentVertex].touched == 0) {
            eAPT->sV[parentVertex].newDelta = tree->vArr[parentVertex].delta -
                                              ((bc_t)tree->vArr[parentVertex].sigma / (bc_t)tree->vArr[startVertex].sigma) *
                                              (bc_t)(tree->vArr[startVertex].delta + 1);

            eAPT->sV[parentVertex].newSigma = tree->vArr[parentVertex].sigma;
            eAPT->sV[parentVertex].touched = -1;
            QueueUp[qUpEnd++] = parentVertex;
        }
    }

    // Starting Multi-Level "BFS" ascent.

    // The ascent continues going up as long as the root has not been reached and that there
    // are elements in the current level of the ascent. The ascent starts in the deepest level
    // of the graph.
    // It was worth noting that in the ascent stage:
    // 1) All previously untouched elements that are touched are marked with "-1".
    // 2) On the way up, it is possible that elements that were touched in the BFS decent, will
    //    touch elements that were not touchded in the decsent and that are below "vertex". These
    //    are elements that do not have shortest paths going through "vertex" ,however, there BC
    //    values have changed due to the changes occuring below them. Because of this, they are
    //    placed in the Multi-level queue.


    while (1) {
        // Removing last element from the queue
        int64_t currElement = -1;

        if (*qEnd < 0 && qUpStart >= qUpEnd) {
            break;
        } else if (*qEnd < 0) {
            currElement = QueueUp[qUpStart++];
        } else if (qEnd >= 0) {
            if (qUpEnd > 0) {
                int64_t qUpStartLevel = -1;

                if (eAPT->sV[QueueUp[qUpStart]].newLevel == 0) {
                    qUpStartLevel = tree->vArr[QueueUp[qUpStart]].level;
                } else {
                    qUpStartLevel = eAPT->sV[QueueUp[qUpStart]].newLevel;
                }

                int64_t qEndLevel = -1;

                if (eAPT->sV[Queue[*qEnd]].newLevel == 0) {
                    qEndLevel = tree->vArr[Queue[*qEnd]].level;
                } else {
                    qEndLevel = eAPT->sV[Queue[*qEnd]].newLevel;
                }

                if (qUpStartLevel < qEndLevel) {
                    currElement = Queue[(*qEnd)--];
                } else {
                    currElement = QueueUp[qUpStart++];
                }
            } else {
                currElement = Queue[(*qEnd)--];
            }
        }

        eAPT->sV[currElement].newEdgesBelow = 0;
        touchedVerticesUp[tvUpEnd++] = currElement;

        int64_t currElementLevel = -1;

        if (eAPT->sV[currElement].newLevel == 0) {
            currElementLevel = tree->vArr[currElement].level;
        } else {
            currElementLevel = eAPT->sV[currElement].newLevel;
        }

        int64_t parentVertexLevel = -1;

        if (eAPT->sV[parentVertex].newLevel == 0) {
            parentVertexLevel = tree->vArr[parentVertex].level;
        } else {
            parentVertexLevel = eAPT->sV[parentVertex].newLevel;
        }

        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement) {
            int64_t k = STINGER_EDGE_DEST;

            int64_t kLevel = -1;

            if (eAPT->sV[k].newLevel == 0) {
                kLevel = tree->vArr[k].level;
            } else {
                kLevel = eAPT->sV[k].newLevel;
            }

            if (kLevel == currElementLevel + 1) {
                if (eAPT->sV[k].touched == 0)
                    eAPT->sV[k].newEdgesBelow = tree->vArr[k].edgesBelow;

                eAPT->sV[currElement].newEdgesBelow += eAPT->sV[k].newEdgesBelow + 1;
            }

            if (kLevel == parentVertexLevel) {
                if (eAPT->sV[parentVertex].touched == 0) {
                    eAPT->sV[parentVertex].newDelta = tree->vArr[parentVertex].delta -
                                                      ((bc_t)tree->vArr[parentVertex].sigma / (bc_t)tree->vArr[startVertex].sigma) *
                                                      (bc_t)(tree->vArr[startVertex].delta + 1);

                    eAPT->sV[parentVertex].newSigma = tree->vArr[parentVertex].sigma;
                    eAPT->sV[parentVertex].touched = -1;
                    QueueUp[qUpEnd++] = parentVertex;
                }
            }

            // Checking that the vertices are in different levels.
            if (kLevel == currElementLevel - 1) {

                // Checking to see if "k" has been touched before.
                if (eAPT->sV[k].touched == 0) {
                    eAPT->sV[k].newDelta = tree->vArr[k].delta;

                    // Marking element as touched in the ascent stage.
                    eAPT->sV[k].newSigma = tree->vArr[k].sigma;
                    eAPT->sV[k].touched = -1;
                    QueueUp[qUpEnd++] = k;
                }

                eAPT->sV[k].newDelta +=
                    ((bc_t)eAPT->sV[k].newSigma / (bc_t)eAPT->sV[currElement].newSigma) *
                    (bc_t)(eAPT->sV[currElement].newDelta + 1);


                // For the elements that are touched in the ascent stage it is necessary to
                // to reduce the values that they previously had.
                // In addition to this, the "parentVertex" that is connected to "vertex", i.e.
                // the vertices of the new edge, needs to increase its betweenness centrality
                // following the new connection, without removing the old delta value.
                if (eAPT->sV[k].touched < 0 && tree->vArr[k].level < tree->vArr[currElement].level) {
                    eAPT->sV[k].newDelta -=
                        ((bc_t)tree->vArr[k].sigma / (bc_t)tree->vArr[currElement].sigma) *
                        (bc_t)(tree->vArr[currElement].delta + 1);
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

    for (int64_t k = 0; k < tvDownEnd; k++) {
        int64_t vertex = touchedVerticesDown[k];
        tree->vArr[vertex].level = eAPT->sV[vertex].newLevel;
    }

    // Handles case where edge deletion creates new connected component.
    if (tree->vArr[startVertex].level == INFINITY_MY) {
        *qStart = 0;
        *qEnd = 1;
        Queue[0] = startVertex;

        eAPT->sV[startVertex].touched = -2;

        while (*qStart != *qEnd) {
            int64_t currElement = Queue[(*qStart)++];

            eAPT->sV[currElement].totalBC -= tree->vArr[currElement].delta;

            tree->vArr[currElement].edgesBelow = 0;
            tree->vArr[currElement].edgesAbove = 0;
            eAPT->sV[currElement].newEdgesAbove = 0;
            eAPT->sV[currElement].newEdgesBelow = 0;
            tree->vArr[currElement].sigma = INFINITY_MY;
            eAPT->sV[currElement].newSigma = 0;
            STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement) {
                int64_t k = STINGER_EDGE_DEST;

                if (eAPT->sV[k].touched != -2) {
                    eAPT->sV[k].touched = -2;
                    touchedVerticesUp[tvUpEnd++] = k;
                    Queue[(*qEnd)++] = k;
                }
            }
            STINGER_FORALL_EDGES_OF_VTX_END();
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
}
