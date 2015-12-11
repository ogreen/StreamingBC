#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#include "streaming_utils.h"
// #include "streaming_utils_brandes.h"
#include "streaming_deletions.h"
#include "timing_util.h"

#include "bcTreeDS.h"

#include "list.h"
#include "dsUtils.h"

#define MAX(X, Y) ((X) > (Y) ? X : Y)

void deleteEdgeWithoutMovement(bcForest* forest, struct stinger *sStinger,
                    uint64_t currRoot, uint64_t startVertex, uint64_t parentVertex,
                    uint64_t deletedPathsToRoot, extraArraysPerThread *eAPT)
{
    bcTree *tree = forest->forest[currRoot];

    uint64_t *QueueDown = eAPT->QueueDown;
    uint64_t *QueueUp = eAPT->QueueUp;

    int64_t NV = forest->NV;

    eAPT->sV[startVertex].newPathsToRoot = tree->vArr[startVertex].pathsToRoot;
    eAPT->sV[startVertex].touched = 1;
    eAPT->sV[startVertex].diffPath = deletedPathsToRoot;
    eAPT->sV[startVertex].newPathsToRoot += eAPT->sV[startVertex].diffPath;
    
    QueueDown[0] = startVertex;
    int64_t qStart = 0, qEnd = 1;

    int64_t deepestLevel = tree->vArr[startVertex].level;
    int64_t initialLevel = tree->vArr[startVertex].level;

    // Starting BFS descent from "startVertex" down to all vertices that have shortest paths through "startVertex" to "currRoot".
    // All elements that will be touched will received a positive value in their touched field.
    // In this implementaiton, "STACKS" are not used for the dependency accumulation stage. A multi-level queue is used instad.
    // Each level in the tree has a queue and a counter specifying the size of the queue. (?)
    // For simplicity, all elements are pushed into both the multi-level queue and into the regular queue for the BFS traversal.
    while (qStart < qEnd)
    {
        uint64_t currElement = QueueDown[qStart];
        int64_t currElementLevel = tree->vArr[currElement].level;
        int64_t currElementTouched = eAPT->sV[currElement].touched;
        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement)
        {
            uint64_t neighbor = STINGER_EDGE_DEST;

            if (currElementLevel + 1 == tree->vArr[neighbor].level)
            {
                if (eAPT->sV[neighbor].touched == 0)
                {
                    deepestLevel = MAX(deepestLevel, tree->vArr[neighbor].level);                
                    eAPT->sV[neighbor].newPathsToRoot = tree->vArr[currElement].pathsToRoot;
                    QueueDown[qEnd++] = neighbor;
                    eAPT->sV[neighbor].touched = currElementLevel + 1;
                    eAPT->sV[neighbor].newPathsToRoot += eAPT->sV[currElement].diffPath;
                    eAPT->sV[neighbor].diffPath = eAPT->sV[currElement].diffPath;

                }
                else if (eAPT->sV[neighbor].touched == currElementLevel + 1)
                {
                    eAPT->sV[neighbor].diffPath += eAPT->sV[currElement].diffPath;
                    eAPT->sV[neighbor].newPathsToRoot -= eAPT->sV[currElement].diffPath;
                }
            }
        }
        STINGER_FORALL_EDGES_OF_VTX_END();

        #if COUNT_TRAVERSALS==1
        eAPT->dynamicTraverseEdgeCounter += stinger_typed_outdegree(sStinger, currElement, 0);
        eAPT->dynamicTraverseVerticeCounter++;
        #endif
        
        qStart++;       
    }

    int64_t QUpStart = 0, QUpEnd = 0;
    int64_t currElement;    
    int64_t qDownEndMarker = --qEnd; 
 
    // Starting the multi-level dependency accumulation / "BFS" ascent.
    // The ascent continues going up as long as the root has not been reached and that there
    // are elements in the current level of the ascent. The ascent starts in the deepest level
    // of the graph.
    // It is worth noting that in the ascent stage:
    // 1) All the previously untouched elements that are touched are marked with "-1".
    // 2) On the way up, it is possible that elements that were touched in the BFS descent will
    //    touch elements that were not touched in the descent and that are below startVertex. These
    //    are elements that do not have shortest paths going through startVertex, though their BC
    //    values still changed due to the changes occuring below them. Because of this, they are
    //    placed in the multi-level queue.
    while (1)
    {
        //printf("qEnd, QUpStart, QUpEnd:  %d, %d, %d\n", qEnd, QUpStart, QUpEnd);
        if (qEnd < 0 && QUpStart >= QUpEnd)
        {
            break;
        }
        else if (qEnd < 0)
        {
            currElement = QueueUp[QUpStart++];
        }
        else if(qEnd >= 0)
        {
            if (QUpEnd > 0)
            {
/*                printf("before\n");
                if(QueueUp==NULL || QueueDown==NULL) 
                    printf("baa\n");
                printf("QUpStart: %d\n", QUpStart);
                printf("QueueUp[QUpStart]: %d\n", QueueUp[QUpStart]);
                printf("qEnd: %d\n", qEnd);
                printf("QueueDown[qEnd]: %d\n", QueueDown[qEnd]);*/
                if (tree->vArr[QueueUp[QUpStart]].level < tree->vArr[QueueDown[qEnd]].level)
                {
                    currElement = QueueDown[qEnd--];
                }
                else 
                {
                    currElement = QueueUp[QUpStart++];
                }
            }
            else
            {
                currElement = QueueDown[qEnd--];
            }
        }

        int64_t currElementLevel = tree->vArr[currElement].level;
        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement)
        {
            uint64_t neighbor = STINGER_EDGE_DEST;

            if (tree->vArr[neighbor].level = currElementLevel - 1)
            {
                if (eAPT->sV[neighbor].touched == 0)
                {
                    eAPT->sV[neighbor].newDelta = tree->vArr[neighbor].delta;
                    eAPT->sV[neighbor].touched = -1;
                    QueueUp[QUpEnd++] = neighbor;
                    eAPT->sV[neighbor].newPathsToRoot = tree->vArr[neighbor].pathsToRoot; // Inspect this line more.                             
                }

                eAPT->sV[neighbor].newDelta += ((bc_t)eAPT->sV[neighbor].newPathsToRoot / 
                                                (bc_t)eAPT->sV[currElement].newPathsToRoot) * 
                                                (bc_t)(eAPT->sV[currElement].newDelta + 1);

                if (eAPT->sV[neighbor].touched == -1 && (neighbor != parentVertex || currElement != startVertex))
                {
                    eAPT->sV[neighbor].newDelta -= ((bc_t)tree->vArr[neighbor].pathsToRoot / 
                                                    (bc_t)tree->vArr[neighbor].pathsToRoot) * 
                                                    ((bc_t)tree->vArr[currElement].delta + 1);
                }
            }
        }
        STINGER_FORALL_EDGES_OF_VTX_END();

        #if COUNT_TRAVERSALS==1
        eAPT->dynamicTraverseEdgeCounter += stinger_typed_outdegree(sStinger, currElement, 0);
        eAPT->dynamicTraverseVerticeCounter++;
        #endif

        if (currElement != currRoot)
        {
            eAPT->sV[currElement].totalBC += eAPT->sV[currElement].newDelta - tree->vArr[currElement].delta;
        }
    }

    for (uint64_t i = 0; i <= qDownEndMarker; i++)
    {
        uint64_t k = QueueDown[i];
        tree->vArr[k].delta = eAPT->sV[k].newDelta;
        tree->vArr[k].pathsToRoot = eAPT->sV[k].newPathsToRoot;
        eAPT->sV[k].diffPath = 0;
        eAPT->sV[k].touched = 0;
        eAPT->sV[k].newDelta = 0.0;
        eAPT->sV[k].newPathsToRoot = 0;
    }
    for (uint64_t i = 0; i < QUpEnd; i++)
    {
        uint64_t k = QueueUp[i];
        tree->vArr[k].delta = eAPT->sV[k].newDelta;
        eAPT->sV[k].diffPath = 0;
        eAPT->sV[k].touched = 0;
        eAPT->sV[k].newDelta = 0.0;
        eAPT->sV[k].newPathsToRoot = 0;
    }
}

void moveDownTreeBrandes(bcForest *forest, struct stinger *sStinger,
                    uint64_t currRoot, uint64_t startVertex, uint64_t parentVertex,
                    uint64_t prevDist, extraArraysPerThread *eAPT)
{
    bcTree *tree = forest->forest[currRoot];
    
    int64_t NV = forest->NV;
    
    uint64_t *QueueDown = eAPT->QueueDown;
    uint64_t *QueueUp = eAPT->QueueUp;
    uint64_t *QueueSame = eAPT->QueueSame;

    int64_t newLevel[NV]; // Might include newLevel in eAPT in the future.

    for (uint64_t i = 0; i < NV; i++)
    {
        newLevel[i] = 0;
        eAPT->sV[i].newPathsToRoot = tree->vArr[i].pathsToRoot;
    }

    QueueDown[0] = startVertex;
    int64_t qDownStart = 0, qDownEnd = 1, qSameStart = 0, qSameEnd = 0;

    int64_t stopLevel = tree->vArr[startVertex].level;
    
    newLevel[startVertex] = INFINITY_MY;
    eAPT->sV[startVertex].newPathsToRoot = INFINITY_MY;
    eAPT->sV[startVertex].newDelta = 0;
    eAPT->sV[startVertex].touched = 1;

    int64_t deepestLevel = stopLevel;
    
    list_ptr *multiLevelQueues = eAPT->multiLevelQueues;

    // Sets newLevel to tree level for all siblings of startVertex and infinity for all other neighbors of startVertex.
    // Puts siblingts of startVertex in QueueSame and every other neighbor in QueueDown.
    STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, startVertex)
    {
        int64_t neighbor = STINGER_EDGE_DEST;

        if (tree->vArr[neighbor].level == tree->vArr[startVertex].level)
        {
            newLevel[neighbor] = tree->vArr[neighbor].level;
            eAPT->sV[neighbor].newPathsToRoot = tree->vArr[neighbor].pathsToRoot;
            QueueSame[qSameStart++] = neighbor;
        }
        else
        {
            newLevel[neighbor] = INFINITY_MY;
            eAPT->sV[neighbor].newPathsToRoot = INFINITY_MY;
            deepestLevel = MAX(deepestLevel, tree->vArr[neighbor].level);
        }

        QueueDown[qDownEnd++] = neighbor;
        eAPT->sV[neighbor].touched = 1;
        eAPT->sV[neighbor].newDelta = 0;
    }   
    STINGER_FORALL_EDGES_OF_VTX_END();

    qDownStart = 1;
    while (qDownStart != qDownEnd)
    {
        int64_t currElement = QueueDown[qDownStart++];
        
        // Adds every vertex on same level as startVertex with at least one child to QueueSame.
        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement)
        {
            int64_t neighbor = STINGER_EDGE_DEST;

            if (eAPT->sV[neighbor].touched == 0)
            {
                if (tree->vArr[neighbor].level == stopLevel && tree->vArr[currElement].level == stopLevel + 1)
                {
                    newLevel[neighbor] = tree->vArr[neighbor].level;
                    eAPT->sV[neighbor].newPathsToRoot = tree->vArr[neighbor].pathsToRoot;
                    QueueSame[qSameEnd++] = neighbor;
                    eAPT->sV[neighbor].touched = 1;
                    QueueDown[qDownEnd++] = neighbor;
                    eAPT->sV[neighbor].newDelta = 0;
                }
                else if (tree->vArr[neighbor].level > stopLevel)
                {
                    newLevel[neighbor] = INFINITY_MY;
                    eAPT->sV[neighbor].newPathsToRoot = INFINITY_MY;
                    eAPT->sV[neighbor].touched = 1;
                    QueueDown[qDownEnd++] = neighbor;
                    eAPT->sV[neighbor].newDelta = 0;
                }
            }
        }
        STINGER_FORALL_EDGES_OF_VTX_END();
    }

    qDownEnd = qSameEnd;
    qDownStart = 0;
    for (int64_t i = 0; i < qDownEnd; i++)
    {
        QueueDown[i] = QueueSame[i];
        eAPT->sV[i].newDelta = 0;
        //newLevel[parentVertex] = tree->vArr[parentVertex].level; Commented out on 12/9/15 5:50pm
        //eAPT->sV[i].newPathsToRoot = tree->vArr[parentVertex].pathsToRoot; Commented out on 12/9/15 5:50pm
        append(multiLevelQueues[stopLevel], makeNode(QueueSame[i]));
    }

    eAPT->sV[parentVertex].touched = -1;
    newLevel[parentVertex] = tree->vArr[parentVertex].level;
    eAPT->sV[parentVertex].newPathsToRoot = tree->vArr[parentVertex].pathsToRoot;
    
    append(multiLevelQueues[stopLevel - 1], makeNode(parentVertex));
    eAPT->sV[parentVertex].newDelta = tree->vArr[parentVertex].delta - 
                            ((bc_t)tree->vArr[parentVertex].pathsToRoot / (bc_t)tree->vArr[startVertex].pathsToRoot) *
                            (bc_t)(tree->vArr[startVertex].delta + 1);
   
    // Starting BFS descent. 
    // While QueueDown is not empty.
    while (qDownStart != qDownEnd)
    {
        uint64_t currElement = QueueDown[qDownStart++];

        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement)
        {
            uint64_t neighbor = STINGER_EDGE_DEST;
            
            // If this is a neighbor and has not been found.
            if (newLevel[neighbor] > newLevel[currElement])
            {   
                // Checking if "neighbor" has been found.
                if (newLevel[neighbor] == INFINITY_MY)
                {
                    newLevel[neighbor] = newLevel[currElement] + 1;
                    QueueDown[qDownEnd++] = neighbor;
                    eAPT->sV[neighbor].newDelta = 0;
                    deepestLevel = MAX(deepestLevel, newLevel[neighbor]);
                    append(multiLevelQueues[newLevel[neighbor]], makeNode(neighbor));       
                    if (eAPT->sV[neighbor].touched == 0)
                        printf("ooooooooooooooooooooooooooooooooooops\n");
                }
                
                if (eAPT->sV[neighbor].newPathsToRoot == INFINITY_MY)
                {
                    // "neighbor" has not been found and therefore its paths to the root are through its parent.
                    eAPT->sV[neighbor].newPathsToRoot = eAPT->sV[currElement].newPathsToRoot;
                }
                else
                {
                    // "neighbor" has been found and has multiple paths to the root through its multiple parents.
                    eAPT->sV[neighbor].newPathsToRoot += eAPT->sV[currElement].newPathsToRoot;
                }        
            }
        }
        STINGER_FORALL_EDGES_OF_VTX_END();

#if COUNT_TRAVERSALS == 1
        eAPT->dynamicTraverseEdgeCounter += stinger_typed_outdegree(sStinger, currElement, 0);
        eAPT->dynamicTraverseVerticeCounter++;
#endif

    }

    for (int64_t i = 0; i < NV; i++)
    {
        if (eAPT->sV[i].touched != 0)
        {
            tree->vArr[i].level = newLevel[i];
        }
    }

    // Starting Multi-Level "BFS" ascent.
    //
    // The ascent continues going up as long as the root has not been reached and that there
    // are elements in the current level of the ascent. The ascent starts in the deepest level
    // of the grpah.
    // It is worth noting that in the ascent stage:
    // 1) All previously untouched elements that are touched are makred with "-1".
    // 2) On the way up, it is possible that elements that were touched in the BFS descent will
    //    touch elements that were not touched in the descent that that are below "startVertex". these
    //    are elements that do not have shortest paths going through "startVertex". However, their BC
    //    vlaues have changed due to the the changes occurring below them. Because of this, they are
    //    placed in the multi-level queue.
   
    while (deepestLevel >= 0 && multiLevelQueues[deepestLevel]->size > 0) 
    {
        node_t *currQueue = getFirst(multiLevelQueues[deepestLevel]);
        while (currQueue != NULL)
        {
            uint64_t currElement = currQueue->id;
            
            STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement)
            {
                uint64_t neighbor = STINGER_EDGE_DEST;

                if (tree->vArr[neighbor].level == tree->vArr[currElement].level - 1)
                {
                    if (eAPT->sV[neighbor].touched == 0)
                    {
                        eAPT->sV[neighbor].newDelta = tree->vArr[neighbor].delta;
                        eAPT->sV[neighbor].touched = -1;
                        append(multiLevelQueues[tree->vArr[neighbor].level], makeNode(neighbor));
                    }

                    eAPT->sV[neighbor].newDelta += 
                                ((bc_t)eAPT->sV[neighbor].newPathsToRoot / (bc_t)eAPT->sV[currElement].newPathsToRoot) *
                                (bc_t)(eAPT->sV[currElement].newDelta + 1);

                    if (eAPT->sV[neighbor].touched < 0)
                    {
                        eAPT->sV[neighbor].newDelta -=
                                ((bc_t)tree->vArr[neighbor].pathsToRoot / (bc_t)tree->vArr[currElement].pathsToRoot) *
                                (bc_t)(tree->vArr[currElement].delta + 1);
                    }
                }
            }
            STINGER_FORALL_EDGES_OF_VTX_END();
            
#if COUNT_TRAVERSALS == 1
            eAPT->dynamicTraverseEdgeCounter += stinger_typed_outdegree(sStinger, currElement, 0);
            eAPT->dynamicTraverseVerticeCounter++;
#endif
            if (currElement != currRoot)
            {
                eAPT->sV[currElement].totalBC += eAPT->sV[currElement].newDelta - tree->vArr[currElement].delta;
            }
            currQueue = currQueue->next;
        }
        deepestLevel--;
    }

    for (int64_t i = 0; i < NV; i++)
    {
        if (eAPT->sV[i].touched > 0)
        {
            tree->vArr[i].pathsToRoot = eAPT->sV[i].newPathsToRoot;
        }
        if (eAPT->sV[i].touched != 0)
        {
            tree->vArr[i].delta = eAPT->sV[i].newDelta;
        }
    }
}


/* 
void moveDownTreeBrandes(bcForeset *forest, struct stinger *sStinger,
                    uint64_t currRoot, uint64_t startVertex, uint64_t parentVertex,
                    uint64_t prevDist, extraArraysPerThreat *eAPT)
{
    bcTree *tree = forest->forest[currRoot];
    
    int64_t NV = forest->NV;
    uint64_t *QueueDown = eAPT->QueueDown;
    uint64_t *QueueUp = eAPT->QueueUp;
    uint64_t *QueueSame = eAPT->QueueSame;

    list_ptr *multiLevelQueues = eAPT->multiLevelQueues;

    for (uint64_t i = 0; i < NV; i++)
    {
        eAPT->sV[i].newPathsToRoot = tree->vArr[i].pathsToRoot;
    }

    eAPT->sV[parentVertex].newPathsToRoot = tree->vArr[parentVertex].pathsToRoot;
    eAPT->sV[startVertex].newPathsToRoot = tree->vArr[startVertex].pathsToRoot;

    int64_t qStart = 0;
    int64_t qEnd = 1;

    QueueDown[0] = startVertex;
    eAPT->sV[startVertex].touched = 1;
    eAPT->sV[startVertex].newPathsToRoot = eAPT->sV[parentVertex].newPathsToRoot;
    eAPT->sV[startVertex].diffPath = eAPT->sV[parentVertex].newPathsToRoot;
    eAPT->sV[startVertex].movementDelta = prevDist;
    eAPT->sV[startVertex].IMoved = 1;

    int64_t deepestLevel = tree->vArr[parentIndex].level + 1;

    while (qStart != qEnd)
    {
        int64_t currElement = QueueDown[qStart];
        int64_t currElementTouched = eAPT->sV[currElement].touched;
        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement)
        {
            uint64_t neighbor = STINGER_EDGE_DEST;

            if (eAPT->sV[neighbor].touched == 0)
            {
                int64_t computedDelta = eAPT->sV[currElement].movementDelta - (tree->vArr[currElement].level - tree->vArr[neighbor].level + 1); // Not sure about this line.
                eAPT->sV[neighbor].touched = currElementTouched + 1;
                
                if (computedDelta > 0)
                {
                    eAPT->sV[neighbor].newPathsToRoot = eAPT->sV[currElement].diffPath;
                    eAPT->sV[neighbor].diffPath = eAPT->sV[currElement].diffPath;
                    eAPT->sV[neighbor].movementDelta = computedDelta;
                    eAPT->sV[neighbor].IMoved = 1;
                    QueueDown[qEnd++] = neighbor;
                }
                else if (computedDelta == 0)
                {
                    eAPT->sV[neighbor].diffPath += eAPT->sV[currElement].diffPath;
                    eAPT->sV[neighbor].newPathsToRoot -= eAPT->sV[currElement].diffPath;
                    eAPT->sV[neighbor].IMoved = 0;
                    QueueDown[qEnd++] = neighbor;    
                }
                else
                {
                    eAPT->sV[neighbor].touched = 0; // Not sure about this line.
                }
            }
            else if (eAPT->sV[neighbor].touched == currElementTouched + 1)
            {
                int64_t computedDelta = eAPT->sV[currElement].movementDelta - (tree->vArr[currElement].level - tree->vArr[neighbor].level + 1);
                if (computedDelta >= 0)
                {
                    eAPT->sV[neighbor].diffpath += eAPT->sV[currElement].diffPath;
                    eAPT->sV[neighbor].newPathsToRoot -= eAPT->sV[currElement].diffPath;
                }
            }
        }
        STINGER_FORALL_EDGES_OF_VTX_END();
        
        #if COUNT_TRAVERSALS==1
        eAPT->dynamicTraverseEdgeCounter += stinger_typed_outdegree(sStinger, currElement, 0);
        eAPT->dynamicTraverseVerticeCounter++;
        #endif

        tree->vArr[currElement].level -= eAPT->sV[currElement].movementDelta;
        deepestLevel = MAX(deepestLevel, tree->vArr[currElement].level);
        qStart++;
    }

    // Starting dependency accumulation ascent.
    qEnd = 0;
    node_t *temp_node;
    for (int lev = tree->vArr[startVertex].level; lev < NV; lev++)
    {
        while(multiLevelQueues[lev]->size > 0)
        {
            temp_node = getFirst(multiLevelQueues[lev]);
            QueueDown[qEnd++] = temp_node->id;
            deleteFirst(multiLevelQueues[lev]);
        }
    }
    qEnd--;

    int64_t qDownEndMarker = qEnd;
    int64_t QUpStart = 0, QDownStart = 0;
    int64_t QSameStart = 0, QSameEnd = 0;
    int64_t currElement = 0;
    int64_t upCounter = 0;
    int64_t depthDown = tree->vArr[QueueDown[qEnd]].level;
    int64_t depthUp = -1, depthSame = -1;
    
    // The ascent continues going up as long as the root has not been reached and that there
    // are elements in the current level of the ascent. The ascent starts in the deepest level
    // of the graph.
    // It's worth noting that in the ascent stage:
    // 1) All previously untouched elements that are touched are marked with "-1".
    // 2) On the way up, it is possible that elements that were touched in the BFS descent will
    //    touch elements that were not touched in the descent that are belwo "startVErtex". These
    //    are elements that do not have shortest paths going through "startVertex", however their BC
    //    values have changed due to the changes occurring below them. Because of this, they are
    //    placed in the multilevel queue.
    // 3) There are vertices that did not move and that one of their neighbors move down (such that
    //    the vertices are no longer in the same level ). Consequently, the number of shortest paths going
    //    through the vertex that did not move was increased. These vertices will be touched as -2
    //    and added to the queue and the ascent will continue from these vertices as well.
    while (1)
    {
        if (depthDown == -1 && depthSame == -1 && depthUp == -1)
            break;

        if (depthDown >= depthSame && depthDown >= depthUp)
        {
            currElement = QueueDown[qEnd--];
            if (qEnd < 0)
                depthDown = -1;
            else
                depthDown = tree->vArr[QueueDown[qEnd]].level;
        }
        else if (depthUp >= depthSame && depthUp >= depthDown)
        {
            currElement = QueueUp[QUpStart++];
            if (QUpStart >= QUpEnd)
                depthSame = -1;
            else
                depthSame = tree->vArr[QueueUp[QUpStart]].level;
        }
        else if (depthSame >= depthUp && depthSame >= depthDown)
        {
            currElement = QueueDown[QDownStart++];
            if (QDownStart >= QDownEnd)
                depthSame = -1;
            else
                depthSame = tree->vArr[QueueDown[QDownStart]].level;
        }
    
        int64_t currElementLevel = tree->vArr[currElement].level;
        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement)
        {
            uint64_t neighbor = STINGER_EDGE_DEST;
            if (tree->vArr[neighbor].level == currElementLevel - 1)
            {
                if(eAPT->sV[neighbor].touched == 0)
                {
                    eAPT->sV[neighbor].newDelta = tree->vArr[neighbor].delta;
                    upCounter++;
                    eAPT->sV[neighbor].touched = -1;
                    QueueUp[QUpEnd++] = neighbor;
                    
                    if (depthUp == -1)
                        depthUp = tree->vArr[QueueUp[QUpStart]].level;
                    
                    // Not sure about this line.        
                    eAPT->sV[neighbor].newPathsToRoot = tree->vArr[neighbor].pathsToRoot;
                    
                }

                eAPT->sV[neighbor].newDelta += 
                    ((bc_t)eAPT->sV[neighbor].newPathsToRoot / (bc_t)eAPT->sV[currElement].newPathsToRoot) *
                    (1 + eAPT->sV[currElement].newDelta);

                if (eAPT->sV[neighbor].touched < 0 && (neighbor != parentVertex || currElement != startVertex))
                {
                    
                }
            }
        }
        STINGER_FORALL_EDGES_OF_VTX_END();

    }
}*/
