

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#include "streaming_utils.h"
//#include "streaming_utils_brandes.h"
#include "timing_util.h"

#include "bcTreeDS.h"

#include "list.h"
#include "dsUtils.h"

// Case 2
void removeEdgeWithoutMovementBrandes(bcForest* forest, struct stinger* sStinger, uint64_t currRoot, 
                        uint64_t startVertex, uint64_t parentVertex, uint64_t deletedPathsFromRoot, 
                        extraArraysPerThread *eAPT)
{
    bcTree* tree = forest->forest[currRoot];

    int64_t NV = forest->NV;
    uint64_t Queue[NV];

    uint64_t queueBFSTREE[NV][NV];
    int64_t levelCounter[NV];

    for(uint64_t k = 0; k < NV; k++)
    {
        eAPT->sV[k].touched = 0;
        eAPT->sV[k].diffPath = 0;
        eAPT->sV[k].newPathsToRoot = tree->vArr[k].pathsToRoot;
        eAPT->sV[k].newDelta = 0.0;
        levelCounter[k] = 0;
        
        for (uint64_t l = 0; l < NV; l++)
            queueBFSTREE[k][l] = 0;
    }


    eAPT->sV[startVertex].touched = 1;
    eAPT->sV[startVertex].newPathsToRoot -= deletedPathsFromRoot;
    eAPT->sV[startVertex].diffPath = deletedPathsFromRoot;

    Queue[0] = startVertex;
    int64_t qStart=0,qEnd=1;

    int64_t deepestLevel = tree->vArr[startVertex].level;
    //printf("currRoot, startVertex_level: %ld, %ld\n", currRoot, tree->vArr[startVertex].level);
    queueBFSTREE[tree->vArr[startVertex].level][levelCounter[tree->vArr[startVertex].level]]=startVertex;
    levelCounter[tree->vArr[startVertex].level]++;


    // Starting BFS decent from "startVertex", down to all the vertices that have shortest paths through "startVertex"
    // All elements that will be touched will receive a positive value in their touched field.
    // In this implementation, "STACKS" are not used for the "moving up" stage. Rather, a multi-level queue is used.
    // Each level in the tree(max depth NV) has a queue and a counter specifiying how deep a specific deepth-queue is.
    // For simplicity, all elements are pushed both into the multi-level queue and into the regular queue which is used
    // for the BFS traversal.
    while(qStart!=qEnd)
    {
        uint64_t currElement = Queue[qStart];

        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
        {
            uint64_t k = STINGER_EDGE_DEST;

	    // if this vertex has not been added yet
            if ((tree->vArr[currElement].level + 1) == (tree->vArr[k].level))
            {
	        if(eAPT->sV[k].touched == 0)
		{
		    // Checking if a "deeper level" has been reached.
		    if(deepestLevel < tree->vArr[k].level)
                        deepestLevel = tree->vArr[k].level;

		    // insert this vertex into the BFS queue
		    Queue[qEnd++] = k;
                    queueBFSTREE[tree->vArr[k].level][levelCounter[tree->vArr[k].level]] = k;
                    levelCounter[tree->vArr[k].level]++;

		    // indicate that it is in the next level of the BFS
		    eAPT->sV[k].touched = eAPT->sV[currElement].touched + 1;
		    // add new paths to root that go through current BFS Vertex
		    eAPT->sV[k].newPathsToRoot -= eAPT->sV[currElement].diffPath;
		    // pass on my new paths to root for its search
		    eAPT->sV[k].diffPath = eAPT->sV[currElement].diffPath;

		}
	        
                // otherwise if it has been touched, but is specifically in the next level
		// of the search (meaning it has more than one edge to the current level)
		else if(eAPT->sV[k].touched == eAPT->sV[currElement].touched + 1)
		{
		    // add new paths to root that go through current BFS Vertex
		    eAPT->sV[k].newPathsToRoot -= eAPT->sV[currElement].diffPath;
		    // pass on my new paths to root for its search
		    eAPT->sV[k].diffPath += eAPT->sV[currElement].diffPath;
		}     
            }

        }
        STINGER_FORALL_EDGES_OF_VTX_END();

	qStart++;
    }

    // The parent vertex needs to be placed in the queue for the dependency accumulation stage.
    // Also, it no longer has a child and so the delta from the child needs to be removed.
    queueBFSTREE[tree->vArr[parentVertex].level][levelCounter[tree->vArr[parentVertex].level]]=parentVertex;
    levelCounter[tree->vArr[parentVertex].level]++;
    eAPT->sV[parentVertex].newDelta = tree->vArr[parentVertex].delta -
                ((bc_t)tree->vArr[parentVertex].pathsToRoot / (bc_t)tree->vArr[startVertex].pathsToRoot) *
                (bc_t)(tree->vArr[startVertex].delta + 1);

    eAPT->sV[parentVertex].touched = -1;

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
    while (deepestLevel >=0 && levelCounter[deepestLevel] > 0)
    {

        int64_t currQueue = 0;

        while (currQueue < levelCounter[deepestLevel])
        {
            // Removing last element from the queue
            uint64_t currElement = queueBFSTREE[deepestLevel][currQueue];
            currQueue++;

            STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
            {
                uint64_t k = STINGER_EDGE_DEST;


                // Checking that the vertices are in different levels.
                if (tree->vArr[k].level == (tree->vArr[currElement].level - 1))
                {
                    // Checking to see if "k" has been touched before.
                    if (eAPT->sV[k].touched == 0)
                    {
                        eAPT->sV[k].newDelta=tree->vArr[k].delta;

                        // Marking element as touched in the ascent stage.
                        eAPT->sV[k].touched=-1;
                        queueBFSTREE[tree->vArr[k].level][levelCounter[tree->vArr[k].level]] = k;
                        levelCounter[tree->vArr[k].level]++;

                    }

                    eAPT->sV[k].newDelta +=
                        ((bc_t)eAPT->sV[k].newPathsToRoot / (bc_t)eAPT->sV[currElement].newPathsToRoot) *
                        (bc_t)(eAPT->sV[currElement].newDelta + 1);



                    // For the elements that are touched in the ascent stage it is necessary to
                    // to reduce the values that they previously had.
                    // In addition to this, the "parentVertex" that is connected to "vertex", i.e.
                    // the vertices of the new edge, needs to increase its betweenness centrality
                    // following the new connection, without removing the old delta value.
                    if (eAPT->sV[k].touched < 0 && ( k != parentVertex  || currElement != startVertex))
                    {
                        eAPT->sV[k].newDelta -=
                            ((bc_t)tree->vArr[k].pathsToRoot / (bc_t)tree->vArr[currElement].pathsToRoot) * 
                            (bc_t)(tree->vArr[currElement].delta + 1);

                    }
                }

            }
            STINGER_FORALL_EDGES_OF_VTX_END();

	    if (currElement != currRoot)
            {
                //forest->totalBC[currElement]+=newDelta[currElement]-tree->delta[currElement];
                eAPT->sV[currElement].totalBC += eAPT->sV[currElement].newDelta - tree->vArr[currElement].delta;
            }

        }
        deepestLevel--;

    }

    for(uint64_t k = 0; k < NV; k++)
    {
        if(eAPT->sV[k].touched != 0)
        {
            tree->vArr[k].delta = eAPT->sV[k].newDelta;
            tree->vArr[k].pathsToRoot = eAPT->sV[k].newPathsToRoot;
        }
    }

    for (uint64_t k = 0; k < NV; k++)
    {
        eAPT->sV[k].diffPath = 0;
        eAPT->sV[k].touched = 0;
        eAPT->sV[k].newDelta = 0.0;
        eAPT->sV[k].newPathsToRoot = 0;
    }
}

// Case 3
void moveDownTreeBrandes(bcForest* forest, struct stinger* sStinger, uint64_t currRoot, 
                uint64_t startVertex, uint64_t parentVertex, extraArraysPerThread *eAPT)
{
    bcTree* tree = forest->forest[currRoot];

    int64_t NV = forest->NV;

    int64_t Queue[NV], topQueue[NV];

    int64_t newLevel[NV];

    int64_t levelCounter[NV];
    uint64_t queueBFSTREE[NV][NV];

    for(uint64_t k = 0; k < NV; k++) 
    {
        eAPT->sV[k].touched = 0;
        newLevel[k]=0;
        eAPT->sV[k].newPathsToRoot=tree->vArr[k].pathsToRoot;
        eAPT->sV[k].newDelta=0.0;
        for (uint64_t l = 0; l < NV; l++)
            queueBFSTREE[k][l] = 0;
    }



    Queue[0] = startVertex;
    int64_t qStart=0,qEnd=1, tqStart=0, tqEnd=0;


    int64_t stopLevel = tree->vArr[startVertex].level;

    newLevel[startVertex] = INFINITY_MY;
    eAPT->sV[startVertex].newPathsToRoot = INFINITY_MY;
    eAPT->sV[startVertex].newDelta = 0.0;
    eAPT->sV[startVertex].touched = 1;

    int64_t deepestLevel=stopLevel;

    STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,startVertex)
    {
        int64_t k = STINGER_EDGE_DEST;

        if (tree->vArr[k].level == tree->vArr[startVertex].level)
        {
            newLevel[k] = tree->vArr[k].level;
            eAPT->sV[k].newPathsToRoot = tree->vArr[k].pathsToRoot;
            topQueue[tqEnd++] = k;
        }
        else
        {
            newLevel[k] = INFINITY_MY;
            eAPT->sV[k].newPathsToRoot = INFINITY_MY;
            if(deepestLevel < tree->vArr[k].level)
                deepestLevel = tree->vArr[k].level;
        }
        Queue[qEnd++] = k;
        eAPT->sV[k].touched = 1;
        eAPT->sV[k].newDelta = 0.0;
    }
    STINGER_FORALL_EDGES_OF_VTX_END();

    qStart = 1;
    while (qStart != qEnd)
    {
        int64_t currElement = Queue[qStart++];

        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
        {
            int64_t k = STINGER_EDGE_DEST;

            if (eAPT->sV[k].touched == 0)
            {

                if (tree->vArr[k].level == stopLevel && tree->vArr[currElement].level == (stopLevel + 1))
                {
                    newLevel[k]=tree->vArr[k].level;
                    eAPT->sV[k].newPathsToRoot=tree->vArr[k].pathsToRoot;
                    topQueue[tqEnd++]=k;
                    eAPT->sV[k].touched=1;
                    Queue[qEnd++]=k;
                    eAPT->sV[k].newDelta=0.0;
                }
                else if (tree->vArr[k].level>stopLevel)
                {
                    newLevel[k] = INFINITY_MY;
                    eAPT->sV[k].newPathsToRoot = INFINITY_MY;
                    eAPT->sV[k].touched = 1;
                    Queue[qEnd++] = k;
                    eAPT->sV[k].newDelta = 0.0;
                }

            }

        }
        STINGER_FORALL_EDGES_OF_VTX_END();

    }


    qEnd = tqEnd;
    qStart = 0;
    for(int64_t t = 0; t < qEnd; t++)
    {
        Queue[t] = topQueue[t];
        eAPT->sV[t].newDelta=0.0;

        queueBFSTREE[stopLevel][levelCounter[stopLevel]++]=topQueue[t];
    }

    eAPT->sV[parentVertex].touched = -1;
    newLevel[parentVertex] = tree->vArr[parentVertex].level;
    eAPT->sV[parentVertex].newPathsToRoot = tree->vArr[parentVertex].pathsToRoot;

    queueBFSTREE[stopLevel-1][levelCounter[stopLevel-1]++] = parentVertex;
    eAPT->sV[parentVertex].newDelta = tree->vArr[parentVertex].delta -
        ((bc_t)tree->vArr[parentVertex].pathsToRoot / (bc_t)tree->vArr[startVertex].pathsToRoot) *
        (bc_t)(tree->vArr[startVertex].delta + 1);





    // While queue is not empty
    while (qStart != qEnd)
    {
        uint64_t currElement = Queue[qStart++];

        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
        {
            uint64_t k = STINGER_EDGE_DEST;


            // If this is a neighbor and has not been found
            if (newLevel[k] > newLevel[currElement])
            {

                // Checking if "k" has been found.
                if (newLevel[k] == INFINITY_MY)
                {
                    newLevel[k] = newLevel[currElement] + 1;
                    Queue[qEnd++] = k;
                    eAPT->sV[k].newDelta = 0.0;

		    if(deepestLevel < newLevel[k])
                        deepestLevel = newLevel[k];

                    queueBFSTREE[newLevel[k]][levelCounter[newLevel[k]]++] = k;
                    if (eAPT->sV[k].touched == 0)
                        printf("AAAAAAAAAAAAAAAAAAHHHHHHHHHHHHHHHHHHHHHHH\n");
                }

                if (eAPT->sV[k].newPathsToRoot == INFINITY_MY)
                {
                    // k has not been found and therefore its paths to the roots are through its parent.
                    eAPT->sV[k].newPathsToRoot = eAPT->sV[currElement].newPathsToRoot;
                }
                else
                {
                    // k has been found and has multiple paths to the root as it has multiple parents.
                    eAPT->sV[k].newPathsToRoot += eAPT->sV[currElement].newPathsToRoot;
                }
            }
        }
        STINGER_FORALL_EDGES_OF_VTX_END();
    }



    for(int64_t k=0; k<NV;k++)
    {
        if(eAPT->sV[k].touched!=0)// && k!=parentVertex)
        {
            tree->vArr[k].level=newLevel[k];
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



    while(deepestLevel >=0 && levelCounter[deepestLevel] > 0)
    {
        int64_t currQueue = 0;

        while(currQueue<levelCounter[deepestLevel])
        {
            // Removing last element from the queue
            uint64_t currElement = queueBFSTREE[deepestLevel][currQueue];
            currQueue++;

            STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
            {
                uint64_t k = STINGER_EDGE_DEST;

                // Checking that the vertices are in different levels.
                if(tree->vArr[k].level == (tree->vArr[currElement].level-1))
                {

                    // Checking to see if "k" has been touched before.
                    if(eAPT->sV[k].touched==0)
                    {
                        eAPT->sV[k].newDelta=tree->vArr[k].delta;

                        // Marking element as touched in the ascent stage.
                        eAPT->sV[k].touched=-1;
                    	queueBFSTREE[tree->vArr[k].level][levelCounter[tree->vArr[k].level]++]=k;
                    }

                    eAPT->sV[k].newDelta +=
                        ((bc_t)eAPT->sV[k].newPathsToRoot/(bc_t)eAPT->sV[currElement].newPathsToRoot)*
                        (bc_t)(eAPT->sV[currElement].newDelta+1);


                    // For the elements that are touched in the ascent stage it is necessary to
                    // to reduce the values that they previously had.
                    // In addition to this, the "parentVertex" that is connected to "vertex", i.e.
                    // the vertices of the new edge, needs to increase its betweenness centrality
                    // following the new connection, without removing the old delta value.
                    if(eAPT->sV[k].touched<0)// && ( k!=parentVertex  || currElement!=startVertex))
                    {
                        eAPT->sV[k].newDelta -=
                        ((bc_t)tree->vArr[k].pathsToRoot/(bc_t)tree->vArr[currElement].pathsToRoot)*
                        (bc_t)(tree->vArr[currElement].delta+1);

                    }

                }

            }
            STINGER_FORALL_EDGES_OF_VTX_END();

    	    if(currElement!=currRoot)
            {
                //forest->totalBC[currElement]+=newDelta[currElement]-tree->delta[currElement];
                eAPT->sV[currElement].totalBC+=eAPT->sV[currElement].newDelta-tree->vArr[currElement].delta;
            }

        }
        deepestLevel--;

    }

    for(int64_t k=0; k<NV;k++)
    {
        if(eAPT->sV[k].touched>0)
        {
            tree->vArr[k].pathsToRoot=eAPT->sV[k].newPathsToRoot;
        }
        if(eAPT->sV[k].touched!=0)
        {
            tree->vArr[k].delta=eAPT->sV[k].newDelta;
        }
    }

    for (uint64_t k = 0; k < NV; k++)
    {
        eAPT->sV[k].diffPath = 0;
        eAPT->sV[k].touched = 0;
        eAPT->sV[k].newDelta = 0.0;
        eAPT->sV[k].newPathsToRoot = 0;
    }
}

