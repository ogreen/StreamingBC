#include "omp.h"
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "main.h"
#include "operation.h"
#include "streaming_utils.h"
#include "streaming_deletions.h"
#include "streaming_insertions.h"
#include "xmalloc.h"
#include "list.h"
#include "dsUtils.h"
#include <sys/time.h>
#include <math.h>
//#define IS_PARALLEL
#include "timer.h"


static struct timeval tvBegin, tvEnd, tvDiff;
int timeval_subtract(struct timeval *result, struct timeval *t2, struct timeval *t1)
{
	long int diff = (t2->tv_usec + 1000000 * t2->tv_sec) - (t1->tv_usec + 1000000 * t1->tv_sec);
	result->tv_sec = diff / 1000000;
	result->tv_usec = diff % 1000000;

	return (diff<0);
}

StreamingExtraInfo insertEdgeBrandes(bcForest* forest, struct stinger* sStinger,
		uint64_t newU, uint64_t newV, uint64_t * rootArrayForApproximation,int64_t NK,extraArraysPerThread** eAPT)
{

	omp_set_num_threads(NT);


	int64_t adjRootArray[NK];
	int64_t moveRootArray[NK];

	for(int64_t tk=0; tk<NK; tk++) {adjRootArray[tk]=0; moveRootArray[tk]=0;}


	uint64_t currRoot = 0;
	uint64_t samelevel = 0, compConn = 0, adjacent=0, movement=0;
        
	for(currRoot = 0; currRoot < NK; currRoot++)
	{
		uint64_t i = rootArrayForApproximation[currRoot];
		int64_t thread = 0;
		extraArraysPerThread* myExtraArrays = eAPT[thread];
		bcTree* tree = forest->forest[i];
		int64_t diff = tree->vArr[newU].level - tree->vArr[newV].level;
		// New edge is connecting vertices in the same level. Nothing needs to be done.
		//        continue;
		if(diff==0)
		{
			eAPT[thread]->samelevelCounter++;
			continue;
		}
		// Newly inserted edge is causing movement in the tee
		if(diff < -1 || diff > 1)
		{
			moveRootArray[ eAPT[thread]->movementCounter++]=i;
		}
		// Newly inserted edge is connecting vertices that were in adjacent levels before insertions
		else if(diff == -1 || diff == 1)
		{
			adjRootArray[ eAPT[thread]->adjacentCounter++]=i;
		}
	}



	for(uint64_t thread=0; thread<1; ++thread){
		//printf("Thread=%d,", thread);
		samelevel += eAPT[thread]->samelevelCounter;
		compConn += eAPT[thread]->compConnCounter;
		adjacent += eAPT[thread]->adjacentCounter;
		movement += eAPT[thread]->movementCounter;

		//		eAPT[thread]->samelevelCounter=0;
		eAPT[thread]->compConnCounter=0;
		eAPT[thread]->adjacentCounter=0;
		eAPT[thread]->movementCounter=0;
	}

	int64_t newRootArray[NK];
	int64_t counter=0, thread=0;
	for(int64_t m=0; m<movement; m++)
	{
		newRootArray[counter++]=moveRootArray[m++];
	}

	for(int64_t a=0; a<adjacent; a++)
	{
		newRootArray[counter++]=adjRootArray[a++];
	}

	double times[NT];
	uint64_t r;

	#pragma omp parallel for schedule(dynamic,1)
	for(r = 0; r < counter; r++)
		//    for(r = 0; r < NK; r++)
	{
		int64_t i=newRootArray[r];
		//            int64_t i=rootArrayForApproximation[r];
		int64_t thread = omp_get_thread_num();
		extraArraysPerThread* myExtraArrays = eAPT[thread];

		bcTree* tree = forest->forest[i];
		int64_t diff = tree->vArr[newU].level - tree->vArr[newV].level;
		// New edge is connecting vertices in the same level. Nothing needs to be done.
		//        		if(diff==0)
		//        		{
		//        			eAPT[thread]->samelevelCounter++;
		//                    continue;
		//        		}
		// Newly inserted edge is causing movement in the tee
		if(diff < -1 || diff > 1)
		{
			if(diff<-1)
			{
				moveUpTreeBrandes(forest, sStinger, i, newV, newU, (-diff) - 1,  myExtraArrays);
			}
			else
			{
				moveUpTreeBrandes(forest,  sStinger, i, newU, newV, (diff) - 1, myExtraArrays);
			}
			eAPT[thread]->movementCounter++;
		}
		// Newly inserted edge is connecting vertices that were in adjacent levels before insertions
		else if(diff == -1 || diff == 1)
		{
			if(diff==-1)
			{
				addEdgeWithoutMovementBrandes(forest, sStinger, i, newV, newU, tree->vArr[newU].pathsToRoot,myExtraArrays);
			}
			else
			{
				addEdgeWithoutMovementBrandes(forest, sStinger, i, newU, newV, tree->vArr[newV].pathsToRoot, myExtraArrays);
			}
			eAPT[thread]->adjacentCounter++;
		}

		#pragma omp barriar

		int64_t tlow=(NV*thread)/NT;
		int64_t thigh=(NV*(thread+1))/NT-1;
                printf("tlow: %d\n", tlow);
                printf("thigh: %d\n", thigh);
                printf("NT: %d\n", NT);
                printf("before_dynamic: totalBC[newU: %d]: %d\n", newU, forest->totalBC[newU]); 
                printf("before_dynamic: totalBC[newV: %d]: %d\n", newV, forest->totalBC[newV]); 
		for(uint64_t v=tlow;v<thigh;v++){
                        //printf("before_dynamic: totalBC[%d]: %d\n", v, forest->totalBC[v]);
			for(uint64_t t=0;t<NT;t++){
				forest->totalBC[v]+=eAPT[t]->sV[v].totalBC;
			}
                        //printf("after_dynamic: totalBC[%d]: %d\n", v, forest->totalBC[v]);
		}
                //printf("eAPT[0]->sV[newU: %d].totalBC: %d\n", newU, eAPT[0]->sV[newU].totalBC);
                printf("after_dynamic: totalBC[newU: %d]: %d\n", newU, forest->totalBC[newU]); 
                printf("after_dynamic: totalBC[newV: %d]: %d\n", newV, forest->totalBC[newV]);

	}



        StreamingExtraInfo returnSEI={0,0,0,0};
        returnSEI.sameLevel= samelevel;
        returnSEI.adjacent = adjacent;
        returnSEI.movement = movement;

        return returnSEI;
}

StreamingExtraInfo deleteEdgeBrandes(bcForest *forest, struct stinger *sStinger, uint64_t oldU, uint64_t oldV,
                            uint64_t *rootArrayForApproximation, int64_t NK, extraArraysPerThread **eAPT)
{
    omp_set_num_threads(NT);

    int64_t adjRootArray[NK];
    int64_t moveRootArray[NK];

    for (int64_t i = 0; i < NK; i++)
    {
        adjRootArray[i] = 0;
        moveRootArray[i] = 0;
    }

    uint64_t currRoot = 0;
    uint64_t samelevel = 0, compConn = 0, adjacent = 0, movement = 0;

    for (currRoot = 0; currRoot < NK; currRoot++)
    {
        uint64_t i = rootArrayForApproximation[currRoot];
        int64_t thread = 0;
        extraArraysPerThread *myExtraArrays = eAPT[thread];
        bcTree *tree = forest->forest[i];
        int64_t diff = tree->vArr[oldU].level - tree->vArr[oldV].level;

        // Case I: Old edge connected vertices in the same level. Nothing needs to be done.
        if (diff == 0)
        {
            eAPT[thread]->samelevelCounter++;
            continue;
        }

        int64_t numParents = 0;
        int64_t childVertex = oldU;
        if (tree->vArr[oldU].level < tree->vArr[oldV].level)
        {
            childVertex = oldV;
        }

        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, childVertex)
        {
            uint64_t neighbor = STINGER_EDGE_DEST;
            if (tree->vArr[neighbor].level + 1 == tree->vArr[childVertex].level)
            {
                numParents++;
            }
        }
        STINGER_FORALL_EDGES_OF_VTX_END();
        
        if (numParents > 1)
        {
            adjRootArray[eAPT[thread]->adjacentCounter++] = i;
        }
        else if (numParents == 1)
        {
            moveRootArray[eAPT[thread]->movementCounter++] = i;
        }
    }

    for (uint64_t thread = 0; thread < 1; thread++)
    {
        samelevel += eAPT[thread]->samelevelCounter;
        compConn += eAPT[thread]->compConnCounter;
        adjacent += eAPT[thread]->adjacentCounter;
        movement += eAPT[thread]->movementCounter;

        eAPT[thread]->compConnCounter = 0;
        eAPT[thread]->adjacentCounter = 0;
        eAPT[thread]->movementCounter = 0;
    }

    int64_t newRootArray[NK];
    int64_t counter = 0, thread = 0;
    for (int64_t m = 0; m < movement; m++)
    {
        newRootArray[counter++] = moveRootArray[m++]; // Why are there two m++'s?
    }

    for (int64_t a = 0; a < adjacent; a++)
    {
        newRootArray[counter++] = adjRootArray[a++];
    }

    double times[NT];
    uint64_t r;

    #pragma omp parallel for schedule(dynamic,1)
    for (r = 0; r < counter; r++)
    {
        int64_t i = newRootArray[r];
        int64_t thread = omp_get_thread_num();
        extraArraysPerThread *myExtraArrays = eAPT[thread];

        bcTree *tree = forest->forest[i];

        int64_t numParents;
        int64_t childVertex = oldU;
        int64_t parentVertex = oldV;
        if (tree->vArr[oldU].level < tree->vArr[oldV].level)
        {
            childVertex = oldV;
            parentVertex = oldU;
        }
        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, childVertex)
        {
            uint64_t neighbor = STINGER_EDGE_DEST;

            if (tree->vArr[neighbor].level + 1 == tree->vArr[childVertex].level)
            {
                numParents++;
            }
        }
        STINGER_FORALL_EDGES_OF_VTX_END();

        if (numParents > 1)
        {
            deleteEdgeWithoutMovement(forest, sStinger, i, childVertex, parentVertex, tree->vArr[parentVertex].pathsToRoot, myExtraArrays);
            eAPT[thread]->adjacentCounter++;
        }
        else
        {
            moveDownTreeBrandes(forest, sStinger, i, childVertex, parentVertex, 1, myExtraArrays);
            eAPT[thread]->movementCounter++;
        }

        int64_t tlow = (NV * thread) / NT;
        int64_t thigh = (NV * (thread + 1)) / NT - 1;
        #pragma omp barriar
        for (uint64_t v = tlow; v < thigh; v++)
        {
            for (uint64_t t = 0; t < NT; t++)
            {
                forest->totalBC[v] += eAPT[t]->sV[v].totalBC;
            }
            printf("dynamic: totalBC[%d]: %d\n", v, forest->totalBC[v]);
        }
    }

    StreamingExtraInfo returnSEI = {0,0,0,0};

    returnSEI.sameLevel = samelevel;
    returnSEI.adjacent = adjacent;
    returnSEI.movement = movement;
    
    return returnSEI;
}
