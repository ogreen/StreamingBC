#include "omp.h" 
#include <stdint.h> 
#include <stdlib.h> 
#include <stdio.h>

#include "streamingbc_aux.h"
//#include "timer.h"

//#define ANCHORED 3
#define PARENT_ANCHORED 3 
#define SIBLING_ANCHORED 4

void addEdgeWithoutMovementBrandes(bcForest* forest, struct stinger* sStinger,
		uint64_t currRoot, uint64_t startVertex, uint64_t parentVertex,
		uint64_t addedPathsToRoot, extraArraysPerThread* eAPT, uint64_t cores) {

        //omp_set_num_threads(cores);
	bcTree* tree = forest->forest[currRoot];

	uint64_t *QueueDown = eAPT->QueueDown;
	uint64_t *QueueUp = eAPT->QueueUp;
        uint64_t *QueueDownBorders = eAPT->QueueSame;

	int64_t NV = forest->NV;

        
	eAPT->sV[startVertex].newPathsToRoot = tree->vArr[startVertex].pathsToRoot;
        eAPT->sV[parentVertex].newEdgesBelow = tree->vArr[parentVertex].edgesBelow;
        eAPT->sV[startVertex].newEdgesBelow = tree->vArr[startVertex].edgesBelow;
        eAPT->sV[parentVertex].newEdgesAbove = tree->vArr[parentVertex].edgesAbove;
        eAPT->sV[startVertex].newEdgesAbove = tree->vArr[startVertex].edgesAbove;
	eAPT->sV[startVertex].touched = 1;
        eAPT->sV[startVertex].newPathsToRoot += addedPathsToRoot;
	eAPT->sV[startVertex].diffPath = addedPathsToRoot;
        eAPT->sV[startVertex].newEdgesAbove += eAPT->sV[parentVertex].newEdgesAbove + 1;
        eAPT->sV[parentVertex].newEdgesBelow += eAPT->sV[startVertex].newEdgesBelow + 1;

	QueueDown[0] = startVertex;
	//int64_t qStart=0,qEnd=1;

        int64_t *qStart = &(eAPT->qStart);
        int64_t *qEnd = &(eAPT->qEnd);
        int64_t *qStart_nxt = &(eAPT->qStart_nxt);
        int64_t *qEnd_nxt = &(eAPT->qEnd_nxt);

        int64_t qDownBIndex = 0;
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
            QueueDownBorders[qDownBIndex++] = *qStart;
            QueueDownBorders[qDownBIndex++] = *qEnd;

            //#pragma omp parallel for schedule(dynamic, 1)
            int64_t thread_nums = cores;
            if ((*qEnd - *qStart) < thread_nums) {
                thread_nums = *qEnd - *qStart;
            }

            #pragma omp parallel num_threads(thread_nums)
            {
                #pragma omp for
                for (int64_t i = *qStart; i < *qEnd; i++) {
                    uint64_t currElement = QueueDown[i];
                    int64_t levelCurrPlusOne = tree->vArr[currElement].level + 1;
                    int64_t touchedCurrPlusOne = eAPT->sV[currElement].touched + 1;

                    STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
                    {
                        uint64_t k = STINGER_EDGE_DEST;
	
                        // if this vertex has not been added yet
                        if(levelCurrPlusOne == (tree->vArr[k].level)){
                                
                            if (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 0, currElement)) {
                                // Checking if a "deeper level" has been reached.
                                if (deepestLevel < tree->vArr[k].level) {
                                        deepestLevel = tree->vArr[k].level;
                                }

                                __sync_fetch_and_add(&(eAPT->sV[k].newEdgesAbove), tree->vArr[k].edgesAbove);
                                __sync_fetch_and_sub(&(eAPT->sV[k].newEdgesAbove), tree->vArr[currElement].edgesAbove);
                                __sync_fetch_and_add(&(eAPT->sV[k].newEdgesAbove), eAPT->sV[currElement].newEdgesAbove);

                                //NEW
                                __sync_fetch_and_add(&(eAPT->sV[k].newPathsToRoot), tree->vArr[k].pathsToRoot);
                                
                                // insert this vertex into the BFS queue
                                QueueDown[__sync_fetch_and_add(qEnd_nxt, 1)] = k;
                                // indicate that it is in the next level of the BFS
                                // add new paths to root that go through current BFS Vertex
                                __sync_fetch_and_add(&(eAPT->sV[k].newPathsToRoot), eAPT->sV[currElement].diffPath);
                                // pass on my new paths to root for its search
                                __sync_fetch_and_add(&(eAPT->sV[k].diffPath), eAPT->sV[currElement].diffPath);
                            }
                            // otherwise if it has been touched, but is specifically in the next level
                            // of the search (meaning it has more than one edge to the current level)
                            else if (eAPT->sV[k].touched != currElement) {
                                __sync_fetch_and_add(&(eAPT->sV[k].newEdgesAbove), -tree->vArr[currElement].edgesAbove);
                                __sync_fetch_and_add(&(eAPT->sV[k].newEdgesAbove), eAPT->sV[currElement].newEdgesAbove);

                                // add new paths to root that go through current BFS Vertex
                                __sync_fetch_and_add(&(eAPT->sV[k].newPathsToRoot), eAPT->sV[currElement].diffPath);
                                // pass on my new paths to root for its search
                                __sync_fetch_and_add(&(eAPT->sV[k].diffPath), eAPT->sV[currElement].diffPath);
                            }
                        }
                    }
                    STINGER_FORALL_EDGES_OF_VTX_END();

#if COUNT_TRAVERSALS==1
                    eAPT->dynamicTraverseEdgeCounter+=stinger_typed_outdegree(sStinger,currElement,0);
                    eAPT->dynamicTraverseVerticeCounter++;
#endif
                }
            }
            *qStart = *qStart_nxt;
            *qEnd = *qEnd_nxt;
            *qStart_nxt = *qEnd;
            *qEnd_nxt = *qStart_nxt;
	}

        #if 0
        printf("currRoot, %ld\n", currRoot);
        for (int64_t i = 0; i < qDownBIndex; i += 2) {
            printf("border: %ld, %ld\n", QueueDownBorders[i], QueueDownBorders[i + 1]);
        }
        #endif

	int64_t QUpStart=0,QUpEnd=0;
	// int64_t myCont=1;
	int64_t currElement;
	(*qEnd)--;

	int64_t qDownEndMarker= *qEnd;
        int64_t qDownEnd = *qEnd;

        *qStart = 0;
        *qEnd = 0;
        *qStart_nxt = 0;
        *qEnd_nxt = 0;

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
	while(!(qDownBIndex <= 0 && *qStart >= *qEnd && *qStart_nxt >= *qEnd_nxt)) {
            
            if (qDownBIndex >= 2) {
            
                int64_t thread_nums = cores;
                if (QueueDownBorders[qDownBIndex - 1] - QueueDownBorders[qDownBIndex - 2] < thread_nums) {
                    thread_nums = QueueDownBorders[qDownBIndex - 1] - QueueDownBorders[qDownBIndex - 2];
                }
                
                #pragma omp parallel num_threads(thread_nums)
                {
                    #pragma omp for
                    for (int64_t i = QueueDownBorders[qDownBIndex - 2]; i < QueueDownBorders[qDownBIndex - 1]; i++) {

                        int64_t currElement = QueueDown[i];
                        //eAPT->sV[currElement].touched = -1;
                        int64_t levelCurrMinusOne = tree->vArr[currElement].level-1;
                        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
                        {
                            uint64_t k = STINGER_EDGE_DEST;

                            if(tree->vArr[k].level == levelCurrMinusOne){
                                // Checking to see if "k" has been touched before.
                                //if(eAPT->sV[k].touched==0){
                                if (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 0, -1)) {
                                    eAPT->sV[k].newDelta += tree->vArr[k].delta;
                                    // Marking element as touched in the ascent stage.
                                    eAPT->sV[k].touched=-1;
                                    //QueueUp[QUpEnd] = k;
                                    //QUpEnd++;
                                    //QueueUp[*qEnd_nxt] = k;
                                    //(*qEnd_nxt)++;
                                    QueueUp[__sync_fetch_and_add(qEnd_nxt, 1)] = k;
                                    //eAPT->sV[k].newPathsToRoot += tree->vArr[k].pathsToRoot;
                                    __sync_fetch_and_add(&(eAPT->sV[k].newPathsToRoot), tree->vArr[k].pathsToRoot);
                                    if (k != parentVertex) {
                                        //eAPT->sV[k].newEdgesBelow += tree->vArr[k].edgesBelow;
                                        __sync_fetch_and_add(&(eAPT->sV[k].newEdgesBelow), tree->vArr[k].edgesBelow);
                                    }
                                }
                                if (k != parentVertex && tree->vArr[k].level <= tree->vArr[parentVertex].level) {
                                    //eAPT->sV[k].newEdgesBelow -= tree->vArr[currElement].edgesBelow;
                                    __sync_fetch_and_add(&(eAPT->sV[k].newEdgesBelow), -tree->vArr[currElement].edgesBelow);
                                    //eAPT->sV[k].newEdgesBelow += eAPT->sV[currElement].newEdgesBelow;
                                    __sync_fetch_and_add(&(eAPT->sV[k].newEdgesBelow), eAPT->sV[currElement].newEdgesBelow);
                                }
                            }
                            

                            if (tree->vArr[k].level == tree->vArr[currElement].level + 1 && eAPT->sV[k].touched != 0) {
                                eAPT->sV[currElement].newDelta +=
                                            ((bc_t)eAPT->sV[currElement].newPathsToRoot/(bc_t)eAPT->sV[k].newPathsToRoot)*
                                            (bc_t)(eAPT->sV[k].newDelta+1);


                                // For the elements that are touched in the ascent stage it is necessary to
                                // to reduce the values that they previously had.
                                // In addition to this, the "parentVertex" that is connected to "vertex", i.e.
                                // the vertices of the new edge, needs to increase its betweenness centrality
                                // following the new connection, without removing the old delta value.
                                if(eAPT->sV[currElement].touched<0 && (currElement != parentVertex || k != startVertex)){
                                        eAPT->sV[currElement].newDelta -=
                                                ((bc_t)tree->vArr[currElement].pathsToRoot/(bc_t)tree->vArr[k].pathsToRoot)*
                                                (bc_t)(tree->vArr[k].delta+1);
                                }
                            }
                        }
                        STINGER_FORALL_EDGES_OF_VTX_END();

#if COUNT_TRAVERSALS==1
                        eAPT->dynamicTraverseEdgeCounter+=stinger_typed_outdegree(sStinger,currElement,0);
                        eAPT->dynamicTraverseVerticeCounter++;
#endif
                        if(currElement!=currRoot){
                            eAPT->sV[currElement].totalBC+=eAPT->sV[currElement].newDelta - tree->vArr[currElement].delta;
                        }
                    }
                }
            }

            qDownBIndex -= 2;

            *qStart = *qStart_nxt;
            *qEnd = *qEnd_nxt;
            *qStart_nxt = *qEnd;
            *qEnd_nxt = *qStart_nxt;
            int64_t thread_nums = cores;
            if (*qEnd - *qStart < cores) {
                thread_nums = *qEnd - *qStart;
            }

            #pragma omp parallel num_threads(thread_nums)
            {
                #pragma omp for
                for (int64_t i = *qStart; i < *qEnd; i++) {
                    int64_t currElement = QueueUp[i];
                    //eAPT->sV[currElement].touched = -1;
                    //(*qStart)++;
                    //__sync_fetch_and_add(qStart, 1);
                    int64_t levelCurrMinusOne = tree->vArr[currElement].level-1;
                    STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
                    {
                        uint64_t k = STINGER_EDGE_DEST;

                        if(tree->vArr[k].level == levelCurrMinusOne){
                                // Checking to see if "k" has been touched before.
                            //if(eAPT->sV[k].touched==0)
                            if (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 0, -1)) {
                                eAPT->sV[k].newDelta += tree->vArr[k].delta;
                                // Marking element as touched in the ascent stage.
                                //eAPT->sV[k].touched=-1;
                                //QueueUp[QUpEnd] = k;
                                //QUpEnd++;
                                //QueueUp[*qEnd_nxt] = k;
                                //(*qEnd_nxt)++;
                                QueueUp[__sync_fetch_and_add(qEnd_nxt, 1)] = k;
                                //eAPT->sV[k].newPathsToRoot += tree->vArr[k].pathsToRoot;
                                __sync_fetch_and_add(&(eAPT->sV[k].newPathsToRoot), tree->vArr[k].pathsToRoot);
                                if (k != parentVertex) {
                                    //eAPT->sV[k].newEdgesBelow += tree->vArr[k].edgesBelow;
                                    __sync_fetch_and_add(&(eAPT->sV[k].newEdgesBelow), tree->vArr[k].edgesBelow);
                                }
                            }

                            if (k != parentVertex && tree->vArr[k].level <= tree->vArr[parentVertex].level) {
                                //eAPT->sV[k].newEdgesBelow -= tree->vArr[currElement].edgesBelow;
                                __sync_fetch_and_add(&(eAPT->sV[k].newEdgesBelow), -tree->vArr[currElement].edgesBelow);
                                //eAPT->sV[k].newEdgesBelow += eAPT->sV[currElement].newEdgesBelow;
                                __sync_fetch_and_add(&(eAPT->sV[k].newEdgesBelow), eAPT->sV[currElement].newEdgesBelow);
                            }

                        }

                        if (tree->vArr[k].level == tree->vArr[currElement].level + 1 && eAPT->sV[k].touched != 0) {
                            eAPT->sV[currElement].newDelta +=
                                        ((bc_t)eAPT->sV[currElement].newPathsToRoot/(bc_t)eAPT->sV[k].newPathsToRoot)*
                                        (bc_t)(eAPT->sV[k].newDelta+1);

                            // For the elements that are touched in the ascent stage it is necessary to
                            // to reduce the values that they previously had.
                            // In addition to this, the "parentVertex" that is connected to "vertex", i.e.
                            // the vertices of the new edge, needs to increase its betweenness centrality
                            // following the new connection, without removing the old delta value.
                            if(eAPT->sV[currElement].touched<0 && (currElement != parentVertex || k != startVertex)){
                                    eAPT->sV[currElement].newDelta -=
                                            ((bc_t)tree->vArr[currElement].pathsToRoot/(bc_t)tree->vArr[k].pathsToRoot)*
                                            (bc_t)(tree->vArr[k].delta+1);
                            }
                        }
                    }
                    STINGER_FORALL_EDGES_OF_VTX_END();

#if COUNT_TRAVERSALS==1
                    eAPT->dynamicTraverseEdgeCounter+=stinger_typed_outdegree(sStinger,currElement,0);
                    eAPT->dynamicTraverseVerticeCounter++;
#endif
                    if(currElement!=currRoot){
                            eAPT->sV[currElement].totalBC+=eAPT->sV[currElement].newDelta - tree->vArr[currElement].delta;
                    }
                }
            }
            *qStart = *qEnd;
	}

	for(uint64_t c = 0; c <= qDownEndMarker; c++){
		uint64_t k=QueueDown[c];
		tree->vArr[k].delta=eAPT->sV[k].newDelta;
		tree->vArr[k].pathsToRoot=eAPT->sV[k].newPathsToRoot;
                tree->vArr[k].edgesAbove = eAPT->sV[k].newEdgesAbove;
                eAPT->sV[k].diffPath=0;
		eAPT->sV[k].touched=0;
		eAPT->sV[k].newDelta=0.0;
		eAPT->sV[k].newPathsToRoot=0;
                eAPT->sV[k].newEdgesAbove = 0;
	}
        eAPT->sV[startVertex].newEdgesAbove = 0;
        eAPT->sV[parentVertex].newEdgesAbove = 0;

	//for(uint64_t c = 0; c < QUpEnd; c++){
        for (uint64_t c = 0; c < *qEnd; c++) {
		uint64_t k=QueueUp[c];
		tree->vArr[k].delta=eAPT->sV[k].newDelta;
		tree->vArr[k].edgesBelow = eAPT->sV[k].newEdgesBelow;
		eAPT->sV[k].diffPath=0;
		eAPT->sV[k].touched=0;
		eAPT->sV[k].newDelta=0.0;
                eAPT->sV[k].newEdgesBelow = 0;
		eAPT->sV[k].newPathsToRoot=0;
                eAPT->sV[k].newEdgesBelow = 0;
	}

        eAPT->sV[parentVertex].newEdgesBelow = 0;
        eAPT->sV[startVertex].newEdgesBelow = 0;

        eAPT->qStart = 0;
        eAPT->qEnd = 0;
        eAPT->qStart_nxt = 0;
        eAPT->qEnd_nxt = 0;

	return;
}


void moveUpTreeBrandes(bcForest* forest, struct stinger* sStinger,
		uint64_t currRoot, uint64_t startVertex, uint64_t parentVertex,
		uint64_t prevDist, extraArraysPerThread* eAPT, uint64_t cores){

	bcTree* tree = forest->forest[currRoot];
	int64_t NV = forest->NV;
	uint64_t *QueueDown=eAPT->QueueDown;
	uint64_t *QueueUp=eAPT->QueueUp;
	uint64_t *QueueSame=eAPT->QueueSame;
        uint64_t *QueueDownBorders = eAPT->Stack;

	list_ptr* multiLevelQueues = eAPT->multiLevelQueues;
        queue_t* queue = eAPT->queue;
        level_node_t* levelIndices = eAPT->levelIndices;

	//NEW
	eAPT->sV[parentVertex].newPathsToRoot = tree->vArr[parentVertex].pathsToRoot;
	eAPT->sV[startVertex].newPathsToRoot = tree->vArr[startVertex].pathsToRoot;

	//int64_t qStart = 0;
	//int64_t qEnd = 1;

        int64_t *qStart = &(eAPT->qStart);
        int64_t *qEnd = &(eAPT->qEnd);
        int64_t *qStart_nxt = &(eAPT->qStart_nxt);
        int64_t *qEnd_nxt = &(eAPT->qEnd_nxt);
        int64_t qDownBIndex = 0;

        int64_t *qStartSame = &(eAPT->qStartSame);
        int64_t *qEndSame = &(eAPT->qEndSame);
        int64_t *qStartSame_nxt = &(eAPT->qStartSame_nxt);
        int64_t *qEndSame_nxt = &(eAPT->qEndSame_nxt);

        *qEnd = 1;
        *qStart_nxt = 1;
        *qEnd_nxt = 1;


	QueueDown[0] = startVertex;
	eAPT->sV[startVertex].touched = 1;
	eAPT->sV[startVertex].newPathsToRoot = eAPT->sV[parentVertex].newPathsToRoot;
	eAPT->sV[startVertex].diffPath = eAPT->sV[parentVertex].newPathsToRoot;

	eAPT->sV[startVertex].movementDelta = prevDist;
	eAPT->sV[startVertex].IMoved = 1;

        eAPT->sV[parentVertex].newEdgesAbove = tree->vArr[parentVertex].edgesAbove;
        eAPT->sV[startVertex].newEdgesAbove = eAPT->sV[parentVertex].newEdgesAbove + 1;
	int64_t deepestLevel = tree->vArr[parentVertex].level+1;

	// Starting BFS decent from "startVertex", down to all the vertices that have shortest paths through "startVertex"
	// All elements that will be touched will receive a positive value in their touched field.
	// In this implementation, "STACKS" are not used for the "moving up" stage. Rather, a multi-level queue is used.
	// Each level in the tree(max depth NV) has a queue and a counter specifiying how deep a specific deepth-queue is.
	// For simplicity, all elements are pushed both into the multi-level queue and into the regular queue which is used
	// for the BFS traversal.
	while(*qStart < *qEnd) {
                QueueDownBorders[qDownBIndex++] = *qStart;
                QueueDownBorders[qDownBIndex++] = *qEnd;

                int64_t thread_nums = cores;
                if ((*qEnd - *qStart) < cores) {
                    thread_nums = *qEnd - *qStart;
                }

                #pragma omp parallel num_threads(thread_nums)
                {
                    #pragma omp for
                    for (int64_t i = *qStart; i < *qEnd; i++) {
                        //int64_t currElement = QueueDown[*qStart];

                        int64_t currElement = QueueDown[i];
                        int64_t touchedCurrPlusOne = eAPT->sV[currElement].touched+1;
                        __sync_fetch_and_add(&(eAPT->sV[currElement].newEdgesAbove), -eAPT->sV[currElement].newEdgesAbove);

                        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
                        {
                                uint64_t k = STINGER_EDGE_DEST;

                                int64_t computedDelta = eAPT->sV[currElement].movementDelta -
                                                    (tree->vArr[currElement].level - tree->vArr[k].level + 1);

                                int64_t newCurrLevel = 0;
                                __sync_fetch_and_add(&newCurrLevel, tree->vArr[currElement].level);
                                __sync_fetch_and_add(&newCurrLevel, -eAPT->sV[currElement].movementDelta);

                                uint64_t newKLevel = 0;
                                __sync_fetch_and_add(&newKLevel, tree->vArr[k].level);
                                __sync_fetch_and_add(&newKLevel, -computedDelta);

                                int go_in = 0;
                                if (computedDelta < 0 && eAPT->sV[k].touched == 0) {
                                    if (computedDelta >= 0 && newKLevel < newCurrLevel) {
                                        __sync_fetch_and_add(&(eAPT->sV[k].newEdgesAbove), tree->vArr[k].edgesAbove);
                                        __sync_fetch_and_add(&(eAPT->sV[currElement].newEdgesAbove), eAPT->sV[k].newEdgesAbove + 1);
                                    }

                                    if (computedDelta < 0 && tree->vArr[k].level < newCurrLevel) {
                                        __sync_fetch_and_add(&(eAPT->sV[currElement].newEdgesAbove), tree->vArr[k].edgesAbove + 1);
                                    }
                                }
                                //if(eAPT->sV[k].touched == 0 ){}
                                else if (go_in || __sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 0, touchedCurrPlusOne)) {
                                        if (computedDelta >= 0 && newKLevel < newCurrLevel) {
                                            __sync_fetch_and_add(&(eAPT->sV[k].newEdgesAbove), tree->vArr[k].edgesAbove);
                                            __sync_fetch_and_add(&(eAPT->sV[currElement].newEdgesAbove), eAPT->sV[k].newEdgesAbove + 1);
                                        }

                                        if (computedDelta < 0 && tree->vArr[k].level < newCurrLevel) {
                                            __sync_fetch_and_add(&(eAPT->sV[currElement].newEdgesAbove), tree->vArr[k].edgesAbove + 1);
                                        }
                                        
                                        // if the adjacent vertex should be moved, put it in the queue
                                        if(computedDelta > 0){
                                            __sync_fetch_and_add(&(eAPT->sV[k].newPathsToRoot), eAPT->sV[currElement].diffPath);
                                            __sync_fetch_and_add(&(eAPT->sV[k].diffPath), eAPT->sV[currElement].diffPath);
                                            __sync_fetch_and_add(&(eAPT->sV[k].movementDelta), computedDelta);
                                            __sync_fetch_and_add(&(eAPT->sV[k].IMoved), 2);
                                            QueueDown[__sync_fetch_and_add(qEnd_nxt, 1)] = k;
                                        }
                                        // Vertex that will not be moved has been found.
                                        else if(computedDelta == 0){
                                            //NEW
                                            printf("here, currRoot, currElement, k: %ld, %ld, %ld\n", currRoot, currElement, k);
                                            __sync_fetch_and_add(&(eAPT->sV[k].newPathsToRoot), tree->vArr[k].pathsToRoot);
                                            __sync_fetch_and_add(&(eAPT->sV[k].newPathsToRoot), eAPT->sV[currElement].diffPath);
                                            __sync_fetch_and_add(&(eAPT->sV[k].diffPath), eAPT->sV[currElement].diffPath);
                                            __sync_fetch_and_add(&(eAPT->sV[k].movementDelta), computedDelta);
                                            __sync_fetch_and_add(&(eAPT->sV[k].IMoved), -eAPT->sV[k].IMoved);
                                            QueueDown[__sync_fetch_and_add(qEnd_nxt, 1)] = k;
                                        }
                                        // Vertex that the number of shortest path to the root does not change has been found.
                                        // This vertex is not marked as it might be touched on the way up.
                                                                                 
                                        // if adjacent and in the next level
                                }
                                //else if(eAPT->sV[k].touched == eAPT->sV[currElement].touched + 1)
                                else if(eAPT->sV[k].touched == touchedCurrPlusOne){
                                        //int64_t computedDelta = eAPT->sV[currElement].movementDelta - (tree->vArr[currElement].level - tree->vArr[k].level + 1);
                                        if (computedDelta >= 0 && newKLevel < newCurrLevel) {
                                            __sync_fetch_and_add(&(eAPT->sV[currElement].newEdgesAbove), eAPT->sV[k].newEdgesAbove + 1);
                                        }

                                        if (computedDelta < 0 && tree->vArr[k].level < newCurrLevel) {
                                            __sync_fetch_and_add(&(eAPT->sV[currElement].newEdgesAbove), eAPT->sV[k].newEdgesAbove + 1);
                                        } 
                                        if(computedDelta >= 0){
                                            __sync_fetch_and_add(&(eAPT->sV[k].newPathsToRoot), eAPT->sV[currElement].diffPath);
                                            __sync_fetch_and_add(&(eAPT->sV[k].diffPath), eAPT->sV[currElement].diffPath);
                                        }

                                        
                                }
                                else if (computedDelta >= 0 && newKLevel < newCurrLevel) {
                                    __sync_fetch_and_add(&(eAPT->sV[currElement].newEdgesAbove), eAPT->sV[k].newEdgesAbove + 1);
                                }
                                else if (computedDelta < 0 && tree->vArr[k].level < newCurrLevel) {
                                    __sync_fetch_and_add(&(eAPT->sV[currElement].newEdgesAbove), eAPT->sV[k].newEdgesAbove + 1);
                                }
                        }
                        STINGER_FORALL_EDGES_OF_VTX_END();
                        

#if COUNT_TRAVERSALS==1
                        eAPT->dynamicTraverseEdgeCounter+=stinger_typed_outdegree(sStinger,currElement,0);
                        eAPT->dynamicTraverseVerticeCounter++;
#endif
                        // move ourself and retire
                        __sync_fetch_and_add(&(tree->vArr[currElement].level), -eAPT->sV[currElement].movementDelta);
                        appendDS2(queue, levelIndices, tree->vArr[currElement].level, currElement, omp_get_thread_num());
                        // Checking if a "deeper level" has been reached.
                        if(deepestLevel<tree->vArr[currElement].level)
                                deepestLevel=tree->vArr[currElement].level;
                    }
                }

                *qStart = *qStart_nxt;
                *qEnd = *qEnd_nxt;
                *qStart_nxt = *qEnd;
                *qEnd_nxt = *qStart_nxt;
	}

    	// Starting Multi-Level "BFS" ascent.
	*qEnd=0;
        queue_node_t* temp_node;
	for(int lev = tree->vArr[startVertex].level; lev<NV; lev++){
                temp_node = getFirstDS(queue, levelIndices, lev);
                while (temp_node != NULL) {
                        QueueDown[(*qEnd)++] = temp_node->data;
                        deleteFirstDS(queue, levelIndices, lev);
                        temp_node = getFirstDS(queue, levelIndices, lev);
		}
	}

        (*qEnd)--;
	//NEW
	int64_t qDownEndMarker= *qEnd;

	int64_t QUpStart=0,QUpEnd=0;
	int64_t QSameStart=0,QSameEnd=0;
	int64_t currElement=0; //dummy initilization - variable will be initialized in function.
	int64_t upCounter=0;
	int64_t depthDown=tree->vArr[QueueDown[*qEnd]].level,depthUp=-1,depthSame=-1;
        
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

        while (!(qDownBIndex <= 0 && *qStart >= *qEnd && *qStart_nxt >= *qEnd_nxt && *qStartSame >= *qEndSame && *qStartSame_nxt >= *qEndSame_nxt)) {

            int operation = -1; // 0 - down, 1 - up, 2 - same
            if(depthDown==-1 && depthSame==-1 && depthUp==-1)
                break;
            if(depthUp>=depthSame && depthUp >=depthDown){
                operation = 1;
                if (*qEnd_nxt > *qStart_nxt)
                    depthUp=-1;
                else {
                    depthUp=tree->vArr[QueueUp[*qStart_nxt]].level;
                }
            }
            else if(depthDown>=depthSame && depthDown >=depthUp){
                operation = 0;
                if (qDownBIndex < 2 || QueueDownBorders[qDownBIndex - 2] > QueueDownBorders[qDownBIndex - 1])
                    depthDown = -1;
                else if (qDownBIndex > 2) {
                    depthDown=tree->vArr[QueueDown[QueueDownBorders[qDownBIndex - 2] - 1]].level;
                }
            }
            
            else if(depthDown<=depthSame && depthUp <=depthSame){
                operation = 2;
                if (*qEndSame_nxt > *qStartSame_nxt)
                    depthSame=-1;
                else
                    depthSame=tree->vArr[QueueSame[*qStartSame_nxt]].level;
            }

            if (operation == 0 && qDownBIndex >= 2) {

                int64_t thread_nums = cores;
                if (QueueDownBorders[qDownBIndex - 1] - QueueDownBorders[qDownBIndex - 2] < thread_nums) {
                    thread_nums = QueueDownBorders[qDownBIndex - 1] - QueueDownBorders[qDownBIndex - 2];
                }

                #pragma omp parallel num_threads(thread_nums)
                {
                    #pragma omp for
                    for (int64_t i = QueueDownBorders[qDownBIndex - 2]; i < QueueDownBorders[qDownBIndex - 1]; i++) {

                        int64_t currElement = QueueDown[i];
                        int64_t levelCurrMinusOne = tree->vArr[currElement].level-1;
                        eAPT->sV[currElement].newEdgesBelow = 0;
                        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
                        {
                                uint64_t k = STINGER_EDGE_DEST;
                                // Checking that the vertices are in different levels.

                                if (tree->vArr[k].level == tree->vArr[currElement].level + 1) {

                                    if (eAPT->sV[k].touched == 0) {
                                        eAPT->sV[currElement].newEdgesBelow += tree->vArr[k].edgesBelow + 1;
                                    } else {
                                        eAPT->sV[currElement].newEdgesBelow += eAPT->sV[k].newEdgesBelow + 1;
                                    }
                                }

                                if(tree->vArr[k].level == levelCurrMinusOne){
                                        // Checking to see if "k" has been touched before.
                                        if (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 0, -1)) {
                                                eAPT->sV[k].newDelta += tree->vArr[k].delta;
                                                upCounter++;
                                                // Marking element as touched in the ascent stage.
                                                eAPT->sV[k].touched=-1;

                                                __sync_bool_compare_and_swap(&depthUp, -1, tree->vArr[k].level);
                                                QueueUp[__sync_fetch_and_add(qEnd_nxt, 1)] = k;

                                                if (k != parentVertex)
                                                    eAPT->sV[k].newPathsToRoot += tree->vArr[k].pathsToRoot;
                                        }
                                }
                                
                                if (tree->vArr[k].level == tree->vArr[currElement].level + 1 && eAPT->sV[k].touched != 0) {
                                    eAPT->sV[currElement].newDelta +=
                                                ((bc_t)eAPT->sV[currElement].newPathsToRoot/(bc_t)eAPT->sV[k].newPathsToRoot)*
                                                (bc_t)(eAPT->sV[k].newDelta+1);
                                    // For the elements that are touched in the ascent stage it is necessary to
                                    // to reduce the values that they previously had.
                                    // In addition to this, the "parentVertex" that is connected to "startVertex", i.e.
                                    // the vertices of the new edge, needs to increase its betweenness centrality
                                    // following the new connection, without removing the old delta value.
                                    if(eAPT->sV[currElement].touched<0 && ( currElement!=parentVertex || k!=startVertex)){
                                        eAPT->sV[currElement].newDelta -=
                                                ((bc_t)tree->vArr[currElement].pathsToRoot/(bc_t)tree->vArr[k].pathsToRoot)*
                                                (bc_t)(tree->vArr[k].delta+1);
                                        
                                    }
                                }

                                        
                                // Vertices that did not move and that one of their neighbors move up(such that
                                // the vertices are now in the same level).
                                if(tree->vArr[k].level == tree->vArr[currElement].level && ((eAPT->sV[currElement].IMoved==1 && eAPT->sV[k].IMoved<0) )){
                                    // Checking to see if "k" has been touched before.
                                    //if(eAPT->sV[k].touched==0){
                                    if (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 0, -1)) {
                                            eAPT->sV[k].newDelta += tree->vArr[k].delta;
                                            
                                            upCounter++;
                                            // Marking element as touched in the ascent stage.
                                            eAPT->sV[k].touched=-2;
                                            __sync_bool_compare_and_swap(&depthSame, -1, tree->vArr[k].level);
                                            QueueSame[__sync_fetch_and_add(qEndSame_nxt, 1)] = k;
                                            eAPT->sV[k].newPathsToRoot += tree->vArr[k].pathsToRoot;
                                    }
                                    // Paths that previosul went through this vertex no longer go through them, thus the
                                    // shortest path count(BC) is reduced.
                                                                        
                                }

                                if(tree->vArr[k].level == tree->vArr[currElement].level && ((eAPT->sV[k].IMoved==1 && eAPT->sV[currElement].IMoved<0) )){
                                    eAPT->sV[currElement].newDelta -=
                                            ((bc_t)tree->vArr[currElement].pathsToRoot/(bc_t)tree->vArr[k].pathsToRoot)*
                                            (bc_t)(tree->vArr[k].delta+1);

                                }
                        }
                        STINGER_FORALL_EDGES_OF_VTX_END();

#if COUNT_TRAVERSALS==1
                        eAPT->dynamicTraverseEdgeCounter+=stinger_typed_outdegree(sStinger,currElement,0);
                        eAPT->dynamicTraverseVerticeCounter++;
#endif
                        if(currElement!=currRoot){
                                eAPT->sV[currElement].totalBC+=eAPT->sV[currElement].newDelta-tree->vArr[currElement].delta;
                        }        
                    }
                }
                qDownBIndex -= 2;

            }


            if (operation == 1) {
                *qStart = *qStart_nxt;
                *qEnd = *qEnd_nxt;
                *qStart_nxt = *qEnd;
                *qEnd_nxt = *qStart_nxt;

                int64_t thread_nums = cores;
                if (*qEnd - *qStart < cores) {
                    thread_nums = *qEnd - *qStart;
                }

                #pragma omp parallel num_threads(thread_nums)
                {
                    #pragma omp for
                    for (int64_t i = *qStart; i < *qEnd; i++) {
                        int64_t currElement = QueueUp[i];
                        int64_t levelCurrMinusOne = tree->vArr[currElement].level-1;
                        eAPT->sV[currElement].newEdgesBelow = 0;
                        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
                        {
                                uint64_t k = STINGER_EDGE_DEST;
                                // Checking that the vertices are in different levels.

                                if (tree->vArr[k].level == tree->vArr[currElement].level + 1) {

                                    if (eAPT->sV[k].touched == 0) {
                                        eAPT->sV[currElement].newEdgesBelow += tree->vArr[k].edgesBelow + 1;
                                    } else {
                                        eAPT->sV[currElement].newEdgesBelow += eAPT->sV[k].newEdgesBelow + 1;
                                    }
                                }

                                if(tree->vArr[k].level == levelCurrMinusOne){
                                        // Checking to see if "k" has been touched before.
                                        if (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 0, -1)) {
                                                eAPT->sV[k].newDelta += tree->vArr[k].delta;
                                                upCounter++;
                                                // Marking element as touched in the ascent stage.
                                                eAPT->sV[k].touched=-1;

                                                //QueueUp[QUpEnd] = k;
                                                __sync_bool_compare_and_swap(&depthUp, -1, tree->vArr[k].level);
                                                QueueUp[__sync_fetch_and_add(qEnd_nxt, 1)] = k;

                                                if (k != parentVertex)
                                                    eAPT->sV[k].newPathsToRoot += tree->vArr[k].pathsToRoot;
                                        }
                                }
                                
                                if (tree->vArr[k].level == tree->vArr[currElement].level + 1 && eAPT->sV[k].touched != 0) {
                                    eAPT->sV[currElement].newDelta +=
                                                ((bc_t)eAPT->sV[currElement].newPathsToRoot/(bc_t)eAPT->sV[k].newPathsToRoot)*
                                                (bc_t)(eAPT->sV[k].newDelta+1);
                                    // For the elements that are touched in the ascent stage it is necessary to
                                    // to reduce the values that they previously had.
                                    // In addition to this, the "parentVertex" that is connected to "startVertex", i.e.
                                    // the vertices of the new edge, needs to increase its betweenness centrality
                                    // following the new connection, without removing the old delta value.
                                    if(eAPT->sV[currElement].touched<0 && ( currElement!=parentVertex || k!=startVertex)){
                                        eAPT->sV[currElement].newDelta -=
                                                ((bc_t)tree->vArr[currElement].pathsToRoot/(bc_t)tree->vArr[k].pathsToRoot)*
                                                (bc_t)(tree->vArr[k].delta+1);
                                        
                                    }
                                }

                                        
                                // Vertices that did not move and that one of their neighbors move up(such that
                                // the vertices are now in the same level).
                                if(tree->vArr[k].level == tree->vArr[currElement].level && ((eAPT->sV[currElement].IMoved==1 && eAPT->sV[k].IMoved<0) )){
                                    // Checking to see if "k" has been touched before.
                                    //if(eAPT->sV[k].touched==0){
                                    if (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 0, -1)) {
                                            eAPT->sV[k].newDelta += tree->vArr[k].delta;
                                            
                                            upCounter++;
                                            // Marking element as touched in the ascent stage.
                                            eAPT->sV[k].touched=-2;
                                            //QueueSame[QSameEnd]=k;
                                            
                                            __sync_bool_compare_and_swap(&depthSame, -1, tree->vArr[k].level);
                                            QueueSame[__sync_fetch_and_add(qEndSame_nxt, 1)]=k;
                                            eAPT->sV[k].newPathsToRoot += tree->vArr[k].pathsToRoot;
                                    }
                                    // Paths that previosul went through this vertex no longer go through them, thus the
                                    // shortest path count(BC) is reduced.
                                                                        
                                }

                                if(tree->vArr[k].level == tree->vArr[currElement].level && ((eAPT->sV[k].IMoved==1 && eAPT->sV[currElement].IMoved<0) )){
                                    eAPT->sV[currElement].newDelta -=
                                            ((bc_t)tree->vArr[currElement].pathsToRoot/(bc_t)tree->vArr[k].pathsToRoot)*
                                            (bc_t)(tree->vArr[k].delta+1);

                                }
                        }
                        STINGER_FORALL_EDGES_OF_VTX_END();

#if COUNT_TRAVERSALS==1
                        eAPT->dynamicTraverseEdgeCounter+=stinger_typed_outdegree(sStinger,currElement,0);
                        eAPT->dynamicTraverseVerticeCounter++;
#endif
                        if(currElement!=currRoot){
                                eAPT->sV[currElement].totalBC+=eAPT->sV[currElement].newDelta-tree->vArr[currElement].delta;
                        } 
                    }
                }
                *qStart = *qEnd;
            }

            if (operation == 2) {
                *qStartSame = *qStartSame_nxt;
                *qEndSame = *qEndSame_nxt;
                *qStartSame_nxt = *qEndSame;
                *qEndSame_nxt = *qStartSame_nxt;

                int64_t thread_nums = cores;
                if (*qEndSame - *qStartSame < cores) {
                    thread_nums = *qEndSame - *qStartSame;
                }

                #pragma omp parallel num_threads(thread_nums)
                {
                    #pragma omp for
                    for (int64_t i = *qStartSame; i < *qEndSame; i++) {
                        int64_t currElement = QueueSame[i];
                        int64_t levelCurrMinusOne = tree->vArr[currElement].level-1;
                        eAPT->sV[currElement].newEdgesBelow = 0;
                        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
                        {
                                uint64_t k = STINGER_EDGE_DEST;
                                // Checking that the vertices are in different levels.

                                if (tree->vArr[k].level == tree->vArr[currElement].level + 1) {

                                    if (eAPT->sV[k].touched == 0) {
                                        eAPT->sV[currElement].newEdgesBelow += tree->vArr[k].edgesBelow + 1;
                                    } else {
                                        eAPT->sV[currElement].newEdgesBelow += eAPT->sV[k].newEdgesBelow + 1;
                                    }
                                }

                                if(tree->vArr[k].level == levelCurrMinusOne){
                                        // Checking to see if "k" has been touched before.
                                        if (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 0, -1)) {
                                                eAPT->sV[k].newDelta += tree->vArr[k].delta;
                                                upCounter++;
                                                // Marking element as touched in the ascent stage.
                                                eAPT->sV[k].touched=-1;

                                                __sync_bool_compare_and_swap(&depthUp, -1, tree->vArr[k].level);
                                                QueueUp[__sync_fetch_and_add(qEnd_nxt, 1)] = k;

                                                if (k != parentVertex)
                                                    eAPT->sV[k].newPathsToRoot += tree->vArr[k].pathsToRoot;
                                        }
                                }
                                
                                if (tree->vArr[k].level == tree->vArr[currElement].level + 1 && eAPT->sV[k].touched != 0) {
                                    eAPT->sV[currElement].newDelta +=
                                                ((bc_t)eAPT->sV[currElement].newPathsToRoot/(bc_t)eAPT->sV[k].newPathsToRoot)*
                                                (bc_t)(eAPT->sV[k].newDelta+1);
                                    // For the elements that are touched in the ascent stage it is necessary to
                                    // to reduce the values that they previously had.
                                    // In addition to this, the "parentVertex" that is connected to "startVertex", i.e.
                                    // the vertices of the new edge, needs to increase its betweenness centrality
                                    // following the new connection, without removing the old delta value.
                                    if(eAPT->sV[currElement].touched<0 && ( currElement!=parentVertex || k!=startVertex)){
                                        eAPT->sV[currElement].newDelta -=
                                                ((bc_t)tree->vArr[currElement].pathsToRoot/(bc_t)tree->vArr[k].pathsToRoot)*
                                                (bc_t)(tree->vArr[k].delta+1);
                                        
                                    }
                                }

                                        
                                // Vertices that did not move and that one of their neighbors move up(such that
                                // the vertices are now in the same level).
                                if(tree->vArr[k].level == tree->vArr[currElement].level && ((eAPT->sV[currElement].IMoved==1 && eAPT->sV[k].IMoved<0) )){
                                    // Checking to see if "k" has been touched before.
                                    //if(eAPT->sV[k].touched==0){
                                    if (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 0, -1)) {
                                            eAPT->sV[k].newDelta += tree->vArr[k].delta;
                                            
                                            upCounter++;
                                            // Marking element as touched in the ascent stage.
                                            eAPT->sV[k].touched=-2;
                                            __sync_bool_compare_and_swap(&depthSame, -1, tree->vArr[k].level);
                                            QueueSame[__sync_fetch_and_add(qEndSame_nxt, 1)]=k;
                                            eAPT->sV[k].newPathsToRoot += tree->vArr[k].pathsToRoot;
                                    }
                                    // Paths that previosul went through this vertex no longer go through them, thus the
                                    // shortest path count(BC) is reduced.
                                                                        
                                }

                                if(tree->vArr[k].level == tree->vArr[currElement].level && ((eAPT->sV[k].IMoved==1 && eAPT->sV[currElement].IMoved<0) )){
                                    eAPT->sV[currElement].newDelta -=
                                            ((bc_t)tree->vArr[currElement].pathsToRoot/(bc_t)tree->vArr[k].pathsToRoot)*
                                            (bc_t)(tree->vArr[k].delta+1);

                                }
                        }
                        STINGER_FORALL_EDGES_OF_VTX_END();

#if COUNT_TRAVERSALS==1
                        eAPT->dynamicTraverseEdgeCounter+=stinger_typed_outdegree(sStinger,currElement,0);
                        eAPT->dynamicTraverseVerticeCounter++;
#endif
                        if(currElement!=currRoot){
                                eAPT->sV[currElement].totalBC+=eAPT->sV[currElement].newDelta-tree->vArr[currElement].delta;
                        }                                             
                    }
                } 
                *qStartSame = *qEndSame;
            }

        }

        for(uint64_t c = 0; c <= qDownEndMarker; c++){
		uint64_t k=QueueDown[c];
		tree->vArr[k].delta=eAPT->sV[k].newDelta;
		tree->vArr[k].pathsToRoot=eAPT->sV[k].newPathsToRoot;
                tree->vArr[k].edgesBelow = eAPT->sV[k].newEdgesBelow;
                tree->vArr[k].edgesAbove = eAPT->sV[k].newEdgesAbove;
		eAPT->sV[k].diffPath=0;
		eAPT->sV[k].touched=0;
		eAPT->sV[k].newDelta=0.0;
		eAPT->sV[k].movementDelta=0;
		eAPT->sV[k].IMoved = -1;
		eAPT->sV[k].newPathsToRoot=0;
                eAPT->sV[k].newEdgesAbove = 0;
                eAPT->sV[k].newEdgesBelow = 0;
	}
        
        eAPT->sV[startVertex].newEdgesAbove = 0;
        eAPT->sV[parentVertex].newEdgesAbove = 0;

        for (uint64_t c = 0; c < *qEndSame; c++) {
		uint64_t k=QueueSame[c];
		tree->vArr[k].delta=eAPT->sV[k].newDelta;
                tree->vArr[k].edgesBelow = eAPT->sV[k].newEdgesBelow;
		eAPT->sV[k].diffPath=0;
		eAPT->sV[k].touched=0;
		eAPT->sV[k].newDelta=0.0;
		eAPT->sV[k].movementDelta=0;
		eAPT->sV[k].IMoved = -1;
		eAPT->sV[k].newPathsToRoot=0;
                eAPT->sV[k].newEdgesBelow = 0;
	}

	for(uint64_t c = 0; c < *qEnd; c++){
		uint64_t k=QueueUp[c];
		tree->vArr[k].delta=eAPT->sV[k].newDelta;
                tree->vArr[k].edgesBelow = eAPT->sV[k].newEdgesBelow;
		eAPT->sV[k].diffPath=0;
		eAPT->sV[k].touched=0;
		eAPT->sV[k].newDelta=0.0;
		eAPT->sV[k].movementDelta=0;
		eAPT->sV[k].IMoved = -1;
		eAPT->sV[k].newPathsToRoot=0;
                eAPT->sV[k].newEdgesBelow = 0;
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

	return;
}

// Case 2 
void removeEdgeWithoutMovementBrandes(bcForest* forest, struct stinger* sStinger, uint64_t currRoot,
                        uint64_t startVertex, uint64_t parentVertex, uint64_t deletedPathsFromRoot,
                        extraArraysPerThread *eAPT, uint64_t cores) {
    bcTree* tree = forest->forest[currRoot];

    int64_t NV = forest->NV;
    int64_t *Queue = eAPT->QueueSame;
    int64_t *QueueDown = eAPT->QueueDown;
    int64_t *QueueUp = eAPT->QueueUp;
    int64_t *QueueDownBorders = eAPT->Stack;
 
    eAPT->sV[startVertex].newEdgesBelow = tree->vArr[startVertex].edgesBelow;
    eAPT->sV[parentVertex].newEdgesBelow = tree->vArr[parentVertex].edgesBelow;
    eAPT->sV[startVertex].newEdgesAbove = tree->vArr[startVertex].edgesAbove;
    eAPT->sV[parentVertex].newEdgesAbove = tree->vArr[parentVertex].edgesAbove;
    eAPT->sV[startVertex].newPathsToRoot = tree->vArr[startVertex].pathsToRoot;
    eAPT->sV[startVertex].touched = 1;
    eAPT->sV[startVertex].newPathsToRoot -= deletedPathsFromRoot;
    eAPT->sV[startVertex].diffPath = deletedPathsFromRoot;
    eAPT->sV[startVertex].newEdgesAbove -= eAPT->sV[parentVertex].newEdgesAbove + 1;
    eAPT->sV[parentVertex].newEdgesBelow -= eAPT->sV[startVertex].newEdgesBelow + 1;

    QueueDown[0] = startVertex;
    
    //int64_t qDownStart=0,qDownEnd=1;
    int64_t *qDownStart = &(eAPT->qStart);
    int64_t *qDownEnd = &(eAPT->qEnd);
    int64_t *qDownStart_nxt = &(eAPT->qStart_nxt);
    int64_t *qDownEnd_nxt = &(eAPT->qEnd_nxt);

    int64_t qDownBIndex = 0;
    *qDownEnd = 1;
    *qDownStart_nxt = 1;
    *qDownEnd_nxt = 1;

    int64_t deepestLevel = tree->vArr[startVertex].level;
    
    queue_t* queue = eAPT->queue;
    level_node_t* levelIndices = eAPT->levelIndices;
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
                uint64_t currElement = QueueDown[i];

                if (currElement != startVertex) {
                    __sync_fetch_and_add(&(eAPT->sV[currElement].newEdgesAbove), tree->vArr[currElement].edgesAbove);
                }

                STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
                {
                    uint64_t k = STINGER_EDGE_DEST;

                    if (currElement != startVertex
                            && tree->vArr[currElement].level - 1 == tree->vArr[k].level
                            && tree->vArr[currElement].level >= tree->vArr[startVertex].level) {

                        if (eAPT->sV[k].touched != 0) {
                            __sync_fetch_and_add(&(eAPT->sV[currElement].newEdgesAbove), -tree->vArr[k].edgesAbove);
                            __sync_fetch_and_add(&(eAPT->sV[currElement].newEdgesAbove), eAPT->sV[k].newEdgesAbove);
                        }
                    }

                    // if this vertex has not been added yet
                    if ((tree->vArr[currElement].level + 1) == (tree->vArr[k].level)) {
                        if (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 0, currElement)) {
                            // Checking if a "deeper level" has been reached.
                            if(deepestLevel < tree->vArr[k].level)
                                deepestLevel = tree->vArr[k].level;

                            // insert this vertex into the BFS queue
                            QueueDown[__sync_fetch_and_add(qDownEnd_nxt, 1)] = k;
                           
                            __sync_fetch_and_add(&(eAPT->sV[k].newPathsToRoot), tree->vArr[k].pathsToRoot);

                            // indicate that it is in the next level of the BFS

                            // add new paths to root that go through current BFS Vertex
                            __sync_fetch_and_add(&(eAPT->sV[k].newPathsToRoot), -eAPT->sV[currElement].diffPath);
                            // pass on my new paths to root for its search
                            __sync_fetch_and_add(&(eAPT->sV[k].diffPath), eAPT->sV[currElement].diffPath);
                        }        
                        // otherwise if it has been touched, but is specifically in the next level
                        // of the search (meaning it has more than one edge to the current level)
                        //else if(eAPT->sV[k].touched == eAPT->sV[currElement].touched + 1){}
                        else if (eAPT->sV[k].touched != currElement) {
                            // add new paths to root that go through current BFS Vertex
                            __sync_fetch_and_add(&(eAPT->sV[k].newPathsToRoot), -eAPT->sV[currElement].diffPath);

                            // pass on my new paths to root for its search
                            __sync_fetch_and_add(&(eAPT->sV[k].diffPath), eAPT->sV[currElement].diffPath);
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
                        __sync_fetch_and_add(&(eAPT->sV[currElement].newEdgesBelow), tree->vArr[currElement].edgesBelow);
                    }
                    STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
                    {
                        uint64_t k = STINGER_EDGE_DEST;


                        if (currElement != parentVertex
                                && tree->vArr[currElement].level <= tree->vArr[parentVertex].level
                                && tree->vArr[k].level > tree->vArr[currElement].level) {


                            if (eAPT->sV[k].touched != 0) {
                                __sync_fetch_and_add(&(eAPT->sV[k].newEdgesBelow), -tree->vArr[k].edgesBelow);
                                __sync_fetch_and_add(&(eAPT->sV[k].newEdgesBelow), eAPT->sV[k].newEdgesBelow);
                            }
                        }

                        
                        if (tree->vArr[k].level == tree->vArr[parentVertex].level && __sync_bool_compare_and_swap(&(eAPT->sV[parentVertex].touched), 0, -1)) {
                            QueueUp[__sync_fetch_and_add(qDownEnd_nxt, 1)] = parentVertex;
                            __sync_fetch_and_add(&(eAPT->sV[parentVertex].newPathsToRoot), tree->vArr[parentVertex].pathsToRoot);
                            eAPT->sV[parentVertex].newDelta = tree->vArr[parentVertex].delta -
                                    ((bc_t)tree->vArr[parentVertex].pathsToRoot / (bc_t)tree->vArr[startVertex].pathsToRoot) *
                                    (bc_t)(tree->vArr[startVertex].delta + 1);


                        }
                        // Checking that the vertices are in different levels.
                        if (tree->vArr[k].level == (tree->vArr[currElement].level - 1)){
                            // Checking to see if "k" has been touched before.
                            if (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 0, -1)) {
                                eAPT->sV[k].newDelta=tree->vArr[k].delta;
                                // Marking element as touched in the ascent stage.
                                __sync_fetch_and_add(&(eAPT->sV[k].newPathsToRoot), tree->vArr[k].pathsToRoot);
                                QueueUp[__sync_fetch_and_add(qDownEnd_nxt, 1)] = k;
                            }
                        }

                        if (tree->vArr[k].level == tree->vArr[currElement].level + 1 && eAPT->sV[k].touched != 0) {

                            eAPT->sV[currElement].newDelta +=
                                ((bc_t)eAPT->sV[currElement].newPathsToRoot / (bc_t)eAPT->sV[k].newPathsToRoot) *
                                (bc_t)(eAPT->sV[k].newDelta + 1);

                            // For the elements that are touched in the ascent stage it is necessary to
                            // to reduce the values that they previously had.
                            // In addition to this, the "parentVertex" that is connected to "vertex", i.e.
                            // the vertices of the new edge, needs to increase its betweenness centrality
                            // following the new connection, without removing the old delta value.
                            if (eAPT->sV[currElement].touched < 0 && ( currElement != parentVertex || k != startVertex)){
                                eAPT->sV[currElement].newDelta -=
                                    ((bc_t)tree->vArr[currElement].pathsToRoot / (bc_t)tree->vArr[k].pathsToRoot) *
                                    (bc_t)(tree->vArr[k].delta + 1);

                            }
                        }
                    }
                    STINGER_FORALL_EDGES_OF_VTX_END();

#if COUNT_TRAVERSALS==1
                    eAPT->dynamicTraverseEdgeCounter += stinger_typed_outdegree(sStinger, currElement, 0);
                    eAPT->dynamicTraverseVerticeCounter++;
#endif
                    if (currElement != currRoot){
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
                    __sync_fetch_and_add(&(eAPT->sV[currElement].newEdgesBelow), tree->vArr[currElement].edgesBelow);
                }
                STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
                {
                    uint64_t k = STINGER_EDGE_DEST;


                    if (currElement != parentVertex
                            && tree->vArr[currElement].level <= tree->vArr[parentVertex].level
                            && tree->vArr[k].level > tree->vArr[currElement].level) {

                        if (eAPT->sV[k].touched != 0) {
                            __sync_fetch_and_add(&(eAPT->sV[currElement].newEdgesBelow), -tree->vArr[k].edgesBelow);
                            __sync_fetch_and_add(&(eAPT->sV[currElement].newEdgesBelow), eAPT->sV[k].newEdgesBelow);
                        }
                    }

                    
                    if (tree->vArr[k].level == tree->vArr[parentVertex].level && __sync_bool_compare_and_swap(&(eAPT->sV[parentVertex].touched), 0, -1)) {
                        QueueUp[__sync_fetch_and_add(qDownEnd_nxt, 1)] = parentVertex;
                        __sync_fetch_and_add(&(eAPT->sV[parentVertex].newPathsToRoot), tree->vArr[parentVertex].pathsToRoot);
                        eAPT->sV[parentVertex].newDelta = tree->vArr[parentVertex].delta -
                                ((bc_t)tree->vArr[parentVertex].pathsToRoot / (bc_t)tree->vArr[startVertex].pathsToRoot) *
                                (bc_t)(tree->vArr[startVertex].delta + 1);


                    }
                    // Checking that the vertices are in different levels.
                    if (tree->vArr[k].level == (tree->vArr[currElement].level - 1)){
                        // Checking to see if "k" has been touched before.
                        if (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 0, -1)) {
                            eAPT->sV[k].newDelta=tree->vArr[k].delta;
                            // Marking element as touched in the ascent stage.
                            __sync_fetch_and_add(&(eAPT->sV[k].newPathsToRoot), tree->vArr[k].pathsToRoot);
                            QueueUp[__sync_fetch_and_add(qDownEnd_nxt, 1)] = k;
                        }
                    }

                    if (tree->vArr[k].level == tree->vArr[currElement].level + 1 && eAPT->sV[k].touched != 0) {

                        eAPT->sV[currElement].newDelta +=
                            ((bc_t)eAPT->sV[currElement].newPathsToRoot / (bc_t)eAPT->sV[k].newPathsToRoot) *
                            (bc_t)(eAPT->sV[k].newDelta + 1);

                        // For the elements that are touched in the ascent stage it is necessary to
                        // to reduce the values that they previously had.
                        // In addition to this, the "parentVertex" that is connected to "vertex", i.e.
                        // the vertices of the new edge, needs to increase its betweenness centrality
                        // following the new connection, without removing the old delta value.
                        if (eAPT->sV[currElement].touched < 0 && ( currElement != parentVertex || k != startVertex)){
                            eAPT->sV[currElement].newDelta -=
                                ((bc_t)tree->vArr[currElement].pathsToRoot / (bc_t)tree->vArr[k].pathsToRoot) *
                                (bc_t)(tree->vArr[k].delta + 1);

                        }
                    }
                }
                STINGER_FORALL_EDGES_OF_VTX_END();

#if COUNT_TRAVERSALS==1
                eAPT->dynamicTraverseEdgeCounter += stinger_typed_outdegree(sStinger, currElement, 0);
                eAPT->dynamicTraverseVerticeCounter++;
#endif
                if (currElement != currRoot){
                    eAPT->sV[currElement].totalBC += eAPT->sV[currElement].newDelta - tree->vArr[currElement].delta;
                }
            }
        }
        *qDownStart = *qDownEnd;
    }
    

    for (uint64_t q = 0; q <= qDownEndMarker; q++) {
        uint64_t k = QueueDown[q];
        if (eAPT->sV[k].touched != 0) {
            tree->vArr[k].delta = eAPT->sV[k].newDelta;
            tree->vArr[k].pathsToRoot = eAPT->sV[k].newPathsToRoot;
        }
        tree->vArr[k].edgesAbove = eAPT->sV[k].newEdgesAbove;
        eAPT->sV[k].newEdgesAbove = 0;
        eAPT->sV[k].newEdgesBelow = 0;
        eAPT->sV[k].diffPath = 0;
        eAPT->sV[k].touched = 0;
        eAPT->sV[k].newDelta = 0.0;
        eAPT->sV[k].newPathsToRoot = 0;
    }

    eAPT->sV[startVertex].newEdgesAbove = 0;
    eAPT->sV[parentVertex].newEdgesAbove = 0;

    for (uint64_t q = 0; q < *qDownEnd; q++) {
        uint64_t k = QueueUp[q];
        if (eAPT->sV[k].touched != 0) {
            tree->vArr[k].delta = eAPT->sV[k].newDelta;
            tree->vArr[k].pathsToRoot = eAPT->sV[k].newPathsToRoot;
        }
        tree->vArr[k].edgesBelow = eAPT->sV[k].newEdgesBelow;
        eAPT->sV[k].newEdgesBelow = 0;
        eAPT->sV[k].newEdgesAbove = 0;
        eAPT->sV[k].diffPath = 0;
        eAPT->sV[k].touched = 0;
        eAPT->sV[k].newDelta = 0.0;
        eAPT->sV[k].newPathsToRoot = 0;
    }
    
    eAPT->sV[startVertex].newEdgesBelow = 0;
    eAPT->sV[parentVertex].newEdgesBelow = 0;

    queue->size = 0;

    eAPT->qStart = 0;
    eAPT->qEnd = 0;
    eAPT->qStart_nxt = 0;
    eAPT->qEnd_nxt = 0;
    
}

void moveDownTreeBrandes(bcForest* forest, struct stinger* sStinger, uint64_t currRoot,
                uint64_t startVertex, uint64_t parentVertex, extraArraysPerThread *eAPT, uint64_t cores) {
    bcTree* tree = forest->forest[currRoot];

    int64_t NV = forest->NV;
    int64_t *Queue = eAPT->QueueDown;
    int64_t *QueueUp = eAPT->QueueUp;
    int64_t *topQueue = eAPT->QueueSame;
    int64_t *QueueDownBorders = eAPT->Stack;
    int64_t *tqBorders = eAPT->tqBorders;

    int64_t *touchedVerticesDown = eAPT->touchedVerticesDown;
    int64_t *touchedVerticesUp = eAPT->touchedVerticesUp;
    queue_t* queue = eAPT->queue;
    level_node_t* levelIndices = eAPT->levelIndices;

    Queue[0] = startVertex;
    //int64_t tqStart=0, tqEnd=0;
    int64_t tvDownEnd = 0, tvUpEnd = 0;
    int64_t stopLevel = tree->vArr[startVertex].level;

    int64_t *qStart = &(eAPT->qStart);
    int64_t *qEnd = &(eAPT->qEnd);
    int64_t *qStart_nxt = &(eAPT->qStart_nxt);
    int64_t *qEnd_nxt = &(eAPT->qEnd_nxt);

    int64_t *tqStart = &(eAPT->tqStart);
    int64_t *tqEnd = &(eAPT->tqEnd);
    int64_t *tqStart_nxt = &(eAPT->tqStart_nxt);
    int64_t *tqEnd_nxt = &(eAPT->tqEnd_nxt);

    int64_t qDownBIndex = 0, tqBIndex = 0;
    *qEnd = 1;
    *qStart_nxt = 1;
    *qEnd_nxt = 1;

    *tqStart = 0;
    *tqEnd = 0;
    *tqStart_nxt = 0;
    *tqEnd_nxt = 0;

    eAPT->sV[startVertex].newLevel = INFINITY_MY;
    eAPT->sV[startVertex].newPathsToRoot = INFINITY_MY;
    eAPT->sV[startVertex].newDelta = 0.0;
    eAPT->sV[startVertex].newEdgesAbove = INFINITY_MY;
    eAPT->sV[startVertex].touched = 1;
    touchedVerticesDown[tvDownEnd++] = startVertex;

    int64_t deepestLevel=stopLevel;

    *qStart = 0;
    while (*qStart != *qEnd) {
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
                //int64_t currElement = Queue[(*qStart)++];
                int64_t currElement = Queue[i];

                STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement)
                {
                    int64_t k = STINGER_EDGE_DEST;
                    
                    if (tree->vArr[k].level == tree->vArr[currElement].level + 1 && __sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 0, -1)) {
                        __sync_fetch_and_add(&(eAPT->sV[k].newEdgesAbove), INFINITY_MY);
                        __sync_fetch_and_add(&(eAPT->sV[k].newLevel), INFINITY_MY);
                        __sync_fetch_and_add(&(eAPT->sV[k].newPathsToRoot), INFINITY_MY);
                        touchedVerticesDown[__sync_fetch_and_add(&tvDownEnd, 1)] = k;
                        Queue[__sync_fetch_and_add(qEnd_nxt, 1)] = k;
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

                STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement)
                {
                    int64_t l = STINGER_EDGE_DEST;

                    if (tree->vArr[l].level == tree->vArr[currElement].level - 1) {
                        if (eAPT->sV[l].touched == 0) {
                            parentOutsideSubtree = l;
                            parentPathsToRoot += tree->vArr[l].pathsToRoot;
                        } else {
                            parentPathsToRoot += tree->vArr[l].pathsToRoot - 1;
                        }
                        parentEdgesAbove += tree->vArr[l].edgesAbove + 1;
                    } else if (tree->vArr[l].level == tree->vArr[currElement].level) {
                        if (eAPT->sV[l].touched == 0) {
                            siblingOutsideSubtree = l;
                            siblingPathsToRoot += tree->vArr[l].pathsToRoot;
                        } else {
                            siblingPathsToRoot += tree->vArr[l].pathsToRoot - 1;
                        }
                        siblingEdgesAbove += tree->vArr[l].edgesAbove + 1;
                    }
                }
                STINGER_FORALL_EDGES_OF_VTX_END();

                if (parentOutsideSubtree) {

                    if (eAPT->sV[currElement].touched == 1 || eAPT->sV[currElement].touched == SIBLING_ANCHORED) {
                        topQueue[__sync_fetch_and_add(tqEnd_nxt, 1)] = currElement;
                        eAPT->sV[currElement].newLevel = tree->vArr[parentOutsideSubtree].level + 1;
                    } 

                    eAPT->sV[currElement].touched = PARENT_ANCHORED;
                    eAPT->sV[currElement].newDelta = 0.0;
                } else if (siblingOutsideSubtree) {

                    if (eAPT->sV[currElement].touched == 1) {
                        topQueue[__sync_fetch_and_add(tqEnd_nxt, 1)] = currElement;
                        eAPT->sV[currElement].newLevel = tree->vArr[siblingOutsideSubtree].level + 1;
                    }

                    if (eAPT->sV[currElement].touched != PARENT_ANCHORED) {
                        eAPT->sV[currElement].touched = SIBLING_ANCHORED;
                    }

                    eAPT->sV[currElement].newDelta = 0.0;
                } 

                STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement)
                {
                    int64_t k = STINGER_EDGE_DEST;
                    
                    if (tree->vArr[k].level == tree->vArr[currElement].level + 1 && (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), -1, -2)
                                                                                 || __sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 1, -2))) {

                        if (eAPT->sV[currElement].touched == PARENT_ANCHORED && __sync_bool_compare_and_swap(&(eAPT->sV[k].touched), -2, PARENT_ANCHORED)
                                                                             || __sync_bool_compare_and_swap(&(eAPT->sV[k].touched), SIBLING_ANCHORED, PARENT_ANCHORED)) {}

                        else if (eAPT->sV[currElement].touched == SIBLING_ANCHORED && __sync_bool_compare_and_swap(&(eAPT->sV[k].touched), -2, SIBLING_ANCHORED)) {}
                         
                        else {
                            __sync_fetch_and_add(&(eAPT->sV[k].touched), 3);
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
    // Phase 3 ------
    //qsort_r(topQueue, tqEnd, sizeof(int64_t), compare_levels, eAPT);
   
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
    // -------
        
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

    // Phase 4 --------
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

                        Queue[__sync_fetch_and_add(qEnd, 1)] = topQueue[i];
                        appendDS(queue, levelIndices, eAPT->sV[topQueue[i]].newLevel, topQueue[i]);
                    }
                }
            }
            tqLevel += 2;
        }
        *qStart_nxt = *qEnd;
        *qEnd_nxt = *qStart_nxt;

        int64_t thread_nums = cores;
        if (*qEnd - *qStart < cores) {
            thread_nums = *qEnd - *qStart;
        }
        #pragma omp parallel num_threads(thread_nums)
        {
            #pragma omp for
            for (int64_t i = *qStart; i < *qEnd; i++) {
                uint64_t currElement = Queue[i];

                eAPT->sV[currElement].newEdgesAbove = 0;
                if (deepestLevel < eAPT->sV[currElement].newLevel) {
                    deepestLevel = eAPT->sV[currElement].newLevel;
                }

                STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
                {
                    int64_t k = STINGER_EDGE_DEST;
                    if (eAPT->sV[k].newLevel > eAPT->sV[currElement].newLevel) {
                        // Checking if "k" has been found.
                        
                        __sync_bool_compare_and_swap(&(eAPT->sV[k].newLevel), INFINITY_MY, eAPT->sV[currElement].newLevel + 1);  
                        
                        if (__sync_bool_compare_and_swap(&(eAPT->sV[k].touched), 1, 5) ||
                            __sync_bool_compare_and_swap(&(eAPT->sV[k].touched), PARENT_ANCHORED, 5) ||
                            __sync_bool_compare_and_swap(&(eAPT->sV[k].touched), SIBLING_ANCHORED, 5)) {

                            Queue[__sync_fetch_and_add(qEnd_nxt, 1)] = k;
                            eAPT->sV[k].newDelta = 0.0;

                            if (deepestLevel < eAPT->sV[k].newLevel)
                                deepestLevel = eAPT->sV[k].newLevel;
                            appendDS(queue, levelIndices, eAPT->sV[k].newLevel, k);
                            if (eAPT->sV[k].touched == 0) {
                                printf("currRoot, currElement, k: %ld, %ld, %ld\n", currRoot, currElement, k);
                                printf("AAAAAAAAAAAAAAAAAAHHHHHHHHHHHHHHHHHHHHHHH\n");
                            }
                        } 
                    }

                    if (eAPT->sV[currElement].newLevel == tree->vArr[k].level + 1 && eAPT->sV[k].touched == 0) {
                        if (__sync_bool_compare_and_swap(&(eAPT->sV[currElement].newPathsToRoot), INFINITY_MY, tree->vArr[k].pathsToRoot)) {}
                        else { __sync_fetch_and_add(&(eAPT->sV[currElement].newPathsToRoot), tree->vArr[k].pathsToRoot); }

                        if (__sync_bool_compare_and_swap(&(eAPT->sV[currElement].newEdgesAbove), INFINITY_MY, tree->vArr[k].edgesAbove + 1)) {}
                        else { __sync_fetch_and_add(&(eAPT->sV[currElement].newEdgesAbove), tree->vArr[k].edgesAbove + 1); }
                    } else if (eAPT->sV[currElement].newLevel == eAPT->sV[k].newLevel + 1 && eAPT->sV[k].touched != 0) {

                        if (__sync_bool_compare_and_swap(&(eAPT->sV[currElement].newPathsToRoot), INFINITY_MY, eAPT->sV[k].newPathsToRoot)) {
                        } else {
                            __sync_fetch_and_add(&(eAPT->sV[currElement].newPathsToRoot), eAPT->sV[k].newPathsToRoot);
                        }
                        
                        if (__sync_bool_compare_and_swap(&(eAPT->sV[currElement].newEdgesAbove), INFINITY_MY, eAPT->sV[k].newEdgesAbove)) {
                        } else {
                            __sync_fetch_and_add(&(eAPT->sV[currElement].newEdgesAbove), eAPT->sV[k].newEdgesAbove + 1);
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
    // ---------------

    /// Phase 5 ----------
    *qEnd = 0;
    //for (int64_t lev = tree->vArr[startVertex].level; lev < NV; lev++) {

    // If it is not a case 4
    if (deepestLevel != INFINITY_MY) {
        for (int64_t lev = stopLevel; lev <= deepestLevel; lev++) {
        //for (int64_t lev = 0; lev < NV; lev++) {
            int64_t index = levelIndices[lev].front;
            int64_t levelEmpty = 1;
            while (index != -1) {
                levelEmpty = 0;
                queue_node_t *temp_node = queue->nodes + index;
                Queue[(*qEnd)++] = temp_node->data;
                index = temp_node->next;
            }
            levelIndices[lev].front = -1;
            levelIndices[lev].back = -1;
            //if (levelEmpty) {
            // break;
            //}
        }
        queue->size = 0;
    }
    // -----------

    // Phase 6 -------------
    // -----------

    //int64_t qUpStart = 0, qUpEnd = 0;
    int64_t *qUpStart = &(eAPT->qStartSame);
    int64_t *qUpEnd = &(eAPT->qEndSame);
    (*qEnd)--;

    //if (tree->vArr[startVertex].level == INFINITY_MY) {
    if (eAPT->sV[startVertex].newLevel == INFINITY_MY) {
        if (eAPT->sV[parentVertex].touched == 0) {
            eAPT->sV[parentVertex].newDelta = tree->vArr[parentVertex].delta -
                    ((bc_t)tree->vArr[parentVertex].pathsToRoot / (bc_t)tree->vArr[startVertex].pathsToRoot) *
                    (bc_t)(tree->vArr[startVertex].delta + 1);

            eAPT->sV[parentVertex].newPathsToRoot = tree->vArr[parentVertex].pathsToRoot;
            eAPT->sV[parentVertex].touched = -1;
            //QueueUp[qUpEnd++] = parentVertex;
            QueueUp[(*qUpEnd)++] = parentVertex;
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

    // Phase 7 --------------

    while (1) {
        // Removing last element from the queue
        int64_t currElement = -1;
        //if (*qEnd < 0 && qUpStart >= qUpEnd) {
        if (*qEnd < 0 && *qUpStart >= *qUpEnd) {
            break;
        } else if (*qEnd < 0) {
            //currElement = QueueUp[qUpStart++];
            currElement = QueueUp[(*qUpStart)++];
        } else if (qEnd >= 0) {
            //if (qUpEnd > 0) {
            if (*qUpEnd > 0) {
                int64_t qUpStartLevel = -1;
                //if (eAPT->sV[QueueUp[qUpStart]].newLevel == 0) {
                if (eAPT->sV[QueueUp[*qUpStart]].newLevel == 0) {
                    //qUpStartLevel = tree->vArr[QueueUp[qUpStart]].level;
                    qUpStartLevel = tree->vArr[QueueUp[*qUpStart]].level;
                } else {
                    //qUpStartLevel = eAPT->sV[QueueUp[qUpStart]].newLevel;
                    qUpStartLevel = eAPT->sV[QueueUp[*qUpStart]].newLevel;
                }

                int64_t qEndLevel = -1;
                if (eAPT->sV[Queue[*qEnd]].newLevel == 0) {
                    qEndLevel = tree->vArr[Queue[*qEnd]].level;
                } else {
                    qEndLevel = eAPT->sV[Queue[*qEnd]].newLevel;
                }

                //if (tree->vArr[QueueUp[qUpStart]].level < tree->vArr[Queue[qEnd]].level) {
                if (qUpStartLevel < qEndLevel) {
                    currElement = Queue[(*qEnd)--];
                } else {
                    //currElement = QueueUp[qUpStart++];
                    currElement = QueueUp[(*qUpStart)++];
                }
            } else {
                currElement = Queue[(*qEnd)--];
            }
        } else {
            printf("Should never be here.\n");
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

        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
        {
            uint64_t k = STINGER_EDGE_DEST;
           
            int64_t kLevel = -1;
            if (eAPT->sV[k].newLevel == 0) {
                kLevel = tree->vArr[k].level;
            } else {
                kLevel = eAPT->sV[k].newLevel;
            }

            //if (tree->vArr[k].level == tree->vArr[currElement].level + 1) {
            if (kLevel == currElementLevel + 1) {
                if (eAPT->sV[k].touched != 0) {
                    eAPT->sV[currElement].newEdgesBelow += eAPT->sV[k].newEdgesBelow + 1;
                } else if (eAPT->sV[k].touched == 0) {
                    eAPT->sV[currElement].newEdgesBelow += tree->vArr[k].edgesBelow + 1;
                }
            }

            //if (tree->vArr[k].level == tree->vArr[parentVertex].level) {
            if (kLevel == parentVertexLevel) {
                if (eAPT->sV[parentVertex].touched == 0) {
                    eAPT->sV[parentVertex].newDelta = tree->vArr[parentVertex].delta -
                            ((bc_t)tree->vArr[parentVertex].pathsToRoot / (bc_t)tree->vArr[startVertex].pathsToRoot) *
                            (bc_t)(tree->vArr[startVertex].delta + 1);

                    eAPT->sV[parentVertex].newPathsToRoot = tree->vArr[parentVertex].pathsToRoot;
                    eAPT->sV[parentVertex].touched = -1;
                    //QueueUp[qUpEnd++] = parentVertex;
                    QueueUp[(*qUpEnd)++] = parentVertex;
                }
            }
            
            // Checking that the vertices are in different levels.
            if (kLevel == currElementLevel - 1) {
            //if(tree->vArr[k].level == (tree->vArr[currElement].level-1)){

                // Checking to see if "k" has been touched before.
                if(eAPT->sV[k].touched==0) { // && tree->vArr[k].level < stopLevel){
                    eAPT->sV[k].newDelta=tree->vArr[k].delta;

                    // Marking element as touched in the ascent stage.
                    eAPT->sV[k].newPathsToRoot = tree->vArr[k].pathsToRoot;
                    eAPT->sV[k].touched=-1;
                    //QueueUp[qUpEnd++] = k;
                    QueueUp[(*qUpEnd)++] = k;
                }                 
                
                eAPT->sV[k].newDelta +=
                    ((bc_t)eAPT->sV[k].newPathsToRoot/(bc_t)eAPT->sV[currElement].newPathsToRoot)*
                    (bc_t)(eAPT->sV[currElement].newDelta+1);


                // For the elements that are touched in the ascent stage it is necessary to
                // to reduce the values that they previously had.
                // In addition to this, the "parentVertex" that is connected to "vertex", i.e.
                // the vertices of the new edge, needs to increase its betweenness centrality
                // following the new connection, without removing the old delta value.
                if(eAPT->sV[k].touched<0 && tree->vArr[k].level < tree->vArr[currElement].level) // && ( k!=parentVertex || currElement!=startVertex))
                {
                    eAPT->sV[k].newDelta -=
                        ((bc_t)tree->vArr[k].pathsToRoot/(bc_t)tree->vArr[currElement].pathsToRoot)*
                        (bc_t)(tree->vArr[currElement].delta+1);
                }
            }
             
        }
        STINGER_FORALL_EDGES_OF_VTX_END();
#if COUNT_TRAVERSALS==1
        eAPT->dynamicTraverseEdgeCounter += stinger_typed_outdegree(sStinger, currElement, 0);
        eAPT->dynamicTraverseVerticeCounter++;
#endif
        if(currElement!=currRoot){
            eAPT->sV[currElement].totalBC+=eAPT->sV[currElement].newDelta-tree->vArr[currElement].delta;
        }
    }
    // --------
    for (int64_t k = 0; k < tvDownEnd; k++) {
        int64_t vertex = touchedVerticesDown[k];
        tree->vArr[vertex].level = eAPT->sV[vertex].newLevel;
    }
    // Handles case where edge deletion creates new connected component.
    // Phase 8 ------
    if (tree->vArr[startVertex].level == INFINITY_MY){
        *qStart = 0;
        *qEnd = 1;
        Queue[0] = startVertex;
        
        eAPT->sV[startVertex].touched = -2;

        while (*qStart != *qEnd){
            int64_t currElement = Queue[(*qStart)++];
            
            eAPT->sV[currElement].totalBC -= tree->vArr[currElement].delta;

            tree->vArr[currElement].edgesBelow = 0;
            tree->vArr[currElement].edgesAbove = 0;
            eAPT->sV[currElement].newEdgesAbove = 0;
            eAPT->sV[currElement].newEdgesBelow = 0;
            tree->vArr[currElement].pathsToRoot = INFINITY_MY;
            eAPT->sV[currElement].newPathsToRoot = 0;
            STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement)
            {
                uint64_t k = STINGER_EDGE_DEST;

                if (eAPT->sV[k].touched != -2) {
                    eAPT->sV[k].touched = -2;
                    touchedVerticesUp[tvUpEnd++] = k;
                    Queue[(*qEnd)++] = k;
                }
            }
            STINGER_FORALL_EDGES_OF_VTX_END();
        }
    }
    // --------

    // Phase 9 -------
    for (int64_t q = 0; q < tvDownEnd; q++) {
        int64_t k = touchedVerticesDown[q];
        if(eAPT->sV[k].touched>0){
            tree->vArr[k].pathsToRoot=eAPT->sV[k].newPathsToRoot;
        }
        if(eAPT->sV[k].touched!=0){
            tree->vArr[k].delta=eAPT->sV[k].newDelta;
        }
        tree->vArr[k].edgesAbove = eAPT->sV[k].newEdgesAbove;
        eAPT->sV[k].diffPath = 0;
        eAPT->sV[k].touched = 0;
        eAPT->sV[k].newDelta = 0.0;
        eAPT->sV[k].newPathsToRoot = 0;
        eAPT->sV[k].newLevel = 0;
        eAPT->sV[k].newEdgesAbove = 0;
    }
    // -----

    // Phase 10 ----
    for (int64_t q = 0; q < tvUpEnd; q++) {
        int64_t k = touchedVerticesUp[q];
        if(eAPT->sV[k].touched>0){
            tree->vArr[k].pathsToRoot=eAPT->sV[k].newPathsToRoot;
        }
        if(eAPT->sV[k].touched!=0){
            tree->vArr[k].delta=eAPT->sV[k].newDelta;
        }
        tree->vArr[k].edgesBelow = eAPT->sV[k].newEdgesBelow;
        eAPT->sV[k].diffPath = 0;
        eAPT->sV[k].touched = 0;
        eAPT->sV[k].newDelta = 0.0;
        eAPT->sV[k].newPathsToRoot = 0;
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
    // -----


