#include "omp.h"
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#include "streaming_utils.h"
#include "streaming_insertions.h"
#include "timing_util.h"

#include "bcTreeDS.h"

#include "list.h"
#include "dsUtils.h"

#include "timer.h"

#ifdef COUNTING_INFO_ON
uint64_t countFullVertex = 0;
uint64_t countFullEdgeUp = 0;
uint64_t countStreamVertexDown = 0;
uint64_t countStreamVertexUp = 0;
uint64_t countStreamEdgeDown = 0;
uint64_t countStreamEdgeUp = 0;

#endif



void bfsBrandes(bcForest* forest, struct stinger* sStinger,extraArraysPerThread* eAPT)
{
	for(int64_t i=0; i<forest->NV; i++)
		forest->totalBC[i]=0.0;

	//  ODED-WARNING - USING OMP ON THIS ALGORITHM W/O ADDING ATOMICS/CRITICAL SECTIONS IS DANERGOUS TO YOUR HEALTH
	//	#pragma omp parallel for
	for(uint64_t i = 0; i < forest->NV; i++) {
		//		VB("\tTree[%ld]...\n",i)
		bfsBrandesPerTree(forest,sStinger,i, forest->totalBC,eAPT);
	}
	//	printf("bfsBrandes completed, forest=%ld\n",forest);
}


void bfsBrandesForApproxCase(bcForest* forest, struct stinger* sStinger, uint64_t* rootArrayForApproximation,
		uint64_t rootArraySizeForApproximation,extraArraysPerThread* eAPT)
{
	//	printf("bfsBrandesForApproxCase started, forest=%ld\n",forest);
	//bfsBrandes(forest, sStinger);

	for(int64_t i=0; i<forest->NV; i++)
		forest->totalBC[i]=0.0;

	for(uint64_t i = 0; i < rootArraySizeForApproximation; i++) {
		bfsBrandesPerTree(forest,sStinger,rootArrayForApproximation[i], forest->totalBC,eAPT);
	}

}


void bfsBrandesForApproxCaseParallel(bcForest* forest, struct stinger* sStinger, uint64_t* rootArrayForApproximation,
		uint64_t NK,extraArraysPerThread** eAPT,int64_t NT)
{
	omp_set_num_threads(NT);

	for(int64_t i=0; i<forest->NV; i++)
		forest->totalBC[i]=0.0;

#pragma omp parallel for schedule(dynamic,1)
	for(uint64_t i = 0; i < NK; i++) {

		int64_t thread = omp_get_thread_num();
		bfsBrandesPerTree(forest,sStinger,rootArrayForApproximation[i], forest->totalBC,eAPT[thread]);
		//            printf("%ld\n",rootArrayForApproximation[i]);
	}

#pragma omp parallel for
	for(uint64_t v=0;v<forest->NV;v++) {
		//printf("before_static: totalBC[%ld]: %lf\n", v, forest->totalBC[v]);
		//printf("%lf\n", forest->totalBC[v]);
		for(uint64_t t=0;t<NT;t++) {
                	forest->totalBC[v]+=eAPT[t]->sV[v].totalBC;
                        eAPT[t]->sV[v].totalBC = 0.0;
                }
		//printf("after_static: totalBC[%ld]: %lf\n", v, forest->totalBC[v]);
		//printf("totalBC[%ld]: %lf\n", v, forest->totalBC[v]);
        }
}


uint64_t bfsBrandesPerTree(bcForest* forest, struct stinger* sStinger,
		uint64_t currRoot, bc_t* totalBC,extraArraysPerThread* eAPT)
{
	bcTree* tree = forest->forest[currRoot];

	for(uint64_t j = 0; j < tree->NV; j++)
	{
		tree->vArr[j].level = INFINITY_MY;
		tree->vArr[j].pathsToRoot = INFINITY_MY;
		tree->vArr[j].delta = 0;
	}
	tree->vArr[currRoot].level = 0;
	tree->vArr[currRoot].pathsToRoot = 1;

	uint64_t* Stack= eAPT->Stack;
	uint64_t* Queue = eAPT->QueueDown;


	Queue[0] = currRoot;
	int64_t qStart=0,qEnd=1;
	int64_t sStart=0;
	int64_t k;
	// While queue is not empty
	while(qStart!=qEnd)
	{
		uint64_t currElement = Queue[qStart];
		Stack[sStart] = currElement;
		sStart++;
		qStart++;

		STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
		{
			k = STINGER_EDGE_DEST;

			// If this is a neighbor and has not been found
			if(tree->vArr[k].level > tree->vArr[currElement].level)
			{
				// Checking if "k" has been found.
				if(tree->vArr[k].level==INFINITY_MY)
				{
					tree->vArr[k].level = tree->vArr[currElement].level+1;
					Queue[qEnd++] = k;
					tree->vArr[k].delta=0;
				}

				if(tree->vArr[k].pathsToRoot == INFINITY_MY)
				{
					// k has not been found and therefore its paths to the roots are through its parent.
					tree->vArr[k].pathsToRoot = tree->vArr[currElement].pathsToRoot;
				}
				else
				{
					// k has been found and has multiple paths to the root as it has multiple parents.
					tree->vArr[k].pathsToRoot += tree->vArr[currElement].pathsToRoot;
				}
			}

		}
		STINGER_FORALL_EDGES_OF_VTX_END();

#if COUNT_TRAVERSALS==1
		eAPT->staticTraverseEdgeCounter+=stinger_typed_outdegree(sStinger,currElement,0);
		eAPT->staticTraverseVerticeCounter++;
#endif		

	}

	// Using Brandes algorithm to compute BC for a specific tree.
	// Essentially, use the stack which the elements are placed in depth-reverse order, to "climb" back
	// up the tree, all the way to the root.


	int64_t sEnd = sStart-1;

	while(sEnd>=0)
	{
		uint64_t currElement = Stack[sEnd];


		STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
		{
			k = STINGER_EDGE_DEST;


			// If this is a neighbor and has not been found
			if((tree->vArr[k].level == (tree->vArr[currElement].level-1)))
			{

				tree->vArr[k].delta +=
					((bc_t)tree->vArr[k].pathsToRoot/(bc_t)tree->vArr[currElement].pathsToRoot)*
					(bc_t)(tree->vArr[currElement].delta+1);


			}
		}
		STINGER_FORALL_EDGES_OF_VTX_END();
#if COUNT_TRAVERSALS==1
		eAPT->staticTraverseEdgeCounter+=stinger_typed_outdegree(sStinger,currElement,0);
		eAPT->staticTraverseVerticeCounter++;
#endif		

		if(currElement!=currRoot)
		{
			//printf("currElement, delta: %ld, %lf\n", currElement, tree->vArr[currElement].delta);
                        //forest->eAPT->sV[currElement].totalBC+=tree->vArr[currElement].delta; 
                        eAPT->sV[currElement].totalBC+=tree->vArr[currElement].delta;
		}

		sEnd--;
	}


	return -1;
}





void addEdgeWithoutMovementBrandes(bcForest* forest, struct stinger* sStinger,
		uint64_t currRoot, uint64_t startVertex, uint64_t parentVertex,
		uint64_t addedPathsToRoot,  extraArraysPerThread* eAPT)
{

	bcTree* tree = forest->forest[currRoot];

	uint64_t *QueueDown=eAPT->QueueDown;
	uint64_t *QueueUp=eAPT->QueueUp;

	int64_t NV = forest->NV;


	eAPT->sV[startVertex].newPathsToRoot = tree->vArr[startVertex].pathsToRoot;

	eAPT->sV[startVertex].touched = 1;
	eAPT->sV[startVertex].newPathsToRoot += addedPathsToRoot;
	eAPT->sV[startVertex].diffPath = addedPathsToRoot;

	QueueDown[0] = startVertex;
	int64_t qStart=0,qEnd=1;

	int64_t deepestLevel = tree->vArr[startVertex].level;
	int64_t intialLevel = tree->vArr[startVertex].level;



	// Starting BFS decent from "startVertex", down to all the vertices that have shortest paths through "startVertex"
	// All elements that will be touched will receive a positive value in their touched field.
	// In this implementation, "STACKS" are not used for the "moving up" stage. Rather, a multi-level queue is used.
	// Each level in the tree(max depth NV) has a queue and a counter specifiying how deep a specific deepth-queue is.
	// For simplicity, all elements are pushed both into the multi-level queue and into the regular queue which is used
	// for the BFS traversal.
	while(qStart<qEnd)
	{
		uint64_t currElement = QueueDown[qStart];
		int64_t levelCurrPlusOne = tree->vArr[currElement].level+1;
		int64_t touchedCurrPlusOne = eAPT->sV[currElement].touched+1;
		STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
		{
			uint64_t k = STINGER_EDGE_DEST;

			// if this vertex has not been added yet
			if(levelCurrPlusOne == (tree->vArr[k].level))
			{
				if(eAPT->sV[k].touched == 0)
				{
					// Checking if a "deeper level" has been reached.
					if(deepestLevel<tree->vArr[k].level)
					{
						deepestLevel=tree->vArr[k].level;
					}

					//NEW
					eAPT->sV[k].newPathsToRoot = tree->vArr[k].pathsToRoot;

					// insert this vertex into the BFS queue
					QueueDown[qEnd++] = k;

					// indicate that it is in the next level of the BFS
					//eAPT->sV[k].touched = eAPT->sV[currElement].touched + 1;
					eAPT->sV[k].touched = touchedCurrPlusOne;
					// add new paths to root that go through current BFS Vertex
					eAPT->sV[k].newPathsToRoot += eAPT->sV[currElement].diffPath;
					// pass on my new paths to root for its search
					eAPT->sV[k].diffPath = eAPT->sV[currElement].diffPath;


				}
				// otherwise if it has been touched, but is specifically in the next level
				// of the search (meaning it has more than one edge to the current level)
				else if(eAPT->sV[k].touched == touchedCurrPlusOne)
					//else if(eAPT->sV[k].touched == eAPT->sV[currElement].touched + 1)
				{
					// add new paths to root that go through current BFS Vertex
					eAPT->sV[k].newPathsToRoot += eAPT->sV[currElement].diffPath;
					// pass on my new paths to root for its search
					eAPT->sV[k].diffPath += eAPT->sV[currElement].diffPath;
				}
			}

		}
		STINGER_FORALL_EDGES_OF_VTX_END();


#if COUNT_TRAVERSALS==1
		eAPT->dynamicTraverseEdgeCounter+=stinger_typed_outdegree(sStinger,currElement,0);
		eAPT->dynamicTraverseVerticeCounter++;
#endif
		qStart++;
	}


	int64_t QUpStart=0,QUpEnd=0;
	//    int64_t myCont=1;
	int64_t currElement;
	qEnd--;

	//NEW 
	int64_t qDownEndMarker= qEnd;

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
	while(1)
	{
		if(qEnd<0 && QUpStart>=QUpEnd)
		{
			break;
		}
		else if(qEnd<0)
		{
			currElement=QueueUp[QUpStart++];
		}
		else if (qEnd>=0)
		{
			if(QUpEnd>0)
			{
				if(tree->vArr[QueueUp[QUpStart]].level < tree->vArr[QueueDown[qEnd]].level)
				{
					currElement=QueueDown[qEnd--];
				}
				else
				{
					currElement=QueueUp[QUpStart++];
				}
			}
			else
				currElement=QueueDown[qEnd--];

		}
		else
		{
			printf("should never be here\n");
		}

		int64_t levelCurrMinusOne = tree->vArr[currElement].level-1;
		STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
		{
			uint64_t k = STINGER_EDGE_DEST;

                        //printf("currRoot, currElement, k: %ld, %ld, %ld\n", currRoot, currElement, k);
    			// Checking that the vertices are in different levels.
			//if(tree->vArr[k].level == (tree->vArr[currElement].level-1))
			if(tree->vArr[k].level == levelCurrMinusOne)
			{
				// Checking to see if "k" has been touched before.
				if(eAPT->sV[k].touched==0)
				{
					eAPT->sV[k].newDelta=tree->vArr[k].delta;

					// Marking element as touched in the ascent stage.
					eAPT->sV[k].touched=-1;
					QueueUp[QUpEnd] = k;
					QUpEnd++;


					//NEW
					eAPT->sV[k].newPathsToRoot = tree->vArr[k].pathsToRoot;

				}

				eAPT->sV[k].newDelta +=
					((bc_t)eAPT->sV[k].newPathsToRoot/(bc_t)eAPT->sV[currElement].newPathsToRoot)*
					(bc_t)(eAPT->sV[currElement].newDelta+1);

                                // For the elements that are touched in the ascent stage it is necessary to
				// to reduce the values that they previously had.
				// In addition to this, the "parentVertex" that is connected to "vertex", i.e.
				// the vertices of the new edge, needs to increase its betweenness centrality
				// following the new connection, without removing the old delta value.
				if(eAPT->sV[k].touched<0 && ( k!=parentVertex  || currElement!=startVertex))
				{
					eAPT->sV[k].newDelta -=
						((bc_t)tree->vArr[k].pathsToRoot/(bc_t)tree->vArr[currElement].pathsToRoot)*
						(bc_t)(tree->vArr[currElement].delta+1);

	                        }
			}

		}
		STINGER_FORALL_EDGES_OF_VTX_END();

#if COUNT_TRAVERSALS==1
		eAPT->dynamicTraverseEdgeCounter+=stinger_typed_outdegree(sStinger,currElement,0);
		eAPT->dynamicTraverseVerticeCounter++;
		//printf("%ld,", eAPT->dynamicTraverseVerticeCounter);
#endif

		if(currElement!=currRoot)
		{
			//forest->eAPT->sV[currElement].totalBC+=eAPT->sV[currElement].newDelta-tree->vArr[currElement].delta;
			eAPT->sV[currElement].totalBC+=eAPT->sV[currElement].newDelta - tree->vArr[currElement].delta;
		}
	}
	/*
	   for(uint64_t k = 0; k < NV; k++)
	   {
	   if(eAPT->sV[k].touched!=0)
	   {
	   tree->vArr[k].delta=eAPT->sV[k].newDelta;
	   tree->vArr[k].pathsToRoot=eAPT->sV[k].newPathsToRoot;
	   }
	   }
	 */
	//New  
	for(uint64_t c = 0; c <= qDownEndMarker; c++)
	{
		uint64_t k=QueueDown[c];
		tree->vArr[k].delta=eAPT->sV[k].newDelta;
		tree->vArr[k].pathsToRoot=eAPT->sV[k].newPathsToRoot;
		eAPT->sV[k].diffPath=0;
		eAPT->sV[k].touched=0;
		eAPT->sV[k].newDelta=0.0;
		eAPT->sV[k].newPathsToRoot=0;
	}

	for(uint64_t c = 0; c < QUpEnd; c++)
	{
		uint64_t k=QueueUp[c];
		tree->vArr[k].delta=eAPT->sV[k].newDelta;
		eAPT->sV[k].diffPath=0;
		eAPT->sV[k].touched=0;
		eAPT->sV[k].newDelta=0.0;
		eAPT->sV[k].newPathsToRoot=0;
	}


	return;
}


void moveUpTreeBrandes(bcForest* forest, struct stinger* sStinger,
		uint64_t currRoot, uint64_t startVertex, uint64_t parentVertex,
		uint64_t prevDist,  extraArraysPerThread* eAPT)
{

	bcTree* tree = forest->forest[currRoot];

	int64_t NV = forest->NV;

	uint64_t *QueueDown=eAPT->QueueDown;
	uint64_t *QueueUp=eAPT->QueueUp;
	uint64_t *QueueSame=eAPT->QueueSame;

	list_ptr* multiLevelQueues = eAPT->multiLevelQueues;
	/*      	for(uint64_t k = 0; k < NV; k++) {

			if(eAPT->sV[k].touched!=0 || eAPT->sV[k].newPathsToRoot!=0 || eAPT->sV[k].movementDelta!=0 || eAPT->sV[k].IMoved!=-1 || eAPT->sV[k].diffPath!=0)
			printf("%ld,%ld\n", k, eAPT->sV[k].touched);
			} 
	 */
	//	for(uint64_t k = 0; k < NV; k++) {

	//		if(eAPT->sV[k].touched!=0)
	//			printf("%ld,%ld\n", k, eAPT->sV[k].touched);
	//		eAPT->sV[k].touched = 0;
	//		eAPT->sV[k].movementDelta=0;
	//		eAPT->sV[k].diffPath = 0;
	//		eAPT->sV[k].IMoved=-1;

	//		eAPT->sV[k].newPathsToRoot=tree->vArr[k].pathsToRoot;
	//		eAPT->sV[k].newDelta=0.0;
	//	}

	/*
	   for(uint64_t k = 0; k < NV; k++) {
	   eAPT->sV[k].newPathsToRoot=tree->vArr[k].pathsToRoot;
	   }
	 */

	//printf("%ld,",toc());
	//return;
	//tic();                 


	//NEW
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


	int64_t deepestLevel = tree->vArr[parentVertex].level+1;

	/*
#if PARENT_INFO_ON
emptyList(tree->parentList[startVertex]);
append(tree->parentList[startVertex],makeNode(parentVertex));
#endif
	 */

	// Starting BFS decent from "startVertex", down to all the vertices that have shortest paths through "startVertex"
	// All elements that will be touched will receive a positive value in their touched field.
	// In this implementation, "STACKS" are not used for the "moving up" stage. Rather, a multi-level queue is used.
	// Each level in the tree(max depth NV) has a queue and a counter specifiying how deep a specific deepth-queue is.
	// For simplicity, all elements are pushed both into the multi-level queue and into the regular queue which is used
	// for the BFS traversal.
	while(qStart != qEnd) {
		int64_t currElement = QueueDown[qStart];

		int64_t touchedCurrPlusOne = eAPT->sV[currElement].touched+1;
		STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
		{
			uint64_t k = STINGER_EDGE_DEST;

			if(eAPT->sV[k].touched == 0 )
			{
				// compute distance the adjacent vertex should be moved
				int64_t computedDelta = eAPT->sV[currElement].movementDelta - 
                                                    (tree->vArr[currElement].level - tree->vArr[k].level + 1);
				//eAPT->sV[k].touched = eAPT->sV[currElement].touched + 1;
				eAPT->sV[k].touched = touchedCurrPlusOne;
				// if the adjacent vertex should be moved, put it in the queue
				if(computedDelta > 0)
				{
					eAPT->sV[k].newPathsToRoot = eAPT->sV[currElement].diffPath;
					eAPT->sV[k].diffPath = eAPT->sV[currElement].diffPath;
					eAPT->sV[k].movementDelta = computedDelta;
					eAPT->sV[k].IMoved=1;
					QueueDown[qEnd++] = k;
				}
				// Vertex that will not be moved has been found.
				else if(computedDelta == 0)
				{      
					//NEW
					eAPT->sV[k].newPathsToRoot = tree->vArr[k].pathsToRoot;

					eAPT->sV[k].newPathsToRoot += eAPT->sV[currElement].diffPath;
					eAPT->sV[k].diffPath += eAPT->sV[currElement].diffPath;
					eAPT->sV[k].movementDelta = computedDelta;
					eAPT->sV[k].IMoved=0;
					QueueDown[qEnd++] = k;

				}
				// Vertex that the number of shortest path to the root does not change has been found.
				// This vertex is not marked as it might be touched on the way up.
				else
				{
					eAPT->sV[k].touched=0;
				}

				// if adjacent and in the next level
			}
			//else if(eAPT->sV[k].touched == eAPT->sV[currElement].touched + 1)
			else if(eAPT->sV[k].touched == touchedCurrPlusOne)
			{
				int64_t computedDelta = eAPT->sV[currElement].movementDelta - (tree->vArr[currElement].level - tree->vArr[k].level + 1);
				if(computedDelta >= 0)
				{
					eAPT->sV[k].newPathsToRoot += eAPT->sV[currElement].diffPath;
					eAPT->sV[k].diffPath += eAPT->sV[currElement].diffPath;
				}
			}

		}
		STINGER_FORALL_EDGES_OF_VTX_END();

#if COUNT_TRAVERSALS==1
		eAPT->dynamicTraverseEdgeCounter+=stinger_typed_outdegree(sStinger,currElement,0);
		eAPT->dynamicTraverseVerticeCounter++;
#endif 


		// move ourself and retire
		tree->vArr[currElement].level -= eAPT->sV[currElement].movementDelta;

		// Placing in multilevel queue only after it has moved in the tree.
		//        queueBFSTREE[tree->vArr[currElement].level][levelCounter[tree->vArr[currElement].level]]=currElement;
		//        levelCounter[tree->vArr[currElement].level]++;

		append(multiLevelQueues[tree->vArr[currElement].level],makeNode(currElement));



		// Checking if a "deeper level" has been reached.
		if(deepestLevel<tree->vArr[currElement].level)
			deepestLevel=tree->vArr[currElement].level;

		qStart++;
	}

	// Starting Multi-Level "BFS" ascent.
	qEnd=0;
	node_t* temp_node;
	for(int lev = tree->vArr[startVertex].level; lev<NV; lev++)
	{
		while (multiLevelQueues[lev]->size >0)
		{
			temp_node = getFirst(multiLevelQueues[lev]);
			QueueDown[qEnd++] = temp_node->id;
			deleteFirst(multiLevelQueues[lev]);

		}

	}
	qEnd--;
	//NEW 
	int64_t qDownEndMarker= qEnd;

	int64_t QUpStart=0,QUpEnd=0;
	int64_t QSameStart=0,QSameEnd=0;
	int64_t currElement=0;  //dummy initilization - variable will be initialized in function.
	int64_t upCounter=0;
	int64_t depthDown=tree->vArr[QueueDown[qEnd]].level,depthUp=-1,depthSame=-1;


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

	while(1)
	{
		if(depthDown==-1 && depthSame==-1 && depthUp==-1)
			break;
		if(depthDown>=depthSame && depthDown >=depthUp)
		{
			currElement=QueueDown[qEnd];
			qEnd--;
			if (qEnd<0)
				depthDown=-1;
			else
				depthDown=tree->vArr[QueueDown[qEnd]].level;
		}
		else if(depthUp>=depthSame && depthUp >=depthDown)
		{
			currElement=QueueUp[QUpStart];
			QUpStart++;
			if (QUpStart>=QUpEnd)
				depthUp=-1;
			else
				depthUp=tree->vArr[QueueUp[QUpStart]].level;
		}
		else if(depthDown<=depthSame && depthUp <=depthSame)
		{
			currElement=QueueSame[QSameStart];
			QSameStart++;
			if (QSameStart>=QSameEnd)
				depthSame=-1;
			else
				depthSame=tree->vArr[QueueSame[QSameStart]].level;
		}


		int64_t levelCurrMinusOne = tree->vArr[currElement].level-1;
		STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
		{
			uint64_t k = STINGER_EDGE_DEST;

			// Checking that the vertices are in different levels.
			//if(tree->vArr[k].level == (tree->vArr[currElement].level-1))
			if(tree->vArr[k].level == levelCurrMinusOne)
			{
				// Checking to see if "k" has been touched before.
				if(eAPT->sV[k].touched==0)
				{
					eAPT->sV[k].newDelta=tree->vArr[k].delta;

					upCounter++;
					// Marking element as touched in the ascent stage.
					eAPT->sV[k].touched=-1;

					QueueUp[QUpEnd] = k;
					if(depthUp==-1)
						depthUp=tree->vArr[QueueUp[QUpStart]].level;
					QUpEnd++;

					//NEWi
					eAPT->sV[k].newPathsToRoot = tree->vArr[k].pathsToRoot;

					//                        queueBFSTREE[tree->vArr[k].level][levelCounter[tree->vArr[k].level]]=k;
					//                        levelCounter[tree->vArr[k].level]++;

				}

				eAPT->sV[k].newDelta +=
					((bc_t)eAPT->sV[k].newPathsToRoot/(bc_t)eAPT->sV[currElement].newPathsToRoot)*
					(bc_t)(eAPT->sV[currElement].newDelta+1);

				// For the elements that are touched in the ascent stage it is necessary to
				// to reduce the values that they previously had.
				// In addition to this, the "parentVertex" that is connected to "startVertex", i.e.
				// the vertices of the new edge, needs to increase its betweenness centrality
				// following the new connection, without removing the old delta value.
				if(eAPT->sV[k].touched<0 && ( k!=parentVertex  || currElement!=startVertex))
				{
					eAPT->sV[k].newDelta -=
						((bc_t)tree->vArr[k].pathsToRoot/(bc_t)tree->vArr[currElement].pathsToRoot)*
						(bc_t)(tree->vArr[currElement].delta+1);

				}

			}
			// Vertices that did not move and that one of their neighbors move up(such that
			//    the vertices are now in the same level).
			else if(tree->vArr[k].level == tree->vArr[currElement].level && ((eAPT->sV[currElement].IMoved==1 && eAPT->sV[k].IMoved<0) ))
			{
				// Checking to see if "k" has been touched before.
				if(eAPT->sV[k].touched==0)
				{
					eAPT->sV[k].newDelta=tree->vArr[k].delta;

					upCounter++;
					// Marking element as touched in the ascent stage.
					eAPT->sV[k].touched=-2;
					QueueSame[QSameEnd]=k;
					if(depthSame==-1)
						depthSame=tree->vArr[QueueSame[QSameStart]].level;
					QSameEnd++;

					//NEW
					eAPT->sV[k].newPathsToRoot = tree->vArr[k].pathsToRoot;


					//                        queueBFSTREE[tree->vArr[k].level][levelCounter[tree->vArr[k].level]]=k;
					//                        levelCounter[tree->vArr[k].level]++;

				}

				// Paths that previosul went through this vertex no longer go through them, thus the
				// shortest path count(BC) is reduced.
				eAPT->sV[k].newDelta -=
					((bc_t)tree->vArr[k].pathsToRoot/(bc_t)tree->vArr[currElement].pathsToRoot)*
					(bc_t)(tree->vArr[currElement].delta+1);


			}

		}
		STINGER_FORALL_EDGES_OF_VTX_END();

#if COUNT_TRAVERSALS==1
		eAPT->dynamicTraverseEdgeCounter+=stinger_typed_outdegree(sStinger,currElement,0);
		eAPT->dynamicTraverseVerticeCounter++;
#endif


		if(currElement!=currRoot)
		{
			//forest->eAPT->sV[currElement].totalBC+=eAPT->sV[currElement].newDelta-tree->vArr[currElement].delta;
			eAPT->sV[currElement].totalBC+=eAPT->sV[currElement].newDelta-tree->vArr[currElement].delta;
		}


	}  



	for(uint64_t c = 0; c <= qDownEndMarker; c++)
	{
		uint64_t k=QueueDown[c];
		tree->vArr[k].delta=eAPT->sV[k].newDelta;
		tree->vArr[k].pathsToRoot=eAPT->sV[k].newPathsToRoot;
		eAPT->sV[k].diffPath=0;
		eAPT->sV[k].touched=0;
		eAPT->sV[k].newDelta=0.0;
		eAPT->sV[k].movementDelta=0;
		eAPT->sV[k].IMoved=-1;
		eAPT->sV[k].newPathsToRoot=0;
	}


	for(uint64_t c = 0; c < QSameEnd; c++)
	{
		uint64_t k=QueueSame[c];
		tree->vArr[k].delta=eAPT->sV[k].newDelta;
		eAPT->sV[k].diffPath=0;
		eAPT->sV[k].touched=0;
		eAPT->sV[k].newDelta=0.0;
		eAPT->sV[k].movementDelta=0;
		eAPT->sV[k].IMoved=-1;
		eAPT->sV[k].newPathsToRoot=0;
	}
	for(uint64_t c = 0; c < QUpEnd; c++)
	{
		uint64_t k=QueueUp[c];
		tree->vArr[k].delta=eAPT->sV[k].newDelta;
		eAPT->sV[k].diffPath=0;
		eAPT->sV[k].touched=0;
		eAPT->sV[k].newDelta=0.0;
		eAPT->sV[k].movementDelta=0;
		eAPT->sV[k].IMoved=-1;
		eAPT->sV[k].newPathsToRoot=0;
	}
	
	return;
}




