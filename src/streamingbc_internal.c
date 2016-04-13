#include "omp.h"
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#include "streamingbc_aux.h"
//#include "timer.h"

#define ANCHORED 3

void addEdgeWithoutMovementBrandes(bcForest* forest, struct stinger* sStinger,
		uint64_t currRoot, uint64_t startVertex, uint64_t parentVertex,
		uint64_t addedPathsToRoot,  extraArraysPerThread* eAPT){

	bcTree* tree = forest->forest[currRoot];
	uint64_t *QueueDown=eAPT->QueueDown;
	uint64_t *QueueUp=eAPT->QueueUp;
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
	int64_t qStart=0,qEnd=1;

	int64_t deepestLevel = tree->vArr[startVertex].level;
	int64_t intialLevel = tree->vArr[startVertex].level;

	// Starting BFS decent from "startVertex", down to all the vertices that have shortest paths through "startVertex"
	// All elements that will be touched will receive a positive value in their touched field.
	// In this implementation, "STACKS" are not used for the "moving up" stage. Rather, a multi-level queue is used.
	// Each level in the tree(max depth NV) has a queue and a counter specifiying how deep a specific deepth-queue is.
	// For simplicity, all elements are pushed both into the multi-level queue and into the regular queue which is used
	// for the BFS traversal.
	while(qStart<qEnd){
		uint64_t currElement = QueueDown[qStart];
		int64_t levelCurrPlusOne = tree->vArr[currElement].level+1;
		int64_t touchedCurrPlusOne = eAPT->sV[currElement].touched+1;

                if (currElement != startVertex)
                    eAPT->sV[currElement].newEdgesAbove = tree->vArr[currElement].edgesAbove;
		STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
		{
			uint64_t k = STINGER_EDGE_DEST;

			// if this vertex has not been added yet
			if (currElement != startVertex 
                                && tree->vArr[currElement].level >= tree->vArr[startVertex].level
                                && tree->vArr[k].level < tree->vArr[currElement].level) {
                            
                            if (eAPT->sV[k].touched == 0)
                                eAPT->sV[k].newEdgesAbove = tree->vArr[k].edgesAbove;

                            eAPT->sV[currElement].newEdgesAbove -= tree->vArr[k].edgesAbove;
                            eAPT->sV[currElement].newEdgesAbove += eAPT->sV[k].newEdgesAbove;
                        }
                        if(levelCurrPlusOne == (tree->vArr[k].level)){
                                
                                if(eAPT->sV[k].touched == 0){
					// Checking if a "deeper level" has been reached.
					if(deepestLevel<tree->vArr[k].level){
						deepestLevel=tree->vArr[k].level;
					}

					//NEW
					eAPT->sV[k].newPathsToRoot = tree->vArr[k].pathsToRoot;
					// insert this vertex into the BFS queue
					QueueDown[qEnd++] = k;
					// indicate that it is in the next level of the BFS
					eAPT->sV[k].touched = touchedCurrPlusOne;
					// add new paths to root that go through current BFS Vertex
					eAPT->sV[k].newPathsToRoot += eAPT->sV[currElement].diffPath;
					// pass on my new paths to root for its search
					eAPT->sV[k].diffPath = eAPT->sV[currElement].diffPath;
				}
				// otherwise if it has been touched, but is specifically in the next level
				// of the search (meaning it has more than one edge to the current level)
				else if(eAPT->sV[k].touched == touchedCurrPlusOne){
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
	while(1){
		if(qEnd<0 && QUpStart>=QUpEnd){
			break;
		}
		else if(qEnd<0){
			currElement=QueueUp[QUpStart++];
		}
		else if (qEnd>=0){
			if(QUpEnd>0){
				if(tree->vArr[QueueUp[QUpStart]].level < tree->vArr[QueueDown[qEnd]].level){
					currElement=QueueDown[qEnd--];
				}
				else{
					currElement=QueueUp[QUpStart++];
				}
			}
			else
				currElement=QueueDown[qEnd--];
		}
		else{
			printf("should never be here\n");
		}
		int64_t levelCurrMinusOne = tree->vArr[currElement].level-1;

                if (currElement != parentVertex)
                    eAPT->sV[currElement].newEdgesBelow = tree->vArr[currElement].edgesBelow;
		STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
		{
			uint64_t k = STINGER_EDGE_DEST;

			if (currElement != parentVertex 
                                && tree->vArr[currElement].level <= tree->vArr[parentVertex].level
                                && tree->vArr[k].level > tree->vArr[currElement].level) {
                            
                            if (eAPT->sV[k].touched == 0)
                                eAPT->sV[k].newEdgesBelow = tree->vArr[k].edgesBelow;
                            
                            eAPT->sV[currElement].newEdgesBelow -= tree->vArr[k].edgesBelow;
                            eAPT->sV[currElement].newEdgesBelow += eAPT->sV[k].newEdgesBelow;
                        }
                        
                        if(tree->vArr[k].level == levelCurrMinusOne){
				// Checking to see if "k" has been touched before.
				if(eAPT->sV[k].touched==0){
					eAPT->sV[k].newDelta=tree->vArr[k].delta;
					// Marking element as touched in the ascent stage.
					eAPT->sV[k].touched=-1;
					QueueUp[QUpEnd] = k;
					QUpEnd++;
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
				if(eAPT->sV[k].touched<0 && ( k!=parentVertex  || currElement!=startVertex)){
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
#endif
		if(currElement!=currRoot){
			eAPT->sV[currElement].totalBC+=eAPT->sV[currElement].newDelta - tree->vArr[currElement].delta;
		}
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
	}

	for(uint64_t c = 0; c < QUpEnd; c++){
		uint64_t k=QueueUp[c];
		tree->vArr[k].delta=eAPT->sV[k].newDelta;
		tree->vArr[k].edgesBelow = eAPT->sV[k].newEdgesBelow;
		eAPT->sV[k].diffPath=0;
		eAPT->sV[k].touched=0;
		eAPT->sV[k].newDelta=0.0;
                eAPT->sV[k].newEdgesBelow = 0;
		eAPT->sV[k].newPathsToRoot=0;
	}
	return;
}


void moveUpTreeBrandes(bcForest* forest, struct stinger* sStinger,
		uint64_t currRoot, uint64_t startVertex, uint64_t parentVertex,
		uint64_t prevDist,  extraArraysPerThread* eAPT){

	bcTree* tree = forest->forest[currRoot];
	int64_t NV = forest->NV;
	uint64_t *QueueDown=eAPT->QueueDown;
	uint64_t *QueueUp=eAPT->QueueUp;
	uint64_t *QueueSame=eAPT->QueueSame;

	list_ptr* multiLevelQueues = eAPT->multiLevelQueues;
        queue_t* queue = eAPT->queue;
        level_node_t* levelIndices = eAPT->levelIndices;

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

        eAPT->sV[startVertex].newEdgesAbove = eAPT->sV[parentVertex].newEdgesAbove + 1;
	int64_t deepestLevel = tree->vArr[parentVertex].level+1;

	// Starting BFS decent from "startVertex", down to all the vertices that have shortest paths through "startVertex"
	// All elements that will be touched will receive a positive value in their touched field.
	// In this implementation, "STACKS" are not used for the "moving up" stage. Rather, a multi-level queue is used.
	// Each level in the tree(max depth NV) has a queue and a counter specifiying how deep a specific deepth-queue is.
	// For simplicity, all elements are pushed both into the multi-level queue and into the regular queue which is used
	// for the BFS traversal.
	while(qStart != qEnd) {
		int64_t currElement = QueueDown[qStart];

		int64_t touchedCurrPlusOne = eAPT->sV[currElement].touched+1;
                eAPT->sV[currElement].newEdgesAbove = 0;
		STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
		{
			uint64_t k = STINGER_EDGE_DEST;

                        int64_t computedDelta = eAPT->sV[currElement].movementDelta -
                                            (tree->vArr[currElement].level - tree->vArr[k].level + 1);

                        int64_t newCurrLevel = tree->vArr[currElement].level - eAPT->sV[currElement].movementDelta;
                        int64_t newKLevel = tree->vArr[k].level - computedDelta;


                        if (computedDelta >= 0 && newKLevel < newCurrLevel) {
                            if (eAPT->sV[k].touched == 0)
                                eAPT->sV[k].newEdgesAbove = tree->vArr[k].edgesAbove;
                            eAPT->sV[currElement].newEdgesAbove += eAPT->sV[k].newEdgesAbove + 1;
                        } else if (computedDelta < 0 && tree->vArr[k].level < newCurrLevel) {
                            if (eAPT->sV[k].touched == 0)
                                eAPT->sV[k].newEdgesAbove = tree->vArr[k].edgesAbove;
                            eAPT->sV[currElement].newEdgesAbove += eAPT->sV[k].newEdgesAbove + 1;
                        }

			if(eAPT->sV[k].touched == 0 ){
				// compute distance the adjacent vertex should be moved
				int64_t computedDelta = eAPT->sV[currElement].movementDelta - 
                                                    (tree->vArr[currElement].level - tree->vArr[k].level + 1);
				eAPT->sV[k].touched = touchedCurrPlusOne;
				// if the adjacent vertex should be moved, put it in the queue
				if(computedDelta > 0){
					eAPT->sV[k].newPathsToRoot = eAPT->sV[currElement].diffPath;
					eAPT->sV[k].diffPath = eAPT->sV[currElement].diffPath;
					eAPT->sV[k].movementDelta = computedDelta;
					eAPT->sV[k].IMoved=1;
					QueueDown[qEnd++] = k;
				}
				// Vertex that will not be moved has been found.
				else if(computedDelta == 0){      
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
				else{
					eAPT->sV[k].touched=0;
				}

				// if adjacent and in the next level
			}
			//else if(eAPT->sV[k].touched == eAPT->sV[currElement].touched + 1)
			else if(eAPT->sV[k].touched == touchedCurrPlusOne){
				int64_t computedDelta = eAPT->sV[currElement].movementDelta - (tree->vArr[currElement].level - tree->vArr[k].level + 1);
				if(computedDelta >= 0){
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
                appendDS(queue, levelIndices, tree->vArr[currElement].level, currElement);  
		// Checking if a "deeper level" has been reached.
		if(deepestLevel<tree->vArr[currElement].level)
			deepestLevel=tree->vArr[currElement].level;

		qStart++;
	}

	// Starting Multi-Level "BFS" ascent.
	qEnd=0;
        queue_node_t* temp_node;
	for(int lev = tree->vArr[startVertex].level; lev<NV; lev++){
                temp_node = getFirstDS(queue, levelIndices, lev);
                while (temp_node != NULL) {
                        QueueDown[qEnd++] = temp_node->data;
                        deleteFirstDS(queue, levelIndices, lev);
                        temp_node = getFirstDS(queue, levelIndices, lev);
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

	while(1){
		if(depthDown==-1 && depthSame==-1 && depthUp==-1)
			break;
		if(depthDown>=depthSame && depthDown >=depthUp){
			currElement=QueueDown[qEnd];
			qEnd--;
			if (qEnd<0)
				depthDown=-1;
			else
				depthDown=tree->vArr[QueueDown[qEnd]].level;
		}
		else if(depthUp>=depthSame && depthUp >=depthDown){
			currElement=QueueUp[QUpStart];
			QUpStart++;
			if (QUpStart>=QUpEnd)
				depthUp=-1;
			else
				depthUp=tree->vArr[QueueUp[QUpStart]].level;
		}
		else if(depthDown<=depthSame && depthUp <=depthSame){
			currElement=QueueSame[QSameStart];
			QSameStart++;
			if (QSameStart>=QSameEnd)
				depthSame=-1;
			else
				depthSame=tree->vArr[QueueSame[QSameStart]].level;
		}


		int64_t levelCurrMinusOne = tree->vArr[currElement].level-1;
                eAPT->sV[currElement].newEdgesBelow = 0;
		STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
		{
			uint64_t k = STINGER_EDGE_DEST;
			// Checking that the vertices are in different levels.

                        if (tree->vArr[k].level == tree->vArr[currElement].level + 1) {

                            if (eAPT->sV[k].touched == 0)
                                eAPT->sV[k].newEdgesBelow = tree->vArr[k].edgesBelow;
                            eAPT->sV[currElement].newEdgesBelow += eAPT->sV[k].newEdgesBelow + 1;
                        }

			if(tree->vArr[k].level == levelCurrMinusOne){
				// Checking to see if "k" has been touched before.
				if(eAPT->sV[k].touched==0){
					eAPT->sV[k].newDelta=tree->vArr[k].delta;
					upCounter++;
					// Marking element as touched in the ascent stage.
					eAPT->sV[k].touched=-1;

					QueueUp[QUpEnd] = k;
					if(depthUp==-1)
						depthUp=tree->vArr[QueueUp[QUpStart]].level;
					QUpEnd++;

					eAPT->sV[k].newPathsToRoot = tree->vArr[k].pathsToRoot;
				}

				eAPT->sV[k].newDelta +=
					((bc_t)eAPT->sV[k].newPathsToRoot/(bc_t)eAPT->sV[currElement].newPathsToRoot)*
					(bc_t)(eAPT->sV[currElement].newDelta+1);

				// For the elements that are touched in the ascent stage it is necessary to
				// to reduce the values that they previously had.
				// In addition to this, the "parentVertex" that is connected to "startVertex", i.e.
				// the vertices of the new edge, needs to increase its betweenness centrality
				// following the new connection, without removing the old delta value.
				if(eAPT->sV[k].touched<0 && ( k!=parentVertex  || currElement!=startVertex)){
					eAPT->sV[k].newDelta -=
						((bc_t)tree->vArr[k].pathsToRoot/(bc_t)tree->vArr[currElement].pathsToRoot)*
						(bc_t)(tree->vArr[currElement].delta+1);
				}
			}
			// Vertices that did not move and that one of their neighbors move up(such that
			//    the vertices are now in the same level).
			else if(tree->vArr[k].level == tree->vArr[currElement].level && ((eAPT->sV[currElement].IMoved==1 && eAPT->sV[k].IMoved<0) )){
				// Checking to see if "k" has been touched before.
				if(eAPT->sV[k].touched==0){
					eAPT->sV[k].newDelta=tree->vArr[k].delta;

					upCounter++;
					// Marking element as touched in the ascent stage.
					eAPT->sV[k].touched=-2;
					QueueSame[QSameEnd]=k;
					if(depthSame==-1)
						depthSame=tree->vArr[QueueSame[QSameStart]].level;
					QSameEnd++;
					eAPT->sV[k].newPathsToRoot = tree->vArr[k].pathsToRoot;
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
		if(currElement!=currRoot){
			eAPT->sV[currElement].totalBC+=eAPT->sV[currElement].newDelta-tree->vArr[currElement].delta;
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
		eAPT->sV[k].IMoved=-1;
		eAPT->sV[k].newPathsToRoot=0;
	}


	for(uint64_t c = 0; c < QSameEnd; c++){
		uint64_t k=QueueSame[c];
		tree->vArr[k].delta=eAPT->sV[k].newDelta;
                tree->vArr[k].edgesBelow = eAPT->sV[k].newEdgesBelow;
		eAPT->sV[k].diffPath=0;
		eAPT->sV[k].touched=0;
		eAPT->sV[k].newDelta=0.0;
		eAPT->sV[k].movementDelta=0;
		eAPT->sV[k].IMoved=-1;
		eAPT->sV[k].newPathsToRoot=0;
	}
	for(uint64_t c = 0; c < QUpEnd; c++){
		uint64_t k=QueueUp[c];
		tree->vArr[k].delta=eAPT->sV[k].newDelta;
                tree->vArr[k].edgesBelow = eAPT->sV[k].newEdgesBelow;
		eAPT->sV[k].diffPath=0;
		eAPT->sV[k].touched=0;
		eAPT->sV[k].newDelta=0.0;
		eAPT->sV[k].movementDelta=0;
		eAPT->sV[k].IMoved=-1;
		eAPT->sV[k].newPathsToRoot=0;
	}
        queue->size = 0;	
	return;
}


// Case 2
void removeEdgeWithoutMovementBrandes(bcForest* forest, struct stinger* sStinger, uint64_t currRoot, 
                        uint64_t startVertex, uint64_t parentVertex, uint64_t deletedPathsFromRoot, 
                        extraArraysPerThread *eAPT) {
    bcTree* tree = forest->forest[currRoot];

    int64_t NV = forest->NV;
    int64_t *Queue = eAPT->QueueSame; 
    int64_t *QueueDown = eAPT->QueueDown;
    int64_t *QueueUp = eAPT->QueueUp;
 
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
    
    int64_t qDownStart=0,qDownEnd=1;
    int64_t deepestLevel = tree->vArr[startVertex].level;
    
    queue_t* queue = eAPT->queue;
    level_node_t* levelIndices = eAPT->levelIndices;
    // Starting BFS decent from "startVertex", down to all the vertices that have shortest paths through "startVertex"
    // All elements that will be touched will receive a positive value in their touched field.
    // In this implementation, "STACKS" are not used for the "moving up" stage. Rather, a multi-level queue is used.
    // Each level in the tree(max depth NV) has a queue and a counter specifiying how deep a specific deepth-queue is.
    // For simplicity, all elements are pushed both into the multi-level queue and into the regular queue which is used
    // for the BFS traversal.
    while(qDownStart!=qDownEnd){
        uint64_t currElement = QueueDown[qDownStart];

        if (currElement != startVertex)
            eAPT->sV[currElement].newEdgesAbove = tree->vArr[currElement].edgesAbove;

        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
        {
            uint64_t k = STINGER_EDGE_DEST;

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
    	        if(eAPT->sV[k].touched == 0) {
        	    // Checking if a "deeper level" has been reached.
        	    if(deepestLevel < tree->vArr[k].level)
                        deepestLevel = tree->vArr[k].level;

        	    // insert this vertex into the BFS queue
        	    QueueDown[qDownEnd++] = k;
                   
                    eAPT->sV[k].newPathsToRoot = tree->vArr[k].pathsToRoot; 
                    // indicate that it is in the next level of the BFS
                    eAPT->sV[k].touched = eAPT->sV[currElement].touched + 1;
                    // add new paths to root that go through current BFS Vertex
                    eAPT->sV[k].newPathsToRoot -= eAPT->sV[currElement].diffPath;
                    // pass on my new paths to root for its search
                    eAPT->sV[k].diffPath = eAPT->sV[currElement].diffPath;
                }        
                // otherwise if it has been touched, but is specifically in the next level
                // of the search (meaning it has more than one edge to the current level)
                else if(eAPT->sV[k].touched == eAPT->sV[currElement].touched + 1){
                    // add new paths to root that go through current BFS Vertex
                    eAPT->sV[k].newPathsToRoot -= eAPT->sV[currElement].diffPath;
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
	qDownStart++;
    }

    // The parent vertex needs to be placed in the queue for the dependency accumulation stage.
    // Also, it no longer has a child and so the delta from the child needs to be removed.
    
    int64_t qUpStart = 0, qUpEnd = 0;
    qDownEnd--;
    int64_t qDownEndMarker = qDownEnd;
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
     
    while (1){
        
        uint64_t currElement = -1;
        if (qDownEnd < 0 && qUpStart >= qUpEnd) {
            break;
        } else if (qDownEnd < 0) {
            currElement = QueueUp[qUpStart++];
        } else if (qDownEnd >= 0) {
            if (qUpEnd > 0) {
                if (tree->vArr[QueueUp[qUpStart]].level < tree->vArr[QueueDown[qDownEnd]].level) {
                    currElement = QueueDown[qDownEnd--];
                } else {
                    currElement = QueueUp[qUpStart++];
                }
            } else {
                currElement = QueueDown[qDownEnd--];
            }
        } else {
            printf("Should never be here.\n");
        }

        if (currElement != parentVertex)
            eAPT->sV[currElement].newEdgesBelow = tree->vArr[currElement].edgesBelow;
        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
        {
            uint64_t k = STINGER_EDGE_DEST;


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
                eAPT->sV[parentVertex].newPathsToRoot = tree->vArr[parentVertex].pathsToRoot;
                eAPT->sV[parentVertex].newDelta = tree->vArr[parentVertex].delta -
                        ((bc_t)tree->vArr[parentVertex].pathsToRoot / (bc_t)tree->vArr[startVertex].pathsToRoot) *
                        (bc_t)(tree->vArr[startVertex].delta + 1);


            }
            // Checking that the vertices are in different levels.
            if (tree->vArr[k].level == (tree->vArr[currElement].level - 1)){
                // Checking to see if "k" has been touched before.
                if (eAPT->sV[k].touched == 0){
                    eAPT->sV[k].newDelta=tree->vArr[k].delta;
                    // Marking element as touched in the ascent stage.
                    eAPT->sV[k].touched=-1;
                    eAPT->sV[k].newPathsToRoot = tree->vArr[k].pathsToRoot;
                    QueueUp[qUpEnd++] = k;
                }

                
                eAPT->sV[k].newDelta +=
                    ((bc_t)eAPT->sV[k].newPathsToRoot / (bc_t)eAPT->sV[currElement].newPathsToRoot) *
                    (bc_t)(eAPT->sV[currElement].newDelta + 1);

                // For the elements that are touched in the ascent stage it is necessary to
                // to reduce the values that they previously had.
                // In addition to this, the "parentVertex" that is connected to "vertex", i.e.
                // the vertices of the new edge, needs to increase its betweenness centrality
                // following the new connection, without removing the old delta value.
                if (eAPT->sV[k].touched < 0 && ( k != parentVertex  || currElement != startVertex)){
                    eAPT->sV[k].newDelta -=
                        ((bc_t)tree->vArr[k].pathsToRoot / (bc_t)tree->vArr[currElement].pathsToRoot) * 
                        (bc_t)(tree->vArr[currElement].delta + 1);
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

    for (uint64_t q = 0; q < qUpEnd; q++) {
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

    queue->size = 0;
}

// Case 3
void moveDownTreeBrandes(bcForest* forest, struct stinger* sStinger, uint64_t currRoot, 
                uint64_t startVertex, uint64_t parentVertex, extraArraysPerThread *eAPT)
{
    bcTree* tree = forest->forest[currRoot];

    int64_t NV = forest->NV;
    int64_t *Queue = eAPT->QueueDown;
    int64_t *QueueUp = eAPT->QueueUp;
    int64_t *topQueue = eAPT->QueueSame;

    int64_t *touchedVerticesDown = eAPT->touchedVerticesDown;
    int64_t *touchedVerticesUp = eAPT->touchedVerticesUp;
    queue_t* queue = eAPT->queue;
    level_node_t* levelIndices = eAPT->levelIndices;

    Queue[0] = startVertex;
    int64_t qStart=0,qEnd=1, tqStart=0, tqEnd=0, tvDownEnd = 0, tvUpEnd = 0;
    int64_t stopLevel = tree->vArr[startVertex].level;

    eAPT->sV[startVertex].newLevel = INFINITY_MY;
    eAPT->sV[startVertex].newPathsToRoot = INFINITY_MY;
    eAPT->sV[startVertex].newDelta = 0.0;
    eAPT->sV[startVertex].newEdgesAbove = tree->vArr[startVertex].edgesAbove;
    eAPT->sV[startVertex].touched = 1;
    touchedVerticesDown[tvDownEnd++] = startVertex;  

    int64_t deepestLevel=stopLevel;

    #if 0
    // Phase 1 ---
    STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,startVertex)
    {
        int64_t k = STINGER_EDGE_DEST;

        if (tree->vArr[k].level == tree->vArr[startVertex].level){
            eAPT->sV[k].newLevel = tree->vArr[k].level;
            eAPT->sV[k].newPathsToRoot = tree->vArr[k].pathsToRoot;
            topQueue[tqEnd++] = k;
        }
        else{
            eAPT->sV[k].newLevel = INFINITY_MY;
            eAPT->sV[k].newPathsToRoot = INFINITY_MY;
            if(deepestLevel < tree->vArr[k].level)
                deepestLevel = tree->vArr[k].level;
       }
        Queue[qEnd++] = k;
        eAPT->sV[k].newEdgesAbove = tree->vArr[k].edgesAbove;
        eAPT->sV[k].touched = 1;
        touchedVerticesDown[tvDownEnd++] = k;  
        eAPT->sV[k].newDelta = 0.0;
    }
    STINGER_FORALL_EDGES_OF_VTX_END();
    // --------------
    #endif
    
    
    #if 0
    int64_t startvSiblingOutsideTree = 0;
    int64_t startvSiblingPathsToRoot = 0;

    STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, startVertex)
    {
        int64_t l = STINGER_EDGE_DEST;
        
        if (tree->vArr[l].level == tree->vArr[startVertex].level && eAPT->sV[l].touched == 0) {
            startvSiblingOutsideTree = l;
            startvSiblingPathsToRoot += tree->vArr[l].pathsToRoot;
        }
    }
    STINGER_FORALL_EDGES_OF_VTX_END();

    if (startvSiblingOutsideTree) {
        eAPT->sV[startVertex].touched = ANCHORED;
        eAPT->sV[startVertex].newLevel = tree->vArr[startvSiblingOutsideTree].level + 1;
        eAPT->sV[startVertex].newPathsToRoot = startvSiblingPathsToRoot;
        topQueue[tqEnd++] = startVertex;
    }
    #endif

    qStart = 0;
    while (qStart != qEnd) {
        int64_t currElement = Queue[qStart++];

        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement)
        {
            int64_t k = STINGER_EDGE_DEST;
            
            if (tree->vArr[k].level == tree->vArr[currElement].level + 1 && eAPT->sV[k].touched == 0) {
                eAPT->sV[k].newEdgesAbove = tree->vArr[k].edgesAbove; 
                eAPT->sV[k].newLevel = INFINITY_MY;
                eAPT->sV[k].newPathsToRoot = INFINITY_MY;
                touchedVerticesDown[tvDownEnd++] = k;
                Queue[qEnd++] = k;
                eAPT->sV[k].newDelta = 0.0;

            }
        }
        STINGER_FORALL_EDGES_OF_VTX_END();

        int64_t parentOutsideSubtree = 0;
        int64_t siblingOutsideSubtree = 0;
        int64_t parentPathsToRoot = 0;
        int64_t siblingPathsToRoot = 0; 

        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement)
        {
            int64_t l = STINGER_EDGE_DEST;
            if (tree->vArr[l].level == tree->vArr[currElement].level - 1 && eAPT->sV[l].touched == 0) {
                parentOutsideSubtree = l;
                parentPathsToRoot += tree->vArr[l].pathsToRoot;
            } else if (tree->vArr[l].level == tree->vArr[currElement].level && eAPT->sV[l].touched == 0) {
                siblingOutsideSubtree = l;
                siblingPathsToRoot += tree->vArr[l].pathsToRoot;
            }
        }
        STINGER_FORALL_EDGES_OF_VTX_END();

        if (parentOutsideSubtree) {
            eAPT->sV[currElement].touched = ANCHORED;
            eAPT->sV[currElement].newLevel = tree->vArr[parentOutsideSubtree].level + 1;
            eAPT->sV[currElement].newPathsToRoot = parentPathsToRoot;
            eAPT->sV[currElement].newEdgesAbove = tree->vArr[parentOutsideSubtree].edgesAbove + 1;
            topQueue[tqEnd++] = currElement;
            eAPT->sV[currElement].newDelta = 0.0;
        } else if (siblingOutsideSubtree) {
            if (eAPT->sV[currElement].touched != ANCHORED) {
                topQueue[tqEnd++] = currElement;
                eAPT->sV[currElement].newLevel = tree->vArr[siblingOutsideSubtree].level + 1;
            }
            eAPT->sV[currElement].touched = ANCHORED;
            if (eAPT->sV[currElement].newPathsToRoot == INFINITY_MY) {
                eAPT->sV[currElement].newPathsToRoot = siblingPathsToRoot;
            } else {
                eAPT->sV[currElement].newPathsToRoot += siblingPathsToRoot;
            }
            eAPT->sV[currElement].newEdgesAbove = tree->vArr[siblingOutsideSubtree].edgesAbove + 1;
            eAPT->sV[currElement].newDelta = 0.0;
        } 

        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement)
        {
            int64_t k = STINGER_EDGE_DEST;
            
            if (tree->vArr[k].level == tree->vArr[currElement].level + 1 && eAPT->sV[k].touched == 0) {
                if (eAPT->sV[currElement].touched == ANCHORED) {
                    eAPT->sV[k].touched = ANCHORED;
                } else {
                    eAPT->sV[k].touched = 1;
                }
            }
        }
        STINGER_FORALL_EDGES_OF_VTX_END();


    }
    #if 0
    qStart = 1; 
    // Phase 2 --------------
    while (qStart != qEnd){
        int64_t currElement = Queue[qStart++];

        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
        {
            int64_t k = STINGER_EDGE_DEST;
               
            if (eAPT->sV[k].touched == 0){
                eAPT->sV[k].newEdgesAbove = tree->vArr[k].edgesAbove;
                if (tree->vArr[k].level == stopLevel && tree->vArr[currElement].level == (stopLevel + 1)){
                    eAPT->sV[k].newLevel = tree->vArr[k].level;
                    eAPT->sV[k].newPathsToRoot=tree->vArr[k].pathsToRoot;
                    topQueue[tqEnd++]=k;
                    eAPT->sV[k].touched=1;
                    touchedVerticesDown[tvDownEnd++] = k;  
                    Queue[qEnd++]=k;
                    eAPT->sV[k].newDelta=0.0;
                }
                else if (tree->vArr[k].level>stopLevel){
                    eAPT->sV[k].newLevel = INFINITY_MY;
                    eAPT->sV[k].newPathsToRoot = INFINITY_MY;
                    eAPT->sV[k].touched = 1;
                    touchedVerticesDown[tvDownEnd++] = k;  
                    Queue[qEnd++] = k;
                    eAPT->sV[k].newDelta = 0.0;
                }
            }
        }
        STINGER_FORALL_EDGES_OF_VTX_END();

#if COUNT_TRAVERSALS==1
        eAPT->dynamicTraverseEdgeCounter += stinger_typed_outdegree(sStinger, currElement, 0);
        eAPT->dynamicTraverseVerticeCounter++;
#endif
    }
    //--------
    #endif

    qEnd = tqEnd;
    qStart = 0;
    // Phase 3 ------
    for(int64_t t = 0; t < qEnd; t++){
        Queue[t] = topQueue[t];
        //eAPT->sV[t].newDelta=0.0;
        //appendDS(queue, levelIndices, stopLevel, topQueue[t]);
        appendDS(queue, levelIndices, eAPT->sV[topQueue[t]].newLevel, topQueue[t]);
    }
    // ------- 

    // While queue is not empty

    // Phase 4 --------
    while (qStart != qEnd){
        uint64_t currElement = Queue[qStart++];

        eAPT->sV[currElement].newEdgesAbove = 0;
        if (deepestLevel < eAPT->sV[currElement].newLevel) {
            deepestLevel = eAPT->sV[currElement].newLevel;
        }

        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
        {
            uint64_t k = STINGER_EDGE_DEST;
            if (eAPT->sV[k].newEdgesAbove == 0) {
                eAPT->sV[k].newEdgesAbove = tree->vArr[k].edgesAbove;
            }
            
            if (eAPT->sV[k].newLevel > eAPT->sV[currElement].newLevel) {
                // Checking if "k" has been found.
                
                if (eAPT->sV[k].newLevel == INFINITY_MY) {
                    eAPT->sV[k].newLevel = eAPT->sV[currElement].newLevel + 1;
                    Queue[qEnd++] = k;
                    eAPT->sV[k].newDelta = 0.0; 

                    if (deepestLevel < eAPT->sV[k].newLevel)
                        deepestLevel = eAPT->sV[k].newLevel;
                    appendDS(queue, levelIndices, eAPT->sV[k].newLevel, k);
                    if (eAPT->sV[k].touched == 0) {
                        printf("currRoot, currElement, k: %ld, %ld, %ld\n", currRoot, currElement, k);
                        printf("AAAAAAAAAAAAAAAAAAHHHHHHHHHHHHHHHHHHHHHHH\n");
                    }
                }
                if (eAPT->sV[k].newPathsToRoot == INFINITY_MY){
                    // k has not been found and therefore its paths to the roots are through its parent.
                    eAPT->sV[k].newPathsToRoot = eAPT->sV[currElement].newPathsToRoot;
                }
                else{
                    // k has been found and has multiple paths to the root as it has multiple parents.
                    eAPT->sV[k].newPathsToRoot += eAPT->sV[currElement].newPathsToRoot;
                    
                }    
            }
           
            if (eAPT->sV[k].touched != 0 && eAPT->sV[k].newLevel < eAPT->sV[currElement].newLevel) {
                eAPT->sV[currElement].newEdgesAbove += eAPT->sV[k].newEdgesAbove + 1;
            } else if (eAPT->sV[k].touched == 0 && tree->vArr[k].level < eAPT->sV[currElement].newLevel) {
                eAPT->sV[currElement].newEdgesAbove += eAPT->sV[k].newEdgesAbove + 1; 
            }

        }


        STINGER_FORALL_EDGES_OF_VTX_END();
#if COUNT_TRAVERSALS==1
        eAPT->dynamicTraverseEdgeCounter += stinger_typed_outdegree(sStinger, currElement, 0);
        eAPT->dynamicTraverseVerticeCounter++;
#endif
    }  
    // --------------- 

    /// Phase 5 ----------
    qEnd = 0;
    //for (int64_t lev = tree->vArr[startVertex].level; lev < NV; lev++) {

    for (int64_t lev = stopLevel; lev <= deepestLevel; lev++) {
    //for (int64_t lev = 0; lev < NV; lev++) {
        int64_t index = levelIndices[lev].front;
        int64_t levelEmpty = 1;
        while (index != -1) {
            levelEmpty = 0;
            queue_node_t *temp_node = queue->nodes + index;
            Queue[qEnd++] = temp_node->data;
            index = temp_node->next;
        }
        levelIndices[lev].front = -1;
        levelIndices[lev].back = -1;
        //if (levelEmpty) {
        //    break;
        //}
    }
    queue->size = 0; 
    // ----------- 

    // Phase 6 -------------
    #if 0
    for (int64_t k = 0; k < tvDownEnd; k++) {
        int64_t vertex = touchedVerticesDown[k];
        tree->vArr[vertex].level = eAPT->sV[vertex].newLevel;
    }
    #endif
    // ----------- 

    int64_t qUpStart = 0, qUpEnd = 0;
    qEnd--;

    if (tree->vArr[startVertex].level == INFINITY_MY) {
        if (eAPT->sV[parentVertex].touched == 0) {
            eAPT->sV[parentVertex].newDelta = tree->vArr[parentVertex].delta -
                    ((bc_t)tree->vArr[parentVertex].pathsToRoot / (bc_t)tree->vArr[startVertex].pathsToRoot) *
                    (bc_t)(tree->vArr[startVertex].delta + 1);

            eAPT->sV[parentVertex].newPathsToRoot = tree->vArr[parentVertex].pathsToRoot;
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

    // Phase 7 -------------- 

    while (1) {
        // Removing last element from the queue
        int64_t currElement = -1;
        if (qEnd < 0 && qUpStart >= qUpEnd) {
            break;
        } else if (qEnd < 0) {
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
                if (eAPT->sV[Queue[qEnd]].newLevel == 0) {
                    qEndLevel = tree->vArr[Queue[qEnd]].level;
                } else {
                    qEndLevel = eAPT->sV[Queue[qEnd]].newLevel;
                }

                //if (tree->vArr[QueueUp[qUpStart]].level < tree->vArr[Queue[qEnd]].level) {
                if (qUpStartLevel < qEndLevel) {
                    currElement = Queue[qEnd--];
                } else {
                    currElement = QueueUp[qUpStart++];
                }
            } else {
                currElement = Queue[qEnd--];
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
                if (eAPT->sV[k].touched == 0)
                    eAPT->sV[k].newEdgesBelow = tree->vArr[k].edgesBelow;
                eAPT->sV[currElement].newEdgesBelow += eAPT->sV[k].newEdgesBelow + 1;
            }

            //if (tree->vArr[k].level == tree->vArr[parentVertex].level) {
            if (kLevel == parentVertexLevel) {
                if (eAPT->sV[parentVertex].touched == 0) {
                    eAPT->sV[parentVertex].newDelta = tree->vArr[parentVertex].delta -
                            ((bc_t)tree->vArr[parentVertex].pathsToRoot / (bc_t)tree->vArr[startVertex].pathsToRoot) *
                            (bc_t)(tree->vArr[startVertex].delta + 1);

                    eAPT->sV[parentVertex].newPathsToRoot = tree->vArr[parentVertex].pathsToRoot;
                    eAPT->sV[parentVertex].touched = -1;
                    QueueUp[qUpEnd++] = parentVertex;
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
                    QueueUp[qUpEnd++] = k;
                }                 
                
                eAPT->sV[k].newDelta +=
                    ((bc_t)eAPT->sV[k].newPathsToRoot/(bc_t)eAPT->sV[currElement].newPathsToRoot)*
                    (bc_t)(eAPT->sV[currElement].newDelta+1);


                // For the elements that are touched in the ascent stage it is necessary to
                // to reduce the values that they previously had.
                // In addition to this, the "parentVertex" that is connected to "vertex", i.e.
                // the vertices of the new edge, needs to increase its betweenness centrality
                // following the new connection, without removing the old delta value.
                if(eAPT->sV[k].touched<0 && tree->vArr[k].level < tree->vArr[currElement].level) // && ( k!=parentVertex  || currElement!=startVertex))
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
        qStart = 0;
        qEnd = 1;
        Queue[0] = startVertex;
        
        eAPT->sV[startVertex].touched = -2;

        while (qStart != qEnd){
            int64_t currElement = Queue[qStart++];
            
            eAPT->sV[currElement].totalBC -= tree->vArr[currElement].delta;

            tree->vArr[currElement].edgesBelow = 0;
            tree->vArr[currElement].edgesAbove = 0;
            eAPT->sV[currElement].newEdgesAbove = 0;
            eAPT->sV[currElement].newEdgesBelow = 0;
            STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement)
            {
                uint64_t k = STINGER_EDGE_DEST;

                if (eAPT->sV[k].touched != -2) {
                    eAPT->sV[k].touched = -2;
                    touchedVerticesUp[tvUpEnd++] = k;
                    Queue[qEnd++] = k;
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
    // ----- 
}

