

#include <alloca.h>
#include <stdio.h>
#include <stdlib.h>
#include "xmalloc.h"

#include "bcTreeDS.h"


/**
* @brief Creates the data structures needed for computing BC.
*       These include the level of each vertex in the BFS tree.
*       The number of shortest paths each vertex has to the root.
*       The delta value computed in the dependency accumulation stage/
*
* @param numVertices The number of vertices in the graph
*
* @return Returns the created data structure.
*/
bcTree* CreateTree(int64_t numVertices)
{
    bcTree* newTree;
    newTree = (bcTree*)xmalloc(sizeof(bcTree));

    newTree->NV = numVertices;
    newTree->vArr = (bcV*)xmalloc(numVertices*sizeof(bcV));


    return newTree;
}

/**
* @brief Destroys the allocated tree.
*
* @param deadTree The tree that needs to be unallocated.
*
* @return None
*/
void DestroyTree(bcTreePtr deadTree)
{
    free(deadTree->vArr);

    free(deadTree);
}

/**
* @brief Creates a tree for each vertex in the rootArray. Thus, a total of rootArraySize trees are created.
*
* @param newForest This is an OUT parameter that contains the forest for approximate case.
* @param numVertices The number of vertices in the graph
* @param rootArray The array containing vertex id's to be considered for approximate case
* @param rootArraySize The number of vertices in the graph to be considered for approximate case
*
* @return None (see parameter newForest)
*/
bcForest* CreateForestForApproxCase(bcForestPtr* newForest, int64_t numVertices, uint64_t* rootArray, uint64_t rootArraySize)
{
	bcForest* temp;//=*newForest;
    temp = (bcForest*)xmalloc(sizeof(bcForest));

    temp->NV = numVertices;
    temp->forest = (bcTreePtr*)xmalloc(numVertices*sizeof(bcTreePtr));
    temp->totalBC = (bc_t*)xmalloc(numVertices*sizeof(bc_t));

    int64_t i;
    for(i=0;i<rootArraySize;i++)
    {
      temp->forest[rootArray[i]] = CreateTree(numVertices);
    }
	return temp;

}



/**
* @brief Destroys the allocated forest for approximate case
*
* @param deadForest The forest that needs to be unallocated.
* @param rootArray The array containing vertex id's to be considered for approximate case
* @param rootArraySize The number of vertices in the graph to be considered for approximate case
*
* @return None
*/
void DestroyForestForApproxCase(bcForestPtr* deadForest, uint64_t* rootArray, uint64_t rootArraySize)
{

	bcForest* temp=*deadForest;
    int64_t i;
    for(i=0;i<(rootArraySize);i++)
    {
        DestroyTree(temp->forest[rootArray[i]]);
    }
   free(temp->totalBC);

   free(temp->forest);

    free(temp);

    return;
}

