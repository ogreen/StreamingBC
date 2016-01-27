#include "omp.h"
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "operation.h"
#include "streaming_deletions.h"
#include "streaming_insertions.h"
#include "xmalloc.h"

#include "streamingbc_aux.h"

#include "dsUtils.h"
#include <math.h>
//#define IS_PARALLEL
#include "timer.h"


StreamingExtraInfo insertEdgeBrandes(bcForest* forest, struct stinger* sStinger,
        uint64_t newU, uint64_t newV, uint64_t * rootArrayForApproximation,int64_t NK, int64_t NV, int64_t NT,
        extraArraysPerThread** eAPT){
    omp_set_num_threads(NT);
    //printf("newU: %ld\n", newU);
    //printf("newV: %ld\n", newV);

    int64_t adjRootArray[NK];
    int64_t moveRootArray[NK];

    for(int64_t tk=0; tk<NK; tk++) {adjRootArray[tk]=0; moveRootArray[tk]=0;}

    uint64_t currRoot = 0;
    uint64_t samelevel = 0, compConn = 0, adjacent=0, movement=0;

    for(currRoot = 0; currRoot < NK; currRoot++){
        uint64_t i = rootArrayForApproximation[currRoot];
        int64_t thread = 0;
        extraArraysPerThread* myExtraArrays = eAPT[thread];
        bcTree* tree = forest->forest[i];
        int64_t diff = tree->vArr[newU].level - tree->vArr[newV].level;
        // New edge is connecting vertices in the same level. Nothing needs to be done.
        if(diff==0){
            eAPT[thread]->samelevelCounter++;
            continue;
        }
        // Newly inserted edge is causing movement in the tee
        if(diff < -1 || diff > 1){
            moveRootArray[ eAPT[thread]->movementCounter++]=i;
        }
        // Newly inserted edge is connecting vertices that were in adjacent levels before insertions
        else if(diff == -1 || diff == 1){
            adjRootArray[ eAPT[thread]->adjacentCounter++]=i;
        }
    }

    for(uint64_t thread=0; thread<1; ++thread){
        samelevel += eAPT[thread]->samelevelCounter;
        compConn += eAPT[thread]->compConnCounter;
        adjacent += eAPT[thread]->adjacentCounter;
        movement += eAPT[thread]->movementCounter;

        eAPT[thread]->compConnCounter=0;
        eAPT[thread]->adjacentCounter=0;
        eAPT[thread]->movementCounter=0;
    }

    int64_t newRootArray[NK];
    int64_t counter=0, thread=0;
    for(int64_t m=0; m<movement; m++){
        newRootArray[counter++]=moveRootArray[m];
    }
    for(int64_t a=0; a<adjacent; a++){
        newRootArray[counter++]=adjRootArray[a];
    }

    double times[NT];
    uint64_t r;

    int64_t edgeCountValues[NV];

    //for (uint64_t k = 0; k < NV; k++)
    //    edgeCountValues[k] = -1;

    #pragma omp parallel for schedule(dynamic,1)
//    for(r = 0; r < counter; r++)
    for(r = 0; r < NK; r++)
    {
//        int64_t i=newRootArray[r];
        int64_t i=rootArrayForApproximation[r];
        int64_t thread = omp_get_thread_num();
        extraArraysPerThread* myExtraArrays = eAPT[thread];

        bcTree* tree = forest->forest[i];
        int64_t diff = tree->vArr[newU].level - tree->vArr[newV].level;
        // New edge is connecting vertices in the same level. Nothing needs to be done.
        //              if(diff==0)
        //              {
        //                  eAPT[thread]->sameevelCounter++;
        //                    continue;
        //              }
        // Newly inserted edge is causing movement in the tee
        if(diff < -1 || diff > 1){
            uint64_t prevEdgeCount   = myExtraArrays->dynamicTraverseEdgeCounter;
            uint64_t prevVertexCount = myExtraArrays->dynamicTraverseVerticeCounter;

            if(diff<-1){
                moveUpTreeBrandes(forest, sStinger, i, newV, newU, (-diff) - 1,  myExtraArrays);
            }
            else{
                moveUpTreeBrandes(forest,  sStinger, i, newU, newV, (diff) - 1, myExtraArrays);
            }

            uint64_t edgeCount   = myExtraArrays->dynamicTraverseEdgeCounter - prevEdgeCount;
            uint64_t vertexCount = myExtraArrays->dynamicTraverseVerticeCounter - prevVertexCount;            
            //edgeCountValues[i] = edgeCount;
            //printf("%ld\n", edgeCount);
            //printf("Case III: edgeCount: %lf\n", ((double)edgeCount / (double)totalEdges));
            eAPT[thread]->movementCounter++;
        }
        // Newly inserted edge is connecting vertices that were in adjacent levels before insertions
        else if(diff == -1 || diff == 1){
            uint64_t prevEdgeCount   = myExtraArrays->dynamicTraverseEdgeCounter;
            uint64_t prevVertexCount = myExtraArrays->dynamicTraverseVerticeCounter;

            if(diff==-1){
                addEdgeWithoutMovementBrandes(forest, sStinger, i, newV, newU, tree->vArr[newU].pathsToRoot,myExtraArrays);
            }
            else{
                addEdgeWithoutMovementBrandes(forest, sStinger, i, newU, newV, tree->vArr[newV].pathsToRoot, myExtraArrays);
            }
             
            uint64_t edgeCount   = myExtraArrays->dynamicTraverseEdgeCounter - prevEdgeCount;
            uint64_t vertexCount = myExtraArrays->dynamicTraverseVerticeCounter - prevVertexCount;
                
            //edgeCountValues[i] = edgeCount;
            //printf("Case II: edgeCount: %lf\n", ((double)edgeCount / (double)totalEdges));
            eAPT[thread]->adjacentCounter++;
        }
//        #pragma omp barriar

    }
    
    /*for (uint64_t u = 0; u < NV; u++)
    {
        int64_t neighborEdgeCount = 0;
        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, u)
        {
            uint64_t v = STINGER_EDGE_DEST;
            
            if (edgeCountValues[v] != -1)   
                neighborEdgeCount += edgeCountValues[v];

        }
        STINGER_FORALL_EDGES_OF_VTX_END();
       
        double averageNeighborEdgeCount = 0.0;
        
        if (stinger_typed_outdegree(sStinger, u, 0) >= 1)
            averageNeighborEdgeCount = neighborEdgeCount / stinger_typed_outdegree(sStinger, u, 0);
        
        if (edgeCountValues[u] != -1 && edgeCountValues[u] < 125000)
            //printf("%ld\n", averageNeighborEdgeCount);
            printf("%ld,%lf\n", edgeCountValues[u], averageNeighborEdgeCount);
            //printf("edgeCount[u], averageNeighborEdgeCount[u]: %ld, %lf\n", edgeCountValues[u], averageNeighborEdgeCount);
    }*/
 
    int64_t tlow=(NV*thread)/NT;
    int64_t thigh=(NV*(thread+1))/NT;
    for(uint64_t v=tlow;v<thigh;v++){
        for(uint64_t t=0;t<NT;t++){
            forest->totalBC[v]+=eAPT[t]->sV[v].totalBC;
            eAPT[t]->sV[v].totalBC = 0.0;
        }
    }

    StreamingExtraInfo returnSEI={0,0,0,0};
    returnSEI.sameLevel= samelevel;
    returnSEI.adjacent = adjacent;
    returnSEI.movement = movement;

return returnSEI;
}

StreamingExtraInfo deleteEdgeBrandes(bcForest *forest, struct stinger *sStinger, uint64_t oldU, uint64_t oldV,
                            uint64_t *rootArrayForApproximation, int64_t NK, int64_t NV, int64_t NT,
                            extraArraysPerThread **eAPT){
    omp_set_num_threads(NT);

    int64_t adjRootArray[NK];
    int64_t moveRootArray[NK];

    for (int64_t i = 0; i < NK; i++){
        adjRootArray[i] = 0;
        moveRootArray[i] = 0;
    }

    //printf("oldU: %ld\n", oldU);
    //printf("oldV: %ld\n", oldV);

    uint64_t currRoot = 0;
    uint64_t samelevel = 0, compConn = 0, adjacent = 0, movement = 0;
    
    uint64_t r;
    int64_t thread = 0;

    int64_t edgeCountValues[NV];
    
    //for (uint64_t k = 0; k < NV; k++)
    //    edgeCountValues[k] = -1;

    #pragma omp parallel for schedule(dynamic,1)
    for (r = 0; r < NK; r++){ 
        int64_t i = rootArrayForApproximation[r];
        int64_t thread = omp_get_thread_num();
        extraArraysPerThread *myExtraArrays = eAPT[thread];

        bcTree *tree = forest->forest[i];

        int64_t extraParents = 0;
        int64_t childVertex = oldU;
        int64_t parentVertex = oldV;
        
        int64_t diff = tree->vArr[oldU].level - tree->vArr[oldV].level;
        if (diff == 0){ 
            eAPT[thread]->samelevelCounter++;
            samelevel++;
            continue;
        }
        
        if (tree->vArr[oldU].level < tree->vArr[oldV].level){
            childVertex = oldV;
            parentVertex = oldU;
        }
        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, childVertex)
        {
            uint64_t neighbor = STINGER_EDGE_DEST;

            if (tree->vArr[neighbor].level + 1 == tree->vArr[childVertex].level)
            {
                extraParents++;
            }
        }
        STINGER_FORALL_EDGES_OF_VTX_END();
        
        if (extraParents >= 1){                 
            //uint64_t prevEdgeCount   = myExtraArrays->dynamicTraverseEdgeCounter;
            //uint64_t prevVertexCount = myExtraArrays->dynamicTraverseVerticeCounter;
            
            removeEdgeWithoutMovementBrandes(forest, sStinger, i, childVertex, parentVertex, 
                        tree->vArr[parentVertex].pathsToRoot, myExtraArrays);

            //uint64_t edgeCount = myExtraArrays->dynamicTraverseEdgeCounter - prevEdgeCount;
            //uint64_t vertexCount = myExtraArrays->dynamicTraverseVerticeCounter - prevVertexCount;


            //printf("Case II: root, edgeCount, vertexCount: %ld, %ld, %ld\n", i, edgeCount, vertexCount);
            //edgeCountValues[i] = edgeCount;
            eAPT[thread]->adjacentCounter++;
            adjacent++;
        }
        else{ 
            //uint64_t prevEdgeCount   = myExtraArrays->dynamicTraverseEdgeCounter;
            //uint64_t prevVertexCount = myExtraArrays->dynamicTraverseVerticeCounter;
            
            moveDownTreeBrandes(forest, sStinger, i, childVertex, parentVertex, myExtraArrays);
            
            //uint64_t edgeCount = myExtraArrays->dynamicTraverseEdgeCounter - prevEdgeCount;
            //uint64_t vertexCount = myExtraArrays->dynamicTraverseVerticeCounter - prevVertexCount;

            //printf("Case III: childVertex, parentVertex: %ld, %ld\n", childVertex, parentVertex);
            //printf("Case III: root, edgeCount, vertexCount: %ld, %ld, %ld\n", i, edgeCount, vertexCount);
            //edgeCountValues[i] = edgeCount; 
            eAPT[thread]->movementCounter++;
            movement++;
        }
    }
    
    /* Edge counting.
    for (uint64_t u = 0; u < NV; u++)
    {
        int64_t neighborEdgeCount = 0;
        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, u)
        {
            uint64_t v = STINGER_EDGE_DEST;
            
            if (edgeCountValues[v] != -1)   
                neighborEdgeCount += edgeCountValues[v];

        }
        STINGER_FORALL_EDGES_OF_VTX_END();
       
        double averageNeighborEdgeCount = 0.0;
        
        if (stinger_typed_outdegree(sStinger, u, 0) >= 1)
            averageNeighborEdgeCount = neighborEdgeCount / stinger_typed_outdegree(sStinger, u, 0);
        
        if (edgeCountValues[u] != -1)
            //printf("%ld\n", averageNeighborEdgeCount);
            printf("%ld,%lf\n", edgeCountValues[u], averageNeighborEdgeCount);
            //printf("edgeCount[u], averageNeighborEdgeCount[u]: %ld, %lf\n", edgeCountValues[u], averageNeighborEdgeCount);
    }*/
 
    
    int64_t tlow = (NV * thread) / NT;
    int64_t thigh = (NV * (thread + 1)) / NT ;
       
    for (uint64_t v = tlow; v < thigh; v++)
    { 
        for (uint64_t t = 0; t < NT; t++)
        {
            forest->totalBC[v] += eAPT[t]->sV[v].totalBC;
            eAPT[t]->sV[v].totalBC = 0.0;
        }
    } 
    
    //printf("movement: %ld\n", movement);    
    StreamingExtraInfo returnSEI = {0,0,0,0};

    returnSEI.sameLevel = samelevel;
    returnSEI.adjacent = adjacent;
    returnSEI.movement = movement;
    
    return returnSEI;
}
