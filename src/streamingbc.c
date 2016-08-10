#include <omp.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <unistd.h>

#include "stinger.h"
#include "streamingbc.h"
#include "streamingbc_aux.h"
#include "timer.h"

bcForest* streamingBCCreateForestExact(uint64_t NV){
    return CreateForestExact(NV);
}

bcForest* streamingBCCreateForestApproximate(uint64_t NV, uint64_t NK, uint64_t *rootArray){
    return CreateForestApproximate(NV, rootArray, NK);
}

extraArraysPerThread** streamingBCCreateAuxilary(int64_t threadCount,int64_t NV){
    return createExtraArraysForThreads(threadCount,NV);
}

void streamingBCInitStaticExact(bcForest* forest, struct stinger*stingerGraph, uint64_t NT,extraArraysPerThread** auxilary){
    BrandesExactParallel(forest,stingerGraph,NT, auxilary);
}

void streamingBCInitStaticApproximate(bcForest* forest, struct stinger*stingerGraph, uint64_t NT,extraArraysPerThread** auxilary,uint64_t NK, uint64_t *rootArray){
    BrandesApproxCaseParallel(forest, stingerGraph, rootArray, NK,auxilary,NT);
}

void streamingBCDeleteAuxilary(extraArraysPerThread** parallelExtra, int64_t threadCount, int64_t NV){
    destroyExtraArraysForThreads(parallelExtra,threadCount,NV);
}


void streamingBCDeleteForestExact(bcForestPtr* deadForest){
    DestroyForestExact(deadForest);
}
void streamingBCDeleteForestApproximate(bcForestPtr* deadForest, uint64_t rootArraySize,uint64_t* rootArray){
    DestroyForestApproximate(deadForest,rootArray, rootArraySize);
}


StreamingExtraInfo insertVertexStreamingBC(bcForest* forest, struct stinger* sStinger, uint64_t src,
                       uint64_t* adjacencyArray,uint64_t adjacencySize, uint64_t * rootArrayForApproximation, 
                       int64_t NK, int64_t NV, int64_t NT, extraArraysPerThread** eAPT){
    StreamingExtraInfo oneSEI,returnsei;
    returnsei.adjacent=0; returnsei.movement=0; returnsei.sameLevel=0;
    for(uint64_t d=0; d<adjacencySize; d++){
        uint64_t dest = adjacencyArray[d];
        stinger_insert_edge(sStinger,0,src,dest,0,0);
        stinger_insert_edge(sStinger,0,dest,src,0,0);
        oneSEI=insertEdgeStreamingBC(forest,sStinger,src,dest,rootArrayForApproximation,NK, NV,NT,eAPT, 1);
        returnsei.adjacent += oneSEI.adjacent; returnsei.movement += oneSEI.movement; returnsei.sameLevel += oneSEI.sameLevel;
    }
    return returnsei;    
}


int64_t compareArrays(const void *arr1, const void *arr2) {
    const int64_t* one = (const int64_t *) arr1;
    const int64_t* two = (const int64_t *) arr2;

    return two[1] - one[1];
}

StreamingExtraInfo insertEdgeStreamingBC(bcForest* forest, struct stinger* sStinger,
        uint64_t newU, uint64_t newV, uint64_t * rootArrayForApproximation,int64_t NK, int64_t NV, int64_t NT,
        extraArraysPerThread** eAPT, uint32_t loadBalancing){
    omp_set_num_threads(NT);

    int64_t workPerVertex[NK][2]; // First column has vertex ids, second col has work values per id.
    uint64_t currRoot = 0;
    uint64_t samelevel = 0, compConn = 0, adjacent=0, movement=0;
    int64_t workIndex = 0;

    for(currRoot = 0; currRoot < NK; currRoot++){
        uint64_t i = rootArrayForApproximation[currRoot];
        int64_t thread = 0;
        extraArraysPerThread* myExtraArrays = eAPT[thread];
        bcTree* tree = forest->forest[i];

        if (loadBalancing == 1)
            workPerVertex[workIndex][0] = i;

        int64_t diff = tree->vArr[newU].level - tree->vArr[newV].level; 

        if (loadBalancing == 1) {
            if (diff < 0) {
                workPerVertex[workIndex++][1] = 2 * tree->vArr[newV].edgesBelow + tree->vArr[newV].edgesAbove;
            } else if (diff > 0) {
                workPerVertex[workIndex++][1] = 2 * tree->vArr[newU].edgesBelow + tree->vArr[newU].edgesAbove;
            } else {
                workPerVertex[workIndex++][1] = 0;
            }
        }
    }

    if (loadBalancing == 1)
        qsort((const int*)&workPerVertex, workIndex, sizeof(int64_t[2]), compareArrays);

    uint64_t r;

    #pragma omp parallel for schedule(dynamic,1)
    for(r = 0; r < NK; r++)
    {
        int64_t i = workPerVertex[r][0];
        if (loadBalancing == 0) {
            i = rootArrayForApproximation[r];
        }

        int64_t thread = omp_get_thread_num();
        extraArraysPerThread* myExtraArrays = eAPT[thread];

        bcTree* tree = forest->forest[i];
        int64_t diff = tree->vArr[newU].level - tree->vArr[newV].level;

        if(diff < -1 || diff > 1) {

            if(diff<-1){
                moveUpTreeBrandes(forest, sStinger, i, newV, newU, (-diff) - 1,  myExtraArrays, (uint64_t)1);
            }
            else{
                moveUpTreeBrandes(forest,  sStinger, i, newU, newV, (diff) - 1, myExtraArrays, (uint64_t)1);
            }
            //printf("%.9lf\n", (double)caseTime); fflush(stdout); 
            eAPT[thread]->movementCounter++;
        }
        // Newly inserted edge is connecting vertices that were in adjacent levels before insertions
        else if(diff == -1 || diff == 1){

            //printf("Case II root: %ld\n", i); fflush(stdout);
            if(diff==-1) {
                addEdgeWithoutMovementBrandes(forest, sStinger, i, newV, newU, tree->vArr[newU].pathsToRoot,myExtraArrays, (uint64_t)1);
            }
            else{
                addEdgeWithoutMovementBrandes(forest, sStinger, i, newU, newV, tree->vArr[newV].pathsToRoot, myExtraArrays, (uint64_t)1);
            }
            //printf("%.9lf\n", (double) caseTime); fflush(stdout);
            eAPT[thread]->adjacentCounter++;
        } 

    }
     
    #pragma omp parallel for
    for(uint64_t v=0;v<NV;v++){
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


StreamingExtraInfo deleteVertexStreamingBC(bcForest* forest, struct stinger* sStinger, uint64_t src,
                       uint64_t* adjacencyArray,uint64_t* adjacencySize, uint64_t * rootArrayForApproximation, 
                       int64_t NK, int64_t NV, int64_t NT, extraArraysPerThread** eAPT){
    StreamingExtraInfo oneSEI,returnsei;
    returnsei.adjacent=0; returnsei.movement=0; returnsei.sameLevel=0;

    int64_t d=0;
    STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,src)
    {
        uint64_t dest = STINGER_EDGE_DEST;
        adjacencyArray[d]=dest;
        stinger_remove_edge(sStinger, 0, src, dest);
        stinger_remove_edge(sStinger, 0, dest, src);

        oneSEI=deleteEdgeStreamingBC(forest,sStinger,src,dest,rootArrayForApproximation,NK, NV,NT,eAPT, 1);
        returnsei.adjacent += oneSEI.adjacent; returnsei.movement += oneSEI.movement; returnsei.sameLevel += oneSEI.sameLevel;
    }
    STINGER_FORALL_EDGES_OF_VTX_END();        
    *adjacencySize=d;
    return returnsei;
}

StreamingExtraInfo deleteEdgeStreamingBC(bcForest *forest, struct stinger *sStinger, uint64_t oldU, uint64_t oldV,
                            uint64_t *rootArrayForApproximation, int64_t NK, int64_t NV, int64_t NT,
                            extraArraysPerThread **eAPT, uint32_t loadBalancing) {
    omp_set_num_threads(NT);

    uint64_t currRoot = 0;
    uint64_t samelevel = 0, compConn = 0, adjacent = 0, movement = 0;
    
    uint64_t r;
    int64_t thread = 0;

    int64_t workPerVertex[NK][2];
    int64_t workIndex = 0;

    if (loadBalancing == 1) {
        for (r = 0; r < NK; r++) {
            int64_t i = rootArrayForApproximation[r];
            bcTree *tree = forest->forest[i];
            int64_t diff = tree->vArr[oldU].level - tree->vArr[oldV].level;
            workPerVertex[workIndex][0] = i;
            if (diff < 0) {
                workPerVertex[workIndex++][1] = 2 * tree->vArr[oldV].edgesBelow + tree->vArr[oldV].edgesAbove;
            } else if (diff > 0) {
                workPerVertex[workIndex++][1] = 2 * tree->vArr[oldU].edgesBelow + tree->vArr[oldU].edgesAbove;
            } else {
                workPerVertex[workIndex++][1] = 0;
            }
        }
    }
    
    if (loadBalancing == 1)
        qsort((const int*)&workPerVertex, workIndex, sizeof(int64_t[2]), compareArrays);
    
    #pragma omp parallel for schedule(dynamic,1)
    for (r = 0; r < NK; r++){ 
        int64_t i = workPerVertex[r][0];
        if (loadBalancing == 0)
            i = rootArrayForApproximation[r];

        int64_t thread = omp_get_thread_num();
        extraArraysPerThread *myExtraArrays = eAPT[thread];

        bcTree *tree = forest->forest[i];

        int64_t extraParents = 0;
        int64_t childVertex = oldU;
        int64_t parentVertex = oldV;
        
        int64_t diff = tree->vArr[oldU].level - tree->vArr[oldV].level;
        
        if (diff == 0){ 
            //printf("%.9lf\n", 0.0f);
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
            removeEdgeWithoutMovementBrandes(forest, sStinger, i, childVertex, parentVertex, 
                        tree->vArr[parentVertex].pathsToRoot, myExtraArrays, (uint64_t)1);
            //printf("%.9lf\n", (double) caseTime); fflush(stdout);
            eAPT[thread]->adjacentCounter++;
            adjacent++;
        }
        else{ 
            moveDownTreeBrandes(forest, sStinger, i, childVertex, parentVertex, myExtraArrays, (uint64_t)2);
            //printf("%.9lf\n",  (double) caseTime); fflush(stdout);
            eAPT[thread]->movementCounter++;
            movement++;
        } 
    }  
  
    int64_t tlow = (NV * thread) / NT;
    int64_t thigh = (NV * (thread + 1)) / NT ;
       
    for (uint64_t v = tlow; v < NV; v++)
    { 
        for (uint64_t t = 0; t < NT; t++)
        {
            forest->totalBC[v] += eAPT[t]->sV[v].totalBC;
            eAPT[t]->sV[v].totalBC = 0.0;
        }
    }
    
    StreamingExtraInfo returnSEI = {0,0,0,0};

    returnSEI.sameLevel = samelevel;
    returnSEI.adjacent = adjacent;
    returnSEI.movement = movement;
    
    return returnSEI;
}

