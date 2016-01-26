#pragma once


#include "stinger.h"
#include "bcTreeDS.h"
#include "streamingbc_aux.h"
#include "csr.h"

uint64_t** createParentArray(csrGraph* graph,uint64_t NV);
void destroyParentArray(uint64_t** parentArray,uint64_t NV);

uint64_t*** createParallelParentArray(csrGraph* graph,uint64_t NV, uint64_t threadCount);
void destroyParallelParentArray(uint64_t*** parallelParentArray,uint64_t NV,uint64_t threadCount);

uint64_t** createParentArrayStinger(struct stinger* s,uint64_t NV);
uint64_t*** createParallelParentArrayStinger(struct stinger* S,uint64_t NV, uint64_t threadCount);


uint64_t** createMultiLevelQueue(uint64_t NV);
void destroyMultiLevelQueue(uint64_t** multiLevelQueue,uint64_t NV);
uint64_t*** createParallelMultiLevelQueue(uint64_t NV, uint64_t threadCount);
void destroyParallelMultiLevelQueue(uint64_t*** parallelMultiLevelQueue,uint64_t NV,uint64_t threadCount);

bcTree** createParallelForest(int64_t threadCount,int64_t NV);
void destroyParallelForest(bcTree** parallelForest, int64_t threadCount);


list_ptr** createParallelList(int64_t threadCount,int64_t NV);
void destroyParallelList(list_ptr** parallelList, int64_t threadCount,int64_t NV);

float** createParallelBetweennessArray(int64_t threadCount,int64_t NV);
void destroyParallelBetweennessArray(float** parallelList, int64_t threadCount);

 typedef struct{

    int64_t diffPath;
    int64_t touched;
    int64_t newPathsToRoot;
    int64_t movementDelta;
    bc_t newDelta;
    bc_t totalBC;

    int64_t IMoved;
}sbcV;


typedef struct {

    sbcV* sV;

    int64_t* QueueDown;
    int64_t* QueueUp;
    int64_t* QueueSame;
    int64_t* Stack;

	uint64_t samelevelCounter;
	uint64_t compConnCounter;
	uint64_t adjacentCounter;
	uint64_t movementCounter;

    uint64_t staticTraverseVerticeCounter;
    uint64_t staticTraverseEdgeCounter;

    uint64_t dynamicTraverseVerticeCounter;
    uint64_t dynamicTraverseEdgeCounter;


    list_ptr* multiLevelQueues;

	uint64_t dummy[8];

} extraArraysPerThread;

extraArraysPerThread* createExtraArraysPerThread(int64_t NV);
void ClearCounters(extraArraysPerThread* eapt);
void destroyExtraArraysPerThread(extraArraysPerThread* eapt,int64_t NV);
extraArraysPerThread** createExtraArraysForThreads(int64_t threadCount,int64_t NV);
void destroyExtraArraysForThreads(extraArraysPerThread** parallelExtra, int64_t threadCount, int64_t NV);


#define COUNT_TRAVERSALS 1
