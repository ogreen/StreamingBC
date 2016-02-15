#pragma once

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include "stinger.h"
#include "xmalloc.h"

#define INFINITY_MY LONG_MAX
#define COUNT_TRAVERSALS 1

typedef double bc_t;

//------------------------------------------------
//------------------------------------------------
// Betweennes centrality data structures


typedef struct{
    int64_t level;
    int64_t pathsToRoot;
    bc_t delta;
}bcV;

typedef struct {
    int64_t NV;
    bcV* vArr;
	double dummy[4];
} bcTree;
typedef bcTree* bcTreePtr;

typedef struct {
	bcTreePtr* forest;
    bc_t* totalBC;
    int64_t NV;
} bcForest;
typedef bcForest* bcForestPtr;

typedef struct
{
    uint64_t movement;
    uint64_t adjacent;
    uint64_t sameLevel;
    uint64_t connectedComponents;
}StreamingExtraInfo;

bcForest* CreateForestExact(int64_t numVertices);
void DestroyForestExact(bcForest** deadForest);

bcForest* CreateForestApproximate(int64_t numVertices, uint64_t* rootArray, uint64_t rootArraySize);
void DestroyForestForApproximate(bcForest** deadForest, uint64_t* rootArray, uint64_t rootArraySize);


//------------------------------------------------
//------------------------------------------------
// List

typedef struct struct_node{
    int64_t id;
    int64_t aux_val;
    struct struct_node *next;
} node_t;

typedef struct{
    node_t *head;
    node_t *tail;
    int64_t size;
} list_t;

typedef list_t* list_ptr;

list_t* makeList(void);
node_t* makeNode(int64_t);
void append(list_t*,node_t*);
node_t *getFirst(list_t*);
void deleteFirst(list_t*);
void printList(list_t*);

void emptyList(list_t*);

void makeArrayOfLists(list_ptr** aL,int64_t numberOfLists);
void destroyArrayOfLists(list_ptr** aL,int64_t numberOfLists);


//------------------------------------------------
//------------------------------------------------
// Auxilary data structures


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
    int64_t newLevel;
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


