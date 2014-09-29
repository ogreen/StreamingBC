#pragma once

#include "stinger.h"
#include "bcTreeDS.h"
#include "list.h"

#include "dsUtils.h"


void bfsBrandes(bcForest* forest, struct stinger* sStinger,extraArraysPerThread* eAPT);

void bfsBrandesForApproxCase(bcForest* forest, struct stinger* sStinger, uint64_t* rootArrayForApproximation, 
		uint64_t rootArraySizeForApproximation,extraArraysPerThread* eAPT);
 void bfsBrandesForApproxCaseParallel(bcForest* forest, struct stinger* sStinger, uint64_t* rootArrayForApproximation,
		uint64_t rootArraySizeForApproximation,extraArraysPerThread** eAPT,int64_t NT);

uint64_t bfsBrandesPerTree(bcForest* forest, struct stinger* sStinger,
                            uint64_t currRoot, bc_t* totalBC,extraArraysPerThread* eAPT);


void addEdgeWithoutMovementBrandes(bcForest* forest, struct stinger* sStinger,
                            uint64_t currRoot, uint64_t startVertex, uint64_t parentVertex,
                            uint64_t addedPathsToRoot,  extraArraysPerThread* eAPT);



void moveUpTreeBrandes(bcForest* forest, struct stinger* sStinger,
                            uint64_t currRoot, uint64_t startVertex, uint64_t parentVertex,
                            uint64_t prevDist,  extraArraysPerThread* eAPT);

void addEdgeConnectingComponents(bcForest* forest, struct stinger* sStinger,
                            uint64_t newU, uint64_t newV, bc_t* totalBC);


void addEdgeConnectingComponentsBrandes(bcForest* forest,  struct stinger* sStinger,
                            uint64_t lowerVertex, uint64_t lowerTreeSize, uint64_t* lowerTreeVertices,
                            uint64_t higherVertex, uint64_t higherTreeSize, uint64_t* higherTreeVertices, bc_t* totalBC);


void bfsBrandesOld(uint64_t NV, uint64_t * graph[], uint64_t * level[],
                    uint64_t * pathsToRoot[],bc_t* deltaB[],bc_t* BCB,struct stinger* sStinger);

uint64_t bfsBrandesPerTreeOld(uint64_t NV, uint64_t * graph[], uint64_t currRoot, uint64_t * level[],
                    uint64_t * pathsToRoot[],bc_t * deltaB[],bc_t* BCB,struct stinger* sStinger);


void moveUpTreeBrandesOld(uint64_t * levels, uint64_t * pathsToRoot, uint64_t ** graph,
	int64_t vertex, int64_t parentVertex, uint64_t dist, uint64_t NV,
	bc_t* deltaB, bc_t* BCB,uint64_t root,
	uint64_t** queueBFSTREE, int64_t* levelCounter,bc_t* newdeltaB,uint64_t* oldPathsToRoot,
	struct stinger* sStinger);



void addEdgeWithoutMovementBrandesOld(uint64_t * level, uint64_t * pathsToRoot, uint64_t ** graph,
	uint64_t vertex, uint64_t deltaTop, uint64_t NV, bc_t* deltaB, bc_t* BCB,uint64_t root,
	uint64_t** queueBFSTREE, int64_t* levelCounter,bc_t* newdeltaB,uint64_t* oldPathsToRoot,uint64_t parentVertex,
	struct stinger* sStinger);


void insertEdgeBrandesOld(uint64_t u, uint64_t v, uint64_t NV, uint64_t * graph[], uint64_t * level[],
                       uint64_t * pathsToRoot[],bc_t* delta[],bc_t* BC,bc_t* newDelta[],uint64_t * oldPathsToRoot[],struct stinger* sStinger);




void SanityCheck(uint64_t NV,uint64_t * graph[],bc_t* BCSanity);


#define FOR_ADJ 0
#define FOR_STINGER 1
#define FOR_CSR 2
#define FOR_TYPE FOR_STINGER
