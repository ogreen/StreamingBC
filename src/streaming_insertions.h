#pragma once

#include "streamingbc_aux.h"
#include "stinger.h"

//----------------------------------------------
//----------------------------------------------
// Static Graph BC


void bfsBrandes(bcForest* forest, struct stinger* sStinger,extraArraysPerThread* eAPT);

void bfsBrandesForApproxCase(bcForest* forest, struct stinger* sStinger, uint64_t* rootArrayForApproximation, 
		uint64_t rootArraySizeForApproximation,extraArraysPerThread* eAPT);
void bfsBrandesForApproxCaseParallel(bcForest* forest, struct stinger* sStinger, uint64_t* rootArrayForApproximation,
		uint64_t rootArraySizeForApproximation,extraArraysPerThread** eAPT,int64_t NT);

uint64_t bfsBrandesPerTree(bcForest* forest, struct stinger* sStinger,
                            uint64_t currRoot, bc_t* totalBC,extraArraysPerThread* eAPT);

//----------------------------------------------
//----------------------------------------------
// Dynamic Graph Insertion


void addEdgeWithoutMovementBrandes(bcForest* forest, struct stinger* sStinger,
                            uint64_t currRoot, uint64_t startVertex, uint64_t parentVertex,
                            uint64_t addedPathsToRoot,  extraArraysPerThread* eAPT);

void moveUpTreeBrandes(bcForest* forest, struct stinger* sStinger,
                            uint64_t currRoot, uint64_t startVertex, uint64_t parentVertex,
                            uint64_t prevDist,  extraArraysPerThread* eAPT);


//----------------------------------------------
//----------------------------------------------
// Dynamic Graph Deletion

void removeEdgeWithoutMovementBrandes(bcForest* forest, struct stinger* sStinger,
                            uint64_t currRoot, uint64_t startVertex, uint64_t parentVertex,
                            uint64_t deletedPathsFromRoot, bc_t* totalBC);


void moveDownTreeBrandes(bcForest* forest, struct stinger* sStinger,
                            uint64_t currRoot, uint64_t startVertex, uint64_t parentVertex, bc_t* totalBC);





