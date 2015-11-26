#pragma once

#include "stinger.h"
#include "bcTreeDS.h"
#include "list.h"

#include "dsUtils.h"

void deleteEdgeWithoutMovement(bcForest *forest, struct stinger *sStinger, uint64_t currRoot, uint64_t startVertex,
                        uint64_t parentVertex, uint64_t deletedPathsToRoot, extraArraysPerThread *eAPT);

void moveDownTreeBrandes(bcForest *forest, struct stinger *sStinger, uint64_t currRoot, uint64_t startVertex,
                        uint64_t parentVertex, uint64_t prevDist, extraArraysPerThread *eAPT);
