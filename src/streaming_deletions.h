#pragma once

#include "streamingbc_aux.h"
#include "stinger.h"
#include "bcTreeDS.h"

#include "dsUtils.h"


void removeEdgeWithoutMovementBrandes(bcForest* forest, struct stinger* sStinger,
                            uint64_t currRoot, uint64_t startVertex, uint64_t parentVertex,
                            uint64_t deletedPathsFromRoot, bc_t* totalBC);


void moveDownTreeBrandes(bcForest* forest, struct stinger* sStinger,
                            uint64_t currRoot, uint64_t startVertex, uint64_t parentVertex, bc_t* totalBC);

