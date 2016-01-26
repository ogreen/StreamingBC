#pragma once


#include "stinger.h"
#include "streamingbc_aux.h"


#include "dsUtils.h"

StreamingExtraInfo insertEdgeBrandes(bcForest* forest, struct stinger* sStinger,
                       uint64_t newU, uint64_t newV, uint64_t * rootArrayForApproximation,int64_t NK, int64_t NV, int64_t NT, extraArraysPerThread** eAPT);

StreamingExtraInfo deleteEdgeBrandes(bcForest *forest, struct stinger *sStinger,
                       uint64_t oldU, uint64_t oldV, uint64_t *rootArrayForApproximation, int64_t NK, int64_t Nv, int64_t NT, extraArraysPerThread **eAPT);

