#pragma once

#include "stinger.h"
#include "bcTreeDS.h"

#include "dsUtils.h"

StreamingExtraInfo insertEdgeBrandes(bcForest* forest, struct stinger* sStinger,
                       uint64_t newU, uint64_t newV, uint64_t * rootArrayForApproximation,int64_t NK,extraArraysPerThread** eAPT);

//StreamingExtraInfo removeEdgeBrandes(bcForest* forest, uint64_t * matGraph[], struct stinger* sStinger,
//                       uint64_t newU, uint64_t newV,uint64_t ** queueBFSTREE, uint64_t * rootArrayForApproximation);
