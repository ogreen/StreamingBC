
#include "omp.h"
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#include "streamingbc_internal.h"
#include "streamingbc_aux.h"

void BrandesExact(bcForest * forest, struct stinger * sStinger, extraArraysPerThread * eAPT)
{
    for (int64_t i = 0; i < forest->NV; i++)
        forest->totalBC[i] = 0.0;

    for (int64_t i = 0; i < forest->NV; i++) {
        //		VB("\tTree[%ld]...\n",i)
        BrandesSingleTree(forest, sStinger, i, forest->totalBC, eAPT);
    }

    //	printf("bfsBrandes completed, forest=%ld\n",forest);
}


void BrandesExactParallel(bcForest * forest, struct stinger * sStinger, int64_t NT, extraArraysPerThread ** eAPT)
{
    for (int64_t i = 0; i < forest->NV; i++)
        forest->totalBC[i] = 0.0;

    #pragma omp parallel for schedule(dynamic,1)

    for (int64_t i = 0; i < forest->NV; i++) {
        int64_t thread = omp_get_thread_num();
        BrandesSingleTree(forest, sStinger, i, forest->totalBC, eAPT[thread]);
    }

    #pragma omp parallel for

    for (int64_t v = 0; v < forest->NV; v++) {
        for (int64_t t = 0; t < NT; t++) {
            forest->totalBC[v] += eAPT[t]->sV[v].totalBC;
            eAPT[t]->sV[v].totalBC = 0.0;
        }
    }
}

void BrandesApproxCase(bcForest * forest, struct stinger * sStinger, int64_t * rootArrayForApproximation,
                       int64_t rootArraySizeForApproximation, extraArraysPerThread * eAPT)
{
    //	printf("bfsBrandesForApproxCase started, forest=%ld\n",forest);
    //bfsBrandes(forest, sStinger);

    for (int64_t i = 0; i < forest->NV; i++)
        forest->totalBC[i] = 0.0;

    for (int64_t i = 0; i < rootArraySizeForApproximation; i++) {
        BrandesSingleTree(forest, sStinger, rootArrayForApproximation[i], forest->totalBC, eAPT);
    }
}

void BrandesApproxCaseParallel(bcForest * forest, struct stinger * sStinger, int64_t * rootArrayForApproximation,
                               int64_t NK, extraArraysPerThread ** eAPT, int64_t NT)
{
    omp_set_num_threads(NT);

    for (int64_t i = 0; i < forest->NV; i++)
        forest->totalBC[i] = 0.0;

    //printf("static threads here: %ld\n", NT);
    #pragma omp parallel for schedule(dynamic,1)

    for (int64_t i = 0; i < NK; i++) {
        int64_t thread = omp_get_thread_num();
        BrandesSingleTree(forest, sStinger, rootArrayForApproximation[i], forest->totalBC, eAPT[thread]);
    }

    #pragma omp parallel for

    for (int64_t v = 0; v < forest->NV; v++) {
        for (int64_t t = 0; t < NT; t++) {
            forest->totalBC[v] += eAPT[t]->sV[v].totalBC;
            eAPT[t]->sV[v].totalBC = 0.0;
        }
    }
}


int64_t BrandesSingleTree(bcForest * forest, struct stinger * sStinger,
                           int64_t currRoot, bc_t * totalBC, extraArraysPerThread * eAPT)
{
    bcTree * tree = forest->forest[currRoot];

    for (int64_t j = 0; j < tree->NV; j++) {
        tree->vArr[j].level = INFINITY_MY;
        tree->vArr[j].sigma = INFINITY_MY;
        tree->vArr[j].delta = 0;
        tree->vArr[j].edgesBelow = 0;
        tree->vArr[j].edgesAbove = 0;
    }

    tree->vArr[currRoot].level = 0;
    tree->vArr[currRoot].sigma = 1;

    int64_t * Stack = eAPT->Stack;
    int64_t * Queue = eAPT->QueueDown;

    Queue[0] = currRoot;
    int64_t qStart = 0, qEnd = 1;
    int64_t sStart = 0;
    int64_t k;

    // While queue is not empty
    while (qStart != qEnd) {
        int64_t currElement = Queue[qStart];
        Stack[sStart] = currElement;
        sStart++;
        qStart++;

        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement) {
            k = STINGER_EDGE_DEST;

            // If this is a neighbor and has not been found
            if (tree->vArr[k].level > tree->vArr[currElement].level) {
                tree->vArr[k].edgesAbove += tree->vArr[currElement].edgesAbove + 1;

                // Checking if "k" has been found.
                if (tree->vArr[k].level == INFINITY_MY) {
                    tree->vArr[k].level = tree->vArr[currElement].level + 1;
                    Queue[qEnd++] = k;
                    tree->vArr[k].delta = 0;
                }

                if (tree->vArr[k].sigma == INFINITY_MY) {
                    // k has not been found and therefore its paths to the roots are through its parent.
                    tree->vArr[k].sigma = tree->vArr[currElement].sigma;
                } else {
                    // k has been found and has multiple paths to the root as it has multiple parents.
                    tree->vArr[k].sigma += tree->vArr[currElement].sigma;
                }
            }
        }
        STINGER_FORALL_EDGES_OF_VTX_END();

#if COUNT_TRAVERSALS==1
        eAPT->staticTraverseEdgeCounter += stinger_typed_outdegree(sStinger, currElement, 0);
        eAPT->staticTraverseVerticeCounter++;
#endif

    }

    // Using Brandes algorithm to compute BC for a specific tree.
    // Essentially, use the stack which the elements are placed in depth-reverse order, to "climb" back
    // up the tree, all the way to the root.
    int64_t sEnd = sStart - 1;

    while (sEnd >= 0) {
        int64_t currElement = Stack[sEnd];
        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger, currElement) {
            k = STINGER_EDGE_DEST;

            // If this is a neighbor and has not been found
            if ((tree->vArr[k].level == (tree->vArr[currElement].level - 1))) {
                tree->vArr[k].delta +=
                    ((bc_t)tree->vArr[k].sigma / (bc_t)tree->vArr[currElement].sigma) *
                    (bc_t)(tree->vArr[currElement].delta + 1);
            } else if (tree->vArr[k].level == tree->vArr[currElement].level + 1) {
                tree->vArr[currElement].edgesBelow += tree->vArr[k].edgesBelow + 1;
            }
        }
        STINGER_FORALL_EDGES_OF_VTX_END();
#if COUNT_TRAVERSALS==1
        eAPT->staticTraverseEdgeCounter += stinger_typed_outdegree(sStinger, currElement, 0);
        eAPT->staticTraverseVerticeCounter++;
#endif

        if (currElement != currRoot) {
            eAPT->sV[currElement].totalBC += tree->vArr[currElement].delta;
        }

        sEnd--;
    }

    return -1;
}
