
#include "dsUtils.h"
#include "xmalloc.h"
#include "streaming_insertions.h"

// Creates the parent array needed for the computation of BC.
// The output of this function is an array of arrays. Each vertex is allocated an array based on it's adjaency.
uint64_t** createParentArray(csrGraph* graph,uint64_t NV)
{
    uint64_t** parentArray = (uint64_t**)malloc(NV*sizeof(uint64_t*));
    for(uint64_t v=0; v<NV;v++)
    {
        uint64_t edgeCount=graph->vertexPointerArray[v+1]-graph->vertexPointerArray[v];
        parentArray[v]=(uint64_t*)malloc(edgeCount*sizeof(uint64_t));
    }

    return parentArray;
}

// Destroys the parent array
void destroyParentArray(uint64_t** parentArray,uint64_t NV)
{

    for(uint64_t v=0; v<NV;v++)
    {
        free(parentArray[v]);
    }

    free(parentArray);
}

// Creates a parent array for each thread/core.
uint64_t*** createParallelParentArray(csrGraph* graph,uint64_t NV, uint64_t threadCount)
{
    uint64_t*** parallelParentArray = (uint64_t***)malloc(threadCount*sizeof(uint64_t**));

    for(uint64_t t=0; t<threadCount;t++)
    {
        parallelParentArray[t]=createParentArray(graph,NV);
        if(parallelParentArray[t]==NULL)
            printf("Failed to allocated memory for parallel parent array\n");
    }

    return parallelParentArray;
}

// Destroys the parent array of each thread/core.
void destroyParallelParentArray(uint64_t*** parallelParentArray,uint64_t NV,uint64_t threadCount)
{

    for(uint64_t t=0; t<threadCount;t++)
    {
        destroyParentArray(parallelParentArray[t],NV);
    }

    free(parallelParentArray);
}


 // Creates the parent array needed for the computation of BC.
// The output of this function is an array of arrays. Each vertex is allocated an array based on it's adjaency.
uint64_t** createParentArrayStinger(struct stinger* S,uint64_t NV)
{
    uint64_t** parentArray = (uint64_t**)malloc(NV*sizeof(uint64_t*));
	for(uint64_t v=0; v<NV;v++)
    {
        uint64_t edgeCount=stinger_outdegree(S,v);
        parentArray[v]=(uint64_t*)malloc(edgeCount*sizeof(uint64_t));
    }

    return parentArray;
}


// Creates a parent array for each thread/core.
uint64_t*** createParallelParentArrayStinger(struct stinger* S,uint64_t NV, uint64_t threadCount)
{
    uint64_t*** parallelParentArray = (uint64_t***)malloc(threadCount*sizeof(uint64_t**));

 //   	printf("thread count %d   %d \n", threadCount,NV);

	for(uint64_t t=0; t<threadCount;t++)
    {
//		printf("thread count %d\n", threadCount);
        parallelParentArray[t]=createParentArrayStinger(S,NV);
        if(parallelParentArray[t]==NULL)
            printf("Failed to allocated memory for parallel parent array\n");
    }

    return parallelParentArray;
}








// Creates the multi-level queue for the computation of BC. In this case O(NV^2) is allocated.
uint64_t** createMultiLevelQueue(uint64_t NV)
{
    uint64_t** multiLevelQueue = (uint64_t**)malloc(NV*sizeof(uint64_t*));

    for(uint64_t v=0; v<NV;v++)
    {
        multiLevelQueue[v]=(uint64_t*)malloc(NV*sizeof(uint64_t));
    }

    return multiLevelQueue;
}

// Destroys the multi level queue
void destroyMultiLevelQueue(uint64_t** multiLevelQueue,uint64_t NV)
{

    for(uint64_t v=0; v<NV;v++)
    {
        free(multiLevelQueue[v]);
    }

    free(multiLevelQueue);
}

// Creates a multi level queue for each thread/core.
uint64_t*** createParallelMultiLevelQueue(uint64_t NV, uint64_t threadCount)
{
    uint64_t*** parallelMultiLevelQueue = (uint64_t***)malloc(threadCount*sizeof(uint64_t**));

    for(uint64_t t=0; t<threadCount;t++)
    {
        parallelMultiLevelQueue[t]=createMultiLevelQueue(NV);
    }

    return parallelMultiLevelQueue;
}
// Destroys the multi level queue of each thread/core.
void destroyParallelMultiLevelQueue(uint64_t*** parallelMultiLevelQueue,uint64_t NV,uint64_t threadCount)
{

    for(uint64_t t=0; t<threadCount;t++)
    {
        destroyMultiLevelQueue(parallelMultiLevelQueue[t],NV);
    }

    free(parallelMultiLevelQueue);
}




// Creates a parent list for each thread/core.
list_ptr** createParallelList(int64_t threadCount,int64_t NV)
{
    list_ptr** parallelList = (list_ptr**)xcalloc(threadCount,sizeof(list_ptr**));
    for(int64_t i=0; i<threadCount;i++)
    {
        makeArrayOfLists(&parallelList[i],NV);
    }

    return parallelList;
}

// Destroys the parent list of each thread/core.
void destroyParallelList(list_ptr** parallelList, int64_t threadCount,int64_t NV)
{
    for(int64_t i=0; i<threadCount;i++)
    {
        destroyArrayOfLists(&parallelList[i],NV);
    }

    free(parallelList);
}


float** createParallelBetweennessArray(int64_t threadCount,int64_t NV)
{
    float** totalBC = (float**)malloc((threadCount)*sizeof(float*));

    for(int64_t i=0;i<threadCount;i++)
    {
        totalBC[i] = malloc(sizeof(float)*NV);
    }
    return totalBC;
}

void destroyParallelBetweennessArray(float** parallelScore, int64_t threadCount)
{
    for(int64_t i=0; i<threadCount;i++)
    {
//        print64_tf("destroying %d\n",i); fflush(stdout);
        free(parallelScore[i]);
    }
//        printf("destroying last\n"); fflush(stdout);
    free(parallelScore);
}

extraArraysPerThread* createExtraArraysPerThread(int64_t NV)
{
    extraArraysPerThread* eapt = (extraArraysPerThread*)malloc(sizeof(extraArraysPerThread));

    eapt->sV = (sbcV*)malloc(NV*sizeof(sbcV));
	
	for(int64_t v=0; v<NV; v++)
	{
        	
		eapt->sV[v].diffPath=0;
		eapt->sV[v].touched=0;
		eapt->sV[v].newPathsToRoot=0;
		eapt->sV[v].movementDelta=0;
		eapt->sV[v].newDelta=0;
		eapt->sV[v].totalBC=0;
		eapt->sV[v].IMoved=-1;
	}

    eapt->QueueDown = (int64_t*)malloc(NV*sizeof(int64_t));
    eapt->QueueUp = (int64_t*)malloc(NV*sizeof(int64_t));
    eapt->QueueSame = (int64_t*)malloc(NV*sizeof(int64_t));
    eapt->Stack = (int64_t*)malloc(NV*sizeof(int64_t));

    makeArrayOfLists(&eapt->multiLevelQueues,NV);

    eapt->samelevelCounter=0;
    eapt->compConnCounter=0;
    eapt->adjacentCounter=0;
    eapt->movementCounter=0;
    eapt->staticTraverseVerticeCounter=0;
    eapt->staticTraverseEdgeCounter=0;
    eapt->dynamicTraverseVerticeCounter=0;
    eapt->dynamicTraverseEdgeCounter=0;

	return eapt;
}

void ClearCounters(extraArraysPerThread* eapt)
{
    eapt->samelevelCounter=0;
    eapt->compConnCounter=0;
    eapt->adjacentCounter=0;
    eapt->movementCounter=0;

    eapt->staticTraverseVerticeCounter=0;
    eapt->staticTraverseEdgeCounter=0;
    eapt->dynamicTraverseVerticeCounter=0;
    eapt->dynamicTraverseEdgeCounter=0;

}
void destroyExtraArraysPerThread(extraArraysPerThread* eapt,int64_t NV)
{
    free(eapt->sV);
    free(eapt->QueueDown);
    free(eapt->QueueUp);
    free(eapt->QueueSame);
    free(eapt->Stack);


	destroyArrayOfLists(&eapt->multiLevelQueues,NV);
    free(eapt);
}

extraArraysPerThread** createExtraArraysForThreads(int64_t threadCount,int64_t NV)
{
    extraArraysPerThread** extraInfoArray = (extraArraysPerThread**)xcalloc(threadCount,sizeof(extraArraysPerThread*));
    for(int64_t i=0; i<threadCount;i++)
    {
        extraInfoArray[i] = createExtraArraysPerThread(NV);
    }

    return extraInfoArray;
}

void destroyExtraArraysForThreads(extraArraysPerThread** parallelExtra, int64_t threadCount, int64_t NV)
{
    for(int64_t i=0; i<threadCount;i++)
    {
        destroyExtraArraysPerThread(parallelExtra[i],NV);
    }

    free(parallelExtra);
}


