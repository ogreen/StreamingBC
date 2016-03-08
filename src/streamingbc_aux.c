
#include "streamingbc_aux.h"

//------------------------------------------
// Forest Creation

/**
* @brief Creates the data structures needed for computing BC.
*       These include the level of each vertex in the BFS tree.
*       The number of shortest paths each vertex has to the root.
*       The delta value computed in the dependency accumulation stage/
*
* @param numVertices The number of vertices in the graph
*
* @return Returns the created data structure.
*/
bcTree* CreateTree(int64_t numVertices){
    bcTree* newTree;
    newTree = (bcTree*)xmalloc(sizeof(bcTree));

    newTree->NV = numVertices;
    newTree->vArr = (bcV*)xmalloc(numVertices*sizeof(bcV));

    return newTree;
}

/**
* @brief Destroys the allocated tree.
*
* @param deadTree The tree that needs to be unallocated.
*
* @return None
*/
void DestroyTree(bcTreePtr deadTree){
    free(deadTree->vArr);
    free(deadTree);
}


/**
* @brief Creates a tree for each vertex in the graph.
*
* @param numVertices The number of vertices in the graph
* @param rootArray The array containing vertex id's to be considered for approximate case
* @param rootArraySize The number of vertices in the graph to be considered for approximate case
*
* @return None (see parameter newForest)
*/
bcForest* CreateForestExact(int64_t numVertices){
    bcForest* temp;
    temp = (bcForest*)xmalloc(sizeof(bcForest));
    temp->NV = numVertices;
    temp->forest = (bcTreePtr*)xmalloc(numVertices*sizeof(bcTreePtr));
    temp->totalBC = (bc_t*)xmalloc(numVertices*sizeof(bc_t));

    int64_t i;
    for(i=0;i<numVertices;i++)
    {
      temp->forest[i] = CreateTree(numVertices);
    }
    return temp;
}

/**
* @brief Destroys the allocated forest for exact computation.
*
* @param deadForest The forest that needs to be unallocated.
* @param rootArray The array containing vertex id's to be considered for approximate case
* @param rootArraySize The number of vertices in the graph to be considered for approximate case
*
* @return None
*/
void DestroyForestExact(bcForestPtr* deadForest){
    bcForest* temp=*deadForest;
    int64_t i;
    for(i=0;i<(temp->NV);i++)
    {
        DestroyTree(temp->forest[i]);
    }
   free(temp->totalBC);
   free(temp->forest);
   free(temp);
    return;
}



/**
* @brief Creates a tree for each vertex in the rootArray. Thus, a total of rootArraySize trees are created.
*
* @param numVertices The number of vertices in the graph
* @param rootArray The array containing vertex id's to be considered for approximate case
* @param rootArraySize The number of vertices in the graph to be considered for approximate case
*
* @return None (see parameter newForest)
*/
bcForest* CreateForestApproximate(int64_t numVertices, uint64_t* rootArray, uint64_t rootArraySize){
	bcForest* temp;;
    temp = (bcForest*)xmalloc(sizeof(bcForest));
    temp->NV = numVertices;
    temp->forest = (bcTreePtr*)xmalloc(numVertices*sizeof(bcTreePtr));
    temp->totalBC = (bc_t*)xmalloc(numVertices*sizeof(bc_t));

    int64_t i;
    for(i=0;i<rootArraySize;i++)
    {
      temp->forest[rootArray[i]] = CreateTree(numVertices);
    }
	return temp;
}

/**
* @brief Destroys the allocated forest for approximate case
*
* @param deadForest The forest that needs to be unallocated.
* @param rootArray The array containing vertex id's to be considered for approximate case
* @param rootArraySize The number of vertices in the graph to be considered for approximate case
*
* @return None
*/
void DestroyForestApproximate(bcForestPtr* deadForest, uint64_t* rootArray, uint64_t rootArraySize){
	bcForest* temp=*deadForest;
    int64_t i;
    for(i=0;i<(rootArraySize);i++)
    {
        DestroyTree(temp->forest[rootArray[i]]);
    }
   free(temp->totalBC);
   free(temp->forest);
   free(temp);
    return;
}


//----------------------------------------------
// Auxilary data structures

// Creates the multi-level queue for the computation of BC. In this case O(NV^2) is allocated.
uint64_t** createMultiLevelQueue(uint64_t NV){
    uint64_t** multiLevelQueue = (uint64_t**)malloc(NV*sizeof(uint64_t*));
    for(uint64_t v=0; v<NV;v++)    {
        multiLevelQueue[v]=(uint64_t*)malloc(NV*sizeof(uint64_t));
    }
    return multiLevelQueue;
}

// Destroys the multi level queue
void destroyMultiLevelQueue(uint64_t** multiLevelQueue,uint64_t NV){
    for(uint64_t v=0; v<NV;v++){
        free(multiLevelQueue[v]);
    }
    free(multiLevelQueue);
}

// Creates a multi level queue for each thread/core.
uint64_t*** createParallelMultiLevelQueue(uint64_t NV, uint64_t threadCount){
    uint64_t*** parallelMultiLevelQueue = (uint64_t***)malloc(threadCount*sizeof(uint64_t**));
    for(uint64_t t=0; t<threadCount;t++){
        parallelMultiLevelQueue[t]=createMultiLevelQueue(NV);
    }
    return parallelMultiLevelQueue;
}
// Destroys the multi level queue of each thread/core.
void destroyParallelMultiLevelQueue(uint64_t*** parallelMultiLevelQueue,uint64_t NV,uint64_t threadCount){
    for(uint64_t t=0; t<threadCount;t++){
        destroyMultiLevelQueue(parallelMultiLevelQueue[t],NV);
    }
    free(parallelMultiLevelQueue);
}

// Creates a parent list for each thread/core.
list_ptr** createParallelList(int64_t threadCount,int64_t NV){
    list_ptr** parallelList = (list_ptr**)xcalloc(threadCount,sizeof(list_ptr**));
    for(int64_t i=0; i<threadCount;i++){
        makeArrayOfLists(&parallelList[i],NV);
    }
    return parallelList;
}

// Destroys the parent list of each thread/core.
void destroyParallelList(list_ptr** parallelList, int64_t threadCount,int64_t NV){
    for(int64_t i=0; i<threadCount;i++){
        destroyArrayOfLists(&parallelList[i],NV);
    }
    free(parallelList);
}


float** createParallelBetweennessArray(int64_t threadCount,int64_t NV){
    float** totalBC = (float**)malloc((threadCount)*sizeof(float*));
    for(int64_t i=0;i<threadCount;i++){
        totalBC[i] = malloc(sizeof(float)*NV);
    }
    return totalBC;
}

void destroyParallelBetweennessArray(float** parallelScore, int64_t threadCount){
    for(int64_t i=0; i<threadCount;i++){
        free(parallelScore[i]);
    }
    free(parallelScore);
}

extraArraysPerThread* createExtraArraysPerThread(int64_t NV){
    extraArraysPerThread* eapt = (extraArraysPerThread*)malloc(sizeof(extraArraysPerThread));
    eapt->sV = (sbcV*)malloc(NV*sizeof(sbcV));
    
    for(int64_t v=0; v<NV; v++){
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
    eapt->touchedVerticesUp = (int64_t*) malloc(NV * sizeof(int64_t));
    eapt->touchedVerticesDown = (int64_t*) malloc(NV * sizeof(int64_t));
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

void ClearCounters(extraArraysPerThread* eapt){
    eapt->samelevelCounter=0;
    eapt->compConnCounter=0;
    eapt->adjacentCounter=0;
    eapt->movementCounter=0;

    eapt->staticTraverseVerticeCounter=0;
    eapt->staticTraverseEdgeCounter=0;
    eapt->dynamicTraverseVerticeCounter=0;
    eapt->dynamicTraverseEdgeCounter=0;

}
void destroyExtraArraysPerThread(extraArraysPerThread* eapt,int64_t NV){
    free(eapt->sV);
    free(eapt->QueueDown);
    free(eapt->QueueUp);
    free(eapt->QueueSame);
    free(eapt->Stack);

    destroyArrayOfLists(&eapt->multiLevelQueues,NV);
    free(eapt);
}

extraArraysPerThread** createExtraArraysForThreads(int64_t threadCount,int64_t NV){
    extraArraysPerThread** extraInfoArray = (extraArraysPerThread**)xcalloc(threadCount,sizeof(extraArraysPerThread*));
    for(int64_t i=0; i<threadCount;i++){
        extraInfoArray[i] = createExtraArraysPerThread(NV);
    }
    return extraInfoArray;
}

void destroyExtraArraysForThreads(extraArraysPerThread** parallelExtra, int64_t threadCount, int64_t NV){
    for(int64_t i=0; i<threadCount;i++){
        destroyExtraArraysPerThread(parallelExtra[i],NV);
    }
    free(parallelExtra);
}





//--------------------------------------------------------------
//--------------------------------------------------------------
// List


list_t* makeList()
{
    list_t *l  = (list_t*)xmalloc(sizeof(list_t));
    l->head=NULL;
    l->tail=NULL;
    l->size=0;
    return l;
}

node_t* makeNode(int64_t id)
{
    node_t *newnode = (node_t*)xmalloc(sizeof(node_t));
    newnode->id = id;
    newnode->next = NULL;
    return newnode;
}

/* note append will append at the last. */
void append(list_t *L, node_t *n)
{
    if(L->size==0)
    {
        L->head=n;
        L->tail=n;
    }
    else
    {
        L->tail->next = n;
        L->tail = n;
    }
    L->size++;

}

node_t* getFirst(list_t *L)
{
    return L->head;
}

void deleteFirst(list_t *L)
{
    node_t *n = L->head;
    L->head = n->next;
    free(n);
    L->size--;
}

void printList(list_t *L)
{
    node_t *n;
    printf("Printing list of size:%ld\n",L->size);
    if(L->size==0) return;
    for(n=L->head; n!=L->tail; n=n->next)
    {
        printf("%ld,",n->id);
    }
    printf("%ld,",n->id);
    printf("\n\n\n\n");
}


void emptyList(list_t* L)
{
    while(L->size >0)
    {
        node_t *n = L->head;
        L->head = n->next;
        free(n);
        L->size--;
    }
}


void makeArrayOfLists(list_ptr** aL,int64_t numberOfLists)
{
    int64_t i;

    *aL = (list_ptr*)xmalloc(numberOfLists*sizeof(list_ptr*));

    for(i=0;i<numberOfLists;i++)
    {
        (*aL)[i] = makeList();
    }
}

void destroyArrayOfLists(list_ptr** aL,int64_t numberOfLists)
{
    int64_t i;
    for(i=0;i<numberOfLists;i++)
    {
        emptyList((*aL)[i]);
        free((*aL)[i]);
    }

    free(*aL);
}

