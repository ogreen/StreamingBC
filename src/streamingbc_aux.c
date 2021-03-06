
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
bcTree * CreateTree(int64_t numVertices)
{
    bcTree * newTree;
    newTree = (bcTree *)xmalloc(sizeof(bcTree));

    newTree->NV = numVertices;
    newTree->vArr = (bcV *)xmalloc(numVertices * sizeof(bcV));

    return newTree;
}

/**
* @brief Destroys the allocated tree.
*
* @param deadTree The tree that needs to be unallocated.
*
* @return None
*/
void DestroyTree(bcTreePtr deadTree)
{
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
bcForest * CreateForestExact(int64_t numVertices)
{
    bcForest * temp;
    temp = (bcForest *)xmalloc(sizeof(bcForest));
    temp->NV = numVertices;
    temp->forest = (bcTreePtr *)xmalloc(numVertices * sizeof(bcTreePtr));
    temp->totalBC = (bc_t *)xmalloc(numVertices * sizeof(bc_t));

    int64_t i;

    for (i = 0; i < numVertices; i++) {
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
void DestroyForestExact(bcForestPtr * deadForest)
{
    bcForest * temp = *deadForest;
    int64_t i;

    for (i = 0; i < (temp->NV); i++) {
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
bcForest * CreateForestApproximate(int64_t numVertices, int64_t * rootArray, int64_t rootArraySize)
{
    bcForest * temp;;
    temp = (bcForest *)xmalloc(sizeof(bcForest));
    temp->NV = numVertices;
    temp->forest = (bcTreePtr *)xmalloc(numVertices * sizeof(bcTreePtr));
    temp->totalBC = (bc_t *)xmalloc(numVertices * sizeof(bc_t));

    int64_t i;

    for (i = 0; i < rootArraySize; i++) {
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
void DestroyForestApproximate(bcForestPtr * deadForest, int64_t * rootArray, int64_t rootArraySize)
{
    bcForest * temp = *deadForest;
    int64_t i;

    for (i = 0; i < (rootArraySize); i++) {
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
int64_t ** createMultiLevelQueue(int64_t NV)
{
    int64_t ** multiLevelQueue = (int64_t **)malloc(NV * sizeof(int64_t *));

    for (int64_t v = 0; v < NV; v++)    {
        multiLevelQueue[v] = (int64_t *)malloc(NV * sizeof(int64_t));
    }

    return multiLevelQueue;
}

// Destroys the multi level queue
void destroyMultiLevelQueue(int64_t ** multiLevelQueue, int64_t NV)
{
    for (int64_t v = 0; v < NV; v++) {
        free(multiLevelQueue[v]);
    }

    free(multiLevelQueue);
}

// Creates a multi level queue for each thread/core.
int64_t ** * createParallelMultiLevelQueue(int64_t NV, int64_t threadCount)
{
    int64_t ** * parallelMultiLevelQueue = (int64_t ** *)malloc(threadCount * sizeof(int64_t **));

    for (int64_t t = 0; t < threadCount; t++) {
        parallelMultiLevelQueue[t] = createMultiLevelQueue(NV);
    }

    return parallelMultiLevelQueue;
}
// Destroys the multi level queue of each thread/core.
void destroyParallelMultiLevelQueue(int64_t ** * parallelMultiLevelQueue, int64_t NV, int64_t threadCount)
{
    for (int64_t t = 0; t < threadCount; t++) {
        destroyMultiLevelQueue(parallelMultiLevelQueue[t], NV);
    }

    free(parallelMultiLevelQueue);
}

// Creates a parent list for each thread/core.
list_ptr ** createParallelList(int64_t threadCount, int64_t NV)
{
    list_ptr ** parallelList = (list_ptr **)xcalloc(threadCount, sizeof(list_ptr **));

    for (int64_t i = 0; i < threadCount; i++) {
        makeArrayOfLists(&parallelList[i], NV);
    }

    return parallelList;
}

// Destroys the parent list of each thread/core.
void destroyParallelList(list_ptr ** parallelList, int64_t threadCount, int64_t NV)
{
    for (int64_t i = 0; i < threadCount; i++) {
        destroyArrayOfLists(&parallelList[i], NV);
    }

    free(parallelList);
}


float ** createParallelBetweennessArray(int64_t threadCount, int64_t NV)
{
    float ** totalBC = (float **)malloc((threadCount) * sizeof(float *));

    for (int64_t i = 0; i < threadCount; i++) {
        totalBC[i] = (float *)malloc(sizeof(float) * NV);
    }

    return totalBC;
}

void destroyParallelBetweennessArray(float ** parallelScore, int64_t threadCount)
{
    for (int64_t i = 0; i < threadCount; i++) {
        free(parallelScore[i]);
    }

    free(parallelScore);
}

extraArraysPerThread * createExtraArraysPerThread(int64_t NV)
{
    extraArraysPerThread * eapt = (extraArraysPerThread *)malloc(sizeof(extraArraysPerThread));
    eapt->sV = (sbcV *)malloc(NV * sizeof(sbcV));

    for (int64_t v = 0; v < NV; v++) {
        eapt->sV[v].diffPath = 0;
        eapt->sV[v].touched = 0;
        eapt->sV[v].newSigma = 0;
        eapt->sV[v].movementDelta = 0;
        eapt->sV[v].newDelta = 0;
        eapt->sV[v].totalBC = 0;
        eapt->sV[v].IMoved = -1;
        eapt->sV[v].newEdgesAbove = 0;
        eapt->sV[v].newEdgesBelow = 0;
        eapt->sV[v].newLevel = 0;
    }

    eapt->QueueDown = (int64_t *)malloc(NV * sizeof(int64_t));
    eapt->QueueUp = (int64_t *)malloc(NV * sizeof(int64_t));
    eapt->QueueSame = (int64_t *)malloc(2 * NV * sizeof(int64_t));
    eapt->Stack = (int64_t *)malloc(2 * NV * sizeof(int64_t));
    eapt->touchedVerticesUp = (int64_t *) malloc(NV * sizeof(int64_t));
    eapt->touchedVerticesDown = (int64_t *) malloc(NV * sizeof(int64_t));
    eapt->tqBorders = (int64_t *) malloc(NV * sizeof(int64_t));
    //makeArrayOfLists(&eapt->multiLevelQueues,NV);
    eapt->queue = (queue_t *) malloc(sizeof(queue_t));
    eapt->queue->size = 0;
    eapt->queue->nodes = (queue_node_t *) malloc(NV * sizeof(queue_node_t));
    eapt->levelIndices = (level_node_t *) malloc(NV * sizeof(level_node_t));

    for (int64_t i = 0; i < NV; i++) {
        eapt->levelIndices[i].front = -1;
        eapt->levelIndices[i].back = -1;
        eapt->queue->nodes[i].next = -2;
        eapt->levelIndices[i].back_flag = -1;
    }

    eapt->samelevelCounter = 0;
    eapt->compConnCounter = 0;
    eapt->adjacentCounter = 0;
    eapt->movementCounter = 0;
    eapt->staticTraverseVerticeCounter = 0;
    eapt->staticTraverseEdgeCounter = 0;
    eapt->dynamicTraverseVerticeCounter = 0;
    eapt->dynamicTraverseEdgeCounter = 0;

    eapt->qStart = 0;
    eapt->qEnd = 0;
    eapt->qStart_nxt = 0;
    eapt->qEnd_nxt = 0;

    eapt->qStartSame = 0;
    eapt->qEndSame = 0;
    eapt->qStartSame_nxt = 0;
    eapt->qEndSame_nxt = 0;

    eapt->tqStart = 0;
    eapt->tqEnd = 0;
    eapt->tqStart_nxt = 0;
    eapt->tqEnd_nxt = 0;
    return eapt;
}

void ClearCounters(extraArraysPerThread * eapt)
{
    eapt->samelevelCounter = 0;
    eapt->compConnCounter = 0;
    eapt->adjacentCounter = 0;
    eapt->movementCounter = 0;

    eapt->staticTraverseVerticeCounter = 0;
    eapt->staticTraverseEdgeCounter = 0;
    eapt->dynamicTraverseVerticeCounter = 0;
    eapt->dynamicTraverseEdgeCounter = 0;

}
void destroyExtraArraysPerThread(extraArraysPerThread * eapt, int64_t NV)
{
    free(eapt->sV);
    free(eapt->QueueDown);
    free(eapt->QueueUp);
    free(eapt->QueueSame);
    free(eapt->Stack);
    free(eapt->queue->nodes);
    free(eapt->queue);
    free(eapt->levelIndices);
    free(eapt->tqBorders);
    //destroyArrayOfLists(&eapt->multiLevelQueues,NV);
    free(eapt);
}

extraArraysPerThread ** createExtraArraysForThreads(int64_t threadCount, int64_t NV)
{
    extraArraysPerThread ** extraInfoArray = (extraArraysPerThread **)xcalloc(threadCount, sizeof(extraArraysPerThread *));

    for (int64_t i = 0; i < threadCount; i++) {
        extraInfoArray[i] = createExtraArraysPerThread(NV);
    }

    return extraInfoArray;
}

void destroyExtraArraysForThreads(extraArraysPerThread ** parallelExtra, int64_t threadCount, int64_t NV)
{
    for (int64_t i = 0; i < threadCount; i++) {
        destroyExtraArraysPerThread(parallelExtra[i], NV);
    }

    free(parallelExtra);
}





//--------------------------------------------------------------
//--------------------------------------------------------------
// List


list_t * makeList()
{
    list_t * l  = (list_t *)xmalloc(sizeof(list_t));
    l->head = NULL;
    l->tail = NULL;
    l->size = 0;
    return l;
}

node_t * makeNode(int64_t id)
{
    node_t * newnode = (node_t *)xmalloc(sizeof(node_t));
    newnode->id = id;
    newnode->next = NULL;
    return newnode;
}

void appendDS(queue_t * queue, level_node_t * levelIndices, int64_t level, int64_t data)
{
    appendDS2(queue, levelIndices, level, data, 0);
}

/* note will be appended at the end. */
void appendDS2(queue_t * queue, level_node_t * levelIndices, int64_t level, int64_t data, int64_t thread_num)
{

    int64_t next = __sync_fetch_and_add(&(queue->size), 1);
    queue_node_t * list = queue->nodes;

    while (1) {
        if (__sync_bool_compare_and_swap(&(levelIndices[level].back_flag), -1, thread_num)) {
            if (levelIndices[level].back == -1) {
                levelIndices[level].front = next;
                levelIndices[level].back = next;
                list[next].next = -1;
            } else {
                list[levelIndices[level].back].next = next;
                levelIndices[level].back = next;
                list[next].next = -1;
            }

            levelIndices[level].back_flag = -1;
            list[next].data = data;
            break;
        }
    }
}

queue_node_t * getFirstDS(queue_t * queue, level_node_t * levelIndices, int64_t level)
{
    int64_t index = levelIndices[level].front;

    if (index == -1) {
        return NULL;
    }

    return &queue->nodes[index];
}

void deleteFirstDS(queue_t * queue, level_node_t * levelIndices, int64_t level)
{
    int64_t front = levelIndices[level].front;
    int64_t back = levelIndices[level].back;

    if (front == back) {
        levelIndices[level].front = -1;
        levelIndices[level].back = -1;
    } else {
        int64_t next = queue->nodes[front].next;
        levelIndices[level].front = next;
        queue->nodes[front].next = -2;
    }
}

void printListDS(queue_t * queue, level_node_t * levelIndices, int64_t level)
{
    int64_t index = levelIndices[level].front;
    printf("DS: ");

    while (index != -1) {
        printf("%ld ", queue->nodes[index].data);
        index = queue->nodes[index].next;
    }

    printf("\n");
}

void append(list_t * L, node_t * n)
{
    if (L->size == 0) {
        L->head = n;
        L->tail = n;
    } else {
        L->tail->next = n;
        L->tail = n;
    }

    L->size++;

}

node_t * getFirst(list_t * L)
{
    return L->head;
}

void deleteFirst(list_t * L)
{
    node_t * n = L->head;
    L->head = n->next;
    free(n);
    L->size--;
}


void compareLists(list_t * L, queue_t * queue, level_node_t * levelIndices, int64_t level)
{
    int64_t index = levelIndices[level].front;

    if (index == -1 && L->size == 0) return;

    if (index == -1 && L->size != 0) {
        printf("Error: DS is empty when LL is of size: %ld\n", L->size);
    }

    if (index != -1 && L->size == 0) {
        printf("Error: LL is empty when DS is not.\n");
    }

    node_t * n = L->head;

    while (index != -1 && n != NULL) {
        int64_t dsElem = queue->nodes[index].data;
        int64_t listElem = n->id;

        if (dsElem != listElem) {
            printf("Error: lists differ at level: %ld, DS: %ld, List; %ld\n", level, dsElem, listElem);
            return;
        }

        n = n->next;
        index = queue->nodes[index].next;
    }
}

void printList(list_t * L)
{
    node_t * n;
    //printf("Printing list of size:%ld\n",L->size);
    printf("LL: ");

    if (L->size == 0) return;

    for (n = L->head; n != L->tail; n = n->next) {
        printf("%ld,", n->id);
    }

    printf("%ld,", n->id);
    printf("\n");
}


void emptyList(list_t * L)
{
    while (L->size > 0) {
        node_t * n = L->head;
        L->head = n->next;
        free(n);
        L->size--;
    }
}


void makeArrayOfLists(list_ptr ** aL, int64_t numberOfLists)
{
    int64_t i;

    *aL = (list_ptr *)xmalloc(numberOfLists * sizeof(list_ptr *));

    for (i = 0; i < numberOfLists; i++) {
        (*aL)[i] = makeList();
    }
}

void destroyArrayOfLists(list_ptr ** aL, int64_t numberOfLists)
{
    int64_t i;

    for (i = 0; i < numberOfLists; i++) {
        emptyList((*aL)[i]);
        free((*aL)[i]);
    }

    free(*aL);
}

