
#include "streaming_utils.h"

#include <alloca.h>
#include <stdio.h>
#include <stdlib.h>

#include "stinger.h"
#include "stinger-traversal.h"
#include "stinger-utils.h"
#include "xmalloc.h"

#include "csr.h"
#include "list.h"

#include "graphGen.h"

void uniformRandGraph(uint64_t NV, uint64_t NE, struct stinger ** stingerGraph, uint64_t randomSeed)
{
    if(randomSeed==0)
        srand(time(NULL));
    else
    srand(randomSeed);

    int64_t* edgeSrc = (int64_t*)xmalloc((NE*2)*sizeof(int64_t));
    int64_t* edgeDest = (int64_t*)xmalloc((NE*2)*sizeof(int64_t));
    int64_t* edgeWeight = (int64_t*)xmalloc((NE*2)*sizeof(int64_t));

	if(edgeSrc == NULL || edgeDest==NULL || edgeWeight==NULL)
		printf("NULL Pointer in random graph creation\n");
    int64_t found;

	for(uint64_t i = 0; i < NE; i++) {
		uint64_t u = rand() % NV;
		uint64_t v = rand() % NV;

        found=0;

		if(!found)
		{
			edgeSrc[2*i]=u;edgeDest[2*i]=v;
			edgeSrc[2*i+1]=v;edgeDest[2*i+1]=u;
            edgeWeight[2*i]=0; edgeWeight[2*i+1]=0;

		}
		else
		{
		    i--;
		}
	}
    *stingerGraph = edge_list_to_stinger(NV,2*NE,edgeSrc,edgeDest,edgeWeight,NULL,NULL,0);
    free(edgeSrc); free(edgeDest); free(edgeWeight);

}


void FillAdjacencyMat(uint64_t NV, uint64_t NE, uint64_t * graph[],struct stinger * stingerGraph, uint64_t randomSeed)
{
    for(int64_t v=0; v<NV;v++)
    {
        STINGER_FORALL_EDGES_OF_VTX_BEGIN(stingerGraph,v)
        {
            graph[STINGER_EDGE_SOURCE][STINGER_EDGE_DEST]=1;
        }
        STINGER_FORALL_EDGES_OF_VTX_END();
    }

}



void uniformRandGraphMin2Components(uint64_t NV, uint64_t NE, uint64_t * graph[], struct stinger ** stingerGraph, uint64_t randomSeed, float percentageOfBigComponent)
{
    if(randomSeed==0)
        srand(time(NULL));
    else
    srand(randomSeed);

    int64_t* edgeSrc = (int64_t*)xmalloc((NE*2)*sizeof(int64_t));
    int64_t* edgeDest = (int64_t*)xmalloc((NE*2)*sizeof(int64_t));
    int64_t* edgeWeight = (int64_t*)xmalloc((NE*2)*sizeof(int64_t));

	if(edgeSrc == NULL || edgeDest==NULL || edgeWeight==NULL)
		printf("NULL Pointer in random graph creation\n");

    // Highest value vertex in the big component
    int64_t maxVertex = NV*percentageOfBigComponent;

    int64_t inSameComponent;

	for(uint64_t i = 0; i < NE; i++) {
		uint64_t u = rand() % NV;
		uint64_t v = rand() % NV;
		inSameComponent = 0;
		if((u<=maxVertex && v<=maxVertex) || (u>maxVertex && v>maxVertex) )
            inSameComponent = 1;

		if(graph[u][v] == 0 && u != v && inSameComponent) {
			graph[u][v] = 1;
			graph[v][u] = 1;

			edgeSrc[2*i]=u;edgeDest[2*i]=v;
			edgeSrc[2*i+1]=v;edgeDest[2*i+1]=u;
            edgeWeight[2*i]=0; edgeWeight[2*i+1]=0;
		} else {
			i--;
		}
	}

    *stingerGraph = edge_list_to_stinger(NV,2*NE,edgeSrc,edgeDest,edgeWeight,NULL,NULL,0);
    free(edgeSrc); free(edgeDest); free(edgeWeight);

}


void rMatRandGraph(uint64_t NV, uint64_t NE, uint64_t * graph[],struct stinger ** stingerGraph, uint64_t randomSeed,int64_t SCALE)
{
    if(randomSeed==0)
        srand(time(NULL));
    else
    srand(randomSeed);

    int64_t* edgeSrc = (int64_t*)xmalloc((NE*2)*sizeof(int64_t));
    int64_t* edgeDest = (int64_t*)xmalloc((NE*2)*sizeof(int64_t));
    int64_t* edgeWeight = (int64_t*)xmalloc((NE*2)*sizeof(int64_t));

	if(edgeSrc == NULL || edgeDest==NULL || edgeWeight==NULL)
		printf("NULL Pointer in random graph creation\n");


    int64_t iout, jout;
    double A=0.55, B=0.15, C=0.15, D=0.15;


	for(uint64_t n = 0; n < NE; n++)
	{
		uint64_t u,v;

        static int64_t xxx = 0;
        size_t rni = 0;
        int64_t i = 0, j = 0;
        int64_t bit = ((int64_t) 1) << (SCALE - 1);


        while (1) {
            const double r = (double)rand()/(double)RAND_MAX;
            if (r > A) {                /* outside quadrant 1 */
              if (r <= A + B)           /* in quadrant 2 */
                j |= bit;
              else if (r <= A + B + C)  /* in quadrant 3 */
                i |= bit;
              else {                    /* in quadrant 4 */
                j |= bit;
                i |= bit;
              }
            }
        if (1 == bit)
          break;

       bit >>= 1;
        }
        /* Iterates SCALE times. */
        iout = i;
        jout = j;

        u=iout; v=jout;


        if(graph[u][v] == 0 && u != v)
        {
            graph[u][v] = 1;
            graph[v][u] = 1;
            edgeSrc[2*n]=u;edgeDest[2*n]=v;
            edgeSrc[2*n+1]=v;edgeDest[2*n+1]=u;
            edgeWeight[2*n]=0; edgeWeight[2*n+1]=0;

        }
        else
        {
            n--;
        }



	}
        printf("BAAH\n\n\n");

    *stingerGraph = edge_list_to_stinger(NV,2*NE,edgeSrc,edgeDest,edgeWeight,NULL,NULL,0);
    free(edgeSrc); free(edgeDest); free(edgeWeight);
}


void CreateRMatEdge(uint64_t* out_u, uint64_t* out_v, uint64_t NV,int64_t SCALE)
{

    int64_t iout, jout;
    {
        double A=0.55, B=0.1, C=0.1, D=0.25;

		uint64_t u,v;

        int64_t i = 0, j = 0;
        int64_t bit = ((int64_t) 1) << (SCALE - 1);

        while (1) {
            const double r = (double)rand()/(double)RAND_MAX;
            if (r > A) {                // outside quadrant 1
              if (r <= A + B)           // in quadrant 2
                j |= bit;
              else if (r <= A + B + C)  // in quadrant 3
                i |= bit;
              else {                    // in quadrant 4
                j |= bit;
                i |= bit;
              }
            }
        if (1 == bit)
          break;

        bit >>= 1;
        }
        // Iterates SCALE times.
        iout = i;
        jout = j;

        *out_u =iout; *out_v=jout;

    }

}

void CreateUniformEdge(uint64_t* out_u, uint64_t* out_v, uint64_t NV,uint64_t* graph[])
{
    while(1)
    {
        *out_u = rand() % NV;*out_v = rand() % NV;

        if(graph[*out_u][*out_v] == 0 && *out_u!=*out_v) {
            return;
        }
    }
}

// This function doesn't allow the creation of duplicate edges.
void CreateCSRGraph(csrGraph* csrG,uint64_t NV, uint64_t NE,uint64_t randomSeed,graphType gt, uint64_t SCALE)
{
// uint64_t * graph[],struct stinger ** stingerGraph, int64_t SCALE)
    if(randomSeed==0)
        srand(time(NULL));
    else
        srand(randomSeed);

    list_ptr* parentList;
    makeArrayOfLists(&parentList,NV);

    uint64_t out_u ,out_v,e;
    int64_t found;
    for(e=0; e<NE;e++)
    {
        found=0;

        switch(gt)
        {
            case GT_ER:
                out_u = rand() % NV;
                out_v = rand() % NV;
                break;
            case GT_RMAT:
                CreateRMatEdge(&out_u,&out_v,NV,SCALE);
                break;
            default:
                break;
        }

            append(parentList[out_u],makeNode(out_v));
            append(parentList[out_v],makeNode(out_u));
        continue;

        node_t* iterator = parentList[out_u]->head;
        while(iterator!=NULL && !found)
        {
            uint64_t k = iterator->id;
            iterator = iterator->next;

            if(k==out_v)
                found=1;
        }
        if(!found)
        {
            append(parentList[out_u],makeNode(out_v));
            append(parentList[out_v],makeNode(out_u));
        }
        else
            e--;
    }


    int64_t edgeCounter=0, edgeCounterVertex;
    csrG->vertexPointerArray[0]=0;
    for(int64_t v=0; v<NV;v++)
    {
        edgeCounterVertex=0;

        node_t* iterator = parentList[v]->head;
        while(iterator!=NULL)
        {
            uint64_t k = iterator->id;
            iterator = iterator->next;

            csrG->edgeArray[edgeCounter]=k;
            edgeCounter++;
            edgeCounterVertex++;
        }
        csrG->vertexPointerArray[v+1]= csrG->vertexPointerArray[v] + edgeCounterVertex;
    }

    destroyArrayOfLists(&parentList,NV);
}
