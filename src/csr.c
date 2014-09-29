#include "csr.h"

// NE number of undirected edges.
csrGraph* CreateCSRFromStinger(struct stinger* stingerGraph,int64_t NV,int64_t NE)
{
    csrGraph* newGraph = (csrGraph*)malloc(sizeof(csrGraph));
    newGraph->NV = NV;
    newGraph->NE = NE*2;

    newGraph->vertexPointerArray=(int64_t*)malloc(sizeof(int64_t)*(NV+1));
    newGraph->edgeArray=(int64_t*)malloc(sizeof(int64_t)*NE);

    int64_t edgeCounter=0;
    newGraph->vertexPointerArray[0]=0;
    for(int64_t v=0; v<NV;v++)
    {
        newGraph->vertexPointerArray[v+1]= newGraph->vertexPointerArray[v]+stinger_outdegree(stingerGraph,v);

        STINGER_FORALL_EDGES_OF_VTX_BEGIN(stingerGraph,v)
        {
            newGraph->edgeArray[edgeCounter]=STINGER_EDGE_DEST;
            edgeCounter++;
        }
        STINGER_FORALL_EDGES_OF_VTX_END();
    }

//    printf("%d,%ld,%d",NE*2,newGraph->vertexPointerArray[NV],edgeCounter);


    return newGraph;
}

// NE number of undirected edges.
void CreateStingerFromCSR(csrGraph* csr,struct stinger** stingerGraph)
{
        *stingerGraph=stinger_new();

    int64_t* offset=(int64_t*)malloc(sizeof(int64_t)*(csr->NV+1));
    int64_t* edges=(int64_t*)malloc(sizeof(int64_t)*(csr->NE));
    int64_t* weight=(int64_t*)malloc(sizeof(int64_t)*(csr->NE));

    for(int64_t v=0; v<=csr->NV;v++)
        offset[v]=csr->vertexPointerArray[v];

    for(int64_t e=0; e<csr->NE;e++)
    {
        edges[e]=csr->edgeArray[e];
        weight[e]=0;
    }


    stinger_set_initial_edges (*stingerGraph /* G */ ,
				csr->NV,
				0,
				offset ,
				edges ,
				weight,
				NULL ,
				NULL ,
				0);

	free(offset); free(edges); free(weight);
    return;
}



// NE number of undirected edges.
csrGraph* CreateEmptyUndirectedCSR(int64_t NV,int64_t NE)
{
    csrGraph* newGraph = (csrGraph*)malloc(sizeof(csrGraph));
    newGraph->NV = NV;
    newGraph->NE = NE*2;

    newGraph->vertexPointerArray=(int64_t*)malloc(sizeof(int64_t)*(newGraph->NV+1));
    newGraph->edgeArray=(int64_t*)malloc(sizeof(int64_t)*(newGraph->NE));

    newGraph->vertexPointerArray[0]=0;
    for(int64_t v=0; v<NV;v++)
    {
        newGraph->vertexPointerArray[v+1]= 0;
    }

    return newGraph;
}

csrGraph* CreateEmptyDirectedCSR(int64_t NV,int64_t NE)
{
    csrGraph* newGraph = (csrGraph*)malloc(sizeof(csrGraph));
    newGraph->NV = NV;
    newGraph->NE = NE;

    newGraph->vertexPointerArray=(int64_t*)malloc(sizeof(int64_t)*(newGraph->NV+1));
    newGraph->edgeArray=(int64_t*)malloc(sizeof(int64_t)*(newGraph->NE));

    newGraph->vertexPointerArray[0]=0;
    for(int64_t v=0; v<NV;v++)
    {
        newGraph->vertexPointerArray[v+1]= 0;
    }

    return newGraph;
}



void FreeCSR(csrGraph* graph)
{
    free(graph->vertexPointerArray);
    free(graph->edgeArray);
    free(graph);
}
