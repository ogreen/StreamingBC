#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <assert.h>

#include <time.h>
#include <math.h>

#include "streamingbc_aux.h"

#include "streaming_insertions.h"
#include "operation.h"

#include "timer.h"
#include "stinger-utils.h"
#include "stinger-traversal.h"
#include "xmalloc.h"

#include "omp.h"

#define PRINT_NEIGHBORS(u) printf("Neighbors of %ld : ",u);{for(uint64_t k = 0; k < NV; k++) { if(graphAfter[u][k] == 1) printf("%ld ,",k);}  printf("\n");}
#define PRINT_NEIGHBORS_GRAPH(u) printf("Neighbors of %ld : ",u);{for(uint64_t k = 0; k < NV; k++) { if(graphAfter[u][k] == 1) printf("%ld ,",k);}  printf("\n");}
#define PRINT_VERTICES(r) {int64_t s; for(s=0;s<NV;s++){printf("%ld,",s);} printf("\n");}
#define PRINT_LEVEL_TREE(r) {int64_t s; for(s=0;s<NV;s++){printf("%ld,",level[r][s]);} printf("\n");}

#define LINE_SIZE 100000
/*int64_t noComputeDiffComp; // Not really used?
testCase graphTestCase; // GraphType ER=0,RMAT=1,real graph=2, connec
int64_t SCALE; // used only for RMAT type graphs
uint64_t** parentArray;
uint64_t*  parentCounter;
// To be used on real graphs only
int64_t* srcVerToDelete;
int64_t* destVerToDelete;
*/

#define COUNT 50
#define INSERTING 0

//int64_t * rootArrayForApproximation;

typedef enum{
    UP_INSERT = 0,
    UP_DELETE,
} updateType;


typedef struct{
    int64_t NV;
    int64_t NE;
    int64_t* edgeArray;
    int64_t* vertexPointerArray;
} csrGraph;


// NE number of undirected edges.
csrGraph* CreateCSRFromStinger(struct stinger* stingerGraph,int64_t NV,int64_t NE);
void CreateStingerFromCSR(csrGraph* csr,struct stinger** stingerGraph);
//csrGraph* CreateEmptyCSR(int64_t NV,int64_t NE);
csrGraph* CreateEmptyUndirectedCSR(int64_t NV,int64_t NE);
csrGraph* CreateEmptyDirectedCSR(int64_t NV,int64_t NE);

void FreeCSR(csrGraph* graph);



void hostParseArgsVitalUpdate(int argc, char** argv, int64_t *NV, int64_t *NE, int64_t *NK, int64_t *NT,
                                        int64_t *randomSeed, int64_t *iterationCount, char *initial_graph_name[1024]);

void CreateRandomEdgeListFromGraph(struct stinger* stingerGraph, int64_t NV, int64_t* insertionArraySrc,
		int64_t* insertionArrayDest, int64_t insertionCount);

double updateEdgeNEW(struct stinger* stingerGraph,StreamingExtraInfo* oneSEI,
		extraArraysPerThread** eAPT_perThread, uint64_t * rootArrayForApproximation, int64_t NK,
		int64_t NV, int64_t NT, bcForest* beforeBCForest,int64_t u_new, int64_t v_new, int64_t *iterationCount) {

	uint64_t iterator;
	float timeFullBeforeMulti = 0, timeFullAfterMulti = 0, timeStreamMulti=0;//timeStream = 0, timeSummation = 0;
	float_t avgComponents = 0;
	fflush(stdout);
	iterationCount=1;
	double  timeFullBefore = 0, timeFullAfter = 0, timeStream=0;//timeStream = 0, timeSummation = 0;

	tic();
	*oneSEI = insertEdgeStreamingBC(beforeBCForest, stingerGraph, u_new, v_new, rootArrayForApproximation,NK,NV,NT,eAPT_perThread);
	timeStreamMulti = toc();

	return timeStreamMulti;
}

double deleteEdgeNEW(struct stinger *stingerGraph, StreamingExtraInfo* oneSEI, extraArraysPerThread **eAPT_perThread,
                uint64_t *rootArrayForApproximation, int64_t NK, int64_t NV, int64_t NT, 
                bcForest *beforeBCForest, int64_t u_old, int64_t v_old, int64_t *iterationCount)
{
    uint64_t iterator;
    float timeFullBeforeMulti = 0, timeFullAfterMulti = 0, timeStreamMulti = 0;
    float_t avgComponents = 0;

    fflush(stdout);
    iterationCount = 1;

    double timeFullBefore = 0, timeFullAfter = 0, timeStream = 0;

    tic();
    *oneSEI = deleteEdgeStreamingBC(beforeBCForest, stingerGraph, u_old, v_old, rootArrayForApproximation, NK, NV, NT,  eAPT_perThread);
    timeStreamMulti = toc();

    return timeStreamMulti;
}

int doubleCompareforMinSort (const void * a, const void * b){
	if ( *(double*)a <  *(double*)b ) return 1;
	if ( *(double*)a == *(double*)b ) return 0;
	if ( *(double*)a >  *(double*)b ) return -1;
}

int main(int argc, char *argv[])
{
	
        int64_t NV;
        int64_t NE;
        int64_t NK;
        int64_t NT;
        int64_t randomSeed;
        char initial_graph_name[1024];
        int64_t iterationCount;
        const int64_t threadArraySize = 1;
        int64_t threadArray[] = {1};//{1,5,10,15,20,25,30,35,40}; 
        
        int64_t insertionArraySrc[COUNT];
        int64_t insertionArrayDest[COUNT];
        int64_t deletionArraySrc[COUNT];
        int64_t deletionArrayDest[COUNT];

        hostParseArgsVitalUpdate(argc, argv, &NV, &NE, &NK, &NT, &randomSeed, &iterationCount, &initial_graph_name);
   
        printf("NV, NE, NK, NT: %ld, %ld, %ld, %ld\n", NV, NE, NK, NT); 
        printf("initial_graph_name: %s\n", initial_graph_name);
        double timingDynamic[threadArraySize][COUNT];
	double timingStatic[threadArraySize];
	double timingDynamicTotal[threadArraySize];

	StreamingExtraInfo seiDynamic[threadArraySize][COUNT];
	StreamingExtraInfo seiDynamicTotal[threadArraySize];

	int64_t staticTraverseVerticeCounter[threadArraySize];
	int64_t staticTraverseEdgeCounter[threadArraySize];
        
	int64_t dynamicTraverseVerticeCounter[threadArraySize][COUNT];
	int64_t dynamicTraverseEdgeCounter[threadArraySize][COUNT];
	int64_t dynamicTraverseVerticeCounterTotal[threadArraySize];
	int64_t dynamicTraverseEdgeCounterTotal[threadArraySize];
	int64_t dynamicTraverseEdgeCounterMax[threadArraySize][COUNT];
/*
        double timingDynamicDeleting[threadArraySize][DEL_COUNT];
        double timingStaticDeleting[threadArraySize];
        double timingDynamicTotalDeleting[threadArraySize];

        StreamingExtraInfo seiDynamicDeleting[threadArraySize][DEL_COUNT];
        StreamingExtraInfo seiDynamicTotalDeleting[threadArraySize];

        int64_t staticTraverseVerticeCounterDeleting[threadArraySize];
        int64_t staticTraverseEdgeCounterDeleting[threadArraySize];

        int64_t dynamicTraverseVerticeCounterDeleting[threadArraySize][DEL_COUNT];
        int64_t dynamicTraverseEdgeCounterDeleting[threadArraySize][DEL_COUNT];
        int64_t dynamicTraverseVerticeCounterTotalDeleting[threadArraySize];
        int64_t dynamicTraverseEdgeCounterTotalDeleting[threadArraySize];
        int64_t dynamicTraverseEdgeCounterMaxDeleting[threadArraySize][DEL_COUNT];
*/
	for (int64_t threadCount=0; threadCount<threadArraySize; threadCount++)
	{
		NT = threadArray[threadCount];
		if(randomSeed==0)
			srand(time(NULL));
		else
			srand(randomSeed);


		FILE * fp = fopen(initial_graph_name, "r");

		struct stinger* stingerGraph;
		char line[LINE_SIZE];


		// Read data from file
		if(1) {
			fgets(line, LINE_SIZE, fp);
			sscanf(line, "%ld %ld", &NV, &NE);

			NV=NV+1;
			NE=NE*2;

			csrGraph* CSRGraph = CreateEmptyDirectedCSR(NV,NE);

			CSRGraph->vertexPointerArray[0]=0;
			CSRGraph->vertexPointerArray[1]=0;
			int64_t counter=0;
			for(int64_t u = 1; fgets(line, LINE_SIZE, fp); u++) {
				uint64_t neigh=0;
				uint64_t v = 0;
				char * ptr = line;
				int read = 0;

				while(sscanf(ptr, "%ld%n", &v, &read) > 0) {
					ptr += read;
					neigh++;
					CSRGraph->edgeArray[counter++]=v;
				}
    
				CSRGraph->vertexPointerArray[u+1]=CSRGraph->vertexPointerArray[u]+neigh;
			}
                        
			CreateStingerFromCSR(CSRGraph,&stingerGraph);

			FreeCSR(CSRGraph);

		}

		fclose(fp);

		uint64_t * rootArrayForApproximation = NULL;
                uint64_t * staticTraverseEdgeCounterRoot = NULL;
		rootArrayForApproximation = (uint64_t*)xmalloc(sizeof(uint64_t)*NK);
                staticTraverseEdgeCounterRoot = (uint64_t *)xmalloc(sizeof(uint64_t) * NV);

		//assert (NK == NV);
        if(NK==NV){
        	//printf("checking this\n");
			for(int64_t vr=0;vr<NK;vr++)
				rootArrayForApproximation[vr] = vr;
		}
		else{
			int64_t* flag = (uint64_t*)xmalloc(sizeof(uint64_t)*NV);
			for(int64_t vv=0;vv<NV;vv++){
				flag[vv]=0;
			}

			for(int64_t vr=0;vr<NK;vr++){
                                int tempV=rand()%NV;
				if(flag[tempV]==1)
				{
					vr--;
					continue;
				}
				rootArrayForApproximation[vr] = tempV;
				flag[tempV]=1;
                                //            printf("%ld,",rootArrayForApproximation[vr]);
			}
			free(flag);

		}
		//    printf("***************\n");

		double timeFullBeforeMulti=0;


		bcForest* beforeBCForest=NULL;
                              
                #if INSERTING==1 
				CreateRandomEdgeListFromGraph(stingerGraph,NV,insertionArraySrc,insertionArrayDest,COUNT);
                #else
                CreateRandomEdgeListFromGraphDeleting(stingerGraph, NV, deletionArraySrc, 
                                            deletionArrayDest, COUNT);
                #endif
               
		if(0)
		{
			//-------Compute static BC
			//-------START

			beforeBCForest=CreateForestForApproxCase(&beforeBCForest, NV, rootArrayForApproximation, NK);

			extraArraysPerThread* eAPT = createExtraArraysPerThread(NV);

			tic();
			BrandesApproxCase(beforeBCForest,stingerGraph, rootArrayForApproximation, NK,eAPT);
			timeFullBeforeMulti=toc();

			destroyExtraArraysPerThread(eAPT,NV);
			//-------END
		}
		else
		{

			//-------Compute static BC - Parallel
			//-------START

			beforeBCForest=CreateForestForApproxCase(&beforeBCForest, NV, rootArrayForApproximation, NK);
			extraArraysPerThread** eAPT_perThread2 = createExtraArraysForThreads(NT,NV);

			tic();
			BrandesApproxCaseParallel(beforeBCForest,stingerGraph, rootArrayForApproximation, NK,eAPT_perThread2,NT);
                        timingStatic[threadCount]=toc();
			staticTraverseVerticeCounter[threadCount]=eAPT_perThread2[0]->staticTraverseVerticeCounter;
			staticTraverseEdgeCounter[threadCount]=eAPT_perThread2[0]->staticTraverseEdgeCounter;

			//   printf("********* %ld %ld  *********", staticTraverseVerticeCounter[threadCount], staticTraverseEdgeCounter[threadCount]);

			/*
			   printf("---------------------\n---------------------\n");
			   for(int64_t t=0;t<NT;t++)
			   printf("%ld %ld %ld\n", t,
			   eAPT_perThread2[t]->staticTraverseVerticeCounter,eAPT_perThread2[t]->staticTraverseEdgeCounter);
			   printf("---------------------\n---------------------\n");
			 */
			destroyExtraArraysForThreads(eAPT_perThread2,NT,NV);

			//-------END
		}
                 

		//-------Compute streaming BC
		//-------START
		StreamingExtraInfo oneSEI;
		StreamingExtraInfo globalSEI = {0,0,0,0};

		extraArraysPerThread** eAPT_perThread = createExtraArraysForThreads(NT,NV);

		double streamingTimeTotal=0.0,minTime=10000000000.0,maxTime=0.0,itTime=0.0;

		timingDynamicTotal[threadCount]=0.0;
		dynamicTraverseVerticeCounterTotal[threadCount]=0;
		dynamicTraverseEdgeCounterTotal[threadCount]=0;

		for(int64_t count=0; count<COUNT; count++)
		{
                        #if INSERTING==1
			int64_t src = insertionArraySrc[count];
			int64_t dest = insertionArrayDest[count];

			stinger_insert_edge(stingerGraph,0,src,dest,0,0);
			stinger_insert_edge(stingerGraph,0,dest,src,0,0);

			timingDynamic[threadCount][count] =updateEdgeNEW(stingerGraph,&oneSEI, eAPT_perThread, rootArrayForApproximation, 
                                                                            NK, NV, NT, beforeBCForest,src, dest, &iterationCount); 
                        #else
                        int64_t src = deletionArraySrc[count];
                        int64_t dest = deletionArrayDest[count];

                        stinger_remove_edge(stingerGraph, 0, src, dest);
                        stinger_remove_edge(stingerGraph, 0, dest, src);
                            
                        timingDynamic[threadCount][count] = deleteEdgeNEW(stingerGraph, &oneSEI, eAPT_perThread, rootArrayForApproximation, 
                                                                            NK, NV, NT, beforeBCForest, src, dest, &iterationCount);
                        #endif

			timingDynamicTotal[threadCount]+=timingDynamic[threadCount][count];

			globalSEI.adjacent += oneSEI.adjacent; globalSEI.movement += oneSEI.movement; globalSEI.sameLevel += oneSEI.sameLevel;
			seiDynamic[threadCount][count].movement = oneSEI.movement; 
			seiDynamic[threadCount][count].adjacent = oneSEI.adjacent; 
			seiDynamic[threadCount][count].sameLevel = oneSEI.sameLevel; 
			oneSEI.adjacent=0; oneSEI.movement=0; oneSEI.sameLevel=0;


			dynamicTraverseVerticeCounter[threadCount][count]=eAPT_perThread[0]->dynamicTraverseVerticeCounter;
			dynamicTraverseEdgeCounter[threadCount][count]=eAPT_perThread[0]->dynamicTraverseEdgeCounter;
			dynamicTraverseVerticeCounterTotal[threadCount]+=dynamicTraverseVerticeCounter[threadCount][count];
			dynamicTraverseEdgeCounterTotal[threadCount]+=dynamicTraverseEdgeCounter[threadCount][count];


			dynamicTraverseEdgeCounterMax[threadCount][count]= eAPT_perThread[0]->dynamicTraverseEdgeCounter;

			for(int t=0;t<NT;t++)                                                      
			{
				if(dynamicTraverseEdgeCounterMax[threadCount][count]<eAPT_perThread[t]->dynamicTraverseEdgeCounter)
					dynamicTraverseEdgeCounterMax[threadCount][count]=eAPT_perThread[t]->dynamicTraverseEdgeCounter;
				ClearCounters(eAPT_perThread[t]);
			}
		}

		seiDynamicTotal[threadCount].movement =globalSEI.movement; 
		seiDynamicTotal[threadCount].adjacent =globalSEI.adjacent; 
		seiDynamicTotal[threadCount].sameLevel =globalSEI.sameLevel; 

		bcForest* afterBCForest=NULL;

		afterBCForest=CreateForestForApproxCase(&afterBCForest, NV, rootArrayForApproximation, NK);
		extraArraysPerThread** eAPT_perThreadAfter = createExtraArraysForThreads(NT,NV);
		BrandesApproxCaseParallel(afterBCForest,stingerGraph, rootArrayForApproximation, NK,eAPT_perThreadAfter,NT);

		for(int a=0; a<beforeBCForest->NV; a++)
                {
                    //printf("before totalBC[%ld]: %lf\n", a, beforeBCForest->totalBC[a]);
                    //printf("after  totalBC[%ld]: %lf\n", a, afterBCForest->totalBC[a]);
                    if(beforeBCForest->totalBC[a]-afterBCForest->totalBC[a] > 0.001 ||
                        afterBCForest->totalBC[a] - beforeBCForest->totalBC[a] > 0.001) {
			printf("Error in computation %d, before: %lf  after: %lf\n", a,beforeBCForest->totalBC[a],afterBCForest->totalBC[a]);
	    	    }
                }
		destroyExtraArraysForThreads(eAPT_perThreadAfter,NT,NV);

                
                /*        	
                //-------Compute static BC - Parallel
                //-------START

                beforeBCForest=CreateForestForApproxCase(&beforeBCForest, NV, rootArrayForApproximation, NK);
                extraArraysPerThread** eAPT_perThread2 = createExtraArraysForThreads(NT,NV);

                tic();
                bfsBrandesForApproxCaseParallel(beforeBCForest,stingerGraph, rootArrayForApproximation, NK,eAPT_perThread2,NT);
                timingStatic[threadCount]=toc();
                staticTraverseVerticeCounter[threadCount]=eAPT_perThread2[0]->staticTraverseVerticeCounter;
                staticTraverseEdgeCounter[threadCount]=eAPT_perThread2[0]->staticTraverseEdgeCounter;

                //   printf("********* %ld %ld  *********", staticTraverseVerticeCounter[threadCount], staticTraverseEdgeCounter[threadCount]);

                *
                   printf("---------------------\n---------------------\n");
                   for(int64_t t=0;t<NT;t++)
                   printf("%ld %ld %ld\n", t,
                   eAPT_perThread2[t]->staticTraverseVerticeCounter,eAPT_perThread2[t]->staticTraverseEdgeCounter);
                   printf("---------------------\n---------------------\n");
                 *
                
                destroyExtraArraysForThreads(eAPT_perThread2,NT,NV);
                */
                destroyExtraArraysForThreads(eAPT_perThread,NT,NV);
		//-------END
                


		DestroyForestForApproxCase(&beforeBCForest,rootArrayForApproximation, NK);

		//-------END

		free(rootArrayForApproximation);

		stinger_free(stingerGraph);

	}
	if(1)
	{

		for(int64_t count=0; count<COUNT;count++)
		{
			printf("%ld,",count);
			printf("%ld,",seiDynamic[0][count].sameLevel);
			printf("%ld,",seiDynamic[0][count].adjacent);  
			printf("%ld,", seiDynamic[0][count].movement); 
			printf("%9lf, ",(double)(timingStatic[0])); // Min speedup 
			for (int64_t threadCount=0; threadCount<threadArraySize; threadCount++)
			{	
				/*				printf("%ld,%lf,%ld, %ld, %ld, %lf, %lf, %lf\n", 
				 *				ins,timingDynamic[0][ins] ,
				 (double)timingStatic[0]/(double)timingDynamic[0][ins],
				 (double)staticTraverseVerticeCounter[0]/(double)dynamicTraverseVerticeCounter[0][ins],
				 (double)staticTraverseEdgeCounter[0]/(double)dynamicTraverseEdgeCounter[0][ins]);
				 printf("%9lf, ",(double)staticTraverseEdgeCounter[0]/(double)dynamicTraverseEdgeCounterMax[threadCount]);
				 */
				printf("%9lf, ",(double)(timingDynamic[threadCount][count])); // Min speedup 
				printf("%9lf, ",(double)(timingDynamic[threadCount][count])/(double)timingStatic[threadCount]); // Min speedup 
				printf("%9lf",(double)(dynamicTraverseEdgeCounterMax[threadCount][count])/(double)(staticTraverseEdgeCounter[threadCount]) );
				printf("%9lf, ",(double)(dynamicTraverseEdgeCounterMax[threadCount][count])/(double)(staticTraverseEdgeCounter[0]) );

				//				printf("%3ld, %3ld, ",ins,threadArray[threadCount]); // Number of threads used
				//				printf("%9lf, ",(double)timingStatic[threadCount]/(double)(timingDynamic[threadCount][ins])); // Min speedup 
				//				printf("%9lf, ",(double)staticTraverseEdgeCounter[threadCount] /(double)dynamicTraverseEdgeCounterMax[threadCount][ins]);
			}
			printf("\n");
		}
	}
	if(0)
	{
		for (int64_t threadCount=0; threadCount<threadArraySize; threadCount++)
		{   
			int64_t threads = threadArray[threadCount];
			qsort(timingDynamic[threadCount],COUNT,sizeof(double), doubleCompareforMinSort);
			printf("%3ld, ",threads); // Number of threads used
			printf("%9lf, ",timingStatic[threadCount]); // Time of static algorithm
			printf("%9lf, ",timingDynamicTotal[threadCount]/COUNT); // Time of dynamic algorithm
			printf("%9lf, ",1.0/(double)timingDynamicTotal[threadCount]*COUNT); // Updates per second
			printf("%9lf, ",timingStatic[0]/timingStatic[threadCount]); // Speedup of static algorithm
			printf("%9lf, ",timingDynamicTotal[0]/timingDynamicTotal[threadCount]); // Speedup of dynamic algorithm


			printf("%9lf, ",timingStatic[0]/(timingDynamic[threadCount][0])); // Min speedup 
			printf("%9lf, ",timingStatic[0]/(timingDynamic[threadCount][COUNT/4])); // Min speedup 
			printf("%9lf, ",timingStatic[0]/(timingDynamic[threadCount][COUNT/2])); // Min speedup 
			printf("%9lf, ",timingStatic[0]/(timingDynamicTotal[threadCount]/COUNT)); // Average speedup
			printf("%9lf, ",timingStatic[0]/(timingDynamic[threadCount][3*COUNT/4])); // Min speedup 
			printf("%9lf, ",timingStatic[0]/(timingDynamic[threadCount][COUNT-1])); // Max Speedup

			printf("%9lf, ",timingStatic[threadCount]/(timingDynamicTotal[threadCount]/COUNT)); //
			printf("%9lf, ",timingStatic[threadCount]/(timingDynamic[threadCount][0])); // Min speedup 
			printf("%9lf, ",timingStatic[threadCount]/(timingDynamic[threadCount][COUNT/4])); // Min speedup 
			printf("%9lf, ",timingStatic[threadCount]/(timingDynamic[threadCount][COUNT/2])); // Min speedup 
			printf("%9lf, ",timingStatic[threadCount]/(timingDynamicTotal[threadCount]/COUNT)); // Average speedup
			printf("%9lf, ",timingStatic[threadCount]/(timingDynamic[threadCount][3*COUNT/4])); // Min speedup 
			printf("%9lf, ",timingStatic[threadCount]/(timingDynamic[threadCount][COUNT-1])); // Max Speedup

			//printf("%9lf, ",(double)staticTraverseEdgeCounter[0]/(double)dynamicTraverseEdgeCounterTotal[0]);

			/*		printf("%9lf, "); //
					printf("%9lf, "); //
					printf("%9lf, "); //
			/*		printf("%9lf, "); //
			printf("%9lf,"); //
			printf("%9lf,"); //
			 */		//printf("%lf,"); //
			//printf("%lf, %lf, %lf, %lf, %lf, %lf, %lf\n",timeFullBeforeMulti,streamingTimeTotal/COUNT,minTime,maxTime, 
			//		(double)globalSEI.sameLevel/(double)COUNT, (double)globalSEI.adjacent/(double)COUNT, (double)globalSEI.movement/(double)COUNT);
			printf("\n");
		} 

	}
}


void CreateRandomEdgeListFromGraph(struct stinger* stingerGraph, int64_t NV, int64_t* insertionArraySrc,
		int64_t* insertionArrayDest, int64_t insertionCount)
{
	int64_t ins=0,src,dest,srcAdj,destInAdj,destCounter;

	while (ins<insertionCount)
	{
		src = rand() % NV;

		srcAdj=stinger_typed_outdegree(stingerGraph,src,0);

		// No adjacenies
		if(srcAdj==0)
			continue;
		destInAdj=rand()%srcAdj;
		destCounter=0;dest=0;
		STINGER_FORALL_EDGES_OF_VTX_BEGIN(stingerGraph,srcAdj)
                {
			dest=STINGER_EDGE_DEST;
		    if (destInAdj==destCounter)
		    {
			break;
		    }

		    destCounter++;
		}
                STINGER_FORALL_EDGES_OF_VTX_END();

		if(src==dest)
			continue;
		if(src> NV || dest >NV)
			printf("bbbbbbbbbbbbbbbbbaaaaaaaaaaaaa %ld\n",ins);
		stinger_remove_edge(stingerGraph,0,src,dest);
		stinger_remove_edge(stingerGraph,0,dest,src);

		insertionArraySrc[ins]=src;
		insertionArrayDest[ins]=dest;
		//printf("%ld %ld %ld\n",src,dest,destInAdj);  fflush(stdout);
		ins++;
	}
}

void CreateRandomEdgeListFromGraphDeleting(struct stinger* stingerGraph, int64_t NV, int64_t* deletionArraySrc,
                int64_t* deletionArrayDest, int64_t deletionCount)
{
    
    int64_t del = 0, src, dest;

    /*
    stinger_insert_edge(stingerGraph, 0, 6945, 16642, 0, 0);
    stinger_insert_edge(stingerGraph, 0, 16642, 6945, 0, 0);

    deletionArraySrc[0] = 6945;
    deletionArrayDest[0] = 16642;
    */
    while (del < deletionCount)
    {
        src = rand() % NV;
        dest = rand() % NV;
        if (src == dest)
            continue;

        if (src == 0 || dest == 0)
            continue;

        int result = stinger_insert_edge(stingerGraph, 0, src, dest, 0, 0);

        if (result < 1)
            continue;

        //printf("Src: %ld, Dest: %ld\n", src, dest);
        stinger_insert_edge(stingerGraph, 0, dest, src, 0, 0);
        deletionArraySrc[del] = src;
        deletionArrayDest[del] = dest;
        del++;
    }
    
     
    /*
    int64_t del = 0, src, dest, srcAdj, destInAdj, destCounter;

    while (del < deletionCount)
    {
        src = rand() % NV;

        srcAdj = stinger_typed_outdegree(stingerGraph, src, 0);

        // No adjacencies
        if (srcAdj == 0)
            continue;

        destInAdj = rand() % srcAdj;
        destCounter = 0;
        dest = 0;
        STINGER_FORALL_EDGES_OF_VTX_BEGIN(stingerGraph, srcAdj)
            dest = STINGER_EDGE_DEST;
            if (destInAdj == destCounter)
                break;
            destCounter++;
        STINGER_FORALL_EDGES_OF_VTX_END();

        if (src == dest)
            continue;
        if (src > NV || dest > NV)
            printf("Oops %ld\n", del);

        deletionArraySrc[del] = src;
        deletionArrayDest[del] = dest;
        del++;
    }*/
}

void hostParseArgsVitalUpdate(int argc, char** argv, int64_t *NV, int64_t *NE, int64_t *NK, int64_t *NT,
                                        int64_t *randomSeed, int64_t *iterationCount, char *initial_graph_name[1024]) {
        updateType opType;
	static struct option long_options[] = {
		{"VCount", optional_argument, 0, 'V'},
		{"EFactor", optional_argument, 0, 'E'},
		{"RandomSeed", required_argument, 0, 'R'},
		{"Iterations", required_argument, 0, 'G'},
		{"GraphType", optional_argument, 0, 'I'},
		{"GraphName",optional_argument,0,'N'},
		{"OperationType",optional_argument,0,'O'},
		{"Threads",optional_argument,0,'T'},
		{"ApproximateVCount",optional_argument,0,'K'},

		{0, 0, 0, 0}
	};

	while(1) {
		int32_t option_index = 0;
		int32_t c = getopt_long(argc, argv, "V:E:R:G:I:N:O:T:K:?h", long_options, &option_index);
		extern char * optarg;
		extern int32_t    optind, opterr, optopt;
		int64_t intout = 0;

		if(-1 == c)
			break;

		switch(c) {
			default:
				printf("Unrecognized option: %c\n\n", c);
			case '?':
			case 'h':
				exit(0);
				break;

			case 'V':
				errno = 0;
				intout = strtol(optarg, NULL, 10);
				if(errno || intout < 0) {
					printf("Error - Number of vertices needs to be positive %s\n", optarg);
					exit(-1);
				}
				*NV = intout;
				break;
			case 'E':
				errno = 0;
				intout = strtol(optarg, NULL, 10);
				if(errno || intout < 0) {
					printf("Error - Number of edges needs to be positive %s\n", optarg);
					exit(-1);
				}
				*NE = intout;
				break;
			case 'R':
				errno = 0;
				intout = strtol(optarg, NULL, 10);
				if(errno || intout < 0) {
					printf("Error - Seed needs to be positive. R=0 uses a random seed based on the time%s\n", optarg);
					exit(-1);
				}
				*randomSeed = intout;
				break;
			case 'I':
				errno = 0;
				intout = strtol(optarg, NULL, 10);
				if(errno || intout < 0) {
					printf("Error - Number of iterations needs to be positive %s\n", optarg);
					exit(-1);
				}
				iterationCount = intout;
				break;

			case 'G':
				errno = 0;
				intout = strtol(optarg, NULL, 10);
				if(errno || intout < 0 || intout > 5) {
					printf("Error - Limited number of graph types %s\n", optarg);
					exit(-1);
				}
				//graphTestCase = intout;

			case 'N':
				strcpy(initial_graph_name,optarg);
				break;
			case 'O':
				errno = 0;
				intout = strtol(optarg, NULL, 10);
				if(errno || intout < 0 || intout > 1) {
					printf("Error - Operation can be either 0 or 1. 0 for insertion and 1 for deletion%s\n", optarg);
					exit(-1);
				}
				opType = intout;
				break;
			case 'T':
				errno = 0;
				intout = strtol(optarg, NULL, 10);
				if(errno || intout < 0) {
					printf("Error - Number of threads needs to be positive %s\n", optarg);
					exit(-1);
				}
				*NT = intout;
				if(*NT<1){
					*NT=1;
				}

				//printf("the number of threads is %d\n",NT);
				break;
			case 'K':
				errno = 0;
				intout = strtol(optarg, NULL, 10);
				if(errno || intout < 0) {
					printf("Error - Number of approximate vertices needs to be positive %s\n", optarg);
					exit(-1);
				}
				*NK = intout;
				break;
		}
	}
//	printf("Vertix count : %d , Edge count : %d\n", NV,NE);
//	printf("NT : %d\n", NT);
}



//------------------------------------------
//------------------------------------------
// CSR Graph representation.
// Used only for graph ingestion. 


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