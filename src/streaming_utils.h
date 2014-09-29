#pragma once

#include <stdint.h>

#include "stinger.h"
#include "csr.h"

#define sigma(j,k) pathsToRoot[j][k]
#define sigmaf(j,k) ((bc_t)pathsToRoot[j][k])
#define d(j,k) level[j][k]

#define INFINITY_MY 1073741824

#define VB_ON 0
#define DB_ON 0
#define INFO_ON 0
#define STREAMING_INFO_ON 0
#define COUNTING_INFO_ON 0

//#define PARENT_INFO_ON 0

#define BEFORE_AND_AFTER 1

#if VB_ON
#define VB(...) printf(__VA_ARGS__); fflush(stdout);
#else
#define VB(...)
#endif

#if DB_ON
#define DB(...) printf(__VA_ARGS__); fflush(stdout);
#else
#define DB(...)
#endif

#if INFO_ON
#define INFO(...) printf(__VA_ARGS__); fflush(stdout);
#else
#define INFO(...)
#endif

#if STREAMING_INFO_ON
#define STREAMING_INFO(...) printf(__VA_ARGS__); fflush(stdout);
#else
#define STREAMING_INFO(...)
#endif

#ifdef COUNTING_INFO_ON
extern uint64_t countFullVertex;
extern uint64_t countFullEdgeUp;
extern uint64_t countStreamVertexDown;
extern uint64_t countStreamVertexUp;
extern uint64_t countStreamEdgeDown;
extern uint64_t countStreamEdgeUp;
#endif


typedef double bc_t;

uint64_t ** squareMat(uint64_t NV);
void freeSquareMat(uint64_t ** sq, uint64_t NV);
bc_t ** squareMatFloat(uint64_t NV);
void freeSquareMatFloat(bc_t ** sq, uint64_t NV);
void printGraph(uint64_t ** graph, uint64_t NV);

void intersectAdjacencyWithLevel(uint64_t * outArray, uint64_t * numFound, uint64_t NV,
	uint64_t * adjacency, uint64_t * levelArray, uint64_t level);
void moveUpTree(uint64_t * levels, uint64_t * pathsToRoot, uint64_t ** graph,
	uint64_t vertex, uint64_t ignoreVertex, uint64_t dist, uint64_t NV);
void addEdgeWithoutMovement(uint64_t * levels, uint64_t * pathsToRoot, uint64_t ** graph,
	uint64_t vertex, uint64_t deltaTop, uint64_t NV);
void copyTreeIn(uint64_t * levelsDest, uint64_t * pathsToRootDest, uint64_t existingVertex, uint64_t newVertex,
	uint64_t *levelsSrc, uint64_t * pathsToRootSrc, uint64_t INF, uint64_t NV);


void MatrixBCComput(uint64_t NV, uint64_t * pathsToRoot[],uint64_t * level[],bc_t* perTreeBC[],bc_t* BC);


int64_t serial_shiloach_vishkin_components (struct stinger *S, int64_t nv, int64_t * component_map);
int64_t serial_shiloach_vishkin_componentsCSR (csrGraph* graph, int64_t * component_map);

void selectRootsInMaxComponents(csrGraph* graph,int64_t NV,int64_t* component_map, int64_t componentCount,
                                int64_t numberOfRoots,int64_t* countedEdges, int64_t* selectedRoots);


void read_GML_graph(struct stinger** G, char* filename);
void readSnapGraph(struct stinger** G, char* filename);


void hostParseArgs(int argc, char** argv);
int64_t updateEdge(int argc, char *argv[]);
