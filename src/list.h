#pragma once

#include <stdint.h>

typedef struct struct_node
{
    int64_t id;
    int64_t aux_val;
    struct struct_node *next;
} node_t;

typedef struct
{
    node_t *head;
    node_t *tail;
    int64_t size;
} list_t;

typedef list_t* list_ptr;


list_t* makeList(void);
node_t* makeNode(int64_t);
void append(list_t*,node_t*);
node_t *getFirst(list_t*);
void deleteFirst(list_t*);
void printList(list_t*);


void emptyList(list_t*);


void makeArrayOfLists(list_ptr** aL,int64_t numberOfLists);
void destroyArrayOfLists(list_ptr** aL,int64_t numberOfLists);

