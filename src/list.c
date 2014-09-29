#include "list.h"

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "xmalloc.h"

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

