#include <stdio.h>
#include <stdlib.h>

//#include "omp.h"

int main(void) {
    
    int i = 5;
    //i = 7;
    __sync_bool_compare_and_swap(&i, i, 414);
    printf("i: %ld\n", i);    
    return 0;
}
