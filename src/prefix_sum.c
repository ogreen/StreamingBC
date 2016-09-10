#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "omp.h"

int main(int argc, char **argv) {
    #if 0
    int64_t a[] = {3, 1, 7, 0, 4, 1, 6, 3};
    int64_t b[8];
    int64_t c[8];
    #endif

    int64_t len = 256;
    int64_t a[len];
    int64_t b[len];
    int64_t c[len];
    
    srand(time(NULL));
    for (int64_t i = 0; i < len; i++) {
        a[i] = (rand() % 100) - 50;
    }
    for (int64_t i = 0; i < len; i++) {
        c[i] = a[i];
    }
    
    // sequential.
    //b[0] = a[0];
    b[0] = 0;
    int64_t n = sizeof(a) / sizeof(int64_t);
    printf("n: %d\n", n);
    for (int64_t i = 1; i < n; i++) {
        b[i] = b[i - 1] + a[i - 1];
        //b[i] = b[i - 1] + a[i];
    }  

    // parallel
    int64_t height = (int64_t)(ceil(log2((double)(n)))) - 1;
    
    // up-sweep
    for (int64_t d = 0; d < height; d++) {
        int64_t incr = (int64_t)pow(2, d + 1);
        #pragma omp parallel num_threads(2)
        {
            #pragma omp for
            for(int64_t i = 0; i < n - 1; i += incr) {
                a[i + incr - 1] += a[i + incr / 2 - 1];
            }
        }
    }

    // down-sweep
    a[n - 1] = 0;
    for (int64_t d = height; d >= 0; d--) {
        int64_t incr = (int64_t)pow(2, d + 1);
        #pragma omp parallel num_threads(2)
        {
            #pragma omp for
            for (int64_t i = 0; i < n - 1; i += incr) {
                int64_t t = a[i + incr / 2 - 1];
                a[i + incr / 2 - 1] = a[i + incr - 1];
                a[i + incr - 1] += t;
            }
        }
    }

    printf("a: ");
    for (int64_t i = 0; i < n; i++) {
        printf("%d ", a[i]);
    }
    printf("\n");
    
    printf("b: ");
    for (int64_t i = 0; i < n; i++) {
        printf("%d ", b[i]);
    }
    printf("\n");

    int correct = 1;
    for (int64_t i = 0; i < n; i++) {
        if (a[i] != b[i]) {
            correct = 0;
            break;
        }
    }
    if (correct) {
        printf("correct\n");
    } else {
        printf("incorrect\n");
    }

    return 0;
}
