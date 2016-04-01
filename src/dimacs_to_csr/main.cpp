 

//#include <iostream>
//#include <numeric>
#include <stdlib.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <inttypes.h>


#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#include "main.h"

using namespace std;

void readGraphDIMACS(char* filePath, int32_t** prmoff, int32_t** prmind, int32_t* prmnv, int32_t* prmne);

int main(const int argc, char *argv[])
{
    cudaSetDevice(0);
	// cudaDeviceProp prop;
	// cudaGetDeviceProperties(&prop, 0);
 //    printf("  Device name: %s\n", prop.name);

    int32_t nv, ne,*off,*adj;

    cout << argv[1] << endl;
    readGraphDIMACS(argv[1],&off,&adj,&nv,&ne);

    cout << "Vertices " << nv << endl;
    cout << "Edges " << ne << endl;

    thrust::device_vector<int32_t> d_adjSizeUsed,d_adjSizeMax,*d_adjArray;

    cout << d_adjSizeUsed.size() << endl;
    cout << d_adjSizeMax.size() << endl;
    allocGPUMemory(nv, ne, off, adj, &d_adjArray, &d_adjSizeUsed, d_adjSizeMax);
    cout << d_adjSizeUsed.size() << endl;
    cout << d_adjSizeMax.size() << endl;


    delete[] d_adjArray;
    free(off);
    free(adj);
    // thrust::host_vector<int32_t> h_off(d_off);
    // thrust::host_vector<int32_t> h_adj(d_adj);

    // cout << h_off.size() << endl;
    // cout << h_adj.size() << endl;


    return 0;
}       

