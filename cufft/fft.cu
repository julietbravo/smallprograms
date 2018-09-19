#include <vector>
#include <iostream>
#include "stdio.h"
#include <cufft.h>

#define cuda_safe_call(err) __cuda_safe_call(err, __FILE__, __LINE__)
inline void __cuda_safe_call(cudaError err, const char *file, const int line)
{
    if (cudaSuccess != err)
        printf("cudaSafeCall() failed at %s:%i : %s\n", file, line, cudaGetErrorString(err));
}

inline int check_cufft(cufftResult err)
{
    if (err == CUFFT_SUCCESS)
    {
        //std::cout << "cufft success!" << std::endl;
        return 0;
    }
    else
    {
        if (err == CUFFT_INVALID_PLAN)
            printf("cuFFT plan error: INVALID PLAN\n");
        else if (err == CUFFT_ALLOC_FAILED)
            printf("cuFFT plan error: ALLOC FAILED\n");
        else if (err == CUFFT_INVALID_TYPE)
            printf("cuFFT plan error: INVALID TYPE\n");
        else if (err == CUFFT_INVALID_VALUE)
            printf("cuFFT plan error: INVALID VALUE\n");
        else if (err == CUFFT_INTERNAL_ERROR)
            printf("cuFFT plan error: INTERNAL ERROR\n");
        else if (err == CUFFT_EXEC_FAILED)
            printf("cuFFT plan error: EXEC FAILED\n");
        else if (err == CUFFT_SETUP_FAILED)
            printf("cuFFT plan error: SETUP FAILED\n");
        else if (err == CUFFT_INVALID_SIZE)
            printf("cuFFT plan error: INVALID SIZE\n");
        else if (err == CUFFT_UNALIGNED_DATA)
            printf("cuFFT plan error: UNALIGNED DATA\n");
        else 
            printf("cuFFT plan error: OTHER\n");

        return 1; 
    }
}

int main()
{
    // Grid et al.
    // ------------------------
    const int nloops = 1000;
    const int itot = 512;
    const int jtot = 512;
    const int ktot = 512;
    const int ncells = itot*jtot*ktot;

    // Field at host
    // ------------------------
    std::vector<double> field(ncells);

    // Create device field & tmp
    // ------------------------
    double* field_g;
    double* tmp_g;
    cuda_safe_call(cudaMalloc((void**)&field_g, ncells*sizeof(double)));
    cuda_safe_call(cudaMalloc((void**)&tmp_g, ncells*sizeof(double)));
    cuda_safe_call(cudaMemcpy(field_g, field.data(), ncells, cudaMemcpyHostToDevice));

    // Create FFT plan
    // ------------------------
    cufftHandle iplanf;
    const int rank = 1;

    // Double input
    int i_ni[]    = {itot};
    int i_istride = 1;
    int i_idist   = itot;

    // Double-complex output
    int o_ni[]    = {itot/2+1};
    int o_istride = 1;
    int o_idist   = itot/2+1;

    check_cufft( cufftPlanMany(&iplanf, rank, i_ni, i_ni, i_istride, i_idist, o_ni, o_istride, o_idist, CUFFT_D2Z, jtot*ktot) ); 

    // Calculate FFTs
    // ------------------------
    for (int i=0; i<nloops; ++i)
    {
        check_cufft( cufftExecD2Z(iplanf, (cufftDoubleReal*)field_g, (cufftDoubleComplex*)tmp_g) );
        cudaThreadSynchronize();
    }

    return 0;
}
