#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <stdlib.h>
#include <cstdio>

__host__ __device__ inline double dg4(const double v1, const double v2, const double v3, const double v4, const double v5, const double v6, const double v7) 
{
    return (1./576.)*(v1+v7) + (-54./576.)*(v2+v6) + (783./576.)*(v3+v5) + (-1460./576.)*v4;
}

/* 
4th order diffusion (2D), similar to CPU implementation MicroHH
*/
void diff_cpu_2d(double * const __restrict__ at, const double * const __restrict__ a,
                 const double dxidxi, const double dyidyi, const double dzidzi,
                 const int istart, const int iend, 
                 const int jstart, const int jend, 
                 const int kstart, const int kend, 
                 const int icells, const int ijcells)
{
    const int ii1 = 1;
    const int ii2 = 2;
    const int ii3 = 3;
    const int jj1 = 1*icells;
    const int jj2 = 2*icells;
    const int jj3 = 3*icells;

    const double visc = 0.1;

    for (int k=kstart; k<kend; ++k)
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*icells + k*ijcells;

                at[ijk] += visc * dg4(a[ijk-ii3], a[ijk-ii2], a[ijk-ii1], a[ijk], a[ijk+ii1], a[ijk+ii2], a[ijk+ii3])*dxidxi
                        +  visc * dg4(a[ijk-jj3], a[ijk-jj2], a[ijk-jj1], a[ijk], a[ijk+jj1], a[ijk+jj2], a[ijk+jj3])*dyidyi;
            }
}

/* 
4th order diffusion (3D), similar to CPU implementation MicroHH
*/
void diff_cpu_3d(double * const __restrict__ at, const double * const __restrict__ a,
                 const double dxidxi, const double dyidyi, const double dzidzi,
                 const int istart, const int iend, 
                 const int jstart, const int jend, 
                 const int kstart, const int kend, 
                 const int icells, const int ijcells)
{
    const int ii1 = 1;
    const int ii2 = 2;
    const int ii3 = 3;
    const int jj1 = 1*icells;
    const int jj2 = 2*icells;
    const int jj3 = 3*icells;
    const int kk1 = 1*ijcells;
    const int kk2 = 2*ijcells;
    const int kk3 = 3*ijcells;

    const double visc = 0.1;

    for (int k=kstart; k<kend; ++k)
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*icells + k*ijcells;

                at[ijk] += visc * dg4(a[ijk-ii3], a[ijk-ii2], a[ijk-ii1], a[ijk], a[ijk+ii1], a[ijk+ii2], a[ijk+ii3])*dxidxi
                        +  visc * dg4(a[ijk-jj3], a[ijk-jj2], a[ijk-jj1], a[ijk], a[ijk+jj1], a[ijk+jj2], a[ijk+jj3])*dyidyi
                        +  visc * dg4(a[ijk-kk3], a[ijk-kk2], a[ijk-kk1], a[ijk], a[ijk+kk1], a[ijk+kk2], a[ijk+kk3])*dzidzi;
            }
}

/* 
4th order diffusion (2D), no shared memory use
*/
__global__ void diff_gpu_2d(double * const __restrict__ at, const double * const __restrict__ a,
                            const double dxidxi, const double dyidyi, const double dzidzi,
                            const int istart, const int iend, 
                            const int jstart, const int jend, 
                            const int kstart, const int kend, 
                            const int icells, const int ijcells)
{
    const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
    const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
    const int k = blockIdx.z + kstart;

    const double visc = 0.1;

    if(i < iend && j < jend && k < kend)
    {
        const int ijk = i + j*icells + k*ijcells;

        const int ii1 = 1;
        const int ii2 = 2;
        const int ii3 = 3;
        const int jj1 = 1*icells;
        const int jj2 = 2*icells;
        const int jj3 = 3*icells;

	at[ijk] += visc * dg4(a[ijk-ii3], a[ijk-ii2], a[ijk-ii1], a[ijk], a[ijk+ii1], a[ijk+ii2], a[ijk+ii3])*dxidxi
	        +  visc * dg4(a[ijk-jj3], a[ijk-jj2], a[ijk-jj1], a[ijk], a[ijk+jj1], a[ijk+jj2], a[ijk+jj3])*dyidyi;
    }
}

/* 
4th order diffusion (3D), no shared memory use
*/
__global__ void diff_gpu_3d(double * const __restrict__ at, const double * const __restrict__ a,
                            const double dxidxi, const double dyidyi, const double dzidzi,
                            const int istart, const int iend, 
                            const int jstart, const int jend, 
                            const int kstart, const int kend, 
                            const int icells, const int ijcells)
{
    const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
    const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
    const int k = blockIdx.z + kstart;

    const double visc = 0.1;

    if(i < iend && j < jend && k < kend)
    {
        const int ijk = i + j*icells + k*ijcells;

        const int ii1 = 1;
        const int ii2 = 2;
        const int ii3 = 3;
        const int jj1 = 1*icells;
        const int jj2 = 2*icells;
        const int jj3 = 3*icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;
        const int kk3 = 3*ijcells;

	at[ijk] += visc * dg4(a[ijk-ii3], a[ijk-ii2], a[ijk-ii1], a[ijk], a[ijk+ii1], a[ijk+ii2], a[ijk+ii3])*dxidxi
	        +  visc * dg4(a[ijk-jj3], a[ijk-jj2], a[ijk-jj1], a[ijk], a[ijk+jj1], a[ijk+jj2], a[ijk+jj3])*dyidyi
	        +  visc * dg4(a[ijk-kk3], a[ijk-kk2], a[ijk-kk1], a[ijk], a[ijk+kk1], a[ijk+kk2], a[ijk+kk3])*dzidzi;
    }
}

/* 
4th order diffusion, 2D smem tile
*/
__global__ void diff_gpu_2d_s2d(double * const __restrict__ at, const double * const __restrict__ a,
                                const double dxidxi, const double dyidyi, const double dzidzi,
                                const int istart, const int iend, 
                                const int jstart, const int jend, 
                                const int kstart, const int kend, 
                                const int icells, const int ijcells, const int ngc)
{
    extern __shared__ double as[]; 

    const int tx = threadIdx.x;
    const int ty = threadIdx.y;
    const int i  = blockIdx.x*blockDim.x + threadIdx.x + istart;
    const int j  = blockIdx.y*blockDim.y + threadIdx.y + jstart;
    const int k  = blockIdx.z + kstart;
    const int blockxpad = blockDim.x+2*ngc;

    const double visc = 0.1;

    if(i < iend && j < jend && k < kend)
    {
        const int ijk  = i + j*icells + k*ijcells; // index in global memory
        const int ijks = (tx+ngc) + (ty+ngc)*blockxpad; // Same location in 2d shared mem

        const int ii1 = 1;
        const int ii2 = 2;
        const int ii3 = 3;
        const int jj3 = 3*icells;

        const int jjs1 = 1*blockxpad;
        const int jjs2 = 2*blockxpad;
        const int jjs3 = 3*blockxpad;

        as[ijks] = a[ijk];

        if(ty < ngc)
            as[ijks-jjs3] = a[ijk-jj3];
        if(ty >= blockDim.y-ngc)
            as[ijks+jjs3] = a[ijk+jj3];

        if(tx < ngc)
            as[ijks-ii3] = a[ijk-ii3];
        if(tx >= blockDim.x-ngc)
            as[ijks+ii3] = a[ijk+ii3];

        __syncthreads();

	at[ijk] += visc * dg4(as[ijks-ii3 ], as[ijks-ii2 ], as[ijks-ii1 ], as[ijks], as[ijks+ii1 ], as[ijks+ii2 ], as[ijks+ii3 ])*dxidxi
	        +  visc * dg4(as[ijks-jjs3], as[ijks-jjs2], as[ijks-jjs1], as[ijks], as[ijks+jjs1], as[ijks+jjs2], as[ijks+jjs3])*dyidyi;
    }
}

/* 
4th order diffusion (3d), 2D smem tile
*/
//__global__ void diff_gpu_3d_s2d(double * const __restrict__ at, const double * const __restrict__ a,
//                                const double dxidxi, const double dyidyi, const double dzidzi,
//                                const int istart, const int iend, 
//                                const int jstart, const int jend, 
//                                const int kstart, const int kend, 
//                                const int icells, const int ijcells, const int ngc)
//{
//    extern __shared__ double as[]; 
//
//    const int tx = threadIdx.x;
//    const int ty = threadIdx.y;
//    const int i  = blockIdx.x*blockDim.x + threadIdx.x + istart;
//    const int j  = blockIdx.y*blockDim.y + threadIdx.y + jstart;
//    const int k  = blockIdx.z + kstart;
//    const int blockxpad = blockDim.x+2*ngc;
//
//    const double visc = 0.1;
//
//    if(i < iend && j < jend && k < kend)
//    {
//        const int ijk  = i + j*icells + k*ijcells; // index in global memory
//        const int ijks = (tx+ngc) + (ty+ngc)*blockxpad; // Same location in 2d shared mem
//
//        const int ii1 = 1;
//        const int ii2 = 2;
//        const int ii3 = 3;
//        const int jj3 = 3*icells;
//        const int kk1 = 1*ijcells;
//        const int kk2 = 2*ijcells;
//        const int kk3 = 3*ijcells;
//
//        const int jjs1 = 1*blockxpad;
//        const int jjs2 = 2*blockxpad;
//        const int jjs3 = 3*blockxpad;
//
//        as[ijks] = a[ijk];
//
//        if(ty < ngc)
//            as[ijks-jjs3] = a[ijk-jj3];
//        if(ty >= blockDim.y-ngc)
//            as[ijks+jjs3] = a[ijk+jj3];
//
//        if(tx < ngc)
//            as[ijks-ii3] = a[ijk-ii3];
//        if(tx >= blockDim.x-ngc)
//            as[ijks+ii3] = a[ijk+ii3];
//
//        __syncthreads();
//
//	at[ijk] += visc * dg4(as[ijks-ii3 ], as[ijks-ii2 ], as[ijks-ii1 ], as[ijks], as[ijks+ii1 ], as[ijks+ii2 ], as[ijks+ii3 ])*dxidxi
//	        +  visc * dg4(as[ijks-jjs3], as[ijks-jjs2], as[ijks-jjs1], as[ijks], as[ijks+jjs1], as[ijks+jjs2], as[ijks+jjs3])*dyidyi
//	        +  visc * dg4(a [ijk-kk3],    a[ijk-kk2],   a [ijk-kk1],   as[ijks], a [ijk+kk1],   a [ijk+kk2],   a [ijk+kk3])*dzidzi;
//    }
//}

/* 
4th order diffusion (3d), 2D smem tile
*/
__global__ void diff_gpu_3d_s2d(double * const __restrict__ at, const double * const __restrict__ a,
                                const double dxidxi, const double dyidyi, const double dzidzi,
                                const int istart, const int iend, 
                                const int jstart, const int jend, 
                                const int kstart, const int kend, 
                                const int icells, const int ijcells, const int ngc)
{
    extern __shared__ double as[]; 

    const int tx = threadIdx.x;
    const int ty = threadIdx.y;
    const int i  = blockIdx.x*blockDim.x + threadIdx.x + istart;
    const int j  = blockIdx.y*blockDim.y + threadIdx.y + jstart;
    //const int k  = blockIdx.z + kstart;
    const int blockxpad = blockDim.x+2*ngc;

    const double visc = 0.1;

    if(i < iend && j < jend)
    {
        const int ii1 = 1;
        const int ii2 = 2;
        const int ii3 = 3;
        const int jj3 = 3*icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;
        const int kk3 = 3*ijcells;

        const int jjs1 = 1*blockxpad;
        const int jjs2 = 2*blockxpad;
        const int jjs3 = 3*blockxpad;

        for(int k=kstart; k<kend; ++k)
        {
            const int ijk  = i + j*icells + k*ijcells; // index in global memory
            const int ijks = (tx+ngc) + (ty+ngc)*blockxpad; // Same location in 2d shared mem

            as[ijks] = a[ijk];

            if(ty < ngc)
                as[ijks-jjs3] = a[ijk-jj3];
            if(ty >= blockDim.y-ngc)
                as[ijks+jjs3] = a[ijk+jj3];

            if(tx < ngc)
                as[ijks-ii3] = a[ijk-ii3];
            if(tx >= blockDim.x-ngc)
                as[ijks+ii3] = a[ijk+ii3];

            __syncthreads();

	    at[ijk] += visc * dg4(as[ijks-ii3 ], as[ijks-ii2 ], as[ijks-ii1 ], as[ijks], as[ijks+ii1 ], as[ijks+ii2 ], as[ijks+ii3 ])*dxidxi
	            +  visc * dg4(as[ijks-jjs3], as[ijks-jjs2], as[ijks-jjs1], as[ijks], as[ijks+jjs1], as[ijks+jjs2], as[ijks+jjs3])*dyidyi
	            +  visc * dg4(a [ijk-kk3],    a[ijk-kk2],   a [ijk-kk1],   as[ijks], a [ijk+kk1],   a [ijk+kk2],   a [ijk+kk3])*dzidzi;
        }
    }
}


/* 
Get max difference between two fields
*/
double maxdiff(const double * const __restrict__ a, const double * const __restrict__ b, const int n)
{
    double maxdiff=0;
    double diff=0;
    for(int i=0; i<n; ++i)
    {
        diff = std::abs(a[i]-b[i]);
        if(diff > maxdiff)
            maxdiff = diff;
    }
    return maxdiff;
}

int main()
{
    //
    // Grid
    //
    const double dxi = 0.1;
    const double dyi = 0.1;
    const double dzi = 0.1;
    
    const int itot = 256;
    const int jtot = 256;
    const int ktot = 256;
    const int gc   = 3;
    const int iter = 40;
    
    //
    // Calculate the required variables.
    //
    const int ncells  = (itot+2*gc)*(jtot+2*gc)*(ktot+2*gc);
    const int istart  = gc;
    const int jstart  = gc;
    const int kstart  = gc;
    const int iend    = itot+gc;
    const int jend    = jtot+gc;
    const int kend    = ktot+gc;
    const int icells  = itot+2*gc;
    const int jcells  = jtot+2*gc;
    const int kcells  = ktot+2*gc;
    const int ijcells = (itot+2*gc)*(jtot+2*gc);
    
    // Padded settings, interior aligned to 128 byte blocks
    const int mo        = 16 - gc;           // Padding at start of array 
    const int pl        = 16-(int)itot%16;   // Elements left in last 128 byte block
    const int icellsp   = itot + pl + (pl < 2*gc)*16;
    const int ijcellsp  = icellsp * jcells;  
    const int ncellsp   = ijcellsp * kcells + mo;
    
    //
    // Prepare fields on HOST
    //
    double *a    = new double[ncells];
    double *at   = new double[ncells];
    double *tmp1 = new double[ncells];
    
    for (int n=0; n<ncells; ++n)
    {
    	a [n]   = 0.001 * (std::rand() % 1000) - 0.5;
    	at[n]   = 0.;
    	tmp1[n] = 0.;
    }
    
    // 
    // Prepare fields on DEVICE
    //
    double *ad, *atd;
    cudaMalloc((void **)&ad,  ncellsp*sizeof(double));
    cudaMalloc((void **)&atd, ncellsp*sizeof(double));
    cudaMemcpy2D(&ad[mo],  icellsp*sizeof(double),  a,  icells*sizeof(double), icells*sizeof(double), jcells*kcells, cudaMemcpyHostToDevice);
    cudaMemcpy2D(&atd[mo], icellsp*sizeof(double),  at, icells*sizeof(double), icells*sizeof(double), jcells*kcells, cudaMemcpyHostToDevice);
    
    //
    // CUDA thread blocks
    //
    const int blocki = 32;
    const int blockj = 16;
    const int gridi  = itot/blocki + (itot%blocki > 0);
    const int gridj  = jtot/blockj + (jtot%blockj > 0);
    dim3 gridGPU  (gridi, gridj, ktot);
    dim3 gridGPU2d(gridi, gridj, 1);
    dim3 blockGPU(blocki, blockj, 1);

    // 
    // Timer stuff
    //
    cudaEvent_t startEvent, stopEvent;
    cudaEventCreate(&startEvent);
    cudaEventCreate(&stopEvent);
    float  dt1, dt2;

    //
    // Execute kernels
    //
    //////////////////// CPU //////////////////////////
    cudaEventRecord(startEvent, 0);
    for(int n=0; n<iter; ++n) // iter+1 since GPU version is warmed up with one call
    {
       //diff_cpu_2d(at,  a,  dxi, dyi, dzi, istart, iend, jstart, jend, kstart, kend, icells, ijcells);
       diff_cpu_3d(at,  a,  dxi, dyi, dzi, istart, iend, jstart, jend, kstart, kend, icells, ijcells);
    }
    cudaEventRecord(stopEvent, 0);
    cudaEventSynchronize(stopEvent);
    cudaEventElapsedTime(&dt1, startEvent, stopEvent);
    printf("CPU; elapsed=%f [ms]\n",dt1);
 
    //////////////////// GPU //////////////////////////
    cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);

    cudaEventRecord(startEvent, 0);
    for(int n=0; n<iter; ++n)
    {
        //diff_gpu_2d<<<gridGPU, blockGPU>>> 
        //         (&atd[mo], &ad[mo], dxi, dyi, dzi, istart, iend, jstart, jend, kstart, kend, icellsp, ijcellsp);

        //diff_gpu_3d<<<gridGPU, blockGPU>>> 
        //         (&atd[mo], &ad[mo], dxi, dyi, dzi, istart, iend, jstart, jend, kstart, kend, icellsp, ijcellsp);

        //diff_gpu_2d_s2d<<<gridGPU, blockGPU, (blocki+2*gc)*(blockj+2*gc)*sizeof(double)>>> 
        //         (&atd[mo], &ad[mo], dxi, dyi, dzi, istart, iend, jstart, jend, kstart, kend, icellsp, ijcellsp, gc);

        diff_gpu_3d_s2d<<<gridGPU2d, blockGPU, (blocki+2*gc)*(blockj+2*gc)*sizeof(double)>>> 
                 (&atd[mo], &ad[mo], dxi, dyi, dzi, istart, iend, jstart, jend, kstart, kend, icellsp, ijcellsp, gc);
    }
    cudaEventRecord(stopEvent, 0);
    cudaEventSynchronize(stopEvent);
    cudaEventElapsedTime(&dt2, startEvent, stopEvent);

    //
    // Copy device field to tmp1 
    //
    cudaMemcpy2D(tmp1, icells*sizeof(double), &atd[mo], icellsp*sizeof(double), icells*sizeof(double), jcells*kcells, cudaMemcpyDeviceToHost);

    printf("GPU; elapsed=%f [ms], speedup=%f, maxdiff=%e \n",dt2,dt1/dt2,maxdiff(at,tmp1,ncells));

    return 0;
}
