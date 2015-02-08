#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <stdlib.h>
#include <cstdio>

// Fourth order interpolation function.
__host__ __device__ inline double interp(const double m2, const double m1, const double p1, const double p2)
{
  return (-1./16)*(m2+p2) + (9./16)*(m1+p1);
}

// Fourth order gradient function.
__host__ __device__ inline double grad(const double m2, const double m1, const double p1, const double p2)
{
  return (1./24.)*(m2-p2) + (27./24.)*(p1-m1);
}

/* 
4th order advection on cpu
*/
void advec_cpu(double * const __restrict__ ut, 
               const double * const __restrict__ u, const double * const __restrict__ v, const double * const __restrict__ w,
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

    for (int k=kstart; k<kend; ++k)
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*icells + k*ijcells;

                ut[ijk] += grad( interp( u[ijk-ii3], u[ijk-ii2], u[ijk-ii1], u[ijk    ] ) * interp( u[ijk-ii3], u[ijk-ii2], u[ijk-ii1], u[ijk    ] ),
                                 interp( u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1] ) * interp( u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1] ),
                                 interp( u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2] ) * interp( u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2] ),
                                 interp( u[ijk    ], u[ijk+ii1], u[ijk+ii2], u[ijk+ii3] ) * interp( u[ijk    ], u[ijk+ii1], u[ijk+ii2], u[ijk+ii3] ))

                         + grad( interp( v[ijk-ii2-jj1], v[ijk-ii1-jj1], v[ijk-jj1], v[ijk+ii1-jj1] ) * interp( u[ijk-jj3], u[ijk-jj2], u[ijk-jj1], u[ijk    ] ),
                                 interp( v[ijk-ii2    ], v[ijk-ii1    ], v[ijk    ], v[ijk+ii1    ] ) * interp( u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1] ),
                                 interp( v[ijk-ii2+jj1], v[ijk-ii1+jj1], v[ijk+jj1], v[ijk+ii1+jj1] ) * interp( u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2] ),
                                 interp( v[ijk-ii2+jj2], v[ijk-ii1+jj2], v[ijk+jj2], v[ijk+ii1+jj2] ) * interp( u[ijk    ], u[ijk+jj1], u[ijk+jj2], u[ijk+jj3] ))

                         + grad( interp( w[ijk-ii2-kk1], w[ijk-ii1-kk1], w[ijk-kk1], w[ijk+ii1-kk1] ) * interp( u[ijk-kk3], u[ijk-kk2], u[ijk-kk1], u[ijk    ] ),
                                 interp( w[ijk-ii2    ], w[ijk-ii1    ], w[ijk    ], w[ijk+ii1    ] ) * interp( u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1] ),
                                 interp( w[ijk-ii2+kk1], w[ijk-ii1+kk1], w[ijk+kk1], w[ijk+ii1+kk1] ) * interp( u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2] ),
                                 interp( w[ijk-ii2+kk2], w[ijk-ii1+kk2], w[ijk+kk2], w[ijk+ii1+kk2] ) * interp( u[ijk    ], u[ijk+kk1], u[ijk+kk2], u[ijk+kk3] ));
            }
}


/* 
4th order advection (3D), no shared memory use
*/
__global__ void advec_gpu(double * const __restrict__ ut, 
                          const double * const __restrict__ u, const double * const __restrict__ v, const double * const __restrict__ w,
                          const int istart, const int iend, 
                          const int jstart, const int jend, 
                          const int kstart, const int kend, 
                          const int icells, const int ijcells)
{
    const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
    const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
    const int k = blockIdx.z + kstart;

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

        ut[ijk] += grad( interp( u[ijk-ii3], u[ijk-ii2], u[ijk-ii1], u[ijk    ] ) * interp( u[ijk-ii3], u[ijk-ii2], u[ijk-ii1], u[ijk    ] ),
                         interp( u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1] ) * interp( u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1] ),
                         interp( u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2] ) * interp( u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2] ),
                         interp( u[ijk    ], u[ijk+ii1], u[ijk+ii2], u[ijk+ii3] ) * interp( u[ijk    ], u[ijk+ii1], u[ijk+ii2], u[ijk+ii3] ))

                 + grad( interp( v[ijk-ii2-jj1], v[ijk-ii1-jj1], v[ijk-jj1], v[ijk+ii1-jj1] ) * interp( u[ijk-jj3], u[ijk-jj2], u[ijk-jj1], u[ijk    ] ),
                         interp( v[ijk-ii2    ], v[ijk-ii1    ], v[ijk    ], v[ijk+ii1    ] ) * interp( u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1] ),
                         interp( v[ijk-ii2+jj1], v[ijk-ii1+jj1], v[ijk+jj1], v[ijk+ii1+jj1] ) * interp( u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2] ),
                         interp( v[ijk-ii2+jj2], v[ijk-ii1+jj2], v[ijk+jj2], v[ijk+ii1+jj2] ) * interp( u[ijk    ], u[ijk+jj1], u[ijk+jj2], u[ijk+jj3] ))

                 + grad( interp( w[ijk-ii2-kk1], w[ijk-ii1-kk1], w[ijk-kk1], w[ijk+ii1-kk1] ) * interp( u[ijk-kk3], u[ijk-kk2], u[ijk-kk1], u[ijk    ] ),
                         interp( w[ijk-ii2    ], w[ijk-ii1    ], w[ijk    ], w[ijk+ii1    ] ) * interp( u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1] ),
                         interp( w[ijk-ii2+kk1], w[ijk-ii1+kk1], w[ijk+kk1], w[ijk+ii1+kk1] ) * interp( u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2] ),
                         interp( w[ijk-ii2+kk2], w[ijk-ii1+kk2], w[ijk+kk2], w[ijk+ii1+kk2] ) * interp( u[ijk    ], u[ijk+kk1], u[ijk+kk2], u[ijk+kk3] ));
    }
}

/* 
4th order advection, smem
*/
__global__ void advec_gpu_smem(double * const __restrict__ ut, 
                               const double * const __restrict__ u, const double * const __restrict__ v, const double * const __restrict__ w, 
                               const int istart, const int iend, 
                               const int jstart, const int jend, 
                               const int kstart, const int kend, 
                               const int icells, const int ijcells, const int ngc)
{
    extern __shared__ double shared[]; 
    const int smem_block = (blockDim.x + 2*ngc) * (blockDim.y + 2*ngc);
    double *us = &shared[0];
    double *vs = &shared[smem_block];

    const int tx = threadIdx.x;
    const int ty = threadIdx.y;
    const int i  = blockIdx.x*blockDim.x + threadIdx.x + istart;
    const int j  = blockIdx.y*blockDim.y + threadIdx.y + jstart;
    const int blockxpad = blockDim.x+2*ngc;

    if(i < iend && j < jend)
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
        const int kk4 = 4*ijcells;

        const int jjs1 = 1*blockxpad;
        const int jjs2 = 2*blockxpad;
        const int jjs3 = 3*blockxpad;

        int ijk;
        const int ijks = (tx+ngc) + (ty+ngc)*blockxpad;

        for(int k=kstart; k<kend; ++k)
        {
            ijk  = i + j*icells + k*ijcells; // index in global memory
 
            us[ijks] = u[ijk];
            vs[ijks] = v[ijk];

            if(ty < ngc)
            {
                us[ijks-jjs3] = u[ijk-jj3];
                vs[ijks-jjs3] = v[ijk-jj3];
            }
            if(ty >= blockDim.y-ngc)
            {
                us[ijks+jjs3] = u[ijk+jj3];
                vs[ijks+jjs3] = v[ijk+jj3];
            }

            if(tx < ngc)
            {
                us[ijks-ngc] = u[ijk-ngc];
                vs[ijks-ngc] = v[ijk-ngc];

                if(ty < ngc)
                {
                    us[ijks-jjs3-ngc] = u[ijk-jj3-ngc];
                    vs[ijks-jjs3-ngc] = v[ijk-jj3-ngc];
                }
                if(ty >= blockDim.y-ngc)
                {
                    us[ijks+jjs3-ngc] = u[ijk+jj3-ngc];
                    vs[ijks+jjs3-ngc] = v[ijk+jj3-ngc];
                }
            }
            if(tx >= blockDim.x-ngc)
            {
                us[ijks+ngc] = u[ijk+ngc];
                vs[ijks+ngc] = v[ijk+ngc];

                if(ty < ngc)
                {
                    us[ijks-jjs3+ngc] = u[ijk-jj3+ngc];
                    vs[ijks-jjs3+ngc] = v[ijk-jj3+ngc];
                }
                if(ty >= blockDim.y-ngc)
                {
                    us[ijks+jjs3+ngc] = u[ijk+jj3+ngc];
                    vs[ijks+jjs3+ngc] = v[ijk+jj3+ngc];
                }
            }

            __syncthreads();

            ut[ijk] += grad( interp( us[ijks-ii3], us[ijks-ii2], us[ijks-ii1], us[ijks    ] ) * interp( us[ijks-ii3], us[ijks-ii2], us[ijks-ii1], us[ijks    ] ),
                             interp( us[ijks-ii2], us[ijks-ii1], us[ijks    ], us[ijks+ii1] ) * interp( us[ijks-ii2], us[ijks-ii1], us[ijks    ], us[ijks+ii1] ),
                             interp( us[ijks-ii1], us[ijks    ], us[ijks+ii1], us[ijks+ii2] ) * interp( us[ijks-ii1], us[ijks    ], us[ijks+ii1], us[ijks+ii2] ),
                             interp( us[ijks    ], us[ijks+ii1], us[ijks+ii2], us[ijks+ii3] ) * interp( us[ijks    ], us[ijks+ii1], us[ijks+ii2], us[ijks+ii3] ))

                     + grad( interp( vs[ijks-ii2-jjs1], vs[ijks-ii1-jjs1], vs[ijks-jjs1], vs[ijks+ii1-jjs1] ) * interp( us[ijks-jjs3], us[ijks-jjs2], us[ijks-jjs1], us[ijks     ] ),
                             interp( vs[ijks-ii2    ],  vs[ijks-ii1     ], vs[ijks    ],  vs[ijks+ii1     ] ) * interp( us[ijks-jjs2], us[ijks-jjs1], us[ijks     ], us[ijks+jjs1] ),
                             interp( vs[ijks-ii2+jjs1], vs[ijks-ii1+jjs1], vs[ijks+jjs1], vs[ijks+ii1+jjs1] ) * interp( us[ijks-jjs1], us[ijks     ], us[ijks+jjs1], us[ijks+jjs2] ),
                             interp( vs[ijks-ii2+jjs2], vs[ijks-ii1+jjs2], vs[ijks+jjs2], vs[ijks+ii1+jjs2] ) * interp( us[ijks     ], us[ijks+jjs1], us[ijks+jjs2], us[ijks+jjs3] ))

                     + grad( interp( w[ijk-ii2-kk1], w[ijk-ii1-kk1], w[ijk-kk1], w[ijk+ii1-kk1] ) * interp( u[ijk-kk3], u[ijk-kk2], u[ijk-kk1], us[ijks  ] ),
                             interp( w[ijk-ii2    ], w[ijk-ii1    ], w[ijk    ], w[ijk+ii1    ] ) * interp( u[ijk-kk2], u[ijk-kk1], us[ijks  ], u[ijk+kk1] ),
                             interp( w[ijk-ii2+kk1], w[ijk-ii1+kk1], w[ijk+kk1], w[ijk+ii1+kk1] ) * interp( u[ijk-kk1], us[ijks  ], u[ijk+kk1], u[ijk+kk2] ),
                             interp( w[ijk-ii2+kk2], w[ijk-ii1+kk2], w[ijk+kk2], w[ijk+ii1+kk2] ) * interp( us[ijks  ], u[ijk+kk1], u[ijk+kk2], u[ijk+kk3] ));
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
    const int itot = 256;
    const int jtot = 256;
    const int ktot = 256;
    const int gc   = 3;
    const int iter = 10;
    
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
    double *u    = new double[ncells];
    double *v    = new double[ncells];
    double *w    = new double[ncells];
    double *ut   = new double[ncells];
    double *tmp1 = new double[ncells];
    
    for (int n=0; n<ncells; ++n)
    {
    	u [n]   = 0.001 * (std::rand() % 1000) - 0.5;
    	v [n]   = 0.001 * (std::rand() % 1000) - 0.5;
    	w [n]   = 0.001 * (std::rand() % 1000) - 0.5;
    	ut[n]   = 0.;
    	tmp1[n] = 0.;
    }
    
    // 
    // Prepare fields on DEVICE
    //
    double *ud, *vd, *wd, *utd;
    cudaMalloc((void **)&ud,  ncellsp*sizeof(double));
    cudaMalloc((void **)&vd,  ncellsp*sizeof(double));
    cudaMalloc((void **)&wd,  ncellsp*sizeof(double));
    cudaMalloc((void **)&utd, ncellsp*sizeof(double));
    cudaMemcpy2D(&ud[mo],  icellsp*sizeof(double),  u,  icells*sizeof(double), icells*sizeof(double), jcells*kcells, cudaMemcpyHostToDevice);
    cudaMemcpy2D(&vd[mo],  icellsp*sizeof(double),  v,  icells*sizeof(double), icells*sizeof(double), jcells*kcells, cudaMemcpyHostToDevice);
    cudaMemcpy2D(&wd[mo],  icellsp*sizeof(double),  w,  icells*sizeof(double), icells*sizeof(double), jcells*kcells, cudaMemcpyHostToDevice);
    cudaMemcpy2D(&utd[mo], icellsp*sizeof(double),  ut, icells*sizeof(double), icells*sizeof(double), jcells*kcells, cudaMemcpyHostToDevice);
    
    //
    // CUDA thread blocks
    //
    const int blocki = 32;
    const int blockj = 8;
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
       advec_cpu(ut, u, v, w,istart, iend, jstart, jend, kstart, kend, icells, ijcells);
    }
    cudaEventRecord(stopEvent, 0);
    cudaEventSynchronize(stopEvent);
    cudaEventElapsedTime(&dt1, startEvent, stopEvent);
    printf("CPU; elapsed=%f [ms]\n",dt1);
 
    ////////////////////// GPU //////////////////////////
    cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);

    cudaEventRecord(startEvent, 0);
    for(int n=0; n<iter; ++n)
    {
        //advec_gpu<<<gridGPU, blockGPU>>> 
        //         (&utd[mo], &ud[mo], &vd[mo], &wd[mo], istart, iend, jstart, jend, kstart, kend, icellsp, ijcellsp);
        advec_gpu_smem<<<gridGPU2d, blockGPU, 2*(blocki+2*gc)*(blockj+2*gc)*sizeof(double)>>> 
                 (&utd[mo], &ud[mo], &vd[mo], &wd[mo], istart, iend, jstart, jend, kstart, kend, icellsp, ijcellsp, gc);
    }
    cudaEventRecord(stopEvent, 0);
    cudaEventSynchronize(stopEvent);
    cudaEventElapsedTime(&dt2, startEvent, stopEvent);

    //
    // Copy device field to tmp1 
    //
    cudaMemcpy2D(tmp1, icells*sizeof(double), &utd[mo], icellsp*sizeof(double), icells*sizeof(double), jcells*kcells, cudaMemcpyDeviceToHost);

    printf("GPU; elapsed=%f [ms], speedup=%f, maxdiff=%e \n",dt2,dt1/dt2,maxdiff(ut,tmp1,ncells));

    return 0;
}
