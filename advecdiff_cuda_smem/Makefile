COMP=/hpc/sw/cuda/9.0.176/bin/nvcc
#FLAGS= -O3 -arch=sm_35 -lineinfo -I/hpc/sw/cuda/9.0.176/include/
FLAGS= -O3 -arch=sm_35 -I/hpc/sw/cuda/9.0.176/include/

#ifeq ($(COMP), nvcc)
  #FLAGS= -O3 -arch=sm_35 -lineinfo -I/opt/cuda/include  -Xptxas="-v" -maxrregcount=36 # NVCC
  #FLAGS= -O3 -arch=sm_20 -use_fast_math -ftz=true -prec-div=false -prec-sqrt=false -lineinfo -I/opt/cuda/include # NVCC
  #FLAGS= -O3 -arch=sm_20 -use_fast_math -lineinfo -I/opt/cuda/include # NVCC
#endif
#ifeq ($(COMP), g++)
#  FLAGS= -O3 -ffast-math -mtune=native # G++
#endif

all:
	$(COMP) -o advec advec4.cu $(FLAGS)
	$(COMP) -o diff diff4.cu $(FLAGS) 

advec:
	$(COMP) -o advec advec4.cu $(FLAGS)

diff:
	$(COMP) -o diff diff4.cu $(FLAGS)

simple:
	$(COMP) -o simple simple.cu $(FLAGS)


