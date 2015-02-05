COMP=nvcc

ifeq ($(COMP), nvcc)
  #FLAGS= -O3 -arch=sm_35 -use_fast_math -I/opt/cuda/include # NVCC
  #FLAGS= -O3 -arch=sm_35 -use_fast_math -lineinfo -I/opt/cuda/include # NVCC
  FLAGS= -O3 -arch=sm_20 -use_fast_math -lineinfo -I/opt/cuda/include # NVCC
endif
ifeq ($(COMP), g++)
  FLAGS= -O3 -ffast-math -mtune=native # G++
endif

all:
	$(COMP) -o advec advec4.cu $(FLAGS)
	$(COMP) -o diff diff4.cu $(FLAGS) 

advec:
	$(COMP) -o advec advec4.cu $(FLAGS)

diff:
	$(COMP) -o diff diff4.cu $(FLAGS)

