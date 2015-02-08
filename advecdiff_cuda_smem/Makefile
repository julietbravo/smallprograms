COMP=/opt/cuda/bin/nvcc

#ifeq ($(COMP), nvcc)
  FLAGS= -O3 -arch=sm_35 -lineinfo -I/opt/cuda/include  -Xptxas="-v" # NVCC
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


