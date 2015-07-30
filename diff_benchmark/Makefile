COMP=g++-mp-4.9
COMPF=gfortran-mp-4.9
#COMP=g++ #=clang on macbook
FLAGS=-Ofast -mtune=native

all:
	$(COMP) -o diff   diff.cxx $(FLAGS)
	$(COMPF) -o diff_f diff.f90 $(FLAGS)


