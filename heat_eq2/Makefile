.DEFAULT_GOAL := gfortran  

include_path = /home/anirbankopty/Softwares/FFTW/fftw-install/include
lib_path = /home/anirbankopty/Softwares/FFTW/fftw-install/lib

# SYNTAX
# target: pre-action
#	action

gcc: gcc_compile run

g+: g++_compile run

gfortran: gfortran_compile run

gcc_compile:
	gcc -I${include_path} -L${lib_path} heat_eq.c -o heat_eq.out -lfftw3 -lm

g++_compile:
	g++ -I${include_path} -L${lib_path} heat_eq.cpp -o heat_eq.out -lfftw3 -lm

gfortran_compile:
	gfortran -L${lib_path} heat_eq.f03 -o heat_eq.out -lfftw3

run:
	./heat_eq.out
	# python3 plot.py

clean:
	rm -rf *.out *.dat