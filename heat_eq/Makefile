include_path = /home/anirbankopty/Softwares/FFTW/fftw-install/include
lib_path = /home/anirbankopty/Softwares/FFTW/fftw-install/lib

# SYNTAX
# target: pre-action
#	action

gfortran: gfortran_compile run

# gcc: gcc_compile run

# gcc_compile:
# 	gcc -I${include_path} -L${lib_path} heat_eq.c -o heat_eq.out -lfftw3 -lm

run:
	./heat_eq.out
	# python3 plot.py

gfortran_compile:
	gfortran -L${lib_path} heat_eq_cmplx.f03 -o heat_eq.out -lfftw3

clean:
	rm -rf *.out *.dat