include_path = /home/anirbankopty/Softwares/FFTW/fftw-install/include
lib_path = /home/anirbankopty/Softwares/FFTW/fftw-install/lib

# SYNTAX
# target: pre-action
#	action

gfortran: gfortran_compile run

gcc: gcc_compile run

gcc_compile:
	gcc -I${include_path} -L${lib_path} FFTW_denoise.c -o FFTW_denoise.out -lfftw3 -lm

run:
	./FFTW_denoise.out
	python3 plot.py

gfortran_compile:
	gfortran -L${lib_path} FFTW_denoise.f03 -o FFTW_denoise.out -lfftw3

clean:
	rm -rf *.out *.dat