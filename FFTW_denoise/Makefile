include_path = /home/anirbankopty/Softwares/FFTW/fftw-install/include
lib_path = /home/anirbankopty/Softwares/FFTW/fftw-install/lib

# SYNTAX
# target: pre-action
#	action

gfortran: gfortran_compile run

gcc: gcc_compile run

gcc_compile:
	gcc -I${include_path} -L${lib_path} FFT_denoise.c -o FFT_denoise.out -lfftw3 -lm

run:
	./FFT_denoise.out
	python3 plot.py

gfortran_compile:
	gfortran -L${lib_path} FFT_denoise.f03 -o FFT_denoise.out -lfftw3

clean:
	rm -rf *.out *.dat