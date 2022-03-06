.DEFAULT_GOAL := gfortran  

include_path = /home/anirbankopty/Softwares/FFTW/fftw-install/include
lib_path = /home/anirbankopty/Softwares/FFTW/fftw-install/lib
LIBS = -lfftw3 -lm

# The files to be executed
SOURCE = heat_eq
TARGET = heat_eq.out

# SYNTAX
# target: pre-action
#	action

gcc: gcc_compile run

g+: g++_compile run

gfortran: gfortran_compile run

gcc_compile:
	gcc -I${include_path} -L${lib_path} $(SOURCE).c -o $(TARGET) $(LIBS)

g++_compile:
	g++ -I${include_path} -L${lib_path} $(SOURCE).cpp -o $(TARGET) $(LIBS)

gfortran_compile:
	gfortran -L${lib_path} $(SOURCE).f03 -o $(TARGET) $(LIBS)

run:
	./heat_eq.out
	python3 heatmap.py &

clean:
	rm -rf *.out *.dat

cleanall: clean
	rm -rf *.png
