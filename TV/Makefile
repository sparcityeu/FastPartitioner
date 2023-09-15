main: tv.cpp
	g++ ../shared/hypergraph_rw.cpp -c  -w -O3 -std=c++11 #-p -g
	g++ ../shared/metric.cpp -c -O3 -std=c++11 #-p -g
	g++ ../shared/utils.cpp -c -fopenmp -O3 -std=c++11 #-p -g 
	g++ tv.cpp -c -O3 -std=c++11
	g++ partition.cpp -fopenmp -c -std=c++11 -O3 #-p -g 
	nvcc -c gpu_partition.cu #-Xptxas -v -g -G -lineinfo
	nvcc -o tv.exe hypergraph_rw.o tv.o partition.o metric.o utils.o gpu_partition.o -Xcompiler -fopenmp -O3 #-g -G
	rm -f *.o

gdb : tv.cpp
	gcc ../shared/hypergraph_rw.cpp  -c  -w -O3 -pg -p -g
	g++ ../shared/metric.cpp -c -w -O3 -pg -p -g
	g++ ../shared/utils.cpp -c -w -O3 -pg -p -g -fopenmp
	g++ tv.cpp -c -w -O3 -pg -p -g
	g++ partition.cpp -c -w -O3 -pg -p -g -fopenmp
	g++ -o tv.exe hypergraph_rw.o metric.o partition.o tv.o utils.o -O3 -pg -p -g -fopenmp
	rm -f *.o
clean:
	rm -f *.o
	rm -f *.~
	rm -f *.bin
	rm -f *.txt
