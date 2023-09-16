main: mgav.cpp
	g++ ../shared/hypergraph_rw.cpp -c  -w -O3 -std=c++11
	g++ ../shared/metric.cpp -c -O3 -std=c++11
	g++ ../shared/utils.cpp -c -fopenmp -O3 -std=c++11
	g++ mgav.cpp -c -O3 -std=c++11
	g++ partition.cpp -fopenmp -c -O3 -std=c++11
	nvcc -c gpu_partition.cu
	nvcc -o mgav.exe hypergraph_rw.o mgav.o partition.o metric.o utils.o -Xcompiler -fopenmp -O3 -std=c++11 -lstdc++
	rm -f *.o

gdb : mgav.cpp
	gcc ../shared/hypergraph_rw.cpp  -c  -w -O3 -g
	g++ ../shared/metric.cpp -c -w -O3   -g
	g++ ../shared/utils.cpp -c -fopenmp -w -O3 -g
	g++ mgav.cpp -c -w -O3   -g
	g++ partition.cpp -fopenmp -c -w -O3 -g
	g++ -o mgav.exe hypergraph_rw.o metric.o partition.o mgav.o utils.o -fopenmp -O3 -g
	rm -f *.o
clean:
	rm -f *.o
	rm -f *.~
	rm -f *.bin
	rm -f *.txt
