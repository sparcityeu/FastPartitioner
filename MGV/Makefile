main: mgv.cpp
	g++ ../shared/hypergraph_rw.cpp -c  -w -O3 -std=c++11
	g++ ../shared/metric.cpp -c -O3 -std=c++11
	g++ ../shared/utils.cpp -c -fopenmp -O3 -std=c++11
	g++ mgv.cpp -c -O3 -std=c++11
	g++ partition.cpp -fopenmp -c -O3 -std=c++11
	nvcc -c gpu_partition.cu -Xptxas -v -arch=sm_70
	nvcc -o mgv.exe hypergraph_rw.o mgv.o partition.o metric.o utils.o gpu_partition.o -Xcompiler -fopenmp -O3 -std=c++11 -lstdc++
	rm -f *.o

gdb : mgv.cpp
	g++ ../shared/hypergraph_rw.cpp  -c  -w -O3 -pg -p -g
	g++ ../shared/metric.cpp -c -w -O3 -pg -p -g
	g++ ../shared/utils.cpp -c -fopenmp -w -O3 -pg -p -g
	g++ mgv.cpp -c -w -O3 -pg -p -g
	g++ partition.cpp -fopenmp -c -w -O3 -pg -p -g
	g++ -o mgv.exe hypergraph_rw.o metric.o partition.o mgv.o utils.o -fopenmp -O3 -pg -p -g
	rm -f *.o
clean:
	rm -f *.o
	rm -f *.~
	rm -f *.bin
	rm -f *.txt
