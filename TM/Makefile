main: tm.cpp
	g++ ../shared/hypergraph_rw.cpp -c  -w -O3 -std=c++11
	g++ ../shared/metric.cpp -c -O3 -std=c++11
	g++ ../shared/utils.cpp -c -fopenmp -O3 -std=c++11
	g++ tm.cpp -c -O3 -std=c++11
	g++ partition.cpp -c -fopenmp -O3 -std=c++11
	nvcc -c gpu_partition.cu  
	nvcc -c gpu_partition_optimized.cu
	nvcc -c gpu_partition_threads.cu -Xptxas -v -arch=sm_52 
	nvcc -o tm.exe hypergraph_rw.o tm.o partition.o metric.o utils.o gpu_partition.o gpu_partition_optimized.o gpu_partition_threads.o -Xcompiler -fopenmp -O3 -arch=sm_52
	rm -f *.o

gdb : tm.cpp
	gcc ../shared/hypergraph_rw.cpp  -c  -w -O3 -pg -p -g
	g++ ../shared/metric.cpp -c -w -O3 -pg -p -g
	g++ ../shared/utils.cpp -c -w -O3 -pg -p -g
	g++ tm.cpp -c -w -O3 -pg -p -g
	g++ partition.cpp -c -w -O3 -pg -p -g
	g++ -o tm.exe hypergraph_rw.o metric.o partition.o tm.o utils.o -O3 -pg -p -g
	rm -f *.o
clean:
	rm -f *.o
	rm -f *.~
	rm -f *.bin
	rm -f *.txt
