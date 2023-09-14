main: cn.cpp
	g++ ../shared/hypergraph_rw.cpp -c -w -O3 -std=c++11
	g++ ../shared/metric.cpp -c -O3 -std=c++11
	g++ ../shared/utils.cpp -c -fopenmp -w -O3 -std=c++11
	g++ cn.cpp -c -O3 -std=c++11
	g++ partition.cpp -c -fopenmp  -O3 -std=c++11
	nvcc -c gpu_partition.cu -Xptxas -v -arch=sm_70
	nvcc -o cn.exe hypergraph_rw.o partition.o cn.o metric.o utils.o gpu_partition.o -Xcompiler -fopenmp -O3 -std=c++11 -g
	rm -f *.o

gdb : cn.cpp
	gcc ../shared/hypergraph_rw.cpp  -c  -w -O3 -pg -p -g
	g++ ../shared/metric.cpp -c -w -O3 -pg -p -g
	g++ ../shared/utils.cpp -c -w -O3 -pg -p -g
	g++ cn.cpp -c -w -O3 -pg -p -g
	g++ partition.cpp -c -w -O3 -pg -p -g
	g++ -o cn.exe hypergraph_rw.o metric.o partition.o cn.o utils.o -O3 -pg -p -g
	rm -f *.o

clean :
	rm -f *.o
	rm -f *.~
	rm -f *.bin
	rm -f *.txt
