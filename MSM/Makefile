main: msm.cpp
	g++ ../shared/hypergraph_rw.cpp -c  -w -O3
	g++ ../shared/metric.cpp -c -O3
	g++ ../shared/utils.cpp -c -fopenmp -O3
	g++ msm.cpp -c -O3
	g++ partition.cpp -c -fopenmp -O3
	g++ -o msm.exe hypergraph_rw.o msm.o partition.o metric.o utils.o -fopenmp -O3 -lstdc++
	rm -f *.o

gdb : msm.cpp
	gcc ../shared/hypergraph_rw.cpp  -c  -w -O3 -pg -p -g
	g++ ../shared/metric.cpp -c -w -O3 -pg -p -g
	g++ ../shared/utils.cpp -c -w -O3 -pg -p -g -fopenmp
	g++ msm.cpp -c -w -O3 -pg -p -g
	g++ partition.cpp -c -w -O3 -pg -p -g -fopenmp
	g++ -o msm.exe hypergraph_rw.o metric.o partition.o msm.o utils.o -O3 -pg -p -g -fopenmp
	rm -f *.o
clean:
	rm -f *.o
	rm -f *.~
	rm -f *.bin
	rm -f *.txt
