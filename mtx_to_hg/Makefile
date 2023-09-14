main: main.cpp
	gcc graphio.c -c -w -O3
	gcc mmio.c -c  -w -O3
	g++ -std=c++11 main.cpp -c -O3
	g++ -std=c++11 -o mtx_to_vertex_list mmio.o graphio.o main.o -O3
	rm -f *.o
