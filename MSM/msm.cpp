#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include "partition.h"
#include "../shared/CLI.hpp"
#include "../shared/utils.h"
#include "../shared/hypergraph_rw.h"

using namespace std;

chrono::duration<double, milli> total(0);

int main(int argc, char *argv[])
{
    CLI::App app("@@@@@@@@@@  HYPERGRAPH PARTITIONER @@@@@@@@@@");
    string f_name, K_CLI = "2", ordering_CLI = "0", refinement_CLI = "0", iter_CLI = "1", num_threads_CLI = "1", gpu_enabled_CLI = "0";

    CLI::Option *opt = app.add_option("-f", f_name, "path to your hypergraph")->required();
    CLI::Option *opt1 = app.add_option("-k", K_CLI, "total partitions");
    CLI::Option *opt2 = app.add_option("-o", ordering_CLI, "0 for ID ordering, 1 for BFS ordering");
    CLI::Option *opt3 = app.add_option("-r", refinement_CLI, "0 for no refinement, 1 for refinement");
    CLI::Option *opt4 = app.add_option("-i", iter_CLI, "refinement iteration number, default is 1");
    CLI::Option *opt5 = app.add_option("-t", num_threads_CLI, "num threads default is 1");
    CLI::Option *opt6 = app.add_option("-g", num_threads_CLI, "gpu refinement, default is 0");

    CLI11_PARSE(app, argc, argv);

    int K = atoi(K_CLI.c_str());                    //number of parts
    int ORDERING_TYPE = atoi(ordering_CLI.c_str()); //from parameters.h
    int REFI_OPT = atoi(refinement_CLI.c_str());
    int REFI_ITER = atoi(iter_CLI.c_str());
    int NT = atoi(num_threads_CLI.c_str());
    int GPE = atoi(num_threads_CLI.c_str());
    //GRAPH VARIABLE INITIALIZATIONS
    unsigned int *row_ptr;        // cumulative net counts from 0'th vertex to n-1'th vertex
    unsigned int *col_ind;        // nonzero values
    unsigned int *net_appearence; // stores the number of times that a net appear in the graph to find the inverse graph faster
    int *edge_wgths;
    int *src_nodes;
    unsigned int nov; // number of vertices in the hyper graph
    unsigned int noe; // number of nets in the hyper graph
    long long pins;   // number of pins(nonzeros)

    //GRAPH READ
    if (!partitionCheck(K) || !orderingCheck(ORDERING_TYPE) || !read_patohlike(f_name, &row_ptr, &col_ind, &net_appearence, nov, noe, pins))
    {
        exit(1);
    }

    assignWeights(f_name, row_ptr, col_ind, net_appearence, nov, noe, pins, &edge_wgths, &src_nodes);

    //Graph read completed
    int *ordering = new int[nov];
    int *distance = new int[nov];
    if (ORDERING_TYPE == 0)
    { //use id's to order
        for (int i = 0; i < nov; i++)
        {
            ordering[i] = i;
        }

	if(GPE){
	  partition_gpu(row_ptr, col_ind, net_appearence, nov, noe, pins, edge_wgths, src_nodes, K, ordering, f_name, REFI_OPT, REFI_ITER);//single threaded
	}
	else{
	if(NT == 1)
	  partition(row_ptr, col_ind, net_appearence, nov, noe, pins, edge_wgths, src_nodes, K, ordering, f_name, REFI_OPT, REFI_ITER);//single threaded
	else
	  partition_par(row_ptr, col_ind, net_appearence, nov, noe, pins, edge_wgths, src_nodes, K, ordering, f_name, REFI_OPT, REFI_ITER,NT);//multi threaded
	
	}
    }
    else
    {
        unsigned int *row_ptr_inv; // cumulative net counts from 0'th vertex to n-1'th vertex
        unsigned int *col_ind_inv; // nonzero values
        int *glob_queue = new int[nov];

        createInverse(row_ptr, col_ind, net_appearence, nov, noe, pins, &row_ptr_inv, &col_ind_inv);
        auto start = chrono::steady_clock::now();
        initialize(ordering, 0, nov, -1);
        initialize(glob_queue, 0, nov, -1);

        int rand_node = rand() % nov;
        int farthest_first = bfsFarthest_queue_par(rand_node, row_ptr, col_ind, row_ptr_inv, col_ind_inv, nov, distance, glob_queue, NT);

        cout << "Farthest node in the first BFS :" << farthest_first << " with distance :" << distance[farthest_first] << endl;

        int farthest_all = bfsFarthest_queue_par(farthest_first, row_ptr, col_ind, row_ptr_inv, col_ind_inv, nov, distance, glob_queue, NT);
        cout << "Farthest node in the second BFS :" << farthest_all << " with distance :" << distance[farthest_all] << endl;

        bfsOrdering_queue_par(farthest_all, row_ptr, col_ind, row_ptr_inv, col_ind_inv, nov, distance, glob_queue, ordering, NT);

        cout << "Starting ordering algortihm." << endl;
        //bfsOrdering(farthest_all, row_ptr, col_ind, row_ptr_inv, col_ind_inv, nov, distance, ordering, NT);
        auto end = chrono::steady_clock::now();
        cout << "Ordering took: " << chrono::duration<double, milli>(end - start).count() / 1000 << " seconds." << endl;
        total += (chrono::duration<double, milli>(end - start));
		if(NT == 1)
        	partition(row_ptr, col_ind, net_appearence, nov, noe, pins, edge_wgths, src_nodes, K, ordering, f_name, REFI_OPT, REFI_ITER);//single threaded
		else
			partition_par(row_ptr, col_ind, net_appearence, nov, noe, pins, edge_wgths, src_nodes, K, ordering, f_name, REFI_OPT, REFI_ITER,NT);//multi threaded
    }
}
