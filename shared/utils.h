#ifndef UTILS_H
#define UTILS_H

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <vector>
#include <chrono>
#include <omp.h>
#include <cstring>
#include <random>

using namespace std;

template<typename T>
void initialize(T * arr, int begin_index, int end_index, T value)
{
    for(int i = begin_index; i < end_index; i++)
    {
        arr[i] = value;
    }
}
void createInverse_patohlike(unsigned int * row_ptr, unsigned int * col_ind, unsigned int * net_appearence, int nov, int noe, int pins, unsigned int ** row_ptr_inv_p, unsigned int ** col_ind_inv_p);
void createInverse(unsigned int * row_ptr, unsigned int * col_ind, unsigned int * net_appearence, int nov, int noe, int pins, unsigned int ** row_ptr_inv_p, unsigned int ** col_ind_inv_p);
void colorCoarse(unsigned int *row_ptr_inv, unsigned int *col_ind_inv, int nov, int noe, int no_color);
int bfsFarthest(int start_node, unsigned int * row_ptr, unsigned int * col_ind, unsigned int * row_ptr_inv, unsigned int * col_ind_inv, int nov, int * distance, int NT);
void bfsOrdering(int start_node, unsigned int * row_ptr, unsigned int * col_ind, unsigned int * row_ptr_inv, unsigned int * col_ind_inv, int nov, int * distance, int * ordering, int NT);
void bfsOrderingSeperated(int start_node, unsigned int * row_ptr, int * col_ind, unsigned int * row_ptr_inv, unsigned int * col_ind_inv, int nov, int noe, int * distance, int * src_nodes, int * ordering_src, int * ordering_targ);
int bfsFarthest_queue(int start_node, unsigned int * row_ptr, unsigned int * col_ind, unsigned int * row_ptr_inv, unsigned int * col_ind_inv, int nov, int * distance);
void bfsOrdering_queue(int start_node, unsigned int * row_ptr, unsigned int * col_ind, unsigned int * row_ptr_inv, unsigned int * col_ind_inv, int nov, int * distance, int * ordering);
int bfsFarthest_queue_par(int start_node, unsigned int * row_ptr, unsigned int * col_ind, unsigned int * row_ptr_inv, unsigned int * col_ind_inv,
                            int nov, int * distance, int * glob_queue, int NT);
void bfsOrdering_queue_par(int start_node, unsigned int * row_ptr, unsigned int * col_ind, unsigned int * row_ptr_inv, unsigned int * col_ind_inv,
                            int nov, int * distance, int * glob_queue, int * ordering, int NT);
bool topDown(unsigned int *row, unsigned int *col, unsigned int * row_ptr_inv, unsigned int * col_ind_inv, int *distance, int &level, int *globalQueue, int &globalLen, int *prefixSum, int **localQueuesList, int * local_lengths);
void writePartVector(string f_name, int nov, int K, int * part_vec, string algo_name);
bool partitionCheck(int K);
bool orderingCheck(int ORDERING_TYPE);
void LSH(unsigned int * row_ptr, unsigned int * col_ind, unsigned int * net_appearence, int nov, int noe, int pins);
#endif
