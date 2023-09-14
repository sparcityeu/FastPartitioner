#ifndef PARTITION_H
#define PARTITION_H
#include <chrono>
#include <iostream>
#include <climits>
#include "../shared/metric.h"
#include "../shared/utils.h"
using namespace std;
extern chrono::duration<double, milli> total;

void partition(unsigned int * row_ptr, unsigned int * col_ind, unsigned int * net_appearence, int nov, int noe,long long pins,int * edge_wgths, int * src_nodes, int K, int * ordering, string & f_name ,int & REFI_OPT, int & REFI_ITER,unsigned int * row_ptr_inv, unsigned int * col_ind_inv);

void thin_partition(unsigned int * row_ptr, unsigned int * col_ind, unsigned int * net_appearence, int nov, int noe,long long pins,int * edge_wgths, int * src_nodes, int K, int * ordering, string & f_name ,int & REFI_OPT, int & REFI_ITER,unsigned int * row_ptr_inv, unsigned int * col_ind_inv, int* part_vec);

void partition_par(unsigned int *row_ptr, unsigned int *col_ind, unsigned int * net_appearence, int nov, int noe, long long pins, int *edge_wgths, int *src_nodes, int K, int *ordering, string f_name, int REFI_OPT, int REFI_ITER, int NT, unsigned int * row_ptr_inv, unsigned int * col_ind_inv);

void thin_partition_par(unsigned int *row_ptr, unsigned int *col_ind, unsigned int * net_appearence, int nov, int noe, long long pins, int *edge_wgths, int *src_nodes, int K, int *ordering, string f_name, int REFI_OPT, int REFI_ITER, int NT, unsigned int * row_ptr_inv, unsigned int * col_ind_inv, int* part_vec);

void partition_gpu(unsigned int *row_ptr, unsigned int *col_ind, unsigned int * net_appearence, int nov, int noe, long long pins, int *edge_wgths, int *src_nodes, int K, int *ordering, string f_name, int REFI_OPT, int REFI_ITER, unsigned int * row_ptr_inv, unsigned int * col_ind_inv);

#endif
