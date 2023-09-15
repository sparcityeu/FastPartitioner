#ifndef METRIC_H
#define METRIC_H
#include <iostream>
#include <vector>
using namespace std;

int CN(int * ep_ends, int nov);
int CN_Refined(int * ep_parts, int * ep_cnt, int * ep_starts, int noe);
int CN_Scratch(int * part_vec, unsigned int * row_ptr, unsigned int * col_ind, unsigned int * net_appearence, int K, int nov, int noe, int pins);
int TV(int * ep_ends, int nov);
int TV_Refined(int * ep_starts, int * ep_cnts, int noe, int K);
int TV_Scratch(int * part_vec, unsigned int * row_ptr, unsigned int * col_ind, unsigned int * net_appearence, int K, int nov, int noe, int pins);
int TM(int * part_messages, int K);
int TM_scratch(int K, int * labels,  int * ep_parts, int * ep_starts, int * ep_ends, int nov,int noe, int * src_nodes);
int TM_Refined(int ** part_con_cnt,int K);
int MGM(int * part_messages, int K);
int MGM_scratch(int K, int * labels,  int * ep_parts, int * ep_starts, int * ep_ends, int nov,int noe, int * src_nodes);
int MGV(int *, int);
int MGV_Scratch(int * part_vec, unsigned int * row_ptr, unsigned int * col_ind, int * src_nodes, int * ep_ends, unsigned int * net_appearence, int K, int nov, int noe, int pins);
int MGAV(int *, int*, int);
int MGAV_Scratch(int * part_vec, unsigned int * row_ptr, unsigned int * col_ind, int * src_nodes, int * ep_ends,  int * ep_parts, int * ep_starts, unsigned int * net_appearence, int K, int nov, int noe, int pins);

#endif
