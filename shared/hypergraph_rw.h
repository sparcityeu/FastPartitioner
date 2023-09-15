#ifndef HGRW_H
#define HGRW_H
#include <chrono>
#include <iostream>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include "utils.h"

using namespace std;

void assignWeights(string file_name,unsigned int * row_ptr, unsigned int * col_ind, unsigned int * net_appearence, unsigned int nov, unsigned int noe, long long pins, int ** edge_wgths, int ** src_nodes);
void assignWeights_patohlike(string file_name,unsigned int * row_ptr, unsigned int * col_ind, unsigned int * net_appearence, unsigned int nov, unsigned int noe, long long pins, int ** edge_wgths, int ** src_nodes);
bool read_hg(string f_name, unsigned int ** row_ptr_p, unsigned int ** col_ind_p, unsigned int ** net_appearence, unsigned int &nov, unsigned int &noe, long long &pins);
bool read_patohlike(string f_name, unsigned int ** row_ptr_p, unsigned int ** col_ind_p, unsigned int ** net_appearence, unsigned int &nov, unsigned int &noe, long long &pins);

bool thinning(string f_name, unsigned int ** row_ptr_p, unsigned int ** col_ind_p, unsigned int ** net_appearence, unsigned int &nov, unsigned int &noe, long long &pins, int THINMEAN);
#endif
