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

void metrics(int* part_vec, unsigned int* col_ptr, unsigned int* row_ind, int K, int noe, int nov, int* srcs){

  //for(int i = 0; i < 2000; i++){
  //cout << part_vec[i] << endl;
  //}
  
  int* p2p = new int[K*K];
  int* sm = new int[K];
  int* rm = new int[K];
  int* sv = new int[K];
  int* rv = new int[K];
  int* con = new int[noe];

  int* marker = new int[K];
  for(int i = 0; i < K; i++){
    marker[i] = -1;
    sm[i] = rm[i] = sv[i] = rv[i] = 0;
  }

  memset(p2p, 0, sizeof(int)*K*K);
  memset(con, 0, sizeof(int)*noe);
  
  for(int i = 0; i < noe; i++){
    int src = srcs[i];
    int spart = part_vec[src];
    for(int j = col_ptr[i]; j < col_ptr[i+1]; j++){

      int v = row_ind[j];
      int part = part_vec[v];
      
      if(marker[part] != i){
	marker[part] = i;
	con[i]++;
	//cout << i << " -- " << con[i] << endl;
	if(part != spart){
	  sv[spart]++;
	  rv[part]++;
	  if(p2p[spart*K+part] == 0){
	    p2p[spart*K+part] = 1;
	    sm[spart]++;
	    rm[part]++;
	  }
      	}
      }
    }
  }
  int tv, cn, tm, msv, msrv, msm;
  tv = cn = tm = msv = msrv = msm = 0;

  for(int i = 0; i < noe; i++){
    if(con[i] > 0) tv += con[i]-1;
    if(con[i] > 1) cn++;
  }

  for(int i = 0; i < K; i++){
    tm += sm[i];
    if(msm < sm[i]) msm = sm[i];
    if(msv < sv[i]) msv = sv[i];
    if(msrv < sv[i] + rv[i]) msrv = sv[i] + rv[i];
  }


  cout << "Metrics Computed ------/////// \n";
  cout << "TV: " << tv << endl;
  cout << "CN: " << cn << endl;
  cout << "TM: " << tm << endl;
  cout << "MGM: " << msm << endl;
  cout << "MGV: " << msv << endl;
  cout << "MGAV: " << msrv << endl;
  
  delete[] sm;
  delete[] sv;
  delete[] rv;
  delete[] rm;
  delete[] p2p;
  delete[] con;
  delete[] marker;  
}

int main(int argc, char *argv[])
{
  CLI::App app("@@@@@@@@@@  HYPERGRAPH PARTITIONER @@@@@@@@@@");
  string f_name, K_CLI = "2", ordering_CLI = "0", refinement_CLI = "0", iter_CLI = "1", num_threads_CLI = "1", gpu_enabled_CLI= "0", thinning_CLI = "0", thinning_min_CLI = "0", thinning_mean_CLI = "0";
    ;
  
  CLI::Option *opt = app.add_option("-f", f_name, "path to your hypergraph")->required();
  CLI::Option *opt1 = app.add_option("-k", K_CLI, "total partitions");
  CLI::Option *opt2 = app.add_option("-o", ordering_CLI, "0 for ID ordering, 1 for BFS ordering");
  CLI::Option *opt3 = app.add_option("-r", refinement_CLI, "0 for no refinement, 1 for refinement");
  CLI::Option *opt4 = app.add_option("-i", iter_CLI, "refinement iteration number, default is 1");
  CLI::Option *opt5 = app.add_option("-t", num_threads_CLI, "num threads default is 1");
  CLI::Option *opt6 = app.add_option("-g", gpu_enabled_CLI, "Use GPU refinement algorithm, default is 0");
  CLI::Option *opt7 = app.add_option("-m", thinning_CLI, "Thinning");
  CLI::Option *opt8 = app.add_option("-n", thinning_min_CLI, "Thinning min");
  CLI::Option *opt9 = app.add_option("-z", thinning_mean_CLI, "Thinning mean");
  
  
  CLI11_PARSE(app, argc, argv);
  
  int K = atoi(K_CLI.c_str());                    //number of parts
  int ORDERING_TYPE = atoi(ordering_CLI.c_str()); //from parameters.h
  int REFI_OPT = atoi(refinement_CLI.c_str());
  int REFI_ITER = atoi(iter_CLI.c_str());
  int NT = atoi(num_threads_CLI.c_str());
  int GPE = atoi(gpu_enabled_CLI.c_str());
  int THIN = atoi(thinning_CLI.c_str());
  int THINMIN = atoi(thinning_min_CLI.c_str());
  int THINMEAN = atoi(thinning_mean_CLI.c_str());
  //GRAPH VARIABLE INITIALIZATIONS
  unsigned int *row_ptr, *col_ptr, *row_ind;        // cumulative net counts from 0'th vertex to n-1'th vertex
  unsigned int *col_ind;       // nonzero values
  unsigned int *net_appearence; // stores the number of times that a net appear in the graph to find the inverse graph faster
  unsigned int *vertex_appearence; // stores the number of times that a vertex appears in the graph to find the inverse graph faster
  int *edge_wgths;
  int *src_nodes;
  unsigned int nov; // number of vertices in the hyper graph
  unsigned int noe; // number of nets in the hyper graph
  long long pins;   // number of pins(nonzeros)

  //
  

  //THIN VERSIONS
  unsigned int *thin_row_ptr, * thin_col_ptr;        // cumulative net counts from 0'th vertex to n-1'th vertex
  unsigned int *thin_col_ind, * thin_row_ind;        // nonzero values
  unsigned int *thin_net_appearence; // stores the number of times that a net appear in the graph to find the inverse graph faster
  unsigned int *thin_vertex_appearence; // stores the number of times that a vertex appears in the graph to find the inverse graph faster
  int *thin_edge_wgths;
  int *thin_src_nodes;
  unsigned int thin_nov; // number of vertices in the hyper graph
  unsigned int thin_noe; // number of nets in the hyper graph
  long long thin_pins;   // number of pins(nonzeros)


  
  
  //GRAPH READ PATOH-LIKE
  if (!partitionCheck(K) || !orderingCheck(ORDERING_TYPE) || !read_patohlike(f_name, &col_ptr, &row_ind, &vertex_appearence, nov, noe, pins))
    {
      exit(1);
    }

  
  row_ptr = new unsigned int[nov+1];
  col_ind = new unsigned int[pins];

  memset(row_ptr, 0, sizeof(unsigned int)*(nov+1));

  for(int i = 0; i < noe; i++){
    for(int j = col_ptr[i]; j < col_ptr[i + 1]; j++){
      row_ptr[row_ind[j]+1]++;
    }
  }

  for(int i = 1; i <= nov; i++){
    row_ptr[i] += row_ptr[i - 1];
  }

  for(int i = 0; i < noe; i++){
    for(int j = col_ptr[i]; j < col_ptr[i + 1]; j++){
      col_ind[row_ptr[row_ind[j]]++] = i;
    }
  }

  for(int i = nov; i > 0; i--){
    row_ptr[i] = row_ptr[i-1];
  }
  row_ptr[0] = 0;
  
  src_nodes = new int[noe];

  
  for(int i = 0; i < noe; i++){
    int degree = col_ptr[i+1] - col_ptr[i];
    if(degree > 0){
      int src_ind = rand() % (col_ptr[i+1] - col_ptr[i]);
      int src = row_ind[col_ptr[i] + src_ind];
      src_nodes[i] = src;
    }else{
      cout << "Empty net is detected! -> " << i << endl;
      src_nodes[i] = -1;
    }
  }
  
  if(THIN){
    
    
    thin_net_appearence = new unsigned int[noe];
        
    thin_nov = nov;
    double mean = pins/(noe + .0f);
    //int* net_perm = new int[noe];
    int* i_net_perm = new int[noe];

    for(int i = 0; i < noe; i++){
      i_net_perm[i] = -1;
      //net_perm[i] = -1;
    }
    thin_noe = 0;

    thin_col_ptr = new unsigned int[noe+1];
    thin_row_ind = new unsigned int[pins];

    thin_col_ptr[0] = 0;
    
    thin_pins = 0;
    for(int i = 0; i < noe; i++){
      int degree = col_ptr[i + 1] - col_ptr[i];
      if(degree < THINMEAN*mean && degree < THINMIN){
	//net_perm[n_count] = i;
	thin_net_appearence[thin_noe] = degree;
	i_net_perm[i] = thin_noe;
	
	for(int j = col_ptr[i]; j < col_ptr[i+1]; j++){
	  thin_row_ind[thin_pins++] = row_ind[j];
	}
	
	thin_noe++;
	thin_col_ptr[thin_noe] = thin_pins;
      }
    }
    
    thin_row_ptr = new unsigned int[nov + 1];
    thin_col_ind = new unsigned int[pins];
    
    thin_pins = 0;
    for(int i = 0; i < nov; i++){
      for(int j = row_ptr[i]; j < row_ptr[i + 1]; j++){
	int net = col_ind[j];

	if(i_net_perm[net] != -1){
	  thin_col_ind[thin_pins++] = i_net_perm[net];
	}
      }
      thin_row_ptr[i+1] = thin_pins;
    }

    thin_row_ptr[0] = 0;

    thin_src_nodes = new int[thin_noe];
    for(int i = 0; i < noe; i++){
      if(i_net_perm[i] != -1){
	thin_src_nodes[i_net_perm[i]] = src_nodes[i];
      }
    }
    
    //thinning(f_name, &thin_row_ptr, &thin_col_ind,&thin_vertex_appearence, thin_nov, thin_noe, thin_pins, THINMEAN);

    cout << "Graph is thinned \n \t" << "Nov = " << thin_nov << "\n\t" << " noe = " << thin_noe <<
      "\n\t pins = " << thin_pins << endl;
      
      }
  
  //assignWeights_patohlike(f_name, row_ptr, col_ind, vertex_appearence, nov, noe, pins, &edge_wgths, &src_nodes);
  //createInverse_patohlike(row_ptr, col_ind, vertex_appearence, nov, noe, pins, &row_ptr_inv, &col_ind_inv);
  //GRAPH READ
  /*if (!partitionCheck(K) || !orderingCheck(ORDERING_TYPE) || !read_hg(f_name, &row_ptr, &col_ind, &net_appearence, nov, noe, pins))
    {
    exit(1);
    }
    
    assignWeights(f_name, row_ptr, col_ind, net_appearence, nov, noe, pins, &edge_wgths, &src_nodes);
    createInverse(row_ptr, col_ind, net_appearence, nov, noe, pins, &row_ptr_inv, &col_ind_inv);*/
  net_appearence = new unsigned int[noe];
  for(int i = 0; i < noe - 1; i++)
    {
      net_appearence[i] = col_ptr[i + 1] - col_ptr[i];
    }
  
  
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
	partition_gpu(row_ptr, col_ind, net_appearence, nov, noe, pins, edge_wgths, src_nodes, K, ordering, f_name, REFI_OPT, REFI_ITER, col_ptr, row_ind);//patohlike and GPU
      }
      else{
      if(NT == 1)
	{	  
	  if(THIN){
	    int* part_vec = new int[nov];
	    
	    thin_partition(thin_row_ptr, thin_col_ind, thin_net_appearence, thin_nov, thin_noe, thin_pins, edge_wgths, thin_src_nodes, K, ordering, f_name, REFI_OPT, REFI_ITER, thin_col_ptr, thin_row_ind, part_vec);//patohlike

	    metrics(part_vec,col_ptr,row_ind, K,noe, nov,src_nodes);
	  }
	  
	  else{
	    partition(row_ptr, col_ind, net_appearence, nov, noe, pins, edge_wgths, src_nodes, K, ordering, f_name, REFI_OPT, REFI_ITER, col_ptr, row_ind);//patohlike
	  }
	  //partition(row_ptr, col_ind, net_appearence, nov, noe, pins, edge_wgths, src_nodes, K, ordering, f_name, REFI_OPT, REFI_ITER, row_ptr_inv, col_ind_inv);//single threaded
	}
      else
	{
	  if(THIN){
	    int* part_vec = new int[nov];
	    
	    thin_partition_par(thin_row_ptr, thin_col_ind, thin_net_appearence, thin_nov, thin_noe, thin_pins, edge_wgths, thin_src_nodes, K, ordering, f_name, REFI_OPT, REFI_ITER, NT, thin_col_ptr, thin_row_ind,part_vec);//patohlike
	    metrics(part_vec,col_ptr,row_ind, K,noe, nov,src_nodes);
	  }
	  else{
	    partition_par(row_ptr, col_ind, net_appearence, nov, noe, pins, edge_wgths, src_nodes, K, ordering, f_name, REFI_OPT, REFI_ITER, NT, col_ptr, row_ind);//patohlike
	  }
	  //partition_par(row_ptr, col_ind, net_appearence, nov, noe, pins, edge_wgths, src_nodes, K, ordering, f_name, REFI_OPT, REFI_ITER, NT, row_ptr_inv, col_ind_inv);//single threaded
	}
      
      }
    }
  else
    {
      //unsigned int *row_ptr_inv; // cumulative net counts from 0'th vertex to n-1'th vertex
      //unsigned int *col_ind_inv; // nonzero values
      int *glob_queue = new int[nov];
      
      //createInverse(row_ptr, col_ind, net_appearence, nov, noe, pins, &row_ptr_inv, &col_ind_inv);
      cout << "Starting ordering algortihm." << endl;	
      auto start = chrono::steady_clock::now();
      initialize(ordering, 0, nov, -1);
      initialize(glob_queue, 0, nov, -1);
      
      int rand_node = rand() % nov;
      //int farthest_first = bfsFarthest_queue_par(rand_node, row_ptr, col_ind, row_ptr_inv, col_ind_inv, nov, distance, glob_queue, NT);
      int farthest_first = bfsFarthest_queue_par(rand_node, row_ptr, col_ind, col_ptr, row_ind, nov, distance, glob_queue, NT);//patohlike
      
      cout << "Farthest node in the first BFS :" << farthest_first << " with distance :" << distance[farthest_first] << endl;
      
      //int farthest_all = bfsFarthest_queue_par(farthest_first, row_ptr, col_ind, row_ptr_inv, col_ind_inv, nov, distance, glob_queue, NT);
      int farthest_all = bfsFarthest_queue_par(farthest_first, row_ptr, col_ind, col_ptr, row_ind, nov, distance, glob_queue, NT);//patohlike
      cout << "Farthest node in the second BFS :" << farthest_all << " with distance :" << distance[farthest_all] << endl;
      
      //bfsOrdering_queue_par(farthest_all, row_ptr, col_ind, row_ptr_inv, col_ind_inv, nov, distance, glob_queue, ordering, NT);
      bfsOrdering_queue_par(farthest_all, row_ptr, col_ind, col_ptr, row_ind, nov, distance, glob_queue, ordering, NT);//patohlike
      
      
      //bfsOrdering(farthest_all, row_ptr, col_ind, row_ptr_inv, col_ind_inv, nov, distance, ordering, NT);
      auto end = chrono::steady_clock::now();
      cout << "Ordering took: " << chrono::duration<double, milli>(end - start).count() / 1000 << " seconds." << endl;
      total += (chrono::duration<double, milli>(end - start));
      if(NT == 1)
	{
	  partition(row_ptr, col_ind, net_appearence, nov, noe, pins, edge_wgths, src_nodes, K, ordering, f_name, REFI_OPT, REFI_ITER, col_ptr, row_ind);//patohlike
	  
	  //partition(row_ptr, col_ind, net_appearence, nov, noe, pins, edge_wgths, src_nodes, K, ordering, f_name, REFI_OPT, REFI_ITER, row_ptr_inv, col_ind_inv);//single threaded
	}
      else
	{
	  partition_par(row_ptr, col_ind, net_appearence, nov, noe, pins, edge_wgths, src_nodes, K, ordering, f_name, REFI_OPT, REFI_ITER, NT, col_ptr, row_ind);//patohlike
	  
	  //partition_par(row_ptr, col_ind, net_appearence, nov, noe, pins, edge_wgths, src_nodes, K, ordering, f_name, REFI_OPT, REFI_ITER, NT, row_ptr_inv, col_ind_inv);//multi threaded
	}
      
    }

  
  
}
