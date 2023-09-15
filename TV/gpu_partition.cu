#include <cuda.h>
#include <stdio.h>

#define WARP_COUNT 8192
#define G_FOCUS_VERTEX 503000


//Error check-----
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
  if (code != cudaSuccess)
    {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
    }
}
//Error check-----

__device__ void GPU_LDG(bool &check_selected_vertex, int &K, int &active_part_count, int *part_sizes,
			int &decision, int *part_vec, int *active_parts, double &C, double *common_scores, int &vertex, int &wid, int &lid){

  //LDG STARTS
  double max_val = 0;
  int max_index = -1;
  
  for(int i = 0; i < active_part_count; i++){
    
    int part_id = active_parts[wid*K+i]; //THIS IS A HUGE IF
    double penalty = 1 - (1.0 * part_sizes[part_id] / C);
    if((wid*K+part_id > WARP_COUNT*K) || wid*K+part_id < 0)
      printf("wtf: wid(%d) # part_id(%d), # K(%d) \n", wid, part_id, K);
    double score = common_scores[wid*K+part_id] * penalty;

    //printf("Penalty(%f) - Score(%f) - C(%f)\n - part_sizes[part_id](%d)", penalty, score, C, part_sizes[part_id], vertex);
    
    if (score > 0)
      {
	if ((max_index == -1)
	    || (score > max_val)
	    || ((score == max_val)
		&& (part_sizes[part_id]<part_sizes[max_index])))
	  {
	    max_val = score;
	    max_index = part_id;
	  }
      }
  }
  
  if (max_index == -1)
    {
      //find min loaded
      int min_load = INT_MAX;
      int min_load_index = -1;
      for (int i = 0; i < K; i++)
	{
	  if (part_sizes[i] < min_load)
	    {
	      min_load = part_sizes[i];
	      min_load_index = i;
	    }
	}
      max_index = min_load_index;
    }
  decision = max_index;
}


//IT IS CLEAR THAT EVERY PART NEED THEIR OWN *ACTIVE_PARTS*
__global__ void gpu_ref_tv_warps(unsigned int* row_ptr, unsigned int* col_ind,
				 unsigned int* net_appearence,
				 int* d_nov, int* d_noe, long long* d_pins, int* d_K,
				 int* edge_wghts, int* src_nodes, int* ordering, int* INmarker,
				 int* INactive_parts, int* part_sizes, int* part_vec,
				 int* ep_starts, int* ep_ends, int* ep_parts, int* ep_cnt,
				 double* scores,
				 int* d_TOTAL_WARPS, double* d_C,
				 int* common_mark, int* common_active,
				 int* d_no_threads, int* best_parts, double* common_scores,
				 int* part_before, int* active_part_counts, int* my_vertex){

  //int tid = blockDim.x * blockIdx.x + threadIdx.x;
  //if(tid == 0){
  //printf("#####BEFORE STARTING#####\n");
  //for(int i = 0; i < *d_nov; i++){
  //printf("Vertex: %d -- Part id: %d \n",i , part_vec[i]);
  //}
  //printf("#####BEFORE STARTING#####\n");
  //}
  //__syncthreads();

  int epoch = 0;
  
  int tid = blockDim.x * blockIdx.x + threadIdx.x;
  int K = *d_K;
  
  int cnt_good_dec = 0;
  bool check_selected_vertex = false; //MADE LOCAL TO THREADS
  int delta_tv_count = 0;
  
  int wid = tid/32;
  int lid = threadIdx.x%32;
  //int active_part_count = 0;

  /*
  if(lid == 0){
    if(wid > 1000){
      printf("I'm warp: %d \n", wid);
    }
  }
  __syncthreads();
  */
    
  int warp_count = *d_TOTAL_WARPS;
  int nov = *d_nov;

  double C = *d_C;
  //if(lid == 0)
  //printf("wid: %d \n", wid);
  
  __syncthreads();
  

      
  for(int cursor = wid; cursor < nov; cursor += warp_count){
    //int vertex = ordering[cursor];
    int vertex = cursor;
    if(lid == 0){
      //printf("Will try wid: %d -- vertex: %d \n", wid, vertex);
      active_part_counts[wid] = 0;
      //active_part_count = __shfl_sync(0xffffffff, active_part_count, 0);
    }
    __syncwarp();

    
    int entry = row_ptr[vertex] + lid;
    int boundary = row_ptr[vertex+1];
    
    int decision = 0;
    int part_id = part_vec[vertex];
    
    for(int i = entry; i < boundary; i+=32){
      
      int edge = col_ind[i];
      int start_index = ep_starts[edge];
      int end_cnt = ep_ends[edge];
      
      for(int j = 0; j < end_cnt; j++){
	int score_part_id = ep_parts[start_index + j];
	//BU SCORE PART ID'YI BIRDEN FAZLA KEZ ALMASI DISINDA CALISIYOR
	//printf("score_part_id: %d -- wid: %d -- MUL: %d \n", score_part_id, wid, score_part_id*wid*K);
	if(wid*K+score_part_id >= K*warp_count || wid*K+score_part_id < 0)
	  continue;
	if(common_mark[wid*K+score_part_id] != vertex){
	  common_scores[wid*K+score_part_id] = 1;
	  common_mark[wid*K+score_part_id] = vertex;
	  //printf("wid(%d), tid(%d), vertex(%d), spi(%d), apc(%d), edge(%d) \n", wid, tid, vertex, score_part_id, active_part_counts[wid], edge);
	  common_active[wid*K+atomicAdd(&active_part_counts[wid], 1)] = score_part_id;
	  //__syncwarp();
	  //active_part_count++;__syncwarp();
	  //active_part_count = __shfl_sync(0xffffffff, active_part_count, lid);
	}
	else{
	  common_scores[wid*K+score_part_id]++;
	}
      }
    }
  
    ////////LDG
    
    if(lid == 0){
      //printf("wid: %d, vertex: %d, apc: %d \n", wid, vertex, active_part_count);
      common_scores[wid*K+part_id]--;
      
      
      double max_val = 0;
      int max_index = -1;
      
      for(int i = 0; i < active_part_counts[wid]; i++){
	int part_id_ldg = common_active[wid*K+i];
	if(part_id_ldg >= 0 && part_id_ldg < K){
	  
	  double penalty = 1 - (1.0 * part_sizes[part_id_ldg] / C);
	  double score = common_scores[wid*K+part_id_ldg] * penalty;
	  
	
	  //if(vertex == G_FOCUS_VERTEX)
	  //printf("Vertex: %d, APC: %d, PartId: %d, PartSize: %d, Penalty: %f, Score: %f\n", vertex, active_part_counts[wid], part_id_ldg, part_sizes[part_id_ldg], penalty, score);
	
	
	  if (score > 0)
	    {
	      if ((max_index == -1)
		  || (score > max_val)
		  || ((score == max_val)
		      && (part_sizes[part_id_ldg]<part_sizes[max_index])))
		{
		  max_val = score;
		  max_index = part_id_ldg;
		}
	    }
	}
	
	if (max_index == -1)
	  {
	    //find min loaded
	    int min_load = INT_MAX;
	    int min_load_index = -1;
	    for (int i = 0; i < K; i++)
	      {
		if (part_sizes[i] < min_load)
		  {
		    min_load = part_sizes[i];
		    min_load_index = i;
		  }
	    }
	    max_index = min_load_index;
	  }
	
	decision = max_index;
	//if(vertex == G_FOCUS_VERTEX)
	//printf("Decision: %d \n", decision);
	//decision = __shfl_sync(0xffffffff, decision, 0);
	//best_parts[my_vertex[wid]] = decision;
	//part_before[vertex] = part_id;
	////////
	////////LDG
      }
    }
  
  
    //decision = __shfl_sync(0xffffffff, decision, 0);
    if(lid == 0){
    if(decision != part_id){
      //printf("dec was: %d part was: %d \n", decision, part_id);
      int leave_gain = 0;
      int arrival_loss = 0;
      bool to_remove = false;
      //CHOOSE IF WE WILL GAIN OR LOSE AFTER PLACEMENT
      for(int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++){
	
	int edge = col_ind[i];
	int start_index = ep_starts[edge];
	int end_cnt = ep_ends[edge]; //COULD REDUCE REGISTER USAGE BY PREDEFINING THEM?
	bool found = false;
	
	for(int j = 0; j < end_cnt; j++){
	  
	  if(ep_parts[start_index + j] == part_id && ep_cnt[start_index + j] == 1){
	    leave_gain++;
	  }
	  else if(ep_parts[start_index + j] == decision){
	    found = true;
	  }
	}
	if(!found){
	  arrival_loss++;
	}
      }
      //IF WE GAIN AFTER PLACEMENT
      if(leave_gain >= arrival_loss){
	//cnt_good_dec++;
	part_vec[vertex] = decision;
	//part_sizes[part_id]--;
	//part_sizes[decision]++;
	atomicAdd(&part_sizes[part_id], -1);
	atomicAdd(&part_sizes[decision], 1);
	//printf("wid: %d removed one from part %d and add to part %d NEWPARTSIZES: -:%d +:%d\n", wid, part_id, decision, part_sizes[part_id], part_sizes[decision]);
	to_remove = true;
	//to_remove = __shfl_sync(0xffffffff, to_remove, 0);
	
	int entry = row_ptr[vertex];// + lid;
	int boundary = row_ptr[vertex+1];
	for(int i = entry; i < boundary; i++){
	  int edge = col_ind[i];
	  int start_ind = ep_starts[edge];
	  int end_cnt = ep_ends[edge];
	  
	  for(int j = 0; j < end_cnt; j++){
	    
	    if(ep_parts[start_ind + j] == part_id){
	    //we find the part that we need to decrement a connection
	      ep_cnt[start_ind + j]--;
	      
	      if(ep_cnt[start_ind + j] == 0){ //all connnections are removed
		
		ep_parts[start_ind + j] = ep_parts[start_ind + end_cnt - 1]; //bring the part in the end to the deleted pos
		ep_cnt[start_ind + j] = ep_cnt[start_ind + end_cnt - 1]; //bring the count in the end to the deleted pos
		ep_parts[start_ind + end_cnt - 1] = -1;
		ep_cnt[start_ind + end_cnt - 1] = 0;
		ep_ends[edge]--;
		end_cnt--;
	      }
	      break;
	    }
	  }
	}
	//ADDING
	for(int i = entry; i < boundary; i++){
	  
	  int edge = col_ind[i];
	  int start_index = ep_starts[edge];
	  int end_cnt = ep_ends[edge];
	  bool found = false;
	  
	  for(int j = 0; j < end_cnt; j++){
	    
	    if(ep_parts[start_index + j] == decision){
	      //We find the part that we need to decrement a connection
	      ep_cnt[start_index + j]++;
	      found = true;
	      //found = __shfl_sync(0xffffffff, found, lid);
	      break;
	    }
	  }
	  if(!found){
	    ep_parts[start_index + end_cnt] = decision;
	    ep_cnt[start_index + end_cnt] = 1;
	    ep_ends[edge]++;
	  }
	}
      }
    }
    //__syncthreads();
    }
    //__syncthreads();
  }
}



extern void partition_gpu_wrapper(unsigned int* row_ptr, unsigned int* col_ind, unsigned int* net_appearence, int nov, int noe, long long pins, int* edge_wgths, int* src_nodes, int K, int* ordering, int* marker, int* active_parts, int* part_sizes, int* ep_starts, int* ep_ends, int* ep_parts, int* part_vec, int* ep_cnt, double* scores, int REFI_ITER, double C){

  //nov = 1000;

  int device = 0;
  cudaSetDevice(device);

  //Device Attributes-----
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, device);

  int no_multiprocessors = prop.multiProcessorCount;
  int max_thread_per_block = prop.maxThreadsPerBlock;

  printf("--Device Properties--\n");
  printf("Device: ");
  printf("%s \n", prop.name);
  printf("Number of multiprocessors: %d \n", no_multiprocessors);
  printf("Max threads per block: %d \n", max_thread_per_block);
  printf("--Device Properties--\n");
  gpuErrchk( cudaDeviceSynchronize() );
  //Device Attributes-----

  int b1 = no_multiprocessors*32;
  int t1 = max_thread_per_block/32;
  b1 = WARP_COUNT;
  t1 = 32;
  //b1 = 1;
  //t1 = 1;
  int no_threads = b1*t1;
  int no_warps = WARP_COUNT;
  printf("b1: %d -- t1: %d\n", b1, t1);
  

  //DEVICE ALLOCATIONS FOR GLOBAL ARRAYS--------------------------------------------------- 
  unsigned int* d_row_ptr;
  unsigned int* d_col_ind;
  unsigned int* d_net_appearence;
  int* d_ep_starts;
  int* d_ep_parts;
  int* d_ep_ends;
  int* d_ep_cnt;
  int* d_part_sizes;
  int* d_part_vec;
  //
  int* d_edge_wgths;
  int* d_src_nodes;
  int* d_ordering;
  int* d_marker;
  int* d_active_parts;
  double* d_scores;
  ////
  int* common_mark = new int[K*WARP_COUNT];
  int* common_active = new int[K*WARP_COUNT];
  double* common_scores = new double[K*WARP_COUNT];
  int* d_common_mark;
  int* d_common_active;
  double* d_common_scores;
  //////
  int* best_parts = new int[nov];
  memset(best_parts, 0, nov*sizeof(int));
  int* d_best_parts;


  ////////
  int* part_before = new int[nov];
  int* active_part_counts = new int[WARP_COUNT];
  int* my_vertex = new int[WARP_COUNT];
  memset(part_before, 0, nov*sizeof(int));
  memset(active_part_counts, 0, WARP_COUNT*sizeof(int));
  memset(my_vertex, 0, WARP_COUNT*sizeof(int));
  int* d_part_before;
  int* d_active_part_counts;
  int* d_my_vertex;
  cudaMalloc(&d_part_before, nov*sizeof(int));
  cudaMemcpy(d_part_before, part_before, nov*sizeof(int), cudaMemcpyHostToDevice);
  cudaMalloc(&d_active_part_counts, WARP_COUNT*sizeof(int));
  cudaMemcpy(d_active_part_counts, active_part_counts, WARP_COUNT*sizeof(int), cudaMemcpyHostToDevice);
  cudaMalloc(&d_my_vertex, WARP_COUNT*sizeof(int));
  cudaMemcpy(d_my_vertex, my_vertex, WARP_COUNT*sizeof(int), cudaMemcpyHostToDevice);
  ////////
  
  for(int w = 0; w < WARP_COUNT; w++){
    for(int i = 0; i < K; i++){
      common_mark[K*w+i] = -1;//marker[i];
      common_active[K*w+i] = -1;//active_parts[i];
      common_scores[K*w+i] = 0;//scores[i];
    }
  }
  //////
    

  cudaMalloc(&d_part_vec, nov*sizeof(int));
  cudaMalloc(&d_row_ptr, (nov+1)*sizeof(unsigned int));
  cudaMalloc(&d_col_ind, pins*sizeof(unsigned int));
  cudaMalloc(&d_net_appearence, noe*sizeof(unsigned int));
  cudaMalloc(&d_ep_starts, (noe+1)*sizeof(int));
  cudaMalloc(&d_ep_parts, pins*sizeof(int));
  cudaMalloc(&d_ep_ends, noe*sizeof(int));
  cudaMalloc(&d_ep_cnt, pins*sizeof(int));
  cudaMalloc(&d_part_sizes, K*sizeof(int));
  //
  cudaMalloc(&d_edge_wgths, noe*sizeof(int));
  cudaMalloc(&d_src_nodes, noe*sizeof(int));
  cudaMalloc(&d_ordering, nov*sizeof(int));
  cudaMalloc(&d_marker, K*sizeof(int));
  cudaMalloc(&d_active_parts, K*sizeof(int));
  cudaMalloc(&d_scores, K*sizeof(double));
  ////
  cudaMalloc(&d_common_mark, K*WARP_COUNT*sizeof(int));
  cudaMalloc(&d_common_active, K*WARP_COUNT*sizeof(int));
  cudaMalloc(&d_common_scores, K*WARP_COUNT*sizeof(double));
  //////
  cudaMalloc(&d_best_parts, nov*sizeof(int));
  
  
  cudaMemcpy(d_part_vec, part_vec, nov*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_row_ptr, row_ptr, (nov+1)*sizeof(unsigned int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_col_ind, col_ind, pins*sizeof(unsigned int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_net_appearence, net_appearence, noe*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_ep_starts, ep_starts, (noe+1)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_ep_parts, ep_parts, pins*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_ep_ends, ep_ends, noe*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_src_nodes, src_nodes, noe*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_ep_cnt, ep_cnt, pins*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_part_sizes, part_sizes, K*sizeof(int), cudaMemcpyHostToDevice);
  //
  cudaMemcpy(d_edge_wgths, edge_wgths, noe*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_src_nodes, src_nodes, noe*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_ordering, ordering, nov*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_marker, marker, K*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_active_parts, active_parts, K*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_scores, scores, K*sizeof(double), cudaMemcpyHostToDevice);
  ////
  cudaMemcpy(d_common_mark, common_mark, K*WARP_COUNT*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_common_active, common_active, K*WARP_COUNT*sizeof(int), cudaMemcpyHostToDevice);
  //////
  cudaMemcpy(d_best_parts, best_parts, nov*sizeof(int), cudaMemcpyHostToDevice);
  
  int* d_nov;
  cudaMalloc((void**)&d_nov, 1*sizeof(int));
  cudaMemcpy(d_nov, &nov, 1*sizeof(int), cudaMemcpyHostToDevice);

  int* d_noe;
  cudaMalloc((void**)&d_noe, 1*sizeof(int));
  cudaMemcpy(d_noe, &noe, 1*sizeof(int), cudaMemcpyHostToDevice);

  int* d_K;
  cudaMalloc((void**)&d_K, 1*sizeof(int));
  cudaMemcpy(d_K, &K, 1*sizeof(int), cudaMemcpyHostToDevice);

  int* d_no_threads;
  cudaMalloc((void**)&d_no_threads, 1*sizeof(int));
  cudaMemcpy(d_no_threads, &no_threads, 1*sizeof(int), cudaMemcpyHostToDevice);

  long long* d_pins;
  cudaMalloc((void**)&d_pins, 1*sizeof(long long));
  cudaMemcpy(d_pins, &pins, 1*sizeof(long long), cudaMemcpyHostToDevice);
  
  
  double* d_C;
  cudaMalloc((void**)&d_C, 1*sizeof(double));
  cudaMemcpy(d_C, &C, 1*sizeof(double), cudaMemcpyHostToDevice);

  int* d_NO_WARPS;
  cudaMalloc((void**)&d_NO_WARPS, 1*sizeof(int));
  cudaMemcpy(d_NO_WARPS, &no_warps, 1*sizeof(int), cudaMemcpyHostToDevice);
  
  //DEVICE ALLOCATIONS FOR GLOBAL ARRAYS---------------------------------------------------

    
  for(int i = 0; i < REFI_ITER; i++){
    gpu_ref_tv_warps<<<b1,t1>>>(d_row_ptr, d_col_ind, d_net_appearence,
				d_nov, d_noe, d_pins, d_K,
				d_edge_wgths, d_src_nodes, d_ordering, d_marker,
				d_active_parts, d_part_sizes, d_part_vec,
				d_ep_starts, d_ep_ends, d_ep_parts, d_ep_cnt, d_scores,
				d_NO_WARPS, d_C, d_common_mark, d_common_active,
				d_no_threads, d_best_parts, d_common_scores,
				d_part_before, d_active_part_counts, d_my_vertex);
    gpuErrchk( cudaDeviceSynchronize() );
  }

  cudaMemcpy(part_vec, d_part_vec, nov*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(ep_starts, d_ep_starts, noe*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(ep_parts, d_ep_parts, pins*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(ep_ends, d_ep_ends, noe*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(ep_cnt, d_ep_cnt, pins*sizeof(int), cudaMemcpyDeviceToHost);

  
  
  double avg_weigth = 0.0;
  for(int i = 0; i < K; i++){
    avg_weigth += part_sizes[i];
  }
  avg_weigth /= K;
  
  double a_imbal = 0;
  for(int i = 0; i < K; i++) {
    a_imbal = max(a_imbal, part_sizes[i]/avg_weigth);
  }
  
  printf("#############################\n");
  printf("BEFORE ACTUAL IMBAL: %f \n", a_imbal);
  printf("#############################\n");
  cudaMemcpy(part_sizes, d_part_sizes, K*sizeof(int), cudaMemcpyDeviceToHost);
  
  avg_weigth = 0.0;
  for(int i = 0; i < K; i++){
    avg_weigth += part_sizes[i];
  }
  avg_weigth /= K;
  
  a_imbal = 0;
  for(int i = 0; i < K; i++) {
    a_imbal = max(a_imbal, part_sizes[i]/avg_weigth);
  }
  printf("#############################\n");
  printf("AFTER ACTUAL IMBAL: %f \n", a_imbal);
  printf("#############################\n");
  
}

  




