#include <cuda.h>
#include <stdio.h>

#define WARP_COUNT 1024

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

__global__ void gpu_ref_mgv_warps(unsigned int *row_ptr,
				  unsigned int *col_ind,
				  unsigned int *net_appearence,
				  int* nov,
				  int* noe,
				  long long* pins,
				  int *edge_wgths,
				  int *src_nodes,
				  int* K,
				  int *ordering,
				  int* send_volumes,
				  double* scores,
				  int* marker,
				  int* active_parts,
				  int* ep_starts,
				  int* ep_parts,
				  int* ep_ends,
				  int* ep_cnt,
				  int* part_sizes,
				  int* part_vec,
				  int* max_part_change_detect,
				  int* gains,
				  int* losses,
				  int* max_parts,
				  int* max_part_marker){
  
  int tid = blockDim.x * blockIdx.x + threadIdx.x;
  int wid = tid/32;
  int lid = threadIdx.x%32;

  __shared__ active_part_count;
  
  __syncthreads();

  for(int cursor = wid; cursor < *nov; cursor += WARP_COUNT){
    
    int decision;
    int vertex = ordering[z];
    if(lid == 0)
      active_part_count = 0;
    int part_id = part_vec[vertex];
    
    
    for(int i = lid; i < K; i += 32){
      gains[i] = 0;
      losses[i] = 0;
    }

    for(int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++){
      int edge = col_ind[i];
      int start_index = ep_starts[edge];
      int end_cnt = ep_ends[edge];

      for(int j = 0; j < end_cnt; j++){
	int score_part_id = ep_starts[start_index + j];
	if(marker[wid*K+score_part_id] != vertex){
	  scores[wid*K+score_part_id] = 1;
	  marker[wid*K+score_part_id] = vertex;
	  active_parts[wid*K+active_part_count] = score_part_id;
	  active_part_count++;
	}
	else{
	  scores[wid*K+score_part_id]++;
	}
	
      }
    }
  
    scores[part_id]--;

    bool check_selected_vertex = true;
    //GET SCORES ENDS
    GPU_LDG(check_selected_vertex, K, active_part_count, part_sizes, decision, part_vec, active_parts,
	    C, scores, vertex);

    if(part_id != decision){
      
      for(int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++){
	int edge = col_ind[i];
	int start_index = ep_starts[edge];
	int end_cnt = ep_ends[edge];
	int net_src_node = src_nodes[edge];
	int net_src_part = part_vec[net_src_node];
	
	if(net_src_node == vertex){
	  gains[part_id] += end_cnt - 1;
	  losses[decision] += end_cnt - 1;
	}
	
	bool found = false;
	
	for(int j = 0; j < end_cnt; j++){
	  if(ep_parts[start_index + j] == decision){
	    //We found that we need to decrement a connection
	    if(net_src_node != vertex){
	      gains[net_src_part]++;
	    }
	    else{
	      losses[decision]--;
	    }
	  }
	  if(ep_parts[start_index + j] == decision){
	    found = true;
	  }
	}

	if(!found){
	  if(net_src_node != vertex){
	    losses[net_src_part]++;
	  }
	  else{
	    losses[decision]++;
	  }
	}
      }

      int best_part = -1;
      bool good_decision = true;

      for(int i = 0; i < max_part_count; i++){
	int max_part_id = max_parts[i];

	if(gains[max_part_id] < losses[max_part_id]){
	  good_decision = false;
	  break;
	}
      }


      for(int i = 0; good_decision && i < K; i++){
	if(max_part_marker[i] == max_changed_count){
	  continue;
	}

	
	if(send_columes[i] + (losses[i] - gains[i]) >= max_send_vol){
	  good_decision = false;
	}
      }

      if(good_decision){
	best_part = decision;
	move_cnt++;
      }

      //IF WE GAIN AFTER PLACEMENT
      if(best_part != -1){
	part_vec[vertex] = decision;
	part_sizes[part_id]--;
	part_sizes[decision]++;
	
	
	//REMOVING
	for(int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++){
	  int edge = col_ind[i];
	  int net = edge;
	  int start_ind = ep_starts[edge];
	  int end_cnt = ep_ends[edge];
	  int net_src_node = src_nodes[edge];
	  int net_src_part = part_vec[net_src_node];

	  if(net_src_node == vertex){
	    send_volumes[part_id] -= end_cnt - 1;
	  }
	  
	  for (int j = 0; j < end_cnt; j++)
	    {
	      
	      if (ep_parts[start_ind + j] == part_id)
		{ //we find the part that we need to decrement a connection
		  ep_cnt[start_ind + j]--;
		  if (ep_cnt[start_ind + j] == 0)
		    {                                                                //all connections are removed
		      ep_parts[start_ind + j] = ep_parts[start_ind + end_cnt - 1]; //bring the part in the end to the deleted pos
		      ep_cnt[start_ind + j] = ep_cnt[start_ind + end_cnt - 1];     //bring the count in the end to the deleted pos
		      ep_parts[start_ind + end_cnt - 1] = -1;
		      ep_cnt[start_ind + end_cnt - 1] = 0;
		      ep_ends[edge]--;
		      end_cnt--;
		      
		      if (net_src_node != vertex)
			{
			  send_volumes[net_src_part]--;
			}
		    }
		  //break; // we don't need to search no longer
		}
	    }
	}
	
	//Update the edge part info
	for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++)
	  {
	    int edge = col_ind[i];
	    int end_cnt = ep_ends[edge];
	    int start_index = ep_starts[edge];
	    bool found = false;
	    int net_src_node = src_nodes[edge];
	    int net_src_part = part_vec[net_src_node];
	    
	    if (net_src_node == vertex)
	      {
		send_volumes[net_src_part] += end_cnt - 1;
	      }
	    
	    for (int j = 0; j < end_cnt; j++)
	      {
		int part_id = ep_parts[start_index + j];
		if (part_id == decision) // no dublicates
		  {
		    found = true;
		    ep_cnt[start_index + j]++;
		    break;
		  }
	      }
	    
	    if (!found)
	      {
		ep_parts[start_index + end_cnt] = decision;
		ep_cnt[start_index + end_cnt] = 1;
		ep_ends[edge]++;
		send_volumes[net_src_part]++;
	      }
	  }
	
	
	max_part_count = 0;
	max_send_vol = -1;
	max_changed_count = 0;
	initialize(max_part_marker, 0, K, -1);
	for (int i = 0; i < K; i++)
	  {
	    if (send_volumes[i] > max_send_vol)
	      {
		max_part_count = 1;
		max_parts[0] = i;
		max_send_vol = send_volumes[i];
		max_changed_count++;
		max_part_marker[i] = max_changed_count;
	      }
	    else if (send_volumes[i] == max_send_vol)
	      {
		max_parts[max_part_count++] = i;
		max_part_marker[i] = max_changed_count;
	      }
	  }
	
      }
    }
  }
}

				  



extern void partition_gpu_wrapper(unsigned int *row_ptr,
				  unsigned int *col_ind,
				  unsigned int *net_appearence,
				  int nov,
				  int noe,
				  long long pins,
				  int *edge_wgths,
				  int *src_nodes,
				  int K,
				  int *ordering,
				  int REFI_OPT,
				  int REFI_ITER,
				  int* send_volumes,
				  int* max_part_change_detect,
				  int* gains,
				  int* losses,
				  int* max_parts,
				  int* max_part_marker,
				  double* scores,
				  int* marker,
				  int* active_parts,
				  int* ep_starts,
				  int* ep_parts,
				  int* ep_ends,
				  int* ep_cnt,
				  int* part_sizes,
				  int* part_vec){

  printf("Sa from the gpu wrapper \n");


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
  t1 = 256;
  //b1 = 1;
  //t1 = 1;
  int no_threads = b1*t1;
  int no_warps = WARP_COUNT;
  printf("b1: %d -- t1: %d\n", b1, t1);


  int* send_volumes = new int[WARP_COUNT*K];
  int* max_part_change_detect = new int[WARP_COUNT*K];
  int* gains = new int[WARP_COUNT*K];
  int* losses = new int[WARP_COUNT*K];
  int* max_parts = new int[WARP_COUNT*K];
  int* max_part_marker = new int[WARP_COUNT*K];

  for(int i = 0; i < WARP_COUNT*K; i++){
    send_volumes[i] = 0;
    max_part_change_detect[i] = 0;
    gains[i] = 0;
    losses[i] = 0;
    max_parts[i] = 0;
    max_part_marker[i] = 0;
  }
  
  int* d_send_volumes;
  int* d_max_part_change_detect;
  int* d_gains;
  int* d_losses;
  int* d_max_parts;
  int* d_max_part_marker;
  
  cudaMalloc(&d_send_volumes, K*sizeof(int));
  cudaMalloc(&d_max_part_change_detect, K*sizeof(int));
  cudaMalloc(&d_gains, K*sizeof(int));
  cudaMalloc(&d_losses, K*sizeof(int));
  cudaMalloc(&d_max_parts, K*sizeof(int));
  cudaMalloc(&d_max_part_marker, K*sizeof(int));
  
  
  //DEVICE ALLOCATIONS FOR GLOBAL ARRAYS--------------------------------------------------- 
  unsigned int* d_row_ptr; //OK
  unsigned int* d_col_ind; //OK
  unsigned int* d_net_appearence; //OK
  int* d_ep_starts; //OK
  int* d_ep_parts; //OK
  int* d_ep_ends; //OK
  int* d_ep_cnt; //OK
  int* d_part_sizes; // OK
  int* d_part_vec; //OK
  //
  int* d_edge_wgths; // OK
  int* d_src_nodes; // OK
  int* d_ordering; // OK
  int* d_marker; // OK
  int* d_active_parts; // OK
  double* d_scores; // OK
  int* d_send_volumes; //OK

  cudaMalloc(&d_part_vec, nov*sizeof(int));
  cudaMalloc(&d_row_ptr, (noe+1)*sizeof(unsigned int));
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
  cudaMalloc(&d_send_volumes, K*sizeof(double));
  ////
  //cudaMalloc(&d_best_parts, nov*sizeof(int));
  
  
  cudaMemcpy(d_part_vec, part_vec, nov*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_row_ptr, row_ptr, (noe+1)*sizeof(unsigned int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_col_ind, col_ind, pins*sizeof(unsigned int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_net_appearence, net_appearence, noe*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_ep_starts, ep_starts, noe*sizeof(int), cudaMemcpyHostToDevice);
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
  cudaMemcpy(d_send_volumes, send_volumes, K*sizeof(int), cudaMemcpyHostToDevice);
  ////
    
  
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
  
  double C = (nov/K) * 1.05;
  
  double* d_C;
  cudaMalloc((void**)&d_C, 1*sizeof(double));
  cudaMemcpy(d_C, &C, 1*sizeof(double), cudaMemcpyHostToDevice);

  //int* d_NO_WARPS;
  //cudaMalloc((void**)&d_NO_WARPS, 1*sizeof(int));
  //cudaMemcpy(d_NO_WARPS, &no_warps, 1*sizeof(int), cudaMemcpyHostToDevice);

  for(int i = 0; i < REFI_ITER; i++){
    gpu_ref_mgv_warps<<<b1,t1>>>(d_row_ptr,
				 d_col_ind,
				 d_net_appearence,
				 d_nov,
				 d_noe,
				 d_pins,
				 d_edge_wgths,
				 d_src_nodes,
				 d_K,
				 d_ordering,
				 d_send_volumes,
				 d_scores,
				 d_marker,
				 d_active_parts,
				 d_ep_starts,
				 d_ep_parts,
				 d_ep_ends,
				 d_ep_cnt,
				 d_part_sizes,
				 d_part_vec,
				 send_volumes,
				 max_part_change_detect,
				 gains,
				 losses,
				 max_parts,
				 max_part_marker
				 );
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
