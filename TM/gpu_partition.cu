#include <cuda.h>
#include <stdio.h>


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

//CREATE CRITICAL REGION HELPERS-------------------------
__device__ volatile int sem = 0;

__device__ void acquire_semaphore(volatile int *lock){
  while (atomicCAS((int *)lock, 0, 1) != 0);
  }

__device__ void release_semaphore(volatile int *lock){
  *lock = 0;
  __threadfence();
  }
//CREATE CRITICAL REGION HELPERS-------------------------


__global__ void gpu_ref(int* d_receivers, int* d_senders,
			int* d_ref_marker_receiver, int* d_ref_marker_sender,
			int* d_part_vec, int* d_nov,
			unsigned int* d_row_ptr, unsigned int* d_col_ind,
			unsigned int* d_net_appearence,int* d_ep_starts, int* d_ep_parts,
			int* d_ep_ends, int* d_src_nodes, int* d_ep_cnt, int* d_gain_loss_1d,
			int* d_part_connections, int* d_K, int* d_best_parts,
			int* d_part_sizes, int* d_part_messages){
  //printf("Nov: %d \n", *d_nov); Referring to int
  int bd = blockDim.x;
  int bi = blockIdx.x;
  int tx = threadIdx.x;
  int tid = blockDim.x * blockIdx.x + threadIdx.x;
  
  int K = *d_K;

  //
  extern int __shared__ arrays[];
  //
  
  __shared__ int v;
  __shared__ int best_part;
  __shared__ int best_improvement;
  __shared__ int part_id;
  //These are normally partition-global, but 1M thread atomic update is impossible(?), so they will be vertex-wise(C0)
  volatile __shared__ int receiver_cnt;
  volatile __shared__ int sender_cnt;
  __shared__ int gain;
  __shared__ int loss;
  __shared__ int my_cnt;
  __shared__ int tot_improvement;
  //
  //__shared__ int improvements[1024];
  //
  //__shared__ int s_ref_marker_receiver[1024];
  //__shared__ int s_ref_marker_sender[1024];
  //__shared__ int s_receivers[1024];
  //__shared__ int s_senders[1024];
  //

  int* improvements = &arrays[0];
  int* s_ref_marker_receiver = &arrays[K];
  int* s_ref_marker_sender = &arrays[2*K];
  int* s_receivers = &arrays[3*K];
  int* s_senders = &arrays[4*K];
    
  if(tx == 0){
    v = bi;
    best_part = -1;
    best_improvement = 0;
    part_id = d_part_vec[v];
    //(C0)
    receiver_cnt = 0;
    sender_cnt = 0;
    gain = 0;
    loss = 0;
    my_cnt = 0;
    tot_improvement = 0;
    
    //printf("Initializing\n");
    for(int i = 0; i < K; i++){
      improvements[i] = -1;
      s_ref_marker_receiver[i] = -1;
      s_ref_marker_sender[i] = -1;
      s_receivers[i] = -1;
      s_senders[i] = -1;
    }
    //printf("Initialized\n");
    
    //(C0)
  }
  //----__syncthreads();
  
  unsigned int entry = d_row_ptr[v] + tx;
  unsigned int boundary = d_row_ptr[v+1];
  
  
  for(int n_ind = entry; n_ind < boundary; n_ind += 32){
    
    //if(n_ind+tx >= boundary)
    //break;
    
    //int net = d_col_ind[n_ind+tx];
    int net = d_col_ind[n_ind];
    int start_index = d_ep_starts[net];
    int end_cnt = d_ep_ends[net];
    
    if(d_src_nodes[net] == v){
      for(int j = 0; j < end_cnt; j++){
	int target_part = d_ep_parts[start_index + j];

	if(target_part < 0)
	  continue;
	
	if(target_part == part_id && d_ep_cnt[start_index + j] > 1){
	  if(s_ref_marker_receiver[part_id] != v){
	    s_receivers[receiver_cnt++] = part_id; //Mind the race condition (C1)
	    //atomicAdd(&receiver_cnt, 1);
	    s_ref_marker_receiver[part_id] = v;
	  }
	}
	
	if(target_part != part_id){
	  if(s_ref_marker_receiver[target_part] != v){
	    s_receivers[receiver_cnt++] = target_part;//(C1) Bu cok feci waw, raw, war bütün hazardlar
	    //atomicAdd(&receiver_cnt, 1);
	    s_ref_marker_receiver[target_part] = v;
	  }
	  d_gain_loss_1d[part_id*K+target_part]++;
	  
	  if(d_part_connections[part_id*K+target_part] == 1 || d_part_connections[part_id*K+target_part] == d_gain_loss_1d[part_id*K+target_part])
	    //gain++; //(C1)
	    atomicAdd(&gain,1);
	}
      } //end_cnt ends
      
      for(int j = 0; j < end_cnt; j++){
	int target_part = d_ep_parts[start_index + j];

	if(target_part < 0)
	  continue;
	
	if(target_part != part_id){
	  d_gain_loss_1d[part_id*K+target_part]--;
	}
      }
    }
    
    else{
      int p_sender = d_part_vec[d_src_nodes[net]];
      if(s_ref_marker_sender[p_sender] != v){
	s_senders[sender_cnt++] = p_sender;
	//atomicAdd(&sender_cnt, 1);
	s_ref_marker_sender[p_sender] = v;
	int end_cnt = d_ep_ends[net];
	for(int j = 0; j < end_cnt; j++){
	  if(d_ep_parts[start_index + j] == part_id && d_ep_cnt[start_index + j] == 1){
	    if(d_part_connections[p_sender*K+part_id] == 1){
	      atomicAdd(&gain,1);
	    }
	  }
	}
      }
    }
    
  }
  //__syncthreads();
  
  if(gain > 0){
    int fbp_start = tx;
    //int fbp_end = *d_K;
    //int part_loss = 0;
    //printf("fbp_start: %d, fbp_end: %d\n", fbp_start, fbp_end);
    int improvement = 0;
    for(int k = fbp_start; k < *d_K; k += 32){
      improvement = gain; //Unlike vertex-global gain, this should be private for each thread
      
      if(k >= *d_K)
	break;
      
      if(k == part_id)
	continue;
      
      
      for(int i = 0; i < receiver_cnt; i++){
	int rec_part = s_receivers[i];

	if(rec_part < 0)
	  continue;
	
	if(rec_part == k)
	  continue;
	if(d_part_connections[(k*K)+rec_part] == 0){
	  improvement--;
	}
      }
      
      for(int i = 0; i < sender_cnt; i++){
	int sender_part = s_senders[i];
	if(sender_part == k)
	  continue;
	if(d_part_connections[(sender_part*K)+k] == 0){
	  improvement--;
	}
      }

      improvements[k] = improvement;
    }
    // __syncthreads();
  }
    if(tx == 0){
      int my_max = 0;
      for(int t = 0; t < K; t++){
	//printf("Will access improvements[%d], d_best_parts[%d]\n", t, bi);
	if(improvements[t] > my_max){
	  my_max = improvements[t];
	  d_best_parts[v] = t;
	}
      }
    }

    //REMOVING
    if(d_best_parts[v] != -1){
      int vertex = v;
      int decision = d_best_parts[v];
      
      if(tx == 0){
	d_part_vec[v] = decision;
	d_part_sizes[part_id]--;
	d_part_sizes[decision]++;
      }
      
      int begin = d_row_ptr[vertex]+tx;
      int end = d_row_ptr[vertex+1];

      for(int i = begin; i < end; i += 32){
	int edge = d_col_ind[i];
	int net = edge;
	int start_ind = d_ep_starts[edge];
	int end_cnt = d_ep_ends[edge];

	for(int j = 0; j < end_cnt; j++){

	  if(d_ep_parts[start_ind + j] == part_id){
	    d_ep_cnt[start_ind +j]--;

	    if(d_ep_cnt[start_ind + j] == 0){
	      d_ep_parts[start_ind + j] = d_ep_parts[start_ind + end_cnt -1];
	      d_ep_cnt[start_ind + j] = d_ep_cnt[start_ind + end_cnt - 1];
	      d_ep_parts[start_ind + end_cnt -1] = -1;
	      d_ep_cnt[start_ind + end_cnt - 1] = 0;
	      d_ep_ends[edge]--;
	    }
	    break;
	  }
	}
      }

      //Update the edge part info
      
      for(int i = begin; i < end; i += 32){
	int edge = d_col_ind[i];
	int net = edge;
	int start_ind = d_ep_starts[edge];
	int end_cnt = d_ep_ends[edge];

	for(int j = 0; j < end_cnt; j++){
	  
	  if(d_ep_parts[start_ind + j] == part_id){
	    d_ep_cnt[start_ind + j]--;

	    if(d_ep_cnt[start_ind + j] == 0){
	      d_ep_parts[start_ind + j] = d_ep_parts[start_ind + end_cnt -1];
	      d_ep_cnt[start_ind + j] = d_ep_cnt[start_ind + end_cnt -1];
	      d_ep_parts[start_ind + end_cnt - 1] = -1;
	      d_ep_cnt[start_ind + end_cnt -1] = 0;
	      d_ep_ends[edge]--;
	    }
	    break;
	  }
	}
      }

      //Update the edge part info

      for(int i = begin; i < end; i+=32){
	int edge = d_col_ind[i];
	int end_cnt = d_ep_ends[edge];
	int start_index = d_ep_starts[edge];
	bool found = false;

	for(int j = 0; j < end_cnt; j++){
	  int part_idLOCAL = d_ep_parts[start_index + j];

	  if(part_idLOCAL == decision){
	    found = true;
	    d_ep_cnt[start_index + j]++;
	    break;
	  }
	}

	if(!found){
	  d_ep_parts[start_index + end_cnt] = decision;
	  d_ep_cnt[start_index + end_cnt] = 1;
	  d_ep_ends[edge]++;
	}
      }

      
      for(int n_ind = begin; n_ind < end; n_ind += 32){
	int net = d_col_ind[n_ind];
	int start_index = d_ep_starts[net];
	int end_cnt = d_ep_ends[net];

	if(d_src_nodes[net] == v){
	  for(int j = 0; j < end_cnt; j++){
	    int target_part = d_ep_parts[start_index + j];

	    if(target_part < 0)
	      continue;

	    if(target_part == decision){

	      d_part_connections[(decision*K)+target_part]++;

	      if(d_part_connections[(decision*K)+target_part] == 1){
		d_part_messages[decision]++;
	      }
	    }

	    if(target_part != part_id && target_part != decision){
	      d_part_connections[(part_id*K)+target_part]--;

	      if(d_part_connections[(part_id*K)+target_part] == 0) {
		d_part_messages[part_id]--;
	      }
	    }

	    if(target_part == decision && d_ep_cnt[start_index + j] > 1){
	      d_part_connections[(part_id*K)+target_part]--;

	      if(d_part_connections[(part_id*K)+target_part] == 0){
		d_part_messages[part_id]--;
	      }
	    }
	  }
	}
	else{
	  int src_net = d_part_vec[d_src_nodes[net]];
	  int end_cnt = d_ep_ends[net];
	  bool found = false;

	  for(int j = 0; j < end_cnt; j++){

	    if(d_ep_parts[start_index + j] == decision && d_ep_cnt[start_index + j] == 1){
	      d_part_connections[(src_net*K)+decision]++;

	      if(d_part_connections[(src_net*K)+decision] == 1){
		d_part_messages[src_net]++;
	      }
	    }

	    if(d_ep_parts[start_index + j] == part_id)
	      found = true;

	    if(!found){
	      d_part_connections[(src_net*K)+part_id]--;

	      if(d_part_connections[(src_net*K)+part_id] == 0){
		d_part_messages[src_net]--;
	      }
	    }
	  }
	}
      }
    }
    
}
							   
      
      
	      
      

    



///////////////
//__syncthreads();
//if(gain > 0)
//printf("Gain: %d, vertex: %d, receiver_cnt: %d, my_net: %d, thread: %d \n", gain, v, receiver_cnt, net, tid);
//__syncthreads();
///////////////
///////////////
//printf("end_cnt: %d \n", end_cnt);
///////////////

//TO DO: best_part and (consequently) add-remove to new parts should be in a buffered manner like in streaming
extern void partition_gpu_wrapper(unsigned int* row_ptr, unsigned int* col_ind, 
				  unsigned int* net_appearence, int nov, int noe, 
				  long long pins, int* edge_wgths, int* src_nodes, 
				  int K, int* ordering, int REFI_OPT, int REFI_ITER,
				  int* part_vec,
				  int* ep_starts, int* ep_parts, int* ep_ends, int* ep_cnt,
				  int** part_connections, int* part_sizes, int* part_messages){

  
  //Device Attributes-----
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);

  int warp_size = prop.warpSize;
  int shared_per_block = (int)prop.sharedMemPerBlock;
  int shared_mem_per_multiprocessor = (int)prop.sharedMemPerMultiprocessor;
  int no_multiprocessors = prop.multiProcessorCount;
  int max_thread_per_block = prop.maxThreadsPerBlock;
  //int max_blocks_per_multiprocessor = prop.maxBlocksPerMultiProcessor;
  int max_threads_per_multiprocessor = prop.maxThreadsPerMultiProcessor;
  int reg_per_multiprocessor = prop.regsPerMultiprocessor;
  
  printf("--Device Properties--\n");
  printf("Warp size: %d \n", warp_size);
  printf("Shared Memory Per block: %d bytes, %d integers\n", shared_per_block, shared_per_block/4);
  printf("Shared Memory Per multiprocessor: %d bytes, %d integers \n", shared_mem_per_multiprocessor, shared_mem_per_multiprocessor/4);
  printf("Number of multiprocessors: %d \n", no_multiprocessors);
  printf("Max threads per block: %d \n", max_thread_per_block);
  //printf("Max blocks per multiprocessor: %d, \n", max_blocks_per_multiprocessor);
  printf("Max threads per multiprocessor: %d \n", max_threads_per_multiprocessor);
  printf("Number of registers per multiprocessor: %d \n", reg_per_multiprocessor);
  printf("--Device Properties--\n");
  gpuErrchk( cudaDeviceSynchronize() );
  //Device Attributes-----

  //Device insights for occupancy-----
  int block_per_sm = shared_mem_per_multiprocessor / shared_per_block;
  
  
  
  //for(int i = 0; i < 100; i++){
  //printf("in wrapper ep_ends[%d]: %d \n", i, ep_ends[i]);
  //}

  //Flatten Part Connections-----
  int* f_part_connections = new int[K*K];
  for(int row = 0; row < K; row++){
    for(int col = 0; col < K; col++ ){
      f_part_connections[(row*K)+col] = part_connections[row][col];
    }
  }

  int* d_part_connections;
  cudaMalloc(&d_part_connections, K*K*sizeof(int));
  cudaMemcpy(d_part_connections, f_part_connections, K*K*sizeof(int), cudaMemcpyHostToDevice);
  //Flatten Part Connections-----
  
  //DEVICE ALLOCATIONS FOR GLOBAL ARRAYS---------------------------------------------------
  int* d_part_vec;
  unsigned int* d_row_ptr;
  unsigned int* d_col_ind;
  unsigned int* d_net_appearence;
  int* d_ep_starts;
  int* d_ep_parts;
  int* d_ep_ends;
  int* d_src_nodes;
  int* d_ep_cnt;
  int* d_part_sizes;
  int* d_part_messages;
  
  cudaMalloc(&d_part_vec, nov*sizeof(int));
  cudaMalloc(&d_row_ptr, (noe+1)*sizeof(unsigned int));
  cudaMalloc(&d_col_ind, pins*sizeof(unsigned int));
  cudaMalloc(&d_net_appearence, noe*sizeof(unsigned int));
  cudaMalloc(&d_ep_starts, noe*sizeof(int));
  cudaMalloc(&d_ep_parts, pins*sizeof(int));
  cudaMalloc(&d_ep_ends, noe*sizeof(int));
  cudaMalloc(&d_src_nodes, noe*sizeof(int));
  cudaMalloc(&d_ep_cnt, pins*sizeof(int));
  cudaMalloc(&d_part_sizes, K*sizeof(int));
  cudaMalloc(&d_part_messages, K* sizeof(int));
  
  
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
  cudaMemcpy(d_part_messages, part_messages, K*sizeof(int), cudaMemcpyHostToDevice);
  
  
  int* d_nov;
  cudaMalloc((void**)&d_nov, 1*sizeof(int));
  cudaMemcpy(d_nov, &nov, 1*sizeof(int), cudaMemcpyHostToDevice);

  int* d_K;
  cudaMalloc((void **)&d_K, 1*sizeof(int));
  cudaMemcpy(d_K, &K, 1*sizeof(int), cudaMemcpyHostToDevice);
  //DEVICE ALLOCATIONS FOR GLOBAL ARRAYS---------------------------------------------------
  
  //LOCAL ALLOCATIONS AND COPIES-----------------------------------------------------------
  int* gain_loss_1d;
  int* receivers;
  int* senders;
  int* ref_marker_receiver;
  int* ref_marker_sender;
  int* best_parts;
  int* d_gain_loss_1d;
  int* d_receivers;
  int* d_senders;
  int* d_ref_marker_receiver;
  int* d_ref_marker_sender;
  int* d_best_parts;

  gain_loss_1d = (int*)malloc(K*K*sizeof(int));
  receivers = (int*)malloc(K*sizeof(int));
  senders = (int*)malloc(K*sizeof(int));
  ref_marker_receiver = (int*)malloc(K*sizeof(int));
  ref_marker_sender = (int*)malloc(K*sizeof(int));
  best_parts = (int*)malloc(nov*sizeof(int));
  
  //Mimicking initialize()
  for(int i = 0; i < K; i++){
    //Not canon
    receivers[i] = -1;
    senders[i] = -1;
    //Not canon
    ref_marker_receiver[i] = -1;
    ref_marker_sender[i] = -1;
  }
  
  for(int i = 0; i < K*K; i++){
    gain_loss_1d[i] = 0;
  }

  for(int i = 0; i < nov; i++){
    best_parts[i] = -1;
  }

  
  cudaMalloc(&d_gain_loss_1d, K*K*sizeof(int));
  cudaMalloc(&d_receivers, K*sizeof(int));
  cudaMalloc(&d_senders, K*sizeof(int));
  cudaMalloc(&d_ref_marker_receiver, K*sizeof(int));
  cudaMalloc(&d_ref_marker_sender, K*sizeof(int));
  cudaMalloc(&d_best_parts, nov*sizeof(int));

  cudaMemcpy(d_gain_loss_1d, gain_loss_1d, K*K*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_receivers, receivers, K*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_senders, senders, K*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_ref_marker_receiver, ref_marker_receiver, K*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_ref_marker_sender, ref_marker_sender, K*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_best_parts, best_parts, nov*sizeof(int), cudaMemcpyHostToDevice);
  //LOCAL ALLOCATIONS AND COPIES------------------------------------------------------------
  

  for(int i = 0; i < REFI_ITER; i++){
    gpu_ref<<<nov,32, 5*K*sizeof(int)>>>(d_receivers,d_senders, d_ref_marker_receiver, d_ref_marker_sender, d_part_vec, d_nov, d_row_ptr, d_col_ind, d_net_appearence, d_ep_starts, d_ep_parts, d_ep_ends, d_src_nodes, d_ep_cnt, d_gain_loss_1d, d_part_connections, d_K, d_best_parts, d_part_sizes, d_part_messages);
  gpuErrchk( cudaDeviceSynchronize() );
  }
  
  //Debug purposes
  //printf("Just before calling the kernel\n");
  ///gpu_ref<<<1000,32>>>(d_receivers,d_senders, d_ref_marker_receiver, d_ref_marker_sender, d_part_vec, d_nov, d_row_ptr, d_col_ind, d_net_appearence, d_ep_starts, d_ep_parts, d_ep_ends, d_src_nodes, d_ep_cnt, d_gain_loss_1d, d_part_connections, d_K);
  
  //gpuErrchk( cudaPeekAtLastError() );
  


  //
  cudaMemcpy(part_vec, d_part_vec, nov*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(ep_starts, d_ep_starts, noe*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(ep_parts, d_ep_parts, pins*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(ep_ends, d_ep_ends, noe*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(src_nodes, d_src_nodes, noe*sizeof(int), cudaMemcpyDeviceToHost);
  //
  

  
  //for(int i = 0; i < nov; i++){
  //printf("Vertex: (%d), Best_part(%d)\n", i, best_parts[i]);
  //if(best_parts[i] != -1){
  //part_vec[i] = best_parts[i];
  //}
  //}
  

  free(f_part_connections);
  free(gain_loss_1d);
  free(receivers);
  free(senders);
  free(ref_marker_receiver);
  free(ref_marker_sender);
  free(best_parts);

  cudaFree(d_part_vec);
  cudaFree(d_row_ptr);
  cudaFree(d_col_ind);
  cudaFree(d_net_appearence);
  cudaFree(d_ep_starts);
  cudaFree(d_ep_parts);
  cudaFree(d_ep_ends);
  cudaFree(d_src_nodes);
  cudaFree(d_ep_cnt);
  cudaFree(d_gain_loss_1d);
  cudaFree(d_receivers);
  cudaFree(d_senders);
  cudaFree(d_ref_marker_receiver);
  cudaFree(d_ref_marker_sender);
  cudaFree(d_best_parts);
  cudaFree(d_part_sizes);
  cudaFree(d_part_messages);

}
  
