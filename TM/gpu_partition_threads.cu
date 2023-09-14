#include <cuda.h>
#include <stdio.h>


#if defined(_WIN64) || defined(__LP64__)
# define PTR_CONSTRAINT "l"
#else
# define PTR_CONSTRAINT "r"
#endif

//Error check-----
#define gpuErrchk3(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
  if (code != cudaSuccess) 
    {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
    }
}
//Error check-----


__global__ void gpu_ref_optimized_threads(int* d_receivers, int* d_senders,
					  int* d_ref_marker_receiver, int* d_ref_marker_sender,
					  int* d_part_vec, int* d_nov,
					  unsigned int* d_row_ptr, unsigned int* d_col_ind,
					  unsigned int* d_net_appearence,int* d_ep_starts, int* d_ep_parts,
					  int* d_ep_ends, int* d_src_nodes, int* d_ep_cnt, int* d_gain_loss_1d,
					  int* d_part_connections, int* d_K, int* d_best_parts,
					  int* d_part_sizes, int* d_part_messages, int* d_no_threads){
  //printf("Nov: %d \n", *d_nov); Referring to int
  
  //int tx = threadIdx.x;
  int tid = blockDim.x * blockIdx.x + threadIdx.x;
  
  int K = *d_K;

  int v;
  int best_part;
  int best_improvement;
  int part_id;
  //These are normally partition-global, but 1M thread atomic update is impossible(?), so they will be vertex-wise(C0)
  int receiver_cnt;
  int sender_cnt;
  int gain;
  int loss;
  //int my_cnt;
  //int tot_improvement;
  
  for(int cursor = tid; cursor < *d_nov; cursor += *d_no_threads){
    
    v = cursor;
    best_part = -1;
    best_improvement = 0;
    part_id = d_part_vec[v];
    
    receiver_cnt = 0;
    sender_cnt = 0;
    gain = 0;
    loss = 0;
    //my_cnt = 0;
    //tot_improvement = 0;
    
    //printf("Thread(%d)  Vertex(%d) \n", tid, v);
    
    int entry = d_row_ptr[v];
    int boundary = d_row_ptr[v+1];
    
    for(int n_ind = entry; n_ind < boundary;n_ind += 1){
      int net = d_col_ind[n_ind];
      int start_index = d_ep_starts[net];
      int end_cnt = d_ep_ends[net];

      if(d_src_nodes[net] == v){
	for(int j = 0; j < end_cnt; j++){
	  int target_part = d_ep_parts[start_index + j];

	  if(target_part < 0){
	    //printf("d_ep_parts[%d]: %d \n", start_index+j, d_ep_parts[start_index+j]);
	    continue;
	  }
	  	  
	  //(EE) This is one of the erroneous lines
	  
	  if(target_part == part_id && d_ep_cnt[start_index + j] > 1 && d_ref_marker_receiver[part_id] != v){
	    //if(d_ref_marker_receiver[part_id] != v){
	      d_receivers[receiver_cnt++] = part_id; //Mind the race condition (C1)
	      d_ref_marker_receiver[part_id] = v;
	      //}
	  }
	  
	  if(target_part != part_id){
	    if(d_ref_marker_receiver[target_part] != v){
	      d_receivers[receiver_cnt++] = target_part;
	      //GTX980 give error at line 95 with k=32
	      d_ref_marker_receiver[target_part] = v;
	    }
	    d_gain_loss_1d[part_id*K+target_part]++;
	    
	    if(d_part_connections[part_id*K+target_part] == 1 ||
	       d_part_connections[part_id*K+target_part] == d_gain_loss_1d[part_id*K+target_part])
	      gain++; 
	    
	  }
	} //end_cnt ends
	
	for(int j = 0; j < end_cnt; j++){
	  int target_part = d_ep_parts[start_index + j];
	  
	  //if(target_part < 0){
	    //printf("d_ep_parts[%d]: %d \n", start_index+j, d_ep_parts[start_index+j]);
	    //continue;
	  //}
	  //(EE)
	  
	  if(target_part != part_id && target_part > 0){
	    d_gain_loss_1d[part_id*K+target_part]--;
	  }
	}
      }
	else{
	  int p_sender = d_part_vec[d_src_nodes[net]];
	  if(d_ref_marker_sender[p_sender] != v){
	    if(sender_cnt < K){
	      d_senders[sender_cnt++] = p_sender;
	    }
	    d_ref_marker_sender[p_sender] = v;
	    int end_cnt = d_ep_ends[net];
	    for(int j = 0; j < end_cnt; j++){
	      if(d_ep_parts[start_index + j] == part_id && d_ep_cnt[start_index + j] == 1 && d_part_connections[p_sender*K+part_id] == 1){
		//if(d_part_connections[p_sender*K+part_id] == 1){
		  gain++;
		  //}
	      }
	    }
	  }
	}
    }
    
    if(gain > 0){
      
      for(int k = 0; k < K; k += 1){
	if(k == part_id)
	  continue;
	
	loss = 0;
	
	for(int i = 0; i < receiver_cnt; i++){
	  int rec_part = d_receivers[i];
	  
	  //if(rec_part < 0){
	  //printf("d_receivers[%d]: %d \n", i, d_receivers[i]);
	  //continue;
	  //}
	  //(EE)
	  
	  //if(rec_part == k || rec_part < 0)
	  //continue;
      
	  
	  if(d_part_connections[(k*K)+rec_part] == 0 && !(rec_part == k || rec_part < 0)){
	    loss++;
	  }
	}
       
	for(int i = 0; i < sender_cnt; i++){
	  int sender_part = d_senders[i];
	  if((sender_part == k || (sender_part*K)+k > (K*K)) && (d_part_connections[(sender_part*K)+k] == 0))
	    loss++;
	}
	
	if((gain > loss) && ((gain - loss) > best_improvement)){
	  best_part = k;
	  best_improvement = gain-loss;
	}
      }
    }
  
    //ADD AND REMOVE
    if(best_part != -1){
      int vertex = v;
      int decision = best_part;
      
      d_part_vec[v] = decision;
      d_part_sizes[part_id]--;
      d_part_sizes[decision]++;
      
      
      //REMOVING
      for(int i = entry; i < boundary; i++){
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
      
      for(int i = entry; i < boundary; i++){
	  int edge = d_col_ind[i];
	  int end_cnt = d_ep_ends[edge];
	  int start_index = d_ep_starts[edge];
	  bool found = false;
	  
	  for(int j = 0; j < end_cnt; j++){
	    int part_idLOCAL = d_ep_parts[start_index + j];
	    
	    if(part_idLOCAL == decision){
	      found = true;
	      d_ep_cnt[start_index + j]++; //CHECKBOI
	      break;
	    }
	  }
	  
	  if(!found){
	    d_ep_parts[start_index + end_cnt] = decision;
	    d_ep_cnt[start_index + end_cnt] = 1;
	    d_ep_ends[edge]++;
	  }
	}
	
	
      for(int n_ind = entry; n_ind < boundary; n_ind++){
	int net = d_col_ind[n_ind];
	int start_index = d_ep_starts[net];
	int end_cnt = d_ep_ends[net];
	
	  if(d_src_nodes[net] == v){
	    for(int j = 0; j < end_cnt; j++){
	      int target_part = d_ep_parts[start_index + j];
	      
	      if(target_part < 0)
		continue;
	      
	      if(target_part != decision){ //CHECK THIS OUT ON OTHERS AS WELL!!! WRONG STATEMENT WERE FOUND
		
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
	    }
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



  
  
extern void partition_gpu_threads_wrapper(unsigned int* row_ptr, unsigned int* col_ind, 
					  unsigned int* net_appearence, int nov, int noe, 
					  long long pins, int* edge_wgths, int* src_nodes, 
					  int K, int* ordering, int REFI_OPT, int REFI_ITER,
					  int* part_vec,
					  int* ep_starts, int* ep_parts, int* ep_ends, int* ep_cnt,
					  int** part_connections, int* part_sizes, int* part_messages){
  
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
  gpuErrchk3( cudaDeviceSynchronize() );
  //Device Attributes-----

  int b1 = no_multiprocessors*32;
  int t1 = max_thread_per_block/2;
  printf("b1: %d -- t1: %d\n", b1, t1);
  int total_threads = b1*t1;

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


  /*
  for(int i = 0; i < K; i++){
    if(part_vec[i] < 0)
      printf("d_part_vec, i: %d, val: %d\n", i, part_vec[i]);
  }

  for(int i = 0; i < noe; i++){
    if(row_ptr[i] < 0)
      printf("d_row_ptr, i: %d, val: %d\n", i, row_ptr[i]);
  }

  for(int i = 0; i < K; i++){
    if(col_ind[i] < 0)
      printf("d_col_ind, i: %d, val: %d\n", i, col_ind[i]);
  }
  */
  
  //for(int i = 0; i < pins; i++){
  //if(ep_parts[i] < 0)
  //printf("ep_parts[i], %d, val: %d\n", i,ep_parts[i]);
  //}

    
  
  int* d_nov;
  cudaMalloc((void**)&d_nov, 1*sizeof(int));
  cudaMemcpy(d_nov, &nov, 1*sizeof(int), cudaMemcpyHostToDevice);

  int* d_K;
  cudaMalloc((void**)&d_K, 1*sizeof(int));
  cudaMemcpy(d_K, &K, 1*sizeof(int), cudaMemcpyHostToDevice);
  
  int* d_no_threads;
  cudaMalloc((void**)&d_no_threads, 1*sizeof(int));
  cudaMemcpy(d_no_threads, &total_threads, 1*sizeof(int), cudaMemcpyHostToDevice);
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
    gpu_ref_optimized_threads<<<b1,t1>>>(d_receivers,d_senders, d_ref_marker_receiver, d_ref_marker_sender, d_part_vec, d_nov, d_row_ptr, d_col_ind, d_net_appearence, d_ep_starts, d_ep_parts, d_ep_ends, d_src_nodes, d_ep_cnt, d_gain_loss_1d, d_part_connections, d_K, d_best_parts, d_part_sizes, d_part_messages, d_no_threads);
    gpuErrchk3( cudaDeviceSynchronize() );
  }
  
  //
  cudaMemcpy(part_vec, d_part_vec, nov*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(ep_starts, d_ep_starts, noe*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(ep_parts, d_ep_parts, pins*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(ep_ends, d_ep_ends, noe*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(src_nodes, d_src_nodes, noe*sizeof(int), cudaMemcpyDeviceToHost);
  //
  
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
  
