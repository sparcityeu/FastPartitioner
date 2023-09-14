#include <cuda.h>
#include <stdio.h>

#define WARP_COUNT 1024
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

/*
  //CHECK SEG BEFOREHAND
  if(tid == 0){
    printf("Sa, I'm thread %d , I will be your checker\n", tid);
    printf("Nov: %d -- Noe: %d -- Pins: %lld -- K: %d \n", nov, noe, pins, K);
    

    int joker;
    
    for(int i = 0; i < nov; i++){
      joker = d_v_list_push[i];
    }
    printf("d_v_list_push OK! \n");
    
    for(int i = 0; i < nov; i++){
      joker = d_part_vec[i];
    }
    printf("d_part_vec OK! \n");
    
    for(int i = 0; i < noe+1; i++){
      joker = d_row_ptr[i];
    }
    printf("d_row_ptr OK! \n");

    for(int i = 0; i < noe+1; i++){
      joker = d_row_ptr_inv[i];
    }
    printf("d_row_ptr_inv OK! \n");

    for(int i = 0; i < pins; i++){
      joker = d_col_ind[i];
    }
    printf("d_col_ind OK! \n");

    for(int i = 0; i < pins; i++){
      joker = d_col_ind_inv[i];
    }
    printf("d_col_ind_inv OK! \n");

    for(int i = 0; i < noe; i++){
      joker = d_net_appearences[i];
    }
    printf("d_net_appearences OK! \n");

    for(int i = 0; i < noe+1; i++){
      joker = d_ep_starts[i];
    }
    printf("d_ep_starts OK! \n");

    for(int i = 0; i < pins; i++){
      joker = d_ep_parts[i];
    }
    printf("d_ep_parts OK! \n");

    for(int i = 0; i < noe; i++){
      joker = d_ep_ends[i];
    }
    printf("d_ep_ends OK! \n");

    for(int i = 0; i < pins; i++){
      joker = d_ep_cnt[i];
    }
    printf("d_ep_cnt OK! \n");

    for(int i = 0; i < K; i++){
      joker = d_part_vec[i];
    }
    printf("d_part_vec OK! \n");
    
  }
  //CHECK SEG BEFOREHAND
  */

__global__ void gpu_ref_cn_threads(unsigned int* d_row_ptr,
				  unsigned int* d_row_ptr_inv,
				   unsigned int* d_col_ind,
				   unsigned int* d_col_ind_inv,
				   unsigned int* d_net_appearences,
				   int* d_ep_starts,
				   int* d_ep_parts,
				   int* d_ep_ends,
				   int* d_ep_cnt,
				   int* d_part_sizes,
				   int* d_part_vec,
				   int* gloval_d_v_list_push,
				   int* d_noe,
				   int* d_nov,
				   long long* d_pins,
				   int* d_K,
				   int* d_b1,
				   int* d_t1){

  int tid = blockDim.x * blockIdx.x + threadIdx.x;

  int nov = *d_nov;
  int noe = *d_noe;
  long long pins = *d_pins;
  int K = *d_K;

  int b1 = *d_b1;
  int t1 = *d_t1;

  int no_threads = b1*t1;
  int d_v_list_push[32];
  

  if(tid == 0)
    printf("Nov: %d -- Noe: %d -- Pins: %lld -- K: %d \n #Blocks: %d -- #Threads: %d -- #Total: %d\n", nov, noe, pins, K, b1, t1, no_threads);
    
  __syncthreads();
  
  int v_list_cntr;
  int rem_cnt = 0;
  int inc_cnt = 0;

  for(int edge = tid; edge < noe; edge += no_threads){

    //Loop reset for variables
    v_list_cntr = 0; //That is susceptible
    rem_cnt = 0;
    inc_cnt = 0;
    //Loop reset for variables

    if(d_ep_ends[edge] == 2){
      v_list_cntr = 0;
      int part1 = -1, part2 = -1;
      int start_index = d_ep_starts[edge];
      int end_cnt = d_ep_ends[edge];
      //to find p1, pw
      //int part1_count = 0, part2_count = 0;
      int part1_count = 0;
      int part2_count = 0;
      
      for(int j = 0; j < end_cnt; j++){
	if(part1 != -1){
	  part2 = d_ep_parts[start_index + j];
	  part2_count = d_ep_cnt[start_index + j];
	  break;
	}
	else{
	  part1 = d_ep_parts[start_index + j];
	  part1_count = d_ep_cnt[start_index + j];
	}
      }
      
      //now we know p1 and p2

      int old_decision, new_decision;

      if(part1_count < part2_count){
	old_decision = part1;
	new_decision = part2;
      }
      else{
	old_decision = part2;
	new_decision = part1;
      }

      bool not_in_cut = false;

      for(int k = d_row_ptr_inv[edge]; !not_in_cut && k < d_row_ptr_inv[edge + 1]; k++){

	int vertex = d_col_ind_inv[k];
	int part_id = d_part_vec[vertex];

	if(edge == 540327)
	  printf("part_id: %d old_decision: %d \n", part_id, old_decision);
	  
	if(part_id == old_decision){
	  d_v_list_push[v_list_cntr++] = vertex; 
	  
	  for(int j = d_row_ptr[vertex]; !not_in_cut && j < d_row_ptr[vertex + 1]; j++){
	    int net_of_curr_vertex = d_col_ind[j];
	    if(d_ep_ends[net_of_curr_vertex] == 1){ //Which is not in cut
	      not_in_cut = true;
	      break;
	    }
	  }
	}
      }

      if(edge == 540327)
	printf("Edge: 540327 v_list_cntr: %d \n", v_list_cntr);
            
      if(!not_in_cut){
	for(int ind = 0; ind < v_list_cntr; ind++){
	  int v_push = d_v_list_push[ind];
	  int part_id = d_part_vec[v_push];
	  d_part_vec[v_push] = new_decision;
	  d_part_sizes[part_id]--;
	  d_part_sizes[new_decision]++;

	  //REMOVAL
	  for(int i = d_row_ptr[v_push]; i < d_row_ptr[v_push + 1]; i++){
	    int edge = d_col_ind[i];
	    int start_ind = d_ep_starts[edge];
	    int end_cnt = d_ep_ends[edge];

	    for(int j = 0; j < end_cnt; j++){
	      if(d_ep_parts[start_ind + j] == part_id)
		{//We find the part that we need to decrement a connection
		  d_ep_cnt[start_ind + j]--;
		  if(d_ep_cnt[start_ind + j] == 0)
		    {//All connections are removed
		      rem_cnt++;
		      d_ep_parts[start_ind + j] = d_ep_parts[start_ind + end_cnt - 1]; //bring the part in the end to the deleted pos
		      d_ep_cnt[start_ind + j] = d_ep_cnt[start_ind + end_cnt - 1]; //bring the count in the end to the deleted pos
		      d_ep_parts[start_ind + end_cnt - 1] = -1;
		      d_ep_cnt[start_ind + end_cnt - 1] = 0;
		      d_ep_ends[edge]--;
		      end_cnt--;
		    }
		  break;
		}
	    }
	  }
	  for(int i = d_row_ptr[v_push]; i < d_row_ptr[v_push + 1]; i++){//Traversing the edge list of the v_push
	    {
	      int edge = d_col_ind[i];
	      int start_index = d_ep_starts[edge];
	      int end_cnt = d_ep_ends[edge];
	      bool found = false;

	      for(int j = 0; j < end_cnt; j++){
		if(d_ep_parts[start_index + j] == new_decision){//We find the part that we need to increment a connection
		  d_ep_cnt[start_index + j]++;
		  found = true;
		  break;
		}
	      }
	      if(!found){
		d_ep_parts[start_index + end_cnt] = new_decision;
		d_ep_cnt[start_index + end_cnt] = 1;
		d_ep_ends[edge]++;
		inc_cnt++;
	      }
	    }
	  }
	}
      }
    }	    
  }
}
  
extern void partition_gpu_wrapper(unsigned int *row_ptr,
				  unsigned int *col_ind,
				  unsigned int *net_appearences,
				  int nov,
				  int noe,
				  long long pins,
				  int *edge_wgths,
				  int *src_nodes,
				  int K,
				  int *ordering,
				  //string f_name,
				  int REFI_OPT,
				  int REFI_ITER,
				  unsigned int *row_ptr_inv,
				  unsigned int *col_ind_inv,
				  int *part_sizes,
				  int* ep_starts,
				  int* ep_parts,
				  int* ep_ends,
				  int* ep_cnt,
				  int* part_vec){



  printf("SA from the GPU wrapper \n");

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
  b1 = 2048;
  t1 = 512;
  //b1 = 1;
  //t1 = 1;
  int no_threads = b1*t1;
  int no_warps = WARP_COUNT;
  printf("b1: %d -- t1: %d\n", b1, t1);
  

  //DEVICE ALLOCATIONS FOR GLOBAL ARRAYS---------------------------------------------------

  //
  int* v_list_push = new int[nov];
  
  for(int i = 0; i < nov; i++){
    v_list_push[i] = 0;
  }
  
  int* d_v_list_push;
  cudaMalloc(&d_v_list_push, nov*sizeof(int));
  cudaMemcpy(d_v_list_push, v_list_push, nov*sizeof(int), cudaMemcpyHostToDevice);
  //
    
  unsigned int* d_row_ptr;//
  unsigned int* d_col_ind;//
  unsigned int* d_row_ptr_inv;//
  unsigned int* d_col_ind_inv;//
  unsigned int* d_net_appearences;//
  int* d_ep_starts;//
  int* d_ep_parts;//
  int* d_ep_ends;//
  int* d_ep_cnt;//
  int* d_part_sizes;//
  int* d_part_vec;//
  //

  cudaMalloc(&d_row_ptr, (noe+1)*sizeof(unsigned int));
  cudaMalloc(&d_row_ptr_inv, (noe+1)*sizeof(unsigned int));

  cudaMalloc(&d_col_ind, pins*sizeof(unsigned int));
  cudaMalloc(&d_col_ind_inv, pins*sizeof(unsigned int));

  cudaMalloc(&d_net_appearences, noe*sizeof(unsigned int));
  
  cudaMalloc(&d_ep_starts, (noe+1)*sizeof(int));                        
  cudaMalloc(&d_ep_parts, pins*sizeof(int));                          
  cudaMalloc(&d_ep_ends, noe*sizeof(int));                         
  cudaMalloc(&d_ep_cnt, pins*sizeof(int));
  cudaMalloc(&d_part_sizes, K*sizeof(int));
  cudaMalloc(&d_part_vec, nov*sizeof(int));

  //

  cudaMemcpy(d_row_ptr, row_ptr, (noe+1)*sizeof(unsigned int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_row_ptr_inv, row_ptr_inv, (noe+1)*sizeof(unsigned int), cudaMemcpyHostToDevice);

  cudaMemcpy(d_col_ind, col_ind, pins*sizeof(unsigned int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_col_ind_inv, col_ind_inv, pins*sizeof(unsigned int), cudaMemcpyHostToDevice);
  
  cudaMemcpy(d_net_appearences, net_appearences, noe*sizeof(unsigned int), cudaMemcpyHostToDevice);
  
  cudaMemcpy(d_ep_starts, ep_starts, noe*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_ep_parts, ep_parts, pins*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_ep_ends, ep_ends, noe*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_ep_cnt, ep_cnt, pins*sizeof(int), cudaMemcpyHostToDevice);


  cudaMemcpy(d_part_sizes, part_sizes, K*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_part_vec, part_vec, K*sizeof(int), cudaMemcpyHostToDevice);

  //DEVICE ALLOCATIONS FOR GLOBAL ARRAYS---------------------------------------------------

  //DEVICE ALLOCATIONS FOR VARIABLES-------------------------------------------------------
  int* d_nov;
  cudaMalloc((void**)&d_nov, 1*sizeof(int));
  cudaMemcpy(d_nov, &nov, 1*sizeof(int), cudaMemcpyHostToDevice);
 
  int* d_noe;
  cudaMalloc((void**)&d_noe, 1*sizeof(int));
  cudaMemcpy(d_noe, &noe, 1*sizeof(int), cudaMemcpyHostToDevice);

  long long* d_pins;
  cudaMalloc((void**)&d_pins, 1*sizeof(long long));
  cudaMemcpy(d_pins, &pins, 1*sizeof(long long), cudaMemcpyHostToDevice);
  
  int* d_K;
  cudaMalloc((void**)&d_K, 1*sizeof(int));
  cudaMemcpy(d_K, &K, 1*sizeof(int), cudaMemcpyHostToDevice);

  int* d_b1;
  cudaMalloc((void**)&d_b1, 1*sizeof(int));
  cudaMemcpy(d_b1, &b1, 1*sizeof(int), cudaMemcpyHostToDevice);

  int* d_t1;
  cudaMalloc((void**)&d_t1, 1*sizeof(int));
  cudaMemcpy(d_t1, &t1, 1*sizeof(int), cudaMemcpyHostToDevice);

  //for(int i = 0; i < noe; i++){
  //if(ep_ends[i] == 2)
  //printf("Edge w/ 2: %d \n", i );
  //}
  
  for(int i = 0; i < REFI_ITER; i++){
    gpu_ref_cn_threads<<<b1,t1>>>(d_row_ptr,
				  d_row_ptr_inv,
				  d_col_ind,
				  d_col_ind_inv,
				  d_net_appearences,
				  d_ep_starts,
				  d_ep_parts,
				  d_ep_ends,
				  d_ep_cnt,
				  d_part_sizes,
				  d_part_vec,
				  d_v_list_push,
				  d_nov,
				  d_noe,
				  d_pins,
				  d_K,
				  d_b1,
				  d_t1);
      gpuErrchk( cudaDeviceSynchronize() );
  }

  cudaMemcpy(part_vec, d_part_vec, nov*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(ep_starts, d_ep_starts, noe*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(ep_parts, d_ep_parts, pins*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(ep_ends, d_ep_ends, noe*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(ep_cnt, d_ep_cnt, pins*sizeof(int), cudaMemcpyDeviceToHost);

  //for(int i = 0; i < K; i++){
  //printf("Part vec[%d]: %d \n", i, part_vec[i]);
  //}

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
  printf("###########################\n");
}

