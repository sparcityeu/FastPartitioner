#include <cuda.h>
#include <stdio.h>

#define WARP_COUNT 1024


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



gpu_ref_msrv_warps(unsigned int* row_ptr, unsigned int* col_ind,
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
  
  
  //int *gains = new int[K];
  //int *losses = new int[K];
  //int *max_parts = new int[K];
  //int *max_part_marker = new int[K];
  //int *max_part_change_detect = new int [K];
    
  for (iter = 0; REFI_OPT && iter < REFI_ITER; iter++)
    {
        initialize(marker, 0, K, -1);
        initialize(max_part_marker, 0, K, -1);
        int max_all_vol = -1;
        int max_part_count = 0;
        int max_changed_count = 0;
        for (int i = 0; i < K; i++)
        {
            int all_volume = send_volumes[i] + received_volumes[i];
            if (all_volume > max_all_vol)
            {
                max_part_count = 1;
                max_parts[0] = i;
                max_all_vol = all_volume;
                max_changed_count++;
                max_part_marker[i] = max_changed_count;
            }
            else if (all_volume == max_all_vol)
            {
                max_parts[max_part_count++] = i;
                max_part_marker[i] = max_changed_count;
            }
        }
        for (int z = 0; z < nov; z++)
        {
            //GET SCORES STARTS
            int decision;
            int vertex = ordering[z];
            active_part_count = 0;
            int part_id = part_vec[vertex];

            initialize(gains, 0, K, 0);
            initialize(losses, 0, K, 0);

            for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++)
            {
                int edge = col_ind[i];
                int start_index = ep_starts[edge];
                int end_cnt = ep_ends[edge];
                for (int j = 0; j < end_cnt; j++)
                {
                    int score_part_id = ep_parts[start_index + j];
                    if (marker[score_part_id] != vertex)
                    {
                        scores[score_part_id] = 1;
                        marker[score_part_id] = vertex;
                        active_parts[active_part_count] = score_part_id;
                        active_part_count++;
                    }
                    else
                    {
                        scores[score_part_id]++;
                    }
                }
            }

            check_selected_vertex = true;

            //for excluding the selected vertex
            if (check_selected_vertex)
            {
                scores[part_id]--;
            }

            //GET SCORES ENDS
            LDG(check_selected_vertex, K, active_part_count, part_sizes, decision, part_vec, active_parts, C, scores, vertex);
            // IF WE GET DIFFERENT PART
            if (part_id != decision)
            {
                // CHOOSE IF WE WILL GAIN OR LOSE AFTER PLACEMENT
                for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++) // Traversing the edge list of the vertex
                {
                    int edge = col_ind[i];
                    int start_index = ep_starts[edge];
                    int end_cnt = ep_ends[edge];
                    int net_src_node = src_nodes[edge];
                    int net_src_part = part_vec[net_src_node];
                    if (net_src_node == vertex)
                    {
                        gains[part_id] += end_cnt - 1;   // source is changing from the prev part so all conns are assumed to be removed
                        losses[decision] += end_cnt - 1; //new conns
                    }

                    bool found = false;
                    for (int j = 0; j < end_cnt; j++)
                    {
                        if (ep_parts[start_index + j] == part_id && ep_cnt[start_index + j] == 1)
                        { //we find the part that we need to decrement a connection
                            if (net_src_node != vertex)
                            {
                                gains[net_src_part]++; //src will send 1 less connection to target
                                gains[part_id]++;      //target part will receive one less connection
                            }
                            else
                            {
                                losses[decision]--;
                            }
                        }
                        if (ep_parts[start_index + j] == decision)
                        {
                            found = true;
                        }
                    }
                    if (!found)
                    {
                        if (net_src_node != vertex)
                        {
                            losses[net_src_part]++; //source part will send one more
                            losses[decision]++;     //decided part will receive one more connection
                        }
                        else
                        {
                            losses[decision]++;
                        }
                    }
                }

                int best_part = -1;
                bool good_decision = true;
                for (int i = 0; i < max_part_count; i++)
                {
                    int max_part_id = max_parts[i];
                    if (gains[max_part_id] <= losses[max_part_id])
                    {
                        good_decision = false;
                        break;
                    }
                }

                for (int i = 0; good_decision && i < K; i++)
                {
                    if (max_part_marker[i] == max_changed_count)
                    {
                        continue;
                    }

                    if (send_volumes[i] + received_volumes[i] + (losses[i] - gains[i]) >= max_all_vol)
                    {
                        good_decision = false;
                    }
                }

                if (good_decision)
                {
                    best_part = decision;
                }

                // IF WE GAIN AFTER PLACEMENT
                if (best_part != -1)
                {
                    part_vec[vertex] = decision;
                    part_sizes[part_id]--;
                    part_sizes[decision]++;

                    //REMOVING
                    for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++)
                    {
                        int edge = col_ind[i];
                        int net = edge;
                        int start_ind = ep_starts[edge];
                        int end_cnt = ep_ends[edge];
                        int net_src_node = src_nodes[edge];
                        int net_src_part = part_vec[net_src_node];

                        if (net_src_node == vertex)
                        {
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
                                        received_volumes[part_id]--;
                                    }
                                }
                                else
                                {
                                    if (net_src_node == vertex)
                                    {
                                        received_volumes[part_id]++;
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
                                if (net_src_node == vertex)
                                {
                                    received_volumes[decision]--;
                                }
                                break;
                            }
                        }

                        if (!found)
                        {
                            ep_parts[start_index + end_cnt] = decision;
                            ep_cnt[start_index + end_cnt] = 1;
                            ep_ends[edge]++;
                            send_volumes[net_src_part]++;
                            if (net_src_node != vertex)
                            {
                                received_volumes[decision]++;
                            }
                        }
                    }

                    bool max_vol_changed_in_refinement = false;
                    for(int i = 0; i < max_part_count; i++)
                    {
                        int prev_max = max_parts[i];
                        if(max_all_vol != send_volumes[prev_max] + received_volumes[prev_max])
                        {
                            max_vol_changed_in_refinement = true;
                        }
                        max_part_change_detect[i] = max_parts[i];
                    }

                    int prev_max_cnt = max_part_count;
                    for (int i = 0; max_vol_changed_in_refinement && i < K; i++)
                    {
                        if (send_volumes[i] + received_volumes[i] > max_all_vol)
                        {
                            max_part_count = 1;
                            max_parts[0] = i;
                            max_all_vol = send_volumes[i] + received_volumes[i];
                            max_changed_count++;
                            max_part_marker[i] = max_changed_count;
                        }
                        else if (send_volumes[i] + received_volumes[i] == max_all_vol)
                        {
                            max_parts[max_part_count++] = i;
                            max_part_marker[i] = max_changed_count;
                        }
                    }
                }
            }
        }
    }
    
}
  

extern void partition_gpu_wrapper(unsigned int* row_ptr, unsigned int* col_ind, unsigned int* net_appearence, int nov, int noe, long long pins, int* edge_wgths, int* src_nodes, int K, int* ordering, int* marker, int* active_parts, int* part_sizes, int* ep_starts, int* ep_ends, int* ep_parts, int* part_vec, int* ep_cnt, double* scores, int REFI_ITER, double C){  

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
  b1 = 128;
  t1 = 256;
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
  ////
  cudaMalloc(&d_common_mark, K*WARP_COUNT*sizeof(int));
  cudaMalloc(&d_common_active, K*WARP_COUNT*sizeof(int));
  cudaMalloc(&d_common_scores, K*WARP_COUNT*sizeof(double));
  //////
  cudaMalloc(&d_best_parts, nov*sizeof(int));
  
  
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
  

