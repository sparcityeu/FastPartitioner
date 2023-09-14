#include "partition.h"

extern void partition_gpu_wrapper(unsigned int* row_ptr, unsigned int* col_ind, unsigned int* net_appearence, int nov, int noe, long long pins, int* edge_wgths, int* src_nodes, int K, int* ordering, int REFI_OPT, int REFI_ITER, int* part_vec, int* ep_starts, int* ep_parts, int* ep_ends, int* ep_cnt, int** part_connections, int* part_sizes, int* part_messages);//GPU 1

extern void partition_gpu_optimized_wrapper(unsigned int* row_ptr, unsigned int* col_ind, unsigned int* net_appearence, int nov, int noe, long long pins, int* edge_wgths, int* src_nodes, int K, int* ordering, int REFI_OPT, int REFI_ITER, int* part_vec, int* ep_starts, int* ep_parts, int* ep_ends, int* ep_cnt, int** part_connections, int* part_sizes, int* part_messages);//GPU 2

extern void partition_gpu_threads_wrapper(unsigned int* row_ptr, unsigned int* col_ind, unsigned int* net_appearence, int nov, int noe, long long pins, int* edge_wgths, int* src_nodes, int K, int* ordering, int REFI_OPT, int REFI_ITER, int* part_vec, int* ep_starts, int* ep_parts, int* ep_ends, int* ep_cnt, int** part_connections, int* part_sizes, int* part_messages);//GPU 3

//void partition_gpu_optimized(unsigned int *row_ptr, unsigned int *col_ind, unsigned int *net_appearence, int nov, int noe, long long pins, int *edge_wgths, int *src_nodes, int K, int *ordering, string f_name, int REFI_OPT, int REFI_ITER); //This is for testing purposes, just for now I'm calling optimized from base just to see effect

void partition(unsigned int *row_ptr, unsigned int *col_ind, unsigned int *net_appearence, int nov, int noe, long long pins,
               int *edge_wgths, int *src_nodes, int K, int *ordering, string f_name, int REFI_OPT, int REFI_ITER)
{
  double eps = 1.05;          //ratio of balancing
  double C = (nov / K) * eps; //maximum number of vertices that we can assign to single part
  //int active_part_count;
  int src_part_cnt;
  int src_update_cnt;
  int *part_vec = new int[nov];
  int *ep_starts = new int[noe]; // cumulatively stores degree of nets
  int *ep_ends = new int[noe];   // how many parts have vertice from that net
  int *ep_parts = new int[pins]; // stores the parts that nets are connected
  int *marker = new int[K];
  int *src_marker = new int[K];
  int *src_node_parts = new int[noe];
  int *part_sizes = new int[K];
  int *part_messages = new int[K];
  int *active_src; // stores the src node parts while partitioning a vertex
  int *active_src_net_info;
  int active_src_edges_cnt = 0;
  int *src_update;
  int *mes_saved = new int[K];
  int *ep_cnt = new int[pins];
  int **part_connections = new int *[K];
  
  for (int i = 0; i < K; i++)
    {
      part_connections[i] = new int[K];
      for (int j = 0; j < K; j++)
        {
	  part_connections[i][j] = 0;
        }
    }
  
  //vertexi al source ise bunun targetlarının partlarının mesaj sayılarını hesapladım ve buna göre en az mesaj gönderene sourceyi koydum
  //eğer target ise elimdeki vertex capacitye göre bak mesaj olan yere koy
  
  //initializations
  auto start = chrono::steady_clock::now();
  initialize(ep_cnt, 0, pins, 0);
  initialize(marker, 0, K, -1);
  initialize(src_marker, 0, K, -1);
  initialize(part_sizes, 0, K, 0);
  initialize(part_messages, 0, K, 0);
  initialize(mes_saved, 0, K, 0);
  initialize(part_vec, 0, nov, -1);
  initialize(src_node_parts, 0, noe, -1);
  initialize(ep_starts, 0, nov, 0);
  initialize(ep_ends, 0, nov, 0);
  initialize(ep_parts, 0, pins, -1);
  //Calculate edge degrees and max vertex degree
  int max_v_deg = 0;
  for (int i = 0; i < nov; i++)
    {
      if (row_ptr[i + 1] - row_ptr[i] > max_v_deg)
        {
	  max_v_deg = row_ptr[i + 1] - row_ptr[i];
        }
    }
  
  for (int i = 0; i < noe; i++)
    {
      ep_starts[i + 1] = net_appearence[i] + ep_starts[i];
    }
  
  active_src = new int[max_v_deg]; // stores the src node parts while partitioning a vertex
  active_src_net_info = new int[max_v_deg];
  src_update = new int[max_v_deg];
  auto end = chrono::steady_clock::now();
  auto diff = end - start;
  cout << "Initialization took: " << chrono::duration<double, milli>(diff).count() / 1000 << " seconds." << endl;
  total += (chrono::duration<double, milli>(end - start));
  //Start Partitioning
  cout << "Partitioning started" << endl;
  start = chrono::steady_clock::now();
  
  for (int z = 0; z < nov; z++)
    {
      int vertex = ordering[z]; // current vertex we want to add to our partition
      int decision;             // it will store the part id of the decision we made for the vertex part
      src_part_cnt = 0;         // number of different parts that the current net set has it's source node in it
      src_update_cnt = 0;       // number of different nets that they have the <vertex> as their source node
      active_src_edges_cnt = 0;
      for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++)
        {
	  int edge = col_ind[i];
	  if (src_nodes[edge] == vertex)
            {
	      src_update[src_update_cnt++] = edge;
            }
	  
	  if (src_node_parts[edge] != -1)
            { // src node is assigned
	      if (src_marker[src_node_parts[edge]] != vertex)
                {
		  src_marker[src_node_parts[edge]] = vertex;
		  active_src[src_part_cnt++] = src_node_parts[edge];
                }
	      active_src_net_info[active_src_edges_cnt++] = edge;
            }
        }
      
      int max_mes_saved = -1; // maximum number of messages saved for a potential decision on a part <k>
      int max_index = -1;     // part id of the maximum message saved part
      
      for (int k = 0; k < K; k++)
        {
	  if (part_sizes[k] < C)
            { // only if the capacity allows
	      int mes_saved = 0;
	      for (int i = 0; i < src_part_cnt; i++)
                {
		  int src_part = active_src[i];
		  if (part_connections[src_part][k] > 0)
                    {
		      mes_saved++;
                    }
                }
	      for (int i = 0; i < src_update_cnt; i++)
                {
		  int edge = src_update[i];
		  int end_cnt = ep_ends[edge];
		  int start_index = ep_starts[edge];
		  for (int j = 0; j < end_cnt; j++)
                    {
		      int part_id = ep_parts[start_index + j];
		      if (part_connections[k][part_id] > 0)
                        {
			  mes_saved++;
                        }
                    }
                }
	      if (mes_saved > max_mes_saved || (mes_saved == max_mes_saved &&
						(part_sizes[k] < part_sizes[max_index])))
                {
		  max_index = k;
		  max_mes_saved = mes_saved;
                }
            }
        }
      
      if (max_index == -1)
        {
	  int minLoad = INT_MAX;
	  int minLoadIndex = -1;
	  for (int i = 0; i < K; i++)
            {
	      if (part_sizes[i] < minLoad)
                {
		  minLoad = part_sizes[i];
		  minLoadIndex = i;
                }
            }
	  max_index = minLoadIndex;
        }
      decision = max_index;
      part_vec[vertex] = decision; //vertex part is decided
      part_sizes[decision]++;
      
      //Update the edge part info
      for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++)
        {
	  int edge = col_ind[i];
	  int end_cnt = ep_ends[edge];
	  int start_index = ep_starts[edge];
	  bool found = false;
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
            }
        }
      //update messages and source parts
      for (int i = 0; i < active_src_edges_cnt; i++)
        {
	  int edge = active_src_net_info[i];
	  int src_part = src_node_parts[edge];
	  
	  if (src_part == decision)
	    continue;
	  if (part_connections[src_part][decision] == 0)
            {
	      part_messages[src_part]++;
	      part_connections[src_part][decision] = 1;
            }
	  else
            {
	      int end_cnt = ep_ends[edge];
	      int start_index = ep_starts[edge];
	      for (int j = 0; j < end_cnt; j++)
                {
		  int part_id = ep_parts[start_index + j];
		  if (part_id == decision)
                    {
		      if (ep_cnt[start_index + j] == 1) //if net caused a new connection
                        {
			  part_connections[src_part][decision]++;
                        }
		      break;
                    }
                }
            }
        }
      
      for (int i = 0; i < src_update_cnt; i++)
        {
	  int edge = src_update[i];
	  src_node_parts[edge] = decision;
	  int end_cnt = ep_ends[edge];
	  int start_index = ep_starts[edge];
	  
            for (int j = 0; j < end_cnt; j++)
	      {
                int part_id = ep_parts[start_index + j];
                if (part_id == decision)
		  continue;
                if (part_connections[decision][part_id] == 0)
		  {
                    part_messages[decision]++;
                    part_connections[decision][part_id] = 1;
		  }
                else
		  {
                    part_connections[decision][part_id]++;
		  }
	      }
        }
    }
  
  end = chrono::steady_clock::now();
  diff = end - start;
  total += (chrono::duration<double, milli>(end - start));
  cout << "Paritioning took : " << chrono::duration<double, milli>(total).count() / 1000 << " seconds." << endl;
  
  int TV_OUT = TV(ep_ends, noe);
  int CN_OUT = CN(ep_ends, noe);
  int TM_OUT = TM_scratch(K, part_vec, ep_parts, ep_starts, ep_ends, nov, noe, src_nodes);
  int MGM_OUT = MGM_scratch(K, part_vec, ep_parts, ep_starts, ep_ends, nov, noe, src_nodes);
  int MGV_OUT = MGV_Scratch(part_vec, row_ptr, col_ind, src_nodes, ep_ends, net_appearence, K, nov, noe, pins);
  int MGAV_OUT = MGAV_Scratch(part_vec, row_ptr, col_ind, src_nodes, ep_ends, ep_parts, ep_starts, net_appearence, K, nov, noe, pins);
  
  cout << "TV : " << TV_OUT << endl;
  cout << "CN : " << CN_OUT << endl;
  cout << "TM : " << TM_OUT << endl;
  cout << "MGM : " << MGM_OUT << endl;
  cout << "MGV : " << MGV_OUT << endl;
  cout << "MGAV : " << MGAV_OUT << endl;


      double avg_weigth2 = 0.0;
      for(int i = 0; i < K; i++){
	avg_weigth2 += part_sizes[i];
      }
      avg_weigth2 /= K;

      double a_imbal2 = 0;
      for(int i = 0; i < K; i++) {
	a_imbal2 = max(a_imbal2, part_sizes[i]/avg_weigth2);
      }

      printf("#############################\n");
      printf("SEQ BEFORE REFINEMENT ACTUAL IMBAL: %f \n", a_imbal2);
      printf("#############################\n");

  bool once = true;
  start = chrono::steady_clock::now();
  //REFINEMENT CODE STARTS
  int ** gain_loss_2d = new int * [K];
  for(int i = 0; i < K; i++)
    {
      gain_loss_2d[i] = new int [K];
      for(int j = 0; j < K; j++)
        {
	  gain_loss_2d[i][j] = 0;
        }
    }

  for (int iter = 0; REFI_OPT && iter < REFI_ITER; iter++)
    {
      auto bizim_start = chrono::steady_clock::now();
      int *receivers = new int[K];
      int *senders = new int[K];
      
      int *ref_marker_receiver = new int[K];
      int *ref_marker_sender = new int[K];
      initialize(ref_marker_receiver, 0, K, -1);
      initialize(ref_marker_sender, 0, K, -1);
      
      int receiver_cnt, sender_cnt, gain, loss, my_cnt = 0, tot_improvement = 0;
      for (int v = 0; v < nov; v++) //This is for warps
        {
	  int best_part = -1, best_improvement = 0;
	  int part_id = part_vec[v];
	  gain = 0;
	  receiver_cnt = 0;
	  sender_cnt = 0;
	  
	  for (int n_ind = row_ptr[v]; n_ind < row_ptr[v + 1]; n_ind++) //This is for threads
	    {
	      int net = col_ind[n_ind];
	      int start_index = ep_starts[net];
	      int end_cnt = ep_ends[net];
	      if (src_nodes[net] == v)
		{
		  for (int j = 0; j < end_cnt; j++)
                    {
		      int target_part = ep_parts[start_index + j];
		      if (target_part == part_id && ep_cnt[start_index + j] > 1) //if there is still a connection after moving vertex v
                        {
			  if (ref_marker_receiver[part_id] != v)
                            {
			      receivers[receiver_cnt++] = part_id;
			      ref_marker_receiver[part_id] = v;
                            }
                        }
		      
		      if (target_part != part_id)
                        {
			  if(ref_marker_receiver[target_part] != v)
                            {
			      receivers[receiver_cnt++] = target_part;
			      ref_marker_receiver[target_part] = v;
                            }
			  
			  gain_loss_2d[part_id][target_part]++;
			  
			  if (part_connections[part_id][target_part] == 1 ||
			      part_connections[part_id][target_part] == gain_loss_2d[part_id][target_part])
                            {
			      gain++;
                            }
                        }
                    }
		  
		  
		  //BURAYA KADAR Bİ GİRELİM

		  //std::cout << "vertex: " << v << " end_cnt: " << end_cnt << std::endl;
		  //reset gain loss 2d
		  for (int j = 0; j < end_cnt; j++)
                    {
		      int target_part = ep_parts[start_index + j];
		      if (target_part != part_id)
                        {
			  gain_loss_2d[part_id][target_part]--;
                        }
                    }
		}
		  
		  else
		    {
		      int p_sender = part_vec[src_nodes[net]];
		      //cout <<"v: " << v << " p_sender: " << p_sender << endl;
		      if (ref_marker_sender[p_sender] != v)
			{
			  senders[sender_cnt++] = p_sender;
			  ref_marker_sender[p_sender] = v;
			  int end_cnt = ep_ends[net];
			  for (int j = 0; j < end_cnt; j++)
			    {
			      if (ep_parts[start_index + j] == part_id && ep_cnt[start_index + j] == 1)
				{
				  if (part_connections[p_sender][part_id] == 1)
				    {
				      gain++;
				    }
				}
			    }
			}
		    }
		}
	      
	  //if(gain > 0)
	  //std::cout << "Gain for v: " << v << " -> " << gain << std::endl;

	  /*
	  if(once){
	    cout << "v: " << v << endl;
	    cout << "Receiver_cnt: " << receiver_cnt << endl;
	    cout << "Sender_cnt: " << sender_cnt << endl;
	    //once = false;
	  }
	  */
	
	  //FIND BETTER PART
	  if (gain > 0)
            {
	      for (int k = 0; k < K; k++)
                {
		  if (part_sizes[k] >= C) {
		    continue;
		  }
		  
		  if (k == part_id)
		    continue;
		  
		  loss = 0;
		  
		  for (int i = 0; i < receiver_cnt; i++)
                    {
		      int rec_part = receivers[i];
		      if (rec_part == k)
			continue;
		      if (part_connections[k][rec_part] == 0)
                        {
			  loss++;
                        }
                    }
		  
		  for (int i = 0; i < sender_cnt; i++)
                    {
		      int sender_part = senders[i];
		      if (sender_part == k)
			continue;
		      if (part_connections[sender_part][k] == 0)
                        {
			  loss++;
                        }
                    }
		  
		  if ((gain > loss) && ((gain - loss) > best_improvement))
                    {
		      best_part = k;
		      best_improvement = gain - loss;
		      //cout << "gain : " << gain << " loss : " << loss << " v : " << v <<  endl;
                    }
                }
            }

	  //if(gain > 0)
	  //std::cout << "Vertex: " << v << " best_part -> " << best_part << std::endl;
	  
	  //ADD AND REMOVE
	  int added_mes = 0, removed_mes = 0;
	  if (best_part != -1)
            {
	      int vertex = v;
	      int decision = best_part;
	      part_vec[v] = decision;
	      part_sizes[part_id]--;
	      part_sizes[decision]++;
	      
	      //REMOVING
	      for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++)
                {
		  int edge = col_ind[i];
		  int net = edge;
		  int start_ind = ep_starts[edge];
		  int end_cnt = ep_ends[edge];
		  for (int j = 0; j < end_cnt; j++)
                    {
		      if (ep_parts[start_ind + j] == part_id) //-------///
                        { //we find the part that we need to decrement a connection
			  ep_cnt[start_ind + j]--;
			  if (ep_cnt[start_ind + j] == 0)
                            { //all connections are removed
			      ep_parts[start_ind + j] = ep_parts[start_ind + end_cnt - 1]; //bring the part in the end to the deleted pos
			      ep_cnt[start_ind + j] = ep_cnt[start_ind + end_cnt - 1];     //bring the count in the end to the deleted pos
			      ep_parts[start_ind + end_cnt - 1] = -1;
			      ep_cnt[start_ind + end_cnt - 1] = 0;
			      ep_ends[edge]--;
                            }
			  break; // we don't need to search no longer
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
		      }
                }
	      
	      
	      for (int n_ind = row_ptr[v]; n_ind < row_ptr[v + 1]; n_ind++)
                {
		  int net = col_ind[n_ind];
		  int start_index = ep_starts[net];
		  int end_cnt = ep_ends[net];
		  if (src_nodes[net] == v)
                    {
		      for (int j = 0; j < end_cnt; j++)
                        {
			  int target_part = ep_parts[start_index + j];

			  if (target_part != decision) //OK!
                            { //we find the part that we need to increment a connection
			      
			      part_connections[decision][target_part]++;
			      if (part_connections[decision][target_part] == 1)
                                {
				  part_messages[decision]++;
                                }
                            }
			  
			  if (target_part != part_id && target_part != decision) //OK!
                            {
			      part_connections[part_id][target_part]--;
			      if (part_connections[part_id][target_part] == 0) //OK!
                                {
				  part_messages[part_id]--;
                                }
                            }

			  if (target_part == decision && ep_cnt[start_index + j] > 1)
                            {
			      part_connections[part_id][target_part]--;
                                if (part_connections[part_id][target_part] == 0)
				  {
                                    part_messages[part_id]--;
				  }
                            }
                        }
                    }
		  else
                    {
		      int src_net = part_vec[src_nodes[net]];
		      int end_cnt = ep_ends[net];
		      bool found = false;
		      for (int j = 0; j < end_cnt; j++)
                        {
			  if (ep_parts[start_index + j] == decision && ep_cnt[start_index + j] == 1)
                            { //we find the part that we need to decrement a connection
			      part_connections[src_net][decision]++;
			      
			      if (part_connections[src_net][decision] == 1)
				{
				  part_messages[src_net]++;
                                }
                            }
			  if (ep_parts[start_index + j] == part_id)
			    found = true;
                        }
		      if (!found)
                        {
			  part_connections[src_net][part_id]--;
			  if (part_connections[src_net][part_id] == 0)
                            {
			      part_messages[src_net]--;
                            }
                        }
                    }
		  
                }
	      //if (best_improvement != (removed_mes - added_mes))
	      //cout << "v :" << v << " Best improvement  " << best_improvement << "Removed - Added :" << removed_mes - added_mes << endl;
            }
        }
    }
  
  if (REFI_OPT)
    {
      end = chrono::steady_clock::now();
      diff = end - start;
      cout << "Refinement took : " << chrono::duration<double, milli>(diff).count() / 1000 << " seconds." << endl;
      total += (chrono::duration<double, milli>(end - start));
      cout << "Total Paritioning took : " << chrono::duration<double, milli>(total).count() / 1000 << " seconds." << endl;
      
      int TV_Refined_out = TV_Refined(ep_starts, ep_cnt, noe, K);
      int CN_Refined_out = CN(ep_ends, noe);
      int TM_Refined_out = TM_scratch(K, part_vec, ep_parts, ep_starts, ep_ends, nov, noe, src_nodes);
      int MGM_Refined_out = MGM_scratch(K, part_vec, ep_parts, ep_starts, ep_ends, nov, noe, src_nodes);
      int MGV_Refined_out = MGV_Scratch(part_vec, row_ptr, col_ind, src_nodes, ep_ends, net_appearence, K, nov, noe, pins);
      int MGAV_Refined_out = MGAV_Scratch(part_vec, row_ptr, col_ind, src_nodes, ep_ends, ep_parts, ep_starts, net_appearence, K, nov, noe, pins);

      /**/
      //for(int i = 0; i < pins; i++){
      //if(ep_parts[i] < 0)
      //printf("NAIVE ep_parts[i], %d, val: %d\n", i,ep_parts[i]);
      //}
      /**/
      
      
      if (REFI_OPT)
        {
	  cout << "After Refinement TV : " << TV_Refined_out << endl;
	  cout << "After Refinement CN : " << CN_Refined_out << endl;
	  cout << "After Refinement TM : " << TM_Refined_out << endl;
	  cout << "After Refinement MGM : " << MGM_Refined_out << endl;
	  cout << "After Refinement MGV : " << MGV_Refined_out << endl;
	  cout << "After Refinement MGAV : " << MGAV_Refined_out << endl;
	  cout << "Improved Percentage TV : " << 100.0 - ((double(TV_Refined_out) / TV_OUT) * 100) << endl;
	  cout << "Improved Percentage CN : " << 100.0 - ((double(CN_Refined_out) / CN_OUT) * 100) << endl;
	  cout << "Improved Percentage TM : " << 100.0 - ((double(TM_Refined_out) / TM_OUT) * 100) << endl;
	  cout << "Improved Percentage MGM : " << 100.0 - ((double(MGM_Refined_out) / MGM_OUT) * 100) << endl;
	  cout << "Improved Percentage MGV : " << 100.0 - ((double(MGV_Refined_out) / MGV_OUT) * 100) << endl;
	  cout << "Improved Percentage MGAV : " << 100.0 - ((double(MGAV_Refined_out) / MGAV_OUT) * 100) << endl;
        }
      else
        {
          cout << "After Refinement TV : " << -1 << endl;
          cout << "After Refinement CN : " << -1 << endl;
          cout << "After Refinement TM : " << -1 << endl;
          cout << "After Refinement MGM : " << -1 << endl;
          cout << "After Refinement MGV : " << -1 << endl;
          cout << "After Refinement MGAV : " << -1 << endl;
          cout << "Improved Percentage TV : " << -1 << endl;
          cout << "Improved Percentage CN : " << -1 << endl;
          cout << "Improved Percentage TM : " << -1<< endl;
          cout << "Improved Percentage MGM : " << -1 << endl;
          cout << "Improved Percentage MGV : " << -1 << endl;
          cout << "Improved Percentage MGAV : " << -1 << endl;
        }
      
      cout << "Iterated for : " << REFI_ITER << " times\n";

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
      printf("SEQ BEFORE REFINEMENT ACTUAL IMBAL: %f \n", a_imbal);
      printf("#############################\n");
    }
}

void partition_par(unsigned int *row_ptr, unsigned int *col_ind, unsigned int *net_appearence, int nov, int noe, long long pins, int *edge_wgths, int *src_nodes, int K, int *ordering, string f_name, int REFI_OPT, int REFI_ITER, int NT)
{
  double eps = 1.05;          //ratio of balancing
  double C = (nov / K) * eps; //maximum number of vertices that we can assign to single part
  //int active_part_count;
  int src_part_cnt;
  int src_update_cnt;
  int *part_vec = new int[nov];
  int *ep_starts = new int[noe]; // cumulatively stores degree of nets
  int *ep_ends = new int[noe];   // how many parts have vertice from that net
  int *ep_parts = new int[pins]; // stores the parts that nets are connected
  int *marker = new int[K];
  int *src_marker = new int[K];
  int *src_node_parts = new int[noe];
  int *part_sizes = new int[K];
  int *part_messages = new int[K];
  int *active_src; // stores the src node parts while partitioning a vertex
  int *active_src_net_info;
  int active_src_edges_cnt = 0;
  int *src_update;
  int *mes_saved = new int[K];
  int *ep_cnt = new int[pins];
  int **part_connections = new int *[K];

  for (int i = 0; i < K; i++)
    {
      part_connections[i] = new int[K];
      for (int j = 0; j < K; j++)
        {
	  part_connections[i][j] = 0;
        }
    }
  
  //vertexi al source ise bunun targetlarının partlarının mesaj sayılarını hesapladım ve buna göre en az mesaj gönderene sourceyi koydum
  //eğer target ise elimdeki vertex capacitye göre bak mesaj olan yere koy
  
  //initializations
  auto start = chrono::steady_clock::now();
  initialize(ep_cnt, 0, pins, 0);
  initialize(marker, 0, K, -1);
  initialize(src_marker, 0, K, -1);
  initialize(part_sizes, 0, K, 0);
  initialize(part_messages, 0, K, 0);
  initialize(mes_saved, 0, K, 0);
  initialize(part_vec, 0, nov, -1);
  initialize(src_node_parts, 0, noe, -1);
  initialize(ep_starts, 0, nov, 0);
  initialize(ep_ends, 0, nov, 0);
  initialize(ep_parts, 0, pins, -1);
  //Calculate edge degrees and max vertex degree
  int max_v_deg = 0;
  for (int i = 0; i < nov; i++)
    {
      if (row_ptr[i + 1] - row_ptr[i] > max_v_deg)
        {
	  max_v_deg = row_ptr[i + 1] - row_ptr[i];
        }
    }
  
  for (int i = 0; i < noe; i++)
    {
      ep_starts[i + 1] = net_appearence[i] + ep_starts[i];
    }
  
  active_src = new int[max_v_deg]; // stores the src node parts while partitioning a vertex
  active_src_net_info = new int[max_v_deg];
  src_update = new int[max_v_deg];
  auto end = chrono::steady_clock::now();
  auto diff = end - start;
  cout << "Initialization took: " << chrono::duration<double, milli>(diff).count() / 1000 << " seconds." << endl;
  total += (chrono::duration<double, milli>(end - start));
  //Start Partitioning
  cout << "Partitioning started" << endl;
  start = chrono::steady_clock::now();
  
  for (int z = 0; z < nov; z++)
    {
      int vertex = ordering[z]; // current vertex we want to add to our partition
      int decision;             // it will store the part id of the decision we made for the vertex part
      src_part_cnt = 0;         // number of different parts that the current net set has it's source node in it
      src_update_cnt = 0;       // number of different nets that they have the <vertex> as their source node
      active_src_edges_cnt = 0;
      for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++)
        {
	  int edge = col_ind[i];
	  if (src_nodes[edge] == vertex)
            {
	      src_update[src_update_cnt++] = edge;
            }
	  
	  if (src_node_parts[edge] != -1)
            { // src node is assigned
	      if (src_marker[src_node_parts[edge]] != vertex)
                {
		  src_marker[src_node_parts[edge]] = vertex;
		  active_src[src_part_cnt++] = src_node_parts[edge];
                }
	      active_src_net_info[active_src_edges_cnt++] = edge;
            }
        }
      
      for (int k = 0; k < K; k++)
        {
	  int mes_saved_loc = 0;
	  if (part_sizes[k] < C)
            { // only if the capacity allows
	      
	      for (int i = 0; i < src_part_cnt; i++)
                {
		  int src_part = active_src[i];
		  if (part_connections[src_part][k] > 0)
                    {
		      mes_saved_loc++;
                    }
                }
	      for (int i = 0; i < src_update_cnt; i++)
                {
		  int edge = src_update[i];
		  int end_cnt = ep_ends[edge];
		  int start_index = ep_starts[edge];
		  for (int j = 0; j < end_cnt; j++)
                    {
		      int part_id = ep_parts[start_index + j];
		      if (part_connections[k][part_id] > 0)
                        {
			  mes_saved_loc++;
                        }
                    }
                }
            }
	  mes_saved[k] = mes_saved_loc;
        }
      
      int max_mes_saved = -1; // maximum number of messages saved for a potential decision on a part <k>
      int max_index = -1;     // part id of the maximum message saved part
      
      for (int k = 0; k < K; k++)
        {
	  if (mes_saved[k] > max_mes_saved || (mes_saved[k] == max_mes_saved &&
					       (part_sizes[k] < part_sizes[max_index])))
            {
	      max_index = k;
	      max_mes_saved = mes_saved[k];
            }
        }
      
      if (max_index == -1)
        {
	  int minLoad = INT_MAX;
	  int minLoadIndex = -1;
	  for (int i = 0; i < K; i++)
            {
	      if (part_sizes[i] < minLoad)
                {
		  minLoad = part_sizes[i];
		  minLoadIndex = i;
                }
            }
	  max_index = minLoadIndex;
        }
      decision = max_index;
      part_vec[vertex] = decision; //vertex part is decided
      part_sizes[decision]++;
      
      //Update the edge part info
      for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++)
        {
	  int edge = col_ind[i];
	  int end_cnt = ep_ends[edge];
	  int start_index = ep_starts[edge];
	  bool found = false;
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
            }
        }
      //update messages and source parts
      for (int i = 0; i < active_src_edges_cnt; i++)
        {
	  int edge = active_src_net_info[i];
	  int src_part = src_node_parts[edge];
	  
	  if (src_part == decision)
	    continue;
	  if (part_connections[src_part][decision] == 0)
            {
	      part_messages[src_part]++;
	      part_connections[src_part][decision] = 1;
            }
	  else
            {
	      int end_cnt = ep_ends[edge];
	      int start_index = ep_starts[edge];
	      for (int j = 0; j < end_cnt; j++)
                {
		  int part_id = ep_parts[start_index + j];
		  if (part_id == decision)
                    {
		      if (ep_cnt[start_index + j] == 1) //if net caused a new connection
                        {
			  part_connections[src_part][decision]++;
                        }
		      break;
                    }
                }
            }
        }
      
      for (int i = 0; i < src_update_cnt; i++)
        {
	  int edge = src_update[i];
	  src_node_parts[edge] = decision;
	  int end_cnt = ep_ends[edge];
	  int start_index = ep_starts[edge];
	  
	  for (int j = 0; j < end_cnt; j++)
            {
	      int part_id = ep_parts[start_index + j];
	      if (part_id == decision)
		continue;
	      if (part_connections[decision][part_id] == 0)
                {
		  part_messages[decision]++;
		  part_connections[decision][part_id] = 1;
                }
	      else
                {
		  part_connections[decision][part_id]++;
                }
            }
        }
    }
  
  end = chrono::steady_clock::now();
  diff = end - start;
  total += (chrono::duration<double, milli>(end - start));
  cout << "Paritioning took : " << chrono::duration<double, milli>(total).count() / 1000 << " seconds." << endl;
  
  int TV_OUT = TV(ep_ends, noe);
  int CN_OUT = CN(ep_ends, noe);
  int TM_OUT = TM_scratch(K, part_vec, ep_parts, ep_starts, ep_ends, nov, noe, src_nodes);
  int MGM_OUT = MGM_scratch(K, part_vec, ep_parts, ep_starts, ep_ends, nov, noe, src_nodes);
  int MGV_OUT = MGV_Scratch(part_vec, row_ptr, col_ind, src_nodes, ep_ends, net_appearence, K, nov, noe, pins);
  int MGAV_OUT = MGAV_Scratch(part_vec, row_ptr, col_ind, src_nodes, ep_ends, ep_parts, ep_starts, net_appearence, K, nov, noe, pins);
  
  cout << "TV : " << TV_OUT << endl;
  cout << "CN : " << CN_OUT << endl;
  cout << "TM : " << TM_OUT << endl;
  cout << "MGM : " << MGM_OUT << endl;
  cout << "MGV : " << MGV_OUT << endl;
  cout << "MGAV : " << MGAV_OUT << endl;

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
      printf("SEQ BEFORE REFINEMENT ACTUAL IMBAL: %f \n", a_imbal);
      printf("#############################\n");
  //REFINEMENT CODE STARTS
  start = chrono::steady_clock::now();
  int **vertices_par = new int *[NT];
  int **decisions_par = new int *[NT];
  int *refined_cnt_glob = new int[NT];
  
  for (int t = 0; t < NT; t++)
    {
      vertices_par[t] = new int[nov];
      decisions_par[t] = new int[nov];
    }
  
  omp_set_dynamic(0);
  omp_set_num_threads(NT);
  bool max_part_changed = true;//a flag to store max part change info
  
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    int *receivers = new int[K];
    int *senders = new int[K];
    int *ref_marker_receiver = new int[K];
    int *ref_marker_sender = new int[K];
    int ** gain_loss_2d = new int * [K];
    
    for(int i = 0; i < K; i++)
      {
	gain_loss_2d[i] = new int [K];
	for(int j = 0; j < K; j++)
	  {
	    gain_loss_2d[i][j] = 0;
	  }
      }
    for (int iter = 0; REFI_OPT && iter < REFI_ITER; iter++)
      {
	
	initialize(ref_marker_receiver, 0, K, -1);
	initialize(ref_marker_sender, 0, K, -1);
	refined_cnt_glob[tid] = 0;
	int refined_cnt_loc = 0;
	int receiver_cnt, sender_cnt, gain, loss;
	
#pragma omp for schedule(dynamic, 256)
	for (int v = 0; v < nov; v++)
	  {
	    int best_part = -1, best_improvement = 0;
	    int part_id = part_vec[v];
	    gain = 0;
	    receiver_cnt = 0;
	    sender_cnt = 0;
	    
	    for (int n_ind = row_ptr[v]; n_ind < row_ptr[v + 1]; n_ind++)
	      {
		int net = col_ind[n_ind];
		int start_index = ep_starts[net];
		int end_cnt = ep_ends[net];
		if (src_nodes[net] == v)
		  {
		    for (int j = 0; j < end_cnt; j++)
		      {
			int target_part = ep_parts[start_index + j];
			if (target_part == part_id && ep_cnt[start_index + j] > 1) //if there is still a connection after moving vertex v
			  {
			    if (ref_marker_receiver[part_id] != v)
			      {
				receivers[receiver_cnt++] = part_id;
				ref_marker_receiver[part_id] = v;
			      }
			  }
			
			if (target_part != part_id)
			  {
			    if(ref_marker_receiver[target_part] != v)
			      {
				receivers[receiver_cnt++] = target_part;
				ref_marker_receiver[target_part] = v;
			      }
			    
			    gain_loss_2d[part_id][target_part]++;
			    
			    if (part_connections[part_id][target_part] == 1 ||
				part_connections[part_id][target_part] == gain_loss_2d[part_id][target_part])
			      {
				gain++;
			      }
			  }
		      }
		    
		    //reset gain loss 2d
		    for (int j = 0; j < end_cnt; j++)
		      {
			int target_part = ep_parts[start_index + j];
			if (target_part != part_id)
			  {
			    gain_loss_2d[part_id][target_part]--;
			  }
		      }
		  }
		else
		  {
		    int p_sender = part_vec[src_nodes[net]];
		    if (ref_marker_sender[p_sender] != v)
		      {
			senders[sender_cnt++] = p_sender;
			ref_marker_sender[p_sender] = v;
			int end_cnt = ep_ends[net];
			for (int j = 0; j < end_cnt; j++)
			  {
			    if (ep_parts[start_index + j] == part_id && ep_cnt[start_index + j] == 1)
			      {
				if (part_connections[p_sender][part_id] == 1)
				  {
				    gain++;
				  }
			      }
			  }
		      }
		  }
	      }
	    
	    //FIND BETTER PART
	    if (gain > 0)
	      {
		for (int k = 0; k < K; k++)
		  {
		    if (k == part_id)
		      continue;

		    if(part_sizes[k] >= C) {
		      continue;
		    }
		    
		    loss = 0;
		    
		    for (int i = 0; i < receiver_cnt; i++)
		      {
			int rec_part = receivers[i];
			if (rec_part == k)
			  continue;
			if (part_connections[k][rec_part] == 0)
			  {
			    loss++;
			  }
		      }
		    
		    for (int i = 0; i < sender_cnt; i++)
		      {
			int sender_part = senders[i];
			if (sender_part == k)
			  continue;
			if (part_connections[sender_part][k] == 0)
			  {
			    loss++;
			  }
		      }
		    
		    if ((gain > loss) && ((gain - loss) > best_improvement))
		      {
			best_part = k;
			best_improvement = gain - loss;
		      }
		  }
	      }
	    
	    if(best_part != -1)
	      {
		vertices_par[tid][refined_cnt_loc] = v;
		decisions_par[tid][refined_cnt_loc] = best_part;
		refined_cnt_loc++;
	      }
	  }
	refined_cnt_glob[tid] = refined_cnt_loc;
	
#pragma omp single
	for (int i = 0; i < NT; i++)
	  {
	    for (int k = 0; k < refined_cnt_glob[i]; k++)
	      {
		int vertex = vertices_par[i][k];
		int v = vertex;
		int decision = decisions_par[i][k];
		if(part_sizes[decision] >= C * (1.02)) {
		  continue;
		}
		
		int part_id = part_vec[vertex];
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
		    for (int j = 0; j < end_cnt; j++)
		      {
			if (ep_parts[start_ind + j] == part_id)
			  { //we find the part that we need to decrement a connection
			    ep_cnt[start_ind + j]--;
			    if (ep_cnt[start_ind + j] == 0)
			      { //all connections are removed
				ep_parts[start_ind + j] = ep_parts[start_ind + end_cnt - 1]; //bring the part in the end to the deleted pos
				ep_cnt[start_ind + j] = ep_cnt[start_ind + end_cnt - 1];     //bring the count in the end to the deleted pos
				ep_parts[start_ind + end_cnt - 1] = -1;
				ep_cnt[start_ind + end_cnt - 1] = 0;
				ep_ends[edge]--;
			      }
			    break; // we don't need to search no longer
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
		      }
		  }
		
		
		for (int n_ind = row_ptr[v]; n_ind < row_ptr[v + 1]; n_ind++)
		  {
		    int net = col_ind[n_ind];
		    int start_index = ep_starts[net];
		    int end_cnt = ep_ends[net];
		    if (src_nodes[net] == v)
		      {
			for (int j = 0; j < end_cnt; j++)
			  {
			    int target_part = ep_parts[start_index + j];
			    
			    if (target_part != decision)
			      { //we find the part that we need to increment a connection
				
				part_connections[decision][target_part]++;
				if (part_connections[decision][target_part] == 1)
				  {
				    part_messages[decision]++;
				  }
			      }
			    
			    if (target_part != part_id && target_part != decision)
			      {
				part_connections[part_id][target_part]--;
				if (part_connections[part_id][target_part] == 0)
				  {
				    part_messages[part_id]--;
				  }
			      }
			    
			    if (target_part == decision && ep_cnt[start_index + j] > 1)
			      {
				part_connections[part_id][target_part]--;
				if (part_connections[part_id][target_part] == 0)
				  {
				    part_messages[part_id]--;
				  }
			      }
			  }
		      }
		    else
		      {
			int src_net = part_vec[src_nodes[net]];
			int end_cnt = ep_ends[net];
			bool found = false;
			for (int j = 0; j < end_cnt; j++)
			  {
			    if (ep_parts[start_index + j] == decision && ep_cnt[start_index + j] == 1)
			      { //we find the part that we need to decrement a connection
				part_connections[src_net][decision]++;
				
				if (part_connections[src_net][decision] == 1)
				  {
				    part_messages[src_net]++;
				  }
			      }
			    if (ep_parts[start_index + j] == part_id)
			      found = true;
			  }
			if (!found)
			  {
			    part_connections[src_net][part_id]--;
			    if (part_connections[src_net][part_id] == 0)
			      {
				part_messages[src_net]--;
			      }
			  }
		      }
		    
		  }
	      }
	  }
      }
    
  } //pragma
  if (REFI_OPT)
    {
      end = chrono::steady_clock::now();
      diff = end - start;
      cout << "Refinement took : " << chrono::duration<double, milli>(diff).count() / 1000 << " seconds." << endl;
      double ref_time = chrono::duration<double, milli>(diff).count() / 1000;
      total += (chrono::duration<double, milli>(end - start));
      cout << "Total Paritioning took : " << chrono::duration<double, milli>(total).count() / 1000 << " seconds." << endl;
      
      int TV_Refined_out = TV_Refined(ep_starts, ep_cnt, noe, K);
      int CN_Refined_out = CN(ep_ends, noe);
      int TM_Refined_out = TM_scratch(K, part_vec, ep_parts, ep_starts, ep_ends, nov, noe, src_nodes);
      int MGM_Refined_out = MGM_scratch(K, part_vec, ep_parts, ep_starts, ep_ends, nov, noe, src_nodes);
      int MGV_Refined_out = MGV_Scratch(part_vec, row_ptr, col_ind, src_nodes, ep_ends, net_appearence, K, nov, noe, pins);
      int MGAV_Refined_out = MGAV_Scratch(part_vec, row_ptr, col_ind, src_nodes, ep_ends, ep_parts, ep_starts, net_appearence, K, nov, noe, pins);
      
      
      double ref_improved;
      if (REFI_OPT)
        {
	  cout << "After Refinement TV : " << TV_Refined_out << endl;
	  cout << "After Refinement CN : " << CN_Refined_out << endl;
	  cout << "After Refinement TM : " << TM_Refined_out << endl;
	  cout << "After Refinement MGM : " << MGM_Refined_out << endl;
	  cout << "After Refinement MGV : " << MGV_Refined_out << endl;
	  cout << "After Refinement MGAV : " << MGAV_Refined_out << endl;
	  cout << "Improved Percentage TV : " << 100.0 - ((double(TV_Refined_out) / TV_OUT) * 100) << endl;
	  cout << "Improved Percentage CN : " << 100.0 - ((double(CN_Refined_out) / CN_OUT) * 100) << endl;
	  cout << "Improved Percentage TM : " << 100.0 - ((double(TM_Refined_out) / TM_OUT) * 100) << endl;
	  ref_improved = 100.0 - ((double(TM_Refined_out) / TM_OUT) * 100);
	  cout << "Improved Percentage MGM : " << 100.0 - ((double(MGM_Refined_out) / MGM_OUT) * 100) << endl;
	  cout << "Improved Percentage MGV : " << 100.0 - ((double(MGV_Refined_out) / MGV_OUT) * 100) << endl;
	  cout << "Improved Percentage MGAV : " << 100.0 - ((double(MGAV_Refined_out) / MGAV_OUT) * 100) << endl;
        }
      else
        {
          cout << "After Refinement TV : " << -1 << endl;
          cout << "After Refinement CN : " << -1 << endl;
          cout << "After Refinement TM : " << -1 << endl;
          cout << "After Refinement MGM : " << -1 << endl;
          cout << "After Refinement MGV : " << -1 << endl;
          cout << "After Refinement MGAV : " << -1 << endl;
          cout << "Improved Percentage TV : " << -1 << endl;
          cout << "Improved Percentage CN : " << -1 << endl;
          cout << "Improved Percentage TM : " << -1<< endl;
          cout << "Improved Percentage MGM : " << -1 << endl;
          cout << "Improved Percentage MGV : " << -1 << endl;
          cout << "Improved Percentage MGAV : " << -1 << endl;
        }
      cout << "Iterated for : " << REFI_ITER << " times\n";
      cout << "ttt: " << ref_improved << " " << ref_time << " " << REFI_ITER << " K: " << K <<endl;
      cout << "#########################################" << endl;

      
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
      printf("SEQ BEFORE REFINEMENT ACTUAL IMBAL: %f \n", a_imbal);
      printf("#############################\n");
}
}


void partition_gpu(unsigned int *row_ptr, unsigned int *col_ind, unsigned int *net_appearence, int nov, int noe, long long pins, int *edge_wgths, int *src_nodes, int K, int *ordering, string f_name, int REFI_OPT, int REFI_ITER){

  cout << "GPUGPUGPUG" << endl;
  /*------------------------------------------------------------------*/
  /*------------------------------------------------------------------*/
  /*AT THIS POINT, COPYING PARTITIONING FROM CPU, WILL CALL GPU REFINE*/
  /*------------------------------------------------------------------*/
  /*------------------------------------------------------------------*/

  double eps = 1.05;          //ratio of balancing
  double C = (nov / K) * eps; //maximum number of vertices that we can assign to single part
  //int active_part_count;
  int src_part_cnt;
  int src_update_cnt;
  int *part_vec = new int[nov];
  int *ep_starts = new int[noe]; // cumulatively stores degree of nets
  int *ep_ends = new int[noe];   // how many parts have vertice from that net
  int *ep_parts = new int[pins]; // stores the parts that nets are connected
  int *marker = new int[K];
  int *src_marker = new int[K];
  int *src_node_parts = new int[noe];
  int *part_sizes = new int[K];
  int *part_messages = new int[K];
  int *active_src; // stores the src node parts while partitioning a vertex
  int *active_src_net_info;
  int active_src_edges_cnt = 0;
  int *src_update;
  int *mes_saved = new int[K];
  int *ep_cnt = new int[pins];
  int **part_connections = new int *[K];

  for (int i = 0; i < K; i++)
    {
      part_connections[i] = new int[K];
      for (int j = 0; j < K; j++)
        {
	  part_connections[i][j] = 0;
        }
    }
  
  //vertexi al source ise bunun targetlarının partlarının mesaj sayılarını hesapladım ve buna göre en az mesaj gönderene sourceyi koydum
  //eğer target ise elimdeki vertex capacitye göre bak mesaj olan yere koy
  
  //initializations
  auto start = chrono::steady_clock::now();
  initialize(ep_cnt, 0, pins, 0);
  initialize(marker, 0, K, -1);
  initialize(src_marker, 0, K, -1);
  initialize(part_sizes, 0, K, 0);
  initialize(part_messages, 0, K, 0);
  initialize(mes_saved, 0, K, 0);
  initialize(part_vec, 0, nov, -1);
  initialize(src_node_parts, 0, noe, -1);
  initialize(ep_starts, 0, nov, 0);
  initialize(ep_ends, 0, nov, 0);
  initialize(ep_parts, 0, pins, -1);
  //Calculate edge degrees and max vertex degree
  int max_v_deg = 0;
  for (int i = 0; i < nov; i++)
    {
      if (row_ptr[i + 1] - row_ptr[i] > max_v_deg)
        {
	  max_v_deg = row_ptr[i + 1] - row_ptr[i];
        }
    }
  
  for (int i = 0; i < noe; i++)
    {
      ep_starts[i + 1] = net_appearence[i] + ep_starts[i];
    }
  
  active_src = new int[max_v_deg]; // stores the src node parts while partitioning a vertex
  active_src_net_info = new int[max_v_deg];
  src_update = new int[max_v_deg];
  auto end = chrono::steady_clock::now();
  auto diff = end - start;
  cout << "Initialization took: " << chrono::duration<double, milli>(diff).count() / 1000 << " seconds." << endl;
  total += (chrono::duration<double, milli>(end - start));
  //Start Partitioning
  cout << "Partitioning started" << endl;
  start = chrono::steady_clock::now();
  
  for (int z = 0; z < nov; z++)
    {
      int vertex = ordering[z]; // current vertex we want to add to our partition
      int decision;             // it will store the part id of the decision we made for the vertex part
      src_part_cnt = 0;         // number of different parts that the current net set has it's source node in it
      src_update_cnt = 0;       // number of different nets that they have the <vertex> as their source node
      active_src_edges_cnt = 0;
      for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++)
        {
	  int edge = col_ind[i];
	  if (src_nodes[edge] == vertex)
            {
	      src_update[src_update_cnt++] = edge;
            }
	  
	  if (src_node_parts[edge] != -1)
            { // src node is assigned
	      if (src_marker[src_node_parts[edge]] != vertex)
                {
		  src_marker[src_node_parts[edge]] = vertex;
		  active_src[src_part_cnt++] = src_node_parts[edge];
                }
	      active_src_net_info[active_src_edges_cnt++] = edge;
            }
        }
      
      for (int k = 0; k < K; k++)
        {
	  int mes_saved_loc = 0;
	  if (part_sizes[k] < C)
            { // only if the capacity allows
	      
	      for (int i = 0; i < src_part_cnt; i++)
                {
		  int src_part = active_src[i];
		  if (part_connections[src_part][k] > 0)
                    {
		      mes_saved_loc++;
                    }
                }
	      for (int i = 0; i < src_update_cnt; i++)
                {
		  int edge = src_update[i];
		  int end_cnt = ep_ends[edge];
		  int start_index = ep_starts[edge];
		  for (int j = 0; j < end_cnt; j++)
                    {
		      int part_id = ep_parts[start_index + j];
		      if (part_connections[k][part_id] > 0)
                        {
			  mes_saved_loc++;
                        }
                    }
                }
            }
	  mes_saved[k] = mes_saved_loc;
        }
      
      int max_mes_saved = -1; // maximum number of messages saved for a potential decision on a part <k>
      int max_index = -1;     // part id of the maximum message saved part
      
      for (int k = 0; k < K; k++)
        {
	  if (mes_saved[k] > max_mes_saved || (mes_saved[k] == max_mes_saved &&
					       (part_sizes[k] < part_sizes[max_index])))
            {
	      max_index = k;
	      max_mes_saved = mes_saved[k];
            }
        }
      
      if (max_index == -1)
        {
	  int minLoad = INT_MAX;
	  int minLoadIndex = -1;
	  for (int i = 0; i < K; i++)
            {
	      if (part_sizes[i] < minLoad)
                {
		  minLoad = part_sizes[i];
		  minLoadIndex = i;
                }
            }
	  max_index = minLoadIndex;
        }
      decision = max_index;
      part_vec[vertex] = decision; //vertex part is decided
      part_sizes[decision]++;
      
      //Update the edge part info
      for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++)
        {
	  int edge = col_ind[i];
	  int end_cnt = ep_ends[edge];
	  int start_index = ep_starts[edge];
	  bool found = false;
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
            }
        }
      //update messages and source parts
      for (int i = 0; i < active_src_edges_cnt; i++)
        {
	  int edge = active_src_net_info[i];
	  int src_part = src_node_parts[edge];
	  
	  if (src_part == decision)
	    continue;
	  if (part_connections[src_part][decision] == 0)
            {
	      part_messages[src_part]++;
	      part_connections[src_part][decision] = 1;
            }
	  else
            {
	      int end_cnt = ep_ends[edge];
	      int start_index = ep_starts[edge];
	      for (int j = 0; j < end_cnt; j++)
                {
		  int part_id = ep_parts[start_index + j];
		  if (part_id == decision)
                    {
		      if (ep_cnt[start_index + j] == 1) //if net caused a new connection
                        {
			  part_connections[src_part][decision]++;
                        }
		      break;
                    }
                }
            }
        }
      
      for (int i = 0; i < src_update_cnt; i++)
        {
	  int edge = src_update[i];
	  src_node_parts[edge] = decision;
	  int end_cnt = ep_ends[edge];
	  int start_index = ep_starts[edge];
	  
	  for (int j = 0; j < end_cnt; j++)
            {
	      int part_id = ep_parts[start_index + j];
	      if (part_id == decision)
		continue;
	      if (part_connections[decision][part_id] == 0)
                {
		  part_messages[decision]++;
		  part_connections[decision][part_id] = 1;
                }
	      else
                {
		  part_connections[decision][part_id]++;
                }
            }
        }
    }
  
  end = chrono::steady_clock::now();
  diff = end - start;
  total += (chrono::duration<double, milli>(end - start));
  cout << "Paritioning took : " << chrono::duration<double, milli>(total).count() / 1000 << " seconds." << endl;
  
  int TV_OUT = TV(ep_ends, noe);
  int CN_OUT = CN(ep_ends, noe);
  int TM_OUT = TM_scratch(K, part_vec, ep_parts, ep_starts, ep_ends, nov, noe, src_nodes);
  int MGM_OUT = MGM_scratch(K, part_vec, ep_parts, ep_starts, ep_ends, nov, noe, src_nodes);
  int MGV_OUT = MGV_Scratch(part_vec, row_ptr, col_ind, src_nodes, ep_ends, net_appearence, K, nov, noe, pins);
  int MGAV_OUT = MGAV_Scratch(part_vec, row_ptr, col_ind, src_nodes, ep_ends, ep_parts, ep_starts, net_appearence, K, nov, noe, pins);
  
  cout << "TV : " << TV_OUT << endl;
  cout << "CN : " << CN_OUT << endl;
  cout << "TM : " << TM_OUT << endl;
  cout << "MGM : " << MGM_OUT << endl;
  cout << "MGV : " << MGV_OUT << endl;
  cout << "MGAV : " << MGAV_OUT << endl;

  /**/
  for(int i = 0; i < pins; i++){
    if(ep_parts[i] < 0)
      printf("HOST ep_parts[i], %d, val: %d\n", i,ep_parts[i]);
  }
  /**/
  
  /*------------------------------------------------------------------*/
  /*------------------------------------------------------------------*/
  /*AT THIS POINT, COPYING PARTITIONING FROM CPU, WILL CALL GPU REFINE*/
  /*------------------------------------------------------------------*/
  /*------------------------------------------------------------------*/

  auto gstart = chrono::steady_clock::now();
  partition_gpu_wrapper(row_ptr, col_ind, net_appearence, nov, noe, pins, edge_wgths, src_nodes, K, ordering, REFI_OPT, REFI_ITER, part_vec, ep_starts, ep_parts, ep_ends, ep_cnt, part_connections, part_sizes, part_messages);//GPU
  auto partial = chrono::steady_clock::now();
  auto partial_diff = end - start;
  auto halfgain = (chrono::duration<double, milli>(partial - gstart));
  cout << "GPU Partial took : " << chrono::duration<double, milli>(halfgain).count() / 1000 << " seconds." << endl;

  int TV_Refined_out = TV_Refined(ep_starts, ep_cnt, noe, K);
  int CN_Refined_out = CN(ep_ends, noe);
  int TM_Refined_out = TM_scratch(K, part_vec, ep_parts, ep_starts, ep_ends, nov, noe, src_nodes);
  int MGM_Refined_out = MGM_scratch(K, part_vec, ep_parts, ep_starts, ep_ends, nov, noe, src_nodes);
  int MGV_Refined_out = MGV_Scratch(part_vec, row_ptr, col_ind, src_nodes, ep_ends, net_appearence, K, nov, noe, pins);
  int MGAV_Refined_out = MGAV_Scratch(part_vec, row_ptr, col_ind, src_nodes, ep_ends, ep_parts, ep_starts, net_appearence, K, nov, noe, pins);
  
  
  
  
  cout << "After Refinement TV : " << TV_Refined_out << endl;
  cout << "After Refinement CN : " << CN_Refined_out << endl;
  cout << "After Refinement TM : " << TM_Refined_out << endl;
  cout << "After Refinement MGM : " << MGM_Refined_out << endl;
  cout << "After Refinement MGV : " << MGV_Refined_out << endl;
  cout << "After Refinement MGAV : " << MGAV_Refined_out << endl;
  cout << "Improved Percentage TV : " << 100.0 - ((double(TV_Refined_out) / TV_OUT) * 100) << endl;
  cout << "Improved Percentage CN : " << 100.0 - ((double(CN_Refined_out) / CN_OUT) * 100) << endl;
  cout << "Improved Percentage TM : " << 100.0 - ((double(TM_Refined_out) / TM_OUT) * 100) << endl;
  cout << "Improved Percentage MGM : " << 100.0 - ((double(MGM_Refined_out) / MGM_OUT) * 100) << endl;
  cout << "Improved Percentage MGV : " << 100.0 - ((double(MGV_Refined_out) / MGV_OUT) * 100) << endl;
  cout << "Improved Percentage MGAV : " << 100.0 - ((double(MGAV_Refined_out) / MGAV_OUT) * 100) << endl;

}

void partition_gpu_optimized(unsigned int *row_ptr, unsigned int *col_ind, unsigned int *net_appearence, int nov, int noe, long long pins, int *edge_wgths, int *src_nodes, int K, int *ordering, string f_name, int REFI_OPT, int REFI_ITER){
  
  /*------------------------------------------------------------------*/
  /*------------------------------------------------------------------*/
  /*AT THIS POINT, COPYING PARTITIONING FROM CPU, WILL CALL GPU REFINE*/
  /*------------------------------------------------------------------*/
  /*------------------------------------------------------------------*/

  double eps = 1.05;          //ratio of balancing
  double C = (nov / K) * eps; //maximum number of vertices that we can assign to single part
  //int active_part_count;
  int src_part_cnt;
  int src_update_cnt;
  int *part_vec = new int[nov];
  int *ep_starts = new int[noe]; // cumulatively stores degree of nets
  int *ep_ends = new int[noe];   // how many parts have vertice from that net
  int *ep_parts = new int[pins]; // stores the parts that nets are connected
  int *marker = new int[K];
  int *src_marker = new int[K];
  int *src_node_parts = new int[noe];
  int *part_sizes = new int[K];
  int *part_messages = new int[K];
  int *active_src; // stores the src node parts while partitioning a vertex
  int *active_src_net_info;
  int active_src_edges_cnt = 0;
  int *src_update;
  int *mes_saved = new int[K];
  int *ep_cnt = new int[pins];
  int **part_connections = new int *[K];

  for (int i = 0; i < K; i++)
    {
      part_connections[i] = new int[K];
      for (int j = 0; j < K; j++)
        {
	  part_connections[i][j] = 0;
        }
    }
  
  //vertexi al source ise bunun targetlarının partlarının mesaj sayılarını hesapladım ve buna göre en az mesaj gönderene sourceyi koydum
  //eğer target ise elimdeki vertex capacitye göre bak mesaj olan yere koy
  
  //initializations
  auto start = chrono::steady_clock::now();
  initialize(ep_cnt, 0, pins, 0);
  initialize(marker, 0, K, -1);
  initialize(src_marker, 0, K, -1);
  initialize(part_sizes, 0, K, 0);
  initialize(part_messages, 0, K, 0);
  initialize(mes_saved, 0, K, 0);
  initialize(part_vec, 0, nov, -1);
  initialize(src_node_parts, 0, noe, -1);
  initialize(ep_starts, 0, nov, 0);
  initialize(ep_ends, 0, nov, 0);
  initialize(ep_parts, 0, pins, -1);
  //Calculate edge degrees and max vertex degree
  int max_v_deg = 0;
  for (int i = 0; i < nov; i++)
    {
      if (row_ptr[i + 1] - row_ptr[i] > max_v_deg)
        {
	  max_v_deg = row_ptr[i + 1] - row_ptr[i];
        }
    }
  
  for (int i = 0; i < noe; i++)
    {
      ep_starts[i + 1] = net_appearence[i] + ep_starts[i];
    }
  
  active_src = new int[max_v_deg]; // stores the src node parts while partitioning a vertex
  active_src_net_info = new int[max_v_deg];
  src_update = new int[max_v_deg];
  auto end = chrono::steady_clock::now();
  auto diff = end - start;
  cout << "Initialization took: " << chrono::duration<double, milli>(diff).count() / 1000 << " seconds." << endl;
  total += (chrono::duration<double, milli>(end - start));
  //Start Partitioning
  cout << "Partitioning started" << endl;
  start = chrono::steady_clock::now();
  
  for (int z = 0; z < nov; z++)
    {
      int vertex = ordering[z]; // current vertex we want to add to our partition
      int decision;             // it will store the part id of the decision we made for the vertex part
      src_part_cnt = 0;         // number of different parts that the current net set has it's source node in it
      src_update_cnt = 0;       // number of different nets that they have the <vertex> as their source node
      active_src_edges_cnt = 0;
      for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++)
        {
	  int edge = col_ind[i];
	  if (src_nodes[edge] == vertex)
            {
	      src_update[src_update_cnt++] = edge;
            }
	  
	  if (src_node_parts[edge] != -1)
            { // src node is assigned
	      if (src_marker[src_node_parts[edge]] != vertex)
                {
		  src_marker[src_node_parts[edge]] = vertex;
		  active_src[src_part_cnt++] = src_node_parts[edge];
                }
	      active_src_net_info[active_src_edges_cnt++] = edge;
            }
        }
      
      for (int k = 0; k < K; k++)
        {
	  int mes_saved_loc = 0;
	  if (part_sizes[k] < C)
            { // only if the capacity allows
	      
	      for (int i = 0; i < src_part_cnt; i++)
                {
		  int src_part = active_src[i];
		  if (part_connections[src_part][k] > 0)
                    {
		      mes_saved_loc++;
                    }
                }
	      for (int i = 0; i < src_update_cnt; i++)
                {
		  int edge = src_update[i];
		  int end_cnt = ep_ends[edge];
		  int start_index = ep_starts[edge];
		  for (int j = 0; j < end_cnt; j++)
                    {
		      int part_id = ep_parts[start_index + j];
		      if (part_connections[k][part_id] > 0)
                        {
			  mes_saved_loc++;
                        }
                    }
                }
            }
	  mes_saved[k] = mes_saved_loc;
        }
      
      int max_mes_saved = -1; // maximum number of messages saved for a potential decision on a part <k>
      int max_index = -1;     // part id of the maximum message saved part
      
      for (int k = 0; k < K; k++)
        {
	  if (mes_saved[k] > max_mes_saved || (mes_saved[k] == max_mes_saved &&
					       (part_sizes[k] < part_sizes[max_index])))
            {
	      max_index = k;
	      max_mes_saved = mes_saved[k];
            }
        }
      
      if (max_index == -1)
        {
	  int minLoad = INT_MAX;
	  int minLoadIndex = -1;
	  for (int i = 0; i < K; i++)
            {
	      if (part_sizes[i] < minLoad)
                {
		  minLoad = part_sizes[i];
		  minLoadIndex = i;
                }
            }
	  max_index = minLoadIndex;
        }
      decision = max_index;
      part_vec[vertex] = decision; //vertex part is decided
      part_sizes[decision]++;
      
      //Update the edge part info
      for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++)
        {
	  int edge = col_ind[i];
	  int end_cnt = ep_ends[edge];
	  int start_index = ep_starts[edge];
	  bool found = false;
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
            }
        }
      //update messages and source parts
      for (int i = 0; i < active_src_edges_cnt; i++)
        {
	  int edge = active_src_net_info[i];
	  int src_part = src_node_parts[edge];
	  
	  if (src_part == decision)
	    continue;
	  if (part_connections[src_part][decision] == 0)
            {
	      part_messages[src_part]++;
	      part_connections[src_part][decision] = 1;
            }
	  else
            {
	      int end_cnt = ep_ends[edge];
	      int start_index = ep_starts[edge];
	      for (int j = 0; j < end_cnt; j++)
                {
		  int part_id = ep_parts[start_index + j];
		  if (part_id == decision)
                    {
		      if (ep_cnt[start_index + j] == 1) //if net caused a new connection
                        {
			  part_connections[src_part][decision]++;
                        }
		      break;
                    }
                }
            }
        }
      
      for (int i = 0; i < src_update_cnt; i++)
        {
	  int edge = src_update[i];
	  src_node_parts[edge] = decision;
	  int end_cnt = ep_ends[edge];
	  int start_index = ep_starts[edge];
	  
	  for (int j = 0; j < end_cnt; j++)
            {
	      int part_id = ep_parts[start_index + j];
	      if (part_id == decision)
		continue;
	      if (part_connections[decision][part_id] == 0)
                {
		  part_messages[decision]++;
		  part_connections[decision][part_id] = 1;
                }
	      else
                {
		  part_connections[decision][part_id]++;
                }
            }
        }
    }
  
  end = chrono::steady_clock::now();
  diff = end - start;
  total += (chrono::duration<double, milli>(end - start));
  cout << "Paritioning took : " << chrono::duration<double, milli>(total).count() / 1000 << " seconds." << endl;

  int TV_OUT = TV(ep_ends, noe);
  int CN_OUT = CN(ep_ends, noe);
  int TM_OUT = TM_scratch(K, part_vec, ep_parts, ep_starts, ep_ends, nov, noe, src_nodes);
  int MGM_OUT = MGM_scratch(K, part_vec, ep_parts, ep_starts, ep_ends, nov, noe, src_nodes);
  int MGV_OUT = MGV_Scratch(part_vec, row_ptr, col_ind, src_nodes, ep_ends, net_appearence, K, nov, noe, pins);
  int MGAV_OUT = MGAV_Scratch(part_vec, row_ptr, col_ind, src_nodes, ep_ends, ep_parts, ep_starts, net_appearence, K, nov, noe, pins);
  
  cout << "TV : " << TV_OUT << endl;
  cout << "CN : " << CN_OUT << endl;
  cout << "TM : " << TM_OUT << endl;
  cout << "MGM : " << MGM_OUT << endl;
  cout << "MGV : " << MGV_OUT << endl;
  cout << "MGAV : " << MGAV_OUT << endl;

  /*------------------------------------------------------------------*/
  /*------------------------------------------------------------------*/
  /*AT THIS POINT, COPYING PARTITIONING FROM CPU, WILL CALL GPU REFINE*/
  /*------------------------------------------------------------------*/
  /*------------------------------------------------------------------*/

  auto gstart = chrono::steady_clock::now();
  cout << "Am I able to call the kernel" << endl;
  partition_gpu_optimized_wrapper(row_ptr, col_ind, net_appearence, nov, noe, pins, edge_wgths, src_nodes, K, ordering, REFI_OPT, REFI_ITER, part_vec, ep_starts, ep_parts, ep_ends, ep_cnt, part_connections, part_sizes, part_messages);//GPU
  cout << "Am I able to return?" << endl;
  auto partial = chrono::steady_clock::now();
  auto partial_diff = end - start;
  auto halfgain = (chrono::duration<double, milli>(partial - gstart));
  cout << "GPU Partial took : " << chrono::duration<double, milli>(halfgain).count() / 1000 << " seconds." << endl;

  int TV_Refined_out = TV_Refined(ep_starts, ep_cnt, noe, K);
  int CN_Refined_out = CN(ep_ends, noe);
  int TM_Refined_out = TM_scratch(K, part_vec, ep_parts, ep_starts, ep_ends, nov, noe, src_nodes);
  int MGM_Refined_out = MGM_scratch(K, part_vec, ep_parts, ep_starts, ep_ends, nov, noe, src_nodes);
  int MGV_Refined_out = MGV_Scratch(part_vec, row_ptr, col_ind, src_nodes, ep_ends, net_appearence, K, nov, noe, pins);
  int MGAV_Refined_out = MGAV_Scratch(part_vec, row_ptr, col_ind, src_nodes, ep_ends, ep_parts, ep_starts, net_appearence, K, nov, noe, pins);
  
  
  
  
  cout << "After Refinement TV : " << TV_Refined_out << endl;
  cout << "After Refinement CN : " << CN_Refined_out << endl;
  cout << "After Refinement TM : " << TM_Refined_out << endl;
  cout << "After Refinement MGM : " << MGM_Refined_out << endl;
  cout << "After Refinement MGV : " << MGV_Refined_out << endl;
  cout << "After Refinement MGAV : " << MGAV_Refined_out << endl;
  cout << "Improved Percentage TV : " << 100.0 - ((double(TV_Refined_out) / TV_OUT) * 100) << endl;
  cout << "Improved Percentage CN : " << 100.0 - ((double(CN_Refined_out) / CN_OUT) * 100) << endl;
  cout << "Improved Percentage TM : " << 100.0 - ((double(TM_Refined_out) / TM_OUT) * 100) << endl;
  cout << "Improved Percentage MGM : " << 100.0 - ((double(MGM_Refined_out) / MGM_OUT) * 100) << endl;
  cout << "Improved Percentage MGV : " << 100.0 - ((double(MGV_Refined_out) / MGV_OUT) * 100) << endl;
  cout << "Improved Percentage MGAV : " << 100.0 - ((double(MGAV_Refined_out) / MGAV_OUT) * 100) << endl; 
}


void partition_gpu_threads(unsigned int *row_ptr, unsigned int *col_ind, unsigned int *net_appearence, int nov, int noe, long long pins, int *edge_wgths, int *src_nodes, int K, int *ordering, string f_name, int REFI_OPT, int REFI_ITER){
  
  /*------------------------------------------------------------------*/
  /*------------------------------------------------------------------*/
  /*AT THIS POINT, COPYING PARTITIONING FROM CPU, WILL CALL GPU REFINE*/
  /*------------------------------------------------------------------*/
  /*------------------------------------------------------------------*/

  double eps = 1.05;          //ratio of balancing
  double C = (nov / K) * eps; //maximum number of vertices that we can assign to single part
  //int active_part_count;
  int src_part_cnt;
  int src_update_cnt;
  int *part_vec = new int[nov];
  int *ep_starts = new int[noe]; // cumulatively stores degree of nets
  int *ep_ends = new int[noe];   // how many parts have vertice from that net
  int *ep_parts = new int[pins]; // stores the parts that nets are connected
  int *marker = new int[K];
  int *src_marker = new int[K];
  int *src_node_parts = new int[noe];
  int *part_sizes = new int[K];
  int *part_messages = new int[K];
  int *active_src; // stores the src node parts while partitioning a vertex
  int *active_src_net_info;
  int active_src_edges_cnt = 0;
  int *src_update;
  int *mes_saved = new int[K];
  int *ep_cnt = new int[pins];
  int **part_connections = new int *[K];

  for (int i = 0; i < K; i++)
    {
      part_connections[i] = new int[K];
      for (int j = 0; j < K; j++)
        {
	  part_connections[i][j] = 0;
        }
    }
  
  //vertexi al source ise bunun targetlarının partlarının mesaj sayılarını hesapladım ve buna göre en az mesaj gönderene sourceyi koydum
  //eğer target ise elimdeki vertex capacitye göre bak mesaj olan yere koy
  
  //initializations
  auto start = chrono::steady_clock::now();
  initialize(ep_cnt, 0, pins, 0);
  initialize(marker, 0, K, -1);
  initialize(src_marker, 0, K, -1);
  initialize(part_sizes, 0, K, 0);
  initialize(part_messages, 0, K, 0);
  initialize(mes_saved, 0, K, 0);
  initialize(part_vec, 0, nov, -1);
  initialize(src_node_parts, 0, noe, -1);
  initialize(ep_starts, 0, nov, 0);
  initialize(ep_ends, 0, nov, 0);
  initialize(ep_parts, 0, pins, -1);
  //Calculate edge degrees and max vertex degree
  int max_v_deg = 0;
  for (int i = 0; i < nov; i++)
    {
      if (row_ptr[i + 1] - row_ptr[i] > max_v_deg)
        {
	  max_v_deg = row_ptr[i + 1] - row_ptr[i];
        }
    }
  
  for (int i = 0; i < noe; i++)
    {
      ep_starts[i + 1] = net_appearence[i] + ep_starts[i];
    }
  
  active_src = new int[max_v_deg]; // stores the src node parts while partitioning a vertex
  active_src_net_info = new int[max_v_deg];
  src_update = new int[max_v_deg];
  auto end = chrono::steady_clock::now();
  auto diff = end - start;
  cout << "Initialization took: " << chrono::duration<double, milli>(diff).count() / 1000 << " seconds." << endl;
  total += (chrono::duration<double, milli>(end - start));
  //Start Partitioning
  cout << "Partitioning started" << endl;
  start = chrono::steady_clock::now();
  
  for (int z = 0; z < nov; z++)
    {
      int vertex = ordering[z]; // current vertex we want to add to our partition
      int decision;             // it will store the part id of the decision we made for the vertex part
      src_part_cnt = 0;         // number of different parts that the current net set has it's source node in it
      src_update_cnt = 0;       // number of different nets that they have the <vertex> as their source node
      active_src_edges_cnt = 0;
      for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++)
        {
	  int edge = col_ind[i];
	  if (src_nodes[edge] == vertex)
            {
	      src_update[src_update_cnt++] = edge;
            }
	  
	  if (src_node_parts[edge] != -1)
            { // src node is assigned
	      if (src_marker[src_node_parts[edge]] != vertex)
                {
		  src_marker[src_node_parts[edge]] = vertex;
		  active_src[src_part_cnt++] = src_node_parts[edge];
                }
	      active_src_net_info[active_src_edges_cnt++] = edge;
            }
        }
      
      for (int k = 0; k < K; k++)
        {
	  int mes_saved_loc = 0;
	  if (part_sizes[k] < C)
            { // only if the capacity allows
	      
	      for (int i = 0; i < src_part_cnt; i++)
                {
		  int src_part = active_src[i];
		  if (part_connections[src_part][k] > 0)
                    {
		      mes_saved_loc++;
                    }
                }
	      for (int i = 0; i < src_update_cnt; i++)
                {
		  int edge = src_update[i];
		  int end_cnt = ep_ends[edge];
		  int start_index = ep_starts[edge];
		  for (int j = 0; j < end_cnt; j++)
                    {
		      int part_id = ep_parts[start_index + j];
		      if (part_connections[k][part_id] > 0)
                        {
			  mes_saved_loc++;
                        }
                    }
                }
            }
	  mes_saved[k] = mes_saved_loc;
        }
      
      int max_mes_saved = -1; // maximum number of messages saved for a potential decision on a part <k>
      int max_index = -1;     // part id of the maximum message saved part
      
      for (int k = 0; k < K; k++)
        {
	  if (mes_saved[k] > max_mes_saved || (mes_saved[k] == max_mes_saved &&
					       (part_sizes[k] < part_sizes[max_index])))
            {
	      max_index = k;
	      max_mes_saved = mes_saved[k];
            }
        }
      
      if (max_index == -1)
        {
	  int minLoad = INT_MAX;
	  int minLoadIndex = -1;
	  for (int i = 0; i < K; i++)
            {
	      if (part_sizes[i] < minLoad)
                {
		  minLoad = part_sizes[i];
		  minLoadIndex = i;
                }
            }
	  max_index = minLoadIndex;
        }
      decision = max_index;
      part_vec[vertex] = decision; //vertex part is decided
      part_sizes[decision]++;
      
      //Update the edge part info
      for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++)
        {
	  int edge = col_ind[i];
	  int end_cnt = ep_ends[edge];
	  int start_index = ep_starts[edge];
	  bool found = false;
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
            }
        }
      //update messages and source parts
      for (int i = 0; i < active_src_edges_cnt; i++)
        {
	  int edge = active_src_net_info[i];
	  int src_part = src_node_parts[edge];
	  
	  if (src_part == decision)
	    continue;
	  if (part_connections[src_part][decision] == 0)
            {
	      part_messages[src_part]++;
	      part_connections[src_part][decision] = 1;
            }
	  else
            {
	      int end_cnt = ep_ends[edge];
	      int start_index = ep_starts[edge];
	      for (int j = 0; j < end_cnt; j++)
                {
		  int part_id = ep_parts[start_index + j];
		  if (part_id == decision)
                    {
		      if (ep_cnt[start_index + j] == 1) //if net caused a new connection
                        {
			  part_connections[src_part][decision]++;
                        }
		      break;
                    }
                }
            }
        }
      
      for (int i = 0; i < src_update_cnt; i++)
        {
	  int edge = src_update[i];
	  src_node_parts[edge] = decision;
	  int end_cnt = ep_ends[edge];
	  int start_index = ep_starts[edge];
	  
	  for (int j = 0; j < end_cnt; j++)
            {
	      int part_id = ep_parts[start_index + j];
	      if (part_id == decision)
		continue;
	      if (part_connections[decision][part_id] == 0)
                {
		  part_messages[decision]++;
		  part_connections[decision][part_id] = 1;
                }
	      else
                {
		  part_connections[decision][part_id]++;
                }
            }
        }
    }
  
  end = chrono::steady_clock::now();
  diff = end - start;
  total += (chrono::duration<double, milli>(end - start));
  cout << "Paritioning took : " << chrono::duration<double, milli>(total).count() / 1000 << " seconds." << endl;

  int TV_OUT = TV(ep_ends, noe);
  int CN_OUT = CN(ep_ends, noe);
  int TM_OUT = TM_scratch(K, part_vec, ep_parts, ep_starts, ep_ends, nov, noe, src_nodes);
  int MGM_OUT = MGM_scratch(K, part_vec, ep_parts, ep_starts, ep_ends, nov, noe, src_nodes);
  int MGV_OUT = MGV_Scratch(part_vec, row_ptr, col_ind, src_nodes, ep_ends, net_appearence, K, nov, noe, pins);
  int MGAV_OUT = MGAV_Scratch(part_vec, row_ptr, col_ind, src_nodes, ep_ends, ep_parts, ep_starts, net_appearence, K, nov, noe, pins);
  
  cout << "TV : " << TV_OUT << endl;
  cout << "CN : " << CN_OUT << endl;
  cout << "TM : " << TM_OUT << endl;
  cout << "MGM : " << MGM_OUT << endl;
  cout << "MGV : " << MGV_OUT << endl;
  cout << "MGAV : " << MGAV_OUT << endl;

  /*------------------------------------------------------------------*/
  /*------------------------------------------------------------------*/
  /*AT THIS POINT, COPYING PARTITIONING FROM CPU, WILL CALL GPU REFINE*/
  /*------------------------------------------------------------------*/
  /*------------------------------------------------------------------*/

  auto gstart = chrono::steady_clock::now();
  partition_gpu_threads_wrapper(row_ptr, col_ind, net_appearence, nov, noe, pins, edge_wgths, src_nodes, K, ordering, REFI_OPT, REFI_ITER, part_vec, ep_starts, ep_parts, ep_ends, ep_cnt, part_connections, part_sizes, part_messages);//GPU
  auto partial = chrono::steady_clock::now();
  auto partial_diff = end - start;
  auto halfgain = (chrono::duration<double, milli>(partial - gstart));
  double gpu_ref = chrono::duration<double, milli>(halfgain).count() / 1000;
  cout << "GPU Partial took : " << chrono::duration<double, milli>(halfgain).count() / 1000 << " seconds." << endl;

  int TV_Refined_out = TV_Refined(ep_starts, ep_cnt, noe, K);
  int CN_Refined_out = CN(ep_ends, noe);
  int TM_Refined_out = TM_scratch(K, part_vec, ep_parts, ep_starts, ep_ends, nov, noe, src_nodes);
  int MGM_Refined_out = MGM_scratch(K, part_vec, ep_parts, ep_starts, ep_ends, nov, noe, src_nodes);
  int MGV_Refined_out = MGV_Scratch(part_vec, row_ptr, col_ind, src_nodes, ep_ends, net_appearence, K, nov, noe, pins);
  int MGAV_Refined_out = MGAV_Scratch(part_vec, row_ptr, col_ind, src_nodes, ep_ends, ep_parts, ep_starts, net_appearence, K, nov, noe, pins);

  double gpu_improved = 100.0 - ((double(TM_Refined_out) / TM_OUT) * 100);
  
  cout << "After Refinement TV : " << TV_Refined_out << endl;
  cout << "After Refinement CN : " << CN_Refined_out << endl;
  cout << "After Refinement TM : " << TM_Refined_out << endl;
  cout << "After Refinement MGM : " << MGM_Refined_out << endl;
  cout << "After Refinement MGV : " << MGV_Refined_out << endl;
  cout << "After Refinement MGAV : " << MGAV_Refined_out << endl;
  cout << "Improved Percentage TV : " << 100.0 - ((double(TV_Refined_out) / TV_OUT) * 100) << endl;
  cout << "Improved Percentage CN : " << 100.0 - ((double(CN_Refined_out) / CN_OUT) * 100) << endl;
  cout << "Improved Percentage TM : " << 100.0 - ((double(TM_Refined_out) / TM_OUT) * 100) << endl;
  cout << "Improved Percentage MGM : " << 100.0 - ((double(MGM_Refined_out) / MGM_OUT) * 100) << endl;
  cout << "Improved Percentage MGV : " << 100.0 - ((double(MGV_Refined_out) / MGV_OUT) * 100) << endl;
  cout << "Improved Percentage MGAV : " << 100.0 - ((double(MGAV_Refined_out) / MGAV_OUT) * 100) << endl;
  cout << "ttt: " << gpu_improved << " " << gpu_ref << " " << REFI_ITER << " K: " << K <<endl;
  cout << "#########################################" << endl;
}
