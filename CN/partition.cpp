#include "partition.h"
#include <omp.h>
#include <signal.h>

void segfault() { raise(SIGSEGV); }

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
				  //string f_name,
				  int REFI_OPT,
				  int REFI_ITER,
				  unsigned int *row_ptr_inv,
				  unsigned int *col_ind_inv,
				  int* part_sizes,
				  int* ep_starts,
				  int* ep_parts,
				  int* ep_ends,
				  int* ep_cnt,
				  int* part_vec);

void LDG(bool &check_selected_vertex, int &K, int &active_part_count, int *part_sizes, int &decision, int *part_vec, int *active_parts, double &C, double *scores, int &vertex)
{
    //LDG STARTS
    double max_val = 0;
    int max_index = -1;

    //for exluding selected vertex
    if (check_selected_vertex)
    {
        part_sizes[part_vec[vertex]]--;
    }

    for (int i = 0; i < active_part_count; i++)
    {
        int part_id = active_parts[i];
        double penalty = 1 - (1.0 * part_sizes[part_id] / C);
        double score = scores[part_id] * penalty;
        if (score > 0)
        {
            if ((max_index == -1) || (score > max_val) || ((score == max_val) && (part_sizes[part_id] < part_sizes[max_index])))
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

    if (check_selected_vertex)
    {
        part_sizes[part_vec[vertex]]++;
    }
    //LDG ENDS
}

void partition(unsigned int *row_ptr, unsigned int *col_ind, unsigned int *net_appearence, int nov, int noe, long long pins, int *edge_wgths, int *src_nodes, int K,
               int *ordering, string &f_name, int &REFI_OPT, int &REFI_ITER, unsigned int *row_ptr_inv, unsigned int *col_ind_inv)
{
    int active_part_count;          // total number of active partitions while partitioning of a vertex
    int *part_vec = new int[nov];   // stores the part id's of the vertices
    int *marker = new int[K];       // stores the activeness information to determine active parts efficiently
    int *active_parts = new int[K]; // stores the active parts
    int *ep_starts = new int[nov];  // degree of network
    int *ep_ends = new int[nov];    // how many parts that network have been split
    int *ep_parts = new int[pins];  // stores the edge-part info
    int *part_sizes = new int[K];   // stores number of verties in each part
    double eps = 1.05;
    double C = (nov / K) * eps;
    double *scores = new double[K];
    int *ep_cnt = new int[pins];

    int partCount = K; 

    //initializations
    auto start = chrono::steady_clock::now();
    initialize(ep_cnt, 0, pins, 0);
    initialize(marker, 0, K, -1);
    initialize(active_parts, 0, K, -1);
    initialize(part_sizes, 0, K, 0);
    initialize(scores, 0, K, 0.0);
    initialize(part_vec, 0, nov, -1);
    initialize(ep_starts, 0, nov, 0); //checked
    initialize(ep_ends, 0, nov, 0);
    initialize(ep_parts, 0, pins, -1);

    //Calculate edge degrees
    for (int i = 0; i < nov; i++)
    {
        for (int k = row_ptr[i]; k < row_ptr[i + 1]; k++)
        {
            int edge = col_ind[k];
            ep_starts[edge + 1]++;
        }
    }
    //Store edge degrees in a cumulative way
    for (int i = 2; i <= nov; i++)
    {
        ep_starts[i] = ep_starts[i - 1] + ep_starts[i];
    }

    auto end = chrono::steady_clock::now();
    auto diff = end - start;
    cout << "Initialization took : " << chrono::duration<double, milli>(diff).count() / 1000 << " seconds." << endl;

    //Start Partitioning
    cout << "Partitioning started " << endl;
    start = chrono::steady_clock::now();
    for (int z = 0; z < nov; z++) //for each vertex
    {
        int vertex = ordering[z];
        int decision; // Stores the part decision for the <vertex>
        active_part_count = 0;

        for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++)
        {
            int edge = col_ind[i];
            if (ep_ends[edge] == 1) //information is valuable if only net is not cut
            {
                int part_id = ep_parts[ep_starts[edge]];
                if (marker[part_id] != vertex)
                {
                    scores[part_id] = 1;
                    marker[part_id] = vertex;
                    active_parts[active_part_count] = part_id;
                    active_part_count++;
                }
                else
                {
                    scores[part_id]++;
                }
            }
        }

        //Apply LDG Approach
        double max_val = 0;
        int max_index = -1;
        for (int i = 0; i < active_part_count; i++)
        {
            int part_id = active_parts[i];
            double penalty = 1 - (1.0 * part_sizes[part_id] / C);
            double score = scores[part_id] * penalty;
            if (score > 0)
            {
                if ((max_index == -1) || (score > max_val) || ((score == max_val) && (part_sizes[part_id] < part_sizes[max_index])))
                {
                    max_val = score;
                    max_index = part_id;
                }
            }
        }

        if (max_index == -1) //if every net is split or no net has been observed
        {
            //find min loaded
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
        part_sizes[decision]++;
        part_vec[vertex] = decision; //vertex part is decided

        //Update the edge part info
        for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++) // Traversing the edge list of the vertex
        {
            int edge = col_ind[i];
            int j;
            bool found = false;
            int end_cnt = ep_ends[edge];
            int start_index = ep_starts[edge];
            for (j = 0; j < end_cnt; j++)
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
    }

    end = chrono::steady_clock::now();
    diff = end - start;
    total += (chrono::duration<double, milli>(end - start));
    cout << "Paritioning took : " << chrono::duration<double, milli>(diff).count() / 1000 << " seconds." << endl;

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

    // REFINEMENT STARTS :))))))))
    start = chrono::steady_clock::now();

    for (int iter = 0; REFI_OPT && iter < REFI_ITER; iter++)
    {
        int *v_list_push = new int[nov];
        int v_list_cntr;
        int rem_cnt = 0;
        int inc_cnt = 0;

        for (int edge = 0; edge < noe; edge++)
        {
            if (ep_ends[edge] == 2)
            {
                v_list_cntr = 0;
                int part1 = -1, part2 = -1;
                int start_index = ep_starts[edge];
                int end_cnt = ep_ends[edge];
                // to find p1, p2
                int part1_count, part2_count;
                for (int j = 0; j < end_cnt; j++)

                {
                    if (part1 != -1)
                    {
                        part2 = ep_parts[start_index + j];
                        part2_count = ep_cnt[start_index + j];
                        break;
                    }
                    else
                    {
                        part1 = ep_parts[start_index + j];
                        part1_count = ep_cnt[start_index + j];
                    }
                }

                // now we know p1 and p2
                int old_decision, new_decision;
                if (part1_count < part2_count)
                {
                    old_decision = part1;
                    new_decision = part2;
                }
                else
                {
                    old_decision = part2;
                    new_decision = part1;
                }
		
                bool not_in_cut = false;
                for (int k = row_ptr_inv[edge]; !not_in_cut && k < row_ptr_inv[edge + 1]; k++)
                {
                    int vertex = col_ind_inv[k];
                    int part_id = part_vec[vertex];
                    if (part_id == old_decision)
                    {
                        v_list_push[v_list_cntr++] = vertex;
                        for (int j = row_ptr[vertex]; !not_in_cut && j < row_ptr[vertex + 1]; j++)
                        {
                            int net_of_curr_vertex = col_ind[j];
                            if (ep_ends[net_of_curr_vertex] == 1)
                            { //Which is not in cut
                                not_in_cut = true;
                                break;
                            }
                        }
                    }
                }

                if (!not_in_cut)
                {
                    for (int ind = 0; ind < v_list_cntr; ind++)
                    {
                        int v_push = v_list_push[ind];
                        int part_id = part_vec[v_push];
                        part_vec[v_push] = new_decision;
                        part_sizes[part_id]--;
                        part_sizes[new_decision]++;

                        //REMOVAL
                        for (int i = row_ptr[v_push]; i < row_ptr[v_push + 1]; i++)
                        {
                            int edge = col_ind[i];
                            int start_ind = ep_starts[edge];
                            int end_cnt = ep_ends[edge];
                            for (int j = 0; j < end_cnt; j++)
                            {
                                if (ep_parts[start_ind + j] == part_id)
                                { //we find the part that we need to decrement a connection
                                    ep_cnt[start_ind + j]--;
                                    if (ep_cnt[start_ind + j] == 0)
                                    { //all connections are removed
                                        rem_cnt++;
                                        ep_parts[start_ind + j] = ep_parts[start_ind + end_cnt - 1]; //bring the part in the end to the deleted pos
                                        ep_cnt[start_ind + j] = ep_cnt[start_ind + end_cnt - 1];     //bring the count in the end to the deleted pos
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
                        for (int i = row_ptr[v_push]; i < row_ptr[v_push + 1]; i++) // Traversing the edge list of the v_push
                        {
                            int edge = col_ind[i];
                            int start_index = ep_starts[edge];
                            int end_cnt = ep_ends[edge];
                            bool found = false;
                            for (int j = 0; j < end_cnt; j++)
                            {
                                if (ep_parts[start_index + j] == new_decision)
                                { //we find the part that we need to increment a connection
                                    ep_cnt[start_index + j]++;
                                    found = true;
                                    break;
                                }
                            }
                            if (!found)
                            {
                                ep_parts[start_index + end_cnt] = new_decision;
                                ep_cnt[start_index + end_cnt] = 1;
                                ep_ends[edge]++;
                                inc_cnt++;
                            }
                        }
                    }
                }
            }
            //cout << "--------  "<< iter  << "  ---------"<<"\n after CN : " << CN(ep_ends, noe) << endl;
            //cout << "after TV : " << TV_Refined(ep_starts,ep_cnt,noe, K) << endl;
            //cout << "after TM : " << TM_scratch(K, part_vec, ep_parts, ep_starts, ep_ends, nov, noe, src_nodes) << endl << endl;
        }
    }
    end = chrono::steady_clock::now();
    diff = end - start;
    total += (chrono::duration<double, milli>(end - start));
    cout << "Total Paritioning took : " << chrono::duration<double, milli>(total).count() / 1000 << " seconds." << endl;
    cout << "Refinement took : " << chrono::duration<double, milli>(diff).count() / 1000 << " seconds." << endl;

    int TV_Refined_out = TV_Refined(ep_starts, ep_cnt, noe, K);
    int CN_Refined_out = CN(ep_ends, noe);
    int TM_Refined_out = TM_scratch(K, part_vec, ep_parts, ep_starts, ep_ends, nov, noe, src_nodes);
    int MGM_Refined_out = MGM_scratch(K, part_vec, ep_parts, ep_starts, ep_ends, nov, noe, src_nodes);
    int MGV_Refined_out = MGV_Scratch(part_vec, row_ptr, col_ind, src_nodes, ep_ends, net_appearence, K, nov, noe, pins);
    int MGAV_Refined_out = MGAV_Scratch(part_vec, row_ptr, col_ind, src_nodes, ep_ends, ep_parts, ep_starts, net_appearence, K, nov, noe, pins);



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

    cout << "Total Paritioning took : " << chrono::duration<double, milli>(total).count() / 1000 << " seconds." << endl;
    string out_name = f_name + ".fastpart." + to_string(K);
    ofstream outfile;
    outfile.open(out_name);
    for(int i = 0; i < nov; i++)
      {
	outfile << part_vec[i] << '\n';
      }
    outfile.close();
    std::cout << "Wrote part_vec to " << out_name << std::endl; 
}

void partition_par(unsigned int *row_ptr, unsigned int *col_ind, unsigned int *net_appearence, int nov, int noe, long long pins,int *edge_wgths, int *src_nodes, int K, int *ordering, string f_name, int REFI_OPT, int REFI_ITER, int NT,unsigned int *row_ptr_inv, unsigned int *col_ind_inv)
{
    int active_part_count;             // total number of active partitions while partitioning of a vertex
    int *marker = new int[K];          // stores the activeness information to determine active parts efficiently
    int *active_parts = new int[K];    // stores the active parts
    int *part_sizes = new int[K];      // stores number of verties in each part
    int *ep_starts = new int[noe + 1]; // degree of network
    int *ep_ends = new int[noe];       // how many parts that network have been split
    int *ep_parts = new int[pins];     // stores the edge-part info
    int *part_vec = new int[nov];      // stores the part id's of the vertices
    int *ep_cnt = new int[pins];
    double eps = 1.05;
    double C = (nov / K) * eps;
    double *scores = new double[K];
    bool check_selected_vertex;
    
    int partCount = K;

    //initializations
    auto start = chrono::steady_clock::now();
    initialize(ep_cnt, 0, pins, 0);
    initialize(marker, 0, K, -1);
    initialize(active_parts, 0, K, -1);
    initialize(part_sizes, 0, K, 0);
    initialize(scores, 0, K, 0.0);
    initialize(part_vec, 0, nov, -1);
    initialize(ep_starts, 0, noe + 1, 0);
    initialize(ep_ends, 0, noe, 0);
    initialize(ep_parts, 0, pins, -1);

    //Calculate edge degrees
    for (int i = 0; i < nov; i++)
    {
        for (int k = row_ptr[i]; k < row_ptr[i + 1]; k++)
        {
            int edge = col_ind[k];
            ep_starts[edge + 1]++;
        }
    }
    int max_edge_deg = -1;
    //Store edge degrees in a cumulative way
    for (int i = 2; i <= noe; i++)
    {
        ep_starts[i] = ep_starts[i - 1] + ep_starts[i];
        if(ep_starts[i - 1] + ep_starts[i] > max_edge_deg)
            max_edge_deg = ep_starts[i - 1] + ep_starts[i];
    }

    auto end = chrono::steady_clock::now();
    auto diff = end - start;
    cout << "Initialization took : " << chrono::duration<double, milli>(diff).count() / 1000 << " seconds." << endl;
    total += (chrono::duration<double, milli>(end - start));

    //Start Partitioning
    cout << "Partitioning started " << endl;
    start = chrono::steady_clock::now();
    for (int z = 0; z < nov; z++) //for each vertex
    {
        int vertex = ordering[z];
        int decision; // Stores the part decision for the <vertex>
        active_part_count = 0;

        for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++)
        {
            int edge = col_ind[i];
            if (ep_ends[edge] == 1) //information is valuable if only net is not cut
            {
                int part_id = ep_parts[ep_starts[edge]];
                if (marker[part_id] != vertex)
                {
                    scores[part_id] = 1;
                    marker[part_id] = vertex;
                    active_parts[active_part_count] = part_id;
                    active_part_count++;
                }
                else
                {
                    scores[part_id]++;
                }
            }
        }

        //Apply LDG Approach
        double max_val = 0;
        int max_index = -1;
        for (int i = 0; i < active_part_count; i++)
        {
            int part_id = active_parts[i];
            double penalty = 1 - (1.0 * part_sizes[part_id] / C);
            double score = scores[part_id] * penalty;
            if (score > 0)
            {
                if ((max_index == -1) || (score > max_val) || ((score == max_val) && (part_sizes[part_id] < part_sizes[max_index])))
                {
                    max_val = score;
                    max_index = part_id;
                }
            }
        }

        if (max_index == -1) //if every net is split or no net has been observed
        {
            //find min loaded
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
        part_sizes[decision]++;
        part_vec[vertex] = decision; //vertex part is decided

        //Update the edge part info
        for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++) // Traversing the edge list of the vertex
        {
            int edge = col_ind[i];
            int j;
            bool found = false;
            int end_cnt = ep_ends[edge];
            int start_index = ep_starts[edge];
            for (j = 0; j < end_cnt; j++)
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

    //REFINEMENT CODE STARTS
    start = chrono::steady_clock::now();
    int iter = -1;

    int **vertices_par = new int *[NT];
    int **decisions_par = new int *[NT];
    int *ref_cnt_glob = new int[NT];
    int *v_flag = new int[nov];
    for (int t = 0; t < NT; t++)
    {
        vertices_par[t] = new int[nov];
        decisions_par[t] = new int[nov];
    }
    omp_set_dynamic(0);
    omp_set_num_threads(NT);

    for (int iter = 0; REFI_OPT && iter < REFI_ITER; iter++)
    {
        int rem_cnt = 0;
        int inc_cnt = 0;
	
#pragma omp parallel
        {
	  int tid = omp_get_thread_num();
	  int ref_cnt = 0;
	  int *v_list_push = new int[max_edge_deg];
#pragma omp for schedule(dynamic, 256)
	  for (int edge = 0; edge < noe; edge += 1)
            {
	      if (ep_ends[edge] == 2)
                {
                    int v_list_cntr = 0;
                    int part1 = -1, part2 = -1;
                    int start_index = ep_starts[edge];
                    int end_cnt = ep_ends[edge];
                    // to find p1, p2
                    int part1_count, part2_count;
                    for (int j = 0; j < end_cnt; j++)
                    {
                        if (part1 != -1)
                        {
                            part2 = ep_parts[start_index + j];
                            part2_count = ep_cnt[start_index + j];
                            break;
                        }
                        else
                        {
                            part1 = ep_parts[start_index + j];
                            part1_count = ep_cnt[start_index + j];
                        }
                    }

                    // now we know p1 and p2
                    int old_decision, new_decision;
                    if (part1_count < part2_count)
                    {
                        old_decision = part1;
                        new_decision = part2;
                    }
                    else
                    {
                        old_decision = part2;
                        new_decision = part1;
                    }

		    
                    bool not_in_cut = false;
                    for (int k = row_ptr_inv[edge]; !not_in_cut && k < row_ptr_inv[edge + 1]; k++)
                    {
		      int vertex = col_ind_inv[k];
		      int part_id = part_vec[vertex];

		      if(edge == 540327)
			cout << "part_id: " << part_id << " old_decision: " << old_decision << endl;
		      
		      if (part_id == old_decision)
                        {
                            v_list_push[v_list_cntr++] = vertex;
                            //v_flag[vertex] = 1;
                            for (int j = row_ptr[vertex]; !not_in_cut && j < row_ptr[vertex + 1]; j++)
                            {
                                int net_of_curr_vertex = col_ind[j];
                                if (ep_ends[net_of_curr_vertex] == 1)
                                { //Which is not in cut
                                    not_in_cut = true;
                                    break;
                                }
                            }
                        }
                    }

		    if(edge == 540327)
		      cout << "Edge: " << 540327 << " v_list_cntr: " << v_list_cntr << " tid: " << tid <<endl;

                    if (!not_in_cut)
                    {
                        for (int ind = 0; ind < v_list_cntr; ind++)
                        {
                            int v_push = v_list_push[ind];
                            vertices_par[tid][ref_cnt] = v_push;
                            decisions_par[tid][ref_cnt] = new_decision;
                            ref_cnt++;
                        }
                    }
                }
            }

            ref_cnt_glob[tid] = ref_cnt;
#pragma omp single
            {

                for (int i = 0; i < NT; i++)
                {
                    //REMOVAL
                    for (int k = 0; k < ref_cnt_glob[i]; k++)
                    {
                        int v_push = vertices_par[i][k];
                        int part_id = part_vec[v_push];
                        int new_decision = decisions_par[i][k];

                        part_vec[v_push] = new_decision;
                        part_sizes[part_id]--;
                        part_sizes[new_decision]++;

                        for (int i = row_ptr[v_push]; i < row_ptr[v_push + 1]; i++)
                        {
                            int edge = col_ind[i];
                            int start_ind = ep_starts[edge];
                            int end_cnt = ep_ends[edge];
                            for (int j = 0; j < end_cnt; j++)
                            {
                                if (ep_parts[start_ind + j] == part_id)
                                { //we find the part that we need to decrement a connection
                                    ep_cnt[start_ind + j]--;
                                    if (ep_cnt[start_ind + j] == 0)
                                    { //all connections are removed
                                        rem_cnt++;
                                        ep_parts[start_ind + j] = ep_parts[start_ind + end_cnt - 1]; //bring the part in the end to the deleted pos
                                        ep_cnt[start_ind + j] = ep_cnt[start_ind + end_cnt - 1];     //bring the count in the end to the deleted pos
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
                        for (int i = row_ptr[v_push]; i < row_ptr[v_push + 1]; i++) // Traversing the edge list of the v_push
                        {
                            int edge = col_ind[i];
                            int start_index = ep_starts[edge];
                            int put_index = -1; // a possible location for newly introduced connection
                            int end_cnt = ep_ends[edge];
                            bool found = false;
                            for (int j = 0; j < end_cnt; j++)
                            {
                                if (ep_parts[start_index + j] == new_decision)
                                { //we find the part that we need to decrement a connection
                                    ep_cnt[start_index + j]++;
                                    found = true;
                                    break;
                                }
                            }
                            if (!found)
                            {
                                ep_parts[start_index + end_cnt] = new_decision;
                                ep_cnt[start_index + end_cnt] = 1;
                                ep_ends[edge]++;
                                inc_cnt++;
                            }
                        }
                    }
                } //nt
            }
        }
    }
    end = chrono::steady_clock::now();
    diff = end - start;
    total += (chrono::duration<double, milli>(end - start));
    cout << "Total Paritioning took : " << chrono::duration<double, milli>(total).count() / 1000 << " seconds." << endl;
    cout << "Refinement took : " << chrono::duration<double, milli>(diff).count() / 1000 << " seconds." << endl;

    int TV_Refined_out = TV_Refined(ep_starts, ep_cnt, noe, K);
    int CN_Refined_out = CN(ep_ends, noe);
    int TM_Refined_out = TM_scratch(K, part_vec, ep_parts, ep_starts, ep_ends, nov, noe, src_nodes);
    int MGM_Refined_out = MGM_scratch(K, part_vec, ep_parts, ep_starts, ep_ends, nov, noe, src_nodes);
    int MGV_Refined_out = MGV_Scratch(part_vec, row_ptr, col_ind, src_nodes, ep_ends, net_appearence, K, nov, noe, pins);
    int MGAV_Refined_out = MGAV_Scratch(part_vec, row_ptr, col_ind, src_nodes, ep_ends, ep_parts, ep_starts, net_appearence, K, nov, noe, pins);



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

    string out_name = f_name + ".fastpart." + to_string(K);
    ofstream outfile;
    outfile.open(out_name);
    for(int i = 0; i < nov; i++)
      {
	outfile << part_vec[i] << '\n';
      }
    outfile.close();
    std::cout << "Wrote part_vec to " << out_name << std::endl; 
    //writePartVector(f_name, nov, K, part_vec, "TV");
}


////////////#############################3

void thin_partition(unsigned int *row_ptr, unsigned int *col_ind, unsigned int *net_appearence, int nov, int noe, long long pins, int *edge_wgths, int *src_nodes, int K,
		    int *ordering, string &f_name, int &REFI_OPT, int &REFI_ITER, unsigned int *row_ptr_inv, unsigned int *col_ind_inv, int* part_vec)
{
    int active_part_count;          // total number of active partitions while partitioning of a vertex
    //int *part_vec = new int[nov];   // stores the part id's of the vertices
    int *marker = new int[K];       // stores the activeness information to determine active parts efficiently
    int *active_parts = new int[K]; // stores the active parts
    int *ep_starts = new int[nov];  // degree of network
    int *ep_ends = new int[nov];    // how many parts that network have been split
    int *ep_parts = new int[pins];  // stores the edge-part info
    int *part_sizes = new int[K];   // stores number of verties in each part
    double eps = 1.05;
    double C = (nov / K) * eps;
    double *scores = new double[K];
    int *ep_cnt = new int[pins];

    int partCount = K; 

    //initializations
    auto start = chrono::steady_clock::now();
    initialize(ep_cnt, 0, pins, 0);
    initialize(marker, 0, K, -1);
    initialize(active_parts, 0, K, -1);
    initialize(part_sizes, 0, K, 0);
    initialize(scores, 0, K, 0.0);
    initialize(part_vec, 0, nov, -1);
    initialize(ep_starts, 0, nov, 0); //checked
    initialize(ep_ends, 0, nov, 0);
    initialize(ep_parts, 0, pins, -1);

    //Calculate edge degrees
    for (int i = 0; i < nov; i++)
    {
        for (int k = row_ptr[i]; k < row_ptr[i + 1]; k++)
        {
            int edge = col_ind[k];
            ep_starts[edge + 1]++;
        }
    }
    //Store edge degrees in a cumulative way
    for (int i = 2; i <= nov; i++)
    {
        ep_starts[i] = ep_starts[i - 1] + ep_starts[i];
    }

    auto end = chrono::steady_clock::now();
    auto diff = end - start;
    cout << "Initialization took : " << chrono::duration<double, milli>(diff).count() / 1000 << " seconds." << endl;

    //Start Partitioning
    cout << "Partitioning started " << endl;
    start = chrono::steady_clock::now();
    for (int z = 0; z < nov; z++) //for each vertex
    {
        int vertex = ordering[z];
        int decision; // Stores the part decision for the <vertex>
        active_part_count = 0;

        for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++)
        {
            int edge = col_ind[i];
            if (ep_ends[edge] == 1) //information is valuable if only net is not cut
            {
                int part_id = ep_parts[ep_starts[edge]];
                if (marker[part_id] != vertex)
                {
                    scores[part_id] = 1;
                    marker[part_id] = vertex;
                    active_parts[active_part_count] = part_id;
                    active_part_count++;
                }
                else
                {
                    scores[part_id]++;
                }
            }
        }

        //Apply LDG Approach
        double max_val = 0;
        int max_index = -1;
        for (int i = 0; i < active_part_count; i++)
        {
            int part_id = active_parts[i];
            double penalty = 1 - (1.0 * part_sizes[part_id] / C);
            double score = scores[part_id] * penalty;
            if (score > 0)
            {
                if ((max_index == -1) || (score > max_val) || ((score == max_val) && (part_sizes[part_id] < part_sizes[max_index])))
                {
                    max_val = score;
                    max_index = part_id;
                }
            }
        }

        if (max_index == -1) //if every net is split or no net has been observed
        {
            //find min loaded
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
        part_sizes[decision]++;
        part_vec[vertex] = decision; //vertex part is decided

        //Update the edge part info
        for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++) // Traversing the edge list of the vertex
        {
            int edge = col_ind[i];
            int j;
            bool found = false;
            int end_cnt = ep_ends[edge];
            int start_index = ep_starts[edge];
            for (j = 0; j < end_cnt; j++)
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
    }

    end = chrono::steady_clock::now();
    diff = end - start;
    total += (chrono::duration<double, milli>(end - start));
    cout << "Paritioning took : " << chrono::duration<double, milli>(diff).count() / 1000 << " seconds." << endl;

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

    // REFINEMENT STARTS :))))))))
    start = chrono::steady_clock::now();

    for (int iter = 0; REFI_OPT && iter < REFI_ITER; iter++)
    {
        int *v_list_push = new int[nov];
        int v_list_cntr;
        int rem_cnt = 0;
        int inc_cnt = 0;

        for (int edge = 0; edge < noe; edge++)
        {
            if (ep_ends[edge] == 2)
            {
                v_list_cntr = 0;
                int part1 = -1, part2 = -1;
                int start_index = ep_starts[edge];
                int end_cnt = ep_ends[edge];
                // to find p1, p2
                int part1_count, part2_count;
                for (int j = 0; j < end_cnt; j++)

                {
                    if (part1 != -1)
                    {
                        part2 = ep_parts[start_index + j];
                        part2_count = ep_cnt[start_index + j];
                        break;
                    }
                    else
                    {
                        part1 = ep_parts[start_index + j];
                        part1_count = ep_cnt[start_index + j];
                    }
                }

                // now we know p1 and p2
                int old_decision, new_decision;
                if (part1_count < part2_count)
                {
                    old_decision = part1;
                    new_decision = part2;
                }
                else
                {
                    old_decision = part2;
                    new_decision = part1;
                }
		
                bool not_in_cut = false;
                for (int k = row_ptr_inv[edge]; !not_in_cut && k < row_ptr_inv[edge + 1]; k++)
                {
                    int vertex = col_ind_inv[k];
                    int part_id = part_vec[vertex];
                    if (part_id == old_decision)
                    {
                        v_list_push[v_list_cntr++] = vertex;
                        for (int j = row_ptr[vertex]; !not_in_cut && j < row_ptr[vertex + 1]; j++)
                        {
                            int net_of_curr_vertex = col_ind[j];
                            if (ep_ends[net_of_curr_vertex] == 1)
                            { //Which is not in cut
                                not_in_cut = true;
                                break;
                            }
                        }
                    }
                }

                if (!not_in_cut)
                {
                    for (int ind = 0; ind < v_list_cntr; ind++)
                    {
                        int v_push = v_list_push[ind];
                        int part_id = part_vec[v_push];
                        part_vec[v_push] = new_decision;
                        part_sizes[part_id]--;
                        part_sizes[new_decision]++;

                        //REMOVAL
                        for (int i = row_ptr[v_push]; i < row_ptr[v_push + 1]; i++)
                        {
                            int edge = col_ind[i];
                            int start_ind = ep_starts[edge];
                            int end_cnt = ep_ends[edge];
                            for (int j = 0; j < end_cnt; j++)
                            {
                                if (ep_parts[start_ind + j] == part_id)
                                { //we find the part that we need to decrement a connection
                                    ep_cnt[start_ind + j]--;
                                    if (ep_cnt[start_ind + j] == 0)
                                    { //all connections are removed
                                        rem_cnt++;
                                        ep_parts[start_ind + j] = ep_parts[start_ind + end_cnt - 1]; //bring the part in the end to the deleted pos
                                        ep_cnt[start_ind + j] = ep_cnt[start_ind + end_cnt - 1];     //bring the count in the end to the deleted pos
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
                        for (int i = row_ptr[v_push]; i < row_ptr[v_push + 1]; i++) // Traversing the edge list of the v_push
                        {
                            int edge = col_ind[i];
                            int start_index = ep_starts[edge];
                            int end_cnt = ep_ends[edge];
                            bool found = false;
                            for (int j = 0; j < end_cnt; j++)
                            {
                                if (ep_parts[start_index + j] == new_decision)
                                { //we find the part that we need to increment a connection
                                    ep_cnt[start_index + j]++;
                                    found = true;
                                    break;
                                }
                            }
                            if (!found)
                            {
                                ep_parts[start_index + end_cnt] = new_decision;
                                ep_cnt[start_index + end_cnt] = 1;
                                ep_ends[edge]++;
                                inc_cnt++;
                            }
                        }
                    }
                }
            }
            //cout << "--------  "<< iter  << "  ---------"<<"\n after CN : " << CN(ep_ends, noe) << endl;
            //cout << "after TV : " << TV_Refined(ep_starts,ep_cnt,noe, K) << endl;
            //cout << "after TM : " << TM_scratch(K, part_vec, ep_parts, ep_starts, ep_ends, nov, noe, src_nodes) << endl << endl;
        }
    }
    end = chrono::steady_clock::now();
    diff = end - start;
    total += (chrono::duration<double, milli>(end - start));
    cout << "Total Paritioning took : " << chrono::duration<double, milli>(total).count() / 1000 << " seconds." << endl;
    cout << "Refinement took : " << chrono::duration<double, milli>(diff).count() / 1000 << " seconds." << endl;

    int TV_Refined_out = TV_Refined(ep_starts, ep_cnt, noe, K);
    int CN_Refined_out = CN(ep_ends, noe);
    int TM_Refined_out = TM_scratch(K, part_vec, ep_parts, ep_starts, ep_ends, nov, noe, src_nodes);
    int MGM_Refined_out = MGM_scratch(K, part_vec, ep_parts, ep_starts, ep_ends, nov, noe, src_nodes);
    int MGV_Refined_out = MGV_Scratch(part_vec, row_ptr, col_ind, src_nodes, ep_ends, net_appearence, K, nov, noe, pins);
    int MGAV_Refined_out = MGAV_Scratch(part_vec, row_ptr, col_ind, src_nodes, ep_ends, ep_parts, ep_starts, net_appearence, K, nov, noe, pins);



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

    cout << "Total Paritioning took : " << chrono::duration<double, milli>(total).count() / 1000 << " seconds." << endl;
    string out_name = f_name + ".fastpart." + to_string(K);
    ofstream outfile;
    outfile.open(out_name);
    for(int i = 0; i < nov; i++)
      {
	outfile << part_vec[i] << '\n';
      }
    outfile.close();
    std::cout << "Wrote part_vec to " << out_name << std::endl; 
}

void thin_partition_par(unsigned int *row_ptr, unsigned int *col_ind, unsigned int *net_appearence, int nov, int noe, long long pins,int *edge_wgths, int *src_nodes, int K, int *ordering, string f_name, int REFI_OPT, int REFI_ITER, int NT,unsigned int *row_ptr_inv, unsigned int *col_ind_inv, int* part_vec)
{
    int active_part_count;             // total number of active partitions while partitioning of a vertex
    int *marker = new int[K];          // stores the activeness information to determine active parts efficiently
    int *active_parts = new int[K];    // stores the active parts
    int *part_sizes = new int[K];      // stores number of verties in each part
    int *ep_starts = new int[noe + 1]; // degree of network
    int *ep_ends = new int[noe];       // how many parts that network have been split
    int *ep_parts = new int[pins];     // stores the edge-part info
    //int *part_vec = new int[nov];      // stores the part id's of the vertices
    int *ep_cnt = new int[pins];
    double eps = 1.05;
    double C = (nov / K) * eps;
    double *scores = new double[K];
    bool check_selected_vertex;
    
    int partCount = K;

    //initializations
    auto start = chrono::steady_clock::now();
    initialize(ep_cnt, 0, pins, 0);
    initialize(marker, 0, K, -1);
    initialize(active_parts, 0, K, -1);
    initialize(part_sizes, 0, K, 0);
    initialize(scores, 0, K, 0.0);
    initialize(part_vec, 0, nov, -1);
    initialize(ep_starts, 0, noe + 1, 0);
    initialize(ep_ends, 0, noe, 0);
    initialize(ep_parts, 0, pins, -1);

    //Calculate edge degrees
    for (int i = 0; i < nov; i++)
    {
        for (int k = row_ptr[i]; k < row_ptr[i + 1]; k++)
        {
            int edge = col_ind[k];
            ep_starts[edge + 1]++;
        }
    }
    int max_edge_deg = -1;
    //Store edge degrees in a cumulative way
    for (int i = 2; i <= noe; i++)
    {
        ep_starts[i] = ep_starts[i - 1] + ep_starts[i];
        if(ep_starts[i - 1] + ep_starts[i] > max_edge_deg)
            max_edge_deg = ep_starts[i - 1] + ep_starts[i];
    }

    auto end = chrono::steady_clock::now();
    auto diff = end - start;
    cout << "Initialization took : " << chrono::duration<double, milli>(diff).count() / 1000 << " seconds." << endl;
    total += (chrono::duration<double, milli>(end - start));

    //Start Partitioning
    cout << "Partitioning started " << endl;
    start = chrono::steady_clock::now();
    for (int z = 0; z < nov; z++) //for each vertex
    {
        int vertex = ordering[z];
        int decision; // Stores the part decision for the <vertex>
        active_part_count = 0;

        for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++)
        {
            int edge = col_ind[i];
            if (ep_ends[edge] == 1) //information is valuable if only net is not cut
            {
                int part_id = ep_parts[ep_starts[edge]];
                if (marker[part_id] != vertex)
                {
                    scores[part_id] = 1;
                    marker[part_id] = vertex;
                    active_parts[active_part_count] = part_id;
                    active_part_count++;
                }
                else
                {
                    scores[part_id]++;
                }
            }
        }

        //Apply LDG Approach
        double max_val = 0;
        int max_index = -1;
        for (int i = 0; i < active_part_count; i++)
        {
            int part_id = active_parts[i];
            double penalty = 1 - (1.0 * part_sizes[part_id] / C);
            double score = scores[part_id] * penalty;
            if (score > 0)
            {
                if ((max_index == -1) || (score > max_val) || ((score == max_val) && (part_sizes[part_id] < part_sizes[max_index])))
                {
                    max_val = score;
                    max_index = part_id;
                }
            }
        }

        if (max_index == -1) //if every net is split or no net has been observed
        {
            //find min loaded
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
        part_sizes[decision]++;
        part_vec[vertex] = decision; //vertex part is decided

        //Update the edge part info
        for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++) // Traversing the edge list of the vertex
        {
            int edge = col_ind[i];
            int j;
            bool found = false;
            int end_cnt = ep_ends[edge];
            int start_index = ep_starts[edge];
            for (j = 0; j < end_cnt; j++)
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

    //REFINEMENT CODE STARTS
    start = chrono::steady_clock::now();
    int iter = -1;

    int **vertices_par = new int *[NT];
    int **decisions_par = new int *[NT];
    int *ref_cnt_glob = new int[NT];
    int *v_flag = new int[nov];
    for (int t = 0; t < NT; t++)
    {
        vertices_par[t] = new int[nov];
        decisions_par[t] = new int[nov];
    }
    omp_set_dynamic(0);
    omp_set_num_threads(NT);

    for (int iter = 0; REFI_OPT && iter < REFI_ITER; iter++)
    {
        int rem_cnt = 0;
        int inc_cnt = 0;
	
#pragma omp parallel
        {
	  int tid = omp_get_thread_num();
	  int ref_cnt = 0;
	  int *v_list_push = new int[max_edge_deg];
#pragma omp for schedule(dynamic, 256)
	  for (int edge = 0; edge < noe; edge += 1)
            {
	      if (ep_ends[edge] == 2)
                {
                    int v_list_cntr = 0;
                    int part1 = -1, part2 = -1;
                    int start_index = ep_starts[edge];
                    int end_cnt = ep_ends[edge];
                    // to find p1, p2
                    int part1_count, part2_count;
                    for (int j = 0; j < end_cnt; j++)
                    {
                        if (part1 != -1)
                        {
                            part2 = ep_parts[start_index + j];
                            part2_count = ep_cnt[start_index + j];
                            break;
                        }
                        else
                        {
                            part1 = ep_parts[start_index + j];
                            part1_count = ep_cnt[start_index + j];
                        }
                    }

                    // now we know p1 and p2
                    int old_decision, new_decision;
                    if (part1_count < part2_count)
                    {
                        old_decision = part1;
                        new_decision = part2;
                    }
                    else
                    {
                        old_decision = part2;
                        new_decision = part1;
                    }

		    
                    bool not_in_cut = false;
                    for (int k = row_ptr_inv[edge]; !not_in_cut && k < row_ptr_inv[edge + 1]; k++)
                    {
		      int vertex = col_ind_inv[k];
		      int part_id = part_vec[vertex];

		      if(edge == 540327)
			cout << "part_id: " << part_id << " old_decision: " << old_decision << endl;
		      
		      if (part_id == old_decision)
                        {
                            v_list_push[v_list_cntr++] = vertex;
                            //v_flag[vertex] = 1;
                            for (int j = row_ptr[vertex]; !not_in_cut && j < row_ptr[vertex + 1]; j++)
                            {
                                int net_of_curr_vertex = col_ind[j];
                                if (ep_ends[net_of_curr_vertex] == 1)
                                { //Which is not in cut
                                    not_in_cut = true;
                                    break;
                                }
                            }
                        }
                    }

		    if(edge == 540327)
		      cout << "Edge: " << 540327 << " v_list_cntr: " << v_list_cntr << " tid: " << tid <<endl;

                    if (!not_in_cut)
                    {
                        for (int ind = 0; ind < v_list_cntr; ind++)
                        {
                            int v_push = v_list_push[ind];
                            vertices_par[tid][ref_cnt] = v_push;
                            decisions_par[tid][ref_cnt] = new_decision;
                            ref_cnt++;
                        }
                    }
                }
            }

            ref_cnt_glob[tid] = ref_cnt;
#pragma omp single
            {

                for (int i = 0; i < NT; i++)
                {
                    //REMOVAL
                    for (int k = 0; k < ref_cnt_glob[i]; k++)
                    {
                        int v_push = vertices_par[i][k];
                        int part_id = part_vec[v_push];
                        int new_decision = decisions_par[i][k];

                        part_vec[v_push] = new_decision;
                        part_sizes[part_id]--;
                        part_sizes[new_decision]++;

                        for (int i = row_ptr[v_push]; i < row_ptr[v_push + 1]; i++)
                        {
                            int edge = col_ind[i];
                            int start_ind = ep_starts[edge];
                            int end_cnt = ep_ends[edge];
                            for (int j = 0; j < end_cnt; j++)
                            {
                                if (ep_parts[start_ind + j] == part_id)
                                { //we find the part that we need to decrement a connection
                                    ep_cnt[start_ind + j]--;
                                    if (ep_cnt[start_ind + j] == 0)
                                    { //all connections are removed
                                        rem_cnt++;
                                        ep_parts[start_ind + j] = ep_parts[start_ind + end_cnt - 1]; //bring the part in the end to the deleted pos
                                        ep_cnt[start_ind + j] = ep_cnt[start_ind + end_cnt - 1];     //bring the count in the end to the deleted pos
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
                        for (int i = row_ptr[v_push]; i < row_ptr[v_push + 1]; i++) // Traversing the edge list of the v_push
                        {
                            int edge = col_ind[i];
                            int start_index = ep_starts[edge];
                            int put_index = -1; // a possible location for newly introduced connection
                            int end_cnt = ep_ends[edge];
                            bool found = false;
                            for (int j = 0; j < end_cnt; j++)
                            {
                                if (ep_parts[start_index + j] == new_decision)
                                { //we find the part that we need to decrement a connection
                                    ep_cnt[start_index + j]++;
                                    found = true;
                                    break;
                                }
                            }
                            if (!found)
                            {
                                ep_parts[start_index + end_cnt] = new_decision;
                                ep_cnt[start_index + end_cnt] = 1;
                                ep_ends[edge]++;
                                inc_cnt++;
                            }
                        }
                    }
                } //nt
            }
        }
    }
    end = chrono::steady_clock::now();
    diff = end - start;
    total += (chrono::duration<double, milli>(end - start));
    cout << "Total Paritioning took : " << chrono::duration<double, milli>(total).count() / 1000 << " seconds." << endl;
    cout << "Refinement took : " << chrono::duration<double, milli>(diff).count() / 1000 << " seconds." << endl;

    int TV_Refined_out = TV_Refined(ep_starts, ep_cnt, noe, K);
    int CN_Refined_out = CN(ep_ends, noe);
    int TM_Refined_out = TM_scratch(K, part_vec, ep_parts, ep_starts, ep_ends, nov, noe, src_nodes);
    int MGM_Refined_out = MGM_scratch(K, part_vec, ep_parts, ep_starts, ep_ends, nov, noe, src_nodes);
    int MGV_Refined_out = MGV_Scratch(part_vec, row_ptr, col_ind, src_nodes, ep_ends, net_appearence, K, nov, noe, pins);
    int MGAV_Refined_out = MGAV_Scratch(part_vec, row_ptr, col_ind, src_nodes, ep_ends, ep_parts, ep_starts, net_appearence, K, nov, noe, pins);



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

    string out_name = f_name + ".fastpart." + to_string(K);
    ofstream outfile;
    outfile.open(out_name);
    for(int i = 0; i < nov; i++)
      {
	outfile << part_vec[i] << '\n';
      }
    outfile.close();
    std::cout << "Wrote part_vec to " << out_name << std::endl; 
    //writePartVector(f_name, nov, K, part_vec, "TV");
}

void partition_gpu(unsigned int *row_ptr, unsigned int *col_ind, unsigned int *net_appearence, int nov, int noe, long long pins,int *edge_wgths, int *src_nodes, int K, int *ordering, string f_name, int REFI_OPT, int REFI_ITER, unsigned int *row_ptr_inv, unsigned int *col_ind_inv)
{
    int active_part_count;          // total number of active partitions while partitioning of a vertex
    int *part_vec = new int[nov];   // stores the part id's of the vertices
    int *marker = new int[K];       // stores the activeness information to determine active parts efficiently
    int *active_parts = new int[K]; // stores the active parts
    int *ep_starts = new int[nov];  // degree of network
    int *ep_ends = new int[nov];    // how many parts that network have been split
    int *ep_parts = new int[pins];  // stores the edge-part info
    int *part_sizes = new int[K];   // stores number of verties in each part
    double eps = 1.05;
    double C = (nov / K) * eps;
    double *scores = new double[K];
    int *ep_cnt = new int[pins];

    int partCount = K; 

    //initializations
    auto start = chrono::steady_clock::now();
    initialize(ep_cnt, 0, pins, 0);
    initialize(marker, 0, K, -1);
    initialize(active_parts, 0, K, -1);
    initialize(part_sizes, 0, K, 0);
    initialize(scores, 0, K, 0.0);
    initialize(part_vec, 0, nov, -1);
    initialize(ep_starts, 0, nov, 0); //checked
    initialize(ep_ends, 0, nov, 0);
    initialize(ep_parts, 0, pins, -1);

    //Calculate edge degrees
    for (int i = 0; i < nov; i++)
    {
        for (int k = row_ptr[i]; k < row_ptr[i + 1]; k++)
        {
            int edge = col_ind[k];
            ep_starts[edge + 1]++;
        }
    }
    //Store edge degrees in a cumulative way
    for (int i = 2; i <= nov; i++)
    {
        ep_starts[i] = ep_starts[i - 1] + ep_starts[i];
    }

    auto end = chrono::steady_clock::now();
    auto diff = end - start;
    cout << "Initialization took : " << chrono::duration<double, milli>(diff).count() / 1000 << " seconds." << endl;

    //Start Partitioning
    cout << "Partitioning started " << endl;
    start = chrono::steady_clock::now();
    for (int z = 0; z < nov; z++) //for each vertex
    {
        int vertex = ordering[z];
        int decision; // Stores the part decision for the <vertex>
        active_part_count = 0;

        for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++)
        {
            int edge = col_ind[i];
            if (ep_ends[edge] == 1) //information is valuable if only net is not cut
            {
                int part_id = ep_parts[ep_starts[edge]];
                if (marker[part_id] != vertex)
                {
                    scores[part_id] = 1;
                    marker[part_id] = vertex;
                    active_parts[active_part_count] = part_id;
                    active_part_count++;
                }
                else
                {
                    scores[part_id]++;
                }
            }
        }

        //Apply LDG Approach
        double max_val = 0;
        int max_index = -1;
        for (int i = 0; i < active_part_count; i++)
        {
            int part_id = active_parts[i];
            double penalty = 1 - (1.0 * part_sizes[part_id] / C);
            double score = scores[part_id] * penalty;
            if (score > 0)
            {
                if ((max_index == -1) || (score > max_val) || ((score == max_val) && (part_sizes[part_id] < part_sizes[max_index])))
                {
                    max_val = score;
                    max_index = part_id;
                }
            }
        }

        if (max_index == -1) //if every net is split or no net has been observed
        {
            //find min loaded
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
        part_sizes[decision]++;
        part_vec[vertex] = decision; //vertex part is decided

        //Update the edge part info
        for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++) // Traversing the edge list of the vertex
        {
            int edge = col_ind[i];
            int j;
            bool found = false;
            int end_cnt = ep_ends[edge];
            int start_index = ep_starts[edge];
            for (j = 0; j < end_cnt; j++)
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
    }

    end = chrono::steady_clock::now();
    diff = end - start;
    total += (chrono::duration<double, milli>(end - start));
    cout << "Paritioning took : " << chrono::duration<double, milli>(diff).count() / 1000 << " seconds." << endl;

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

    // REFINEMENT STARTS :))))))))
    start = chrono::steady_clock::now();
    /////

    partition_gpu_wrapper(row_ptr, col_ind, net_appearence, nov, noe, pins, edge_wgths, src_nodes, K, ordering, REFI_OPT, REFI_ITER, row_ptr_inv, col_ind_inv, part_sizes, ep_starts, ep_parts, ep_ends,
			  ep_cnt, part_vec);
    /////
    end = chrono::steady_clock::now();
    diff = end - start;
    total += (chrono::duration<double, milli>(end - start));
    cout << "Total Paritioning took : " << chrono::duration<double, milli>(total).count() / 1000 << " seconds." << endl;
    cout << "Refinement took : " << chrono::duration<double, milli>(diff).count() / 1000 << " seconds." << endl;

    int TV_Refined_out = TV_Refined(ep_starts, ep_cnt, noe, K);
    int CN_Refined_out = CN(ep_ends, noe);
    int TM_Refined_out = TM_scratch(K, part_vec, ep_parts, ep_starts, ep_ends, nov, noe, src_nodes);
    int MGM_Refined_out = MGM_scratch(K, part_vec, ep_parts, ep_starts, ep_ends, nov, noe, src_nodes);
    int MGV_Refined_out = MGV_Scratch(part_vec, row_ptr, col_ind, src_nodes, ep_ends, net_appearence, K, nov, noe, pins);
    int MGAV_Refined_out = MGAV_Scratch(part_vec, row_ptr, col_ind, src_nodes, ep_ends, ep_parts, ep_starts, net_appearence, K, nov, noe, pins);



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

    cout << "Total Paritioning took : " << chrono::duration<double, milli>(total).count() / 1000 << " seconds." << endl;
    string out_name = f_name + ".fastpart." + to_string(K);
    ofstream outfile;
    outfile.open(out_name);
    for(int i = 0; i < nov; i++)
      {
	outfile << part_vec[i] << '\n';
      }
    outfile.close();
    std::cout << "Wrote part_vec to " << out_name << std::endl; 
}
