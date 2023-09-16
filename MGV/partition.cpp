#include "partition.h"
#include <omp.h>

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
				  int REFI_ITER, //Rest emerged from partitioning
				  int* send_volumes,
				  double* scores,
				  int* marker,
				  int* active_parts,
				  int* ep_starts,
				  int* ep_parts,
				  int* ep_ends,
				  int* ep_cnt,
				  int* part_sizes,
				  int* part_vec);


void LDG(bool &check_selected_vertex, int &K, int &active_part_count, int *part_sizes, int &decision, int *part_vec, int *active_parts, double &C, double *scores, int &vertex)
{
    //LDG STARTS
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


void partition(unsigned int *row_ptr, unsigned int *col_ind, unsigned int *net_appearence, int nov, int noe, long long pins,
               int *edge_wgths, int *src_nodes, int K, int *ordering, string f_name, int REFI_OPT, int REFI_ITER)
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
    int *send_volumes = new int[K];
    double eps = 1.05;
    double C = (nov / K) * eps;
    double *scores = new double[K];
    bool check_selected_vertex;

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
    initialize(send_volumes, 0, K, 0);

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
    for (int i = 2; i <= noe; i++)
    {
        ep_starts[i] = ep_starts[i - 1] + ep_starts[i];
    }

    auto end = chrono::steady_clock::now();
    auto diff = end - start;
    cout << "Initialization took : " << chrono::duration<double, milli>(diff).count() / 1000 << " seconds." << endl;
    total += (chrono::duration<double, milli>(end - start));

    //Start Partitioning
    cout << "Partitioning started " << endl;
    start = chrono::steady_clock::now();
    for (int z = 0; z < nov; z++)
    {
        int vertex = ordering[z];
        int decision;
        active_part_count = 0;

        //Get Scores
        for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++)
        {
            int edge = col_ind[i];
            for (int j = 0; j < ep_ends[edge]; j++)
            {
                int start_index = ep_starts[edge];
                int part_id = ep_parts[start_index + j];
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
        check_selected_vertex = false;
        //Apply LDG Approach
        LDG(check_selected_vertex, K, active_part_count, part_sizes, decision, part_vec, active_parts, C, scores, vertex);
        part_sizes[decision]++;
        part_vec[vertex] = decision; // vertex part is decided
        // Update the edge part info
        for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++) // Traversing the edge list of the vertex
        {
            int edge = col_ind[i];
            bool found = false;
            int end_cnt = ep_ends[edge];
            int start_index = ep_starts[edge];
            int net_src_part = part_vec[src_nodes[edge]];
            int net_src_node = src_nodes[edge];

            for (int j = 0; j < end_cnt; j++)
            {
                int part_id = ep_parts[start_index + j];
                if (net_src_node == vertex && part_id != net_src_part)
                    send_volumes[net_src_part]++;

                if (part_id == decision) // no dublicates
                {
                    found = true;
                    ep_cnt[start_index + j]++;
                    //break;
                }
            }

            if (!found)
            {
                ep_parts[start_index + end_cnt] = decision;
                ep_cnt[start_index + end_cnt] = 1;
                ep_ends[edge]++;

                //increase send_volumes
                if (src_nodes[edge] != vertex && net_src_part != -1)
                    send_volumes[net_src_part]++;
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
    int *max_part_change_detect = new int [K];
    int *gains = new int[K];
    int *losses = new int[K];
    int *max_parts = new int[K];
    int *max_part_marker = new int[K];
    initialize(gains, 0, K, 0);
    initialize(losses, 0, K, 0);
    initialize(max_part_marker, 0, K, -1);
    int move_cnt = 0;
    for (iter = 0; REFI_OPT && iter < REFI_ITER; iter++)
    {
        initialize(marker, 0, K, -1);
        int max_send_vol = -1;
        int max_part_count = 0;
        int max_changed_count = 0;
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
                                gains[net_src_part]++;
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
                            losses[net_src_part]++;
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
                    if (gains[max_part_id] < losses[max_part_id])
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

                    if (send_volumes[i] + (losses[i] - gains[i]) >= max_send_vol)
                    {
                        good_decision = false;
                    }
                }

                if (good_decision)
                {
                    best_part = decision;
                    move_cnt++;
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
                        }/////////////////////////

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
    cout << "Move cnt :" << move_cnt<<endl;
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

    cout << "Iterated for : " << iter << " times\n";
    writePartVector(f_name, nov, K, part_vec, "TV");
}



void partition_par(unsigned int *row_ptr, unsigned int *col_ind, unsigned int *net_appearence, int nov, int noe, long long pins,
                   int *edge_wgths, int *src_nodes, int K, int *ordering, string f_name, int REFI_OPT, int REFI_ITER, int NT)
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
    int *send_volumes = new int[K];
    double eps = 1.05;
    double C = (nov / K) * eps;
    double *scores = new double[K];
    bool check_selected_vertex;

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
    initialize(send_volumes, 0, K, 0);

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
    for (int i = 2; i <= noe; i++)
    {
        ep_starts[i] = ep_starts[i - 1] + ep_starts[i];
    }

    auto end = chrono::steady_clock::now();
    auto diff = end - start;
    cout << "Initialization took : " << chrono::duration<double, milli>(diff).count() / 1000 << " seconds." << endl;
    total += (chrono::duration<double, milli>(end - start));

    //Start Partitioning
    cout << "Partitioning started " << endl;
    start = chrono::steady_clock::now();
    for (int z = 0; z < nov; z++)
    {
        int vertex = ordering[z];
        int decision;
        active_part_count = 0;

        //Get Scores
        for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++)
        {
            int edge = col_ind[i];
            for (int j = 0; j < ep_ends[edge]; j++)
            {
                int start_index = ep_starts[edge];
                int part_id = ep_parts[start_index + j];
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
        check_selected_vertex = false;
        //Apply LDG Approach
        LDG(check_selected_vertex, K, active_part_count, part_sizes, decision, part_vec, active_parts, C, scores, vertex);
        part_sizes[decision]++;
        part_vec[vertex] = decision; // vertex part is decided
        // Update the edge part info
        for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++) // Traversing the edge list of the vertex
        {
            int edge = col_ind[i];
            bool found = false;
            int end_cnt = ep_ends[edge];
            int start_index = ep_starts[edge];
            int net_src_part = part_vec[src_nodes[edge]];
            int net_src_node = src_nodes[edge];

            for (int j = 0; j < end_cnt; j++)
            {
                int part_id = ep_parts[start_index + j];
                if (net_src_node == vertex && part_id != net_src_part)
                    send_volumes[net_src_part]++;

                if (part_id == decision) // no dublicates
                {
                    found = true;
                    ep_cnt[start_index + j]++;
                    //break;
                }
            }

            if (!found)
            {
                ep_parts[start_index + end_cnt] = decision;
                ep_cnt[start_index + end_cnt] = 1;
                ep_ends[edge]++;

                //increase send_volumes
                if (src_nodes[edge] != vertex && net_src_part != -1)
                    send_volumes[net_src_part]++;
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
    int **vertices_par = new int * [NT];
    int **decisions_par = new int * [NT];
    int *refined_cnt_glob = new int[NT];
    int ** gains_glob = new int * [NT];
    int ** losses_glob = new int * [NT];
    int ** markers_glob = new int * [NT];
    int ** active_part_glob = new int * [NT];
    double ** scores_glob = new double * [NT];

    int *max_parts = new int[K];
    int *max_part_marker = new int[K];
    int *max_part_change_detect = new int [K];

    omp_set_dynamic(0);
    omp_set_num_threads(NT);
    double single_time = 0;
    int move_cnt = 0;

#pragma omp parallel
{
    int tid = omp_get_thread_num();

    gains_glob[tid] = new int[K];
    losses_glob[tid] = new int[K];
    markers_glob[tid] = new int[K];
    active_part_glob[tid] = new int[K];
    scores_glob[tid] = new double[K];
    vertices_par[tid] = new int[nov];
    decisions_par[tid] = new int[nov];
}
start = chrono::steady_clock::now();

for (int iter = 0; REFI_OPT && iter < REFI_ITER; iter++)
{
    int max_send_vol = -1;
    int max_part_count = 0;
    int max_changed_count = 0;
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

#pragma omp parallel
        {
            int tid = omp_get_thread_num();

            int * gains = gains_glob[tid];
            int * losses = losses_glob[tid];
            int * marker = markers_glob[tid];
            int * active_parts = active_part_glob[tid];
            int * my_vertices_par = vertices_par[tid];
            int * my_decisions_par = decisions_par[tid];
            double * scores =  scores_glob[tid];
            int active_part_count;
            initialize(marker, 0, K, -1);

            const int chunk_size = 256;
            int batch_size = chunk_size * NT;


            for(int local_iter_cnt = 0; local_iter_cnt < (nov / batch_size) + 1; local_iter_cnt++)
            {
                int start_point = local_iter_cnt * batch_size + chunk_size * tid;
                int refined_cnt_loc = 0;

                for(int z = start_point; z < start_point + chunk_size && z < nov; z++)
                {
                    //GET SCORES STARTS
                    int decision;
                    int vertex = z;
                    active_part_count = 0;
                    int part_id = part_vec[vertex];
                    for(int j = 0; j < K; j++)
                    {
                        gains[j] = 0;
                        losses[j] = 0;
                    }

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
                                        gains[net_src_part]++;
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
                                    losses[net_src_part]++;
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
                            if (gains[max_part_id] < losses[max_part_id])
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

                            if (send_volumes[i] + (losses[i] - gains[i]) >= max_send_vol)
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
                            my_vertices_par[refined_cnt_loc] = vertex;
                            my_decisions_par[refined_cnt_loc] = decision;
                            refined_cnt_loc++;
                        }
                    }
                }
                #pragma omp barrier
                refined_cnt_glob[tid] = refined_cnt_loc;

                #pragma omp single
                {
                    bool max_changed_in_refinement = false;
                    bool max_send_changed_in_refinement = false;
                    for (int i = 0; !max_changed_in_refinement && i < NT; i++)
                    {
                        for (int k = 0; k < refined_cnt_glob[i]; k++)
                        {
                            int vertex = vertices_par[i][k];
                            int decision = decisions_par[i][k];
                            int part_id = part_vec[vertex];
                            part_vec[vertex] = decision;
                            part_sizes[part_id]--;
                            part_sizes[decision]++;
                            //cout << "init " << endl;
                            //REMOVING
                            for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++)
                            {
                                int edge = col_ind[i];
                                int net = edge;
                                int start_ind = ep_starts[edge];
                                int end_cnt = ep_ends[edge];
                                int net_src_node = src_nodes[edge];
                                int net_src_part = part_vec[net_src_node];
                                //cout << "init 1" << endl;
                                if (net_src_node == vertex)
                                {
                                    send_volumes[part_id] -= end_cnt - 1;
                                }
                                //cout << "src check 1" << endl;

                                for (int j = 0; j < end_cnt; j++)
                                {

                                    if (ep_parts[start_ind + j] == part_id)
                                    { //we find the part that we need to decrement a connection
                                        //cout << "gind part" << endl;
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
                                        //cout << "removed part " << endl;
                                        //break; // we don't need to search no longer
                                    }
                                }
                            }
                            //cout << "rem " << endl;

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


                            for(int i = 0; i < max_part_count; i++)
                            {
                                max_part_change_detect[i] = max_parts[i];
                            }
                                int prev_max_cnt = max_part_count;
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


                                for(int i = 0; i < max_part_count; i++)
                                {
                                    bool found = false;
                                    for(int j = 0; j < prev_max_cnt; j++)
                                    {
                                        if(max_part_change_detect[j] == max_parts[i])
                                        {
                                            found = true;
                                            break;
                                        }
                                    }
                                    if(!found)
                                        max_changed_in_refinement = true;
                                }


                            if(max_changed_in_refinement)
                                break;
                        }
                    }
                    //cout << "end" << endl;
                }
        }
    }//par
    }

    cout << "Move cnt :" << move_cnt<<endl;
    cout << "single took : " << single_time << endl;
    end = chrono::steady_clock::now();
    diff = end - start;
    total += (chrono::duration<double, milli>(end - start));
    cout << "Total Paritioning took : " << chrono::duration<double, milli>(total).count() / 1000 << " seconds." << endl;
    cout << "Refinement took : " << chrono::duration<double, milli>(diff).count() / 1000 << " seconds." << endl;

    int sum = 0;
    for(int k = 0; k < K; k++)
    {
        sum += send_volumes[k];
    }
    cout << "sum " << sum << endl;

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

    cout << "Iterated for : " << REFI_ITER << " times\n";
    writePartVector(f_name, nov, K, part_vec, "TV");
}

void partition_gpu(unsigned int *row_ptr, unsigned int *col_ind, unsigned int *net_appearence, int nov, int noe, long long pins, int *edge_wgths, int *src_nodes, int K, int *ordering, string f_name, int REFI_OPT, int REFI_ITER)
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
    int *send_volumes = new int[K];
    double eps = 1.05;
    double C = (nov / K) * eps;
    double *scores = new double[K];
    bool check_selected_vertex;

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
    initialize(send_volumes, 0, K, 0);

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
    for (int i = 2; i <= noe; i++)
    {
        ep_starts[i] = ep_starts[i - 1] + ep_starts[i];
    }

    auto end = chrono::steady_clock::now();
    auto diff = end - start;
    cout << "Initialization took : " << chrono::duration<double, milli>(diff).count() / 1000 << " seconds." << endl;
    total += (chrono::duration<double, milli>(end - start));

    //Start Partitioning
    cout << "Partitioning started " << endl;
    start = chrono::steady_clock::now();
    for (int z = 0; z < nov; z++)
    {
        int vertex = ordering[z];
        int decision;
        active_part_count = 0;

        //Get Scores
        for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++)
        {
            int edge = col_ind[i];
            for (int j = 0; j < ep_ends[edge]; j++)
            {
                int start_index = ep_starts[edge];
                int part_id = ep_parts[start_index + j];
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
        check_selected_vertex = false;
        //Apply LDG Approach
        LDG(check_selected_vertex, K, active_part_count, part_sizes, decision, part_vec, active_parts, C, scores, vertex);
        part_sizes[decision]++;
        part_vec[vertex] = decision; // vertex part is decided
        // Update the edge part info
        for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++) // Traversing the edge list of the vertex
        {
            int edge = col_ind[i];
            bool found = false;
            int end_cnt = ep_ends[edge];
            int start_index = ep_starts[edge];
            int net_src_part = part_vec[src_nodes[edge]];
            int net_src_node = src_nodes[edge];

            for (int j = 0; j < end_cnt; j++)
            {
                int part_id = ep_parts[start_index + j];
                if (net_src_node == vertex && part_id != net_src_part)
                    send_volumes[net_src_part]++;

                if (part_id == decision) // no dublicates
                {
                    found = true;
                    ep_cnt[start_index + j]++;
                    //break;
                }
            }

            if (!found)
            {
                ep_parts[start_index + end_cnt] = decision;
                ep_cnt[start_index + end_cnt] = 1;
                ep_ends[edge]++;

                //increase send_volumes
                if (src_nodes[edge] != vertex && net_src_part != -1)
                    send_volumes[net_src_part]++;
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
    
    partition_gpu_wrapper(row_ptr, col_ind, net_appearence, nov, noe, pins, edge_wgths, src_nodes, K, ordering, REFI_OPT, REFI_ITER, send_volumes, scores, marker, active_parts, ep_starts, ep_parts, ep_ends, ep_cnt, part_sizes, part_vec);
      
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

    writePartVector(f_name, nov, K, part_vec, "TV");
}

