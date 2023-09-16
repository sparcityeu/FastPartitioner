#include "partition.h"
void checkMe(unsigned int *row_ptr, unsigned int *col_ind, int nov, int noe, long long pins, int *ep_starts, int *ep_ends, int *ep_parts,
             int *ep_cnt, int *part_vec)
{
    int *t_ep_ends = new int[noe]; // how many parts have vertice from that net
    int *t_ep_cnt = new int[pins];
    int *t_ep_parts = new int[pins]; // stores the parts that nets are connected
    initialize(t_ep_ends, 0, nov, 0);
    initialize(t_ep_parts, 0, pins, -1);
    initialize(t_ep_cnt, 0, pins, 0);

    for (int vertex = 0; vertex < nov; vertex++)
    {
        //Update the edge part info
        int decision = part_vec[vertex];
        for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++)
        {
            int edge = col_ind[i];
            int end_cnt = t_ep_ends[edge];
            int start_index = ep_starts[edge];
            bool found = false;
            for (int j = 0; j < end_cnt; j++)
            {
                int part_id = t_ep_parts[start_index + j];
                if (part_id == decision) // no dublicates
                {
                    found = true;
                    t_ep_cnt[start_index + j]++;
                    break;
                }
            }

            if (!found)
            {
                t_ep_parts[start_index + end_cnt] = decision;
                t_ep_cnt[start_index + end_cnt] = 1;
                t_ep_ends[edge]++;
            }
        }
    }

    for (int i = 0; i < noe; i++)
    {
        int edge = i;
        int end_cnt = ep_ends[edge];
        int start_index = ep_starts[edge];
        for (int j = 0; j < end_cnt; j++)
        {
            int part_id = ep_parts[start_index + j];
            int t_end_cnt = t_ep_ends[i];
            bool found = false;
            for (int z = 0; z < t_end_cnt; z++)
            {
                if (part_id == t_ep_parts[start_index + z])
                {
                    if (ep_cnt[start_index + j] != t_ep_cnt[start_index + z])
                        cout << "cnt er" << endl;
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                cout << "PART ERROR" << endl;
            }
        }
    }
}
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
        part_connections[i][i] = -1;
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
    cout << "Paritioning took : " << chrono::duration<double, milli>(diff).count() / 1000 << " seconds." << endl;
    total += (chrono::duration<double, milli>(end - start));

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

    int *gains = new int[K];
    int *losses = new int[K];
    int *max_parts = new int[K];
    int *max_part_marker = new int[K];
    int **gain_loss_2d = new int *[K];
    for (int i = 0; i < K; i++)
    {
        gain_loss_2d[i] = new int[K];
        for (int j = 0; j < K; j++)
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
        for (int v = 0; v < nov; v++)
        {
            int best_part = -1, best_improvement = -1;
            int part_id = part_vec[v];
            initialize(gains, 0, K, 0);
            initialize(max_part_marker, 0, K, -1);
            bool has_gain = false;
            receiver_cnt = 0;
            sender_cnt = 0;
            //decide maximum parts
            int max_mes = -1;
            int max_part_count = 0;
            int max_changed_count = 0;
            for (int i = 0; i < K; i++)
            {
                if (part_messages[i] > max_mes)
                {
                    max_part_count = 1;
                    max_parts[0] = i;
                    max_mes = part_messages[i];
                    max_changed_count++;
                    max_part_marker[i] = max_changed_count;
                }
                else if (part_messages[i] == max_mes)
                {
                    max_parts[max_part_count++] = i;
                    max_part_marker[i] = max_changed_count;
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
                            if (ref_marker_receiver[target_part] != v)
                            {
                                receivers[receiver_cnt++] = target_part;
                                ref_marker_receiver[target_part] = v;
                            }

                            gain_loss_2d[part_id][target_part]++;

                            if (part_connections[part_id][target_part] == 1 ||
                                part_connections[part_id][target_part] == gain_loss_2d[part_id][target_part])
                            {
                                gains[part_id]++;
                                has_gain = true;
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
                                    gains[p_sender]++;
                                    has_gain = true;
                                }
                            }
                        }
                    }
                }
            }

            //FIND BETTER PART
            if (has_gain)
            {
                for (int k = 0; k < K; k++)
                {
                    if (k == part_id)
                        continue;

                    initialize(losses, 0, K, 0);
                    for (int i = 0; i < receiver_cnt; i++)
                    {
                        int rec_part = receivers[i];
                        if (part_connections[k][rec_part] == 0)
                        {
                            losses[k]++;
                        }
                    }

                    for (int i = 0; i < sender_cnt; i++)
                    {
                        int sender_part = senders[i];
                        if (part_connections[sender_part][k] == 0)
                        {
                            losses[sender_part]++;
                        }
                    }

                    bool good_decision = true;
                    int improvement = 0;
                    for (int i = 0; i < max_part_count; i++)
                    {
                        int max_part_id = max_parts[i];
                        if (gains[max_part_id] <= losses[max_part_id])
                        {
                            good_decision = false;
                            break;
                        }
                        improvement += gains[max_part_id] - losses[max_part_id];
                    }

                    for (int i = 0; good_decision && i < K; i++)
                    {
                        if (max_part_marker[i] == max_changed_count)
                        {
                            continue;
                        }

                        if (part_messages[i] + (losses[i] - gains[i]) >= max_mes)
                        {

                            good_decision = false;
                        }
                    }

                    if (good_decision && (improvement > best_improvement))
                    {
                        best_part = k;
                        best_improvement = improvement;
                    }
                }
            }

            //ADD AND REMOVE
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
        auto bizim_end = chrono::steady_clock::now();
        auto end2 = chrono::steady_clock::now();
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
    cout << "Iterated for : " << REFI_ITER << " times\n";
}

void partition_par(unsigned int *row_ptr, unsigned int *col_ind, unsigned int *net_appearence, int nov, int noe, long long pins,
                   int *edge_wgths, int *src_nodes, int K, int *ordering, string f_name, int REFI_OPT, int REFI_ITER, int NT)
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
        part_connections[i][i] = -1;
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
    cout << "Paritioning took: " << chrono::duration<double, milli>(diff).count() / 1000 << " seconds." << endl;
    total += (chrono::duration<double, milli>(end - start));

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
    int **vertices_par = new int *[NT];
    int **decisions_par = new int *[NT];
    int *refined_cnt_glob = new int[NT];
    int ** gains_glob = new int * [NT];
    int ** losses_glob = new int * [NT];
    int *** gain_loss_2d_glob = new int ** [NT];
    int ** ref_marker_receiver_glob = new int * [NT];
    int ** ref_marker_sender_glob = new int * [NT];
    int ** receivers_glob = new int * [NT];
    int ** senders_glob = new int * [NT];
    int *max_parts = new int[K];
    int *max_part_marker = new int[K];
    int *max_part_change_detect = new int [K];

    omp_set_dynamic(0);
    omp_set_num_threads(NT);

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();

        gains_glob[tid] = new int[K];
        losses_glob[tid] = new int[K];
        vertices_par[tid] = new int[nov];
        decisions_par[tid] = new int[nov];
        gain_loss_2d_glob[tid] = new int *[K];
        ref_marker_sender_glob[tid] = new int [K];
        ref_marker_receiver_glob[tid] = new int [K];
        senders_glob[tid] = new int [K];
        receivers_glob[tid] = new int [K];

        for (int i = 0; i < K; i++)
        {
            gain_loss_2d_glob[tid][i] = new int [K];
            for(int j = 0; j < K; j++)
            {
                gain_loss_2d_glob[tid][i][j] = 0;
            }
        }
    }

for (int iter = 0; REFI_OPT && iter < REFI_ITER; iter++)
    {
        int max_mes = -1;
        int max_part_count = 0;
        int max_changed_count = 0;
        initialize(max_part_marker, 0, K, -1);

        for (int i = 0; i < K; i++)
        {
            if (part_messages[i] > max_mes)
            {
                max_part_count = 1;
                max_parts[0] = i;
                max_mes = part_messages[i];
                max_changed_count++;
                max_part_marker[i] = max_changed_count;
            }
            else if (part_messages[i] == max_mes)
            {
                max_parts[max_part_count++] = i;
                max_part_marker[i] = max_changed_count;
            }
        }

        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            int refined_cnt_loc = 0;
            int receiver_cnt, sender_cnt, gain, loss;
            int *gains = gains_glob[tid];
            int *losses = losses_glob[tid];
            int *ref_marker_receiver = ref_marker_receiver_glob[tid];
            int *ref_marker_sender = ref_marker_sender_glob[tid];
            int *receivers = receivers_glob[tid];
            int *senders = senders_glob[tid];
            int **gain_loss_2d = gain_loss_2d_glob[tid];

            initialize(ref_marker_receiver, 0, K, -1);
            initialize(ref_marker_sender, 0, K, -1);

            const int chunk_size = 128;
            int batch_size = chunk_size * NT;

            for (int local_iter_cnt = 0; local_iter_cnt < (nov / batch_size) + 1; local_iter_cnt++)
            {
                int start_point = local_iter_cnt * batch_size + chunk_size * tid;
                refined_cnt_loc = 0;
                for (int z = start_point; z < start_point + chunk_size && z < nov; z++)
                {
                    //GET SCORES STARTS
                    int vertex = z;
                    int v = vertex;
                    int decision;
                    int part_id = part_vec[vertex];
                    int best_part = -1, best_improvement = -1;
                    initialize(gains, 0, K, 0);
                    bool has_gain = false;
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
                                    if (ref_marker_receiver[target_part] != v)
                                    {
                                        receivers[receiver_cnt++] = target_part;
                                        ref_marker_receiver[target_part] = v;
                                    }

                                    gain_loss_2d[part_id][target_part]++;

                                    if (part_connections[part_id][target_part] == 1 ||
                                        part_connections[part_id][target_part] == gain_loss_2d[part_id][target_part])
                                    {
                                        gains[part_id]++;
                                        has_gain = true;
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
                                            gains[p_sender]++;
                                            has_gain = true;
                                        }
                                    }
                                }
                            }
                        }
                    }

                    //FIND BETTER PART
                    if (has_gain)
                    {
                        for (int k = 0; k < K; k++)
                        {
                            if (k == part_id)
                                continue;

                            initialize(losses, 0, K, 0);
                            for (int i = 0; i < receiver_cnt; i++)
                            {
                                int rec_part = receivers[i];
                                if (part_connections[k][rec_part] == 0)
                                {
                                    losses[k]++;
                                }
                            }

                            for (int i = 0; i < sender_cnt; i++)
                            {
                                int sender_part = senders[i];
                                if (part_connections[sender_part][k] == 0)
                                {
                                    losses[sender_part]++;
                                }
                            }

                            bool good_decision = true;
                            int improvement = 0;
                            for (int i = 0; i < max_part_count; i++)
                            {
                                int max_part_id = max_parts[i];
                                if (gains[max_part_id] <= losses[max_part_id])
                                {
                                    good_decision = false;
                                    break;
                                }
                                improvement += gains[max_part_id] - losses[max_part_id];
                            }

                            for (int i = 0; good_decision && i < K; i++)
                            {
                                if (max_part_marker[i] == max_changed_count)
                                {
                                    continue;
                                }

                                if (part_messages[i] + (losses[i] - gains[i]) >= max_mes)
                                {

                                    good_decision = false;
                                }
                            }

                            if (good_decision && (improvement > best_improvement))
                            {
                                best_part = k;
                                best_improvement = improvement;
                            }
                        }
                    }

                    //ADD AND REMOVE
                    if (best_part != -1)
                    {
                        vertices_par[tid][refined_cnt_loc] = v;
                        decisions_par[tid][refined_cnt_loc] = best_part;
                        refined_cnt_loc++;
                    }
                }

                #pragma omp barrier

                refined_cnt_glob[tid] = refined_cnt_loc;

                #pragma omp single
                {
                    bool max_changed_in_refinement = false;
                    bool max_send_changed_in_refinement = false;
                    for (int i = 0; i < NT; i++)
                    {
                        for (int k = 0; k < refined_cnt_glob[i]; k++)
                        {
                            int vertex = vertices_par[i][k];
                            int v = vertex;
                            int decision = decisions_par[i][k];
                            int part_id = part_vec[vertex];
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
                                        }
                                        break; // we don't need to search no longer
                                    }
                                }
                            }
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

                            for(int i = 0; i < max_part_count; i++)
                            {
                                max_part_change_detect[i] = max_parts[i];
                            }
                                int prev_max_cnt = max_part_count;
                                max_part_count = 0;
                                max_mes = -1;
                                max_changed_count = 0;
                                initialize(max_part_marker, 0, K, -1);

                                for (int i = 0; i < K; i++)
                                {
                                    if (part_messages[i] > max_mes)
                                    {
                                        max_part_count = 1;
                                        max_parts[0] = i;
                                        max_mes = part_messages[i];
                                        max_changed_count++;
                                        max_part_marker[i] = max_changed_count;
                                    }
                                    else if (part_messages[i] == max_mes)
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
                }

            }

            auto bizim_end = chrono::steady_clock::now();
            auto end2 = chrono::steady_clock::now();
        }
    } //pragma

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

    cout << "Iterated for : " << REFI_ITER << " times\n";
}
