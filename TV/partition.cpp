#include "partition.h"
#include <omp.h>

#define C_FOCUS_VERTEX -1 //503000

extern void partition_gpu_wrapper(unsigned int* row_ptr,
				  unsigned int* col_ind,
				  unsigned int* net_appearence,
				  int nov,
				  int noe,
				  long long pins,
				  int* edge_wghts,
				  int* src_nodes,
				  int K,
				  int* ordering,
				  int* marker,
				  int* active_parts,
				  int* part_sizes,
				  int* ep_starts,
				  int* ep_ends,
				  int* ep_parts,
				  int* part_vec,
				  int* ep_cnt,
				  double* scores,
				  int REFI_ITER,
				  double C);

void LDG(bool &check_selected_vertex, int &K, int &active_part_count, int *part_sizes, int &decision, int *part_vec, int *active_parts, double &C, double *scores, int &vertex, int &nov)
{
  //CORRECTION
  int max_pw = (nov/K) + ((nov/K)*0.1);
  //289 for copapers 2048 parts
  //std::cout << "max_pw: " << max_pw << std::endl;
    
    //LDG STARTS
    double max_val = 0;
    int max_index = -1;

    for (int i = 0; i < active_part_count; i++)
    {
        int part_id = active_parts[i];
	double penalty = 1 - (1.0 * part_sizes[part_id] / C); //bu aslında coef gibi, daha kucuk partlari favor etmeye calisiyor
        double score = scores[part_id] * penalty;
	
	if(vertex == C_FOCUS_VERTEX){
	  printf("Vertex: %d, APC: %d, PartId: %d, PartSize: %d, Penalty: %f, Score: %f\n", vertex, active_part_count, part_id, part_sizes[part_id], penalty, score);
	} 
	
        if (score > 0)
        {
	  if ((max_index == -1) || (score > max_val) || ((score == max_val) && (part_sizes[part_id] < part_sizes[max_index]))) //NEGATIF GAIN'E IZIN VERMIYOR
            {
                max_val = score;
                max_index = part_id;
            }
        }
    }


    //BALANCE YAPARKEN NEGATIF MOVE'A IZIN VERIYOR
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

    //std::cout << "Dec pw: " << part_sizes[decision] << "\n";
    
    if(part_sizes[decision] > max_pw){
      std::cout << "Max pw " << max_pw << " can not be exceeded. Changed " << decision << " to " << part_vec[vertex] << std::endl;
      decision = part_vec[vertex];
    }
    
    if(vertex == C_FOCUS_VERTEX)
      std::cout << "Decision: " << decision << " pw: " << part_sizes[decision] << std::endl;
}

void LDG2(bool &check_selected_vertex, int &K, int &active_part_count, int *part_sizes, int &decision, int *part_vec, int *active_parts, double &C, double *scores, int &vertex, int refining, int &nov)
{
    //LDG STARTS
    double max_val = 0;
    int max_index = -1;

    int max_pw = (nov/K) + ((nov/K)*0.1);
    int min_pw = (nov/K) - ((nov/K)*0.1);

    for (int i = 0; i < active_part_count; i++)
    {
        int part_id = active_parts[i];
	//if(refining)
	//cout << "IN LDG -- VERTEX: " << vertex << " -- PartID: " << part_id << std::endl;
        double penalty = 1 - (1.0 * part_sizes[part_id] / C);
        double score = scores[part_id] * penalty;
	//if(part_sizes[part_id] > 500)
	//std::cout << "Penalty: " << penalty << " Score: " << score << " C: " << C <<
	//" part_sizes[part_id]: " << part_sizes[part_id] << std::endl;
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
    
    if(part_sizes[decision] > max_pw){
      //std::cout << "Max pw " << max_pw << " can not be exceeded. Changed " << decision << " to " << part_vec[vertex] << std::endl;
      decision = part_vec[vertex];
    }

    if(part_sizes[decision] < min_pw && part_vec[vertex] != -1){
      //std::cout << "Max pw " << min_pw << " can not be exceeded. Changed " << decision << " to " << part_vec[vertex] << std::endl;
      decision = part_vec[vertex];
    }
    
}

void partition(unsigned int *row_ptr, unsigned int *col_ind, unsigned int *net_appearence, int nov, int noe, long long pins, int *edge_wgths, int *src_nodes, int K, int *ordering, string f_name, int REFI_OPT, int REFI_ITER)
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
	  
	  if(edge < 0 || edge >= noe){
	    cout << "EDGE: " << edge << " noe: " << noe <<endl;
	    continue;
	  }
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
        LDG2(check_selected_vertex, K, active_part_count, part_sizes, decision, part_vec, active_parts, C, scores, vertex, 0, nov);
        part_sizes[decision]++;
        part_vec[vertex] = decision; // vertex part is decided
        // Update the edge part info
        for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++) // Traversing the edge list of the vertex
        {
            int edge = col_ind[i];
            bool found = false;
            int end_cnt = ep_ends[edge];
            int start_index = ep_starts[edge];
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

    int old = nov;
    //nov = 1000;

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
    int iter = -1;
    int cnt_good_dec = 0;
    for (iter = 0; REFI_OPT && iter < REFI_ITER; iter++)
    {
        int vertex_cnt = 0; //how many vertex is being alone in a partition
        int delta_tv_count = 0;
        for (int z = 0; z < nov; z++)
        {
            //GET SCORES STARTS
            int vertex = ordering[z];
            int decision;
            active_part_count = 0;
            int part_id = part_vec[vertex];
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
            LDG2(check_selected_vertex, K, active_part_count, part_sizes, decision, part_vec, active_parts, C, scores, vertex, 1, nov);
            // IF WE GET DIFFERENT PART

            if (part_id != decision)
            {
	      //std::cout << "vertex: " << vertex << " Part id: " << part_id << " Decision: " << decision << std::endl;
                int leave_gain = 0, arrival_loss = 0;
                // CHOOSE IF WE WILL GAIN OR LOSE AFTER PLACEMENT
                for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++) // Traversing the edge list of the vertex
                {
                    int edge = col_ind[i];
                    int start_index = ep_starts[edge];
                    int end_cnt = ep_ends[edge];
                    bool found = false;
                    for (int j = 0; j < end_cnt; j++)
                    {
                        if (ep_parts[start_index + j] == part_id && ep_cnt[start_index + j] == 1)
                        { //we find the part that we need to decrement a connection
                            leave_gain++;
                        }
                        else if (ep_parts[start_index + j] == decision)
                        {
                            found = true;
                        }
                    }
                    if (!found)
                    {
                        arrival_loss++;
                    }
                }
                // IF WE GAIN AFTER PLACEMENT
                if (leave_gain >= arrival_loss)
                {
                    cnt_good_dec++;
                    part_vec[vertex] = decision;
                    part_sizes[part_id]--;
                    part_sizes[decision]++;

                    //REMOVAL
                    for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++)
                    {
                        int edge = col_ind[i];
                        int start_ind = ep_starts[edge];
                        int end_cnt = ep_ends[edge];
                        int j = 0;
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
                                }
                                break; // we don't need to search no longer
                            }
                        }
                    }
                    //ADDING
                    for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++) // Traversing the edge list of the vertex
                    {

                        int edge = col_ind[i];
                        int start_index = ep_starts[edge];
                        int end_cnt = ep_ends[edge];
                        bool found = false;
                        for (int j = 0; j < end_cnt; j++)
                        {
                            if (ep_parts[start_index + j] == decision)
                            { //we find the part that we need to decrement a connection
                                ep_cnt[start_index + j]++;
                                found = true;
                                break;
                            }
                        }
                        if (!found)
                        {
                            ep_parts[start_index + end_cnt] = decision;
                            ep_cnt[start_index + end_cnt] = 1;
                            ep_ends[edge]++;
                            delta_tv_count++;
                        }
                    }
                }
            }
        }
    }
    //    cout << "good dec " << cnt_good_dec << endl;
    nov = old;
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
    printf("SEQ AFTER REFINEMENT ACTUAL IMBAL: %f \n", a_imbal);
    printf("#############################\n");

    //for(int i = 0; i < K; i++){
    //std::cout << "part: " << i << " size: " << part_sizes[i] << std::endl;
    //}
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
    double eps = 1.05;
    double C = (nov / K) * eps;
    double *scores = new double[K];
    bool check_selected_vertex;

    int max_pw = (nov/K) + ((nov/K)*0.1);

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
        LDG2(check_selected_vertex, K, active_part_count, part_sizes, decision, part_vec, active_parts, C, scores, vertex, 0, nov);

	if(part_sizes[decision] > max_pw)
	  continue;

	if(part_sizes[decision] > max_pw)
	  continue;
	  
        part_sizes[decision]++;
        part_vec[vertex] = decision; // vertex part is decided
        // Update the edge part info
        for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++) // Traversing the edge list of the vertex
        {
            int edge = col_ind[i];
            bool found = false;
            int end_cnt = ep_ends[edge];
            int start_index = ep_starts[edge];
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
    printf("PARALLEL BEFORE REFINEMENT ACTUAL IMBAL: %f \n", a_imbal);
    printf("#############################\n");

    //REFINEMENT CODE STARTS

    start = chrono::steady_clock::now();

    int **vertices_par = new int *[NT];
    int **decisions_par = new int *[NT];
    int *refined_cnt_glob = new int[NT];
    int **markers_glob = new int *[NT];
    int **active_part_glob = new int *[NT];
    double **scores_glob = new double *[NT];
    int* decisions = new int[nov]; //ADDITION
    int* part_before = new int[nov]; //ADDITION
    int cnt_good_dec = 0;
    omp_set_dynamic(0);
    omp_set_num_threads(NT);
    double single_time = 0;
    int move_cnt = 0;

    double ps = omp_get_wtime();
    //  random_shuffle (permutation.begin(), permutation.end());
    double pe = omp_get_wtime();
    cout << "shuffle took : " << pe - ps << endl;

#pragma omp parallel
    {
        int tid = omp_get_thread_num();

        markers_glob[tid] = new int[K];
        active_part_glob[tid] = new int[K];
        scores_glob[tid] = new double[K];
        vertices_par[tid] = new int[nov];
        decisions_par[tid] = new int[nov];
    }

    for (int iter = 0; REFI_OPT && iter < REFI_ITER; iter++)
    {
        //random_shuffle (permutation.begin(), permutation.end());
#pragma omp parallel
        {
            int tid = omp_get_thread_num();

            int *marker = markers_glob[tid];
            double *scores = scores_glob[tid];
            int *active_parts = active_part_glob[tid];
            int *my_vertices_par = vertices_par[tid];
            int *my_decisions_par = decisions_par[tid];
            int refined_cnt_loc;
            initialize(marker, 0, K, -1);
            int active_part_count;

            const int chunk_size = 8;
            int batch_size = chunk_size * NT;
	    
            for (int local_iter_cnt = 0; local_iter_cnt < (nov / batch_size) + 1; local_iter_cnt++)
	      {
                int start_point = local_iter_cnt * batch_size + chunk_size * tid;
                refined_cnt_loc = 0;
		
                for (int z = start_point; z < start_point + chunk_size && z < nov; z++)
		  {//BURADA BIR GARIPLIK VAR ORDER
		  //GET SCORES STARTS
		  int vertex = ordering[z];
		  int decision;
		  int part_id = part_vec[vertex];
		  active_part_count = 0;
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
                                active_parts[active_part_count++] = score_part_id;
			      }
                            else
			      {
                                scores[score_part_id]++;
			      }
                        }
                    }

		  //
		  check_selected_vertex = true;
		  
		  //for excluding the selected vertex
		  if (check_selected_vertex)
                    {
		      scores[part_id]--;
                    }
		  
		  LDG2(check_selected_vertex, K, active_part_count, part_sizes, decision, part_vec, active_parts, C, scores, vertex, 0, nov);
		  
		  decisions[vertex] = decision;
		  part_before[vertex] = part_id;
		  // IF WE GET DIFFERENT PART
		  if (part_id != decision)
                    {
		      int leave_gain = 0, arrival_loss = 0;
		      // CHOOSE IF WE WILL GAIN OR LOSE AFTER PLACEMENT
		      for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++) // Traversing the edge list of the vertex
                        {
			  int edge = col_ind[i];
			  int start_index = ep_starts[edge];
			  bool found = false;
			  int end_cnt = ep_ends[edge];
			  for (int j = 0; j < end_cnt; j++)
                            {
			      int target_part = ep_parts[start_index + j];
			      if (target_part == part_id && ep_cnt[start_index + j] == 1)
                                { //we find the part that we need to decrement a connection
				  leave_gain++;
				}
			      else if (target_part == decision)
				    {
				      found = true;
				    }
				}
			  if (!found)
                            {
			      arrival_loss++;
                            }
			    }
			  // IF WE GAIN AFTER PLACEMENT
			  if (leave_gain >= arrival_loss)
			    {
			      my_vertices_par[refined_cnt_loc] = vertex;
			      my_decisions_par[refined_cnt_loc] = decision;
			      refined_cnt_loc++;
			      //std::cout << "Part_sizes[part_id] was: " << part_sizes[part_vec[vertex]] << std::endl;
			      //std::cout << "Part_sizes[decision] was: " << part_sizes[decision] << std::endl;
			    }
			}
		    }
#pragma omp barrier
		  
                refined_cnt_glob[tid] = refined_cnt_loc;
		
#pragma omp single
                {
		  for (int i = 0; i < NT; i++)
                    {
		      for (int k = 0; k < refined_cnt_glob[i]; k++)
                        {
			  int vertex = vertices_par[i][k];
			  //int decision = decisions_par[i][k];
			  int decision = decisions[vertex];
			  //std::cout << "Decision: " << decision << std::endl;
			  int part_id = part_vec[vertex];
			  //****part_vec[vertex] = decision;
			  //part_vec[vertex] = decision;
			  //std::cout << "Part id: " << part_id << " Decision: " << decision << std::endl;
			  //std::cout << "Part_sizes[part_id] was: " << part_sizes[part_id] << std::endl;
			  //std::cout << "Part_sizes[decision] was: " << part_sizes[decision] << std::endl;
			  
			  //try{
			  //part_sizes[decision]++;}
			  //catch(...){
			  //std::cout << "catch: " << decision << std::endl;}
			  
			  if(part_sizes[decision] > max_pw){
			    std::cout << "Weird: " << part_sizes[decision] << std::endl;
			    break;
			  }
			  else{
			    part_vec[vertex] = decision;
			    part_sizes[part_id]--;
			  }
			  //std::cout << "Part_sizes[part_id] become: " << part_sizes[part_id] << std::endl;
			  //std::cout << "Part_sizes[decision] become: " << part_sizes[decision] << std::endl;
			  
			  
			  //REMOVAL
			  for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++)
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
                                        {                                                                //all connections are removed
					  ep_parts[start_ind + j] = ep_parts[start_ind + end_cnt - 1]; //bring the part in the end to the deleted pos
					  ep_cnt[start_ind + j] = ep_cnt[start_ind + end_cnt - 1];     //bring the count in the end to the deleted pos
					  ep_parts[start_ind + end_cnt - 1] = -1;
					  ep_cnt[start_ind + end_cnt - 1] = 0;
					  ep_ends[edge]--;
					  end_cnt--;
                                        }
				      break; // we don't need to search no longer
                                    }
                                }
                            }
			  //ADDING
			  for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++) // Traversing the edge list of the vertex
                            {
			      int edge = col_ind[i];
			      int start_index = ep_starts[edge];
			      bool found = false;
			      int end_cnt = ep_ends[edge];
			      for (int j = 0; j < end_cnt; j++)
                                {
				  if (ep_parts[start_index + j] == decision)
                                    { //we find the part that we need to decrement a connection
				      ep_cnt[start_index + j]++;
				      found = true;
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
                    }
                }
	      }
        }
    }
    // cout << "good dec " << cnt_good_dec << endl;
    //    cout << "single took : " << single_time << endl;
    
    end = chrono::steady_clock::now();
    diff = end - start;
    total += (chrono::duration<double, milli>(end - start));
    cout << "Total Partitioning took : " << chrono::duration<double, milli>(total).count() / 1000 << " seconds." << endl;
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
    
    //for(int i = 0; i < K; i++){
    //std::cout << "part: " << i << " size: " << part_sizes[i] << std::endl;
    //}

    //for (int i = 0; i < NT; i++)
    //{
    //for (int k = 0; k < refined_cnt_glob[i]; k++)
    //{
    //int vertext = vertices_par[i][k];
    //int decisiont = decisions_par[i][k];
    //std::cout << "Vertex: " << vertext << " ## ##" << " Decision: " << decisiont << std::endl;
    //}
    //}

    //for(int i = 0; i < nov; i++){
    //std::cout << "Vertex: " << i << " Decision: " << decisions[i] << " Before: " << part_before[i] <<std::endl;
    //}

    
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
    printf("PARALLEL AFTER REFINEMENT ACTUAL IMBAL: %f \n", a_imbal);
    printf("#############################\n");
    
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
        LDG2(check_selected_vertex, K, active_part_count, part_sizes, decision, part_vec, active_parts, C, scores, vertex, 0, nov);
        part_sizes[decision]++;
        part_vec[vertex] = decision; // vertex part is decided
        // Update the edge part info
        for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++) // Traversing the edge list of the vertex
        {
            int edge = col_ind[i];
            bool found = false;
            int end_cnt = ep_ends[edge];
            int start_index = ep_starts[edge];
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
    }

    end = chrono::steady_clock::now();
    diff = end - start;
    total += (chrono::duration<double, milli>(end - start));
    cout << "Partitioning took : " << chrono::duration<double, milli>(total).count() / 1000 << " seconds." << endl;

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
    
    /*-----------------------------------------------------*/
    /*-----------------------------------------------------*/
    /*TO THIS POINT, COPYING FROM CPU, WILL CALL GPU REFINE*/
    /*-----------------------------------------------------*/
    /*-----------------------------------------------------*/

    start = chrono::steady_clock::now();
    partition_gpu_wrapper(row_ptr, col_ind, net_appearence, nov, noe, pins, edge_wgths, src_nodes, K, ordering, marker, active_parts, part_sizes, ep_starts, ep_ends, ep_parts, part_vec, ep_cnt, scores, REFI_ITER, C);
    //cout << "good dec " << cnt_good_dec << endl;
    end = chrono::steady_clock::now();
    diff = end - start;
    total += (chrono::duration<double, milli>(end - start));
    cout << "Total Partitoning took : " << chrono::duration<double, milli>(total).count() / 1000 << " seconds." << endl;
    cout << "GPU Refinement took : " << chrono::duration<double, milli>(diff).count() / 1000 << " seconds." << endl;

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
    //for(int i = 0; i < K; i++){
    //std::cout << "part: " << i << " size: " << part_sizes[i] << std::endl;
    //}
    //double avg_weigth = 0.0;
    //for(int i = 0; i < K; i++){
    //avg_weigth += part_sizes[i];
    //}
    //avg_weigth /= K;

    //double a_imbal = 0;
    //for(int i = 0; i < K; i++) {
    //a_imbal = max(a_imbal, part_sizes[i]/avg_weigth);
    //}
    //std::cout << "ACTUAL IMBAL: " << a_imbal << std::endl;
}
