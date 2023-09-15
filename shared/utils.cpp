#include "utils.h"

using namespace std;

void createInverse_patohlike(unsigned int *row_ptr, unsigned int *col_ind, unsigned int *vertex_appearence, int nov, int noe, int pins, unsigned int **row_ptr_inv_p, unsigned int **col_ind_inv_p)
{
    auto start = chrono::steady_clock::now();
    //cout << "Creating inverse graph." << endl;
    *row_ptr_inv_p = new unsigned int[nov + 1];
    *col_ind_inv_p = new unsigned int[pins];
    (*row_ptr_inv_p)[0] = 0;

    for (int i = 1; i < nov + 1; i++)
    {
        (*row_ptr_inv_p)[i] = (*row_ptr_inv_p)[i - 1] + vertex_appearence[i - 1];
    }
    int *v_cnt = new int[nov];
    initialize(v_cnt, 0, noe, 0);
    for (int i = 0; i < noe; i++)
    {
        for (int k = row_ptr[i]; k < row_ptr[i + 1]; k++)
        {
            int vertex = col_ind[k];
            (*col_ind_inv_p)[(*row_ptr_inv_p)[vertex] + v_cnt[vertex]] = i;
            v_cnt[vertex]++;
        }
    }
    auto end = chrono::steady_clock::now();
    //cout << "Creating inverse graph took: " << chrono::duration<double, milli>(end - start).count() / 1000 << " seconds." << endl;
}

void createInverse(unsigned int *row_ptr, unsigned int *col_ind, unsigned int *net_appearence, int nov, int noe, int pins, unsigned int **row_ptr_inv_p, unsigned int **col_ind_inv_p)
{
    auto start = chrono::steady_clock::now();
    //cout << "Creating inverse graph." << endl;
    *row_ptr_inv_p = new unsigned int[noe + 1];
    *col_ind_inv_p = new unsigned int[pins];
    (*row_ptr_inv_p)[0] = 0;
    for (int i = 1; i < noe + 1; i++)
    {
        (*row_ptr_inv_p)[i] = (*row_ptr_inv_p)[i - 1] + net_appearence[i - 1];
    }
    int *net_cnt = new int[noe];
    initialize(net_cnt, 0, noe, 0);
    for (int i = 0; i < nov; i++)
    {
        for (int k = row_ptr[i]; k < row_ptr[i + 1]; k++)
        {
            int net = col_ind[k];
            (*col_ind_inv_p)[(*row_ptr_inv_p)[net] + net_cnt[net]] = i;
            net_cnt[net]++;
        }
    }
    auto end = chrono::steady_clock::now();
    //cout << "Creating inverse graph took: " << chrono::duration<double, milli>(end - start).count() / 1000 << " seconds." << endl;
}

void colorCoarse(unsigned int *row_ptr_inv, unsigned int *col_ind_inv, int nov, int noe, int no_color)
{
  int* vertex_color = new int[nov];
  int* color_tracker = new int[no_color];
  std::default_random_engine generator;
  std::uniform_int_distribution<int> distribution(0, no_color - 1);
  for(int i = 0; i < nov; i++)
    {
      int color = distribution(generator);
      vertex_color[i] = color;
    }
  
  int max_color;
  for(int i = 0; i < noe; i++)
    {
      max_color = 0;
      for(int k = row_ptr_inv[i]; k < row_ptr_inv[i + 1]; k++)
	{
	  int vertex = col_ind_inv[k];
	  color_tracker[vertex_color[vertex]] += 1;
	  if(color_tracker[vertex_color[vertex]] > color_tracker[max_color])
	    max_color = vertex_color[vertex]; 
	}
      for(int k = row_ptr_inv[i]; k < row_ptr_inv[i + 1]; k++)
	{
	  int vertex = col_ind_inv[k];
	  vertex_color[vertex] = max_color;
	}
      for(int j = 0; j < no_color; j++)
	color_tracker[j] = 0; 
	
    }

  for(int i = 0; i < nov; i++)
    {
      color_tracker[vertex_color[i]] += 1;
    }
  
  for(int j = 0; j < no_color; j++)
    {
      std::cout << "Color " << j << ": " << color_tracker[j] << std::endl; 
    }

  delete[] color_tracker;
}

int bfsFarthest(int start_node, unsigned int *row_ptr, unsigned int *col_ind, unsigned int *row_ptr_inv, unsigned int *col_ind_inv,
                int nov, int *distance, int NT)
{
    //Apply BFS
    for (int i = 0; i < nov; i++)
        distance[i] = -1;     // all nodes are unreachable in the beginning
    distance[start_node] = 0; //starting node is reachable from itself
    int step = 0;             //number of steps = breadth level
    bool changed = true;      //if we did not change anything in distance array we can stop
    while (changed)
    {
        changed = false;
#pragma omp parallel for num_threads(NT)
        for (int i = 0; i < nov; i++)
        {
            if (distance[i] == step)
            {
                for (int k = row_ptr[i]; k < row_ptr[i + 1]; k++)
                {
                    int net = col_ind[k];
                    for (int j = row_ptr_inv[net]; j < row_ptr_inv[net + 1]; j++)
                    {
                        int pin = col_ind_inv[j];
                        if (distance[pin] == -1)
                        {
                            distance[pin] = step + 1;
                            changed = true;
                        }
                    }
                }
            }
        }
        step++;
        //cout << step << endl;
    }

    //Determine the farthest node and apply bfs on that node one more time
    int farthest = start_node; //we are sure that the node itself is reachable by distance 0
    int farthest_distance = 0;
    bool unreachable = false; //if some distance is -1 it means that node is unreachable
    for (int i = 0; i < nov; i++)
    {
        if (distance[i] > farthest_distance)
        {
            farthest = i;
            farthest_distance = distance[i];
        }
        else if (distance[i] == -1)
        {
            unreachable = true;
        }
    }
    if (unreachable)
    {
        cout << "There are unreachable nodes." << endl;
    }
    return farthest;
}

void bfsOrdering(int start_node, unsigned int *row_ptr, unsigned int *col_ind, unsigned int *row_ptr_inv, unsigned int *col_ind_inv,
                 int nov, int *distance, int *ordering, int NT)
{
    ordering[0] = start_node;
    int ordering_cnt = 1;
    //Apply BFS
    for (int i = 0; i < nov; i++)
        distance[i] = -1;     // all nodes are unreachable in the beginning
    distance[start_node] = 0; //starting node is reachable from itself
    int step = 0;             //number of steps = breadth level
    bool changed = true;      //if we did not change anything in distance array we can stop
    int *v_ordered_check = new int[nov];
    initialize(v_ordered_check, 0, nov, 0);
    while (changed)
    {
        changed = false;
#pragma omp parallel for num_threads(NT)
        for (int i = 0; i < nov; i++)
        {
            if (distance[i] == step)
            {
                for (int k = row_ptr[i]; k < row_ptr[i + 1]; k++)
                {
                    int net = col_ind[k];
                    for (int j = row_ptr_inv[net]; j < row_ptr_inv[net + 1]; j++)
                    {
                        int pin = col_ind_inv[j];
                        if (distance[pin] == -1)
                        {
                            distance[pin] = step + 1;
                            changed = true;
                        }
                    }
                }
            }
        }
        step++;
    }

    int *step_cnt = new int[nov];
    initialize(step_cnt, 0, nov, 0);
    for (int k = 0; k < nov; k++)
    {
        step_cnt[distance[k]]++;
    }
    int max_dist = step;

    for (int i = 1; i < max_dist; i++)
    {
        int cnt = 0;
        for (int k = 0; k < nov; k++)
        {
            if (distance[k] == i)
            {
                {
                    ordering[ordering_cnt++] = k;
                    cnt++;
                }
                if (cnt == step_cnt[i])
                    break;
            }
        }
    }

    cout << ordering_cnt << "ORD" << endl;
    if (ordering_cnt != nov)
    {
        cout << "There are unreachable nodes." << endl;
        for (int i = 0; i < nov; i++)
        {
            if (distance[i] == -1)
                ordering[ordering_cnt++] = i;
        }
    }
}

int bfsFarthest_queue(int start_node, unsigned int *row_ptr, unsigned int *col_ind, unsigned int *row_ptr_inv, unsigned int *col_ind_inv, int nov, int *distance)
{
    int *queue = new int[nov];
    int front = 0, rear = 0;
    //Apply BFS
    for (int i = 0; i < nov; i++)
    {
        distance[i] = -1; // all nodes are unreachable in the beginning
        queue[i] = -1;
    }
    distance[start_node] = 0; //starting node is reachable from itself
    queue[front] = start_node;
    while (front <= rear)
    {
        int v = queue[front];
        int dist = distance[v];
        for (int k = row_ptr[v]; k < row_ptr[v + 1]; k++)
        {
            int net = col_ind[k];
            for (int j = row_ptr_inv[net]; j < row_ptr_inv[net + 1]; j++)
            {
                int pin = col_ind_inv[j];
                if (distance[pin] == -1)
                {
                    distance[pin] = dist + 1;
                    queue[++rear] = pin;
                }
            }
        }
        front++;
    }

    //Determine the farthest node and apply bfs on that node one more time
    int farthest = start_node; //we are sure that the node itself is reachable by distance 0
    int farthest_distance = 0;
    bool unreachable = false; //if some distance is -1 it means that node is unreachable
    for (int i = 0; i < nov; i++)
    {
        if (distance[i] > farthest_distance)
        {
            farthest = i;
            farthest_distance = distance[i];
        }
        else if (distance[i] == -1)
        {
            unreachable = true;
        }
    }
    if (unreachable)
    {
        cout << "There are unreachable nodes." << endl;
    }
    return farthest;
}

void bfsOrdering_queue(int start_node, unsigned int *row_ptr, unsigned int *col_ind, unsigned int *row_ptr_inv, unsigned int *col_ind_inv, int nov, int *distance, int *ordering)
{
    ordering[0] = start_node;
    int ordering_cnt = 1;

    int *queue = new int[nov];
    int front = 0, rear = 0;
    //Apply BFS
    for (int i = 0; i < nov; i++)
    {
        distance[i] = -1; // all nodes are unreachable in the beginning
        queue[i] = -1;
    }
    distance[start_node] = 0; //starting node is reachable from itself
    queue[front] = start_node;
    while (front <= rear)
    {
        int v = queue[front];
        int dist = distance[v];
        for (int k = row_ptr[v]; k < row_ptr[v + 1]; k++)
        {
            int net = col_ind[k];
            for (int j = row_ptr_inv[net]; j < row_ptr_inv[net + 1]; j++)
            {
                int pin = col_ind_inv[j];
                if (distance[pin] == -1)
                {
                    distance[pin] = dist + 1;
                    queue[++rear] = pin;
                    ordering[ordering_cnt++] = pin;
                }
            }
        }
        front++;
    }

    if (ordering_cnt != nov)
    {
        cout << "There are unreachable nodes." << endl;
        for (int i = 0; i < nov; i++)
        {
            if (distance[i] == -1)
                ordering[ordering_cnt++] = i;
        }
    }
}

void bfsOrderingSeperated(int start_node, unsigned int *row_ptr, unsigned int *col_ind, unsigned int *row_ptr_inv, unsigned int *col_ind_inv, int nov, int noe, int *distance, int *src_nodes, int *ordering_src, int *ordering_targ)
{
    int *src_status = new int[nov];
    for (int i = 0; i < nov; i++)
    {
        src_status[i] = -1;
    }

    int cnt_src = 0;
    for (int i = 0; i < noe; i++)
    {
        if (src_status[src_nodes[i]] == -1)
            cnt_src++;
        src_status[src_nodes[i]] = 1;
    }
    cout << "CNT SRCX" << cnt_src << endl;
    int ordering_src_cnt = 0;
    int ordering_targ_cnt = 0;
    if (src_status[start_node] == 1)
    {
        ordering_src[0] = start_node;
        ordering_src_cnt++;
    }
    else
    {
        ordering_targ[0] = start_node;
        ordering_targ_cnt++;
    }

    //Apply BFS
    for (int i = 0; i < nov; i++)
        distance[i] = -1;     // all nodes are unreachable in the beginning
    distance[start_node] = 0; //starting node is reachable from itself

    int step = 0;        //number of steps = breadth level
    bool changed = true; //if we did not change anything in distance array we can stop
    while (changed)
    {
        //cout << "girdik" << endl;
        changed = false;
        //int to_be_ordered_size = 0; //number of elements that are needed to be ordered at each level

        for (int i = 0; i < nov; i++)
        {
            if (distance[i] == step)
            {
                for (int k = row_ptr[i]; k < row_ptr[i + 1]; k++)
                {
                    int net = col_ind[k];
                    for (int j = row_ptr_inv[net]; j < row_ptr_inv[net + 1]; j++)
                    {
                        int pin = col_ind_inv[j];
                        if (distance[pin] == -1)
                        {
                            if (src_status[pin] == 1)
                                ordering_src[ordering_src_cnt++] = pin;
                            else
                                ordering_targ[ordering_targ_cnt++] = pin;

                            distance[pin] = step + 1;
                            changed = true;
                            //to_be_ordered [to_be_ordered_size++] = pin; //elements at each level that are needed to be ordered
                        }
                    }
                }
            }
        }
        step++;
    }

    cout << ordering_src_cnt << "--" << ordering_targ_cnt << endl;
    if (ordering_src_cnt + ordering_targ_cnt != nov)
    {
        cout << "There are unreachable nodes." << endl;
        for (int i = 0; i < nov; i++)
        {
            if (distance[i] == -1)
            {
                if (src_status[i] == 1)
                    ordering_src[ordering_src_cnt++] = i;
                else
                    ordering_targ[ordering_targ_cnt++] = i;
            }
        }
        cout << ordering_src_cnt << "--" << ordering_targ_cnt << endl;
        if (ordering_src_cnt + ordering_targ_cnt != nov)
        {
            cout << "STILL THERE ARE UNR" << endl;
        }
    }
}

int bfsFarthest_queue_par(int start_node, unsigned int *row_ptr, unsigned int *col_ind, unsigned int *row_ptr_inv, unsigned int *col_ind_inv, int nov, int *distance, int *glob_queue, int NT)
{
    omp_set_dynamic(0);
    omp_set_num_threads(NT);

    int *prefixSum = new int[NT + 1];
    int globalLen = 1;
    int *local_lengths = new int[NT];
    int **localQueuesList = new int *[NT];
    int level = 0;
    bool improvement = true;
    int farthest = -1, farthest_distance = -1;

    for (int i = 0; i < nov; i++)
    {
        distance[i] = glob_queue[i] = -1;
    }
    for (int i = 0; i < NT; i++)
    {
        localQueuesList[i] = new int[nov];
    }
    distance[start_node] = 0;
    glob_queue[0] = start_node;
    prefixSum[0] = 0;

    do
    {
        improvement = topDown(row_ptr, col_ind, row_ptr_inv, col_ind_inv, distance, level, glob_queue, globalLen, prefixSum, localQueuesList, local_lengths);
    } while (improvement);

    for (int i = 0; i < NT; i++)
    {
        delete[] localQueuesList[i];
    }

    for (int i = 0; i < nov; i++)
    {
        if (distance[i] > farthest_distance)
        {
            farthest = i;
            farthest_distance = distance[i];
        }
    }

    delete[] localQueuesList;
    delete[] prefixSum;
    return farthest;
}

void bfsOrdering_queue_par(int start_node, unsigned int *row_ptr, unsigned int *col_ind, unsigned int *row_ptr_inv, unsigned int *col_ind_inv, int nov, int *distance, int *glob_queue, int *ordering, int NT)
{
    omp_set_dynamic(0);
    omp_set_num_threads(NT);
    int *v_ordered_check = new int[nov];
    int *prefixSum = new int[NT + 1];
    int ordering_cnt = 1;
    int globalLen = 1;
    int *local_lengths = new int[NT];
    int **localQueuesList = new int *[NT];
    int level = 0;
    bool improvement = true;

    for (int i = 0; i < nov; i++)
    {
        distance[i] = glob_queue[i] = ordering[i] = -1;
        v_ordered_check[i] = 0;
    }
    for (int i = 0; i < NT; i++)
    {
        localQueuesList[i] = new int[nov];
    }

    ordering[0] = start_node;
    v_ordered_check[start_node] = 1;
    distance[start_node] = 0;
    glob_queue[0] = start_node;
    prefixSum[0] = 0;

    do
    {
        improvement = topDown(row_ptr, col_ind, row_ptr_inv, col_ind_inv, distance, level, glob_queue, globalLen, prefixSum, localQueuesList, local_lengths);
        for (int t = 0; t < NT; t++)
        {
            int loc_len = local_lengths[t];
            for (int k = 0; k < loc_len; k++)
            {
                int v = localQueuesList[t][k];
                if (!v_ordered_check[v])
                {
                    ordering[ordering_cnt++] = v;
                    v_ordered_check[v] = 1;
                }
            }
        }
    } while (improvement);

    for (int i = 0; i < NT; i++)
    {
        delete[] localQueuesList[i];
    }

    delete[] localQueuesList;
    delete[] prefixSum;
}

bool topDown(unsigned int *row, unsigned int *col, unsigned int *row_inv, unsigned int *col_inv, int *distance, int &level, int *globalQueue, int &globalLen, int *prefixSum, int **localQueuesList, int *local_lengths)
{
    bool improvement = false;

#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int *localQueue = localQueuesList[tid];
        int localLen = 0;

#pragma omp for reduction(|| \
                          : improvement) schedule(guided, 32)
        for (int i = 0; i < globalLen; i++)
        {
            int v = globalQueue[i];
            for (int j = row[v]; j < row[v + 1]; j++)
            {
                int e = col[j];
                for (int k = row_inv[e]; k < row_inv[e + 1]; k++)
                {
                    int u = col_inv[k];
                    if (distance[u] < 0)
                    {
                        distance[u] = level + 1;
                        localQueue[localLen++] = u;
                        improvement = true;
                    }
                }
            }
        }

        local_lengths[tid] = localLen;
        prefixSum[tid + 1] = localLen;
#pragma omp barrier

#pragma omp single
        {
            for (int i = 0; i < omp_get_num_threads(); i++)
            {
                prefixSum[i + 1] += prefixSum[i];
            }
        }

        memcpy(globalQueue + prefixSum[tid], localQueue, sizeof(int) * (prefixSum[tid + 1] - prefixSum[tid]));
        globalLen = prefixSum[omp_get_num_threads()];
    }

    if (improvement)
    {
        level++;
    }
    return improvement;
}

void writePartVector(string f_name, int nov, int K, int *part_vec, string algo_name)
{
    string folder = "/home/fatih/new_tensors";
    string only_hg = f_name.substr(f_name.rfind("/") + 1);
    string path = folder + only_hg + "_fastpart" + "." + to_string(K) + "." + algo_name;
    cout << "Path to part vector :  " << path << endl;
    ofstream out(path.c_str());
    cout << "Writing part vector to file.\n";
    if (out.fail())
        cout << "cant open the file";

    for (int i = 0; i < nov; i++)
    {
        out << part_vec[i] << "\n";
    }
    cout << "Writing complete.\n";
}

bool partitionCheck(int K)
{
    if (K <= 1)
    {
        cout << "You gave an incorrect partition input. You need to give a partition size bigger than 1." << endl;
        return false;
    }
    return true;
}

bool orderingCheck(int ORDERING_TYPE)
{
    if (ORDERING_TYPE != 1 && ORDERING_TYPE != 0)
    {
        cout << "You gave an incorrect ordering input. You need to give 1 as parameter for ordering and 0 for not ordering." << endl;
        return false;
    }
    return true;
}
