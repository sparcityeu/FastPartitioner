#include "metric.h"
#include "utils.h"
int CN(int * ep_ends, int noe)
{
	int CN = 0;
	for(int i = 0; i < noe; i++)
	{
		if(ep_ends[i] > 1)
		{
			CN++;
		}
	}
	return CN;
}

int CN_Refined(int * ep_parts, int * ep_cnt, int * ep_starts, int noe)
{
	int CN = 0;
	for(int e = 0; e < noe; e++)
	{
		int distinct = 0;
		for(int j = ep_starts[e]; j < ep_starts[e + 1]; j++)
		{
			if(ep_parts[j] >= 0 && ep_cnt[j] > 0){
				distinct++;
			}
		}
		if(distinct > 1)
			CN++;
	}
	return CN;
}
int CN_Scratch(int * part_vec, unsigned int * row_ptr, unsigned int * col_ind, unsigned int * net_appearence, int K, int nov, int noe, int pins)
{
    unsigned int *row_ptr_inv;    // cumulative net counts from 0'th vertex to n-1'th vertex
    unsigned int *col_ind_inv;    // nonzero values
    createInverse(row_ptr, col_ind, net_appearence, nov, noe, pins, &row_ptr_inv, &col_ind_inv);
    int * distinct = new int [K];
    int CN = 0;
    for(int i = 0; i < noe; i++){
        initialize(distinct, 0, K, 0);
        int distinct_cnt = 0;
        for(int j = row_ptr_inv[i]; j < row_ptr_inv[i + 1]; j++){
            int v = col_ind_inv[j];
            int part_id = part_vec[v];
            if(distinct[part_id] == 0){
                distinct_cnt++;
                distinct[part_id] = 1;
            }
        }
        if(distinct_cnt > 1){
            CN += 1;
        }
    }
    return CN;
}
int TV(int * ep_ends, int noe)
{
	int TV = 0;
	for(int i = 0; i < noe; i++)
	{
		if(ep_ends[i] != 0)
		{
			TV += ep_ends[i] - 1;
		}
	}
	return TV;
}

int TV_Refined(int * ep_starts, int * ep_cnts, int noe, int K)
{
	int TV = 0;
	for(int e = 0; e < noe; e++)
	{
		int distinct = 0;
		for(int j = ep_starts[e]; j < ep_starts[e + 1]; j++)
		{
			if(ep_cnts[j] > 0){
				distinct++;
			}
		}
		if(distinct > 0)
			TV += distinct - 1;
	}
	return TV;
}

int TV_Scratch(int * part_vec, unsigned int * row_ptr, unsigned int * col_ind, unsigned int * net_appearence, int K, int nov, int noe, int pins)
{
    unsigned int *row_ptr_inv;    // cumulative net counts from 0'th vertex to n-1'th vertex
    unsigned int *col_ind_inv;    // nonzero values
    createInverse(row_ptr, col_ind, net_appearence, nov, noe, pins, &row_ptr_inv, &col_ind_inv);
    int * distinct = new int [K];
    int TV = 0;
    for(int i = 0; i < noe; i++){
        initialize(distinct, 0, K, 0);
        int distinct_cnt = 0;
        for(int j = row_ptr_inv[i]; j < row_ptr_inv[i + 1]; j++){
            int v = col_ind_inv[j];
            int part_id = part_vec[v];
            if(distinct[part_id] == 0){
                distinct_cnt++;
                distinct[part_id] = 1;
            }
        }
        if(distinct_cnt > 1){
            TV += distinct_cnt - 1;
        }
    }
    return TV;
}

int TM(int * part_messages, int K)
{
	int TM = 0;
	for(int i = 0; i < K; i++)
	{
		TM += part_messages[i];
	}
	return TM;
}
int TM_scratch(int K, int * labels,  int * ep_parts, int * ep_starts, int * ep_ends, int nov, int noe, int * src_nodes)
{
	//calculate tm from partition itself
	int ** part_connections = new int*[K];
	int * part_messages = new int[K];
	int TM = 0;
	vector<vector<int> > part_vec(K);
	vector<vector<int> > vertex_src_status(nov);

	int distinct_src_node_cnt = 0;
	for(int i = 0; i < noe; i++)
	  {
	    int src = src_nodes[i];
	    
	    if(src < 0 || src >= nov)
	      cout << "src: " << src << " " << nov << endl;
		
	    if(vertex_src_status[src_nodes[i]].size() == 0)
			distinct_src_node_cnt++;
		vertex_src_status[src_nodes[i]].push_back(i);
	}
	for(int i = 0; i < K; i++)
	{
		part_connections[i] = new int [K];
		part_messages[i] = 0;
		for(int j = 0; j < K; j++)
		{
			if(i != j)
				part_connections[i][j] = -1;                 // no connections in the begining
			else
				part_connections[i][j] = 1;                 // all parts are connected with itself
		}
	}
	for(int i = 0; i < nov; i++)
	{
		part_vec[labels[i]].push_back(i);
	}
	for(int k = 0; k < K; k++)
	{
		for(int v_ind = 0; v_ind < part_vec[k].size(); v_ind++)
		{
			int vertex = part_vec[k][v_ind];
			for(int edge_ind = 0; edge_ind < vertex_src_status[vertex].size(); edge_ind++)
			{
				int edge = vertex_src_status[vertex][edge_ind];
				int end_cnt = ep_ends[edge];
				int start_index = ep_starts[edge];
				for(int j = 0; j < end_cnt; j++)
				{
					int part_id = ep_parts[start_index + j];
					if(k != part_id && part_connections[k][part_id] == -1)
					{
						part_connections[k][part_id] = 1;
						part_messages[k]++;
					}
				}
			}
		}
	}
	for(int i = 0; i < K; i++)
	{
		TM += part_messages[i];
	}
	return TM;
}

int TM_Refined(int ** part_con_cnt,int K){
	int TM = 0;
	for(int i = 0; i < K; i++){
		for(int j = 0; i !=j && j < K; j++){
			if(part_con_cnt[i][j] > 1){
				TM++;
			}
			if(part_con_cnt[i][j] < 1){
				cout << "no way" << endl;
			}
		}
	}
	return TM;
}
int MGM(int * part_messages, int K)
{
	int MGM = 0;
	//cout << endl ;
	for(int i = 0; i < K; i++)
	{
		//cout << part_messages[i] << " - ";
		if(part_messages[i] > MGM)
			MGM = part_messages[i];
	}
	//cout << endl;
	return MGM;
}

int MGM_scratch(int K, int * labels,  int * ep_parts, int * ep_starts, int * ep_ends, int nov,int noe, int * src_nodes)
{
	int ** part_connections = new int*[K];
	int * part_messages = new int[K];
	int TM = 0;
	vector<vector<int> > part_vec(K);
	vector<vector<int> > vertex_src_status(nov);

	int distinct_src_node_cnt = 0;
	for(int i = 0; i < noe; i++)
	{
		if(vertex_src_status[src_nodes[i]].size() == 0)
			distinct_src_node_cnt++;
		vertex_src_status[src_nodes[i]].push_back(i);
	}

	for(int i = 0; i < K; i++)
	{
		part_connections[i] = new int [K];
		part_messages[i] = 0;
		for(int j = 0; j < K; j++)
		{
			if(i != j)
				part_connections[i][j] = -1;                 // no connections in the begining
			else
				part_connections[i][j] = 1;                 // all parts are connected with itself
		}
	}

	for(int i = 0; i < nov; i++)
	{
		part_vec[labels[i]].push_back(i);
	}

	for(int k = 0; k < K; k++)
	{
		for(int v_ind = 0; v_ind < part_vec[k].size(); v_ind++)
		{
			int vertex = part_vec[k][v_ind];
			for(int edge_ind = 0; edge_ind < vertex_src_status[vertex].size(); edge_ind++)
			{
				int edge = vertex_src_status[vertex][edge_ind];
				int end_cnt = ep_ends[edge];
				int start_index = ep_starts[edge];
				for(int j = 0; j < end_cnt; j++)
				{
					int part_id = ep_parts[start_index + j];
					if(k != part_id && part_connections[k][part_id] == -1)
					{
						part_connections[k][part_id] = 1;
						part_messages[k]++;
					}
				}
			}
		}
	}

	int max = -1;
	int max_part = -1;
	for(int i = 0; i < K; i++)
	{
		if(part_messages[i] > max)
		{
			max_part = i;
			max = part_messages[i];
		}
	}
	return max;
}

int MGV(int * send_volumes, int K)
{
	int max = -1;
	int max_part = -1;
	for(int i = 0; i < K; i++)
	{
		if(send_volumes[i] > max)
		{
			max_part = i;
			max = send_volumes[i];
		}
	}
	return max;
}

int MGV_Scratch(int * part_vec, unsigned int * row_ptr, unsigned int * col_ind, int * src_nodes, int * ep_ends,
	 unsigned int * net_appearence, int K, int nov, int noe, int pins)
{

    unsigned int *row_ptr_inv;    // cumulative net counts from 0'th vertex to n-1'th vertex
    unsigned int *col_ind_inv;    // nonzero values
    createInverse(row_ptr, col_ind, net_appearence, nov, noe, pins, &row_ptr_inv, &col_ind_inv);
	int * send_volumes = new int [K];
	initialize(send_volumes, 0, K, 0);
    for(int i = 0; i < noe; i++){
        for(int j = row_ptr_inv[i]; j < row_ptr_inv[i + 1]; j++){
            int v = col_ind_inv[j];
			if(src_nodes[i] == v)
			{
				int sender_part = part_vec[v];
				send_volumes[sender_part] += ep_ends[i] - 1;
			}
        }
    }

	int max = -1;
	int max_part = -1;
	int sum = 0;
	for(int i = 0; i < K; i++)
	{
		if(send_volumes[i] > max)
		{
			max_part = i;
			max = send_volumes[i];
		}
		sum += send_volumes[i];
	}

	return max;
}

int MGAV(int * send, int* received, int K)
{

	int max = -1;
	int max_part = -1;
	for(int i = 0; i < K; i++)
	{
		if(send[i] + received[i] > max)
		{
			max_part = i;
			max = send[i] + received[i];
		}
	}
	return max;
}

int MGAV_Scratch(int * part_vec, unsigned int * row_ptr, unsigned int * col_ind, int * src_nodes, int * ep_ends, int * ep_parts, int * ep_starts,
	 unsigned int * net_appearence, int K, int nov, int noe, int pins)
{
	unsigned int *row_ptr_inv;    // cumulative net counts from 0'th vertex to n-1'th vertex
	unsigned int *col_ind_inv;    // nonzero values
	createInverse(row_ptr, col_ind, net_appearence, nov, noe, pins, &row_ptr_inv, &col_ind_inv);
	int * send = new int [K];
	int * received = new int [K];
	initialize(received, 0, K, 0);
	initialize(send, 0, K, 0);

	for(int i = 0; i < noe; i++){
	  int src_node_part = part_vec[src_nodes[i]];
		send[src_node_part] += ep_ends[i] - 1;
		int end_cnt = ep_ends[i];
		int start_index = ep_starts[i];
		for(int j = 0; j < end_cnt; j++)
		{
			int part_id = ep_parts[start_index + j];
			if(part_id != src_node_part)
			{
				received[part_id]++;
			}
		}
	}

	int max = -1;
	int max_part = -1;
	for(int i = 0; i < K; i++)
	{
		if(send[i] + received[i] > max)
		{
			max_part = i;
			max = send[i] + received[i];
		}
	}

	return max;
}
