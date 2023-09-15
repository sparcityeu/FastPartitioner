#include "hypergraph_rw.h"

using namespace std;

bool thinning(string f_name, unsigned int ** row_ptr_p, unsigned int ** col_ind_p,unsigned int ** vertex_appearence, unsigned int &nov, unsigned int &noe, long long &pins, int min){

  
  cout << "Thinning func " << endl;
  int MIN = min;
  
//GRAPH VARIABLE INITIALIZATIONS
  int index_base;
    
  auto start = chrono::steady_clock::now();
      
  //Read From Original .hg file
  ifstream mean_inp(f_name);
  ifstream hg_inp(f_name);
  hg_inp >> index_base >> nov >> noe >> pins;
  cout << "Index base = " << index_base << " |V| = " << nov << " |N| = " << noe << " |pins| = " << pins << endl;
  (*row_ptr_p) = new unsigned int[noe + 1];
  (*row_ptr_p)[0] = 0;
  (*col_ind_p) = new unsigned int[pins];
  (*vertex_appearence) = new unsigned int[nov];
  initialize(*vertex_appearence, 0, nov, (unsigned int)0);
  string line;
  int v_count = 0;
  int n_count = 1;

  int* degrees = new int[noe];
  
  
  int tot_degree = 0;
  int en = 0;

  
  getline(mean_inp, line);
  while(getline(mean_inp, line)){
  istringstream iss(line);
      unsigned int vertex;
      int v_count_local = 0;
      
      while(iss >> vertex){
	tot_degree++;
	v_count_local++;
      }
      degrees[en] = v_count_local;
      en++;
  }
  

  double mean = tot_degree / noe;
  
  cout << "Mean: " << mean << endl;
  
  en = -1;
  int skipped = 0;
  int new_pin = 0;
  getline(hg_inp, line);//read the empty remaining part
  while(getline(hg_inp, line))
    {
      en++;
      cout << "degree: " << degrees[en] << " en: " << en << " line: " << line <<endl;
      if(degrees[en] < 50*mean && degrees[en] < MIN){
	cout << "degree: " << degrees[en] << " en: " << en << endl;
	istringstream iss(line);
	unsigned int vertex;
	int v_count_local = 0;
	while(iss >> vertex)
	{
	  new_pin++;
	  if(vertex == nov && index_base != 1)
	    cout <<"vertex error " << vertex << endl;
	  else if(vertex > nov)
	    cout <<"vertex error " << vertex << endl;
	  (*col_ind_p)[v_count++] = vertex;
	  (*vertex_appearence)[vertex]++;
	  v_count_local++;
	}
      (*row_ptr_p)[n_count] = (*row_ptr_p)[n_count - 1] + v_count_local;
      n_count++;
      }
      else{
	skipped++;
      }
      noe -= skipped;
      pins = new_pin;

      return true;
    }
auto end = chrono::steady_clock::now();
cout << "Reading hypergraph took: " << chrono::duration <double, milli> (end - start).count()/1000 << " seconds."<< endl;
}



bool read_patohlike(string f_name, unsigned int ** col_ptr_p, unsigned int ** row_ind_p,unsigned int ** vertex_appearence, unsigned int &nov, unsigned int &noe, long long &pins)
{
  cout << "Reading func " << endl;
  
  //GRAPH VARIABLE INITIALIZATIONS
  string bin_name = f_name + ".plbin";
  int index_base;
  ifstream bin_inp(bin_name, ios::binary);
  if(bin_inp.fail())
    {
      cout << "No binary file detected." << endl;
      auto start = chrono::steady_clock::now();
      
      //Read From Original .hg file
      ifstream hg_inp(f_name);
      hg_inp >> index_base >> nov >> noe >> pins;
      cout << "Index base = " << index_base << " |V| = " << nov << " |N| = " << noe << " |pins| = " << pins << endl;
      (*col_ptr_p) = new unsigned int[noe + 1];
      (*col_ptr_p)[0] = 0;
      (*row_ind_p) = new unsigned int[pins];
      (*vertex_appearence) = new unsigned int[nov];
      initialize(*vertex_appearence, 0, nov, (unsigned int)0);
      string line;
      int v_count = 0;
      int n_count = 1;

      int edge_salt = 0;

      
      if(index_base == 0){
	edge_salt = 0;
      }
      else{
	edge_salt = -1;
      }
      
      
      if(!hg_inp.fail())
        {
	  getline(hg_inp, line);//read the empty remaining part
	  while(getline(hg_inp, line) && n_count <= noe)
            {
	      istringstream iss(line);
	      bool empty_line = true;
	      unsigned int vertex;
	      int v_count_local = 0;
	      while(iss >> vertex)
                {
		  empty_line = false;
		  if(vertex == nov && index_base != 1)
		    cout <<"vertex error " << vertex << endl;
		  else if(vertex > nov)
		    cout <<"vertex error " << vertex << endl;
		  (*row_ind_p)[v_count++] = vertex + edge_salt;
		  (*vertex_appearence)[vertex + edge_salt]++;
		  v_count_local++;
                }
	      if(empty_line == false){
		(*col_ptr_p)[n_count] = (*col_ptr_p)[n_count - 1] + v_count_local;
		n_count++;
	      }
            }
	  
	  if(n_count - 1 != noe)
            {
	      cout << "Number of nets in the hypergraph does not match with the hyper parameters." << endl;
	      cout << "Hyperparameter in the file :" << noe << " != " << n_count << " Number of nets in the file."<< endl;
	      cout << "Setting noe to n count " << endl;
	      noe = n_count - 1;
            }
	  if(v_count != pins)
            {
	      cout << "Number of nonzeros in the hypergraph does not match with the hyper parameters." << endl;
	      cout << "Hyperparameter in the file :" << pins << " != " << v_count << " Number of nonzeros in the file."<< endl;
	      cout << "Exiting..." << endl;
	      return false;
            }
        }
      else
        {
	  cout << "File" << f_name << " does not exist." << endl;
	  cout << "Exiting." << endl;
	  return false;
        }
      auto end = chrono::steady_clock::now();
      cout << "Reading hypergraph took: " << chrono::duration <double, milli> (end - start).count()/1000 << " seconds."<< endl;
      start = chrono::steady_clock::now();
      
      //Write to binary file
      cout << "Creating a binary file with name : " << bin_name << endl;
      ofstream bin_out(bin_name, ios :: binary);
      bin_out.write((char *)&noe, sizeof(unsigned int)); bin_out.write((char *)&nov, sizeof(unsigned int));bin_out.write((char *)&pins, sizeof(long long));
      for(int i = 0; i < noe; i++)
        {
	  unsigned int n_reach = (*col_ptr_p)[i+1];
	  bin_out.write((char *)&n_reach, sizeof(unsigned int));
	  for(int k = (*col_ptr_p)[i]; k < (*col_ptr_p)[i + 1]; k++)
            {
	      unsigned int vertex = (*row_ind_p)[k];
	      bin_out.write((char *)&vertex, sizeof(unsigned int));
            }
        }
      bin_out.flush();
      bin_out.close();
      end = chrono::steady_clock::now();
      cout << "Creating binary file took: " << chrono::duration <double, milli> (end - start).count()/1000 << " seconds."<< endl;
    }
  else
    {
      cout << "Binary file detected." << endl;
      cout << "Reading hypergraph from file :" << bin_name << endl;
      auto start = chrono::steady_clock::now();
      
      bin_inp.read((char *)&noe, sizeof(unsigned int));
      bin_inp.read((char *)&nov, sizeof(unsigned int));
      bin_inp.read((char *)&pins, sizeof(long long));
      cout << "|N| = " << noe << " |V| = " << nov << " |pins| = " << pins << endl;
      
      (*col_ptr_p) = new unsigned int[noe + 1];
      (*row_ind_p) = new unsigned int[pins];
      (*vertex_appearence) = new unsigned int[nov];//
      initialize(*vertex_appearence, 0, nov, (unsigned int)0);
      int v_count = 0;
      int n_count = 1;
      unsigned int prev_col_ptr_val = 0,col_ptr_val, row_ind_val;
      (*col_ptr_p)[0] = 0;
      while(!bin_inp.eof() && n_count <= noe)
        {
	  bin_inp.read((char*)&col_ptr_val, sizeof(unsigned int));
	  (*col_ptr_p)[n_count] = col_ptr_val;
	  int v_cnt_local = 0;
	  while(v_cnt_local < col_ptr_val - prev_col_ptr_val)
            {
	      bin_inp.read((char*)&row_ind_val, sizeof(unsigned int));
	      (*row_ind_p)[v_count++] = row_ind_val;//
	      (*vertex_appearence)[row_ind_val]++;//
	      v_cnt_local++;
            }
	  n_count++;
	  prev_col_ptr_val = col_ptr_val;
        }
      auto end = chrono::steady_clock::now();
      cout << "Reading binary file took: " << chrono::duration <double, milli> (end - start).count()/1000 << " seconds."<< endl;
    }
  return true;
}


bool read_hg(string f_name, unsigned int ** row_ptr_p, unsigned int ** col_ind_p,unsigned int ** net_appearence, unsigned int &nov, unsigned int &noe, long long &pins)
{
    //GRAPH VARIABLE INITIALIZATIONS
    string bin_name = f_name + ".bin";

    ifstream bin_inp(bin_name, ios::binary);
    if(bin_inp.fail())
    {
        cout << "No binary file detected." << endl;
        auto start = chrono::steady_clock::now();

        //Read From Original .hg file
        ifstream hg_inp(f_name);
        hg_inp >> nov >> noe >> pins;
        cout << "|V| = " << nov << " |N| = " << noe << " |pins| = " << pins << endl;
        (*row_ptr_p) = new unsigned int[nov + 1];
        (*row_ptr_p)[0] = 0;
        (*col_ind_p) = new unsigned int[pins];
        (*net_appearence) = new unsigned int[noe];
        initialize(*net_appearence, 0, noe, (unsigned int)0);
        string line;
        int v_count = 1;
        int n_count = 0;
        if(!hg_inp.fail())
        {
            getline(hg_inp, line);//read the empty remaining part
            while(getline(hg_inp, line) && v_count <= nov)
            {
                istringstream iss(line);
                unsigned int net;
                int n_count_local = 0;
                while(iss >> net)
                {
                    if(net >= noe)
                    cout <<"net error " << endl;
                    (*col_ind_p)[n_count++] = net;
                    (*net_appearence)[net]++;
                    n_count_local++;
                }
                (*row_ptr_p)[v_count] = (*row_ptr_p)[v_count - 1] + n_count_local;
                v_count++;
            }
            if(v_count - 1 != nov)
            {
                cout << "Number of vertices in the hypergraph does not match with the hyper parameters." << endl;
                cout << "Hyperparameter in the file :" << nov << " != " << v_count << " Number of vertices in the file."<< endl;
                cout << "Exiting..." << endl;
                return false;
            }
            if(n_count != pins)
            {
                cout << "Number of nonzeros in the hypergraph does not match with the hyper parameters." << endl;
                cout << "Hyperparameter in the file :" << pins << " != " << n_count << " Number of nonzeros in the file."<< endl;
                cout << "Exiting..." << endl;
                return false;
            }
        }
        else
        {
            cout << "File" << f_name << " does not exist." << endl;
            cout << "Exiting." << endl;
            return false;
        }
        auto end = chrono::steady_clock::now();
        cout << "Reading hypergraph took: " << chrono::duration <double, milli> (end - start).count()/1000 << " seconds."<< endl;
        start = chrono::steady_clock::now();

        //Write to binary file
        cout << "Creating a binary file with name : " << bin_name << endl;
        ofstream bin_out(bin_name, ios :: binary);
        bin_out.write((char *)&nov, sizeof(unsigned int)); bin_out.write((char *)&noe, sizeof(unsigned int));bin_out.write((char *)&pins, sizeof(long long));
        for(int i = 0; i < nov; i++)
        {
            unsigned int v_deg = (*row_ptr_p)[i+1];
            bin_out.write((char *)&v_deg, sizeof(unsigned int));
            for(int k = (*row_ptr_p)[i]; k < (*row_ptr_p)[i + 1]; k++)
            {
                unsigned int net = (*col_ind_p)[k];
                bin_out.write((char *)&net, sizeof(unsigned int));
            }
        }
        bin_out.flush();
        bin_out.close();
        end = chrono::steady_clock::now();
        cout << "Creating binary file took: " << chrono::duration <double, milli> (end - start).count()/1000 << " seconds."<< endl;
    }
    else
    {
        cout << "Binary file detected." << endl;
        cout << "Reading hypergraph from file :" << bin_name << endl;
        auto start = chrono::steady_clock::now();

        bin_inp.read((char *)&nov, sizeof(unsigned int));
        bin_inp.read((char *)&noe, sizeof(unsigned int));
        bin_inp.read((char *)&pins, sizeof(long long));
        cout << "|V| = " << nov << " |N| = " << noe << " |pins| = " << pins << endl;

        (*row_ptr_p) = new unsigned int[nov + 1];
        (*col_ind_p) = new unsigned int[pins];
        (*net_appearence) = new unsigned int[noe];
        initialize(*net_appearence, 0, noe, (unsigned int)0);
        int v_count = 1;
        int n_count = 0;
        unsigned int prev_row_ptr_val = 0,row_ptr_val, col_ind_val;
        (*row_ptr_p)[0] = 0;
        while(!bin_inp.eof() && v_count <= nov)
        {
            bin_inp.read((char*)&row_ptr_val, sizeof(unsigned int));
            (*row_ptr_p)[v_count] = row_ptr_val;
            int n_cnt_local = 0;
            while(n_cnt_local < row_ptr_val - prev_row_ptr_val)
            {
                bin_inp.read((char*)&col_ind_val, sizeof(unsigned int));
                (*col_ind_p)[n_count++] = col_ind_val;
                (*net_appearence)[col_ind_val]++;
                n_cnt_local++;
            }
            v_count++;
            prev_row_ptr_val = row_ptr_val;
        }
        auto end = chrono::steady_clock::now();
        cout << "Reading binary file took: " << chrono::duration <double, milli> (end - start).count()/1000 << " seconds."<< endl;
    }
    return true;
}

void assignWeights_patohlike(string file_name, unsigned int * row_ptr, unsigned int * col_ind, unsigned int * vertex_appearence, unsigned int nov, unsigned int noe, long long pins, int ** edge_wgths, int ** src_nodes)
{
  unsigned int *row_ptr_inv;// cumulative net counts from 0'th vertex to n-1'th vertex
  unsigned int *col_ind_inv;// nonzero values

  std::cout << "Assigning weights-> noe: " << noe << std::endl;
  
  *edge_wgths = new int[noe];
  *src_nodes= new int[noe];
  createInverse_patohlike(row_ptr, col_ind, vertex_appearence, nov, noe, pins, &row_ptr_inv, &col_ind_inv);
  //initialize result arays
  for(int i = 0; i < noe; i++)
    {
    (*edge_wgths)[i] = -1;
    (*src_nodes)[i] = -1;
    }
  
  int dotLoc = file_name.rfind(".");
  string baseName = file_name.substr(0, dotLoc);
  string outName = baseName + ".ws";
  ifstream fCheck(outName.c_str());
  if(!fCheck.good())
    {
      ofstream out(outName);
      //Determine source nodes
      int edge_cnt = 0;
      for(int i = 0; i <  noe; i++)
	{
	  int deg = row_ptr[i + 1] - row_ptr[i];
	  if(deg == 0)
	    deg = 1;
	  int src_ind = rand() % deg;
	  int src = col_ind[row_ptr[i] + src_ind];
	  (*src_nodes)[i] = src;
	}
      
      for(int i = 0; i < noe; i++)
	{
	  int wght = rand() % 10;
	  (*edge_wgths)[i] = wght;
	  out << wght << " " << (*src_nodes)[i] << endl;
	}
      out.close();
    }
  else
    {
      cout << "Weight source file exists." << endl;
      string line;
      int edge_cnt = 0;
      while(getline(fCheck, line))
	{
	  istringstream iss(line);
	  int wght, src;
	  iss >> wght >> src;
	  (*edge_wgths)[edge_cnt] = wght;
	  (*src_nodes)[edge_cnt] = src;
	  edge_cnt++;
	}
    }
  cout << "Weights of the edges and their source nodes have been assigned." << endl;
}


void assignWeights(string file_name, unsigned int * row_ptr, unsigned int * col_ind, unsigned int * net_appearence, unsigned int nov, unsigned int noe, long long pins, int ** edge_wgths, int ** src_nodes)
{
  unsigned int *row_ptr_inv;// cumulative net counts from 0'th vertex to n-1'th vertex
  unsigned int *col_ind_inv;// nonzero values

  std::cout << "Assigning weights-> noe: " << noe << std::endl;
  
  *edge_wgths = new int[noe];
  *src_nodes= new int[noe];
  createInverse(row_ptr, col_ind, net_appearence, nov, noe, pins, &row_ptr_inv, &col_ind_inv);
  //initialize result arays
  for(int i = 0; i < noe; i++)
    {
    (*edge_wgths)[i] = -1;
    (*src_nodes)[i] = -1;
    }
  
  int dotLoc = file_name.rfind(".");
  string baseName = file_name.substr(0, dotLoc);
  string outName = baseName + ".ws";
  ifstream fCheck(outName.c_str());
  if(!fCheck.good())
    {
      ofstream out(outName);
      //Determine source nodes
      int edge_cnt = 0;
      for(int i = 0; i <  noe; i++)
	{
	  int deg = row_ptr_inv[i + 1] - row_ptr_inv[i];
	  if(deg == 0)
	    deg = 1;
	  int src_ind = rand() % deg;
	  int src = col_ind_inv[row_ptr_inv[i] + src_ind];
	  (*src_nodes)[i] = src;
	}
      
      for(int i = 0; i < noe; i++)
	{
	  int wght = rand() % 10;
	  (*edge_wgths)[i] = wght;
	  out << wght << " " << (*src_nodes)[i] << endl;
	}
      out.close();
    }
  else
    {
      cout << "Weight source file exists." << endl;
      string line;
      int edge_cnt = 0;
      while(getline(fCheck, line))
	{
	  istringstream iss(line);
	  int wght, src;
	  iss >> wght >> src;
	  (*edge_wgths)[edge_cnt] = wght;
	  (*src_nodes)[edge_cnt] = src;
	  edge_cnt++;
	}
    }
  cout << "Weights of the edges and their source nodes have been assigned." << endl;
}
