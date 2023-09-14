extern "C"
{
    #include "graphio.h"
    #include "graph.h"
}

#include<fstream>
#include<iostream>
#include<string>

char gfile[2048];

using namespace std;

int main(int argc, char *argv[])
{
    //GRAPH VARIABLE INITIALIZATIONS
    etype *row_ptr, *row_ptr_inv; // cumulative net counts from 0'th vertex to n-1'th vertex
    vtype *col_ind, *col_ind_inv; // nonzero values
    ewtype *ewghts; // edge weights
    vwtype *vwghts; // vertex weights
    vtype nov, noe; // number of vertices in the hyper graph
    int * edge_wgths;
    int * src_nodes;
    char * fname  = argv[1];
    bool symmetric = false;
    string file_name = string(fname);

    //GRAPH READ
    string inv_name = string(fname) + ".inv";
    strcpy(gfile, fname);
    if(read_graph(gfile, &row_ptr, &col_ind, &ewghts, &vwghts, &nov, 0) == -1)
    {
        printf("error in graph read\n");
        exit(1);
    }

    string hg_name = file_name.substr(0, file_name.find(".mtx")) + ".hg";
    cout << "Creating file with name : " << hg_name << endl;
    ofstream out(hg_name);
    int * diag_check = new int[nov];
    for(int i = 0; i < nov; i++) diag_check[i] = 0;

    out << nov << " " << nov << " " << row_ptr[nov] << endl;
    for(int i = 0; i < nov; i++)
    {
        for(int k = row_ptr[i]; k < row_ptr[i+1]; k++)
        {
            int net = col_ind[k];
            out << net << " ";
        }
        out << endl;
    }
    out.flush();
    out.close();
    return 0;
}
