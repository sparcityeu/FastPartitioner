#include <iostream>
#include <string>
#include <fstream>
#include <sstream>



int main(int argc, char** argv)
{
  std::string fname = argv[1];
  std::ifstream infile(fname);

  int max_vertex = 0;
  int line_vertex = 0;
  int nnz = 0;
  int no_nets = 0;

  /*
  while (!infile.eof()) {
  std::cout << "Current max: " << max_vertex << std::endl;
    std::cout << "Current nnz: " << nnz << std::endl;
    infile >> line_vertex;
    std::cout << "Line vertex: " << line_vertex << std::endl;
    nnz++;

    if(line_vertex > max_vertex)
      max_vertex = line_vertex;
  }
  */
 
  std::string line;
  while (std::getline(infile, line))
    {
      //std::cout << "Current max: " << max_vertex << std::endl;
      //std::cout << "Current nnz: " << nnz << std::endl;
      //std::cout << "Line vertex: " << line_vertex << std::endl;
      no_nets += 1;
      std::istringstream iss(line);

      while(iss >> line_vertex){
	nnz++;

	if(line_vertex > max_vertex)
	  max_vertex = line_vertex;
      }
    }
  
  //std::cout << "Max vertex: " << max_vertex << std::endl;
  //std::cout << "NNz: " << nnz-2 << std::endl;
  //This will print header that PaToH wants
  //-1 and -2 is for overcounted first line
  std::cout << "1 " << max_vertex << " " << no_nets-1 << " " << nnz-2;  
}
