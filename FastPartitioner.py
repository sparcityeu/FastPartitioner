import sys
import subprocess
import argparse


def inform(return_code, stdout, stderr):

    if(return_code == 0):
        print("Success!")
        print("Output:")
        print(stdout.decode())
    else:
        print("Command failed with return code:", return_code)
        print("Output:")
        print(stdout.decode())
        print("Error output:")
        print(stderr.decode())

def run_command(command):

    command_str = ""
    for item in command:
        command_str += item
        command_str += " "
    print("command:", command_str)
    
    # Run the command and wait for it to complete
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()  # Wait for the process to finish and capture its output
    
    return_code = process.returncode
    inform(return_code, stdout, stderr)

    return stdout

        
def convert_mtx_to_hg(mtx_file):

    print("FastPartitioner::convert_mtx_to_hg")
    
    
    command = ["./mtx_to_vertex_list", mtx_file]
    hg_file = mtx_file.replace(".mtx", ".hg")

    print("out file:", hg_file)
    
    stdout = run_command(command)
    
    ##Make header patoh compatible
    command = ["./patohread_help", hg_file]
    stdout = run_command(command)
    new_header = stdout.decode()
    
    reader = open(hg_file, "r")
    lines = reader.readlines()
    reader.close()
    lines[0] = new_header + "\n"

    writer = open(hg_file, "w")
    for line in lines:
        writer.write(line)
    
    return hg_file

'''
def convert_tns_to_hg():


def reorder_tensor():


def cpd():


def partition_with_patoh():


def partition_with_splatt():


def partition_with_fastpart():


def partition_with_streaming():
'''


def get_args_for_fastpart(args):

    if(args.metric == None):
        print("Target metric must be stated for FastPartitioner, exiting..")
        exit(1)
    if(args.partitions == None):
        print("Number of partitions must be stated for FastPartitioner, exiting..")
        exit(1)
        
    exe = "./" + args.metric + ".exe"
        
    return [exe, "-f", args.hypergraph_file, "-k", args.partitions, "-o", args.ordering, "-r", args.refinement,
            "-i", args.iterations, "-t", args.threads, "-g", args.gpu, "-m", args.sparse, "-n", args.sparse_min, "-z", args.sparse_mean]



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="FastPartitioner, automated partition and reorder tool for matrix and tensors.")
    parser.add_argument("function", help="Functionality to execute")
    parser.add_argument("-im", "--matrix-file", help = "Input matrix file")
    parser.add_argument("-it", "--tensor-file", help = "Input tensor file")
    parser.add_argument("-ih", "--hypergraph-file", help = "Input hypergraph file")
    parser.add_argument("-m", "--metric", help = "FastPartitioner target metric")
    parser.add_argument("-k", "--partitions", help = "Number of partitions")
    parser.add_argument("-o", "--ordering", default = "0", help = "Ordering algorithm for partitioning for FastPartitioner, 0 for no ordering, 1 for BFS ordering, default is 0")
    parser.add_argument("-r", "--refinement", default = "0", help = "Refinement option for FastPartitioner, default is 0")
    parser.add_argument("-i", "--iterations", default = "1", help = "Iterations of refinement for FastPartitioner, default is 1")
    parser.add_argument("-t", "--threads", default = "1", help = "Number of threads for FastPartitioner, default is 1")
    parser.add_argument("-g", "--gpu", default = "0", help = "Use GPU refinement algorithm for FastPartitioner, default is 0")
    parser.add_argument("-s", "--sparse", default = "0", help = "sparsification on hypergraphs for FastPartitioner, default is 0")
    parser.add_argument("-smin", "--sparse-min", default = "0", help = "sparsification on on hypergraphs for FastPartitioner using minimum degree, default is 0")
    parser.add_argument("-smean", "--sparse-mean", default = "0", help = "sparsification on on hypergraphs for FastPartitioner using mean degree, default is 0")
    ##Parse the command-line arguments
    args = parser.parse_args()
    function = args.function
    

    def switch_case(function, args):

        if(function == "mtx_to_hg"):
            print("FastPartitioner::mtx_to_hg")
            mtx_file = args.matrix_file
            hg_file = convert_mtx_to_hg(mtx_file)

        if(function == "fastpart"):
            print("FastPartitioner::FastPartitioner")
            command = get_args_for_fastpart(args)
            run_command(command)

    switch_case(function, args)
