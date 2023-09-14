# FastPartitioner

Table of Contents
----------------
* [Hypergraph Definition](#hypergraph)
* [Hypergraph Partitioning](#partitioning)
* [Partitioning Metrics](#metrics)
   * [Cut Net](#met_cn)
   * [Total Volume(k-1, connectivity - 1)](#met_tv)
   * [Total Send Message](#met_tm)
   * [Maximum Send Message](#met_msm)
   * [Maximum Send Volume](#met_msv)
   * [Maximum Send Received Volume](#met_msrv)
* [Requirements Before Build](#reqs)
* [How to Build FastPartitioner](#build)


<a name="hypergraph"> <font size="6">Hypergraph Definition</font></a>
----------------
A *hypergraph* **H = (V, N)** is defined as a set of vertices **V** and a set of nets
(hyperedges) **N** among those vertices. A net n ∈ N is a subset of vertices and
the vertices in n are called its **pins**.

![Hypergraph](readme_fig/HPGRAPH.png)

<div>
    <p> <br /><br /> Figure 1: A hypergraph representation with 4 nets (edges) 
    <b>(N<sub>1</sub>, 
    N<sub>2</sub>, 
    N<sub>3</sub>, 
    N<sub>4</sub>)</b>
    and 8 vertices 
    <b>(v<sub>1</sub>, v<sub>2</sub>, v<sub>3</sub>, v<sub>4</sub>, v<sub>5</sub>,v<sub>6</sub>,v<sub>7</sub>,v<sub>8</sub>)</b>. 
    Pin sets of the nets are given below.
    <br />
      N<sub>1</sub> = {v<sub>1</sub>, v<sub>2</sub>, v<sub>3</sub>}
      <br/>
      N<sub>2</sub> = {v<sub>1</sub>, v<sub>6</sub>,
      v<sub>7</sub>, v<sub>8</sub>} <br />
      N<sub>3</sub>= {v<sub>3</sub>, v<sub>4</sub>, v<sub>5</sub>}<br />
      N<sub>4</sub>= {v<sub>6</sub>,v<sub>7</sub>}<br />
      <br /><br />
  </p>
</div>

<a name="partitioning"> <font  size="6">Hypergraph partitioning</font></a>
----------------
A K-way partition of a hypergraph H is a partition of its vertex set, which
is denoted as &prod; = {V<sub>1</sub>, V<sub>2</sub>, . . . , V<sub>k</sub>}, where <br />
• parts are pairwise disjoint, i.e., V<sub>k</sub> &cap; V<sub>l</sub> = &empty; for all 1 &le; k &le; l &le; K<br />
• each part V<sub>k</sub> is a nonempty subset of V, i.e., V<sub>k</sub> &sube; V and V<sub>k</sub> &ne; &empty; for 1  &le;  k &le;  K,<br />
• the union of K parts is equal to V, <br />

Let W<sub>k</sub> denote the total vertex weight in V<sub>k</sub>, that is W<sub>k</sub> = &sum; <sub>v&isin;V<sub>k</sub></sub> w[v] and W<sub>avg</sub> denote the weight of each part when the total vertex weight is equally distributed, that is W<sub>avg</sub> = &sum; <sub>v&isin;V<sub>k</sub></sub> w[v] / K. If each part V<sub>k</sub> &isin; Π satisfies the balance criterion  <br />
W<sub>k</sub> &le; W<sub>avg</sub>(1 + ε), for k = 1, 2, . . . , K <br />
We say that Π is ε-balanced where ε is called the maximum allowed imbalance
ratio.


![Hypergraph_partition](readme_fig/HPGRAPH_partitioned.png)
<div>
    <p> <br /><br /><br /><br />Figure 2: A partitioning result of a hypergraph with 4 nets <b>(N<sub>1</sub>, N<sub>2</sub>, N<sub>3</sub>, N<sub>4</sub>)</b>  and 8 vertices <b>(v<sub>1</sub>, v<sub>2</sub>, .., v<sub>8</sub>)</b> into <b>4</b> parts. Resulting partitions are given below.<br />
     P<sub>1</sub> = {v<sub>1</sub>, v<sub>2</sub>}<br />
     P<sub>2</sub> = {v<sub>3</sub>, v<sub>4</sub>, v<sub>5</sub>}<br />
     P<sub>3</sub>= {v<sub>6</sub>, v<sub>7</sub>, v<sub>8</sub>}<br />
     <br />
     <br /><br /><br /><br /> <br />
  </p>
</div>


  </p>
</div>

<a name="metrics"> <font size="6">**Partitioning Metrics** </font></a>
----------------
Quality of the hypergraph partitioninig problem is evaluated based on different metrics. These metrics are generally based on the connectivity of the hypernets after the partitioning. One can find detailed information about the metrics that we are interested in this project below.

<a name="met_cn"> <font size="5">**Cut Net** </font></a>

Cut net is the total number of hypernets that have connectivity values more than 1(Basically connected to more than 1 partition).

<a name="met_tv"> <font size="5">**Total Volume(k-1, connectivity - 1)** </font></a>

Total volume is the most popular metric in the literature. It's calculated over a sum of connectivity values of the edges with substracting 1 from them. It's also called k-1 or connectivity - 1 for that reason.

<a name="met_tm"> <font size="5">**Total Message** </font></a>

Total message is a metric that based on  possible communication requirements between the partitions. Each net has a special vertex in it's vertex set called the <b>source node</b> and all the other nodes of it called <b>target nodes</b>. We say a <b>message</b> is sent from partition of a source node of a net to it's partitions of it's target nodes. Then total message is the total number of different communication pairs among the partitions.

<a name="met_msm"> <font size="5">**Maximum Send Message** </font></a>

Maximum send message is the number of messages that is send by the part which sends the most messages. Algorithms for this metric focuses on removing the bottleneck which caused by some processor by being responsible for most of the communication between the parts.

<a name="met_msv"> <font size="5">**Maximum Send Message** </font></a>

Maximum send volume is the cost that is caused by the part which has the source nodes of nets that causes the most connections between the parts when they summed . Algorithms for this metric focuses on removing the bottleneck which caused by some processor by being responsible for most of the communication over the nets.

<a name="met_msrv"> <font size="5">**Maximum Send Message** </font></a>

Maximum send received volume is the cost that is caused by the part which has the source nodes of nets that causes the most connections between the parts when they summed and in addition to that it's target nodes which causes communication with other parts. Algorithms for this metric focuses on removing the bottleneck which caused by some processor by being responsible for most of the communication over the nets.

<a name="build"> <font size="7">**Build** </font></a>
FastPartitioner itself can partition and reorder hypergraphs and tensors. Moreover, it can automate experimenting process leveraging other tools. Those tools are Splatt and PaToH. You need to install those tools and their dependencies to fully employ FastPartitioner capabilities. Please note that **FastPartitioner assumes user have a CUDA capable GPU**.

- Install BLAS and LAPACK             
`sudo dnf install blas-devel`       
`sudo dnf install lapack-devel`


- Clone Splatt under `FastPartitioner/`     
https://github.com/ShadenSmith/splatt

- Download binary distribution of PaToH, unzip and put `patoh` executable under `FastPartitioner/`
https://faculty.cc.gatech.edu/~umit/software.html       



