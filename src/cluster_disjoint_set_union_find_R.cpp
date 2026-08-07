#include <Rcpp.h>
#include <vector>
#include <numeric>
#include <unordered_map>

using namespace Rcpp;

// creates the disjoint sets from a list of pairs to be clustered
class DisjointSet {
  public:
    std::vector<int> parent;
    std::vector<int> rank;
  
  DisjointSet(int n) {
    // Reserve memory all at once to prevent reallocations
    parent.resize(n);
    rank.resize(n, 0);
    std::iota(parent.begin(), parent.end(), 0);
  }
  
  // Inline find operation with path compression for raw performance
  inline int find(int i) {
    int root = i;
    while (root != parent[root]) {
      root = parent[root];
    }
    // Path compression loop (iterative to avoid recursion stack overhead)
    int curr = i;
    while (curr != root) {
      int nxt = parent[curr];
      parent[curr] = root;
      curr = nxt;
    }
    return root;
  }
  
  // Inline union operation by rank
  inline void unite(int i, int j) {
    int root_i = find(i);
    int root_j = find(j);
    if (root_i != root_j) {
      if (rank[root_i] < rank[root_j]) {
        parent[root_i] = root_j;
      } else if (rank[root_i] > rank[root_j]) {
        parent[root_j] = root_i;
      } else {
        parent[root_j] = root_i;
        rank[root_i]++;
      }
    }
  }
};


//' Cluster Node Pairs into Disjoint Sets
//' 
//' @param total_nodes Total number of nodes (0-indexed integers).
//' @param from_nodes std::vector<int> of starting nodes for pairs.
//' @param to_nodes std::vector<int> of ending nodes for pairs.
//' @return A std::vector<int> where each element is a disjoint cluster and contains the node IDs of that cluster.
//' @export

// [[Rcpp::export]]
std::vector<std::vector<int>> get_clusters_from_pairs(int total_nodes, 
                                                     std::vector<int> from_nodes, 
                                                     std::vector<int> to_nodes) 
{
  DisjointSet dsu(total_nodes);
  
  // 1. Process edges quickly using raw pointers for speed
  int num_edges = from_nodes.size();
  
  if (num_edges < to_nodes.size()) {
    Rcpp::stop("The provided list of node pairs do not have equal length. Missing nodes in the to_nodes list."
                 " Aborting clustering pairs to disjoint sets.");
  }
  
  for (int i = 0; i < num_edges; ++i) {
    if (from_nodes[i] >= total_nodes || to_nodes[i] >= total_nodes) {
      Rcpp::stop("The provided list of node pairs contain a node number greater "
                   "than the total number of nodes - outside the provided valid range. "
                   " Aborting clustering pairs to disjoint sets.");
    }
    dsu.unite(from_nodes[i], to_nodes[i]);
  }
  
  // 2. Linear time grouping
  int unique_clusters_count = 0;
  
  // Map the actual DSU root to a contiguous cluster index (0, 1, 2...)
  std::vector<int> root_to_cluster_idx(total_nodes, -1);
  // extract the final list of clusters
  std::vector<std::vector<int>> out_list(total_nodes);
  
  for (int i = 0; i < total_nodes; ++i) {
    int root = dsu.find(i);
    if (root_to_cluster_idx[root] == -1) {
      root_to_cluster_idx[root] = unique_clusters_count++;
    }
    out_list[root_to_cluster_idx[root]].push_back(i);
  }
  
  out_list.resize(unique_clusters_count);
  
  return out_list;
}


// Rcpp::sourceCpp(file.path("src/disjoint_set_union_find_R.cpp"))
// total_nodes=24
// from = rep(c(0,1,2,3),5)
// to = seq(4, 23)
// from = rep(seq(0, 49),50)
// to = seq(0, 2499)
// total_nodes=2500