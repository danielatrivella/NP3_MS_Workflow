#include <Rcpp.h>
#include <vector>
#include <numeric>
#include <unordered_map>

using namespace Rcpp;

// creates the disjoint sets from a list of pairs to be clustered
class DisjointSet {
  std::vector<int> parent;
  std::vector<int> rank;
  
  public:
    DisjointSet(int n) {
      parent.resize(n);
      rank.resize(n, 0);
      std::iota(parent.begin(), parent.end(), 0);
    }
  
  int find(int i) {
    if (parent[i] == i)
      return i;
    return parent[i] = find(parent[i]); // Path compression
  }
  
  void unite(int i, int j) {
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
std::vector<std::vector<int>> get_clusters_from_pairs(int total_nodes, std::vector<int> from_nodes, std::vector<int> to_nodes) {
  DisjointSet dsu(total_nodes);
  
  // Process edges/pairs
  int num_edges = from_nodes.size();
  for (int i = 0; i < num_edges; i++) {
    if (from_nodes[i] >= total_nodes || to_nodes[i] >= total_nodes) {
      Rcpp::stop("The provided list of node pairs contain a node number greater "
        "than the total number of nodes - outside the provided valid range. "
        " Aborting clustering pairs to disjoint sets.");
    }
    dsu.unite(from_nodes[i], to_nodes[i]);
  }
  
  // Group nodes by their root representative
  std::unordered_map<int, std::vector<int>> clusters_map;
  for (int i = 0; i < total_nodes; i++) {
    int root = dsu.find(i);
    clusters_map[root].push_back(i);
  }
  
  // extract the final list of clusters
  std::vector<std::vector<int>> out_list;

  for (auto it = clusters_map.begin(); it != clusters_map.end(); it++) {
    int root = it->first;
    const std::vector<int>& nodes = it->second;

    std::vector<int> cluster_nodes(nodes.begin(), nodes.end());
    out_list.push_back(cluster_nodes);
  }
  
  return out_list;
}
