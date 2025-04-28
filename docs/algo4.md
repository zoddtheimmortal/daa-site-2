# CoreEXACT Algorithm

Yixiang Fang†⋆, Kaiqiang Yu‡, Reynold Cheng‡, Laks V.S. Lakshmanan§, Xuemin Lin†⋆

The following C++ implementations are based on [Efficient Algorithms for Densest Subgraph Discovery](https://www.vldb.org/pvldb/vol12/p1719-fang.pdf).

## Pseudo Code

```bash
Algorithm 4: The algorithm: CoreExact.

Input: G(V, E), Ψ(V0, E0);
Output: The CDS D(V0, E0);

perform core decomposition using Algorithm 3;
locate the (k*, Ψ)-core using pruning criteria;
initialize C ← ∅, D ← ∅, U ← ∅, t ← ℓ*, ρ' ← ℓ*, kmax;
put all the connected components of (k*, Ψ)-core into C;
for each connected component C'(V0C, E0C) ∈ C do
  if |V0C| ≠ then C(V0C, E0C) ← C ∩ (V0C, E0C), Ψ-core;
  build a flow network F(VF, EF, EF') by lines 5-15 of Algorithm 1;
  find minimum s-t cut (S, T) from F(VF, EF, EF');
  if S≠∅ then continue;
  while t - |Σ| ≥ tmin do
    α ← |Σ|;
    build F(VF, EF, EF') by lines 5-15 of Algorithm 1;
    find minimum s-t cut (S, T) from F(VF, EF, EF');
    if S≠∅ then
      t ← α;
    else if α > |T| then remove some vertices from C;
      t ← α;
      U ← S\{s};
    if G[U], Ψ) > ρ(D, U) then D ← G[U];
 return D;

```

## Usage

Save any of the following code as `algo1.cpp`.

Run the following code in a terminal.

```bash
g++ -O3 algo4.cpp
./a.out
```

We use the `-O3 flag` as it enables aggressive optimizations that can significantly improve the performance of the program

Output for the code will be saved in `path_for_output_file`.

Terminal will display errors, debugging and progress statements.

## Implementation

The below implementation is a one-to-one implementation of the research paper. No additional optimizations have been added.

### Code

```cpp
#include <iostream>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <chrono>
#include <unordered_map>
#include <queue>
#include <climits>
#include <cmath>
#include <string>

using namespace std;
using namespace std::chrono;

using Graph = vector<vector<int>>;
using Flow = vector<vector<pair<int,int>>>;

#include <vector>
#include <algorithm>

typedef vector<vector<int>> Graph;
vector<int> core_decomp;
int n,m;

void findAllKCliques(const Graph& G, int k, vector<int>& currentClique, int start, vector<vector<int>>& result) {
    if (currentClique.size() == static_cast<size_t>(k)) {
        result.push_back(currentClique);
        return;
    }

    for (int v = start; v < G.size(); v++) {
        bool isConnectedToAll = true;
        for (int u : currentClique) {
            bool edgeExists = false;
            for (int neighbor : G[u]) {
                if (neighbor == v) {
                    edgeExists = true;
                    break;
                }
            }

            if (!edgeExists) {
                isConnectedToAll = false;
                break;
            }
        }

        if (isConnectedToAll) {
            currentClique.push_back(v);
            findAllKCliques(G, k, currentClique, v + 1, result);
            currentClique.pop_back();
        }
    }
}

vector<vector<int>> findAllKCliquesWrapper(const Graph& G, int k) {
    vector<vector<int>> result;
    vector<int> currentClique;
    findAllKCliques(G, k, currentClique, 0, result);
    return result;
}

Graph createGraph(int n, const vector<pair<int, int>>& edges) {
    Graph G(n);
    for (const auto& edge : edges) {
        G[edge.first].push_back(edge.second);
        G[edge.second].push_back(edge.first);
    }
    for (auto& neighbors : G) {
        sort(neighbors.begin(), neighbors.end());
    }
    return G;
}

void addEdge(Flow& G, int u, int v, int c) {
    G[u].push_back({v, c});
}

void printFlowNetwork(const Flow& network, int n, int num_lambda, ofstream& outfile) {
    outfile << "\nFlow Network:" << endl;
    outfile << "s (0) -> Vertices (1 to " << n << "):" << endl;
    for (const auto& edge : network[0]) {
        if (edge.first >= 1 && edge.first <= n) {
            outfile << "  s -> " << (edge.first - 1) << " (capacity: " << edge.second << ")" << endl;
        }
    }

    int t = network.size() - 1;
    outfile << "Vertices -> t (" << t << "):" << endl;
    for (int i = 1; i <= n; i++) {
        for (const auto& edge : network[i]) {
            if (edge.first == t) {
                outfile << "  " << (i - 1) << " -> t (capacity: " << edge.second << ")" << endl;
            }
        }
    }

    outfile << "Vertices -> (k-1)-cliques:" << endl;
    for (int i = 1; i <= n; i++) {
        for (const auto& edge : network[i]) {
            if (edge.first > n && edge.first < t) {
                outfile << "  " << (i - 1) << " -> lambda_" << (edge.first - n - 1) << " (capacity: " << edge.second << ")" << endl;
            }
        }
    }

    outfile << "(k-1)-cliques -> Vertices:" << endl;
    for (int i = n + 1; i < t; i++) {
        for (const auto& edge : network[i]) {
            if (edge.first >= 1 && edge.first <= n) {
                outfile << "  lambda_" << (i - n - 1) << " -> " << (edge.first - 1) << " (capacity: " << edge.second << ")" << endl;
            }
        }
    }
}

vector<vector<int>> k_cliques;

vector<vector<int>> CLIQUES(const Graph& G, int k) {
    k_cliques = findAllKCliquesWrapper(G, k);
    return k_cliques;
}

void printCliquesCount(const vector<vector<int>>& cliques, ofstream& outfile) {
    unordered_map<int, int> size_count;

    for (const auto& clique : cliques) {
        size_count[clique.size()]++;
    }

    vector<pair<int, int>> sorted_sizes(size_count.begin(), size_count.end());
    sort(sorted_sizes.begin(), sorted_sizes.end());

    for (const auto& pair : sorted_sizes) {
        outfile << "Size " << pair.first << ": " << pair.second << endl;
    }
}

vector<vector<int>> printKCliques(const vector<vector<int>>& k_cliques, ofstream& outfile, int k) {
    vector<vector<int>> k_size_cliques;
    outfile << "\nAll " << k << "-cliques (" << k_cliques.size() << " in total):\n";
    for (const auto& clique : k_cliques) {
        vector<int> temp;
        // outfile << "[";
        for (size_t i = 0; i < clique.size(); ++i) {
            temp.push_back(clique[i]);
            // outfile << clique[i];
            // if (i != clique.size() - 1) {
            //     outfile << ", ";
            // }
        }
        k_size_cliques.push_back(temp);
        // outfile << "]\n";
    }
    return k_size_cliques;
}




bool readGraphFromFile(const string& filename, int& n, vector<pair<int, int>>& edges) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return false;
    }

    string line;
    int m;

    if (getline(file, line)) {
        istringstream iss(line);
        if (!(iss >> n >> m)) {
            cerr << "Error: Invalid format for vertices and edges count" << endl;
            return false;
        }
    } else {
        cerr << "Error: Empty file" << endl;
        return false;
    }

    vector<pair<int, int>> rawEdges;
    unordered_set<int> uniqueVertices;

    while (getline(file, line)) {
        istringstream iss(line);
        int u, v;
        if (iss >> u >> v) {
            rawEdges.emplace_back(u, v);
            uniqueVertices.insert(u);
            uniqueVertices.insert(v);
        } else {
            cerr << "Warning: Invalid edge format in line: " << line << endl;
        }
    }

    if (uniqueVertices.size() != static_cast<size_t>(n)) {
        cerr << "Warning: Number of unique vertices (" << uniqueVertices.size()
             << ") doesn't match the specified count (" << n << ")" << endl;
        n = uniqueVertices.size();
    }

    unordered_map<int, int> vertexToIndex;
    int index = 0;
    for (int v : uniqueVertices) {
        vertexToIndex[v] = index++;
    }

    edges.clear();
    for (const auto& edge : rawEdges) {
        int u = vertexToIndex[edge.first];
        int v = vertexToIndex[edge.second];
        edges.emplace_back(u, v);
    }

    file.close();
    return true;
}

void find_clique_degrees(const vector<vector<int>>& kSizeCliques, vector<int>& clique_degree) {
    for (const auto& clique : kSizeCliques) {
        for (int v : clique) {
            clique_degree[v]++;
        }
    }
}

int find_maximum_clique_degree(const vector<int>& clique_degree) {
    int max_degree = 0;
    for (int degree : clique_degree) {
        if (degree > max_degree) {
            max_degree = degree;
        }
    }
    return max_degree;
}

bool bfs(const Flow& graph, int s, int t, vector<int>& parent) {
    vector<bool> visited(graph.size(), false);
    queue<int> q;
    q.push(s);
    visited[s] = true;
    parent[s] = -1;

    while (!q.empty()) {
        int u = q.front();
        q.pop();

        for (const auto& edge : graph[u]) {
            int v = edge.first;
            int capacity = edge.second;
            if (!visited[v] && capacity > 0) {
                q.push(v);
                parent[v] = u;
                visited[v] = true;
                if (v == t) return true;
            }
        }
    }
    return false;
}

int fordFulkerson(Flow& graph, int s, int t, vector<bool>& minCut, ofstream& outfile) {
    int V = graph.size();
    Flow residualGraph = graph;

    for (int u = 0; u < V; u++) {
        for (const auto& edge : graph[u]) {
            int v = edge.first;
            bool reverseExists = false;
            for (const auto& revEdge : graph[v]) {
                if (revEdge.first == u) {
                    reverseExists = true;
                    break;
                }
            }
            if (!reverseExists) {
                residualGraph[v].push_back({u, 0});
            }
        }
    }

    vector<int> parent(V);
    int maxFlow = 0;

    while (bfs(residualGraph, s, t, parent)) {
        int pathFlow = INT_MAX;
        for (int v = t; v != s; v = parent[v]) {
            int u = parent[v];
            for (const auto& edge : residualGraph[u]) {
                if (edge.first == v) {
                    pathFlow = min(pathFlow, edge.second);
                    break;
                }
            }
        }

        for (int v = t; v != s; v = parent[v]) {
            int u = parent[v];

            for (auto& edge : residualGraph[u]) {
                if (edge.first == v) {
                    edge.second -= pathFlow;
                    break;
                }
            }

            for (auto& edge : residualGraph[v]) {
                if (edge.first == u) {
                    edge.second += pathFlow;
                    break;
                }
            }
        }
        maxFlow += pathFlow;
    }

    minCut.assign(V, false);
    queue<int> q;
    q.push(s);
    minCut[s] = true;

    while (!q.empty()) {
        int u = q.front();
        q.pop();

        for (const auto& edge : residualGraph[u]) {
            int v = edge.first;
            int capacity = edge.second;
            if (!minCut[v] && capacity > 0) {
                minCut[v] = true;
                q.push(v);
            }
        }
    }

    // outfile << "\nResidual Graph after Max Flow:" << endl;
    // for (int u = 0; u < V; u++) {
    //     if (!residualGraph[u].empty()) {
    //         outfile << u << " -> ";
    //         for (const auto& edge : residualGraph[u]) {
    //             outfile << "(" << edge.first << "," << edge.second << ") ";
    //         }
    //         outfile << endl;
    //     }
    // }

    return maxFlow;
}

bool isValidKClique(const Graph& G, const vector<int>& potentialClique) {
    for (size_t i = 0; i < potentialClique.size(); i++) {
        for (size_t j = i+1; j < potentialClique.size(); j++) {
            int u = potentialClique[i];
            int v = potentialClique[j];
            if (find(G[u].begin(), G[u].end(), v) == G[u].end()) {
                return false;
            }
        }
    }
    return true;
}

int findMaxCore() {
    return *max_element(core_decomp.begin(), core_decomp.end());
}

vector<int> getVerticesInCore(int k_core) {
    vector<int> core_vertices;
    for(int v = 0; v < core_decomp.size(); v++) {
        if(core_decomp[v] >= k_core) {
            core_vertices.push_back(v);
        }
    }
    return core_vertices;
}

void findDensestSubgraph(Graph& G, vector<vector<int>> kSizeCliques, int k, int n, ofstream& outfile) {
    double l = 0;
    vector<int> clique_degree(n, 0);
    find_clique_degrees(kSizeCliques, clique_degree);

    int k_max = findMaxCore();
    int psi_size = k;

    // Get initial core vertices
    int initial_core = ceil(l);
    vector<int> V_core = getVerticesInCore(initial_core);

    vector<int> V(G.size());
    for (int i = 0; i < G.size(); i++) {
        V[i] = i;
    }

    double u = (double)find_maximum_clique_degree(clique_degree);
    vector<vector<int>> LAMBDA = findAllKCliquesWrapper(G, k-1);
    vector<int> dense_vertices;

    outfile << "Initial parameters:" << endl;
    outfile << "l: " << l << ", u: " << u << endl;
    outfile << "Number of (k-1)-cliques: " << LAMBDA.size() << endl;
    outfile << endl;

    double epsilon = 1.0 / (n * (n - 1));
    while ((u - l) >= epsilon) {
        double alpha = (l + u) / 2.0;
        outfile << "\nCurrent alpha: " << alpha << endl;

        int s = 0;
        int t = 1 + n + LAMBDA.size();
        Flow Vf(t + 1);

        for (int v : V) {
            addEdge(Vf, s, v + 1, clique_degree[v]);
        }

        for (int v : V) {
            addEdge(Vf, v + 1, t, (alpha * k));
        }

        vector<vector<int>> completions(LAMBDA.size());
        for (int lambda_idx = 0; lambda_idx < LAMBDA.size(); lambda_idx++) {
            const auto& lambda = LAMBDA[lambda_idx];
            for (int v : V) {
                if (find(lambda.begin(), lambda.end(), v) != lambda.end()) {
                    continue;
                }

                vector<int> potential_clique = lambda;
                potential_clique.push_back(v);
                sort(potential_clique.begin(), potential_clique.end());

                if (isValidKClique(G, potential_clique)) {
                    completions[lambda_idx].push_back(v);
                }
            }
        }

        const int INF = 1000000;

        for (int lambda_idx = 0; lambda_idx < LAMBDA.size(); lambda_idx++) {
            int lambda_node = n + 1 + lambda_idx;
            for (int v : completions[lambda_idx]) {
                addEdge(Vf, v + 1, lambda_node, 1);
            }
        }

        for (int lambda_idx = 0; lambda_idx < LAMBDA.size(); lambda_idx++) {
            int lambda_node = n + 1 + lambda_idx;
            for (int v : LAMBDA[lambda_idx]) {
                addEdge(Vf, lambda_node, v + 1, INF);
            }
        }


        vector<bool> minCut;
        int maxFlow = fordFulkerson(Vf, s, t, minCut, outfile);
        outfile << "Max Flow: " << maxFlow << endl;

        outfile << "Min cut (nodes in S side): ";
        for (int i = 0; i < minCut.size(); i++) {
        }

        outfile << "Vertices in source side: ";
        bool any_vertex_in_cut = false;
        for (int v : V) {
            if (minCut[v + 1]) {
                any_vertex_in_cut = true;
            }
        }

        bool s_has_more = false;
        vector<int> current_dense;
        for (int v : V) {
            if (minCut[v + 1]) {  // Correct offset for core vertices
                current_dense.push_back(v);
                s_has_more = true;
            }
        }

        double new_l=l;
        double new_u=u;
        double new_best_density = 0;
        if (s_has_more) {
            new_l = alpha;
            double current_den = static_cast<double>(current_dense.size())
                               / (V.size() * (V.size() - 1));
            if (current_den > new_best_density) {
                new_best_density = current_den;
            }
        } else {
            new_u = alpha;
        }

        if (!any_vertex_in_cut) {
            outfile << "No vertices in source side, updating u to " << alpha << endl;
            u = alpha;
        } else {
            outfile << "Some vertices in source side, updating l to " << alpha << endl;
            l = alpha;

            dense_vertices.clear();
            for (int v : V) {
                if (minCut[v + 1]) {
                    dense_vertices.push_back(v);
                }
            }
        }

        outfile << "Current gap: " << (u - l) << endl;
    }

    Graph D(n);
    if (!dense_vertices.empty()) {
        for (int v : dense_vertices) {
            for (int u : G[v]) {
                if (find(dense_vertices.begin(), dense_vertices.end(), u) != dense_vertices.end()) {
                    D[v].push_back(u);
                }
            }
        }
    }

    outfile << "\nDense subgraph found with " << dense_vertices.size() << " vertices:" << endl;
    if (dense_vertices.empty()) {
        outfile << "No dense subgraph found." << endl;
    } else {
        outfile << "Vertices in dense subgraph: ";
        for (int v : dense_vertices) {
            outfile << v << " ";
        }
        outfile << endl;

        outfile << "Edges in dense subgraph:" << endl;
        for (int v : dense_vertices) {
            outfile << "Vertex " << v << " connected to: ";
            for (int u : D[v]) {
                outfile << u << " ";
            }
            outfile << endl;
        }

        int edge_count = 0;
        for (int v : dense_vertices) {
            edge_count += D[v].size();
        }
        edge_count /= 2;

        double density = 0;
        if (dense_vertices.size() > 1) {
            density = (2.0 * edge_count) / (dense_vertices.size() * (dense_vertices.size() - 1));
        }

        outfile << "Number of edges: " << edge_count << endl;

        int clique_count = 0;
        for (const auto& clique : kSizeCliques) {
            bool in_dense = true;
            for (int v : clique) {
                if (find(dense_vertices.begin(), dense_vertices.end(), v) == dense_vertices.end()) {
                    in_dense = false;
                    break;
                }
            }
            if (in_dense) {
                clique_count++;
            }
        }

        outfile << "Number of k-cliques in dense subgraph: " << clique_count << endl;
        if (dense_vertices.size() > 0) {
            double avg_cliques = (double)clique_count / dense_vertices.size();
            outfile << "Average k-cliques per vertex (Density): " << avg_cliques << endl;
        }
    }
}



vector<int> clique_instances(int v){
    vector<int> clique_indexes;
    for (int i = 0; i < k_cliques.size(); i++) {
        if (find(k_cliques[i].begin(), k_cliques[i].end(), v) != k_cliques[i].end()) {
            clique_indexes.push_back(i);
        }
    }
    return clique_indexes;
}

void performCoreDecomp(Graph& G) {
    vector<int> degree(G.size());
    find_clique_degrees(k_cliques, degree);
    vector<int> vertices(G.size());
    for (int i = 0; i < G.size(); i++) {
        vertices[i] = i;
    }

    sort(vertices.begin(), vertices.end(), [&degree](int a, int b) {
        return degree[a] < degree[b];
    });

    while(!vertices.empty()) {
        int v = vertices[0];
        core_decomp[v] = degree[v];

        vector<int> instances = clique_instances(v);
        for (int& x : instances) {
            for (int& y : k_cliques[x]) {
                if (y != v && degree[y] > degree[v]) {
                    degree[y]--;
                }
            }
        }

        vertices.erase(vertices.begin());

        sort(vertices.begin(), vertices.end(), [&degree](int a, int b) {
            return degree[a] < degree[b];
        });
    }

}


int main() {
    vector<string> input_files = {
        "1as_caida_processed.txt",
    };

    for (const string& input_file : input_files) {
        for (int k = 3; k<=3; k++) {
            n=m=0;
            string output_file = "final_output_" + input_file + "_" + to_string(k) + ".txt";

            vector<pair<int, int>> edges;

            if (!readGraphFromFile(input_file, n, edges)) {
                cerr << "Error reading graph from file: " << input_file << endl;
                continue;
            }

            ofstream outfile(output_file);
            if (!outfile.is_open()) {
                cerr << "Error: Could not open " << output_file << " for writing" << endl;
                continue;
            }

            outfile << "Graph with " << n << " vertices and " << edges.size() << " edges" << endl;
            outfile << "Looking for cliques of size k = " << k << endl;


            Graph G = createGraph(n, edges);

            auto start_time = high_resolution_clock::now();
            vector<vector<int>> cliques = CLIQUES(G, k);

            core_decomp.resize(n, 0);
            performCoreDecomp(G);

            printCliquesCount(cliques, outfile);
            vector<vector<int>> kSizeCliques = printKCliques(k_cliques, outfile, k);


            cout << "Found " << k_cliques.size() << " " << k << "-cliques in " << input_file << "." << endl;

            findDensestSubgraph(G, kSizeCliques, k, n, outfile);

            auto end_time = high_resolution_clock::now();
            outfile << "Execution time: " << duration_cast<milliseconds>(end_time - start_time).count() << " milliseconds" << endl;

            outfile.close();
            cout << "Results written to " << output_file << endl;
        }
    }

    return 0;
}
```

### Results

#### Netscience Dataset

For k=2

```bash
Graph with 1461 vertices and 2742 edges
Looking for cliques of size k = 2
Size 2: 2742

All 2-cliques (2742 in total):
Initial parameters:
l: 0, u: 34
Number of (k-1)-cliques: 1461


Current alpha: 17
Max Flow: 5484
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 17
Current gap: 17

Current alpha: 8.5
Max Flow: 5444
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 8.5
Current gap: 8.5

Current alpha: 12.75
Max Flow: 5484
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 12.75
Current gap: 4.25

Current alpha: 10.625
Max Flow: 5484
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 10.625
Current gap: 2.125

Current alpha: 9.5625
Max Flow: 5484
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 9.5625
Current gap: 1.0625

Current alpha: 9.03125
Max Flow: 5464
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 9.03125
Current gap: 0.53125

Current alpha: 9.29688
Max Flow: 5464
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 9.29688
Current gap: 0.265625

Current alpha: 9.42969
Max Flow: 5464
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 9.42969
Current gap: 0.132812

Current alpha: 9.49609
Max Flow: 5464
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 9.49609
Current gap: 0.0664062

Current alpha: 9.5293
Max Flow: 5484
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 9.5293
Current gap: 0.0332031

Current alpha: 9.5127
Max Flow: 5484
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 9.5127
Current gap: 0.0166016

Current alpha: 9.50439
Max Flow: 5484
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 9.50439
Current gap: 0.00830078

Current alpha: 9.50024
Max Flow: 5484
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 9.50024
Current gap: 0.00415039

Current alpha: 9.49817
Max Flow: 5464
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 9.49817
Current gap: 0.0020752

Current alpha: 9.49921
Max Flow: 5464
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 9.49921
Current gap: 0.0010376

Current alpha: 9.49973
Max Flow: 5464
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 9.49973
Current gap: 0.000518799

Current alpha: 9.49998
Max Flow: 5464
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 9.49998
Current gap: 0.000259399

Current alpha: 9.50011
Max Flow: 5484
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 9.50011
Current gap: 0.0001297

Current alpha: 9.50005
Max Flow: 5484
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 9.50005
Current gap: 6.48499e-05

Current alpha: 9.50002
Max Flow: 5484
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 9.50002
Current gap: 3.24249e-05

Current alpha: 9.5
Max Flow: 5484
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 9.5
Current gap: 1.62125e-05

Current alpha: 9.49999
Max Flow: 5464
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 9.49999
Current gap: 8.10623e-06

Current alpha: 9.5
Max Flow: 5464
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 9.5
Current gap: 4.05312e-06

Current alpha: 9.5
Max Flow: 5464
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 9.5
Current gap: 2.02656e-06

Current alpha: 9.5
Max Flow: 5464
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 9.5
Current gap: 1.01328e-06

Current alpha: 9.5
Max Flow: 5484
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 9.5
Current gap: 5.06639e-07

Current alpha: 9.5
Max Flow: 5484
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 9.5
Current gap: 2.5332e-07

Dense subgraph found with 20 vertices:
Vertices in dense subgraph: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Edges in dense subgraph:
Vertex 527 connected to: 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 529 connected to: 527 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 823 connected to: 527 529 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 1030 connected to: 527 529 823 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 1032 connected to: 527 529 823 1030 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 1034 connected to: 527 529 823 1030 1032 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 1036 connected to: 527 529 823 1030 1032 1034 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 1038 connected to: 527 529 823 1030 1032 1034 1036 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 1040 connected to: 527 529 823 1030 1032 1034 1036 1038 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 1042 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 1044 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 1046 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 1051 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1053 1055 1057 1059 1061 1068 1287
Vertex 1053 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1055 1057 1059 1061 1068 1287
Vertex 1055 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1057 1059 1061 1068 1287
Vertex 1057 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1059 1061 1068 1287
Vertex 1059 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1061 1068 1287
Vertex 1061 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1068 1287
Vertex 1068 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1287
Vertex 1287 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068
Number of edges: 190
Number of k-cliques in dense subgraph: 190
Average k-cliques per vertex: 9.5
Execution time: 1958 milliseconds


```

For k=3

```bash
Graph with 1461 vertices and 2742 edges
Looking for cliques of size k = 3
Size 3: 3764

All 3-cliques (3764 in total):
Initial parameters:
l: 0, u: 173
Number of (k-1)-cliques: 2742


Current alpha: 86.5
Max Flow: 11292
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 86.5
Current gap: 86.5

Current alpha: 43.25
Max Flow: 10452
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 43.25
Current gap: 43.25

Current alpha: 64.875
Max Flow: 11292
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 64.875
Current gap: 21.625

Current alpha: 54.0625
Max Flow: 11112
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 54.0625
Current gap: 10.8125

Current alpha: 59.4688
Max Flow: 11292
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 59.4688
Current gap: 5.40625

Current alpha: 56.7656
Max Flow: 11272
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 56.7656
Current gap: 2.70312

Current alpha: 58.1172
Max Flow: 11292
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 58.1172
Current gap: 1.35156

Current alpha: 57.4414
Max Flow: 11292
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 57.4414
Current gap: 0.675781

Current alpha: 57.1035
Max Flow: 11292
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 57.1035
Current gap: 0.337891

Current alpha: 56.9346
Max Flow: 11272
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 56.9346
Current gap: 0.168945

Current alpha: 57.019
Max Flow: 11292
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 57.019
Current gap: 0.0844727

Current alpha: 56.9768
Max Flow: 11272
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 56.9768
Current gap: 0.0422363

Current alpha: 56.9979
Max Flow: 11272
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 56.9979
Current gap: 0.0211182

Current alpha: 57.0085
Max Flow: 11292
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 57.0085
Current gap: 0.0105591

Current alpha: 57.0032
Max Flow: 11292
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 57.0032
Current gap: 0.00527954

Current alpha: 57.0006
Max Flow: 11292
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 57.0006
Current gap: 0.00263977

Current alpha: 56.9992
Max Flow: 11272
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 56.9992
Current gap: 0.00131989

Current alpha: 56.9999
Max Flow: 11272
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 56.9999
Current gap: 0.000659943

Current alpha: 57.0002
Max Flow: 11292
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 57.0002
Current gap: 0.000329971

Current alpha: 57.0001
Max Flow: 11292
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 57.0001
Current gap: 0.000164986

Current alpha: 57
Max Flow: 11272
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 57
Current gap: 8.24928e-05

Current alpha: 57
Max Flow: 11292
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 57
Current gap: 4.12464e-05

Current alpha: 57
Max Flow: 11292
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 57
Current gap: 2.06232e-05

Current alpha: 57
Max Flow: 11272
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 57
Current gap: 1.03116e-05

Current alpha: 57
Max Flow: 11292
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 57
Current gap: 5.1558e-06

Current alpha: 57
Max Flow: 11292
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 57
Current gap: 2.5779e-06

Current alpha: 57
Max Flow: 11272
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 57
Current gap: 1.28895e-06

Current alpha: 57
Max Flow: 11272
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 57
Current gap: 6.44475e-07

Current alpha: 57
Max Flow: 11272
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 57
Current gap: 3.22238e-07

Dense subgraph found with 20 vertices:
Vertices in dense subgraph: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Edges in dense subgraph:
Vertex 527 connected to: 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 529 connected to: 527 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 823 connected to: 527 529 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 1030 connected to: 527 529 823 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 1032 connected to: 527 529 823 1030 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 1034 connected to: 527 529 823 1030 1032 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 1036 connected to: 527 529 823 1030 1032 1034 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 1038 connected to: 527 529 823 1030 1032 1034 1036 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 1040 connected to: 527 529 823 1030 1032 1034 1036 1038 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 1042 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 1044 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 1046 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 1051 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1053 1055 1057 1059 1061 1068 1287
Vertex 1053 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1055 1057 1059 1061 1068 1287
Vertex 1055 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1057 1059 1061 1068 1287
Vertex 1057 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1059 1061 1068 1287
Vertex 1059 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1061 1068 1287
Vertex 1061 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1068 1287
Vertex 1068 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1287
Vertex 1287 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068
Number of edges: 190
Number of k-cliques in dense subgraph: 1140
Average k-cliques per vertex: 57
Execution time: 4204 milliseconds


```

For k=4

```bash
Graph with 1461 vertices and 2742 edges
Looking for cliques of size k = 4
Size 4: 7159

All 4-cliques (7159 in total):
Initial parameters:
l: 0, u: 970
Number of (k-1)-cliques: 3764


Current alpha: 485
Max Flow: 28636
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 485
Current gap: 485

Current alpha: 242.5
Max Flow: 28636
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 242.5
Current gap: 242.5

Current alpha: 121.25
Max Flow: 18956
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 121.25
Current gap: 121.25

Current alpha: 181.875
Max Flow: 23796
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 181.875
Current gap: 60.625

Current alpha: 212.188
Max Flow: 26216
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 212.188
Current gap: 30.3125

Current alpha: 227.344
Max Flow: 27436
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 227.344
Current gap: 15.1562

Current alpha: 234.922
Max Flow: 28036
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 234.922
Current gap: 7.57812

Current alpha: 238.711
Max Flow: 28336
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 238.711
Current gap: 3.78906

Current alpha: 240.605
Max Flow: 28496
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 240.605
Current gap: 1.89453

Current alpha: 241.553
Max Flow: 28576
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 241.553
Current gap: 0.947266

Current alpha: 242.026
Max Flow: 28616
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 242.026
Current gap: 0.473633

Current alpha: 242.263
Max Flow: 28636
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 242.263
Current gap: 0.236816

Current alpha: 242.145
Max Flow: 28616
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 242.145
Current gap: 0.118408

Current alpha: 242.204
Max Flow: 28616
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 242.204
Current gap: 0.0592041

Current alpha: 242.234
Max Flow: 28616
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 242.234
Current gap: 0.0296021

Current alpha: 242.248
Max Flow: 28616
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 242.248
Current gap: 0.014801

Current alpha: 242.256
Max Flow: 28636
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 242.256
Current gap: 0.00740051

Current alpha: 242.252
Max Flow: 28636
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 242.252
Current gap: 0.00370026

Current alpha: 242.25
Max Flow: 28636
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 242.25
Current gap: 0.00185013

Current alpha: 242.249
Max Flow: 28616
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 242.249
Current gap: 0.000925064

Current alpha: 242.25
Max Flow: 28616
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 242.25
Current gap: 0.000462532

Current alpha: 242.25
Max Flow: 28636
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 242.25
Current gap: 0.000231266

Current alpha: 242.25
Max Flow: 28616
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 242.25
Current gap: 0.000115633

Current alpha: 242.25
Max Flow: 28616
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 242.25
Current gap: 5.78165e-05

Current alpha: 242.25
Max Flow: 28616
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 242.25
Current gap: 2.89083e-05

Current alpha: 242.25
Max Flow: 28616
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 242.25
Current gap: 1.44541e-05

Current alpha: 242.25
Max Flow: 28616
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 242.25
Current gap: 7.22706e-06

Current alpha: 242.25
Max Flow: 28616
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 242.25
Current gap: 3.61353e-06

Current alpha: 242.25
Max Flow: 28616
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 242.25
Current gap: 1.80677e-06

Current alpha: 242.25
Max Flow: 28636
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 242.25
Current gap: 9.03383e-07

Current alpha: 242.25
Max Flow: 28636
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 242.25
Current gap: 4.51691e-07

Dense subgraph found with 20 vertices:
Vertices in dense subgraph: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Edges in dense subgraph:
Vertex 527 connected to: 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 529 connected to: 527 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 823 connected to: 527 529 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 1030 connected to: 527 529 823 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 1032 connected to: 527 529 823 1030 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 1034 connected to: 527 529 823 1030 1032 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 1036 connected to: 527 529 823 1030 1032 1034 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 1038 connected to: 527 529 823 1030 1032 1034 1036 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 1040 connected to: 527 529 823 1030 1032 1034 1036 1038 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 1042 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1044 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 1044 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1046 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 1046 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1051 1053 1055 1057 1059 1061 1068 1287
Vertex 1051 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1053 1055 1057 1059 1061 1068 1287
Vertex 1053 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1055 1057 1059 1061 1068 1287
Vertex 1055 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1057 1059 1061 1068 1287
Vertex 1057 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1059 1061 1068 1287
Vertex 1059 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1061 1068 1287
Vertex 1061 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1068 1287
Vertex 1068 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1287
Vertex 1287 connected to: 527 529 823 1030 1032 1034 1036 1038 1040 1042 1044 1046 1051 1053 1055 1057 1059 1061 1068
Number of edges: 190
Number of k-cliques in dense subgraph: 4845
Average k-cliques per vertex: 242.25
Execution time: 6724 milliseconds

```
