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

#### AS-733 Dataset

For k=2

```bash
Graph with 6474 vertices and 13895 edges
Looking for cliques of size k = 2
Size 2: 12572

All 2-cliques (12572 in total):
Initial parameters:
l: 0, u: 1458
Number of (k-1)-cliques: 6474


Current alpha: 729
Max Flow: 25144
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 729
Current gap: 729

Current alpha: 364.5
Max Flow: 25144
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 364.5
Current gap: 364.5

Current alpha: 182.25
Max Flow: 25144
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 182.25
Current gap: 182.25

Current alpha: 91.125
Max Flow: 25144
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 91.125
Current gap: 91.125

Current alpha: 45.5625
Max Flow: 25144
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 45.5625
Current gap: 45.5625

Current alpha: 22.7812
Max Flow: 25144
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 22.7812
Current gap: 22.7812

Current alpha: 11.3906
Max Flow: 25144
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 11.3906
Current gap: 11.3906

Current alpha: 5.69531
Max Flow: 24547
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 5.69531
Current gap: 5.69531

Current alpha: 8.54297
Max Flow: 25114
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 8.54297
Current gap: 2.84766

Current alpha: 9.9668
Max Flow: 25144
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 9.9668
Current gap: 1.42383

Current alpha: 9.25488
Max Flow: 25144
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 9.25488
Current gap: 0.711914

Current alpha: 8.89893
Max Flow: 25114
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 8.89893
Current gap: 0.355957

Current alpha: 9.0769
Max Flow: 25144
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 9.0769
Current gap: 0.177979

Current alpha: 8.98792
Max Flow: 25114
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 8.98792
Current gap: 0.0889893

Current alpha: 9.03241
Max Flow: 25144
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 9.03241
Current gap: 0.0444946

Current alpha: 9.01016
Max Flow: 25144
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 9.01016
Current gap: 0.0222473

Current alpha: 8.99904
Max Flow: 25114
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 8.99904
Current gap: 0.0111237

Current alpha: 9.0046
Max Flow: 25144
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 9.0046
Current gap: 0.00556183

Current alpha: 9.00182
Max Flow: 25144
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 9.00182
Current gap: 0.00278091

Current alpha: 9.00043
Max Flow: 25144
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 9.00043
Current gap: 0.00139046

Current alpha: 8.99973
Max Flow: 25114
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 8.99973
Current gap: 0.000695229

Current alpha: 9.00008
Max Flow: 25144
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 9.00008
Current gap: 0.000347614

Current alpha: 8.99991
Max Flow: 25114
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 8.99991
Current gap: 0.000173807

Current alpha: 8.99999
Max Flow: 25114
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 8.99999
Current gap: 8.69036e-05

Current alpha: 9.00004
Max Flow: 25144
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 9.00004
Current gap: 4.34518e-05

Current alpha: 9.00002
Max Flow: 25144
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 9.00002
Current gap: 2.17259e-05

Current alpha: 9.00001
Max Flow: 25144
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 9.00001
Current gap: 1.08629e-05

Current alpha: 9
Max Flow: 25144
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 9
Current gap: 5.43147e-06

Current alpha: 9
Max Flow: 25114
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 9
Current gap: 2.71574e-06

Current alpha: 9
Max Flow: 25114
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 9
Current gap: 1.35787e-06

Current alpha: 9
Max Flow: 25114
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 9
Current gap: 6.78934e-07

Current alpha: 9
Max Flow: 25114
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 9
Current gap: 3.39467e-07

Current alpha: 9
Max Flow: 25114
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 9
Current gap: 1.69734e-07

Current alpha: 9
Max Flow: 25114
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 9
Current gap: 8.48668e-08

Current alpha: 9
Max Flow: 25144
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 9
Current gap: 4.24334e-08

Current alpha: 9
Max Flow: 25144
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 9
Current gap: 2.12167e-08

Dense subgraph found with 40 vertices:
Vertices in dense subgraph: 1911 3116 3150 3155 3159 3176 3179 3248 3256 3258 3265 3275 3281 3287 3289 3292 3305 3321 3325 3330 3334 3338 3343 3348 3351 3354 3357 3362 3407 3419 3437 3444 3451 3470 3491 3498 3539 3618 3637 3692
Edges in dense subgraph:
Vertex 1911 connected to: 1911 1911 3330 3334 3338 3343 3348 3351 3354 3357 3419 3437
Vertex 3116 connected to: 3256 3321 3325 3330 3334 3338 3343 3348 3351 3354 3357 3362 3419 3498
Vertex 3150 connected to: 3150 3150 3176 3258 3321 3325 3330 3343 3348 3351 3354
Vertex 3155 connected to: 3321 3325 3330 3334 3338 3343 3348 3351 3354 3357 3419 3437
Vertex 3159 connected to: 3159 3159 3258 3265 3325 3330 3334 3338 3348 3354 3357
Vertex 3176 connected to: 3150 3258 3265 3321 3325 3330 3334 3338 3343 3348 3357
Vertex 3179 connected to: 3321 3325 3330 3334 3338 3343 3348 3351 3357
Vertex 3248 connected to: 3258 3321 3325 3330 3334 3343 3348 3351 3354 3357
Vertex 3256 connected to: 3116 3256 3256 3265 3321 3325 3330 3334 3338 3343 3348 3351 3354 3357
Vertex 3258 connected to: 3150 3159 3176 3248 3258 3258 3281 3287 3289 3325 3334 3343 3444 3451 3470 3491 3618 3637
Vertex 3265 connected to: 3159 3176 3256 3275 3287 3305 3321 3325 3330 3334 3343 3351 3357 3419 3451 3470 3637
Vertex 3275 connected to: 3265 3305 3321 3325 3330 3334 3338 3343 3351 3354 3357 3362 3407 3419 3437 3618 3692
Vertex 3281 connected to: 3258 3289 3321 3330 3334 3338 3343 3348 3351 3354 3357 3362 3419 3437 3692
Vertex 3287 connected to: 3258 3265 3321 3325 3330 3334 3338 3343 3348 3351 3354 3357 3362 3437
Vertex 3289 connected to: 3258 3281 3321 3325 3330 3334 3343 3348 3351 3354 3357 3618
Vertex 3292 connected to: 3321 3325 3330 3334 3338 3343 3348 3351 3357 3419 3437 3444
Vertex 3305 connected to: 3265 3275 3305 3305 3321 3325 3330 3334 3338 3343 3348 3351 3357 3419
Vertex 3321 connected to: 3116 3150 3155 3176 3179 3248 3256 3265 3275 3281 3287 3289 3292 3305 3325 3330 3334 3338 3343 3348 3351 3354 3357 3362 3407 3419 3437 3444 3451 3470 3491 3498 3539 3618 3637 3692
Vertex 3325 connected to: 3116 3150 3155 3159 3176 3179 3248 3256 3258 3265 3275 3287 3289 3292 3305 3321 3325 3325 3330 3334 3338 3343 3348 3351 3354 3357 3362 3407 3419 3437 3444 3451 3470 3491 3539 3618 3637 3692
Vertex 3330 connected to: 1911 3116 3150 3155 3159 3176 3179 3248 3256 3265 3275 3281 3287 3289 3292 3305 3321 3325 3334 3338 3343 3348 3351 3354 3357 3407 3419 3437 3444 3451 3470 3491 3498 3539 3618 3637 3692
Vertex 3334 connected to: 1911 3116 3155 3159 3176 3179 3248 3256 3258 3265 3275 3281 3287 3289 3292 3305 3321 3325 3330 3338 3343 3348 3351 3354 3357 3407 3419 3437 3444 3451 3470 3491 3539 3618 3637 3692
Vertex 3338 connected to: 1911 3116 3155 3159 3176 3179 3256 3275 3281 3287 3292 3305 3321 3325 3330 3334 3343 3348 3354 3362 3407 3419 3437 3451 3470 3491 3498 3539 3637 3692
Vertex 3343 connected to: 1911 3116 3150 3155 3176 3179 3248 3256 3258 3265 3275 3281 3287 3289 3292 3305 3321 3325 3330 3334 3338 3348 3351 3354 3357 3362 3407 3419 3437 3451 3470 3491 3498 3539 3637 3692
Vertex 3348 connected to: 1911 3116 3150 3155 3159 3176 3179 3248 3256 3281 3287 3289 3292 3305 3321 3325 3330 3334 3338 3343 3348 3348 3351 3354 3357 3362 3407 3419 3437 3444 3451 3470 3491 3498 3539 3618 3637 3692
Vertex 3351 connected to: 1911 3116 3150 3155 3179 3248 3256 3265 3275 3281 3287 3289 3292 3305 3321 3325 3330 3334 3343 3348 3357 3407 3437 3444 3451 3470 3491 3498 3539 3618 3637
Vertex 3354 connected to: 1911 3116 3150 3155 3159 3248 3256 3275 3281 3287 3289 3321 3325 3330 3334 3338 3343 3348 3354 3354 3357 3362 3437 3470 3491 3539 3618 3637
Vertex 3357 connected to: 1911 3116 3155 3159 3176 3179 3248 3256 3265 3275 3281 3287 3289 3292 3305 3321 3325 3330 3334 3343 3348 3351 3354 3357 3357 3362 3407 3419 3437 3444 3451 3470 3491 3498 3539 3637 3692
Vertex 3362 connected to: 3116 3275 3281 3287 3321 3325 3338 3343 3348 3354 3357 3362 3362 3419 3491 3498
Vertex 3407 connected to: 3275 3321 3325 3330 3334 3338 3343 3348 3351 3357 3407 3407
Vertex 3419 connected to: 1911 3116 3155 3265 3275 3281 3292 3305 3321 3325 3330 3334 3338 3343 3348 3357 3362 3437
Vertex 3437 connected to: 1911 3155 3275 3281 3287 3292 3321 3325 3330 3334 3338 3343 3348 3351 3354 3357 3419 3437 3437 3444 3618
Vertex 3444 connected to: 3258 3292 3321 3325 3330 3334 3348 3351 3357 3437 3692
Vertex 3451 connected to: 3258 3265 3321 3325 3330 3334 3338 3343 3348 3351 3357 3451 3451
Vertex 3470 connected to: 3258 3265 3321 3325 3330 3334 3338 3343 3348 3351 3354 3357
Vertex 3491 connected to: 3258 3321 3325 3330 3334 3338 3343 3348 3351 3354 3357 3362
Vertex 3498 connected to: 3116 3321 3330 3338 3343 3348 3351 3357 3362
Vertex 3539 connected to: 3321 3325 3330 3334 3338 3343 3348 3351 3354 3357 3539 3539
Vertex 3618 connected to: 3258 3275 3289 3321 3325 3330 3334 3348 3351 3354 3437 3618 3618
Vertex 3637 connected to: 3258 3265 3321 3325 3330 3334 3338 3343 3348 3351 3354 3357
Vertex 3692 connected to: 3275 3281 3321 3325 3330 3334 3338 3343 3348 3357 3444
Number of edges: 371
Number of k-cliques in dense subgraph: 355
Average k-cliques per vertex: 8.875
Execution time: 70489 milliseconds


```

For k=3

```bash
Graph with 6474 vertices and 13895 edges
Looking for cliques of size k = 3
Size 3: 6584

All 3-cliques (6584 in total):
Initial parameters:
l: 0, u: 2047
Number of (k-1)-cliques: 12572


Current alpha: 1023.5
Max Flow: 19752
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 1023.5
Current gap: 1023.5

Current alpha: 511.75
Max Flow: 19752
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 511.75
Current gap: 511.75

Current alpha: 255.875
Max Flow: 19752
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 255.875
Current gap: 255.875

Current alpha: 127.938
Max Flow: 19752
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 127.938
Current gap: 127.938

Current alpha: 63.9688
Max Flow: 19752
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 63.9688
Current gap: 63.9688

Current alpha: 31.9844
Max Flow: 19300
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 31.9844
Current gap: 31.9844

Current alpha: 47.9766
Max Flow: 19752
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 47.9766
Current gap: 15.9922

Current alpha: 39.9805
Max Flow: 19752
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 39.9805
Current gap: 7.99609

Current alpha: 35.9824
Max Flow: 19728
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 35.9824
Current gap: 3.99805

Current alpha: 37.9814
Max Flow: 19752
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 37.9814
Current gap: 1.99902

Current alpha: 36.9819
Max Flow: 19752
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 36.9819
Current gap: 0.999512

Current alpha: 36.4822
Max Flow: 19752
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 36.4822
Current gap: 0.499756

Current alpha: 36.2323
Max Flow: 19752
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 36.2323
Current gap: 0.249878

Current alpha: 36.1074
Max Flow: 19752
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 36.1074
Current gap: 0.124939

Current alpha: 36.0449
Max Flow: 19752
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 36.0449
Current gap: 0.0624695

Current alpha: 36.0137
Max Flow: 19752
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 36.0137
Current gap: 0.0312347

Current alpha: 35.998
Max Flow: 19728
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 35.998
Current gap: 0.0156174

Current alpha: 36.0058
Max Flow: 19752
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 36.0058
Current gap: 0.00780869

Current alpha: 36.0019
Max Flow: 19752
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 36.0019
Current gap: 0.00390434

Current alpha: 36
Max Flow: 19728
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 36
Current gap: 0.00195217

Current alpha: 36.001
Max Flow: 19752
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 36.001
Current gap: 0.000976086

Current alpha: 36.0005
Max Flow: 19752
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 36.0005
Current gap: 0.000488043

Current alpha: 36.0002
Max Flow: 19752
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 36.0002
Current gap: 0.000244021

Current alpha: 36.0001
Max Flow: 19752
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 36.0001
Current gap: 0.000122011

Current alpha: 36.0001
Max Flow: 19752
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 36.0001
Current gap: 6.10054e-05

Current alpha: 36
Max Flow: 19752
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 36
Current gap: 3.05027e-05

Current alpha: 36
Max Flow: 19752
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 36
Current gap: 1.52513e-05

Current alpha: 36
Max Flow: 19728
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 36
Current gap: 7.62567e-06

Current alpha: 36
Max Flow: 19752
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 36
Current gap: 3.81283e-06

Current alpha: 36
Max Flow: 19752
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 36
Current gap: 1.90642e-06

Current alpha: 36
Max Flow: 19728
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 36
Current gap: 9.53209e-07

Current alpha: 36
Max Flow: 19752
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 36
Current gap: 4.76604e-07

Current alpha: 36
Max Flow: 19752
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 36
Current gap: 2.38302e-07

Current alpha: 36
Max Flow: 19752
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 36
Current gap: 1.19151e-07

Current alpha: 36
Max Flow: 19752
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 36
Current gap: 5.95755e-08

Current alpha: 36
Max Flow: 19752
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 36
Current gap: 2.97878e-08

Current alpha: 36
Max Flow: 19752
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 36
Current gap: 1.48939e-08

Dense subgraph found with 33 vertices:
Vertices in dense subgraph: 1911 3116 3155 3256 3265 3275 3281 3287 3289 3292 3305 3321 3325 3330 3334 3338 3343 3348 3351 3354 3357 3362 3407 3419 3437 3444 3451 3470 3491 3539 3618 3637 3692
Edges in dense subgraph:
Vertex 1911 connected to: 1911 1911 3330 3334 3338 3343 3348 3351 3354 3357 3419 3437
Vertex 3116 connected to: 3256 3321 3325 3330 3334 3338 3343 3348 3351 3354 3357 3362 3419
Vertex 3155 connected to: 3321 3325 3330 3334 3338 3343 3348 3351 3354 3357 3419 3437
Vertex 3256 connected to: 3116 3256 3256 3265 3321 3325 3330 3334 3338 3343 3348 3351 3354 3357
Vertex 3265 connected to: 3256 3275 3287 3305 3321 3325 3330 3334 3343 3351 3357 3419 3451 3470 3637
Vertex 3275 connected to: 3265 3305 3321 3325 3330 3334 3338 3343 3351 3354 3357 3362 3407 3419 3437 3618 3692
Vertex 3281 connected to: 3289 3321 3330 3334 3338 3343 3348 3351 3354 3357 3362 3419 3437 3692
Vertex 3287 connected to: 3265 3321 3325 3330 3334 3338 3343 3348 3351 3354 3357 3362 3437
Vertex 3289 connected to: 3281 3321 3325 3330 3334 3343 3348 3351 3354 3357 3618
Vertex 3292 connected to: 3321 3325 3330 3334 3338 3343 3348 3351 3357 3419 3437 3444
Vertex 3305 connected to: 3265 3275 3305 3305 3321 3325 3330 3334 3338 3343 3348 3351 3357 3419
Vertex 3321 connected to: 3116 3155 3256 3265 3275 3281 3287 3289 3292 3305 3325 3330 3334 3338 3343 3348 3351 3354 3357 3362 3407 3419 3437 3444 3451 3470 3491 3539 3618 3637 3692
Vertex 3325 connected to: 3116 3155 3256 3265 3275 3287 3289 3292 3305 3321 3325 3325 3330 3334 3338 3343 3348 3351 3354 3357 3362 3407 3419 3437 3444 3451 3470 3491 3539 3618 3637 3692
Vertex 3330 connected to: 1911 3116 3155 3256 3265 3275 3281 3287 3289 3292 3305 3321 3325 3334 3338 3343 3348 3351 3354 3357 3407 3419 3437 3444 3451 3470 3491 3539 3618 3637 3692
Vertex 3334 connected to: 1911 3116 3155 3256 3265 3275 3281 3287 3289 3292 3305 3321 3325 3330 3338 3343 3348 3351 3354 3357 3407 3419 3437 3444 3451 3470 3491 3539 3618 3637 3692
Vertex 3338 connected to: 1911 3116 3155 3256 3275 3281 3287 3292 3305 3321 3325 3330 3334 3343 3348 3354 3362 3407 3419 3437 3451 3470 3491 3539 3637 3692
Vertex 3343 connected to: 1911 3116 3155 3256 3265 3275 3281 3287 3289 3292 3305 3321 3325 3330 3334 3338 3348 3351 3354 3357 3362 3407 3419 3437 3451 3470 3491 3539 3637 3692
Vertex 3348 connected to: 1911 3116 3155 3256 3281 3287 3289 3292 3305 3321 3325 3330 3334 3338 3343 3348 3348 3351 3354 3357 3362 3407 3419 3437 3444 3451 3470 3491 3539 3618 3637 3692
Vertex 3351 connected to: 1911 3116 3155 3256 3265 3275 3281 3287 3289 3292 3305 3321 3325 3330 3334 3343 3348 3357 3407 3437 3444 3451 3470 3491 3539 3618 3637
Vertex 3354 connected to: 1911 3116 3155 3256 3275 3281 3287 3289 3321 3325 3330 3334 3338 3343 3348 3354 3354 3357 3362 3437 3470 3491 3539 3618 3637
Vertex 3357 connected to: 1911 3116 3155 3256 3265 3275 3281 3287 3289 3292 3305 3321 3325 3330 3334 3343 3348 3351 3354 3357 3357 3362 3407 3419 3437 3444 3451 3470 3491 3539 3637 3692
Vertex 3362 connected to: 3116 3275 3281 3287 3321 3325 3338 3343 3348 3354 3357 3362 3362 3419 3491
Vertex 3407 connected to: 3275 3321 3325 3330 3334 3338 3343 3348 3351 3357 3407 3407
Vertex 3419 connected to: 1911 3116 3155 3265 3275 3281 3292 3305 3321 3325 3330 3334 3338 3343 3348 3357 3362 3437
Vertex 3437 connected to: 1911 3155 3275 3281 3287 3292 3321 3325 3330 3334 3338 3343 3348 3351 3354 3357 3419 3437 3437 3444 3618
Vertex 3444 connected to: 3292 3321 3325 3330 3334 3348 3351 3357 3437 3692
Vertex 3451 connected to: 3265 3321 3325 3330 3334 3338 3343 3348 3351 3357 3451 3451
Vertex 3470 connected to: 3265 3321 3325 3330 3334 3338 3343 3348 3351 3354 3357
Vertex 3491 connected to: 3321 3325 3330 3334 3338 3343 3348 3351 3354 3357 3362
Vertex 3539 connected to: 3321 3325 3330 3334 3338 3343 3348 3351 3354 3357 3539 3539
Vertex 3618 connected to: 3275 3289 3321 3325 3330 3334 3348 3351 3354 3437 3618 3618
Vertex 3637 connected to: 3265 3321 3325 3330 3334 3338 3343 3348 3351 3354 3357
Vertex 3692 connected to: 3275 3281 3321 3325 3330 3334 3338 3343 3348 3357 3444
Number of edges: 300
Number of k-cliques in dense subgraph: 1185
Average k-cliques per vertex: 35.9091
Execution time: 176099 milliseconds

```

For k=4

```bash
Graph with 6474 vertices and 13895 edges
Looking for cliques of size k = 4
Size 4: 5636

All 4-cliques (5636 in total):
Initial parameters:
l: 0, u: 2059
Number of (k-1)-cliques: 6584


Current alpha: 1029.5
Max Flow: 22544
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 1029.5
Current gap: 1029.5

Current alpha: 514.75
Max Flow: 22544
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 514.75
Current gap: 514.75

Current alpha: 257.375
Max Flow: 22544
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 257.375
Current gap: 257.375

Current alpha: 128.688
Max Flow: 22544
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 128.688
Current gap: 128.688

Current alpha: 64.3438
Max Flow: 19708
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 64.3438
Current gap: 64.3438

Current alpha: 96.5156
Max Flow: 22544
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 96.5156
Current gap: 32.1719

Current alpha: 80.4297
Max Flow: 21901
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 80.4297
Current gap: 16.0859

Current alpha: 88.4727
Max Flow: 22544
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 88.4727
Current gap: 8.04297

Current alpha: 84.4512
Max Flow: 22429
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 84.4512
Current gap: 4.02148

Current alpha: 86.4619
Max Flow: 22544
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 86.4619
Current gap: 2.01074

Current alpha: 85.4565
Max Flow: 22544
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 85.4565
Current gap: 1.00537

Current alpha: 84.9539
Max Flow: 22495
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 84.9539
Current gap: 0.502686

Current alpha: 85.2052
Max Flow: 22528
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 85.2052
Current gap: 0.251343

Current alpha: 85.3309
Max Flow: 22544
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 85.3309
Current gap: 0.125671

Current alpha: 85.268
Max Flow: 22544
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 85.268
Current gap: 0.0628357

Current alpha: 85.2366
Max Flow: 22528
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 85.2366
Current gap: 0.0314178

Current alpha: 85.2523
Max Flow: 22544
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 85.2523
Current gap: 0.0157089

Current alpha: 85.2445
Max Flow: 22528
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 85.2445
Current gap: 0.00785446

Current alpha: 85.2484
Max Flow: 22528
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 85.2484
Current gap: 0.00392723

Current alpha: 85.2504
Max Flow: 22544
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 85.2504
Current gap: 0.00196362

Current alpha: 85.2494
Max Flow: 22528
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 85.2494
Current gap: 0.000981808

Current alpha: 85.2499
Max Flow: 22528
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 85.2499
Current gap: 0.000490904

Current alpha: 85.2501
Max Flow: 22544
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 85.2501
Current gap: 0.000245452

Current alpha: 85.25
Max Flow: 22528
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 85.25
Current gap: 0.000122726

Current alpha: 85.2501
Max Flow: 22544
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 85.2501
Current gap: 6.1363e-05

Current alpha: 85.25
Max Flow: 22544
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 85.25
Current gap: 3.06815e-05

Current alpha: 85.25
Max Flow: 22544
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 85.25
Current gap: 1.53407e-05

Current alpha: 85.25
Max Flow: 22544
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 85.25
Current gap: 7.67037e-06

Current alpha: 85.25
Max Flow: 22528
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 85.25
Current gap: 3.83519e-06

Current alpha: 85.25
Max Flow: 22544
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 85.25
Current gap: 1.91759e-06

Current alpha: 85.25
Max Flow: 22528
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 85.25
Current gap: 9.58797e-07

Current alpha: 85.25
Max Flow: 22544
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 85.25
Current gap: 4.79398e-07

Current alpha: 85.25
Max Flow: 22544
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 85.25
Current gap: 2.39699e-07

Current alpha: 85.25
Max Flow: 22544
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 85.25
Current gap: 1.1985e-07

Current alpha: 85.25
Max Flow: 22544
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 85.25
Current gap: 5.99248e-08

Current alpha: 85.25
Max Flow: 22528
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 85.25
Current gap: 2.99624e-08

Current alpha: 85.25
Max Flow: 22544
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 85.25
Current gap: 1.49812e-08

Dense subgraph found with 32 vertices:
Vertices in dense subgraph: 3116 3155 3256 3265 3275 3281 3287 3289 3292 3305 3321 3325 3330 3334 3338 3343 3348 3351 3354 3357 3362 3407 3419 3437 3444 3451 3470 3491 3539 3618 3637 3692
Edges in dense subgraph:
Vertex 3116 connected to: 3256 3321 3325 3330 3334 3338 3343 3348 3351 3354 3357 3362 3419
Vertex 3155 connected to: 3321 3325 3330 3334 3338 3343 3348 3351 3354 3357 3419 3437
Vertex 3256 connected to: 3116 3256 3256 3265 3321 3325 3330 3334 3338 3343 3348 3351 3354 3357
Vertex 3265 connected to: 3256 3275 3287 3305 3321 3325 3330 3334 3343 3351 3357 3419 3451 3470 3637
Vertex 3275 connected to: 3265 3305 3321 3325 3330 3334 3338 3343 3351 3354 3357 3362 3407 3419 3437 3618 3692
Vertex 3281 connected to: 3289 3321 3330 3334 3338 3343 3348 3351 3354 3357 3362 3419 3437 3692
Vertex 3287 connected to: 3265 3321 3325 3330 3334 3338 3343 3348 3351 3354 3357 3362 3437
Vertex 3289 connected to: 3281 3321 3325 3330 3334 3343 3348 3351 3354 3357 3618
Vertex 3292 connected to: 3321 3325 3330 3334 3338 3343 3348 3351 3357 3419 3437 3444
Vertex 3305 connected to: 3265 3275 3305 3305 3321 3325 3330 3334 3338 3343 3348 3351 3357 3419
Vertex 3321 connected to: 3116 3155 3256 3265 3275 3281 3287 3289 3292 3305 3325 3330 3334 3338 3343 3348 3351 3354 3357 3362 3407 3419 3437 3444 3451 3470 3491 3539 3618 3637 3692
Vertex 3325 connected to: 3116 3155 3256 3265 3275 3287 3289 3292 3305 3321 3325 3325 3330 3334 3338 3343 3348 3351 3354 3357 3362 3407 3419 3437 3444 3451 3470 3491 3539 3618 3637 3692
Vertex 3330 connected to: 3116 3155 3256 3265 3275 3281 3287 3289 3292 3305 3321 3325 3334 3338 3343 3348 3351 3354 3357 3407 3419 3437 3444 3451 3470 3491 3539 3618 3637 3692
Vertex 3334 connected to: 3116 3155 3256 3265 3275 3281 3287 3289 3292 3305 3321 3325 3330 3338 3343 3348 3351 3354 3357 3407 3419 3437 3444 3451 3470 3491 3539 3618 3637 3692
Vertex 3338 connected to: 3116 3155 3256 3275 3281 3287 3292 3305 3321 3325 3330 3334 3343 3348 3354 3362 3407 3419 3437 3451 3470 3491 3539 3637 3692
Vertex 3343 connected to: 3116 3155 3256 3265 3275 3281 3287 3289 3292 3305 3321 3325 3330 3334 3338 3348 3351 3354 3357 3362 3407 3419 3437 3451 3470 3491 3539 3637 3692
Vertex 3348 connected to: 3116 3155 3256 3281 3287 3289 3292 3305 3321 3325 3330 3334 3338 3343 3348 3348 3351 3354 3357 3362 3407 3419 3437 3444 3451 3470 3491 3539 3618 3637 3692
Vertex 3351 connected to: 3116 3155 3256 3265 3275 3281 3287 3289 3292 3305 3321 3325 3330 3334 3343 3348 3357 3407 3437 3444 3451 3470 3491 3539 3618 3637
Vertex 3354 connected to: 3116 3155 3256 3275 3281 3287 3289 3321 3325 3330 3334 3338 3343 3348 3354 3354 3357 3362 3437 3470 3491 3539 3618 3637
Vertex 3357 connected to: 3116 3155 3256 3265 3275 3281 3287 3289 3292 3305 3321 3325 3330 3334 3343 3348 3351 3354 3357 3357 3362 3407 3419 3437 3444 3451 3470 3491 3539 3637 3692
Vertex 3362 connected to: 3116 3275 3281 3287 3321 3325 3338 3343 3348 3354 3357 3362 3362 3419 3491
Vertex 3407 connected to: 3275 3321 3325 3330 3334 3338 3343 3348 3351 3357 3407 3407
Vertex 3419 connected to: 3116 3155 3265 3275 3281 3292 3305 3321 3325 3330 3334 3338 3343 3348 3357 3362 3437
Vertex 3437 connected to: 3155 3275 3281 3287 3292 3321 3325 3330 3334 3338 3343 3348 3351 3354 3357 3419 3437 3437 3444 3618
Vertex 3444 connected to: 3292 3321 3325 3330 3334 3348 3351 3357 3437 3692
Vertex 3451 connected to: 3265 3321 3325 3330 3334 3338 3343 3348 3351 3357 3451 3451
Vertex 3470 connected to: 3265 3321 3325 3330 3334 3338 3343 3348 3351 3354 3357
Vertex 3491 connected to: 3321 3325 3330 3334 3338 3343 3348 3351 3354 3357 3362
Vertex 3539 connected to: 3321 3325 3330 3334 3338 3343 3348 3351 3354 3357 3539 3539
Vertex 3618 connected to: 3275 3289 3321 3325 3330 3334 3348 3351 3354 3437 3618 3618
Vertex 3637 connected to: 3265 3321 3325 3330 3334 3338 3343 3348 3351 3354 3357
Vertex 3692 connected to: 3275 3281 3321 3325 3330 3334 3338 3343 3348 3357 3444
Number of edges: 289
Number of k-cliques in dense subgraph: 2724
Average k-cliques per vertex: 85.125
Execution time: 157578 milliseconds

```

#### Hepth Dataset

For k=2

```bash
Graph with 9877 vertices and 51971 edges
Looking for cliques of size k = 2
Size 2: 25973

All 2-cliques (25973 in total):
Initial parameters:
l: 0, u: 65
Number of (k-1)-cliques: 9877


Current alpha: 32.5
Max Flow: 51946
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 32.5
Current gap: 32.5

Current alpha: 16.25
Max Flow: 51946
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 16.25
Current gap: 16.25

Current alpha: 8.125
Max Flow: 51176
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 8.125
Current gap: 8.125

Current alpha: 12.1875
Max Flow: 51722
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 12.1875
Current gap: 4.0625

Current alpha: 14.2188
Max Flow: 51850
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 14.2188
Current gap: 2.03125

Current alpha: 15.2344
Max Flow: 51914
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 15.2344
Current gap: 1.01562

Current alpha: 15.7422
Max Flow: 51946
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 15.7422
Current gap: 0.507812

Current alpha: 15.4883
Max Flow: 51914
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 15.4883
Current gap: 0.253906

Current alpha: 15.6152
Max Flow: 51946
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 15.6152
Current gap: 0.126953

Current alpha: 15.5518
Max Flow: 51946
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 15.5518
Current gap: 0.0634766

Current alpha: 15.52
Max Flow: 51946
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 15.52
Current gap: 0.0317383

Current alpha: 15.5042
Max Flow: 51946
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 15.5042
Current gap: 0.0158691

Current alpha: 15.4962
Max Flow: 51914
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 15.4962
Current gap: 0.00793457

Current alpha: 15.5002
Max Flow: 51946
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 15.5002
Current gap: 0.00396729

Current alpha: 15.4982
Max Flow: 51914
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 15.4982
Current gap: 0.00198364

Current alpha: 15.4992
Max Flow: 51914
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 15.4992
Current gap: 0.000991821

Current alpha: 15.4997
Max Flow: 51914
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 15.4997
Current gap: 0.000495911

Current alpha: 15.4999
Max Flow: 51914
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 15.4999
Current gap: 0.000247955

Current alpha: 15.5001
Max Flow: 51946
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 15.5001
Current gap: 0.000123978

Current alpha: 15.5
Max Flow: 51914
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 15.5
Current gap: 6.19888e-05

Current alpha: 15.5
Max Flow: 51946
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 15.5
Current gap: 3.09944e-05

Current alpha: 15.5
Max Flow: 51946
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 15.5
Current gap: 1.54972e-05

Current alpha: 15.5
Max Flow: 51946
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 15.5
Current gap: 7.7486e-06

Current alpha: 15.5
Max Flow: 51946
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 15.5
Current gap: 3.8743e-06

Current alpha: 15.5
Max Flow: 51914
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 15.5
Current gap: 1.93715e-06

Current alpha: 15.5
Max Flow: 51946
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 15.5
Current gap: 9.68575e-07

Current alpha: 15.5
Max Flow: 51914
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 15.5
Current gap: 4.84288e-07

Current alpha: 15.5
Max Flow: 51914
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 15.5
Current gap: 2.42144e-07

Current alpha: 15.5
Max Flow: 51914
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 15.5
Current gap: 1.21072e-07

Current alpha: 15.5
Max Flow: 51914
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 15.5
Current gap: 6.0536e-08

Current alpha: 15.5
Max Flow: 51946
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 15.5
Current gap: 3.0268e-08

Current alpha: 15.5
Max Flow: 51914
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 15.5
Current gap: 1.5134e-08

Current alpha: 15.5
Max Flow: 51946
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 15.5
Current gap: 7.567e-09

Dense subgraph found with 32 vertices:
Vertices in dense subgraph: 3834 4243 4249 4477 5154 5365 5399 5686 6170 6635 7068 7775 7955 8422 8536 8602 8618 8670 9021 9067 9379 9428 9429 9430 9433 9434 9439 9440 9441 9443 9444 9445
Edges in dense subgraph:
Vertex 3834 connected to: 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 4243 connected to: 3834 3834 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 4249 connected to: 3834 3834 4243 4243 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 4477 connected to: 3834 3834 4243 4243 4249 4249 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 5154 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 5365 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 5399 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 5686 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 6170 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 6635 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 7068 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 7775 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 7955 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 8422 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 8536 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 8602 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 8618 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 8670 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 9021 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 9067 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 9379 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 9428 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 9429 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 9430 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 9433 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 9434 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 9439 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 9440 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 9441 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9443 9443 9444 9444 9445 9445
Vertex 9443 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9444 9444 9445 9445
Vertex 9444 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9445 9445
Vertex 9445 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444
Number of edges: 992
Number of k-cliques in dense subgraph: 496
Average k-cliques per vertex: 15.5
Execution time: 121067 milliseconds

```

For k=3

```bash
Graph with 9877 vertices and 51971 edges
Looking for cliques of size k = 3
Size 3: 28339

All 3-cliques (28339 in total):
Initial parameters:
l: 0, u: 505
Number of (k-1)-cliques: 25973


Current alpha: 252.5
Max Flow: 85017
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 252.5
Current gap: 252.5

Current alpha: 126.25
Max Flow: 82233
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 126.25
Current gap: 126.25

Current alpha: 189.375
Max Flow: 85017
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 189.375
Current gap: 63.125

Current alpha: 157.812
Max Flow: 85017
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 157.812
Current gap: 31.5625

Current alpha: 142.031
Max Flow: 83769
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 142.031
Current gap: 15.7812

Current alpha: 149.922
Max Flow: 84505
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 149.922
Current gap: 7.89062

Current alpha: 153.867
Max Flow: 84889
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 153.867
Current gap: 3.94531

Current alpha: 155.84
Max Flow: 85017
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 155.84
Current gap: 1.97266

Current alpha: 154.854
Max Flow: 84985
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 154.854
Current gap: 0.986328

Current alpha: 155.347
Max Flow: 85017
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 155.347
Current gap: 0.493164

Current alpha: 155.1
Max Flow: 85017
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 155.1
Current gap: 0.246582

Current alpha: 154.977
Max Flow: 84985
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 154.977
Current gap: 0.123291

Current alpha: 155.038
Max Flow: 85017
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 155.038
Current gap: 0.0616455

Current alpha: 155.008
Max Flow: 85017
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 155.008
Current gap: 0.0308228

Current alpha: 154.992
Max Flow: 84985
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 154.992
Current gap: 0.0154114

Current alpha: 155
Max Flow: 84985
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 155
Current gap: 0.00770569

Current alpha: 155.004
Max Flow: 85017
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 155.004
Current gap: 0.00385284

Current alpha: 155.002
Max Flow: 85017
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 155.002
Current gap: 0.00192642

Current alpha: 155.001
Max Flow: 85017
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 155.001
Current gap: 0.000963211

Current alpha: 155
Max Flow: 85017
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 155
Current gap: 0.000481606

Current alpha: 155
Max Flow: 85017
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 155
Current gap: 0.000240803

Current alpha: 155
Max Flow: 85017
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 155
Current gap: 0.000120401

Current alpha: 155
Max Flow: 84985
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 155
Current gap: 6.02007e-05

Current alpha: 155
Max Flow: 85017
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 155
Current gap: 3.01003e-05

Current alpha: 155
Max Flow: 84985
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 155
Current gap: 1.50502e-05

Current alpha: 155
Max Flow: 85017
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 155
Current gap: 7.52509e-06

Current alpha: 155
Max Flow: 85017
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 155
Current gap: 3.76254e-06

Current alpha: 155
Max Flow: 85017
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 155
Current gap: 1.88127e-06

Current alpha: 155
Max Flow: 84985
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 155
Current gap: 9.40636e-07

Current alpha: 155
Max Flow: 85017
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 155
Current gap: 4.70318e-07

Current alpha: 155
Max Flow: 85017
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 155
Current gap: 2.35159e-07

Current alpha: 155
Max Flow: 85017
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 155
Current gap: 1.17579e-07

Current alpha: 155
Max Flow: 84985
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 155
Current gap: 5.87897e-08

Current alpha: 155
Max Flow: 84985
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 155
Current gap: 2.93949e-08

Current alpha: 155
Max Flow: 85017
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 155
Current gap: 1.46974e-08

Current alpha: 155
Max Flow: 84985
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 155
Current gap: 7.34872e-09

Dense subgraph found with 32 vertices:
Vertices in dense subgraph: 3834 4243 4249 4477 5154 5365 5399 5686 6170 6635 7068 7775 7955 8422 8536 8602 8618 8670 9021 9067 9379 9428 9429 9430 9433 9434 9439 9440 9441 9443 9444 9445
Edges in dense subgraph:
Vertex 3834 connected to: 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 4243 connected to: 3834 3834 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 4249 connected to: 3834 3834 4243 4243 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 4477 connected to: 3834 3834 4243 4243 4249 4249 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 5154 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 5365 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 5399 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 5686 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 6170 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 6635 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 7068 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 7775 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 7955 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 8422 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 8536 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 8602 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 8618 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 8670 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 9021 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 9067 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 9379 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 9428 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 9429 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 9430 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 9433 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 9434 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 9439 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9440 9440 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 9440 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9441 9441 9443 9443 9444 9444 9445 9445
Vertex 9441 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9443 9443 9444 9444 9445 9445
Vertex 9443 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9444 9444 9445 9445
Vertex 9444 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9445 9445
Vertex 9445 connected to: 3834 3834 4243 4243 4249 4249 4477 4477 5154 5154 5365 5365 5399 5399 5686 5686 6170 6170 6635 6635 7068 7068 7775 7775 7955 7955 8422 8422 8536 8536 8602 8602 8618 8618 8670 8670 9021 9021 9067 9067 9379 9379 9428 9428 9429 9429 9430 9430 9433 9433 9434 9434 9439 9439 9440 9440 9441 9441 9443 9443 9444 9444
Number of edges: 992
Number of k-cliques in dense subgraph: 4960
Average k-cliques per vertex: 155
Execution time: 347503 milliseconds

```

#### Caida Dataset

For k=2

```bash
Current alpha: 1314
Max Flow: 106762
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 1314
Current gap: 1314

Current alpha: 657
Max Flow: 106762
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 657
Current gap: 657

Current alpha: 328.5
Max Flow: 106762
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 328.5
Current gap: 328.5

Current alpha: 164.25
Max Flow: 106762
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 164.25
Current gap: 164.25

Current alpha: 82.125
Max Flow: 106762
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 82.125
Current gap: 82.125

Current alpha: 41.0625
Max Flow: 106762
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 41.0625
Current gap: 41.0625

Current alpha: 20.5312
Max Flow: 106762
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 20.5312
Current gap: 20.5312

Current alpha: 10.2656
Max Flow: 104678
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 10.2656
Current gap: 10.2656

Current alpha: 15.3984
Max Flow: 106238
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 15.3984
Current gap: 5.13281

Current alpha: 17.9648
Max Flow: 106756
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 17.9648
Current gap: 2.56641

Current alpha: 19.248
Max Flow: 106762
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 19.248
Current gap: 1.2832

Current alpha: 18.6064
Max Flow: 106762
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 18.6064
Current gap: 0.641602

Current alpha: 18.2856
Max Flow: 106762
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 18.2856
Current gap: 0.320801

Current alpha: 18.1252
Max Flow: 106762
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 18.1252
Current gap: 0.1604

Current alpha: 18.045
Max Flow: 106762
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 18.045
Current gap: 0.0802002

Current alpha: 18.0049
Max Flow: 106762
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 18.0049
Current gap: 0.0401001

Current alpha: 17.9849
Max Flow: 106756
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 17.9849
Current gap: 0.02005

Current alpha: 17.9949
Max Flow: 106756
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 17.9949
Current gap: 0.010025

Current alpha: 17.9999
Max Flow: 106756
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 17.9999
Current gap: 0.00501251

Current alpha: 18.0024
Max Flow: 106762
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 18.0024
Current gap: 0.00250626

Current alpha: 18.0012
Max Flow: 106762
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 18.0012
Current gap: 0.00125313

Current alpha: 18.0006
Max Flow: 106762
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 18.0006
Current gap: 0.000626564

Current alpha: 18.0002
Max Flow: 106762
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 18.0002
Current gap: 0.000313282

Current alpha: 18.0001
Max Flow: 106762
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 18.0001
Current gap: 0.000156641

Current alpha: 18
Max Flow: 106762
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 18
Current gap: 7.83205e-005

Current alpha: 18
Max Flow: 106756
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 18
Current gap: 3.91603e-005

Current alpha: 18
Max Flow: 106756
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 18
Current gap: 1.95801e-005

Current alpha: 18
Max Flow: 106756
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 18
Current gap: 9.79006e-006

Current alpha: 18
Max Flow: 106762
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 18
Current gap: 4.89503e-006

Current alpha: 18
Max Flow: 106762
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 18
Current gap: 2.44752e-006

Current alpha: 18
Max Flow: 106762
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 18
Current gap: 1.22376e-006

Current alpha: 18
Max Flow: 106762
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 18
Current gap: 6.11879e-007

Current alpha: 18
Max Flow: 106762
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 18
Current gap: 3.05939e-007

Current alpha: 18
Max Flow: 106762
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 18
Current gap: 1.5297e-007

Current alpha: 18
Max Flow: 106756
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 18
Current gap: 7.64849e-008

Current alpha: 18
Max Flow: 106756
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 18
Current gap: 3.82424e-008

Current alpha: 18
Max Flow: 106756
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 18
Current gap: 1.91212e-008

Current alpha: 18
Max Flow: 106762
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 18
Current gap: 9.56061e-009

Current alpha: 18
Max Flow: 106762
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 18
Current gap: 4.7803e-009

Current alpha: 18
Max Flow: 106762
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 18
Current gap: 2.39015e-009

Current alpha: 18
Max Flow: 106762
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 18
Current gap: 1.19508e-009


Dense subgraph found with 88 vertices:
Vertices in dense subgraph: 9522 9582 9903 10164 10282 10665 11201 11237 11301 11861 13390 13458 13558 13575 13886 14473 14511 14518 14546 14592 14642 14927 14939 14943 15167 15282 15455 15485 15508 15535 15618 15623 15698 15699 15737 15754 15833 15839 15970 16003 16044 16099 16152 16170 16276 16343 16353 16397 16408 16480 16486 16586 16594 16768 16985 17196 17366 17440 17538 17601 17830 18002 18064 18134 18144 18145 18328 18401 18407 18413 18498 18551 18807 18864 19029 19101 19493 19522 19620 19760 20854 22045 22257 22957 24750 25145 25296 25621
Edges in dense subgraph:
Vertex 9522 connected to: 9582 9582 10665 10665 11237 11237 11301 11301 13390 13390 13558 13558 13886 13886 14473 14473 14642 14642 14939 14939 14943 14943 15167 15167 15508 15508 15535 15535 15623 15623 15754 15754 15833 15833 15839 15839 16044 16044 16099 16099 16170 16170 16276 16276 16397 16397 16408 16408 16480 16480 16486 16486 16586 16586 16768 16768 16985 16985 17196 17196 17366 17366 17601 17601 17830 17830 18002 18002 18401 18401 18407 18407 18413 18413 18498 18498 18551 18551 18807 18807 19522 19522 19620 19620 20854 20854 22257 22257 25145 25145 25296 25296 25621 25621
Vertex 9582 connected to: 9522 9522 9903 9903 10164 10164 10282 10282 11237 11237 13558 13558 13886 13886 14473 14473 14518 14518 14546 14546 14592 14592 14927 14927 14939 14939 14943 14943 15282 15282 15455 15455 15485 15485 15623 15623 15698 15698 15699 15699 15737 15737 15754 15754 16003 16003 16044 16044 16099 16099 16408 16408 16586 16586 16768 16768 17601 17601 18145 18145 18401 18401 18407 18407 18413 18413 18498 18498 18551 18551 18864 18864 19101 19101 19620 19620 20854 20854 22045 22045 22957 22957
Vertex 9903 connected to: 9582 9582 10665 10665 11237 11237 13390 13390 14642 14642 14939 14939 14943 14943 15623 15623 15698 15698 15754 15754 16044 16044 16170 16170 16408 16408 16586 16586 17538 17538 18145 18145 19029 19029 19101 19101 19522 19522 19620 19620 19760 19760 20854 20854
Vertex 10164 connected to: 9582 9582 10665 10665 11861 11861 13390 13390 13458 13458 13558 13558 13886 13886 14473 14473 14943 14943 15485 15485 15618 15618 15623 15623 15737 15737 15833 15833 15970 15970 16044 16044 16170 16170 16276 16276 16408 16408 16594 16594 16768 16768 16985 16985 17601 17601 18328 18328 18407 18407 18498 18498 18551 18551 19029 19029 19760 19760 20854 20854
Vertex 10282 connected to: 9582 9582 10665 10665 13390 13390 13458 13458 14473 14473 14943 14943 15508 15508 15623 15623 15737 15737 15833 15833 16044 16044 16343 16343 16353 16353 16397 16397 16408 16408 16586 16586 16594 16594 16985 16985 17601 17601 20854 20854
Vertex 10665 connected to: 9522 9522 9903 9903 10164 10164 10282 10282 11201 11201 11237 11237 11861 11861 13390 13390 13558 13558 13575 13575 14473 14473 14546 14546 14592 14592 14642 14642 14939 14939 14943 14943 15167 15167 15485 15485 15508 15508 15535 15535 15623 15623 15698 15698 15737 15737 15754 15754 15833 15833 15839 15839 15970 15970 16044 16044 16099 16099 16170 16170 16397 16397 16486 16486 16586 16586 16594 16594 16985 16985 17440 17440 17538 17538 18002 18002 18134 18134 18144 18144 18145 18145 18328 18328 18401 18401 18498 18498 18864 18864 19522 19522 19620 19620 22045 22045 22257 22257 22957 22957 24750 24750 25296 25296
Vertex 11201 connected to: 10665 10665 13390 13390 13458 13458 13886 13886 14943 14943 15623 15623 15970 15970 16003 16003 16044 16044 16099 16099 16343 16343 16408 16408 16985 16985 17538 17538 17601 17601 18498 18498 19029 19029 19522 19522 20854 20854
Vertex 11237 connected to: 9522 9522 9582 9582 9903 9903 10665 10665 11861 11861 13558 13558 13575 13575 14473 14473 14592 14592 14927 14927 14939 14939 15167 15167 15508 15508 15535 15535 15623 15623 15698 15698 15699 15699 15737 15737 15754 15754 15833 15833 16099 16099 16170 16170 16276 16276 16397 16397 16408 16408 16985 16985 17196 17196 17366 17366 17538 17538 17601 17601 18002 18002 18064 18064 18145 18145 18328 18328 18401 18401 18407 18407 18864 18864 19522 19522 19620 19620 19760 19760 22045 22045 22257 22257 22957 22957 25621 25621
Vertex 11301 connected to: 9522 9522 13390 13390 13458 13458 14473 14473 14939 14939 14943 14943 15167 15167 15535 15535 15833 15833 15970 15970 16003 16003 16044 16044 16099 16099 16170 16170 16586 16586 17601 17601 18864 18864 19029 19029 19760 19760
Vertex 11861 connected to: 10164 10164 10665 10665 11237 11237 13390 13390 13886 13886 14473 14473 15167 15167 15485 15485 15508 15508 15623 15623 15698 15698 16003 16003 16099 16099 16170 16170 16343 16343 16353 16353 16397 16397 17196 17196 17538 17538 17601 17601 19522 19522 20854 20854 22257 22257
Vertex 13390 connected to: 9522 9522 9903 9903 10164 10164 10282 10282 10665 10665 11201 11201 11301 11301 11861 11861 13458 13458 14511 14511 14518 14518 14642 14642 14927 14927 14939 14939 15282 15282 15455 15455 15485 15485 15535 15535 15618 15618 15623 15623 15698 15698 15699 15699 15737 15737 15754 15754 15833 15833 15839 15839 15970 15970 16003 16003 16044 16044 16152 16152 16276 16276 16353 16353 16408 16408 16480 16480 16486 16486 16594 16594 16768 16768 16985 16985 17538 17538 17830 17830 18002 18002 18064 18064 18134 18134 18144 18144 18413 18413 18498 18498 18807 18807 18864 18864 19029 19029 19493 19493 19522 19522 19760 19760 20854 20854 22045 22045 22257 22257 22957 22957 25145 25145 25296 25296 25621 25621
Vertex 13458 connected to: 10164 10164 10282 10282 11201 11201 11301 11301 13390 13390 14473 14473 14511 14511 14642 14642 15282 15282 15455 15455 15485 15485 15535 15535 15698 15698 15699 15699 15833 15833 15839 15839 15970 15970 16044 16044 16353 16353 16408 16408 16480 16480 16594 16594 16768 16768 17366 17366 17830 17830 18002 18002 18064 18064 18134 18134 18144 18144 18413 18413 18807 18807 18864 18864 19029 19029 22257 22257 25296 25296 25621 25621
Vertex 13558 connected to: 9522 9522 9582 9582 10164 10164 10665 10665 11237 11237 13575 13575 14473 14473 14546 14546 14642 14642 14939 14939 14943 14943 15167 15167 15508 15508 15535 15535 15698 15698 15699 15699 15754 15754 15833 15833 16044 16044 16099 16099 16170 16170 16397 16397 16408 16408 16586 16586 16985 16985 17366 17366 17440 17440 17601 17601 18145 18145 18328 18328 18401 18401 18407 18407 18498 18498 18551 18551 19522 19522 19620 19620 20854 20854
Vertex 13575 connected to: 10665 10665 11237 11237 13558 13558 14473 14473 14939 14939 14943 14943 15485 15485 15535 15535 15623 15623 15698 15698 15737 15737 16044 16044 16099 16099 16408 16408 16586 16586 18145 18145 18498 18498 19620 19620 20854 20854
Vertex 13886 connected to: 9522 9522 9582 9582 10164 10164 11201 11201 11861 11861 14546 14546 14927 14927 15282 15282 15455 15455 15618 15618 15623 15623 15698 15698 15699 15699 15737 15737 16044 16044 16353 16353 16594 16594 16768 16768 17440 17440 18064 18064 18145 18145 18407 18407 18413 18413 18498 18498 18551 18551 18864 18864 20854 20854
Vertex 14473 connected to: 9522 9522 9582 9582 10164 10164 10282 10282 10665 10665 11237 11237 11301 11301 11861 11861 13458 13458 13558 13558 13575 13575 14518 14518 14927 14927 14939 14939 14943 14943 15167 15167 15455 15455 15508 15508 15535 15535 15698 15698 15699 15699 15737 15737 15754 15754 16003 16003 16044 16044 16099 16099 16170 16170 16276 16276 16397 16397 16480 16480 16586 16586 16985 16985 17366 17366 17601 17601 17830 17830 18002 18002 18134 18134 18144 18144 18145 18145 18328 18328 18401 18401 18407 18407 18498 18498 18864 18864 19493 19493 19522 19522 19620 19620 20854 20854 22257 22257 22957 22957 25145 25145 25296 25296
Vertex 14511 connected to: 13390 13390 13458 13458 14943 14943 15167 15167 15508 15508 15623 15623 15698 15698 15754 15754 15833 15833 15970 15970 16044 16044 16099 16099 16170 16170 16408 16408 16586 16586 16985 16985 17196 17196 17538 17538 17601 17601 19029 19029 19522 19522 19760 19760 24750 24750
Vertex 14518 connected to: 9582 9582 13390 13390 14473 14473 14939 14939 14943 14943 15167 15167 15535 15535 15623 15623 15754 15754 15833 15833 15970 15970 16044 16044 16099 16099 16170 16170 16353 16353 16408 16408 16586 16586 16985 16985 17196 17196 17538 17538 18498 18498 19101 19101 19760 19760 20854 20854 22257 22257
Vertex 14546 connected to: 9582 9582 10665 10665 13558 13558 13886 13886 14642 14642 15167 15167 15508 15508 15535 15535 15623 15623 15698 15698 15737 15737 15754 15754 16044 16044 16099 16099 16170 16170 16353 16353 16408 16408 16985 16985 17601 17601 18498 18498 20854 20854
Vertex 14592 connected to: 9582 9582 10665 10665 11237 11237 14943 14943 15508 15508 15623 15623 15737 15737 15833 15833 15839 15839 16044 16044 16170 16170 16343 16343 16353 16353 16397 16397 16586 16586 17538 17538 17601 17601 19101 19101 19620 19620 19760 19760 20854 20854
Vertex 14642 connected to: 9522 9522 9903 9903 10665 10665 13390 13390 13458 13458 13558 13558 14546 14546 14939 14939 15167 15167 15485 15485 15535 15535 15623 15623 15698 15698 15833 15833 15970 15970 16003 16003 16044 16044 16099 16099 16343 16343 16353 16353 16586 16586 16985 16985 17196 17196 17440 17440 17601 17601 17830 17830 18401 18401 19029 19029 20854 20854 22045 22045 24750 24750
Vertex 14927 connected to: 9582 9582 11237 11237 13390 13390 13886 13886 14473 14473 15508 15508 15754 15754 15833 15833 15970 15970 16003 16003 16044 16044 16170 16170 16353 16353 16397 16397 16408 16408 16586 16586 17538 17538 17601 17601 19029 19029 19101 19101 19522 19522 19760 19760 20854 20854 24750 24750
Vertex 14939 connected to: 9522 9522 9582 9582 9903 9903 10665 10665 11237 11237 11301 11301 13390 13390 13558 13558 13575 13575 14473 14473 14518 14518 14642 14642 14943 14943 15167 15167 15485 15485 15508 15508 15535 15535 15698 15698 15699 15699 15737 15737 15754 15754 15833 15833 15839 15839 15970 15970 16044 16044 16099 16099 16170 16170 16276 16276 16353 16353 16397 16397 16408 16408 16486 16486 16586 16586 16985 16985 17196 17196 17366 17366 17601 17601 17830 17830 18002 18002 18064 18064 18134 18134 18144 18144 18145 18145 18401 18401 18413 18413 18864 18864 19493 19493 19522 19522 19620 19620 19760 19760 22257 22257 25145 25145 25296 25296 25621 25621
Vertex 14943 connected to: 9522 9522 9582 9582 9903 9903 10164 10164 10282 10282 10665 10665 11201 11201 11301 11301 13558 13558 13575 13575 14473 14473 14511 14511 14518 14518 14592 14592 14939 14939 15167 15167 15485 15485 15508 15508 15535 15535 15698 15698 15699 15699 15737 15737 15754 15754 15833 15833 15839 15839 16044 16044 16099 16099 16170 16170 16276 16276 16353 16353 16397 16397 16408 16408 16586 16586 16985 16985 17196 17196 17366 17366 17601 17601 18401 18401 18407 18407 18413 18413 18551 18551 18864 18864 20854 20854 22045 22045 22257 22257 25296 25296 25621 25621
Vertex 15167 connected to: 9522 9522 10665 10665 11237 11237 11301 11301 11861 11861 13558 13558 14473 14473 14511 14511 14518 14518 14546 14546 14642 14642 14939 14939 14943 14943 15508 15508 15535 15535 15698 15698 15754 15754 15839 15839 16099 16099 16170 16170 16397 16397 16408 16408 16586 16586 16594 16594 16768 16768 16985 16985 17196 17196 17366 17366 17601 17601 17830 17830 18144 18144 18145 18145 18401 18401 18413 18413 18807 18807 18864 18864 19522 19522 19620 19620 20854 20854 22257 22257 25145 25145 25621 25621
Vertex 15282 connected to: 9582 9582 13390 13390 13458 13458 13886 13886 15508 15508 15623 15623 15833 15833 15970 15970 16044 16044 16170 16170 16343 16343 16397 16397 16408 16408 17538 17538 17601 17601 18498 18498 19029 19029 19101 19101 19522 19522 19760 19760 20854 20854 24750 24750
Vertex 15455 connected to: 9582 9582 13390 13390 13458 13458 13886 13886 14473 14473 15508 15508 15623 15623 15833 15833 15970 15970 16044 16044 16170 16170 16343 16343 16586 16586 17538 17538 17601 17601 18498 18498 19029 19029 19101 19101 19620 19620 19760 19760 20854 20854 24750 24750
Vertex 15485 connected to: 9582 9582 10164 10164 10665 10665 11861 11861 13390 13390 13458 13458 13575 13575 14642 14642 14939 14939 14943 14943 15508 15508 15535 15535 15618 15618 15623 15623 15833 15833 15970 15970 16003 16003 16044 16044 16099 16099 16276 16276 16343 16343 16353 16353 16408 16408 16480 16480 16586 16586 17538 17538 17601 17601 17830 17830 18134 18134 18145 18145 18328 18328 18498 18498 19029 19029 19101 19101 19522 19522 19620 19620 19760 19760 20854 20854 22257 22257 24750 24750
Vertex 15508 connected to: 9522 9522 10282 10282 10665 10665 11237 11237 11861 11861 13558 13558 14473 14473 14511 14511 14546 14546 14592 14592 14927 14927 14939 14939 14943 14943 15167 15167 15282 15282 15455 15455 15485 15485 15535 15535 15618 15618 15698 15698 15754 15754 15970 15970 16044 16044 16099 16099 16152 16152 16170 16170 16397 16397 16586 16586 16985 16985 17440 17440 18144 18144 18145 18145 18401 18401 18498 18498 18551 18551 18807 18807 18864 18864 19493 19493 19522 19522 19620 19620 22257 22257 22957 22957 25621 25621
Vertex 15535 connected to: 9522 9522 10665 10665 11237 11237 11301 11301 13390 13390 13458 13458 13558 13558 13575 13575 14473 14473 14518 14518 14546 14546 14642 14642 14939 14939 14943 14943 15167 15167 15485 15485 15508 15508 15698 15698 15699 15699 15754 15754 15833 15833 15839 15839 15970 15970 16044 16044 16099 16099 16170 16170 16353 16353 16397 16397 16408 16408 16480 16480 16486 16486 16586 16586 16768 16768 16985 16985 17196 17196 17366 17366 17440 17440 17601 17601 18002 18002 18134 18134 18144 18144 18145 18145 18328 18328 18401 18401 18407 18407 18413 18413 18498 18498 18551 18551 18864 18864 19029 19029 19101 19101 19493 19493 19522 19522 19620 19620 20854 20854 22045 22045 22257 22257 22957 22957 25296 25296 25621 25621
Vertex 15618 connected to: 10164 10164 13390 13390 13886 13886 15485 15485 15508 15508 15623 15623 15698 15698 15737 15737 15970 15970 16003 16003 16170 16170 16343 16343 16397 16397 16408 16408 17196 17196 17601 17601 18498 18498 19029 19029 20854 20854
Vertex 15623 connected to: 9522 9522 9582 9582 9903 9903 10164 10164 10282 10282 10665 10665 11201 11201 11237 11237 11861 11861 13390 13390 13575 13575 13886 13886 14511 14511 14518 14518 14546 14546 14592 14592 14642 14642 15282 15282 15455 15455 15485 15485 15618 15618 15698 15698 15699 15699 15737 15737 15754 15754 15839 15839 16003 16003 16044 16044 16099 16099 16152 16152 16276 16276 16343 16343 16353 16353 16408 16408 16480 16480 16486 16486 16594 16594 16768 16768 17196 17196 17366 17366 17538 17538 17601 17601 17830 17830 18002 18002 18064 18064 18134 18134 18401 18401 18407 18407 18413 18413 18498 18498 18551 18551 18807 18807 18864 18864 19029 19029 19101 19101 19493 19493 19760 19760 20854 20854 22045 22045 22957 22957
Vertex 15698 connected to: 9582 9582 9903 9903 10665 10665 11237 11237 11861 11861 13390 13390 13458 13458 13558 13558 13575 13575 13886 13886 14473 14473 14511 14511 14546 14546 14642 14642 14939 14939 14943 14943 15167 15167 15508 15508 15535 15535 15618 15618 15623 15623 15737 15737 15833 15833 15839 15839 15970 15970 16044 16044 16099 16099 16152 16152 16170 16170 16276 16276 16343 16343 16353 16353 16397 16397 16408 16408 16480 16480 16486 16486 16586 16586 16985 16985 17196 17196 17538 17538 17601 17601 17830 17830 18002 18002 18134 18134 18145 18145 18328 18328 18401 18401 18407 18407 18498 18498 19029 19029 19101 19101 19522 19522 19760 19760 20854 20854 24750 24750 25145 25145 25296 25296 25621 25621
Vertex 15699 connected to: 9582 9582 11237 11237 13390 13390 13458 13458 13558 13558 13886 13886 14473 14473 14939 14939 14943 14943 15535 15535 15623 15623 15754 15754 15833 15833 15970 15970 16044 16044 16099 16099 16170 16170 16343 16343 16408 16408 16586 16586 17538 17538 17601 17601 18498 18498 19029 19029 19101 19101 19522 19522 19620 19620 19760 19760 20854 20854
Vertex 15737 connected to: 9582 9582 10164 10164 10282 10282 10665 10665 11237 11237 13390 13390 13575 13575 13886 13886 14473 14473 14546 14546 14592 14592 14939 14939 14943 14943 15618 15618 15623 15623 15698 15698 15754 15754 15839 15839 16003 16003 16044 16044 16152 16152 16343 16343 16397 16397 16408 16408 16586 16586 16985 16985 17196 17196 17538 17538 17601 17601 17830 17830 18145 18145 18401 18401 18407 18407 18498 18498 19029 19029 19101 19101 19522 19522 19620 19620 19760 19760 20854 20854 22045 22045 25296 25296
Vertex 15754 connected to: 9522 9522 9582 9582 9903 9903 10665 10665 11237 11237 13390 13390 13558 13558 14473 14473 14511 14511 14518 14518 14546 14546 14927 14927 14939 14939 14943 14943 15167 15167 15508 15508 15535 15535 15623 15623 15699 15699 15737 15737 15833 15833 15970 15970 16003 16003 16044 16044 16099 16099 16152 16152 16170 16170 16343 16343 16397 16397 16586 16586 16594 16594 16985 16985 17196 17196 17538 17538 17601 17601 18002 18002 18145 18145 18498 18498 19029 19029 19101 19101 19493 19493 19522 19522 19620 19620 19760 19760 20854 20854 22257 22257 25145 25145
Vertex 15833 connected to: 9522 9522 10164 10164 10282 10282 10665 10665 11237 11237 11301 11301 13390 13390 13458 13458 13558 13558 14511 14511 14518 14518 14592 14592 14642 14642 14927 14927 14939 14939 14943 14943 15282 15282 15455 15455 15485 15485 15535 15535 15698 15698 15699 15699 15754 15754 15839 15839 15970 15970 16003 16003 16044 16044 16099 16099 16152 16152 16353 16353 16480 16480 16486 16486 16586 16586 16594 16594 16768 16768 16985 16985 17196 17196 17366 17366 17601 17601 17830 17830 18002 18002 18064 18064 18134 18134 18144 18144 18328 18328 18407 18407 18413 18413 18807 18807 18864 18864 19029 19029 19620 19620 19760 19760 22257 22257 25145 25145 25296 25296 25621 25621
Vertex 15839 connected to: 9522 9522 10665 10665 13390 13390 13458 13458 14592 14592 14939 14939 14943 14943 15167 15167 15535 15535 15623 15623 15698 15698 15737 15737 15833 15833 15970 15970 16044 16044 16099 16099 16170 16170 16353 16353 16408 16408 16586 16586 16985 16985 18145 18145 18498 18498 19029 19029 19101 19101 19760 19760 20854 20854 22257 22257 25145 25145
Vertex 15970 connected to: 10164 10164 10665 10665 11201 11201 11301 11301 13390 13390 13458 13458 14511 14511 14518 14518 14642 14642 14927 14927 14939 14939 15282 15282 15455 15455 15485 15485 15508 15508 15535 15535 15618 15618 15698 15698 15699 15699 15754 15754 15833 15833 15839 15839 16003 16003 16044 16044 16099 16099 16152 16152 16353 16353 16480 16480 16486 16486 16768 16768 17196 17196 17366 17366 17830 17830 18002 18002 18064 18064 18134 18134 18144 18144 18145 18145 18328 18328 18407 18407 18413 18413 18807 18807 19029 19029 19760 19760 22257 22257 24750 24750 25145 25145 25296 25296 25621 25621
Vertex 16003 connected to: 9582 9582 11201 11201 11301 11301 11861 11861 13390 13390 14473 14473 14642 14642 14927 14927 15485 15485 15618 15618 15623 15623 15737 15737 15754 15754 15833 15833 15970 15970 16044 16044 16152 16152 16276 16276 16343 16343 16353 16353 16480 16480 16486 16486 16768 16768 16985 16985 17538 17538 18002 18002 18144 18144 18328 18328 18401 18401 18413 18413 18498 18498 19029 19029 19101 19101 19493 19493 19522 19522 19760 19760 20854 20854 24750 24750 25145 25145 25621 25621
Vertex 16044 connected to: 9522 9522 9582 9582 9903 9903 10164 10164 10282 10282 10665 10665 11201 11201 11301 11301 13390 13390 13458 13458 13558 13558 13575 13575 13886 13886 14473 14473 14511 14511 14518 14518 14546 14546 14592 14592 14642 14642 14927 14927 14939 14939 14943 14943 15282 15282 15455 15455 15485 15485 15508 15508 15535 15535 15623 15623 15698 15698 15699 15699 15737 15737 15754 15754 15833 15833 15839 15839 15970 15970 16003 16003 16099 16099 16170 16170 16276 16276 16343 16343 16353 16353 16408 16408 16480 16480 16486 16486 16586 16586 16594 16594 16768 16768 16985 16985 17196 17196 17366 17366 17440 17440 17601 17601 17830 17830 18002 18002 18064 18064 18134 18134 18144 18144 18145 18145 18328 18328 18401 18401 18407 18407 18413 18413 18498 18498 18551 18551 18807 18807 18864 18864 19029 19029 19101 19101 19760 19760 20854 20854 22045 22045 22957 22957 25145 25145 25296 25296 25621 25621
Vertex 16099 connected to: 9522 9522 9582 9582 10665 10665 11201 11201 11237 11237 11301 11301 11861 11861 13558 13558 13575 13575 14473 14473 14511 14511 14518 14518 14546 14546 14642 14642 14939 14939 14943 14943 15167 15167 15485 15485 15508 15508 15535 15535 15623 15623 15698 15698 15699 15699 15754 15754 15833 15833 15839 15839 15970 15970 16044 16044 16152 16152 16170 16170 16276 16276 16397 16397 16408 16408 16486 16486 16586 16586 16594 16594 16985 16985 17196 17196 17366 17366 17440 17440 17601 17601 17830 17830 18002 18002 18144 18144 18145 18145 18401 18401 18407 18407 18413 18413 18498 18498 18551 18551 18864 18864 19101 19101 19522 19522 19620 19620 20854 20854 22257 22257 25296 25296 25621 25621
Vertex 16152 connected to: 13390 13390 15508 15508 15623 15623 15698 15698 15737 15737 15754 15754 15833 15833 15970 15970 16003 16003 16099 16099 16353 16353 16397 16397 17196 17196 18498 18498 19029 19029 19760 19760 20854 20854 22257 22257 24750 24750
Vertex 16170 connected to: 9522 9522 9903 9903 10164 10164 10665 10665 11237 11237 11301 11301 11861 11861 13558 13558 14473 14473 14511 14511 14518 14518 14546 14546 14592 14592 14927 14927 14939 14939 14943 14943 15167 15167 15282 15282 15455 15455 15508 15508 15535 15535 15618 15618 15698 15698 15699 15699 15754 15754 15839 15839 16044 16044 16099 16099 16397 16397 16408 16408 16480 16480 16486 16486 16586 16586 16594 16594 16985 16985 17196 17196 17366 17366 17440 17440 17538 17538 17601 17601 17830 17830 18002 18002 18064 18064 18134 18134 18144 18144 18145 18145 18407 18407 18413 18413 18498 18498 18551 18551 18807 18807 18864 18864 19029 19029 19493 19493 19522 19522 19620 19620 20854 20854 22257 22257 24750 24750 25145 25145 25296 25296 25621 25621
Vertex 16276 connected to: 9522 9522 10164 10164 11237 11237 13390 13390 14473 14473 14939 14939 14943 14943 15485 15485 15623 15623 15698 15698 16003 16003 16044 16044 16099 16099 16985 16985 17196 17196 17366 17366 18145 18145 18407 18407 18498 18498 19029 19029 19760 19760 20854 20854 24750 24750
Vertex 16343 connected to: 10282 10282 11201 11201 11861 11861 14592 14592 14642 14642 15282 15282 15455 15455 15485 15485 15618 15618 15623 15623 15698 15698 15699 15699 15737 15737 15754 15754 16003 16003 16044 16044 16353 16353 16408 16408 16594 16594 16768 16768 17538 17538 17601 17601 18064 18064 18145 18145 18407 18407 18413 18413 18498 18498 18551 18551 18864 18864 22045 22045 22957 22957
Vertex 16353 connected to: 10282 10282 11861 11861 13390 13390 13458 13458 13886 13886 14518 14518 14546 14546 14592 14592 14642 14642 14927 14927 14939 14939 14943 14943 15485 15485 15535 15535 15623 15623 15698 15698 15833 15833 15839 15839 15970 15970 16003 16003 16044 16044 16152 16152 16343 16343 16486 16486 16586 16586 17196 17196 17440 17440 17538 17538 17601 17601 17830 17830 18002 18002 18328 18328 18401 18401 18498 18498 18551 18551 19029 19029 19101 19101 19493 19493 19522 19522 19760 19760 20854 20854 22957 22957 24750 24750 25296 25296 25621 25621
Vertex 16397 connected to: 9522 9522 10282 10282 10665 10665 11237 11237 11861 11861 13558 13558 14473 14473 14592 14592 14927 14927 14939 14939 14943 14943 15167 15167 15282 15282 15508 15508 15535 15535 15618 15618 15698 15698 15737 15737 15754 15754 16099 16099 16152 16152 16170 16170 16408 16408 16586 16586 16768 16768 16985 16985 17366 17366 18002 18002 18145 18145 18401 18401 18413 18413 18498 18498 18864 18864 19493 19493 19522 19522 19620 19620 19760 19760 22957 22957 25621 25621
Vertex 16408 connected to: 9522 9522 9582 9582 9903 9903 10164 10164 10282 10282 11201 11201 11237 11237 13390 13390 13458 13458 13558 13558 13575 13575 14511 14511 14518 14518 14546 14546 14927 14927 14939 14939 14943 14943 15167 15167 15282 15282 15485 15485 15535 15535 15618 15618 15623 15623 15698 15698 15699 15699 15737 15737 15839 15839 16044 16044 16099 16099 16170 16170 16343 16343 16397 16397 16586 16586 16594 16594 16985 16985 17196 17196 17366 17366 17601 17601 17830 17830 18145 18145 18401 18401 18407 18407 18498 18498 18551 18551 19029 19029 19101 19101 19493 19493 19522 19522 19620 19620 20854 20854 22045 22045 25296 25296
Vertex 16480 connected to: 9522 9522 13390 13390 13458 13458 14473 14473 15485 15485 15535 15535 15623 15623 15698 15698 15833 15833 15970 15970 16003 16003 16044 16044 16170 16170 16985 16985 17196 17196 18498 18498 19029 19029 19760 19760 20854 20854 22257 22257 24750 24750
Vertex 16486 connected to: 9522 9522 10665 10665 13390 13390 14939 14939 15535 15535 15623 15623 15698 15698 15833 15833 15970 15970 16003 16003 16044 16044 16099 16099 16170 16170 16353 16353 18407 18407 18498 18498 19029 19029 19760 19760 20854 20854 24750 24750
Vertex 16586 connected to: 9522 9522 9582 9582 9903 9903 10282 10282 10665 10665 11301 11301 13558 13558 13575 13575 14473 14473 14511 14511 14518 14518 14592 14592 14642 14642 14927 14927 14939 14939 14943 14943 15167 15167 15455 15455 15485 15485 15508 15508 15535 15535 15698 15698 15699 15699 15737 15737 15754 15754 15833 15833 15839 15839 16044 16044 16099 16099 16170 16170 16353 16353 16397 16397 16408 16408 16594 16594 16768 16768 16985 16985 17196 17196 17366 17366 18134 18134 18144 18144 18401 18401 18407 18407 18413 18413 18551 18551 18864 18864 19101 19101 19620 19620 19760 19760 20854 20854 22045 22045 22257 22257 25145 25145 25296 25296 25621 25621
Vertex 16594 connected to: 10164 10164 10282 10282 10665 10665 13390 13390 13458 13458 13886 13886 15167 15167 15623 15623 15754 15754 15833 15833 16044 16044 16099 16099 16170 16170 16343 16343 16408 16408 16586 16586 17538 17538 18401 18401 18498 18498 19029 19029 20854 20854 24750 24750
Vertex 16768 connected to: 9522 9522 9582 9582 10164 10164 13390 13390 13458 13458 13886 13886 15167 15167 15535 15535 15623 15623 15833 15833 15970 15970 16003 16003 16044 16044 16343 16343 16397 16397 16586 16586 16985 16985 17538 17538 17601 17601 18498 18498 19029 19029 19101 19101 19760 19760 20854 20854 24750 24750
Vertex 16985 connected to: 9522 9522 10164 10164 10282 10282 10665 10665 11201 11201 11237 11237 13390 13390 13558 13558 14473 14473 14511 14511 14518 14518 14546 14546 14642 14642 14939 14939 14943 14943 15167 15167 15508 15508 15535 15535 15698 15698 15737 15737 15754 15754 15833 15833 15839 15839 16003 16003 16044 16044 16099 16099 16170 16170 16276 16276 16397 16397 16408 16408 16480 16480 16586 16586 16768 16768 17196 17196 17366 17366 17440 17440 18002 18002 18134 18134 18145 18145 18401 18401 18413 18413 18498 18498 18807 18807 18864 18864 19029 19029 19493 19493 19522 19522 19620 19620 20854 20854 22045 22045 22257 22257 25145 25145 25296 25296 25621 25621
Vertex 17196 connected to: 9522 9522 11237 11237 11861 11861 14511 14511 14518 14518 14642 14642 14939 14939 14943 14943 15167 15167 15535 15535 15618 15618 15623 15623 15698 15698 15737 15737 15754 15754 15833 15833 15970 15970 16044 16044 16099 16099 16152 16152 16170 16170 16276 16276 16353 16353 16408 16408 16480 16480 16586 16586 16985 16985 17366 17366 17538 17538 17601 17601 17830 17830 18064 18064 18498 18498 18551 18551 19101 19101 19493 19493 19522 19522 19760 19760 22257 22257 25145 25145
Vertex 17366 connected to: 9522 9522 11237 11237 13458 13458 13558 13558 14473 14473 14939 14939 14943 14943 15167 15167 15535 15535 15623 15623 15833 15833 15970 15970 16044 16044 16099 16099 16170 16170 16276 16276 16397 16397 16408 16408 16586 16586 16985 16985 17196 17196 17538 17538 17601 17601 18328 18328 18498 18498 18807 18807 19029 19029 19101 19101 19522 19522 19620 19620 19760 19760 20854 20854 22257 22257
Vertex 17440 connected to: 10665 10665 13558 13558 13886 13886 14642 14642 15508 15508 15535 15535 16044 16044 16099 16099 16170 16170 16353 16353 16985 16985 17538 17538 17601 17601 18145 18145 19101 19101 19522 19522 19760 19760 20854 20854
Vertex 17538 connected to: 9903 9903 10665 10665 11201 11201 11237 11237 11861 11861 13390 13390 14511 14511 14518 14518 14592 14592 14927 14927 15282 15282 15455 15455 15485 15485 15623 15623 15698 15698 15699 15699 15737 15737 15754 15754 16003 16003 16170 16170 16343 16343 16353 16353 16594 16594 16768 16768 17196 17196 17366 17366 17440 17440 18064 18064 18328 18328 18407 18407 18413 18413 18498 18498 18551 18551 18807 18807 18864 18864 19029 19029 19493 19493 19760 19760 20854 20854 22045 22045 22957 22957 25145 25145
Vertex 17601 connected to: 9522 9522 9582 9582 10164 10164 10282 10282 11201 11201 11237 11237 11301 11301 11861 11861 13558 13558 14473 14473 14511 14511 14546 14546 14592 14592 14642 14642 14927 14927 14939 14939 14943 14943 15167 15167 15282 15282 15455 15455 15485 15485 15535 15535 15618 15618 15623 15623 15698 15698 15699 15699 15737 15737 15754 15754 15833 15833 16044 16044 16099 16099 16170 16170 16343 16343 16353 16353 16408 16408 16768 16768 17196 17196 17366 17366 17440 17440 18002 18002 18064 18064 18145 18145 18328 18328 18401 18401 18407 18407 18413 18413 18498 18498 18551 18551 18807 18807 18864 18864 19620 19620 20854 20854 22045 22045 22957 22957
Vertex 17830 connected to: 9522 9522 13390 13390 13458 13458 14473 14473 14642 14642 14939 14939 15167 15167 15485 15485 15623 15623 15698 15698 15737 15737 15833 15833 15970 15970 16044 16044 16099 16099 16170 16170 16353 16353 16408 16408 17196 17196 18498 18498 19029 19029 19493 19493 19760 19760 20854 20854
Vertex 18002 connected to: 9522 9522 10665 10665 11237 11237 13390 13390 13458 13458 14473 14473 14939 14939 15535 15535 15623 15623 15698 15698 15754 15754 15833 15833 15970 15970 16003 16003 16044 16044 16099 16099 16170 16170 16353 16353 16397 16397 16985 16985 17601 17601 18145 18145 18498 18498 19620 19620 20854 20854 22257 22257 25145 25145
Vertex 18064 connected to: 11237 11237 13390 13390 13458 13458 13886 13886 14939 14939 15623 15623 15833 15833 15970 15970 16044 16044 16170 16170 16343 16343 17196 17196 17538 17538 17601 17601 18498 18498 19029 19029 19101 19101 19522 19522 19760 19760 20854 20854
Vertex 18134 connected to: 10665 10665 13390 13390 13458 13458 14473 14473 14939 14939 15485 15485 15535 15535 15623 15623 15698 15698 15833 15833 15970 15970 16044 16044 16170 16170 16586 16586 16985 16985 19029 19029 19760 19760 20854 20854
Vertex 18144 connected to: 10665 10665 13390 13390 13458 13458 14473 14473 14939 14939 15167 15167 15508 15508 15535 15535 15833 15833 15970 15970 16003 16003 16044 16044 16099 16099 16170 16170 16586 16586 19029 19029 19760 19760 22257 22257
Vertex 18145 connected to: 9582 9582 9903 9903 10665 10665 11237 11237 13558 13558 13575 13575 13886 13886 14473 14473 14939 14939 15167 15167 15485 15485 15508 15508 15535 15535 15698 15698 15737 15737 15754 15754 15839 15839 15970 15970 16044 16044 16099 16099 16170 16170 16276 16276 16343 16343 16397 16397 16408 16408 16985 16985 17440 17440 17601 17601 18002 18002 18401 18401 18498 18498 18551 18551 18864 18864 19522 19522 19620 19620 20854 20854 22957 22957
Vertex 18328 connected to: 10164 10164 10665 10665 11237 11237 13558 13558 14473 14473 15485 15485 15535 15535 15698 15698 15833 15833 15970 15970 16003 16003 16044 16044 16353 16353 17366 17366 17538 17538 17601 17601 18498 18498 19522 19522 19760 19760 20854 20854 22257 22257
Vertex 18401 connected to: 9522 9522 9582 9582 10665 10665 11237 11237 13558 13558 14473 14473 14642 14642 14939 14939 14943 14943 15167 15167 15508 15508 15535 15535 15623 15623 15698 15698 15737 15737 16003 16003 16044 16044 16099 16099 16353 16353 16397 16397 16408 16408 16586 16586 16594 16594 16985 16985 17601 17601 18145 18145 18498 18498 20854 20854 22045 22045
Vertex 18407 connected to: 9522 9522 9582 9582 10164 10164 11237 11237 13558 13558 13886 13886 14473 14473 14943 14943 15535 15535 15623 15623 15698 15698 15737 15737 15833 15833 15970 15970 16044 16044 16099 16099 16170 16170 16276 16276 16343 16343 16408 16408 16486 16486 16586 16586 17538 17538 17601 17601 18498 18498 19522 19522 20854 20854
Vertex 18413 connected to: 9522 9522 9582 9582 13390 13390 13458 13458 13886 13886 14939 14939 14943 14943 15167 15167 15535 15535 15623 15623 15833 15833 15970 15970 16003 16003 16044 16044 16099 16099 16170 16170 16343 16343 16397 16397 16586 16586 16985 16985 17538 17538 17601 17601 18498 18498 19029 19029 19101 19101 19522 19522 19760 19760 20854 20854
Vertex 18498 connected to: 9522 9522 9582 9582 10164 10164 10665 10665 11201 11201 13390 13390 13558 13558 13575 13575 13886 13886 14473 14473 14518 14518 14546 14546 15282 15282 15455 15455 15485 15485 15508 15508 15535 15535 15618 15618 15623 15623 15698 15698 15699 15699 15737 15737 15754 15754 15839 15839 16003 16003 16044 16044 16099 16099 16152 16152 16170 16170 16276 16276 16343 16343 16353 16353 16397 16397 16408 16408 16480 16480 16486 16486 16594 16594 16768 16768 16985 16985 17196 17196 17366 17366 17538 17538 17601 17601 17830 17830 18002 18002 18064 18064 18145 18145 18328 18328 18401 18401 18407 18407 18413 18413 18551 18551 18807 18807 18864 18864 19029 19029 19522 19522 19620 19620 19760 19760 20854 20854 22045 22045 25296 25296
Vertex 18551 connected to: 9522 9522 9582 9582 10164 10164 13558 13558 13886 13886 14943 14943 15508 15508 15535 15535 15623 15623 16044 16044 16099 16099 16170 16170 16343 16343 16353 16353 16408 16408 16586 16586 17196 17196 17538 17538 17601 17601 18145 18145 18498 18498 19522 19522 19620 19620 20854 20854
Vertex 18807 connected to: 9522 9522 13390 13390 13458 13458 15167 15167 15508 15508 15623 15623 15833 15833 15970 15970 16044 16044 16170 16170 16985 16985 17366 17366 17538 17538 17601 17601 18498 18498 19029 19029 19760 19760 20854 20854 22257 22257 24750 24750
Vertex 18864 connected to: 9582 9582 10665 10665 11237 11237 11301 11301 13390 13390 13458 13458 13886 13886 14473 14473 14939 14939 14943 14943 15167 15167 15508 15508 15535 15535 15623 15623 15833 15833 16044 16044 16099 16099 16170 16170 16343 16343 16397 16397 16586 16586 16985 16985 17538 17538 17601 17601 18145 18145 18498 18498 19101 19101 19522 19522 19760 19760
Vertex 19029 connected to: 9903 9903 10164 10164 11201 11201 11301 11301 13390 13390 13458 13458 14511 14511 14642 14642 14927 14927 15282 15282 15455 15455 15485 15485 15535 15535 15618 15618 15623 15623 15698 15698 15699 15699 15737 15737 15754 15754 15833 15833 15839 15839 15970 15970 16003 16003 16044 16044 16152 16152 16170 16170 16276 16276 16353 16353 16408 16408 16480 16480 16486 16486 16594 16594 16768 16768 16985 16985 17366 17366 17538 17538 17830 17830 18064 18064 18134 18134 18144 18144 18413 18413 18498 18498 18807 18807 19493 19493 19760 19760 22257 22257 22957 22957 24750 24750 25145 25145 25621 25621
Vertex 19101 connected to: 9582 9582 9903 9903 14518 14518 14592 14592 14927 14927 15282 15282 15455 15455 15485 15485 15535 15535 15623 15623 15698 15698 15699 15699 15737 15737 15754 15754 15839 15839 16003 16003 16044 16044 16099 16099 16353 16353 16408 16408 16586 16586 16768 16768 17196 17196 17366 17366 17440 17440 18064 18064 18413 18413 18864 18864 19493 19493 19760 19760 22045 22045
Vertex 19493 connected to: 13390 13390 14473 14473 14939 14939 15508 15508 15535 15535 15623 15623 15754 15754 16003 16003 16170 16170 16353 16353 16397 16397 16408 16408 16985 16985 17196 17196 17538 17538 17830 17830 19029 19029 19101 19101 19620 19620 19760 19760 20854 20854
Vertex 19522 connected to: 9522 9522 9903 9903 10665 10665 11201 11201 11237 11237 11861 11861 13390 13390 13558 13558 14473 14473 14511 14511 14927 14927 14939 14939 15167 15167 15282 15282 15485 15485 15508 15508 15535 15535 15698 15698 15699 15699 15737 15737 15754 15754 16003 16003 16099 16099 16170 16170 16353 16353 16397 16397 16408 16408 16985 16985 17196 17196 17366 17366 17440 17440 18064 18064 18145 18145 18328 18328 18407 18407 18413 18413 18498 18498 18551 18551 18864 18864 19620 19620 19760 19760 20854 20854 22045 22045 22957 22957
Vertex 19620 connected to: 9522 9522 9582 9582 9903 9903 10665 10665 11237 11237 13558 13558 13575 13575 14473 14473 14592 14592 14939 14939 15167 15167 15455 15455 15485 15485 15508 15508 15535 15535 15699 15699 15737 15737 15754 15754 15833 15833 16099 16099 16170 16170 16397 16397 16408 16408 16586 16586 16985 16985 17366 17366 17601 17601 18002 18002 18145 18145 18498 18498 18551 18551 19493 19493 19522 19522 22257 22257
Vertex 19760 connected to: 9903 9903 10164 10164 11237 11237 11301 11301 13390 13390 14511 14511 14518 14518 14592 14592 14927 14927 14939 14939 15282 15282 15455 15455 15485 15485 15623 15623 15698 15698 15699 15699 15737 15737 15754 15754 15833 15833 15839 15839 15970 15970 16003 16003 16044 16044 16152 16152 16276 16276 16353 16353 16397 16397 16480 16480 16486 16486 16586 16586 16768 16768 17196 17196 17366 17366 17440 17440 17538 17538 17830 17830 18064 18064 18134 18134 18144 18144 18328 18328 18413 18413 18498 18498 18807 18807 18864 18864 19029 19029 19101 19101 19493 19493 19522 19522 20854 20854 22045 22045 22957 22957 24750 24750 25145 25145 25296 25296 25621 25621
Vertex 20854 connected to: 9522 9522 9582 9582 9903 9903 10164 10164 10282 10282 11201 11201 11861 11861 13390 13390 13558 13558 13575 13575 13886 13886 14473 14473 14518 14518 14546 14546 14592 14592 14642 14642 14927 14927 14943 14943 15167 15167 15282 15282 15455 15455 15485 15485 15535 15535 15618 15618 15623 15623 15698 15698 15699 15699 15737 15737 15754 15754 15839 15839 16003 16003 16044 16044 16099 16099 16152 16152 16170 16170 16276 16276 16353 16353 16408 16408 16480 16480 16486 16486 16586 16586 16594 16594 16768 16768 16985 16985 17366 17366 17440 17440 17538 17538 17601 17601 17830 17830 18002 18002 18064 18064 18134 18134 18145 18145 18328 18328 18401 18401 18407 18407 18413 18413 18498 18498 18551 18551 18807 18807 19493 19493 19522 19522 19760 19760 22045 22045 22957 22957
Vertex 22045 connected to: 9582 9582 10665 10665 11237 11237 13390 13390 14642 14642 14943 14943 15535 15535 15623 15623 15737 15737 16044 16044 16343 16343 16408 16408 16586 16586 16985 16985 17538 17538 17601 17601 18401 18401 18498 18498 19101 19101 19522 19522 19760 19760 20854 20854
Vertex 22257 connected to: 9522 9522 10665 10665 11237 11237 11861 11861 13390 13390 13458 13458 14473 14473 14518 14518 14939 14939 14943 14943 15167 15167 15485 15485 15508 15508 15535 15535 15754 15754 15833 15833 15839 15839 15970 15970 16099 16099 16152 16152 16170 16170 16480 16480 16586 16586 16985 16985 17196 17196 17366 17366 18002 18002 18144 18144 18328 18328 18807 18807 19029 19029 19620 19620 25145 25145 25296 25296 25621 25621
Vertex 22957 connected to: 9582 9582 10665 10665 11237 11237 13390 13390 14473 14473 15508 15508 15535 15535 15623 15623 16044 16044 16343 16343 16353 16353 16397 16397 17538 17538 17601 17601 18145 18145 19029 19029 19522 19522 19760 19760 20854 20854
Vertex 24750 connected to: 10665 10665 14511 14511 14642 14642 14927 14927 15282 15282 15455 15455 15485 15485 15698 15698 15970 15970 16003 16003 16152 16152 16170 16170 16276 16276 16353 16353 16480 16480 16486 16486 16594 16594 16768 16768 18807 18807 19029 19029 19760 19760 25145 25145
Vertex 25145 connected to: 9522 9522 13390 13390 14473 14473 14939 14939 15167 15167 15698 15698 15754 15754 15833 15833 15839 15839 15970 15970 16003 16003 16044 16044 16170 16170 16586 16586 16985 16985 17196 17196 17538 17538 18002 18002 19029 19029 19760 19760 22257 22257 24750 24750
Vertex 25296 connected to: 9522 9522 10665 10665 13390 13390 13458 13458 14473 14473 14939 14939 14943 14943 15535 15535 15698 15698 15737 15737 15833 15833 15970 15970 16044 16044 16099 16099 16170 16170 16353 16353 16408 16408 16586 16586 16985 16985 18498 18498 19760 19760 22257 22257
Vertex 25621 connected to: 9522 9522 11237 11237 13390 13390 13458 13458 14939 14939 14943 14943 15167 15167 15508 15508 15535 15535 15698 15698 15833 15833 15970 15970 16003 16003 16044 16044 16099 16099 16170 16170 16353 16353 16397 16397 16586 16586 16985 16985 19029 19029 19760 19760 22257 22257
Number of edges: 3086
Number of k-cliques in dense subgraph: 1543
Average k-cliques per vertex (rho): 17.5341
Execution time: 3817135 milliseconds
```
