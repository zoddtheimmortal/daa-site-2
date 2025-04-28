# EXACT Algorithm

Yixiang Fang†⋆, Kaiqiang Yu‡, Reynold Cheng‡, Laks V.S. Lakshmanan§, Xuemin Lin†⋆

The following C++ implementations are based on [Efficient Algorithms for Densest Subgraph Discovery](https://www.vldb.org/pvldb/vol12/p1719-fang.pdf).

## Pseudo Code

```bash
Algorithm 1: The algorithm: Exact.

Input: G(V, E), Ψ(VΨ, EΨ)
Output: The CDS D(VD, ED)

 initialize l ← 0, u ← max degG(v, Ψ) for v ∈ V
 initialize Λ ← all the instances of (h-1)-clique in G, D ← ∅
 while u - l ≥ 1 / (n(n-1)) do
     α ← (l + u) / 2
     VF ← {s} ∪ V ∪ U ∪ {t}   // build a flow network
     for each vertex v ∈ V do
         add an edge s→v with capacity degG(v, Ψ)
         add an edge v→t with capacity α|VΨ|
     for each (h-1)-clique ψ ∈ Λ do
         for each vertex v ∈ ψ do
             add an edge ψ→v with capacity +∞
     for each (h-1)-clique ψ ∈ Λ do
         for each vertex v ∈ V do
             if ψ and v form an h-clique then
                 add an edge v→ψ with capacity 1
     find minimum s-t cut (S, T) from the flow network F(VF, EF)
     if S = {s} then u ← α
     else l ← α, D ← the subgraph induced by S \ {s}
 return D


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
        outfile << "[";
        for (size_t i = 0; i < clique.size(); ++i) {
            temp.push_back(clique[i]);
            outfile << clique[i];
            if (i != clique.size() - 1) {
                outfile << ", ";
            }
        }
        k_size_cliques.push_back(temp);
        outfile << "]\n";
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

    outfile << "\nResidual Graph after Max Flow:" << endl;
    for (int u = 0; u < V; u++) {
        if (!residualGraph[u].empty()) {
            outfile << u << " -> ";
            for (const auto& edge : residualGraph[u]) {
                outfile << "(" << edge.first << "," << edge.second << ") ";
            }
            outfile << endl;
        }
    }

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

void findDensestSubgraph(Graph& G, vector<vector<int>> kSizeCliques, int k, int n, ofstream& outfile) {
    double l = 0;
    vector<int> clique_degree(n, 0);
    find_clique_degrees(kSizeCliques, clique_degree);

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
        outfile << "Density: " << density << endl;

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
            outfile << "Average k-cliques per vertex: " << avg_cliques << endl;
        }
    }
}

int main() {

    vector<string> input_files = {
        "1ca_hepth.txt",
        "1as_caida_processed.txt",

    };

    for (const string& input_file : input_files) {
        for (int k = 2; k <= 3; k++) {
            if(k==2 && input_file=="1ca_hepth.txt") continue;

            string output_file = "11final_output_" + input_file + "_" + to_string(k) + ".txt";

            vector<pair<int, int>> edges;
            int n;

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

#### Yeast Dataset

For k=2

```bash
Current alpha: 28
Max Flow: 3896
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 28
Current gap: 28

Current alpha: 14
Max Flow: 3896
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 14
Current gap: 14

Current alpha: 7
Max Flow: 3896
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 7
Current gap: 7

Current alpha: 3.5
Max Flow: 3896
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 3.5
Current gap: 3.5

Current alpha: 1.75
Max Flow: 3484
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 1.75
Current gap: 1.75

Current alpha: 2.625
Max Flow: 3893
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 2.625
Current gap: 0.875

Current alpha: 3.0625
Max Flow: 3896
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 3.0625
Current gap: 0.4375

Current alpha: 2.84375
Max Flow: 3893
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 2.84375
Current gap: 0.21875

Current alpha: 2.95312
Max Flow: 3893
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 2.95312
Current gap: 0.109375

Current alpha: 3.00781
Max Flow: 3896
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 3.00781
Current gap: 0.0546875

Current alpha: 2.98047
Max Flow: 3893
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 2.98047
Current gap: 0.0273438

Current alpha: 2.99414
Max Flow: 3893
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 2.99414
Current gap: 0.0136719

Current alpha: 3.00098
Max Flow: 3896
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 3.00098
Current gap: 0.00683594

Current alpha: 2.99756
Max Flow: 3893
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 2.99756
Current gap: 0.00341797

Current alpha: 2.99927
Max Flow: 3893
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 2.99927
Current gap: 0.00170898

Current alpha: 3.00012
Max Flow: 3896
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 3.00012
Current gap: 0.000854492

Current alpha: 2.99969
Max Flow: 3893
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 2.99969
Current gap: 0.000427246

Current alpha: 2.99991
Max Flow: 3893
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 2.99991
Current gap: 0.000213623

Current alpha: 3.00002
Max Flow: 3896
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 3.00002
Current gap: 0.000106812

Current alpha: 2.99996
Max Flow: 3893
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 2.99996
Current gap: 5.34058e-005

Current alpha: 2.99999
Max Flow: 3893
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 2.99999
Current gap: 2.67029e-005

Current alpha: 3
Max Flow: 3896
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 3
Current gap: 1.33514e-005

Current alpha: 3
Max Flow: 3893
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 3
Current gap: 6.67572e-006

Current alpha: 3
Max Flow: 3893
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 3
Current gap: 3.33786e-006

Current alpha: 3
Max Flow: 3896
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 3
Current gap: 1.66893e-006

Current alpha: 3
Max Flow: 3893
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 3
Current gap: 8.34465e-007

Current alpha: 3
Max Flow: 3893
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 3
Current gap: 4.17233e-007


Dense subgraph found with 7 vertices:
Vertices in dense subgraph: 805 912 960 1075 1091 1099 1100
Edges in dense subgraph:
Vertex 805 connected to: 912 960 1091 1099 1100
Vertex 912 connected to: 805 960 1075 1091 1099 1100
Vertex 960 connected to: 805 912 1075 1091 1099 1100
Vertex 1075 connected to: 912 960 1099 1100
Vertex 1091 connected to: 805 912 960 1099 1100
Vertex 1099 connected to: 805 912 960 1075 1091 1100
Vertex 1100 connected to: 805 912 960 1075 1091 1099
Number of edges: 19
Density: 0.904762
Number of k-cliques in dense subgraph: 19
Average k-cliques per vertex: 2.71429
Execution time: 50157 milliseconds

```

For k=3

```bash

Current alpha: 9
Max Flow: 618
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 9
Current gap: 9

Current alpha: 4.5
Max Flow: 618
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 4.5
Current gap: 4.5

Current alpha: 2.25
Max Flow: 579
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 2.25
Current gap: 2.25

Current alpha: 3.375
Max Flow: 610
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 3.375
Current gap: 1.125

Current alpha: 3.9375
Max Flow: 617
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 3.9375
Current gap: 0.5625

Current alpha: 4.21875
Max Flow: 618
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 4.21875
Current gap: 0.28125

Current alpha: 4.07812
Max Flow: 618
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 4.07812
Current gap: 0.140625

Current alpha: 4.00781
Max Flow: 618
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 4.00781
Current gap: 0.0703125

Current alpha: 3.97266
Max Flow: 617
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 3.97266
Current gap: 0.0351562

Current alpha: 3.99023
Max Flow: 617
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 3.99023
Current gap: 0.0175781

Current alpha: 3.99902
Max Flow: 617
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 3.99902
Current gap: 0.00878906

Current alpha: 4.00342
Max Flow: 618
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 4.00342
Current gap: 0.00439453

Current alpha: 4.00122
Max Flow: 618
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 4.00122
Current gap: 0.00219727

Current alpha: 4.00012
Max Flow: 618
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 4.00012
Current gap: 0.00109863

Current alpha: 3.99957
Max Flow: 617
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 3.99957
Current gap: 0.000549316

Current alpha: 3.99985
Max Flow: 617
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 3.99985
Current gap: 0.000274658

Current alpha: 3.99998
Max Flow: 617
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 3.99998
Current gap: 0.000137329

Current alpha: 4.00005
Max Flow: 618
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 4.00005
Current gap: 6.86646e-005

Current alpha: 4.00002
Max Flow: 618
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 4.00002
Current gap: 3.43323e-005

Current alpha: 4
Max Flow: 618
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 4
Current gap: 1.71661e-005

Current alpha: 3.99999
Max Flow: 617
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 3.99999
Current gap: 8.58307e-006

Current alpha: 4
Max Flow: 617
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 4
Current gap: 4.29153e-006

Current alpha: 4
Max Flow: 617
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 4
Current gap: 2.14577e-006

Current alpha: 4
Max Flow: 618
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 4
Current gap: 1.07288e-006

Current alpha: 4
Max Flow: 618
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 4
Current gap: 5.36442e-007

Current alpha: 4
Max Flow: 618
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 4
Current gap: 2.68221e-007


Dense subgraph found with 7 vertices:
Vertices in dense subgraph: 805 912 960 1075 1091 1099 1100
Edges in dense subgraph:
Vertex 805 connected to: 912 960 1091 1099 1100
Vertex 912 connected to: 805 960 1075 1091 1099 1100
Vertex 960 connected to: 805 912 1075 1091 1099 1100
Vertex 1075 connected to: 912 960 1099 1100
Vertex 1091 connected to: 805 912 960 1099 1100
Vertex 1099 connected to: 805 912 960 1075 1091 1100
Vertex 1100 connected to: 805 912 960 1075 1091 1099
Number of edges: 19
Density: 0.904762
Number of k-cliques in dense subgraph: 26
Average k-cliques per vertex: 3.71429
Execution time: 60837 milliseconds


```

k=4
```bash

Current alpha: 6.5
Max Flow: 156
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 6.5
Current gap: 6.5

Current alpha: 3.25
Max Flow: 156
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 3.25
Current gap: 3.25

Current alpha: 1.625
Max Flow: 122
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 1.625
Current gap: 1.625

Current alpha: 2.4375
Max Flow: 143
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 2.4375
Current gap: 0.8125

Current alpha: 2.84375
Max Flow: 156
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 2.84375
Current gap: 0.40625

Current alpha: 2.64062
Max Flow: 150
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 2.64062
Current gap: 0.203125

Current alpha: 2.74219
Max Flow: 150
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 2.74219
Current gap: 0.101562

Current alpha: 2.79297
Max Flow: 156
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 2.79297
Current gap: 0.0507812

Current alpha: 2.76758
Max Flow: 156
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 2.76758
Current gap: 0.0253906

Current alpha: 2.75488
Max Flow: 156
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 2.75488
Current gap: 0.0126953

Current alpha: 2.74854
Max Flow: 150
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 2.74854
Current gap: 0.00634766

Current alpha: 2.75171
Max Flow: 156
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 2.75171
Current gap: 0.00317383

Current alpha: 2.75012
Max Flow: 156
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 2.75012
Current gap: 0.00158691

Current alpha: 2.74933
Max Flow: 150
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 2.74933
Current gap: 0.000793457

Current alpha: 2.74973
Max Flow: 150
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 2.74973
Current gap: 0.000396729

Current alpha: 2.74992
Max Flow: 150
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 2.74992
Current gap: 0.000198364

Current alpha: 2.75002
Max Flow: 156
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 2.75002
Current gap: 9.91821e-005

Current alpha: 2.74997
Max Flow: 150
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 2.74997
Current gap: 4.95911e-005

Current alpha: 2.75
Max Flow: 150
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 2.75
Current gap: 2.47955e-005

Current alpha: 2.75001
Max Flow: 156
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 2.75001
Current gap: 1.23978e-005

Current alpha: 2.75
Max Flow: 156
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 2.75
Current gap: 6.19888e-006

Current alpha: 2.75
Max Flow: 156
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 2.75
Current gap: 3.09944e-006

Current alpha: 2.75
Max Flow: 150
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 2.75
Current gap: 1.54972e-006

Current alpha: 2.75
Max Flow: 156
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 2.75
Current gap: 7.7486e-007

Current alpha: 2.75
Max Flow: 156
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 2.75
Current gap: 3.8743e-007


Dense subgraph found with 7 vertices:
Vertices in dense subgraph: 805 912 960 1075 1091 1099 1100
Edges in dense subgraph:
Vertex 805 connected to: 912 960 1091 1099 1100
Vertex 912 connected to: 805 960 1075 1091 1099 1100
Vertex 960 connected to: 805 912 1075 1091 1099 1100
Vertex 1075 connected to: 912 960 1099 1100
Vertex 1091 connected to: 805 912 960 1099 1100
Vertex 1099 connected to: 805 912 960 1075 1091 1100
Vertex 1100 connected to: 805 912 960 1075 1091 1099
Number of edges: 19
Density: 0.904762
Number of k-cliques in dense subgraph: 19
Average k-cliques per vertex: 2.71429
Execution time: 7435 milliseconds


```

k=5

```bash

Current alpha: 3
Max Flow: 40
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 3
Current gap: 3

Current alpha: 1.5
Max Flow: 40
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 1.5
Current gap: 1.5

Current alpha: 0.75
Max Flow: 26
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 0.75
Current gap: 0.75

Current alpha: 1.125
Max Flow: 40
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 1.125
Current gap: 0.375

Current alpha: 0.9375
Max Flow: 33
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 0.9375
Current gap: 0.1875

Current alpha: 1.03125
Max Flow: 40
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 1.03125
Current gap: 0.09375

Current alpha: 0.984375
Max Flow: 33
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 0.984375
Current gap: 0.046875

Current alpha: 1.00781
Max Flow: 40
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 1.00781
Current gap: 0.0234375

Current alpha: 0.996094
Max Flow: 33
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 0.996094
Current gap: 0.0117188

Current alpha: 1.00195
Max Flow: 40
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 1.00195
Current gap: 0.00585938

Current alpha: 0.999023
Max Flow: 33
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 0.999023
Current gap: 0.00292969

Current alpha: 1.00049
Max Flow: 40
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 1.00049
Current gap: 0.00146484

Current alpha: 0.999756
Max Flow: 33
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 0.999756
Current gap: 0.000732422

Current alpha: 1.00012
Max Flow: 40
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 1.00012
Current gap: 0.000366211

Current alpha: 0.999939
Max Flow: 33
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 0.999939
Current gap: 0.000183105

Current alpha: 1.00003
Max Flow: 40
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 1.00003
Current gap: 9.15527e-005

Current alpha: 0.999985
Max Flow: 33
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 0.999985
Current gap: 4.57764e-005

Current alpha: 1.00001
Max Flow: 40
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 1.00001
Current gap: 2.28882e-005

Current alpha: 0.999996
Max Flow: 33
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 0.999996
Current gap: 1.14441e-005

Current alpha: 1
Max Flow: 40
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 1
Current gap: 5.72205e-006

Current alpha: 0.999999
Max Flow: 33
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 0.999999
Current gap: 2.86102e-006

Current alpha: 1
Max Flow: 40
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 1
Current gap: 1.43051e-006

Current alpha: 1
Max Flow: 33
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 1
Current gap: 7.15256e-007

Current alpha: 1
Max Flow: 40
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 1
Current gap: 3.57628e-007


Dense subgraph found with 7 vertices:
Vertices in dense subgraph: 805 912 960 1075 1091 1099 1100
Edges in dense subgraph:
Vertex 805 connected to: 912 960 1091 1099 1100
Vertex 912 connected to: 805 960 1075 1091 1099 1100
Vertex 960 connected to: 805 912 1075 1091 1099 1100
Vertex 1075 connected to: 912 960 1099 1100
Vertex 1091 connected to: 805 912 960 1099 1100
Vertex 1099 connected to: 805 912 960 1075 1091 1100
Vertex 1100 connected to: 805 912 960 1075 1091 1099
Number of edges: 19
Density: 0.904762
Number of k-cliques in dense subgraph: 7
Average k-cliques per vertex: 1
Execution time: 2040 milliseconds

```

#### Netscience Dataset

For k=2

```bash
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
Current gap: 6.48499e-005

Current alpha: 9.50002
Max Flow: 5484
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 9.50002
Current gap: 3.24249e-005

Current alpha: 9.5
Max Flow: 5484
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 9.5
Current gap: 1.62125e-005

Current alpha: 9.49999
Max Flow: 5464
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 9.49999
Current gap: 8.10623e-006

Current alpha: 9.5
Max Flow: 5464
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 9.5
Current gap: 4.05312e-006

Current alpha: 9.5
Max Flow: 5464
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 9.5
Current gap: 2.02656e-006

Current alpha: 9.5
Max Flow: 5464
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 9.5
Current gap: 1.01328e-006

Current alpha: 9.5
Max Flow: 5484
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 9.5
Current gap: 5.06639e-007

Current alpha: 9.5
Max Flow: 5484
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 9.5
Current gap: 2.5332e-007


Dense subgraph found with 20 vertices:
Vertices in dense subgraph: 550 551 552 553 554 555 556 557 558 559 560 561 562 563 564 565 566 567 568 569
Edges in dense subgraph:
Vertex 550 connected to: 551 552 553 554 555 556 557 558 559 560 561 562 563 564 565 566 567 568 569
Vertex 551 connected to: 550 552 553 554 555 556 557 558 559 560 561 562 563 564 565 566 567 568 569
Vertex 552 connected to: 550 551 553 554 555 556 557 558 559 560 561 562 563 564 565 566 567 568 569
Vertex 553 connected to: 550 551 552 554 555 556 557 558 559 560 561 562 563 564 565 566 567 568 569
Vertex 554 connected to: 550 551 552 553 555 556 557 558 559 560 561 562 563 564 565 566 567 568 569
Vertex 555 connected to: 550 551 552 553 554 556 557 558 559 560 561 562 563 564 565 566 567 568 569
Vertex 556 connected to: 550 551 552 553 554 555 557 558 559 560 561 562 563 564 565 566 567 568 569
Vertex 557 connected to: 550 551 552 553 554 555 556 558 559 560 561 562 563 564 565 566 567 568 569
Vertex 558 connected to: 550 551 552 553 554 555 556 557 559 560 561 562 563 564 565 566 567 568 569
Vertex 559 connected to: 550 551 552 553 554 555 556 557 558 560 561 562 563 564 565 566 567 568 569
Vertex 560 connected to: 550 551 552 553 554 555 556 557 558 559 561 562 563 564 565 566 567 568 569
Vertex 561 connected to: 550 551 552 553 554 555 556 557 558 559 560 562 563 564 565 566 567 568 569
Vertex 562 connected to: 550 551 552 553 554 555 556 557 558 559 560 561 563 564 565 566 567 568 569
Vertex 563 connected to: 550 551 552 553 554 555 556 557 558 559 560 561 562 564 565 566 567 568 569
Vertex 564 connected to: 550 551 552 553 554 555 556 557 558 559 560 561 562 563 565 566 567 568 569
Vertex 565 connected to: 550 551 552 553 554 555 556 557 558 559 560 561 562 563 564 566 567 568 569
Vertex 566 connected to: 550 551 552 553 554 555 556 557 558 559 560 561 562 563 564 565 567 568 569
Vertex 567 connected to: 550 551 552 553 554 555 556 557 558 559 560 561 562 563 564 565 566 568 569
Vertex 568 connected to: 550 551 552 553 554 555 556 557 558 559 560 561 562 563 564 565 566 567 569
Vertex 569 connected to: 550 551 552 553 554 555 556 557 558 559 560 561 562 563 564 565 566 567 568
Number of edges: 190
Number of k-cliques in dense subgraph: 190
rho: 9.5
Execution time: 38098 milliseconds



```

For k=3

```bash
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
Current gap: 8.24928e-005

Current alpha: 57
Max Flow: 11292
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 57
Current gap: 4.12464e-005

Current alpha: 57
Max Flow: 11292
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 57
Current gap: 2.06232e-005

Current alpha: 57
Max Flow: 11272
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 57
Current gap: 1.03116e-005

Current alpha: 57
Max Flow: 11292
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 57
Current gap: 5.1558e-006

Current alpha: 57
Max Flow: 11292
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 57
Current gap: 2.5779e-006

Current alpha: 57
Max Flow: 11272
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 57
Current gap: 1.28895e-006

Current alpha: 57
Max Flow: 11272
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 57
Current gap: 6.44475e-007

Current alpha: 57
Max Flow: 11272
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 57
Current gap: 3.22238e-007


Dense subgraph found with 20 vertices:
Vertices in dense subgraph: 550 551 552 553 554 555 556 557 558 559 560 561 562 563 564 565 566 567 568 569
Edges in dense subgraph:
Vertex 550 connected to: 551 552 553 554 555 556 557 558 559 560 561 562 563 564 565 566 567 568 569
Vertex 551 connected to: 550 552 553 554 555 556 557 558 559 560 561 562 563 564 565 566 567 568 569
Vertex 552 connected to: 550 551 553 554 555 556 557 558 559 560 561 562 563 564 565 566 567 568 569
Vertex 553 connected to: 550 551 552 554 555 556 557 558 559 560 561 562 563 564 565 566 567 568 569
Vertex 554 connected to: 550 551 552 553 555 556 557 558 559 560 561 562 563 564 565 566 567 568 569
Vertex 555 connected to: 550 551 552 553 554 556 557 558 559 560 561 562 563 564 565 566 567 568 569
Vertex 556 connected to: 550 551 552 553 554 555 557 558 559 560 561 562 563 564 565 566 567 568 569
Vertex 557 connected to: 550 551 552 553 554 555 556 558 559 560 561 562 563 564 565 566 567 568 569
Vertex 558 connected to: 550 551 552 553 554 555 556 557 559 560 561 562 563 564 565 566 567 568 569
Vertex 559 connected to: 550 551 552 553 554 555 556 557 558 560 561 562 563 564 565 566 567 568 569
Vertex 560 connected to: 550 551 552 553 554 555 556 557 558 559 561 562 563 564 565 566 567 568 569
Vertex 561 connected to: 550 551 552 553 554 555 556 557 558 559 560 562 563 564 565 566 567 568 569
Vertex 562 connected to: 550 551 552 553 554 555 556 557 558 559 560 561 563 564 565 566 567 568 569
Vertex 563 connected to: 550 551 552 553 554 555 556 557 558 559 560 561 562 564 565 566 567 568 569
Vertex 564 connected to: 550 551 552 553 554 555 556 557 558 559 560 561 562 563 565 566 567 568 569
Vertex 565 connected to: 550 551 552 553 554 555 556 557 558 559 560 561 562 563 564 566 567 568 569
Vertex 566 connected to: 550 551 552 553 554 555 556 557 558 559 560 561 562 563 564 565 567 568 569
Vertex 567 connected to: 550 551 552 553 554 555 556 557 558 559 560 561 562 563 564 565 566 568 569
Vertex 568 connected to: 550 551 552 553 554 555 556 557 558 559 560 561 562 563 564 565 566 567 569
Vertex 569 connected to: 550 551 552 553 554 555 556 557 558 559 560 561 562 563 564 565 566 567 568
Number of edges: 190
Number of k-cliques in dense subgraph: 1140
rho: 57
Execution time: 87882 milliseconds



```

For k=4

```bash
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
Current gap: 5.78165e-005

Current alpha: 242.25
Max Flow: 28616
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 242.25
Current gap: 2.89083e-005

Current alpha: 242.25
Max Flow: 28616
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 242.25
Current gap: 1.44541e-005

Current alpha: 242.25
Max Flow: 28616
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 242.25
Current gap: 7.22706e-006

Current alpha: 242.25
Max Flow: 28616
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 242.25
Current gap: 3.61353e-006

Current alpha: 242.25
Max Flow: 28616
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 242.25
Current gap: 1.80677e-006

Current alpha: 242.25
Max Flow: 28636
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 242.25
Current gap: 9.03383e-007

Current alpha: 242.25
Max Flow: 28636
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 242.25
Current gap: 4.51691e-007


Dense subgraph found with 20 vertices:
Vertices in dense subgraph: 550 551 552 553 554 555 556 557 558 559 560 561 562 563 564 565 566 567 568 569
Edges in dense subgraph:
Vertex 550 connected to: 551 552 553 554 555 556 557 558 559 560 561 562 563 564 565 566 567 568 569
Vertex 551 connected to: 550 552 553 554 555 556 557 558 559 560 561 562 563 564 565 566 567 568 569
Vertex 552 connected to: 550 551 553 554 555 556 557 558 559 560 561 562 563 564 565 566 567 568 569
Vertex 553 connected to: 550 551 552 554 555 556 557 558 559 560 561 562 563 564 565 566 567 568 569
Vertex 554 connected to: 550 551 552 553 555 556 557 558 559 560 561 562 563 564 565 566 567 568 569
Vertex 555 connected to: 550 551 552 553 554 556 557 558 559 560 561 562 563 564 565 566 567 568 569
Vertex 556 connected to: 550 551 552 553 554 555 557 558 559 560 561 562 563 564 565 566 567 568 569
Vertex 557 connected to: 550 551 552 553 554 555 556 558 559 560 561 562 563 564 565 566 567 568 569
Vertex 558 connected to: 550 551 552 553 554 555 556 557 559 560 561 562 563 564 565 566 567 568 569
Vertex 559 connected to: 550 551 552 553 554 555 556 557 558 560 561 562 563 564 565 566 567 568 569
Vertex 560 connected to: 550 551 552 553 554 555 556 557 558 559 561 562 563 564 565 566 567 568 569
Vertex 561 connected to: 550 551 552 553 554 555 556 557 558 559 560 562 563 564 565 566 567 568 569
Vertex 562 connected to: 550 551 552 553 554 555 556 557 558 559 560 561 563 564 565 566 567 568 569
Vertex 563 connected to: 550 551 552 553 554 555 556 557 558 559 560 561 562 564 565 566 567 568 569
Vertex 564 connected to: 550 551 552 553 554 555 556 557 558 559 560 561 562 563 565 566 567 568 569
Vertex 565 connected to: 550 551 552 553 554 555 556 557 558 559 560 561 562 563 564 566 567 568 569
Vertex 566 connected to: 550 551 552 553 554 555 556 557 558 559 560 561 562 563 564 565 567 568 569
Vertex 567 connected to: 550 551 552 553 554 555 556 557 558 559 560 561 562 563 564 565 566 568 569
Vertex 568 connected to: 550 551 552 553 554 555 556 557 558 559 560 561 562 563 564 565 566 567 569
Vertex 569 connected to: 550 551 552 553 554 555 556 557 558 559 560 561 562 563 564 565 566 567 568
Number of edges: 190
Number of k-cliques in dense subgraph: 4845
rho: 242.25
Execution time: 144461 milliseconds

```

#### As733(latest) Dataset

For k=2

```bash

Current alpha: 729
Max Flow: 25144
Min cut (nodes in S side): 0
Current gap: 729

Current alpha: 364.5
Max Flow: 25144
Min cut (nodes in S side): 0
Current gap: 364.5

Current alpha: 182.25
Max Flow: 25144
Min cut (nodes in S side): 0
Current gap: 182.25

Current alpha: 91.125
Max Flow: 25144
Min cut (nodes in S side): 0
Current gap: 91.125

Current alpha: 45.5625
Max Flow: 25144
Min cut (nodes in S side): 0
Current gap: 45.5625

Current alpha: 22.7812
Max Flow: 25144
Min cut (nodes in S side): 0
Current gap: 22.7812

Current alpha: 11.3906
Max Flow: 25144
Min cut (nodes in S side): 0
Current gap: 11.3906

Current alpha: 5.69531
Max Flow: 24547
Min cut (nodes in S side): 0 1 2 3 4 5 6 7 8 9 10 11 13 16 17 19 20 21 22 23 25 27 29 32 39 42 46 48 49 50 52 53 57 58 59 60 61 62 63 64 66 68 77 86 87 93 94 101 102 106 127 130 132 134 135 139 149 168 171 172 174 177 181 184 197 198 200 208 218 243 248 249 252 310 329 353 373 418 438 470 476 477 551 552 553 580 601 632 654 663 699 717 772 802 818 857 858 859 867 872 873 875 877 893 900 902 903 904 905 915 921 923 924 927 929 931 933 938 939 941 942 947 996 1033 1041 1046 1080 1113 1118 1120 1121 1312 1482 1673 1815 2182 2188 2195 2205 2219 2593 2959 2989 3255 3363 4335 4388 4472 6475 6476 6477 6478 6479 6480 6481 6482 6483 6484 6485 6487 6490 6491 6493 6494 6495 6496 6497 6499 6501 6503 6506 6513 6516 6520 6522 6523 6524 6526 6527 6531 6532 6533 6534 6535 6536 6537 6538 6540 6542 6551 6560 6561 6567 6568 6575 6576 6580 6601 6604 6606 6608 6609 6613 6623 6642 6645 6646 6648 6651 6655 6658 6671 6672 6674 6682 6692 6717 6722 6723 6726 6784 6803 6827 6847 6892 6912 6944 6950 6951 7025 7026 7027 7054 7075 7106 7128 7137 7173 7191 7246 7276 7292 7331 7332 7333 7341 7346 7347 7349 7351 7367 7374 7376 7377 7378 7379 7389 7395 7397 7398 7401 7403 7405 7407 7412 7413 7415 7416 7421 7470 7507 7515 7520 7554 7587 7592 7594 7595 7786 7956 8147 8289 8656 8662 8669 8679 8693 9067 9433 9463 9729 9837 10809 10862 10946
Current gap: 5.69531

Current alpha: 8.54297
Max Flow: 25114
Min cut (nodes in S side): 0 1 2 3 4 5 6 7 8 9 10 11 17 21 22 23 25 27 29 39 42 48 50 52 57 61 64 86 87 93 94 102 132 134 149 174 181 197 601 632 654 6475 6476 6477 6478 6479 6480 6481 6482 6483 6484 6485 6491 6495 6496 6497 6499 6501 6503 6513 6516 6522 6524 6526 6531 6535 6538 6560 6561 6567 6568 6576 6606 6608 6623 6648 6655 6671 7075 7106 7128
Current gap: 2.84766

Current alpha: 9.9668
Max Flow: 25144
Min cut (nodes in S side): 0
Current gap: 1.42383

Current alpha: 9.25488
Max Flow: 25144
Min cut (nodes in S side): 0
Current gap: 0.711914

Current alpha: 8.89893
Max Flow: 25114
Min cut (nodes in S side): 0 1 2 3 4 5 6 7 8 9 10 11 17 21 22 23 25 27 29 39 42 48 50 52 57 61 64 86 87 93 94 102 132 134 149 174 181 197 601 632 654 6475 6476 6477 6478 6479 6480 6481 6482 6483 6484 6485 6491 6495 6496 6497 6499 6501 6503 6513 6516 6522 6524 6526 6531 6535 6538 6560 6561 6567 6568 6576 6606 6608 6623 6648 6655 6671 7075 7106 7128
Current gap: 0.355957

Current alpha: 9.0769
Max Flow: 25144
Min cut (nodes in S side): 0
Current gap: 0.177979

Current alpha: 8.98792
Max Flow: 25114
Min cut (nodes in S side): 0 1 2 3 4 5 6 7 8 9 10 11 17 21 22 23 25 27 29 39 42 48 50 52 57 61 64 86 87 93 94 102 132 134 149 174 181 197 601 632 654 6475 6476 6477 6478 6479 6480 6481 6482 6483 6484 6485 6491 6495 6496 6497 6499 6501 6503 6513 6516 6522 6524 6526 6531 6535 6538 6560 6561 6567 6568 6576 6606 6608 6623 6648 6655 6671 7075 7106 7128
Current gap: 0.0889893

Current alpha: 9.03241
Max Flow: 25144
Min cut (nodes in S side): 0
Current gap: 0.0444946

Current alpha: 9.01016
Max Flow: 25144
Min cut (nodes in S side): 0
Current gap: 0.0222473

Current alpha: 8.99904
Max Flow: 25114
Min cut (nodes in S side): 0 1 2 3 4 5 6 7 8 9 10 11 17 21 22 23 25 27 29 39 42 48 50 52 57 61 64 86 87 93 94 102 132 134 149 174 181 197 601 632 654 6475 6476 6477 6478 6479 6480 6481 6482 6483 6484 6485 6491 6495 6496 6497 6499 6501 6503 6513 6516 6522 6524 6526 6531 6535 6538 6560 6561 6567 6568 6576 6606 6608 6623 6648 6655 6671 7075 7106 7128
Current gap: 0.0111237

Current alpha: 9.0046
Max Flow: 25144
Min cut (nodes in S side): 0
Current gap: 0.00556183

Current alpha: 9.00182
Max Flow: 25144
Min cut (nodes in S side): 0
Current gap: 0.00278091

Current alpha: 9.00043
Max Flow: 25144
Min cut (nodes in S side): 0
Current gap: 0.00139046

Current alpha: 8.99973
Max Flow: 25114
Min cut (nodes in S side): 0 1 2 3 4 5 6 7 8 9 10 11 17 21 22 23 25 27 29 39 42 48 50 52 57 61 64 86 87 93 94 102 132 134 149 174 181 197 601 632 654 6475 6476 6477 6478 6479 6480 6481 6482 6483 6484 6485 6491 6495 6496 6497 6499 6501 6503 6513 6516 6522 6524 6526 6531 6535 6538 6560 6561 6567 6568 6576 6606 6608 6623 6648 6655 6671 7075 7106 7128
Current gap: 0.000695229

Current alpha: 9.00008
Max Flow: 25144
Min cut (nodes in S side): 0
Current gap: 0.000347614

Current alpha: 8.99991
Max Flow: 25114
Min cut (nodes in S side): 0 1 2 3 4 5 6 7 8 9 10 11 17 21 22 23 25 27 29 39 42 48 50 52 57 61 64 86 87 93 94 102 132 134 149 174 181 197 601 632 654 6475 6476 6477 6478 6479 6480 6481 6482 6483 6484 6485 6491 6495 6496 6497 6499 6501 6503 6513 6516 6522 6524 6526 6531 6535 6538 6560 6561 6567 6568 6576 6606 6608 6623 6648 6655 6671 7075 7106 7128
Current gap: 0.000173807

Current alpha: 8.99999
Max Flow: 25114
Min cut (nodes in S side): 0 1 2 3 4 5 6 7 8 9 10 11 17 21 22 23 25 27 29 39 42 48 50 52 57 61 64 86 87 93 94 102 132 134 149 174 181 197 601 632 654 6475 6476 6477 6478 6479 6480 6481 6482 6483 6484 6485 6491 6495 6496 6497 6499 6501 6503 6513 6516 6522 6524 6526 6531 6535 6538 6560 6561 6567 6568 6576 6606 6608 6623 6648 6655 6671 7075 7106 7128
Current gap: 8.69036e-005

Current alpha: 9.00004
Max Flow: 25144
Min cut (nodes in S side): 0
Current gap: 4.34518e-005

Current alpha: 9.00002
Max Flow: 25144
Min cut (nodes in S side): 0
Current gap: 2.17259e-005

Current alpha: 9.00001
Max Flow: 25144
Min cut (nodes in S side): 0
Current gap: 1.08629e-005

Current alpha: 9
Max Flow: 25144
Min cut (nodes in S side): 0
Current gap: 5.43147e-006

Current alpha: 9
Max Flow: 25114
Min cut (nodes in S side): 0 1 2 3 4 5 6 7 8 9 10 11 17 21 22 23 25 27 29 39 42 48 50 52 57 61 64 86 87 93 94 102 132 134 149 174 181 197 601 632 654 6475 6476 6477 6478 6479 6480 6481 6482 6483 6484 6485 6491 6495 6496 6497 6499 6501 6503 6513 6516 6522 6524 6526 6531 6535 6538 6560 6561 6567 6568 6576 6606 6608 6623 6648 6655 6671 7075 7106 7128
Current gap: 2.71574e-006

Current alpha: 9
Max Flow: 25114
Min cut (nodes in S side): 0 1 2 3 4 5 6 7 8 9 10 11 17 21 22 23 25 27 29 39 42 48 50 52 57 61 64 86 87 93 94 102 132 134 149 174 181 197 601 632 654 6475 6476 6477 6478 6479 6480 6481 6482 6483 6484 6485 6491 6495 6496 6497 6499 6501 6503 6513 6516 6522 6524 6526 6531 6535 6538 6560 6561 6567 6568 6576 6606 6608 6623 6648 6655 6671 7075 7106 7128
Current gap: 1.35787e-006

Current alpha: 9
Max Flow: 25114
Min cut (nodes in S side): 0 1 2 3 4 5 6 7 8 9 10 11 17 21 22 23 25 27 29 39 42 48 50 52 57 61 64 86 87 93 94 102 132 134 149 174 181 197 601 632 654 6475 6476 6477 6478 6479 6480 6481 6482 6483 6484 6485 6491 6495 6496 6497 6499 6501 6503 6513 6516 6522 6524 6526 6531 6535 6538 6560 6561 6567 6568 6576 6606 6608 6623 6648 6655 6671 7075 7106 7128
Current gap: 6.78934e-007

Current alpha: 9
Max Flow: 25114
Min cut (nodes in S side): 0 1 2 3 4 5 6 7 8 9 10 11 17 21 22 23 25 27 29 39 42 48 50 52 57 61 64 86 87 93 94 102 132 134 149 174 181 197 601 632 654 6475 6476 6477 6478 6479 6480 6481 6482 6483 6484 6485 6491 6495 6496 6497 6499 6501 6503 6513 6516 6522 6524 6526 6531 6535 6538 6560 6561 6567 6568 6576 6606 6608 6623 6648 6655 6671 7075 7106 7128
Current gap: 3.39467e-007

Current alpha: 9
Max Flow: 25114
Min cut (nodes in S side): 0 1 2 3 4 5 6 7 8 9 10 11 17 21 22 23 25 27 29 39 42 48 50 52 57 61 64 86 87 93 94 102 132 134 149 174 181 197 601 632 654 6475 6476 6477 6478 6479 6480 6481 6482 6483 6484 6485 6491 6495 6496 6497 6499 6501 6503 6513 6516 6522 6524 6526 6531 6535 6538 6560 6561 6567 6568 6576 6606 6608 6623 6648 6655 6671 7075 7106 7128
Current gap: 1.69734e-007

Current alpha: 9
Max Flow: 25114
Min cut (nodes in S side): 0 1 2 3 4 5 6 7 8 9 10 11 17 21 22 23 25 27 29 39 42 48 50 52 57 61 64 86 87 93 94 102 132 134 149 174 181 197 601 632 654 6475 6476 6477 6478 6479 6480 6481 6482 6483 6484 6485 6491 6495 6496 6497 6499 6501 6503 6513 6516 6522 6524 6526 6531 6535 6538 6560 6561 6567 6568 6576 6606 6608 6623 6648 6655 6671 7075 7106 7128
Current gap: 8.48668e-008

Current alpha: 9
Max Flow: 25144
Min cut (nodes in S side): 0
Current gap: 4.24334e-008

Current alpha: 9
Max Flow: 25144
Min cut (nodes in S side): 0
Current gap: 2.12167e-008


Dense subgraph found with 40 vertices:
Vertices in dense subgraph: 0 1 2 3 4 5 6 7 8 9 10 16 20 21 22 24 26 28 38 41 47 49 51 56 60 63 85 86 92 93 101 131 133 148 173 180 196 600 631 653
Edges in dense subgraph:
Vertex 0 connected to: 1 2 3 4 5 6 7 8 9 10 16 20 21 22 24 26 28 38 41 47 49 51 56 60 63 85 86 92 93 101 131 133 148 173 180 196
Vertex 1 connected to: 0 1 1 2 3 4 5 6 7 8 9 10 16 20 21 22 26 28 38 41 47 49 51 56 60 63 85 86 92 93 101 131 148 173 180 196 600 631
Vertex 2 connected to: 0 1 3 4 5 6 7 8 9 16 20 21 22 24 26 28 38 41 47 49 51 56 60 63 85 86 92 93 101 131 133 148 173 180 196 631 653
Vertex 3 connected to: 0 1 2 4 5 6 7 8 9 16 20 21 22 24 26 28 38 41 47 49 51 56 60 63 85 86 92 101 131 148 173 180 196 600 631 653
Vertex 4 connected to: 0 1 2 3 5 6 8 10 16 20 22 24 26 38 41 47 51 56 60 85 86 92 101 131 133 148 180 196 631 653
Vertex 5 connected to: 0 1 2 3 4 6 7 8 9 10 16 20 21 22 24 26 28 38 41 47 51 56 60 63 85 86 92 93 101 131 133 148 180 196 600 653
Vertex 6 connected to: 0 1 2 3 4 5 6 6 7 8 9 10 16 20 21 22 24 38 41 47 49 51 56 60 63 85 86 92 93 101 131 133 148 173 180 196 631 653
Vertex 7 connected to: 0 1 2 3 5 6 9 16 20 21 22 24 26 28 38 47 49 51 56 60 63 85 92 93 101 131 133 148 173 180 653
Vertex 8 connected to: 0 1 2 3 4 5 6 8 8 9 10 21 22 24 26 47 56 60 63 92 93 101 131 148 173 180 631 653
Vertex 9 connected to: 0 1 2 3 5 6 7 8 9 9 10 16 20 21 22 24 26 28 38 41 47 49 51 56 60 63 85 86 92 101 131 133 148 180 196 631 653
Vertex 10 connected to: 0 1 4 5 6 8 9 10 10 22 24 26 41 101 131 133
Vertex 16 connected to: 0 1 2 3 4 5 6 7 9 16 16 26 28 41
Vertex 20 connected to: 0 1 2 3 4 5 6 7 9 41 47 49
Vertex 21 connected to: 0 1 2 3 5 6 7 8 9 24 173 600
Vertex 22 connected to: 0 1 2 3 4 5 6 7 8 9 10 28 47 600
Vertex 24 connected to: 0 2 3 4 5 6 7 8 9 10 21 41 47 196 600
Vertex 26 connected to: 0 1 2 3 4 5 7 8 9 10 16 28 38 41 47 173 196
Vertex 28 connected to: 0 1 2 3 5 7 9 16 22 26 41 51 56 60 86 180 631
Vertex 38 connected to: 0 1 2 3 4 5 6 7 9 26 38 38
Vertex 41 connected to: 0 1 2 3 4 5 6 9 10 16 20 24 26 28 47 92 101 653
Vertex 47 connected to: 0 1 2 3 4 5 6 7 8 9 20 22 24 26 41 47 47 49 92 173 653
Vertex 49 connected to: 0 1 2 3 6 7 9 20 47 196 600
Vertex 51 connected to: 0 1 2 3 4 5 6 7 9 28 51 51 600
Vertex 56 connected to: 0 1 2 3 4 5 6 7 8 9 28 600
Vertex 60 connected to: 0 1 2 3 4 5 6 7 8 9 28 60 60 101
Vertex 63 connected to: 0 1 2 3 5 6 7 8 9 600
Vertex 85 connected to: 0 1 2 3 4 5 6 7 9
Vertex 86 connected to: 0 1 2 3 4 5 6 9 28 93 600
Vertex 92 connected to: 0 1 2 3 4 5 6 7 8 9 41 47
Vertex 93 connected to: 0 1 2 5 6 7 8 86 93 93 600
Vertex 101 connected to: 0 1 2 3 4 5 6 7 8 9 10 41 60 133
Vertex 131 connected to: 0 1 2 3 4 5 6 7 8 9 10 600
Vertex 133 connected to: 0 2 4 5 6 7 9 10 101
Vertex 148 connected to: 0 1 2 3 4 5 6 7 8 9 148 148
Vertex 173 connected to: 0 1 2 3 6 7 8 21 26 47 173 173 600
Vertex 180 connected to: 0 1 2 3 4 5 6 7 8 9 28 600
Vertex 196 connected to: 0 1 2 3 4 5 6 9 24 26 49
Vertex 600 connected to: 1 3 5 21 22 24 49 51 56 63 86 93 131 173 180 600 600 631
Vertex 631 connected to: 1 2 3 4 6 8 9 28 600 631 631
Vertex 653 connected to: 2 3 4 5 6 7 8 9 41 47 653 653
Number of edges: 371
Density: 0.475641
Number of k-cliques in dense subgraph: 355
Average k-cliques per vertex: 8.875


```

For k=3

```bash

Current alpha: 1023.5
Max Flow: 19752
Min cut (nodes in S side): 0
Current gap: 1023.5

Current alpha: 511.75
Max Flow: 19752
Min cut (nodes in S side): 0
Current gap: 511.75

Current alpha: 255.875
Max Flow: 19752
Min cut (nodes in S side): 0
Current gap: 255.875

Current alpha: 127.938
Max Flow: 19752
Min cut (nodes in S side): 0
Current gap: 127.938

Current alpha: 63.9688
Max Flow: 19752
Min cut (nodes in S side): 0
Current gap: 63.9688

Current alpha: 31.9844
Max Flow: 19300
Min cut (nodes in S side): 0 1 2 3 4 5 6 7 8 9 10 11 17 21 22 23 25 27 29 39 42 48 50 52 57 61 64 86 87 93 102 132 134 149 174 181 197 601 654 6475 6476 6477 6478 6479 6480 6481 6482 6483 6484 6490 6494 6495 6496 6498 6500 6502 6512 6515 6521 6523 6525 6530 6534 6537 6559 6560 6566 6575 6605 6607 6622 6647 6654 6670 6853 6854 6855 6856 6857 6858 6859 6860 6861 6863 6866 6867 6868 6871 6872 6880 6882 6886 6888 6889 6891 6894 6897 6915 6916 6921 6924 6936 6945 6952 6955 6964 7086 8310 8311 8312 8313 8314 8315 8316 8318 8319 8320 8321 8322 8323 8324 8327 8328 8329 8331 8332 8333 8335 8337 8340 8341 8342 8345 8347 8348 8351 8353 8354 8356 8392 8594 8595 8596 8597 8598 8599 8601 8603 8604 8605 8606 8607 8609 8612 8613 8615 8617 8619 8622 8623 8624 8627 8628 8629 8631 8633 8636 8637 8639 8641 8657 8659 8719 8720 8721 8722 8724 8726 8727 8728 8729 8730 8731 8733 8735 8736 8737 8738 8739 8740 8741 8743 8744 8745 8746 8747 8755 8796 8797 8798 8799 8800 8802 8804 8805 8806 8807 8808 8810 8813 8814 8816 8817 8820 8821 8822 8824 8825 8826 8828 8831 8832 8834 8835 8837 8868 8870 8931 8932 8933 8934 8935 8938 8939 8940 8941 8945 8946 8948 8950 8951 8952 8955 8957 8964 8965 8966 8968 8974 8975 8977 8979 8981 8983 9029 9616 9620 9621 9622 9623 9625 9626 9627 9632 9634 9636 9637 9639 9641 9643 9649 9650 9652 9655 9656 9658 9661 9663 9696 10011 10012 10013 10014 10015 10016 10019 10020 10021 10022 10024 10026 10027 10028 10030 10031 10044 10113 10115 10118 10119 10120 10121 10123 10125 10130 10131 10133 10135 10136 10138 10141 10143 10150 10151 10155 10156 10166 10167 10170 10175 10178 10236 10855 10856 10857 10858 10859 10860 10861 10934 10935 10936 11010 11011 11012 11064 11067 11077 11167 11168 11193 11408 11409 11411 11414 11499 11500 11501 11502 11507 11510 11631 11632 11633 11635 11637 11640 11938 11939 11940 11990 12173 12174 12176 12189 12268 12270 12348 12456 12536 12727 12989 13042 13150 13228 13242
Current gap: 31.9844

Current alpha: 47.9766
Max Flow: 19752
Min cut (nodes in S side): 0
Current gap: 15.9922

Current alpha: 39.9805
Max Flow: 19752
Min cut (nodes in S side): 0
Current gap: 7.99609

Current alpha: 35.9824
Max Flow: 19728
Min cut (nodes in S side): 0 1 2 3 4 5 6 7 8 9 10 11 17 21 22 23 25 27 29 39 42 48 50 52 57 61 93 102 132 149 174 181 197 654 6475 6476 6477 6478 6479 6480 6481 6482 6483 6484 6490 6494 6495 6496 6498 6500 6502 6512 6515 6521 6523 6525 6530 6534 6566 6575 6605 6622 6647 6654 6670 6853 6854 6855 6856 6857 6858 6859 6860 6861 6863 6866 6867 6868 6871 6872 6880 6882 6886 6888 6889 6891 6894 6921 6924 6936 6945 6952 6955 6964 8310 8311 8312 8313 8314 8315 8316 8318 8319 8320 8321 8322 8323 8324 8327 8328 8329 8331 8332 8333 8335 8342 8345 8347 8351 8353 8354 8356 8392 8594 8595 8596 8597 8598 8599 8601 8603 8604 8605 8606 8607 8609 8612 8613 8615 8617 8619 8622 8623 8629 8631 8633 8636 8637 8639 8641 8659 8719 8720 8721 8722 8724 8726 8727 8728 8729 8730 8731 8733 8735 8736 8737 8740 8741 8743 8745 8746 8747 8755 8796 8797 8798 8799 8800 8802 8804 8805 8806 8807 8808 8810 8813 8814 8816 8817 8820 8821 8826 8828 8831 8834 8835 8837 8870 8931 8932 8933 8934 8935 8938 8939 8940 8941 8945 8946 8948 8950 8951 8952 8955 8966 8968 8974 8977 8979 8981 8983 9029 9616 9620 9621 9622 9623 9625 9626 9627 9632 9634 9636 9637 9639 9641 9650 9652 9655 9658 9661 9663 9696 10011 10012 10013 10014 10015 10016 10019 10020 10021 10024 10026 10027 10028 10030 10031 10044 10113 10115 10118 10119 10120 10121 10123 10125 10130 10131 10133 10135 10136 10138 10141 10155 10156 10166 10170 10175 10178 10236 10855 10856 10857 10858 10859 10860 10934 10935 10936 11010 11011 11012 11064 11067 11167 11168 11408 11409 11411 11499 11500 11501 11502 11507 11510 11631 11632 11633 11635 11640 11938 11939 11940 11990 12173 12174 12176 12189 12268 12536
Current gap: 3.99805

Current alpha: 37.9814
Max Flow: 19752
Min cut (nodes in S side): 0
Current gap: 1.99902

Current alpha: 36.9819
Max Flow: 19752
Min cut (nodes in S side): 0
Current gap: 0.999512

Current alpha: 36.4822
Max Flow: 19752
Min cut (nodes in S side): 0
Current gap: 0.499756

Current alpha: 36.2323
Max Flow: 19752
Min cut (nodes in S side): 0
Current gap: 0.249878

Current alpha: 36.1074
Max Flow: 19752
Min cut (nodes in S side): 0
Current gap: 0.124939

Current alpha: 36.0449
Max Flow: 19752
Min cut (nodes in S side): 0
Current gap: 0.0624695

Current alpha: 36.0137
Max Flow: 19752
Min cut (nodes in S side): 0
Current gap: 0.0312347

Current alpha: 35.998
Max Flow: 19728
Min cut (nodes in S side): 0 1 2 3 4 5 6 7 8 9 10 11 17 21 22 23 25 27 29 39 42 48 50 52 57 61 93 102 132 149 174 181 197 654 6475 6476 6477 6478 6479 6480 6481 6482 6483 6484 6490 6494 6495 6496 6498 6500 6502 6512 6515 6521 6523 6525 6530 6534 6566 6575 6605 6622 6647 6654 6670 6853 6854 6855 6856 6857 6858 6859 6860 6861 6863 6866 6867 6868 6871 6872 6880 6882 6886 6888 6889 6891 6894 6921 6924 6936 6945 6952 6955 6964 8310 8311 8312 8313 8314 8315 8316 8318 8319 8320 8321 8322 8323 8324 8327 8328 8329 8331 8332 8333 8335 8342 8345 8347 8351 8353 8354 8356 8392 8594 8595 8596 8597 8598 8599 8601 8603 8604 8605 8606 8607 8609 8612 8613 8615 8617 8619 8622 8623 8629 8631 8633 8636 8637 8639 8641 8659 8719 8720 8721 8722 8724 8726 8727 8728 8729 8730 8731 8733 8735 8736 8737 8740 8741 8743 8745 8746 8747 8755 8796 8797 8798 8799 8800 8802 8804 8805 8806 8807 8808 8810 8813 8814 8816 8817 8820 8821 8826 8828 8831 8834 8835 8837 8870 8931 8932 8933 8934 8935 8938 8939 8940 8941 8945 8946 8948 8950 8951 8952 8955 8966 8968 8974 8977 8979 8981 8983 9029 9616 9620 9621 9622 9623 9625 9626 9627 9632 9634 9636 9637 9639 9641 9650 9652 9655 9658 9661 9663 9696 10011 10012 10013 10014 10015 10016 10019 10020 10021 10024 10026 10027 10028 10030 10031 10044 10113 10115 10118 10119 10120 10121 10123 10125 10130 10131 10133 10135 10136 10138 10141 10155 10156 10166 10170 10175 10178 10236 10855 10856 10857 10858 10859 10860 10934 10935 10936 11010 11011 11012 11064 11067 11167 11168 11408 11409 11411 11499 11500 11501 11502 11507 11510 11631 11632 11633 11635 11640 11938 11939 11940 11990 12173 12174 12176 12189 12268 12536
Current gap: 0.0156174

Current alpha: 36.0058
Max Flow: 19752
Min cut (nodes in S side): 0
Current gap: 0.00780869

Current alpha: 36.0019
Max Flow: 19752
Min cut (nodes in S side): 0
Current gap: 0.00390434

Current alpha: 36
Max Flow: 19728
Min cut (nodes in S side): 0 1 2 3 4 5 6 7 8 9 10 11 17 21 22 23 25 27 29 39 42 48 50 52 57 61 93 102 132 149 174 181 197 654 6475 6476 6477 6478 6479 6480 6481 6482 6483 6484 6490 6494 6495 6496 6498 6500 6502 6512 6515 6521 6523 6525 6530 6534 6566 6575 6605 6622 6647 6654 6670 6853 6854 6855 6856 6857 6858 6859 6860 6861 6863 6866 6867 6868 6871 6872 6880 6882 6886 6888 6889 6891 6894 6921 6924 6936 6945 6952 6955 6964 8310 8311 8312 8313 8314 8315 8316 8318 8319 8320 8321 8322 8323 8324 8327 8328 8329 8331 8332 8333 8335 8342 8345 8347 8351 8353 8354 8356 8392 8594 8595 8596 8597 8598 8599 8601 8603 8604 8605 8606 8607 8609 8612 8613 8615 8617 8619 8622 8623 8629 8631 8633 8636 8637 8639 8641 8659 8719 8720 8721 8722 8724 8726 8727 8728 8729 8730 8731 8733 8735 8736 8737 8740 8741 8743 8745 8746 8747 8755 8796 8797 8798 8799 8800 8802 8804 8805 8806 8807 8808 8810 8813 8814 8816 8817 8820 8821 8826 8828 8831 8834 8835 8837 8870 8931 8932 8933 8934 8935 8938 8939 8940 8941 8945 8946 8948 8950 8951 8952 8955 8966 8968 8974 8977 8979 8981 8983 9029 9616 9620 9621 9622 9623 9625 9626 9627 9632 9634 9636 9637 9639 9641 9650 9652 9655 9658 9661 9663 9696 10011 10012 10013 10014 10015 10016 10019 10020 10021 10024 10026 10027 10028 10030 10031 10044 10113 10115 10118 10119 10120 10121 10123 10125 10130 10131 10133 10135 10136 10138 10141 10155 10156 10166 10170 10175 10178 10236 10855 10856 10857 10858 10859 10860 10934 10935 10936 11010 11011 11012 11064 11067 11167 11168 11408 11409 11411 11499 11500 11501 11502 11507 11510 11631 11632 11633 11635 11640 11938 11939 11940 11990 12173 12174 12176 12189 12268 12536
Current gap: 0.00195217

Current alpha: 36.001
Max Flow: 19752
Min cut (nodes in S side): 0
Current gap: 0.000976086

Current alpha: 36.0005
Max Flow: 19752
Min cut (nodes in S side): 0
Current gap: 0.000488043

Current alpha: 36.0002
Max Flow: 19752
Min cut (nodes in S side): 0
Current gap: 0.000244021

Current alpha: 36.0001
Max Flow: 19752
Min cut (nodes in S side): 0
Current gap: 0.000122011

Current alpha: 36.0001
Max Flow: 19752
Min cut (nodes in S side): 0
Current gap: 6.10054e-005

Current alpha: 36
Max Flow: 19752
Min cut (nodes in S side): 0
Current gap: 3.05027e-005

Current alpha: 36
Max Flow: 19752
Min cut (nodes in S side): 0
Current gap: 1.52513e-005

Current alpha: 36
Max Flow: 19728
Min cut (nodes in S side): 0 1 2 3 4 5 6 7 8 9 10 11 17 21 22 23 25 27 29 39 42 48 50 52 57 61 93 102 132 149 174 181 197 654 6475 6476 6477 6478 6479 6480 6481 6482 6483 6484 6490 6494 6495 6496 6498 6500 6502 6512 6515 6521 6523 6525 6530 6534 6566 6575 6605 6622 6647 6654 6670 6853 6854 6855 6856 6857 6858 6859 6860 6861 6863 6866 6867 6868 6871 6872 6880 6882 6886 6888 6889 6891 6894 6921 6924 6936 6945 6952 6955 6964 8310 8311 8312 8313 8314 8315 8316 8318 8319 8320 8321 8322 8323 8324 8327 8328 8329 8331 8332 8333 8335 8342 8345 8347 8351 8353 8354 8356 8392 8594 8595 8596 8597 8598 8599 8601 8603 8604 8605 8606 8607 8609 8612 8613 8615 8617 8619 8622 8623 8629 8631 8633 8636 8637 8639 8641 8659 8719 8720 8721 8722 8724 8726 8727 8728 8729 8730 8731 8733 8735 8736 8737 8740 8741 8743 8745 8746 8747 8755 8796 8797 8798 8799 8800 8802 8804 8805 8806 8807 8808 8810 8813 8814 8816 8817 8820 8821 8826 8828 8831 8834 8835 8837 8870 8931 8932 8933 8934 8935 8938 8939 8940 8941 8945 8946 8948 8950 8951 8952 8955 8966 8968 8974 8977 8979 8981 8983 9029 9616 9620 9621 9622 9623 9625 9626 9627 9632 9634 9636 9637 9639 9641 9650 9652 9655 9658 9661 9663 9696 10011 10012 10013 10014 10015 10016 10019 10020 10021 10024 10026 10027 10028 10030 10031 10044 10113 10115 10118 10119 10120 10121 10123 10125 10130 10131 10133 10135 10136 10138 10141 10155 10156 10166 10170 10175 10178 10236 10855 10856 10857 10858 10859 10860 10934 10935 10936 11010 11011 11012 11064 11067 11167 11168 11408 11409 11411 11499 11500 11501 11502 11507 11510 11631 11632 11633 11635 11640 11938 11939 11940 11990 12173 12174 12176 12189 12268 12536
Current gap: 7.62567e-006

Current alpha: 36
Max Flow: 19752
Min cut (nodes in S side): 0
Current gap: 3.81283e-006

Current alpha: 36
Max Flow: 19752
Min cut (nodes in S side): 0
Current gap: 1.90642e-006

Current alpha: 36
Max Flow: 19728
Min cut (nodes in S side): 0 1 2 3 4 5 6 7 8 9 10 11 17 21 22 23 25 27 29 39 42 48 50 52 57 61 93 102 132 149 174 181 197 654 6475 6476 6477 6478 6479 6480 6481 6482 6483 6484 6490 6494 6495 6496 6498 6500 6502 6512 6515 6521 6523 6525 6530 6534 6566 6575 6605 6622 6647 6654 6670 6853 6854 6855 6856 6857 6858 6859 6860 6861 6863 6866 6867 6868 6871 6872 6880 6882 6886 6888 6889 6891 6894 6921 6924 6936 6945 6952 6955 6964 8310 8311 8312 8313 8314 8315 8316 8318 8319 8320 8321 8322 8323 8324 8327 8328 8329 8331 8332 8333 8335 8342 8345 8347 8351 8353 8354 8356 8392 8594 8595 8596 8597 8598 8599 8601 8603 8604 8605 8606 8607 8609 8612 8613 8615 8617 8619 8622 8623 8629 8631 8633 8636 8637 8639 8641 8659 8719 8720 8721 8722 8724 8726 8727 8728 8729 8730 8731 8733 8735 8736 8737 8740 8741 8743 8745 8746 8747 8755 8796 8797 8798 8799 8800 8802 8804 8805 8806 8807 8808 8810 8813 8814 8816 8817 8820 8821 8826 8828 8831 8834 8835 8837 8870 8931 8932 8933 8934 8935 8938 8939 8940 8941 8945 8946 8948 8950 8951 8952 8955 8966 8968 8974 8977 8979 8981 8983 9029 9616 9620 9621 9622 9623 9625 9626 9627 9632 9634 9636 9637 9639 9641 9650 9652 9655 9658 9661 9663 9696 10011 10012 10013 10014 10015 10016 10019 10020 10021 10024 10026 10027 10028 10030 10031 10044 10113 10115 10118 10119 10120 10121 10123 10125 10130 10131 10133 10135 10136 10138 10141 10155 10156 10166 10170 10175 10178 10236 10855 10856 10857 10858 10859 10860 10934 10935 10936 11010 11011 11012 11064 11067 11167 11168 11408 11409 11411 11499 11500 11501 11502 11507 11510 11631 11632 11633 11635 11640 11938 11939 11940 11990 12173 12174 12176 12189 12268 12536
Current gap: 9.53209e-007

Current alpha: 36
Max Flow: 19752
Min cut (nodes in S side): 0
Current gap: 4.76604e-007

Current alpha: 36
Max Flow: 19752
Min cut (nodes in S side): 0
Current gap: 2.38302e-007

Current alpha: 36
Max Flow: 19752
Min cut (nodes in S side): 0
Current gap: 1.19151e-007

Current alpha: 36
Max Flow: 19752
Min cut (nodes in S side): 0
Current gap: 5.95755e-008

Current alpha: 36
Max Flow: 19752
Min cut (nodes in S side): 0
Current gap: 2.97878e-008

Current alpha: 36
Max Flow: 19752
Min cut (nodes in S side): 0
Current gap: 1.48939e-008


Dense subgraph found with 33 vertices:
Vertices in dense subgraph: 0 1 2 3 4 5 6 7 8 9 10 16 20 21 22 24 26 28 38 41 47 49 51 56 60 92 101 131 148 173 180 196 653
Edges in dense subgraph:
Vertex 0 connected to: 1 2 3 4 5 6 7 8 9 10 16 20 21 22 24 26 28 38 41 47 49 51 56 60 92 101 131 148 173 180 196
Vertex 1 connected to: 0 1 1 2 3 4 5 6 7 8 9 10 16 20 21 22 26 28 38 41 47 49 51 56 60 92 101 131 148 173 180 196
Vertex 2 connected to: 0 1 3 4 5 6 7 8 9 16 20 21 22 24 26 28 38 41 47 49 51 56 60 92 101 131 148 173 180 196 653
Vertex 3 connected to: 0 1 2 4 5 6 7 8 9 16 20 21 22 24 26 28 38 41 47 49 51 56 60 92 101 131 148 173 180 196 653
Vertex 4 connected to: 0 1 2 3 5 6 8 10 16 20 22 24 26 38 41 47 51 56 60 92 101 131 148 180 196 653
Vertex 5 connected to: 0 1 2 3 4 6 7 8 9 10 16 20 21 22 24 26 28 38 41 47 51 56 60 92 101 131 148 180 196 653
Vertex 6 connected to: 0 1 2 3 4 5 6 6 7 8 9 10 16 20 21 22 24 38 41 47 49 51 56 60 92 101 131 148 173 180 196 653
Vertex 7 connected to: 0 1 2 3 5 6 9 16 20 21 22 24 26 28 38 47 49 51 56 60 92 101 131 148 173 180 653
Vertex 8 connected to: 0 1 2 3 4 5 6 8 8 9 10 21 22 24 26 47 56 60 92 101 131 148 173 180 653
Vertex 9 connected to: 0 1 2 3 5 6 7 8 9 9 10 16 20 21 22 24 26 28 38 41 47 49 51 56 60 92 101 131 148 180 196 653
Vertex 10 connected to: 0 1 4 5 6 8 9 10 10 22 24 26 41 101 131
Vertex 16 connected to: 0 1 2 3 4 5 6 7 9 16 16 26 28 41
Vertex 20 connected to: 0 1 2 3 4 5 6 7 9 41 47 49
Vertex 21 connected to: 0 1 2 3 5 6 7 8 9 24 173
Vertex 22 connected to: 0 1 2 3 4 5 6 7 8 9 10 28 47
Vertex 24 connected to: 0 2 3 4 5 6 7 8 9 10 21 41 47 196
Vertex 26 connected to: 0 1 2 3 4 5 7 8 9 10 16 28 38 41 47 173 196
Vertex 28 connected to: 0 1 2 3 5 7 9 16 22 26 41 51 56 60 180
Vertex 38 connected to: 0 1 2 3 4 5 6 7 9 26 38 38
Vertex 41 connected to: 0 1 2 3 4 5 6 9 10 16 20 24 26 28 47 92 101 653
Vertex 47 connected to: 0 1 2 3 4 5 6 7 8 9 20 22 24 26 41 47 47 49 92 173 653
Vertex 49 connected to: 0 1 2 3 6 7 9 20 47 196
Vertex 51 connected to: 0 1 2 3 4 5 6 7 9 28 51 51
Vertex 56 connected to: 0 1 2 3 4 5 6 7 8 9 28
Vertex 60 connected to: 0 1 2 3 4 5 6 7 8 9 28 60 60 101
Vertex 92 connected to: 0 1 2 3 4 5 6 7 8 9 41 47
Vertex 101 connected to: 0 1 2 3 4 5 6 7 8 9 10 41 60
Vertex 131 connected to: 0 1 2 3 4 5 6 7 8 9 10
Vertex 148 connected to: 0 1 2 3 4 5 6 7 8 9 148 148
Vertex 173 connected to: 0 1 2 3 6 7 8 21 26 47 173 173
Vertex 180 connected to: 0 1 2 3 4 5 6 7 8 9 28
Vertex 196 connected to: 0 1 2 3 4 5 6 9 24 26 49
Vertex 653 connected to: 2 3 4 5 6 7 8 9 41 47 653 653
Number of edges: 300
Density: 0.568182
Number of k-cliques in dense subgraph: 1185
Average k-cliques per vertex: 35.9091


```

For k=4

```bash

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
Current gap: 6.1363e-005

Current alpha: 85.25
Max Flow: 22544
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 85.25
Current gap: 3.06815e-005

Current alpha: 85.25
Max Flow: 22544
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 85.25
Current gap: 1.53407e-005

Current alpha: 85.25
Max Flow: 22544
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 85.25
Current gap: 7.67037e-006

Current alpha: 85.25
Max Flow: 22528
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 85.25
Current gap: 3.83519e-006

Current alpha: 85.25
Max Flow: 22544
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 85.25
Current gap: 1.91759e-006

Current alpha: 85.25
Max Flow: 22528
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 85.25
Current gap: 9.58797e-007

Current alpha: 85.25
Max Flow: 22544
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 85.25
Current gap: 4.79398e-007

Current alpha: 85.25
Max Flow: 22544
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 85.25
Current gap: 2.39699e-007

Current alpha: 85.25
Max Flow: 22544
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 85.25
Current gap: 1.1985e-007

Current alpha: 85.25
Max Flow: 22544
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 85.25
Current gap: 5.99248e-008

Current alpha: 85.25
Max Flow: 22528
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 85.25
Current gap: 2.99624e-008

Current alpha: 85.25
Max Flow: 22544
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 85.25
Current gap: 1.49812e-008


Dense subgraph found with 32 vertices:
Vertices in dense subgraph: 0 1 2 3 4 5 6 7 8 9 10 16 20 21 22 24 26 28 38 41 47 49 51 56 60 92 101 131 148 173 180 196
Edges in dense subgraph:
Vertex 0 connected to: 1 2 3 4 5 6 7 8 9 10 16 20 21 22 24 26 28 38 41 47 49 51 56 60 92 101 131 148 173 180 196
Vertex 1 connected to: 0 1 1 2 3 4 5 6 7 8 9 10 16 20 21 22 26 28 38 41 47 49 51 56 60 92 101 131 148 173 180 196
Vertex 2 connected to: 0 1 3 4 5 6 7 8 9 16 20 21 22 24 26 28 38 41 47 49 51 56 60 92 101 131 148 173 180 196
Vertex 3 connected to: 0 1 2 4 5 6 7 8 9 16 20 21 22 24 26 28 38 41 47 49 51 56 60 92 101 131 148 173 180 196
Vertex 4 connected to: 0 1 2 3 5 6 8 10 16 20 22 24 26 38 41 47 51 56 60 92 101 131 148 180 196
Vertex 5 connected to: 0 1 2 3 4 6 7 8 9 10 16 20 21 22 24 26 28 38 41 47 51 56 60 92 101 131 148 180 196
Vertex 6 connected to: 0 1 2 3 4 5 6 6 7 8 9 10 16 20 21 22 24 38 41 47 49 51 56 60 92 101 131 148 173 180 196
Vertex 7 connected to: 0 1 2 3 5 6 9 16 20 21 22 24 26 28 38 47 49 51 56 60 92 101 131 148 173 180
Vertex 8 connected to: 0 1 2 3 4 5 6 8 8 9 10 21 22 24 26 47 56 60 92 101 131 148 173 180
Vertex 9 connected to: 0 1 2 3 5 6 7 8 9 9 10 16 20 21 22 24 26 28 38 41 47 49 51 56 60 92 101 131 148 180 196
Vertex 10 connected to: 0 1 4 5 6 8 9 10 10 22 24 26 41 101 131
Vertex 16 connected to: 0 1 2 3 4 5 6 7 9 16 16 26 28 41
Vertex 20 connected to: 0 1 2 3 4 5 6 7 9 41 47 49
Vertex 21 connected to: 0 1 2 3 5 6 7 8 9 24 173
Vertex 22 connected to: 0 1 2 3 4 5 6 7 8 9 10 28 47
Vertex 24 connected to: 0 2 3 4 5 6 7 8 9 10 21 41 47 196
Vertex 26 connected to: 0 1 2 3 4 5 7 8 9 10 16 28 38 41 47 173 196
Vertex 28 connected to: 0 1 2 3 5 7 9 16 22 26 41 51 56 60 180
Vertex 38 connected to: 0 1 2 3 4 5 6 7 9 26 38 38
Vertex 41 connected to: 0 1 2 3 4 5 6 9 10 16 20 24 26 28 47 92 101
Vertex 47 connected to: 0 1 2 3 4 5 6 7 8 9 20 22 24 26 41 47 47 49 92 173
Vertex 49 connected to: 0 1 2 3 6 7 9 20 47 196
Vertex 51 connected to: 0 1 2 3 4 5 6 7 9 28 51 51
Vertex 56 connected to: 0 1 2 3 4 5 6 7 8 9 28
Vertex 60 connected to: 0 1 2 3 4 5 6 7 8 9 28 60 60 101
Vertex 92 connected to: 0 1 2 3 4 5 6 7 8 9 41 47
Vertex 101 connected to: 0 1 2 3 4 5 6 7 8 9 10 41 60
Vertex 131 connected to: 0 1 2 3 4 5 6 7 8 9 10
Vertex 148 connected to: 0 1 2 3 4 5 6 7 8 9 148 148
Vertex 173 connected to: 0 1 2 3 6 7 8 21 26 47 173 173
Vertex 180 connected to: 0 1 2 3 4 5 6 7 8 9 28
Vertex 196 connected to: 0 1 2 3 4 5 6 9 24 26 49
Number of edges: 289
Density: 0.582661
Number of k-cliques in dense subgraph: 2724
Average k-cliques per vertex: 85.125


```

#### As733 (Closest to specified)

For k=2

```bash

Current alpha: 207.5
Max Flow: 6264
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 207.5
Current gap: 207.5

Current alpha: 103.75
Max Flow: 6264
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 103.75
Current gap: 103.75

Current alpha: 51.875
Max Flow: 6264
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 51.875
Current gap: 51.875

Current alpha: 25.9375
Max Flow: 6264
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 25.9375
Current gap: 25.9375

Current alpha: 12.9688
Max Flow: 6264
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 12.9688
Current gap: 12.9688

Current alpha: 6.48438
Max Flow: 6084
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 6.48438
Current gap: 6.48438

Current alpha: 9.72656
Max Flow: 6264
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 9.72656
Current gap: 3.24219

Current alpha: 8.10547
Max Flow: 6256
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 8.10547
Current gap: 1.62109

Current alpha: 8.91602
Max Flow: 6264
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 8.91602
Current gap: 0.810547

Current alpha: 8.51074
Max Flow: 6264
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 8.51074
Current gap: 0.405273

Current alpha: 8.30811
Max Flow: 6256
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 8.30811
Current gap: 0.202637

Current alpha: 8.40942
Max Flow: 6256
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 8.40942
Current gap: 0.101318

Current alpha: 8.46008
Max Flow: 6256
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 8.46008
Current gap: 0.0506592

Current alpha: 8.48541
Max Flow: 6256
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 8.48541
Current gap: 0.0253296

Current alpha: 8.49808
Max Flow: 6256
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 8.49808
Current gap: 0.0126648

Current alpha: 8.50441
Max Flow: 6264
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 8.50441
Current gap: 0.0063324

Current alpha: 8.50124
Max Flow: 6264
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 8.50124
Current gap: 0.0031662

Current alpha: 8.49966
Max Flow: 6256
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 8.49966
Current gap: 0.0015831

Current alpha: 8.50045
Max Flow: 6264
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 8.50045
Current gap: 0.00079155

Current alpha: 8.50006
Max Flow: 6264
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 8.50006
Current gap: 0.000395775

Current alpha: 8.49986
Max Flow: 6256
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 8.49986
Current gap: 0.000197887

Current alpha: 8.49996
Max Flow: 6256
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 8.49996
Current gap: 9.89437e-005

Current alpha: 8.50001
Max Flow: 6264
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 8.50001
Current gap: 4.94719e-005

Current alpha: 8.49998
Max Flow: 6256
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 8.49998
Current gap: 2.47359e-005

Current alpha: 8.49999
Max Flow: 6256
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 8.49999
Current gap: 1.2368e-005

Current alpha: 8.5
Max Flow: 6264
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 8.5
Current gap: 6.18398e-006

Current alpha: 8.5
Max Flow: 6256
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 8.5
Current gap: 3.09199e-006

Current alpha: 8.5
Max Flow: 6256
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 8.5
Current gap: 1.546e-006

Current alpha: 8.5
Max Flow: 6256
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 8.5
Current gap: 7.72998e-007

Current alpha: 8.5
Max Flow: 6264
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 8.5
Current gap: 3.86499e-007


Dense subgraph found with 30 vertices:
Vertices in dense subgraph: 486 499 605 660 691 726 735 760 780 786 788 801 833 837 846 852 854 865 870 879 888 900 924 940 951 960 973 1034 1273 1361
Edges in dense subgraph:
Vertex 486 connected to: 486 486 660 660 691 691 735 735 786 786 837 837 870 870 888 888 924 924 1273 1273
Vertex 499 connected to: 605 605 691 691 735 735 786 786 837 837 865 865 870 870 879 879 888 888 1273 1273
Vertex 605 connected to: 499 499 660 660 691 691 726 726 735 735 788 788 801 801 846 846 852 852 865 865 870 870 879 879 888 888 900 900 924 924 951 951 960 960 973 973 1034 1034 1273 1273 1361 1361
Vertex 660 connected to: 486 486 605 605 660 660 726 726 735 735 780 780 786 786 833 833 852 852 854 854 865 865 870 870 879 879 888 888 900 900 940 940 951 951 973 973 1034 1034 1273 1273 1361 1361
Vertex 691 connected to: 486 486 499 499 605 605 780 780 786 786 788 788 801 801 833 833 846 846 852 852 854 854 865 865 870 870 879 879 888 888 900 900 951 951 960 960 1273 1273 1361 1361
Vertex 726 connected to: 605 605 660 660 735 735 760 760 786 786 837 837 865 865 870 870 879 879 888 888 1273 1273
Vertex 735 connected to: 486 486 499 499 605 605 660 660 726 726 760 760 780 780 786 786 801 801 833 833 837 837 846 846 852 852 865 865 870 870 879 879 888 888 900 900 924 924 940 940 951 951 960 960 973 973 1034 1034 1273 1273 1361 1361
Vertex 760 connected to: 726 726 735 735 786 786 837 837 846 846 865 865 879 879 888 888 900 900 924 924 1273 1273
Vertex 780 connected to: 660 660 691 691 735 735 786 786 837 837 865 865 870 870 879 879 888 888 1034 1034 1273 1273
Vertex 786 connected to: 486 486 499 499 660 660 691 691 726 726 735 735 760 760 780 780 786 786 788 788 801 801 833 833 837 837 846 846 852 852 854 854 865 865 870 870 879 879 888 888 900 900 924 924 940 940 951 951 960 960 973 973 1034 1034 1273 1273 1361 1361
Vertex 788 connected to: 605 605 691 691 786 786 837 837 852 852 854 854 870 870 879 879 1273 1273
Vertex 801 connected to: 605 605 691 691 735 735 786 786 837 837 865 865 870 870 879 879 888 888 1273 1273
Vertex 833 connected to: 660 660 691 691 735 735 786 786 837 837 865 865 870 870 879 879 888 888 1273 1273
Vertex 837 connected to: 486 486 499 499 726 726 735 735 760 760 780 780 786 786 788 788 801 801 833 833 846 846 852 852 865 865 870 870 879 879 888 888 900 900 940 940 951 951 960 960 973 973 1273 1273 1361 1361
Vertex 846 connected to: 605 605 691 691 735 735 760 760 786 786 837 837 865 865 870 870 879 879 888 888 924 924 1273 1273 1361 1361
Vertex 852 connected to: 605 605 660 660 691 691 735 735 786 786 788 788 837 837 854 854 865 865 870 870 879 879 888 888 924 924 1273 1273
Vertex 854 connected to: 660 660 691 691 786 786 788 788 852 852 865 865 870 870 879 879 888 888 900 900 973 973 1034 1034 1361 1361
Vertex 865 connected to: 499 499 605 605 660 660 691 691 726 726 735 735 760 760 780 780 786 786 801 801 833 833 837 837 846 846 852 852 854 854 865 865 870 870 879 879 888 888 900 900 924 924 940 940 951 951 960 960 973 973 1273 1273 1361 1361
Vertex 870 connected to: 486 486 499 499 605 605 660 660 691 691 726 726 735 735 780 780 786 786 788 788 801 801 833 833 837 837 846 846 852 852 854 854 865 865 870 870 879 879 888 888 900 900 924 924 940 940 951 951 960 960 973 973 1034 1034 1273 1273
Vertex 879 connected to: 499 499 605 605 660 660 691 691 726 726 735 735 760 760 780 780 786 786 788 788 801 801 833 833 837 837 846 846 852 852 854 854 865 865 870 870 888 888 900 900 924 924 940 940 951 951 960 960 973 973 1034 1034 1273 1273 1361 1361
Vertex 888 connected to: 486 486 499 499 605 605 660 660 691 691 726 726 735 735 760 760 780 780 786 786 801 801 833 833 837 837 846 846 852 852 854 854 865 865 870 870 879 879 900 900 924 924 940 940 951 951 960 960 973 973 1034 1034 1273 1273 1361 1361
Vertex 900 connected to: 605 605 660 660 691 691 735 735 760 760 786 786 837 837 854 854 865 865 870 870 879 879 888 888 1273 1273
Vertex 924 connected to: 486 486 605 605 735 735 760 760 786 786 846 846 852 852 865 865 870 870 879 879 888 888 1273 1273 1361 1361
Vertex 940 connected to: 660 660 735 735 786 786 837 837 865 865 870 870 879 879 888 888 1273 1273
Vertex 951 connected to: 605 605 660 660 691 691 735 735 786 786 837 837 865 865 870 870 879 879 888 888 951 951 1273 1273
Vertex 960 connected to: 605 605 691 691 735 735 786 786 837 837 865 865 870 870 879 879 888 888 960 960 1273 1273 1361 1361
Vertex 973 connected to: 605 605 660 660 735 735 786 786 837 837 854 854 865 865 870 870 879 879 888 888 1273 1273
Vertex 1034 connected to: 605 605 660 660 735 735 780 780 786 786 854 854 870 870 879 879 888 888 1273 1273
Vertex 1273 connected to: 486 486 499 499 605 605 660 660 691 691 726 726 735 735 760 760 780 780 786 786 788 788 801 801 833 833 837 837 846 846 852 852 865 865 870 870 879 879 888 888 900 900 924 924 940 940 951 951 960 960 973 973 1034 1034 1361 1361
Vertex 1361 connected to: 605 605 660 660 691 691 735 735 786 786 837 837 846 846 854 854 865 865 879 879 888 888 924 924 960 960 1273 1273
Number of edges: 495
Density: 1.13793
Number of k-cliques in dense subgraph: 244
Average k-cliques per vertex: 8.13333
Execution time: 58623 milliseconds


```

For k=3

```bash

Current alpha: 320
Max Flow: 7509
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 320
Current gap: 320

Current alpha: 160
Max Flow: 7509
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 160
Current gap: 160

Current alpha: 80
Max Flow: 7509
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 80
Current gap: 80

Current alpha: 40
Max Flow: 7509
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 40
Current gap: 40

Current alpha: 20
Max Flow: 6417
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 20
Current gap: 20

Current alpha: 30
Max Flow: 7413
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 30
Current gap: 10

Current alpha: 35
Max Flow: 7509
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 35
Current gap: 5

Current alpha: 32.5
Max Flow: 7509
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 32.5
Current gap: 2.5

Current alpha: 31.25
Max Flow: 7497
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 31.25
Current gap: 1.25

Current alpha: 31.875
Max Flow: 7509
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 31.875
Current gap: 0.625

Current alpha: 31.5625
Max Flow: 7509
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 31.5625
Current gap: 0.3125

Current alpha: 31.4062
Max Flow: 7509
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 31.4062
Current gap: 0.15625

Current alpha: 31.3281
Max Flow: 7497
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 31.3281
Current gap: 0.078125

Current alpha: 31.3672
Max Flow: 7509
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 31.3672
Current gap: 0.0390625

Current alpha: 31.3477
Max Flow: 7509
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 31.3477
Current gap: 0.0195312

Current alpha: 31.3379
Max Flow: 7509
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 31.3379
Current gap: 0.00976562

Current alpha: 31.333
Max Flow: 7497
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 31.333
Current gap: 0.00488281

Current alpha: 31.3354
Max Flow: 7509
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 31.3354
Current gap: 0.00244141

Current alpha: 31.3342
Max Flow: 7509
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 31.3342
Current gap: 0.0012207

Current alpha: 31.3336
Max Flow: 7509
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 31.3336
Current gap: 0.000610352

Current alpha: 31.3333
Max Flow: 7497
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 31.3333
Current gap: 0.000305176

Current alpha: 31.3335
Max Flow: 7509
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 31.3335
Current gap: 0.000152588

Current alpha: 31.3334
Max Flow: 7509
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 31.3334
Current gap: 7.62939e-005

Current alpha: 31.3334
Max Flow: 7509
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 31.3334
Current gap: 3.8147e-005

Current alpha: 31.3333
Max Flow: 7497
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 31.3333
Current gap: 1.90735e-005

Current alpha: 31.3333
Max Flow: 7509
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 31.3333
Current gap: 9.53674e-006

Current alpha: 31.3333
Max Flow: 7509
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 31.3333
Current gap: 4.76837e-006

Current alpha: 31.3333
Max Flow: 7509
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 31.3333
Current gap: 2.38419e-006

Current alpha: 31.3333
Max Flow: 7497
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 31.3333
Current gap: 1.19209e-006

Current alpha: 31.3333
Max Flow: 7509
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 31.3333
Current gap: 5.96046e-007

Current alpha: 31.3333
Max Flow: 7509
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 31.3333
Current gap: 2.98023e-007


Dense subgraph found with 28 vertices:
Vertices in dense subgraph: 499 605 660 691 726 735 760 780 786 801 833 837 846 852 854 865 870 879 888 900 924 940 951 960 973 1034 1273 1361
Edges in dense subgraph:
Vertex 499 connected to: 605 605 691 691 735 735 786 786 837 837 865 865 870 870 879 879 888 888 1273 1273
Vertex 605 connected to: 499 499 660 660 691 691 726 726 735 735 801 801 846 846 852 852 865 865 870 870 879 879 888 888 900 900 924 924 951 951 960 960 973 973 1034 1034 1273 1273 1361 1361
Vertex 660 connected to: 605 605 660 660 726 726 735 735 780 780 786 786 833 833 852 852 854 854 865 865 870 870 879 879 888 888 900 900 940 940 951 951 973 973 1034 1034 1273 1273 1361 1361
Vertex 691 connected to: 499 499 605 605 780 780 786 786 801 801 833 833 846 846 852 852 854 854 865 865 870 870 879 879 888 888 900 900 951 951 960 960 1273 1273 1361 1361
Vertex 726 connected to: 605 605 660 660 735 735 760 760 786 786 837 837 865 865 870 870 879 879 888 888 1273 1273
Vertex 735 connected to: 499 499 605 605 660 660 726 726 760 760 780 780 786 786 801 801 833 833 837 837 846 846 852 852 865 865 870 870 879 879 888 888 900 900 924 924 940 940 951 951 960 960 973 973 1034 1034 1273 1273 1361 1361
Vertex 760 connected to: 726 726 735 735 786 786 837 837 846 846 865 865 879 879 888 888 900 900 924 924 1273 1273
Vertex 780 connected to: 660 660 691 691 735 735 786 786 837 837 865 865 870 870 879 879 888 888 1034 1034 1273 1273
Vertex 786 connected to: 499 499 660 660 691 691 726 726 735 735 760 760 780 780 786 786 801 801 833 833 837 837 846 846 852 852 854 854 865 865 870 870 879 879 888 888 900 900 924 924 940 940 951 951 960 960 973 973 1034 1034 1273 1273 1361 1361
Vertex 801 connected to: 605 605 691 691 735 735 786 786 837 837 865 865 870 870 879 879 888 888 1273 1273
Vertex 833 connected to: 660 660 691 691 735 735 786 786 837 837 865 865 870 870 879 879 888 888 1273 1273
Vertex 837 connected to: 499 499 726 726 735 735 760 760 780 780 786 786 801 801 833 833 846 846 852 852 865 865 870 870 879 879 888 888 900 900 940 940 951 951 960 960 973 973 1273 1273 1361 1361
Vertex 846 connected to: 605 605 691 691 735 735 760 760 786 786 837 837 865 865 870 870 879 879 888 888 924 924 1273 1273 1361 1361
Vertex 852 connected to: 605 605 660 660 691 691 735 735 786 786 837 837 854 854 865 865 870 870 879 879 888 888 924 924 1273 1273
Vertex 854 connected to: 660 660 691 691 786 786 852 852 865 865 870 870 879 879 888 888 900 900 973 973 1034 1034 1361 1361
Vertex 865 connected to: 499 499 605 605 660 660 691 691 726 726 735 735 760 760 780 780 786 786 801 801 833 833 837 837 846 846 852 852 854 854 865 865 870 870 879 879 888 888 900 900 924 924 940 940 951 951 960 960 973 973 1273 1273 1361 1361
Vertex 870 connected to: 499 499 605 605 660 660 691 691 726 726 735 735 780 780 786 786 801 801 833 833 837 837 846 846 852 852 854 854 865 865 870 870 879 879 888 888 900 900 924 924 940 940 951 951 960 960 973 973 1034 1034 1273 1273
Vertex 879 connected to: 499 499 605 605 660 660 691 691 726 726 735 735 760 760 780 780 786 786 801 801 833 833 837 837 846 846 852 852 854 854 865 865 870 870 888 888 900 900 924 924 940 940 951 951 960 960 973 973 1034 1034 1273 1273 1361 1361
Vertex 888 connected to: 499 499 605 605 660 660 691 691 726 726 735 735 760 760 780 780 786 786 801 801 833 833 837 837 846 846 852 852 854 854 865 865 870 870 879 879 900 900 924 924 940 940 951 951 960 960 973 973 1034 1034 1273 1273 1361 1361
Vertex 900 connected to: 605 605 660 660 691 691 735 735 760 760 786 786 837 837 854 854 865 865 870 870 879 879 888 888 1273 1273
Vertex 924 connected to: 605 605 735 735 760 760 786 786 846 846 852 852 865 865 870 870 879 879 888 888 1273 1273 1361 1361
Vertex 940 connected to: 660 660 735 735 786 786 837 837 865 865 870 870 879 879 888 888 1273 1273
Vertex 951 connected to: 605 605 660 660 691 691 735 735 786 786 837 837 865 865 870 870 879 879 888 888 951 951 1273 1273
Vertex 960 connected to: 605 605 691 691 735 735 786 786 837 837 865 865 870 870 879 879 888 888 960 960 1273 1273 1361 1361
Vertex 973 connected to: 605 605 660 660 735 735 786 786 837 837 854 854 865 865 870 870 879 879 888 888 1273 1273
Vertex 1034 connected to: 605 605 660 660 735 735 780 780 786 786 854 854 870 870 879 879 888 888 1273 1273
Vertex 1273 connected to: 499 499 605 605 660 660 691 691 726 726 735 735 760 760 780 780 786 786 801 801 833 833 837 837 846 846 852 852 865 865 870 870 879 879 888 888 900 900 924 924 940 940 951 951 960 960 973 973 1034 1034 1361 1361
Vertex 1361 connected to: 605 605 660 660 691 691 735 735 786 786 837 837 846 846 854 854 865 865 879 879 888 888 924 924 960 960 1273 1273
Number of edges: 458
Density: 1.21164
Number of k-cliques in dense subgraph: 872
Average k-cliques per vertex: 31.1429
Execution time: 192001 milliseconds

```

#### CaHepth Dataset

For k=2

```bash
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
Current gap: 6.19888e-005

Current alpha: 15.5
Max Flow: 51946
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 15.5
Current gap: 3.09944e-005

Current alpha: 15.5
Max Flow: 51946
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 15.5
Current gap: 1.54972e-005

Current alpha: 15.5
Max Flow: 51946
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 15.5
Current gap: 7.7486e-006

Current alpha: 15.5
Max Flow: 51946
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 15.5
Current gap: 3.8743e-006

Current alpha: 15.5
Max Flow: 51914
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 15.5
Current gap: 1.93715e-006

Current alpha: 15.5
Max Flow: 51946
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 15.5
Current gap: 9.68575e-007

Current alpha: 15.5
Max Flow: 51914
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 15.5
Current gap: 4.84288e-007

Current alpha: 15.5
Max Flow: 51914
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 15.5
Current gap: 2.42144e-007

Current alpha: 15.5
Max Flow: 51914
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 15.5
Current gap: 1.21072e-007

Current alpha: 15.5
Max Flow: 51914
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 15.5
Current gap: 6.0536e-008

Current alpha: 15.5
Max Flow: 51946
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 15.5
Current gap: 3.0268e-008

Current alpha: 15.5
Max Flow: 51914
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 15.5
Current gap: 1.5134e-008

Current alpha: 15.5
Max Flow: 51946
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 15.5
Current gap: 7.567e-009


Dense subgraph found with 32 vertices:
Vertices in dense subgraph: 1800 2260 3601 4321 5158 5539 5858 6063 6095 6205 7471 7472 7473 7475 7478 7480 7481 7482 7483 7484 7485 7487 7488 7489 7491 7496 7497 7552 8113 8284 8726 9301
Edges in dense subgraph:
Vertex 1800 connected to: 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 2260 connected to: 1800 1800 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 3601 connected to: 1800 1800 2260 2260 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 4321 connected to: 1800 1800 2260 2260 3601 3601 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 5158 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 5539 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 5858 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 6063 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 6095 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 6205 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7471 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7472 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7473 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7475 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7478 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7480 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7481 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7482 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7483 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7484 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7485 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7487 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7488 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7489 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7491 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7496 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7497 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7552 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 8113 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8284 8284 8726 8726 9301 9301
Vertex 8284 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8726 8726 9301 9301
Vertex 8726 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 9301 9301
Vertex 9301 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726
Number of edges: 992
Number of k-cliques in dense subgraph: 496
rho: 15.5
Execution time: 2193118 milliseconds

```

For k=3

```bash
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
Current gap: 6.02007e-005

Current alpha: 155
Max Flow: 85017
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 155
Current gap: 3.01003e-005

Current alpha: 155
Max Flow: 84985
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 155
Current gap: 1.50502e-005

Current alpha: 155
Max Flow: 85017
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 155
Current gap: 7.52509e-006

Current alpha: 155
Max Flow: 85017
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 155
Current gap: 3.76254e-006

Current alpha: 155
Max Flow: 85017
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 155
Current gap: 1.88127e-006

Current alpha: 155
Max Flow: 84985
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 155
Current gap: 9.40636e-007

Current alpha: 155
Max Flow: 85017
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 155
Current gap: 4.70318e-007

Current alpha: 155
Max Flow: 85017
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 155
Current gap: 2.35159e-007

Current alpha: 155
Max Flow: 85017
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 155
Current gap: 1.17579e-007

Current alpha: 155
Max Flow: 84985
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 155
Current gap: 5.87897e-008

Current alpha: 155
Max Flow: 84985
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 155
Current gap: 2.93949e-008

Current alpha: 155
Max Flow: 85017
Min cut (nodes in S side): Vertices in source side: No vertices in source side, updating u to 155
Current gap: 1.46974e-008

Current alpha: 155
Max Flow: 84985
Min cut (nodes in S side): Vertices in source side: Some vertices in source side, updating l to 155
Current gap: 7.34872e-009


Dense subgraph found with 32 vertices:
Vertices in dense subgraph: 1800 2260 3601 4321 5158 5539 5858 6063 6095 6205 7471 7472 7473 7475 7478 7480 7481 7482 7483 7484 7485 7487 7488 7489 7491 7496 7497 7552 8113 8284 8726 9301
Edges in dense subgraph:
Vertex 1800 connected to: 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 2260 connected to: 1800 1800 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 3601 connected to: 1800 1800 2260 2260 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 4321 connected to: 1800 1800 2260 2260 3601 3601 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 5158 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 5539 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 5858 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 6063 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 6095 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 6205 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7471 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7472 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7473 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7475 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7478 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7480 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7481 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7482 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7483 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7484 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7485 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7487 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7488 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7489 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7491 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7496 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7497 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7552 7552 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 7552 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 8113 8113 8284 8284 8726 8726 9301 9301
Vertex 8113 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8284 8284 8726 8726 9301 9301
Vertex 8284 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8726 8726 9301 9301
Vertex 8726 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 9301 9301
Vertex 9301 connected to: 1800 1800 2260 2260 3601 3601 4321 4321 5158 5158 5539 5539 5858 5858 6063 6063 6095 6095 6205 6205 7471 7471 7472 7472 7473 7473 7475 7475 7478 7478 7480 7480 7481 7481 7482 7482 7483 7483 7484 7484 7485 7485 7487 7487 7488 7488 7489 7489 7491 7491 7496 7496 7497 7497 7552 7552 8113 8113 8284 8284 8726 8726
Number of edges: 992
Number of k-cliques in dense subgraph: 4960
rho: 155
Execution time: 6632691 milliseconds

```


#### asCaida Dataset

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
rho: 17.5341
Execution time: 25074950 milliseconds


```
