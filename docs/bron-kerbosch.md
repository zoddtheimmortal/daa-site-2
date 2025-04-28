# Listing All Maximal Cliques in Sparse Graphs in Near-optimal Time

David Eppstein, Maarten Loffler, and Darren Strash

## Pseudocode

The following is the pseudocode of the algorithm

### Bron-Kerbosch without Pivot

```py

def bron_kerbosch(P, R, X):
    if not P and not X:
        print("Maximal Clique:", R)
        return

    for v in list(P):  # Iterate over a copy to avoid modification issues
        bron_kerbosch(P & neighbors(v), R | {v}, X & neighbors(v))
        P.remove(v)
        X.add(v)


```

### Bron-Kerbosch with Pivoting

```py
def bron_kerbosch_pivot(P, R, X):
    if not P and not X:
        print("Maximal Clique:", R)
        return

    u = select_pivot(P | X)  # Choose pivot vertex u to maximize |P ∩ Γ(u)|
    for v in list(P - neighbors(u)):  # Consider vertices not in Γ(u)
        bron_kerbosch_pivot(P & neighbors(v), R | {v}, X & neighbors(v))
        P.remove(v)
        X.add(v)
```

### Bron-Kerbosch with Degeneracy Ordering

```py
def bron_kerbosch_degeneracy(V, E):
    ordered_vertices = degeneracy_ordering(V, E)  # Compute degeneracy ordering

    for i, v in enumerate(ordered_vertices):
        P = {u for u in neighbors(v) if ordered_vertices.index(u) > i}
        X = {u for u in neighbors(v) if ordered_vertices.index(u) < i}
        bron_kerbosch_pivot(P, {v}, X)
```

## Recursive Approach

The following C++ implementation is based on [Efficient Enumeration of Maximal Cliques in Sparse Graphs.](https://arxiv.org/abs/1006.5440)

```cpp
#include <bits/stdc++.h>
#include<iostream>
#include<fstream>
#include<sstream>
#include <chrono>
using namespace std::chrono;
using namespace std;
map<int,int> maximal_cliques;



struct pair_hash {
    template <class T1, class T2>
    size_t operator()(const pair<T1, T2>& p) const {
        return hash<T1>()(p.first) ^ hash<T2>()(p.second);
    }
};

unordered_set<int> find_union(const unordered_set<int>& R, int node) {
    unordered_set<int> result = R;
    result.insert(node);
    return result;
}

// Function to find union of two sets
unordered_set<int> find_union_set(const unordered_set<int>& R, const unordered_set<int>& X) {
    unordered_set<int> result = R;
    for (auto it : X) {
        result.insert(it);
    }
    return result;
}

// Function to find intersection of a set and neighbors
unordered_set<int> find_intersect(const unordered_set<int>& R, const vector<int>& neighbours) {
    unordered_set<int> result;
    for (auto it : neighbours) {
        if (R.find(it) != R.end())
            result.insert(it);
    }
    return result;
}
// Find set difference
unordered_set<int> find_difference(const unordered_set<int>& P, const vector<int>& neighbours) {
    unordered_set<int> difference = P;
    for (auto it : neighbours) {
        difference.erase(it);
    }
    return difference;
}


// Forward declaration
void bron_kerbosch_pivot(set<int> P, set<int> R, set<int> X,
                        vector<vector<int>>&edges);

// Find pivot to optimize Bron-Kerbosch
int find_pivot(const unordered_set<int>& P,const unordered_set<int>&R, const unordered_set<int>& X, const vector<vector<int>>& edges) {
    // Create union of P and X
    unordered_set<int> union_set;
    union_set.insert(P.begin(), P.end());
    union_set.insert(X.begin(), X.end());

    int max_connections = -1;
    int pivot = -1;

    // Find the vertex with the most connections to vertices in P
    for (int u : union_set) {
        int connections = 0;

        // Count how many vertices in P are adjacent to u
        for (int neighbor : edges[u]) {
            if (P.find(neighbor) != P.end()) {
                connections++;
            }
        }

        // Update pivot if current vertex has more connections
        if (connections > max_connections) {
            max_connections = connections;
            pivot = u;
        }
    }

    return pivot;
}


// Bron-Kerbosch with pivot
void bron_kerbosch_pivot(unordered_set<int> P, unordered_set<int> R, unordered_set<int> X, vector<vector<int>>&edges) {

    unordered_set<int>P_union_X=find_union_set(P,X);
    if(P_union_X.empty())
    {
        int cliqueSize=R.size();
        maximal_cliques[cliqueSize]++;
        return ;
    }

    // Use pivot to reduce recursive calls
    int pivot = find_pivot(P,R, X, edges);
    unordered_set<int> p_diff_neighbours = find_difference(P, edges[pivot]);

    for (auto node : p_diff_neighbours) {
        unordered_set<int> new_P = find_intersect(P, edges[node]);
        unordered_set<int> new_R = find_union(R, node);
        unordered_set<int> new_X = find_intersect(X, edges[node]);

        bron_kerbosch_pivot(new_P, new_R, new_X, edges);

        P.erase(node);
        X.insert(node);
    }
}


// Generate degeneracy ordering
void find_order(vector<int>& ans, vector<int>& degree,
                vector<vector<int>>&edges,
                vector<int>&vis) {
    for(int k=0;k<edges.size();k++)
    {
        int u=-1;
        int minDeg=INT_MAX;
        for(int i=0;i<edges.size();i++)
        {
            if(!vis[i]&&degree[i]<minDeg)
            {
                minDeg=degree[i];
                u=i;
            }
        }
        if(u==-1)
            break;
        vis[u]=1;
        ans.push_back(u);
        for(int w:edges[u])
        {
            if(!vis[w])
                degree[w]--;

        }
    }
    reverse(ans.begin(),ans.end());

}

// Modified Bron-Kerbosch with degeneracy ordering
void bron_kerbosch_modified(
                           vector<vector<int>>&edges,
                           const vector<int>& order,vector<int>&pos) {
    int nodes=edges.size();
    for(int i=0;i<nodes;i++)
    {
        int node_ini=order[i];
        unordered_set<int>P;
        for(int w:edges[node_ini])
        {
            if(pos[w]>pos[node_ini]) P.insert(w);

        }
        unordered_set<int> X;
        for(int w:edges[node_ini]){
            if(pos[w]<pos[node_ini]) X.insert(w);
        }
        unordered_set<int> R;
        R.insert(node_ini);
        bron_kerbosch_pivot(P,R,X,edges);
    }



}

// Main function to find maximal cliques
void  bron_kerbosch(vector<vector<int>>&edges) {
    int nodes=edges.size();
    vector<int>degree(nodes,0);
    vector<int>vis(nodes,0);
    set<int> all_nodes;

    for(int i=0;i<nodes;i++)
    {
        degree[i]=edges[i].size();
    }

    vector<int> deg_order;
    vector<int>pos(nodes);
    find_order(deg_order, degree, edges, vis);
    for(int i=0;i<nodes;i++) pos[deg_order[i]] = i;

    bron_kerbosch_modified(edges, deg_order,pos);


}


int main(int argc, char* argv[]) {
    ifstream infile(argv[1]);
    ofstream outfile(argv[2]);
    string line;
    unordered_set<pair<int,int>,pair_hash> edgeList;
    set<int>nodes;
    int maxVertex = 0;
    printf("Reading dataset...........\n");
    while(getline(infile, line)){
        if(line.empty() || line[0] == '#')
            continue;
        istringstream iss(line);
        int u, v;
        if(iss >> u >> v){
            if(u!=v){
            edgeList.insert({u, v});
            edgeList.insert({v,u});
            nodes.insert(u);
            nodes.insert(v);
            }
        }
    }
    infile.close();
    printf("Mapping down edges to create edges vector....\n");
    map<int, int> nodeMapping;
    int newId = 0;
    for (int node : nodes) {
        nodeMapping[node] = newId++;
    }

    // Update edge list with new mapped values
    vector<vector<int>> newEdgeList(newId);
    for (auto& edge : edgeList) {
        int newU = nodeMapping[edge.first];
        int newV = nodeMapping[edge.second];
        newEdgeList[newU].push_back(newV);
    }

    printf("Entered the bron_kerbosch function......\n");
    auto start = high_resolution_clock::now();
    bron_kerbosch(newEdgeList);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    printf("Execution completed, printing results........\n");
    outfile<<"# For the dataset:"<<argv[1]<<endl;
    outfile<<"Execution time(ms):   "<<duration.count()<<endl;
    long long int count=0;
    for(auto it :maximal_cliques)
    {
        count+=it.second;
        outfile<<"No of cliques for Size:"<<it.first<<" are :"<<it.second<<endl;
    }
    outfile<<"Total number of maximal cliques:  "<<count<<endl;
    printf("Program executed successfully!!!\n");
    outfile.close();


    return 0;
}
```

### Usage

Save the above code as `bron-kerbosch.cpp`

Run the following code in a terminal.

```bash
g++ -O3 bron-kerbosch.cpp
./a.out <input_file_path> <output_file_path>
```

Output for the code will be saved in `<output_file_path>`.
Terminal will display any error, debugging and progress statements.

### Results

#### Wiki vote dataset

The algorithm takes 2.21 mins to run on the Wiki Vote dataset.

```bash
Execution time(ms):   132813
No of cliques for Size:2 are :8655
No of cliques for Size:3 are :13718
No of cliques for Size:4 are :27292
No of cliques for Size:5 are :48416
No of cliques for Size:6 are :68872
No of cliques for Size:7 are :83266
No of cliques for Size:8 are :76732
No of cliques for Size:9 are :54456
No of cliques for Size:10 are :35470
No of cliques for Size:11 are :21736
No of cliques for Size:12 are :11640
No of cliques for Size:13 are :5449
No of cliques for Size:14 are :2329
No of cliques for Size:15 are :740
No of cliques for Size:16 are :208
No of cliques for Size:17 are :23


```

#### Email-Enron dataset

The algorithm takes 2.59 mins to run on the Email Enron dataset.

```bash

Execution time(ms):   155667
No of cliques for Size:2 are :14070
No of cliques for Size:3 are :7077
No of cliques for Size:4 are :13319
No of cliques for Size:5 are :18143
No of cliques for Size:6 are :22715
No of cliques for Size:7 are :25896
No of cliques for Size:8 are :24766
No of cliques for Size:9 are :22884
No of cliques for Size:10 are :21393
No of cliques for Size:11 are :17833
No of cliques for Size:12 are :15181
No of cliques for Size:13 are :11487
No of cliques for Size:14 are :7417
No of cliques for Size:15 are :3157
No of cliques for Size:16 are :1178
No of cliques for Size:17 are :286
No of cliques for Size:18 are :41
No of cliques for Size:19 are :10
No of cliques for Size:20 are :6
```

### Issues

#### Stack Overflow

1.Since the algorithm uses recursion, larger graphs may run into stack overflow due to deep recursion. To avoid this, we increase the stack space to 512 Mb.

We do this using the following code:

VirtualAlloc is best for allocating large stack space dynamically.

```cpp
#include <windows.h>
#include <iostream>

void increase_stack_size(SIZE_T stack_size = 512 * 1024 * 1024) { // 512 MB
    LPVOID stack = VirtualAlloc(nullptr, stack_size, MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);

    if (stack == nullptr) {
        std::cerr << "Error increasing stack size: " << GetLastError() << std::endl;
    } else {
        std::cout << "Stack size increased to " << (stack_size / (1024 * 1024)) << " MB" << std::endl;
    }
}

```

This ensures that we don't run out of stack space, even for larger graphs.

2.Also due to deep recursions on larger dataset the Algorithm may run into SEGMENTATION fault , Thus the code requires more optimizations.

## Optimized Approach

(Using inline function and Mapping down nodes to continuous set of nodes stored in vectors)

### Code

```cpp
#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
using namespace std::chrono;
using namespace std;

#define ll long long
#define FREQ 1000

ll clique_count = 0;
vector<long long> clique_counts;
long long max_clique_size = 0;

vector<int> deg_order;
vector<int> pos;

typedef vector<int> VertexSet;


inline bool is_neighbor(int u, int v, const vector<vector<int>>& edges) {
    const vector<int>& neighbors = edges[u];
    return binary_search(neighbors.begin(), neighbors.end(), v);
}


inline void set_intersection_with_adj(const vector<int>& adj, const VertexSet& vec, VertexSet& result) {
    for (int v : adj) {
        if (binary_search(vec.begin(), vec.end(), v)) {
            result.push_back(v);
        }
    }
}


inline int find_pivot_sorted(const VertexSet& P, const VertexSet& X,
                           const vector<vector<int>>& edges) {
    int best_pivot = -1;
    int max_connections = -1;


    for (int u : X) {
        int connections = 0;
        for (int neighbor : edges[u]) {
            if (binary_search(P.begin(), P.end(), neighbor)) {
                connections++;
            }
        }

        if (connections > max_connections) {
            max_connections = connections;
            best_pivot = u;
        }
    }


    for (int u : P) {
        int connections = 0;
        for (int neighbor : edges[u]) {
            if (binary_search(P.begin(), P.end(), neighbor)) {
                connections++;
            }
        }

        if (connections > max_connections) {
            max_connections = connections;
            best_pivot = u;
        }
    }

    return best_pivot == -1 ? (P.empty() ? X[0] : P[0]) : best_pivot;
}


inline void record_clique(size_t size) {
    if (size >= clique_counts.size()) {
        clique_counts.resize(size + 100, 0);
    }
    clique_counts[size]++;
    max_clique_size = max(max_clique_size, (long long)size);
    clique_count++;
    if(clique_count % FREQ == 0) {
        cerr << "Found " << clique_count << " cliques so far." << endl;
    }
}


inline void bron_kerbosch_pivot(VertexSet& P, VertexSet& R, VertexSet& X,
                           const vector<vector<int>>& edges) {
    if (P.empty() && X.empty()) {
        record_clique(R.size());
        return;
    }


    int pivot = find_pivot_sorted(P, X, edges);


    VertexSet P_minus_N_pivot;
    P_minus_N_pivot.reserve(P.size());


    vector<bool> is_pivot_neighbor(edges.size(), false);
    for (int neighbor : edges[pivot]) {
        is_pivot_neighbor[neighbor] = true;
    }

    for (int v : P) {
        if (!is_pivot_neighbor[v]) {
            P_minus_N_pivot.push_back(v);
        }
    }


    VertexSet new_P, new_X, new_R;
    new_P.reserve(P.size());
    new_X.reserve(X.size());
    new_R.reserve(R.size() + 1);

    for (int v : P_minus_N_pivot) {

        new_R = R;
        new_R.push_back(v);


        new_P.clear();
        set_intersection_with_adj(edges[v], P, new_P);
        sort(new_P.begin(), new_P.end());


        new_X.clear();
        set_intersection_with_adj(edges[v], X, new_X);
        sort(new_X.begin(), new_X.end());


        bron_kerbosch_pivot(new_P, new_R, new_X, edges);


        auto it = lower_bound(P.begin(), P.end(), v);
        if (it != P.end() && *it == v) {
            P.erase(it);
        }

        X.push_back(v);
        sort(X.begin(), X.end());
    }
}


void find_order(vector<int>& ans, vector<int>& degree,
               vector<vector<int>>& edges,
               vector<int>& vis) {
    int n = edges.size();
    int max_degree = 0;


    for (int i = 0; i < n; i++) {
        max_degree = max(max_degree, degree[i]);
    }


    vector<list<int>> buckets(max_degree + 1);
    vector<list<int>::iterator> pos(n);


    for (int i = 0; i < n; i++) {
        buckets[degree[i]].push_back(i);
        pos[i] = --buckets[degree[i]].end();
    }

    for (int k = 0; k < n; k++) {

        int d = 0;
        while (d <= max_degree && buckets[d].empty()) d++;

        if (d > max_degree) break;


        int node = buckets[d].front();
        buckets[d].pop_front();

        vis[node] = 1;
        ans.push_back(node);


        for (int neighbor : edges[node]) {
            if (!vis[neighbor]) {

                buckets[degree[neighbor]].erase(pos[neighbor]);


                degree[neighbor]--;


                buckets[degree[neighbor]].push_back(neighbor);
                pos[neighbor] = --buckets[degree[neighbor]].end();
            }
        }
    }

    reverse(ans.begin(), ans.end());
}


inline void bron_kerbosch_modified(const vector<vector<int>>& edges,
                                 const vector<int>& order, const vector<int>& pos) {
    int n = edges.size();


    VertexSet P, X, R;
    P.reserve(n);
    X.reserve(n);
    R.reserve(n);

    for (int i = 0; i < n; i++) {
        int v = order[i];
        int v_pos = pos[v];


        P.clear();
        X.clear();
        R.clear();

        R.push_back(v);


        for (int neighbor : edges[v]) {
            int neighbor_pos = pos[neighbor];
            if (neighbor_pos > v_pos) {
                P.push_back(neighbor);
            } else if (neighbor_pos < v_pos) {
                X.push_back(neighbor);
            }
        }


        sort(P.begin(), P.end());
        sort(X.begin(), X.end());

        bron_kerbosch_pivot(P, R, X, edges);
    }
}

inline void bron_kerbosch(const vector<vector<int>>& edges) {
    bron_kerbosch_modified(edges, deg_order, pos);
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <input_file> <output_file>" << endl;
        return 1;
    }

    ifstream infile(argv[1]);
    if (!infile) {
        cerr << "Error: Cannot open input file " << argv[1] << endl;
        return 1;
    }

    ofstream outfile(argv[2]);
    if (!outfile) {
        cerr << "Error: Cannot open output file " << argv[2] << endl;
        return 1;
    }

    string line;
    vector<pair<int, int>> rawEdges;
    int maxNodeId = -1;

    cerr << "Reading input file..." << endl;
    while (getline(infile, line)) {
        if (line.empty() || line[0] == '#')
            continue;

        istringstream iss(line);
        int u, v;
        if (iss >> u >> v) {
            if (u != v) {
                rawEdges.push_back({u, v});
                rawEdges.push_back({v, u});
                maxNodeId = max(maxNodeId, max(u, v));
            }
        }
    }
    infile.close();

    cerr << "Finished reading input file." << endl;
    cerr << "Number of edges: " << rawEdges.size() / 2 << endl;

    if (maxNodeId == -1) {
        outfile << "# For the dataset: " << argv[1] << endl;
        outfile << "Execution time(ms): 0" << endl;
        outfile << "Total number of maximal cliques: 0" << endl;
        outfile.close();
        return 0;
    }


    vector<bool> nodeExists(maxNodeId + 1, false);
    for (auto& edge : rawEdges) {
        nodeExists[edge.first] = true;
        nodeExists[edge.second] = true;
    }


    int nodeCount = 0;
    vector<int> nodeMapping(maxNodeId + 1, -1);
    for (int i = 0; i <= maxNodeId; i++) {
        if (nodeExists[i]) {
            nodeMapping[i] = nodeCount++;
        }
    }

    cerr << "Graph size: " << nodeCount << " nodes" << endl;


    vector<vector<int>> newEdgeList(nodeCount);
    for (auto& edge : rawEdges) {
        int newU = nodeMapping[edge.first];
        int newV = nodeMapping[edge.second];
        newEdgeList[newU].push_back(newV);
    }


    for (int i = 0; i < nodeCount; i++) {
        auto& neighbors = newEdgeList[i];
        if (!neighbors.empty()) {

            sort(neighbors.begin(), neighbors.end());
            neighbors.erase(unique(neighbors.begin(), neighbors.end()), neighbors.end());
        }
    }

    int v = nodeCount;
    vector<int> degree(v, 0);
    vector<int> vis(v, 0);
    pos.resize(v);
    deg_order.clear();

    for (int i = 0; i < v; i++) {
        degree[i] = newEdgeList[i].size();
    }


    clique_counts.resize(1000, 0);

    find_order(deg_order, degree, newEdgeList, vis);
    for (int i = 0; i < v; i++) pos[deg_order[i]] = i;

    cerr << "Starting Bron-Kerbosch algorithm..." << endl;
    auto start = high_resolution_clock::now();
    bron_kerbosch(newEdgeList);
    auto stop = high_resolution_clock::now();
    cerr << "Algorithm completed." << endl;

    auto duration = duration_cast<milliseconds>(stop - start);

    outfile << "# For the dataset: " << argv[1] << endl;
    outfile << "Execution time(ms): " << duration.count() << endl;

    long long total_count = 0;
    for (int i = 1; i <= max_clique_size; i++) {
        if (i < clique_counts.size() && clique_counts[i] > 0) {
            outfile << "Size " << i << ": " << clique_counts[i] << endl;
            total_count += clique_counts[i];
        }
    }

    outfile << "Total number of maximal cliques: " << total_count << endl;
    outfile << "Maximal clique size: " << max_clique_size << endl;
    outfile.close();

    cerr << "Results written to output file." << endl;
    return 0;
}
```

### Optimizations

1. We used Inline functions into our optimized code , this avoids unnecessary copy creation of large sets and edge-list
2. We passed parameters by reference hence optimizing large number of function calling.
3. We Mapped down the random nodes number into a continuous set of nodes , so that we can avoid unnecssary complexities of using unordered-map for storing edges and instead using vector with O(1) constant T.C. for fetching adjacency list for specific nodes.

### Usage

Save the above code as `bron-opti.cpp`

Run the following code in a terminal.

```bash
g++ -O3 bron-opti.cpp
./a.out <input_file_path> <output_file_path>
```

Output for the code will be saved in `<output_file_path>`.
Terminal will display any error, debugging and progress statements.

### Results

#### Wiki vote dataset

The algorithm takes 0.13 mins to run on the Wiki Vote dataset.

```bash
Execution time(ms): 7546
Size 2: 8655
Size 3: 13718
Size 4: 27292
Size 5: 48416
Size 6: 68872
Size 7: 83266
Size 8: 76732
Size 9: 54456
Size 10: 35470
Size 11: 21736
Size 12: 11640
Size 13: 5449
Size 14: 2329
Size 15: 740
Size 16: 208
Size 17: 23
Total number of maximal cliques: 459002
Maximal clique size: 17

```

#### Email-Enron dataset

The algorithm takes 0.1 mins to run on the Email Enron dataset.

```bash

Execution time(ms): 6208
Size 2: 14070
Size 3: 7077
Size 4: 13319
Size 5: 18143
Size 6: 22715
Size 7: 25896
Size 8: 24766
Size 9: 22884
Size 10: 21393
Size 11: 17833
Size 12: 15181
Size 13: 11487
Size 14: 7417
Size 15: 3157
Size 16: 1178
Size 17: 286
Size 18: 41
Size 19: 10
Size 20: 6
Total number of maximal cliques: 226859
Maximal clique size: 20



```

#### Skitter dataset

The algorithm takes 63.51 mins to run on the Email Enron dataset.

```bash
Execution time(ms): 3810505
Size 2: 2319807
Size 3: 3171609
Size 4: 1823321
Size 5: 939336
Size 6: 684873
Size 7: 598284
Size 8: 588889
Size 9: 608937
Size 10: 665661
Size 11: 728098
Size 12: 798073
Size 13: 877282
Size 14: 945194
Size 15: 980831
Size 16: 939987
Size 17: 839330
Size 18: 729601
Size 19: 639413
Size 20: 600192
Size 21: 611976
Size 22: 640890
Size 23: 673924
Size 24: 706753
Size 25: 753633
Size 26: 818353
Size 27: 892719
Size 28: 955212
Size 29: 999860
Size 30: 1034106
Size 31: 1055653
Size 32: 1017560
Size 33: 946717
Size 34: 878552
Size 35: 809485
Size 36: 744634
Size 37: 663650
Size 38: 583922
Size 39: 520239
Size 40: 474301
Size 41: 420796
Size 42: 367879
Size 43: 321829
Size 44: 275995
Size 45: 222461
Size 46: 158352
Size 47: 99522
Size 48: 62437
Size 49: 39822
Size 50: 30011
Size 51: 25637
Size 52: 17707
Size 53: 9514
Size 54: 3737
Size 55: 2042
Size 56: 1080
Size 57: 546
Size 58: 449
Size 59: 447
Size 60: 405
Size 61: 283
Size 62: 242
Size 63: 146
Size 64: 84
Size 65: 49
Size 66: 22
Size 67: 4
Total number of maximal cliques: 37322355
Maximal clique size: 67
```

### Improvement

1. For Wiki-Vote dataset we can observe a time reduction of 94.12%  by using optimized approach
2. For Email-Enron dataset we can observe a time reduction of 96.14% by using optimized approach
3. For skitter dataser we can observe that we are able to find results without running into SEGMENTATION faults
