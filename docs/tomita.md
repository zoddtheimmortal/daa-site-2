# The worst-case time complexity for generating all maximal cliques and computational experiments

Etsuji Tomita, Akira Tanaka, Haruhisa Takahashi

## Pseudocode

The following is the pseudocode of the algorithm

```py
procedure CLIQUES(G)
    /* Graph G = (V, E) */

begin
    /* global variable Q is to constitute a clique */
    Q := ∅;
    EXPAND(V, V)
end of CLIQUES

procedure EXPAND(SUBG, CAND)
begin
    if SUBG = ∅
        then print ("clique,")
        /* to represent that Q is a maximal clique */
    else
        u := a vertex in SUBG that maximizes | CAND ∩ Γ(u) |;
        /* let EXT_u = CAND − Γ(u); */
        /* FINI := ∅; */

        while CAND − Γ(u) ≠ ∅
        do
            q := a vertex in (CAND − Γ(u));
            print (q, ",");
            /* to represent the next statement */
            Q := Q ∪ {q};
            SUBG_q := SUBG ∩ Γ(q);
            CAND_q := CAND ∩ Γ(q);
            EXPAND(SUBG_q, CAND_q);
            CAND := CAND − {q}; /* FINI := FINI ∪ {q}; */
            print ("back,");
            /* to represent the next statement */
            Q := Q − {q}
        od
    fi
end of EXPAND
```

## Recursive Approach

The following C++ implementation is based on [Computational Techniques for Maximum Clique Problems.](https://doi.org/10.1016/j.tcs.2006.06.015)

```cpp
#include <iostream>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <chrono>
#include <unordered_map>

using namespace std;
using namespace std::chrono;

using Graph = vector<vector<int>>;

vector<int> Q;
vector<vector<int>> maximal_cliques;

vector<int> intersection(const vector<int>& a, const vector<int>& b) {
    vector<int> result;
    set_intersection(a.begin(), a.end(), b.begin(), b.end(), back_inserter(result));
    return result;
}

vector<int> difference(vector<int> a, const vector<int>& b) {
    vector<int> result;
    set_difference(a.begin(), a.end(), b.begin(), b.end(), back_inserter(result));
    return result;
}

void EXPAND(const Graph& G, vector<int>& cand, vector<int>& SUBG) {
    if (SUBG.empty()){
        if(Q.size() <= 1) {
            return;
        }
        maximal_cliques.push_back(Q);
        return;
    }

    vector<int> fini = difference(SUBG, cand);
    int pivot = -1, max_size = -1;

    for (int u : SUBG) {
        vector<int> intersection_result = intersection(cand, G[u]);
        if ((int)intersection_result.size() > max_size) {
            max_size = intersection_result.size();
            pivot = u;
        }
    }

    vector<int> ext = difference(cand, G[pivot]);

    for (int q : ext) {
        Q.push_back(q);
        vector<int> cand_q = intersection(cand, G[q]);
        vector<int> subg_q = intersection(SUBG, G[q]);

        EXPAND(G, cand_q, subg_q);

        Q.pop_back();
        cand.erase(remove(cand.begin(), cand.end(), q), cand.end());
        fini.push_back(q);
    }
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

vector<vector<int>> CLIQUES(const Graph& G) {
    Q.clear();
    maximal_cliques.clear();

    vector<int> cand(G.size()), subg(G.size());
    for (int i = 0; i < G.size(); i++) {
        cand[i] = i;
        subg[i] = i;
    }

    EXPAND(G, cand, subg);
    return maximal_cliques;
}

void printCliquesCount(vector<vector<int>> cliques, ofstream& outfile) {
    unordered_map<int, int> size_count;

    for (const auto& clique : cliques) {
        size_count[clique.size()]++;
    }

    vector<pair<int, int>> sorted_sizes(size_count.begin(), size_count.end());
    sort(sorted_sizes.begin(), sorted_sizes.end());

    int totalcnt=0;
    for (const auto& pair : sorted_sizes) {
        outfile << "Size " << pair.first << ": " << pair.second << endl;
        totalcnt+=pair.second;
    }
    outfile<<"Total number of maximal cliques: "<<totalcnt<<endl;
}



bool readGraphFromFile(const string& filename, int& n, vector<pair<int, int>>& edges) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return false;
    }

    string line;
    if (getline(file, line)) {
        istringstream iss(line);
        int m;
        if (!(iss >> n >> m)) {
            cerr << "Error: Invalid format for vertices and edges count" << endl;
            return false;
        }
    } else {
        cerr << "Error: Empty file" << endl;
        return false;
    }

    edges.clear();
    while (getline(file, line)) {
        istringstream iss(line);
        int u, v;
        if (iss >> u >> v) {
            if (u < 0 || u >= n || v < 0 || v >= n) {
                cerr << "Warning: Edge (" << u << ", " << v << ") contains invalid vertex index" << endl;
                continue;
            }
            edges.emplace_back(u, v);
        } else {
            cerr << "Warning: Invalid edge format in line: " << line << endl;
        }
    }
    file.close();
    return true;
}

int main(int argc, char* argv[]) {
    int n;
    vector<pair<int, int>> edges;

    if(argc != 3) {
        cerr << "Usage: " << argv[0] << " <input.txt> <output.txt>" << endl;
        return 1;
    }

    if (!readGraphFromFile(argv[1], n, edges)) {
        cerr << "Failed to read graph from input.txt" << endl;
        return 1;
    }

    ofstream outfile(argv[2]);
    if (!outfile.is_open()) {
        cerr << "Error: Could not open output.txt for writing" << endl;
        return 1;
    }

    Graph G = createGraph(n, edges);

    auto start_time = high_resolution_clock::now();
    vector<vector<int>> cliques = CLIQUES(G);
    auto end_time = high_resolution_clock::now();

    printCliquesCount(cliques, outfile);
    outfile << "Execution time: " << duration_cast<milliseconds>(end_time - start_time).count() << " milliseconds" << endl;

    outfile.close();
    return 0;
}

```

### Usage

Save the above code as `tomita.cpp`

Run the following code in a terminal.

`sudo` permissions are required to increase the recursion stack limit.

```bash
g++ -O3 tomita.cpp
./a.out <input_file_path> <output_file_path>
```

Output for the code will be saved in `output.txt`.
Terminal will display any error, debugging and progress statements.

### Issues

#### Stack Overflow

Since the algorithm uses recursion, larger graphs may run into stack overflow due to deep recursion. To avoid this, we increase the stack space to 512 Mb.

We do this using the `sys/resource.h` library with the following code:

```cpp
void increase_stack_size(rlim_t stack_size = 512 * 1024 * 1024) {  // 512 MB stack
    struct rlimit rl;

    int result = getrlimit(RLIMIT_STACK, &rl);
    if (result != 0) {
        cerr << "Error getting stack limit: " << strerror(errno) << endl;
        return;
    }

    if (rl.rlim_cur < stack_size) {
        rl.rlim_cur = stack_size;
        if (rl.rlim_max < rl.rlim_cur) {
            rl.rlim_max = rl.rlim_cur;  // Also increase hard limit if necessary
        }

        result = setrlimit(RLIMIT_STACK, &rl);
        if (result != 0) {
            cerr << "Error setting stack limit: " << strerror(errno) << endl;
        } else {
            cerr << "Stack size increased to " << (stack_size / (1024 * 1024)) << " MB" << endl;
        }
    }
}

```

This ensures that we don't run out of stack space, even for larger graphs.

### Time Complexity

This algorithm incorporates pruning techniques similar to those used in the Bron–Kerbosch algorithm.
The maximal cliques produced are organized in a tree-like structure.  
The worst-case time complexity of this algorithm is _O(3<sup>n/3</sup>)_ for a graph with _n_ vertices.  
This complexity is optimal with respect to _n_, as there can be up to _O(3<sup>n/3</sup>)_ maximal cliques in an _n_-vertex graph.

## Optimized Approach

```cpp
#include <iostream>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <chrono>
#include <unordered_map>

using namespace std;
using namespace std::chrono;

using Graph = vector<vector<int>>;

#define FREQ 1000
#define ll long long

vector<int> Q;
unordered_map<int,int> distro;

ll max_clique_size = 0;
ll clique_count = 0;

vector<int> intersection(const vector<int>& a, const vector<int>& b) {
    vector<int> result;
    set_intersection(a.begin(), a.end(), b.begin(), b.end(), back_inserter(result));
    return result;
}

vector<int> difference(vector<int> a, const vector<int>& b) {
    vector<int> result;
    set_difference(a.begin(), a.end(), b.begin(), b.end(), back_inserter(result));
    return result;
}

void EXPAND(const Graph& G, vector<int>& cand, vector<int>& SUBG) {
    if (SUBG.empty()){
        if(Q.size() <= 1) {
            return;
        }

        clique_count++;
        max_clique_size = max(max_clique_size, (ll)Q.size());
        distro[Q.size()]++;

        if(clique_count % FREQ == 0) {
            cerr << "Found Cliques: " << clique_count << endl;
        }
        return;
    }

    vector<int> fini = difference(SUBG, cand);
    int pivot = -1, max_size = -1;

    for (int u : SUBG) {
        vector<int> intersection_result = intersection(cand, G[u]);
        if ((int)intersection_result.size() > max_size) {
            max_size = intersection_result.size();
            pivot = u;
        }
    }

    vector<int> ext = difference(cand, G[pivot]);

    for (int q : ext) {
        Q.push_back(q);
        vector<int> cand_q = intersection(cand, G[q]);
        vector<int> subg_q = intersection(SUBG, G[q]);

        EXPAND(G, cand_q, subg_q);

        Q.pop_back();
        cand.erase(remove(cand.begin(), cand.end(), q), cand.end());
        fini.push_back(q);
    }
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

void CLIQUES(const Graph& G) {
    Q.clear();

    vector<int> cand(G.size()), subg(G.size());
    for (int i = 0; i < G.size(); i++) {
        cand[i] = i;
        subg[i] = i;
    }

    EXPAND(G, cand, subg);
}

bool readGraphFromFile(const string& filename, int& n, vector<pair<int, int>>& edges) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return false;
    }

    string line;
    if (getline(file, line)) {
        istringstream iss(line);
        int m;
        if (!(iss >> n >> m)) {
            cerr << "Error: Invalid format for vertices and edges count" << endl;
            return false;
        }
    } else {
        cerr << "Error: Empty file" << endl;
        return false;
    }

    edges.clear();
    while (getline(file, line)) {
        istringstream iss(line);
        int u, v;
        if (iss >> u >> v) {
            if (u < 0 || u >= n || v < 0 || v >= n) {
                cerr << "Warning: Edge (" << u << ", " << v << ") contains invalid vertex index" << endl;
                continue;
            }
            edges.emplace_back(u, v);
        } else {
            cerr << "Warning: Invalid edge format in line: " << line << endl;
        }
    }
    file.close();
    return true;
}

int main(int argc, char* argv[]) {
    int n;
    vector<pair<int, int>> edges;

    if(argc != 3) {
        cerr << "Usage: " << argv[0] << " <input.txt> <output.txt>" << endl;
        return 1;
    }

    if (!readGraphFromFile(argv[1], n, edges)) {
        cerr << "Failed to read graph from input.txt" << endl;
        return 1;
    }

    ofstream outfile(argv[2]);
    if (!outfile.is_open()) {
        cerr << "Error: Could not open output.txt for writing" << endl;
        return 1;
    }

    Graph G = createGraph(n, edges);

    auto start_time = high_resolution_clock::now();
    CLIQUES(G);
    auto end_time = high_resolution_clock::now();

    outfile << "Execution time: " << duration_cast<milliseconds>(end_time - start_time).count() << " milliseconds" << endl;
    outfile << "Maximal clique size: " << max_clique_size << endl;
    outfile << "Clique count: " << clique_count << endl;
    outfile << "Clique distribution:" << endl;

    for (const auto& entry : distro) {
        outfile << "Size " << entry.first << ": " << entry.second << endl;
    }

    outfile.close();
    return 0;
}
```

### Usage

Save the above code as `tomita.cpp`

Run the following code in a terminal.

`sudo` permissions are required to increase the recursion stack limit.

```bash
g++ -O3 tomita.cpp
./a.out <input_file_path> <output_file_path>
```

Output for the code will be saved in `output.txt`.
Terminal will display any error, debugging and progress statements.

## Results

### Wikipedia Vote Network Dataset

```bash
Execution time: 1215 milliseconds
Maximal clique size: 17
Clique count: 459002
Clique distribution:
Size 17: 23
Size 16: 208
Size 15: 740
Size 5: 48416
Size 6: 68872
Size 4: 27292
Size 7: 83266
Size 8: 76732
Size 3: 13718
Size 9: 54456
Size 10: 35470
Size 2: 8655
Size 11: 21736
Size 12: 11640
Size 13: 5449
Size 14: 2329
```

### Enron Email Network Dataset

```bash
Execution time: 1727 milliseconds
Maximal clique size: 20
Clique count: 226859
Clique distribution:
Size 20: 6
Size 18: 41
Size 19: 10
Size 17: 286
Size 16: 1178
Size 15: 3157
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
```

### As Skitter Network Dataset

```bash
Execution time: 1820871 milliseconds
Maximal clique size: 67
Clique count: 37322355
Clique distribution:
Size 59: 447
Size 60: 405
Size 61: 283
Size 64: 84
Size 62: 242
Size 66: 22
Size 67: 4
Size 24: 706753
Size 23: 673924
Size 25: 753633
Size 27: 892719
Size 26: 818353
Size 28: 955212
Size 21: 611976
Size 22: 640890
Size 20: 600192
Size 19: 639413
Size 18: 729601
Size 17: 839330
Size 16: 939987
Size 15: 980831
Size 2: 2319807
Size 3: 3171609
Size 33: 946717
Size 4: 1823321
Size 63: 146
Size 34: 878552
Size 5: 939336
Size 7: 598284
Size 6: 684873
Size 65: 49
Size 10: 665661
Size 8: 588889
Size 9: 608937
Size 12: 798073
Size 11: 728098
Size 13: 877282
Size 14: 945194
Size 30: 1034106
Size 31: 1055653
Size 29: 999860
Size 32: 1017560
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
Size 46: 158352
Size 45: 222461
Size 47: 99522
Size 48: 62437
Size 49: 39822
Size 53: 9514
Size 52: 17707
Size 51: 25637
Size 50: 30011
Size 54: 3737
Size 56: 1080
Size 55: 2042
Size 57: 546
Size 58: 449
```
