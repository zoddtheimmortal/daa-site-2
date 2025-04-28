# Arboricity and Subgraph Listing Algorithms

Norishige Chiba and Takao Nishizeki

The following C++ implementations are based on [Chiba & Nishizeki (1985)](https://doi.org/10.1137/0214017).

## Pseudo Code

```py
procedure UPDATE (i, C);
begin
  if i = n + 1 then
    print out a new clique C
  else
    begin
      if C - N(i) ≠ ∅ then UPDATE (i + 1, C);

      { Prepare for tests }
      { Compute T[y] = |N(y) ∩ C ∩ N(i)| for y ∈ V - C - {i} }
      for each vertex x ∈ C ∩ N(i)
        do for each vertex y ∈ N(x) - C - {i}
          do T[y] := T[y] + 1;

      { Compute S[y] = |N(y) ∩ (C - N(i))| for y ∈ V - C }
      for each vertex x ∈ C - N(i)
        do for each vertex y ∈ N(x) - C
          do S[y] := S[y] + 1;

      FLAG := true;

      { Maximality test }
      if there exists a vertex y ∈ N(i) - C such that y < i and T[y] = |C ∩ N(i)|
        then FLAG := false; { (C ∩ N(i)) ∪ {i} is not a clique of G }

      { Lexicographic test }
      { C ∩ N(i) corresponds to C_i in Lemma 6 }
      sort all the vertices in C - N(i) in ascending order j₁ < j₂ < ... < jₚ,
      where p = |C - N(i)|;

      for k := 1 to p
        do for each vertex y ∈ N(jₖ) - C such that y < i and T[y] = |C ∩ N(i)|
          do
            if y = jₖ then
              S[y] := S[y] - 1 { Alter S[y] to S(y) }
            	else
		   if (jₖ is the first vertex which satisfies y < jₖ)
                     then { S[y] := S(jₖ); }
                     if (S[y] + k = p) and (y ≥ jₖ₋₁) and (j0 = 0)
                        then FLAG := false; { C is not lexicographically largest }

       if C ∩ N(i) ≠ ∅
        then for each vertex y ∉ C ∪ {i} such that y < i, T[y] = |C ∩ N(i)| and S[y] = 0
          do begin
            { Access y from the adjacency list of a vertex in C ∩ N(i) }
            if jₚ < y then FLAG := false; { C is not lexicographically largest }
          end;
       else if jₚ < i - 1 then FLAG := false; { C is not lexicographically largest }

      { Reinitialize S and T }
      for each vertex x ∈ C ∩ N(i)
        do for each vertex y ∈ N(x) - C - {i}
          do T[y] := 0;

      for each vertex x ∈ C - N(i)
        do for each vertex y ∈ N(x) - C
          do S[y] := 0;

      { FLAG is true if and only if (C ∩ N(i)) ∪ {i} is a clique of G, and C is the
      lexicographically largest clique of G_i containing C ∩ N(i). }

      if FLAG then
        begin
          SAVE := C - N(i);
          C := (C ∩ N(i)) ∪ {i};
          UPDATE (i + 1, C);
          C := (C - {i}) ∪ SAVE;
        end;
    end;
end;

begin { of CLIQUE }
  number the vertices of a given graph G in such a way that
  d(1) ≤ d(2) ≤ ... ≤ d(n);

  for i := 1 to n { Initialize S and T }
    do begin
      S[i] := 0;
      T[i] := 0;
    end;

  C := {1};
  UPDATE(2, C);
end { of CLIQUE };

```

## Usage

Save any of the following code as `chiba.cpp`.

Run the following code in a terminal.

`sudo` permissions are required to increase the recursion stack limit.

```bash
g++ -O3 chiba.cpp
sudo ./a.out <path_for_input_file>  <path_for_output_file>
```

We use the `-O3 flag` as it enables aggressive optimizations that can significantly improve the performance of the program

Output for the code will be saved in `path_for_output_file`.

Terminal will display errors, debugging and progress statements.

## Naive Approach

The Naive Approach is a one-to-one implementation of the research paper. No additional optimizations have been added.

### Code

```cpp
#include <bits/stdc++.h>
#include <sys/resource.h>
using namespace std;

using ll = long long;
using vi = vector<int>;
using vl = vector<ll>;

#define all(x) x.begin(), x.end()
#define fr(i, a, b) for (ll i = a; i < (b); ++i)
#define rf(i, a, b) for (ll i = b; i >=(a); i--)
#define nL "\n"
#define fast_io ios_base::sync_with_stdio(false);cin.tie(nullptr)

#define FREQ 1000

ll n,m;
vector<set<int>> adj;
vector<int> S,T;

ll max_clique_size = 0;
ll clique_count = 0;
unordered_map<ll,ll> distribution;

void UPDATE(int i, set<int>&C){
    if(i>n) return;
    if(i==n){
        clique_count++;
        distribution[C.size()]++;
        max_clique_size=max(max_clique_size,(ll)C.size());

        if(clique_count%FREQ==0){
            cerr<<"Found: "<<clique_count<<", Max Clique Size: "<<max_clique_size<<nL;
        }
        return;
    }

    set<int> intersection;
    set_intersection(C.begin(),C.end(),adj[i].begin(),adj[i].end(),
    inserter(intersection, intersection.begin()));

    set<int> difference;
    set_difference(C.begin(),C.end(),adj[i].begin(),adj[i].end(),
    inserter(difference, difference.begin()));

    if(!difference.empty()){
        UPDATE(i+1,C);
    }

    for(int x:intersection){
        for(int y:adj[x]){
            if(C.find(y)==C.end() && y!=i){
                T[y]++;
            }
        }
    }

    for(int x:difference){
        for(int y:adj[x]){
            if(C.find(y)==C.end()){
                S[y]++;
            }
        }
    }

    bool FLAG=true;

    for(int y:adj[i]){
        if(C.find(y)==C.end() && y<i && T[y]==intersection.size()){
            FLAG=false;
        }
    }

    vector<int> diff_idx(difference.begin(), difference.end());
    sort(diff_idx.begin(), diff_idx.end());

    int p=diff_idx.size();
    fr(k,1,p+1){
        int jk=diff_idx[k-1];
        bool first_handled=false;

        for(int y:adj[jk]){

            if(C.find(y)==C.end()&&y<i&&T[y]==intersection.size()){
                if(y>=jk){
                    S[y]--;
                }
                else{
                    if(!first_handled){
                        if(S[y]+k-1==p&&((k-1==0&&y>=0)||(y>=diff_idx[k-2]))){
                            FLAG=false;
                            first_handled=true;
                        }

                    }
                }
            }
        }
    }

    int jp=(!diff_idx.empty())?diff_idx.back():0;
    if(!intersection.empty()){
        fr(y,0,i){
            if(C.find(y)==C.end()&&y!=i&&
            T[y]==intersection.size()&&S[y]==0){
                if(jp<y){
                    FLAG=false;
                    break;
                }
            }
        }
    }
    else if(jp<i-1){
        FLAG=false;
    }

    for(int x:intersection){
        for(int y:adj[x]){
            if(C.find(y)==C.end()&&y!=i){
                T[y]=0;
            }
        }
    }

    for(int x:difference){
        for(int y:adj[x]){
            if(C.find(y)==C.end()){
                S[y]=0;
            }
        }
    }

    if(FLAG){
        set<int> SAVE;
        for(int x:difference){
            SAVE.insert(x);
        }

        C.clear();
        for(int x:intersection){
            C.insert(x);
        }
        C.insert(i);

        UPDATE(i+1,C);

        C.erase(i);
        for(int x:SAVE){
            C.insert(x);
        }
    }
}

void CLIQUE(){
    int start_vertex = 0;
    while(adj[start_vertex].empty()){
        start_vertex++;
    }

    S.assign(n,0);
    T.assign(n,0);

    set<int> C = {start_vertex};
    UPDATE(start_vertex+1,C);
}

void solve(){
    cin>>n>>m;
    adj.resize(n);

    cerr<<"Input Size: "<<n<<", Edges: "<<m<<nL;
    cerr<<"Reading Input..."<<nL;
    fr(i,0,m){
        ll f,s;
        cin>>f>>s;
        adj[f].insert(s);
        adj[s].insert(f);
    }

    vector<pair<int,int>> degree;
    fr(i,0,n){
        degree.push_back({adj[i].size(),i});
    }
    sort(all(degree));
    unordered_map<int,int> new_index;
    fr(i,0,n){
        new_index[degree[i].second] = i;
    }
    vector<set<int>> new_adj(n);
    fr(i,0,n){
        for(int x : adj[i]){
            new_adj[new_index[i]].insert(new_index[x]);
        }
    }
    adj = new_adj;

    cerr<<"Input Read Complete."<<nL;
    cerr<<"Starting Clique Search..."<<nL;
    auto start_time=chrono::high_resolution_clock::now();
    CLIQUE();
    auto end_time=chrono::high_resolution_clock::now();
    cerr<<"Clique Search Complete."<<nL;

    auto duration = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    cout<<"Time taken: "<<duration.count()/1000<<" ms"<<nL;
    cout<<"Max Clique Size: "<<max_clique_size<<nL;
    cout<<"Clique Count: "<<clique_count<<nL;
    cout<<"Distribution: "<<nL;
    for(auto& [size, count] : distribution){
        cout<<"Size "<<size<<": "<<count<<nL;
    }
}

int main(int argc, char* argv[]){
    fast_io;

    if(argc < 3) {
        cerr << "Usage: " << argv[0] << " <input_file> <output_file>" << endl;
        return 1;
    }

    ifstream infile(argv[1]);
    if(!infile) {
        cerr << "Error: Could not open input file " << argv[1] << endl;
        return 1;
    }

    ofstream outfile(argv[2]);
    if(!outfile) {
        cerr << "Error: Could not open output file " << argv[2] << endl;
        return 1;
    }

    cin.rdbuf(infile.rdbuf());
    cout.rdbuf(outfile.rdbuf());

    ll t=1;


    while(t--){
        solve();
    }

    infile.close();
    outfile.close();

    return 0;
}

```

### Results

#### Wiki Vote Dataset

The algorithm takes 12.63 mins to run on the Wiki Vote dataset.
Since the optimized version runs faster than this version, this approach was not run for other datasets.

```bash
Time taken: 758266 ms
Max Clique Size: 17
Clique Count: 459002
Distribution:
Size 17: 23
Size 16: 208
Size 15: 740
Size 2: 8655
Size 4: 27292
Size 5: 48416
Size 7: 83266
Size 6: 68872
Size 3: 13718
Size 8: 76732
Size 9: 54456
Size 10: 35470
Size 11: 21736
Size 12: 11640
Size 13: 5449
Size 14: 2329

```

### Issues

-   Due to being unoptimized, the algo runs additional code that which is not needed for every recursive call. This contributes to the slow performance of the code (mainly for larger datasets).
-   Since this is a recursive approach, the algorithm may run into SEGMENTATION fault during deep recursions.

## Optimized Approach

### Code

```cpp
#include <bits/stdc++.h>
#include <sys/resource.h>
using namespace std;

using ll = long long;
using vi = vector<int>;
using vl = vector<ll>;

#define all(x) x.begin(), x.end()
#define fr(i, a, b) for (ll i = a; i < (b); ++i)
#define rf(i, a, b) for (ll i = b; i >=(a); i--)
#define nL "\n"
#define fast_io ios_base::sync_with_stdio(false);cin.tie(nullptr)

#define FREQ 1000

ll n,m;
vector<set<int>> adj;
vector<int> S,T;

ll max_clique_size = 0;
ll clique_count = 0;
unordered_map<ll,ll> distribution;

void increase_stack_size(rlim_t stack_size = 512 * 1024 * 1024) {
    struct rlimit rl;

    int result = getrlimit(RLIMIT_STACK, &rl);
    if (result != 0) {
        cerr << "Error getting stack limit: " << strerror(errno) << endl;
        return;
    }

    if (rl.rlim_cur < stack_size) {
        rl.rlim_cur = stack_size;
        if (rl.rlim_max < rl.rlim_cur) {
            rl.rlim_max = rl.rlim_cur;
        }

        result = setrlimit(RLIMIT_STACK, &rl);
        if (result != 0) {
            cerr << "Error setting stack limit: " << strerror(errno) << endl;
        } else {
            cerr << "Stack size increased to " << (stack_size / (1024 * 1024)) << " MB" << endl;
        }
    }
}

void UPDATE(int i, set<int>&C){
    if(i>n) return;
    if(i==n){
        clique_count++;
        distribution[C.size()]++;
        max_clique_size=max(max_clique_size,(ll)C.size());

        if(clique_count%FREQ==0){
            cerr<<"Found: "<<clique_count<<", Max Clique Size: "<<max_clique_size<<nL;
        }
        return;
    }

    vector<int> intersection;
    vector<int> difference;

    vector<int> adj_vec(adj[i].begin(), adj[i].end());

    int p1 = 0,p2 = 0;
    while (p1 < C.size() && p2 < adj_vec.size()) {
        int val1 = *next(C.begin(), p1);
        int val2 = adj_vec[p2];

        if (val1 == val2) {
            intersection.push_back(val1);
            p1++;
            p2++;
        } else if (val1 < val2) {
            difference.push_back(val1);
            p1++;
        } else {
            p2++;
        }
    }

    while (p1 < C.size()) {
        difference.push_back(*next(C.begin(), p1));
        p1++;
    }

    if(!difference.empty()){
        UPDATE(i+1,C);
    }

    for(int x:intersection){
        for(int y:adj[x]){
            if(C.find(y)==C.end() && y!=i){
                T[y]++;
            }
        }
    }

    for(int x:difference){
        for(int y:adj[x]){
            if(C.find(y)==C.end()){
                S[y]++;
            }
        }
    }

    bool FLAG=true;

    for(int y:adj[i]){
        if(C.find(y)==C.end() && y<i && T[y]==intersection.size()){
            FLAG=false;
            break;
        }
    }

    if(FLAG){
        vector<int> diff_idx(difference.begin(), difference.end());
        sort(diff_idx.begin(), diff_idx.end());

        int p=diff_idx.size();
        fr(k,1,p+1){
            int jk=diff_idx[k-1];
            bool first_handled=false;

            for(int y:adj[jk]){

                if(C.find(y)==C.end()&&y<i&&T[y]==intersection.size()){
                    if(y>=jk){
                        S[y]--;
                    }
                    else{
                        if(!first_handled){
                            if(S[y]+k-1==p&&((k-1==0&&y>=0)||(y>=diff_idx[k-2]))){
                                FLAG=false;
                                first_handled=true;
                            }

                        }
                    }
                }
            }
        }

        if(FLAG){
            int jp=(!diff_idx.empty())?diff_idx.back():0;
            if(!intersection.empty()){
                fr(y,0,i){
                    if(C.find(y)==C.end()&&y!=i&&
                    T[y]==intersection.size()&&S[y]==0){
                        if(jp<y){
                            FLAG=false;
                            break;
                        }
                    }
                }
            }
            else if(jp<i-1){
                FLAG=false;
            }
        }
    }

    for(int x:intersection){
        for(int y:adj[x]){
            if(C.find(y)==C.end()&&y!=i){
                T[y]=0;
            }
        }
    }

    for(int x:difference){
        for(int y:adj[x]){
            if(C.find(y)==C.end()){
                S[y]=0;
            }
        }
    }

    if(FLAG){
        set<int> SAVE(difference.begin(), difference.end());

        C.clear();
        C.insert(intersection.begin(), intersection.end());
        C.insert(i);

        UPDATE(i+1,C);

        C.erase(i);
        C.insert(SAVE.begin(), SAVE.end());
    }
}

void CLIQUE(){
    int start_vertex = 0;
    while(adj[start_vertex].empty()){
        start_vertex++;
    }

    S.assign(n,0);
    T.assign(n,0);

    set<int> C = {start_vertex};
    UPDATE(start_vertex+1,C);
}

void solve(){
    cin>>n>>m;
    adj.resize(n);

    cerr<<"Input Size: "<<n<<", Edges: "<<m<<nL;
    cerr<<"Reading Input..."<<nL;
    fr(i,0,m){
        ll f,s;
        cin>>f>>s;
        adj[f].insert(s);
        adj[s].insert(f);
    }

    vector<pair<int,int>> degree;
    fr(i,0,n){
        degree.push_back({adj[i].size(),i});
    }
    sort(all(degree));
    unordered_map<int,int> new_index;
    fr(i,0,n){
        new_index[degree[i].second] = i;
    }
    vector<set<int>> new_adj(n);
    fr(i,0,n){
        for(int x : adj[i]){
            new_adj[new_index[i]].insert(new_index[x]);
        }
    }
    adj = new_adj;


    cerr<<"Input Read Complete."<<nL;
    cerr<<"Starting Clique Search..."<<nL;
    auto start_time=chrono::high_resolution_clock::now();
    CLIQUE();
    auto end_time=chrono::high_resolution_clock::now();
    cerr<<"Clique Search Complete."<<nL;

    auto duration = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    cout<<"Time taken: "<<duration.count()/1000<<" ms"<<nL;
    cout<<"Max Clique Size: "<<max_clique_size<<nL;
    cout<<"Clique Count: "<<clique_count<<nL;
    cout<<"Distribution: "<<nL;
    for(auto& [size, count] : distribution){
        cout<<"Size "<<size<<": "<<count<<nL;
    }
}

int main(int argc, char* argv[]){
    fast_io;
    increase_stack_size();

    if(argc < 3) {
        cerr << "Usage: " << argv[0] << " <input_file> <output_file>" << endl;
        return 1;
    }

    ifstream infile(argv[1]);
    if(!infile) {
        cerr << "Error: Could not open input file " << argv[1] << endl;
        return 1;
    }

    ofstream outfile(argv[2]);
    if(!outfile) {
        cerr << "Error: Could not open output file " << argv[2] << endl;
        return 1;
    }

    cin.rdbuf(infile.rdbuf());
    cout.rdbuf(outfile.rdbuf());

    ll t=1;


    while(t--){
        solve();
    }

    infile.close();
    outfile.close();

    return 0;
}
```

### Optimizations

#### FLAG Check

Since only Step 4, 6 and 7 from [Chiba & Nishizeki (Page 221)](https://doi.org/10.1137/0214017) can modify the `FLAG`, we add if conditions after Step 4 and 6 to determine if going into Step 6 and 7 respectively are necessary.
In Step 4 and 7, we add break conditions after `FLAG` is set to false, as further iterations won't change the state of the `FLAG`.
This reduces computational costs needed for the extra steps, making our algorithm run faster.

```cpp
// Step 4: First FLAG check with early exit
for(int y:adj[i]){
    if(C.find(y)==C.end() && y<i && T[y]==intersection.size()){
        FLAG=false;
        break;  // Early exit once FLAG is false
    }
}

// Only proceed to Step 6 if FLAG is still true
if(FLAG){
    // Step 6 processing
    // ...

    // Only proceed to Step 7 if FLAG is still true after Step 6
    if(FLAG){
        // Step 7 processing
        // ...
    }
}
```

#### Optimized Intersection and Difference Calculation

The intersection and difference calculations are optimized using vectors instead of sets for better performance.

```cpp
vector<int> intersection;
vector<int> difference;
vector<int> adj_vec(adj[i].begin(), adj[i].end());
int p1 = 0, p2 = 0;
while (p1 < C.size() && p2 < adj_vec.size()) {
    int val1 = *next(C.begin(), p1);
    int val2 = adj_vec[p2];
    if (val1 == val2) {
        intersection.push_back(val1);
        p1++;
        p2++;
    } else if (val1 < val2) {
        difference.push_back(val1);
        p1++;
    } else {
        p2++;
    }
}
while (p1 < C.size()) {
    difference.push_back(*next(C.begin(), p1));
    p1++;
}
```

#### Increased Stack Size

In the Optimized Approach, the stack size is increased to handle deeper recursion.

```cpp
void increase_stack_size(rlim_t stack_size = 512 * 1024 * 1024) {
    struct rlimit rl;
    int result = getrlimit(RLIMIT_STACK, &rl);
    if (result != 0) {
        cerr << "Error getting stack limit: " << strerror(errno) << endl;
        return;
    }
    if (rl.rlim_cur < stack_size) {
        rl.rlim_cur = stack_size;
        if (rl.rlim_max < rl.rlim_cur) {
            rl.rlim_max = rl.rlim_cur;
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

### Results

#### Wiki Vote Dataset

The algorithm takes 6.98 mins to run on the Wiki Vote dataset.
This shows a 44% decrease in run time.

```bash
Time taken: 419446 ms
Max Clique Size: 17
Clique Count: 459002
Distribution:
Size 17: 23
Size 16: 208
Size 15: 740
Size 2: 8655
Size 4: 27292
Size 5: 48416
Size 7: 83266
Size 6: 68872
Size 3: 13718
Size 8: 76732
Size 9: 54456
Size 10: 35470
Size 11: 21736
Size 12: 11640
Size 13: 5449
Size 14: 2329
```

#### Email Enron Dataset

The algorithm takes 9.3 mins to run on the Email Enron dataset.
This shows a 44% decrease in run time.

```bash
Time taken: 558357 ms
Max Clique Size: 20
Clique Count: 226859
Distribution:
Size 20: 6
Size 19: 10
Size 18: 41
Size 17: 286
Size 16: 1178
Size 14: 7417
Size 2: 14070
Size 15: 3157
Size 4: 13319
Size 3: 7077
Size 5: 18143
Size 7: 25896
Size 6: 22715
Size 8: 24766
Size 9: 22884
Size 10: 21393
Size 13: 11487
Size 11: 17833
Size 12: 15181

```

#### Skitter Dataset

The algorithm runs into `SEGMENTATION` fault after a few minutes, even when the stack size is increased to 1GB.

### Issues

-   Since this is also a recursive approach, the algorithm may run out of stack space and give `SEGMENTATION` fault for deep recursion. Although we can avoid this by increasing stack space, it's better to use an iterative approach.

## Iterative Approach

We simulate the recursion by using a stack and a custom struct which holds the same parameters as our recursive function.

### Code

```cpp
#include <bits/stdc++.h>
#include <sys/resource.h>
using namespace std;

using ll = long long;
using vi = vector<int>;
using vl = vector<ll>;

#define all(x) x.begin(), x.end()
#define fr(i, a, b) for (ll i = a; i < (b); ++i)
#define rf(i, a, b) for (ll i = b; i >=(a); i--)
#define nL "\n"
#define fast_io ios_base::sync_with_stdio(false);cin.tie(nullptr)

#define FREQ 1000

ll n,m;
vector<set<int>> adj;
vector<int> S,T;

ll max_clique_size = 0;
ll clique_count = 0;
unordered_map<ll,ll> distribution;

void UPDATE(int i, set<int>& C) {
    struct StackFrame {
        int i;
        set<int> C;
        vector<int> intersection;
        vector<int> difference;
        vector<int> diff_idx;
        bool FLAG;
        int state;
        set<int> SAVE;
    };

    stack<StackFrame> callStack;
    callStack.push({i, C, {}, {}, {}, true, 0, {}});

    while (!callStack.empty()) {
        auto& frame = callStack.top();

        if (frame.i > n) {
            callStack.pop();
            continue;
        }

        if (frame.i == n) {
            clique_count++;
            distribution[frame.C.size()]++;
            max_clique_size = max(max_clique_size, (ll)frame.C.size());

            if (clique_count % FREQ == 0) {
                cerr << "Found: " << clique_count << ", Max Clique Size: " << max_clique_size << nL;
            }

            callStack.pop();
            continue;
        }

        if (frame.state == 0) {
            vector<int> adj_vec(adj[frame.i].begin(), adj[frame.i].end());

            int p1 = 0, p2 = 0;
            while (p1 < frame.C.size() && p2 < adj_vec.size()) {
                int val1 = *next(frame.C.begin(), p1);
                int val2 = adj_vec[p2];

                if (val1 == val2) {
                    frame.intersection.push_back(val1);
                    p1++;
                    p2++;
                } else if (val1 < val2) {
                    frame.difference.push_back(val1);
                    p1++;
                } else {
                    p2++;
                }
            }

            while (p1 < frame.C.size()) {
                frame.difference.push_back(*next(frame.C.begin(), p1));
                p1++;
            }

            if (!frame.difference.empty()) {

                frame.state = 1;


                callStack.push({frame.i + 1, frame.C, {}, {}, {}, true, 0, {}});
                continue;
            } else {
                frame.state = 1;
            }
        }


        if (frame.state == 1) {

            for (int x : frame.intersection) {
                for (int y : adj[x]) {
                    if (frame.C.find(y) == frame.C.end() && y != frame.i) {
                        T[y]++;
                    }
                }
            }

            for (int x : frame.difference) {
                for (int y : adj[x]) {
                    if (frame.C.find(y) == frame.C.end()) {
                        S[y]++;
                    }
                }
            }

            frame.FLAG = true;


            for (int y : adj[frame.i]) {
                if (frame.C.find(y) == frame.C.end() && y < frame.i && T[y] == frame.intersection.size()) {
                    frame.FLAG = false;
                    break;
                }
            }


            if (frame.FLAG) {
                frame.diff_idx = vector<int>(frame.difference.begin(), frame.difference.end());
                sort(frame.diff_idx.begin(), frame.diff_idx.end());

                int p = frame.diff_idx.size();
                for (int k = 1; k <= p; k++) {
                    int jk = frame.diff_idx[k-1];
                    bool first_handled = false;

                    for (int y : adj[jk]) {
                        if (frame.C.find(y) == frame.C.end() && y < frame.i && T[y] == frame.intersection.size()) {
                            if (y >= jk) {
                                S[y]--;
                            } else {
                                if (!first_handled) {
                                    if (S[y] + k - 1 == p && ((k - 1 == 0 && y >= 0) || (y >= frame.diff_idx[k-2]))) {
                                        frame.FLAG = false;
                                        first_handled = true;
                                    }
                                }
                            }
                        }
                    }
                }


                if (frame.FLAG) {
                    int jp = (!frame.diff_idx.empty()) ? frame.diff_idx.back() : 0;
                    if (!frame.intersection.empty()) {
                        for (int y = 0; y < frame.i; y++) {
                            if (frame.C.find(y) == frame.C.end() && y != frame.i &&
                            T[y] == frame.intersection.size() && S[y] == 0) {
                                if (jp < y) {
                                    frame.FLAG = false;
                                    break;
                                }
                            }
                        }
                    } else if (jp < frame.i - 1) {
                        frame.FLAG = false;
                    }
                }
            }


            for (int x : frame.intersection) {
                for (int y : adj[x]) {
                    if (frame.C.find(y) == frame.C.end() && y != frame.i) {
                        T[y] = 0;
                    }
                }
            }

            for (int x : frame.difference) {
                for (int y : adj[x]) {
                    if (frame.C.find(y) == frame.C.end()) {
                        S[y] = 0;
                    }
                }
            }


            if (frame.FLAG) {
                frame.SAVE = set<int>(frame.difference.begin(), frame.difference.end());

                set<int> newC;
                newC.insert(frame.intersection.begin(), frame.intersection.end());
                newC.insert(frame.i);

                frame.state = 2;


                callStack.push({frame.i + 1, newC, {}, {}, {}, true, 0, {}});
                continue;
            } else {

                callStack.pop();
                continue;
            }
        }


        if (frame.state == 2) {

            if (!callStack.empty()) {
                auto& prevFrame = callStack.top();
                prevFrame.C.erase(prevFrame.i);
                prevFrame.C.insert(prevFrame.SAVE.begin(), prevFrame.SAVE.end());
            }

            callStack.pop();
        }
    }


    if (!callStack.empty()) {
        C = callStack.top().C;
    }
}

void CLIQUE(){
    int start_vertex = 0;
    while(adj[start_vertex].empty()){
        start_vertex++;
    }

    S.assign(n,0);
    T.assign(n,0);

    set<int> C = {start_vertex};
    UPDATE(start_vertex+1,C);
}

void solve(){
    cin>>n>>m;
    adj.resize(n);

    cerr<<"Input Size: "<<n<<", Edges: "<<m<<nL;
    cerr<<"Reading Input..."<<nL;
    fr(i,0,m){
        ll f,s;
        cin>>f>>s;
        adj[f].insert(s);
        adj[s].insert(f);
    }

    vector<pair<int,int>> degree;
    fr(i,0,n){
        degree.push_back({adj[i].size(),i});
    }
    sort(all(degree));
    unordered_map<int,int> new_index;
    fr(i,0,n){
        new_index[degree[i].second] = i;
    }
    vector<set<int>> new_adj(n);
    fr(i,0,n){
        for(int x : adj[i]){
            new_adj[new_index[i]].insert(new_index[x]);
        }
    }
    adj = new_adj;


    cerr<<"Input Read Complete."<<nL;
    cerr<<"Starting Clique Search..."<<nL;
    auto start_time=chrono::high_resolution_clock::now();
    CLIQUE();
    auto end_time=chrono::high_resolution_clock::now();
    cerr<<"Clique Search Complete."<<nL;

    auto duration = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    cout<<"Time taken: "<<duration.count()/1000<<" ms"<<nL;
    cout<<"Max Clique Size: "<<max_clique_size<<nL;
    cout<<"Clique Count: "<<clique_count<<nL;
    cout<<"Distribution: "<<nL;
    for(auto& [size, count] : distribution){
        cout<<"Size "<<size<<": "<<count<<nL;
    }
}

int main(int argc, char* argv[]){
    fast_io;

    if(argc < 3) {
        cerr << "Usage: " << argv[0] << " <input_file> <output_file>" << endl;
        return 1;
    }

    ifstream infile(argv[1]);
    if(!infile) {
        cerr << "Error: Could not open input file " << argv[1] << endl;
        return 1;
    }

    ofstream outfile(argv[2]);
    if(!outfile) {
        cerr << "Error: Could not open output file " << argv[2] << endl;
        return 1;
    }

    cin.rdbuf(infile.rdbuf());
    cout.rdbuf(outfile.rdbuf());

    ll t=1;


    while(t--){
        solve();
    }

    infile.close();
    outfile.close();

    return 0;
}
```

Here, the state variable in StackFrame holds info regarding the sequence of the call. If the call is a first call, we process it similar to a normal recursive call. But if it's not, we only run the part where the recursive function would run after it comes back from a recursive call.

### Results

Results are similar to the [Optimized Approach](#optimized-approach) for Wiki and Email Dataset.

For Skitter, the algorithm ran for 6hrs and found 52k cliques.

This gives us an `estimated run time` of **4269hrs** for **37M cliques**. Thus, the algorithm was terminated.
