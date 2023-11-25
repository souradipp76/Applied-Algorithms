#include <iostream>
#include <vector>
#include <string>
#include <stack>
#include <queue>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <unordered_map>
#include <cmath>
#include <cfloat>

#define INF 1000000
using namespace std;

void readBinaryFile(const string filename) {
    FILE* f = fopen(filename.c_str(), "rb");
    float m;
    fread(&m, sizeof(float), 1, f);
    printf("%f\n", m);
    fclose(f);
}

void writeBinaryFile(const string filename, float m) {
    FILE* f = fopen(filename.c_str(), "wb");
    fwrite(&m, sizeof(double), 1, f);
    fclose(f);
}

vector<string> readFile(const string filename){
    vector<string> lines;
    ifstream file(filename);
    string s;
    if ( file.is_open() ) {
        while ( file ) {
            getline(file, s);
            lines.push_back(s);
        }
    }
    file.close();
    lines.pop_back();
    return lines;
}

void writeFile(const string filename, vector<string> lines)
{
    ofstream file(filename);
    if (file.is_open()) {
        for(auto line: lines){
            file << line;
        }
    }
    file.close();
}

vector<string> splitLine(const string s){
    vector<string> words;
    stringstream ss(s);  
    string word;
    while (ss >> word) {
        words.push_back(word);
    }
    ss.clear();
    return words;
}

void findSmallestCycleMeanInPath(vector<vector<double>> &D, vector<vector<int>> &P, 
        int k, int vs, double &minCycleMean, vector<int> &level, vector<int> &minCycle) {
    vector<int> level_stack;
    int v = vs;
    for(int j = k; j>=0 ;j--){
        if(level[v] > -1){
            double l = (D[level[v]][v] - D[j][v])/(level[v] - j);
            int cycleLength = level[v] - j;
            if(l<minCycleMean) {
                minCycleMean = l;
                minCycle.clear();
                int c = level_stack.size();
                for(int i=c-1;i>=c-cycleLength;i--){
                    minCycle.push_back(level_stack[i]);
                    //cout<<cycle[i]<<"->";
                }
            }
        }

        level[v] = j;
        level_stack.push_back(v);
        v = P[j][v];
    }

    int l = level_stack.size();
    for(int i = 0;i<l;i++){
        level[level_stack[i]] = -1;
    }
}

double findSmallestCycleMean(vector<vector<double>> &D, vector<vector<int>> &P, int k, int n, vector<int> &level, vector<int> &minCycle) {
    double minCycleMean = INF;
    for(int v=0; v<n; v++){
        if(D[k][v] != INF){
            findSmallestCycleMeanInPath(D, P, k, v, minCycleMean, level, minCycle);
        }
    }
    return minCycleMean;
}

bool isPower2(int k){
    double x = log2(k);
    return ceil(x) == floor(x);
}

vector<pair<double, pair<int, int>>> calculatePi(vector<vector<double>> &D, vector<vector<int>> &P, 
        int k, int n, double minCycleMean){
    vector<pair<double, pair<int, int>>> pi(n);
    for(int v = 0; v <n; v++){
        pair<int, int> p = {-1, 0};
        pi[v] = {INF, p};
        for(int j = 0; j <= k; j++){
            if(pi[v].first > D[j][v] - j*minCycleMean){
                pi[v].first = D[j][v] - j*minCycleMean;
                pi[v].second.first = P[j][v];
                pi[v].second.second = j;
            }
        }
    }
    return pi;
}

bool checkPi(vector<vector<pair<int, double>>> &adj, vector<pair<double, pair<int, int>>> &pi, double &minCycleMean){
    int n = adj.size();
    for(int i=0;i<n;i++){
        int len = adj[i].size();
        for(int j = 0;j<len;j++){
            int v = adj[i][j].first;
            long double x = pi[i].first + adj[i][j].second - minCycleMean + 1e-6;
            if(pi[v].first > x){
                return false;
            }
        }
    }
    return true;
}

bool checkEarlyTermination(vector<vector<pair<int, double>>> &adj, vector<vector<double>> &D, vector<vector<int>> &P,  
        int k, vector<int> &level, double &minMean, vector<int> &minCycle){
    int n = adj.size();
    if(isPower2(k) || k == n-1) {
        // cout<<"here"<<endl;
        double minCycleMean = findSmallestCycleMean(D, P, k, n, level, minCycle);
        // cout<<minCycleMean<<endl;
        minMean = minCycleMean;

        vector<pair<double, pair<int, int>>> pi = calculatePi(D, P, k, n, minMean);
        // cout<<"Pi"<<endl;

        bool done = checkPi(adj, pi, minMean);
        // cout<<"Done:"<<done<<endl;
        if(done) {
            return true;
        }
    }
    // cout<<"here2"<<endl;
    return false;
}   

void init(vector<vector<double>> &D, vector<vector<int>> &P, int n){
    vector<double> d(n, INF);
    vector<int> p(n, -1);
    D.push_back(d);
    P.push_back(p);
}

double mcm(vector<vector<pair<int, double>>> &adj, vector<int> &minCycle) {
    int n = adj.size();
    double minMean = INF;
    vector<vector<double>> D;
    vector<vector<int>> P;
    vector<int> level(n, -1);

    init(D, P, n);
    D[0][0] = 0;

    for(int k = 1;k<n+1;k++) {
        // cout<<"######## k = "<<k<<" ########"<<endl;
        init(D, P, n);
        for(long u=0;u<n;u++) {
            long len = adj[u].size();
            for(long j=0;j<len;j++){
                long v = adj[u][j].first;
                double w = adj[u][j].second;
                long double newPath = D[k-1][u] == INF ? INF : D[k-1][u] + w;
                if(D[k][v] > newPath) {
                    D[k][v] = D[k-1][u]+w;
                    P[k][v] = u;
                }
            }
        }

        // for(int i=0;i<=k;i++) {
        //     for(int j=0;j<n;j++){
        //         cout<<D[i][j]<<":";
        //     }
        //     cout<<endl;
        // }

        // for(int i=0;i<=k;i++) {
        //     for(int j=0;j<n;j++){
        //         cout<<P[i][j]<<":";
        //     }
        //     cout<<endl;
        // }
        
        bool terminate = checkEarlyTermination(adj, D, P, k, level, minMean, minCycle);
        // cout<<"T:"<<terminate<<endl;

        if(terminate){
            break;
        }
    }

    return minMean;
}


vector<vector<pair<int, double>>> constructGraph(vector<string> &gr, long &n, long &e){
    vector<vector<pair<int, double>>> adj;
    for(auto line: gr){
        vector<string> parts = splitLine(line);
        if(parts[0] == "V") {
            n = stoi(parts[1]);
            for(long i=0;i<n;i++){
                vector<pair<int, double>> p;
                adj.push_back(p);
            }
        } else if (parts[0] == "E"){
            e = stol(parts[1]);
        } else {
            int y = stoi(parts[0]);
            int x = stoi(parts[1]);
            double w = stod(parts[2]);
            adj[x-1].push_back({y-1, w});
        }
    }

    return adj;
}


int main(int argc, char *argv[])
{
    vector<string> graph = readFile(argv[1]);
    long n = 0, e = 0;
    vector<vector<pair<int, double>>> adj = constructGraph(graph, n, e);
    cout<<n<<":"<<e<<endl;

    vector<int> minMeanCycle;
    double minCycleMean = mcm(adj, minMeanCycle);
    cout << minCycleMean <<endl;
    writeBinaryFile(argv[2], minCycleMean);

    int c = minMeanCycle.size();
    string res = "";
    for(int i=c-1;i>=0;i--){
        string delimeter = (i == 0) ? "\n" : " ";
        res += to_string(minMeanCycle[i]+1) + delimeter;
    }
    cout<<res;

    vector<string> lines;
    lines.push_back(res);
    writeFile(argv[3], lines);
   
    return EXIT_SUCCESS;
}