#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cmath>

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
        int k, int v, double &minCycleMean, vector<int> &level, vector<int> &minCycle) {
    vector<int> path;
    for(int j = k; j>=0 ;j--){
        if(level[v] > -1){
            double lambda = (D[level[v]][v] - D[j][v])/(level[v] - j);
            int c = level[v] - j;
            if(lambda < minCycleMean) {
                minCycleMean = lambda;
                minCycle.clear();
                int l = path.size();
                for(int i = l-1; i >= l-c; i--){
                    minCycle.push_back(path[i]);
                }
            }
        }

        level[v] = j;
        path.push_back(v);
        v = P[j][v];
    }

    int l = path.size();
    for(int i = 0; i<l; i++){
        level[path[i]] = -1;
    }
}

double findSmallestCycleMean(vector<vector<double>> &D, vector<vector<int>> &P, int k, int n, vector<int> &level, vector<int> &minCycle) {
    double minCycleMean = INFINITY;
    for(int v=0; v<n; v++){
        if(!isinf(D[k][v])){
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
        pi[v] = {INFINITY, p};
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

double getTolerance(double minCycleMean){
    if(!isinf(minCycleMean)) {
        if(abs(minCycleMean) < 1.0e-15 ) {
            return 1.0e-15;
        } else {
            return abs(minCycleMean)*1.0e-6;
        }   
    }

    return 1.0e-6;
}

bool checkPi(vector<vector<pair<int, double>>> &adj, vector<pair<double, pair<int, int>>> &pi, double &minCycleMean){
    int n = adj.size();
    for(int i=0;i<n;i++){
        int len = adj[i].size();
        for(int j = 0;j<len;j++){
            int v = adj[i][j].first;
            long double x = pi[i].first + adj[i][j].second - minCycleMean + getTolerance(minCycleMean);
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
        double minCycleMean = findSmallestCycleMean(D, P, k, n, level, minCycle);
        // cout<<minCycleMean<<endl;
        minMean = minCycleMean;

        vector<pair<double, pair<int, int>>> pi = calculatePi(D, P, k, n, minMean);

        bool done = checkPi(adj, pi, minMean);
        // cout<<"Done:"<<done<<endl;
        if(done) {
            return true;
        }
    }
    return false;
}   

void init(vector<vector<double>> &D, vector<vector<int>> &P, int n){
    vector<double> d(n, INFINITY);
    vector<int> p(n, -1);
    D.push_back(d);
    P.push_back(p);
}

double mcm(vector<vector<pair<int, double>>> &adj, vector<int> &minCycle) {
    int n = adj.size();
    double minMean = INFINITY;
    vector<vector<double>> D;
    vector<vector<int>> P;
    vector<int> level(n, -1);

    init(D, P, n);
    D[0][0] = 0;

    for(int k = 1;k<n+1;k++) {
        init(D, P, n);
        for(long u=0;u<n;u++) {
            long len = adj[u].size();
            for(long j=0;j<len;j++){
                long v = adj[u][j].first;
                double w = adj[u][j].second;
                long double newDist = isinf(D[k-1][u]) ? INFINITY : D[k-1][u] + w;
                if(D[k][v] > newDist) {
                    D[k][v] = D[k-1][u]+w;
                    P[k][v] = u;
                }
            }
        }

        // cout<<"######## k = "<<k<<" ########"<<endl;
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
            cout<<"Terminated at K = "<<k<<endl;
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
            n = stol(parts[1]);
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

// double fRand(double fMin, double fMax)
// {
//     double f = (double)rand() / RAND_MAX;
//     return fMin + f * (fMax - fMin);
// }

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

    int cycleSize = minMeanCycle.size();
    string cycleStr = "";
    for(int i=cycleSize-1;i>=0;i--){
        int node = minMeanCycle[i]+1;
        string delimeter = (i == 0) ? "\n" : " ";
        cycleStr += to_string(node) + delimeter;
    }
    cout<<cycleStr;
    vector<string> lines;
    lines.push_back(cycleStr);
    writeFile(argv[3], lines);

    ////////// Testing with random graph ///////////////
    // srand(10);
    // vector<string> test_gr;
    // test_gr.push_back("V 50000");
    // test_gr.push_back("E 200000");
    // for(int i =0;i<200000;i++){
    //     int x = 1 + (rand() % 50000);
    //     int y = 1 + (rand() % 50000);
    //     double c = fRand(-100.0, 100.0);
    //     test_gr.push_back(to_string(y)+" "+to_string(x)+ " "+ to_string(c));
    // }

    // vector<vector<pair<int, double>>> adj2 = constructGraph(test_gr, n, e);
    // cout<<n<<":"<<e<<endl;

    // vector<int> minMeanCycle2;
    // double minCycleMean2 = mcm(adj2, minMeanCycle2);
    // cout << minCycleMean2 <<endl;

    // int s = minMeanCycle2.size();
    // string res = "";
    // for(int i=s-1;i>=0;i--){
    //     int node = minMeanCycle2[i]+1;
    //     string delimeter = (i == 0) ? "\n" : " ";
    //     res += to_string(node) + delimeter;
    // }
    // cout<<res;
   
    return EXIT_SUCCESS;
}