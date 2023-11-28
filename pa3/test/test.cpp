#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cmath>

using namespace std;

void writeFile(const string filename, vector<string> lines)
{
    ofstream file(filename);
    if (file.is_open()) {
        for(auto line: lines){
            file << line<<endl;
        }
    }
    file.close();
}

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

int main(int argc, char *argv[])
{
    srand(42);
    vector<string> test_gr;
    int n = stoi(argv[1]);
    int e = stoi(argv[2]);
    test_gr.push_back("V "+to_string(n));
    test_gr.push_back("E "+to_string(e));
    for(int i =0;i<e;i++){
        int x = 1 + (rand() % n);
        int y = 1 + (rand() % n);
        if(x == y) {
            continue;
        }
        double c = fRand(-100.0, 100.0);
        test_gr.push_back(to_string(y)+" "+to_string(x)+ " "+ to_string(c));
    }

    writeFile(argv[3], test_gr);
    return EXIT_SUCCESS;
}