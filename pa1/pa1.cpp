#include<bits/stdc++.h>

vector<string> readFile(const string filename)
{
    vector<string> lines
    ifstream file(filename);
    string s;
    while (getline(file, s))
        lines.push_back(s);

    return lines;
}

vector<string> splitLine(const string s){
    vector<string> words;
    stringstream ss(s);  
    string word;
    while (ss >> word) {
        words.push_back(word)
    }
    return words;
}

class TreeeNode {
    private:
        int id;
        int capacitance;
        TreeeNode* left;
        TreeeNode* right;
    public:
    
    TreeNode(int id, int cap) : id(x), capacitance(cap), left(nullptr), right(nullptr) {}
    TreeNode(int id, int cap, TreeNode *left, TreeNode *right) : id(x), capacitance(cap), left(left), right(right) {}
}



int main(string[] argv) {

    vector<string> wire = readFile(argv[2]);
    vector<string> rc = splitLine(wire[0]);
    float r = rc[0], c = rc[1];

    cout<<r<< ":" <, c <<endl;
}