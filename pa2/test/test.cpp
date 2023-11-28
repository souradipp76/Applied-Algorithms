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
#include <climits>

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

class TreeNode {
    public:
        string nodeId = "";
        vector<pair<int, int>> rectangles;
        vector<pair<int, int>> indices;
        int X = 0;
        int Y = 0;
        char cutLine = '-';
        TreeNode* left = NULL;
        TreeNode* right = NULL;
    
        TreeNode(string id, vector<pair<int, int>> rect, vector<pair<int, int>> ind) : nodeId(id), rectangles(rect), indices(ind) {}
        TreeNode(string id, char cut, TreeNode* l, TreeNode* r) : nodeId(id), cutLine(cut), left(l), right(r) {}
};


void deleteTree(TreeNode* root) {
    if(root == NULL) {
        return;
    }

    deleteTree(root->left);
    deleteTree(root->right);

    delete root;
}

void postOrder(TreeNode* root, vector<string> &lines) {
    if(root == NULL) {
        return;
    }

    postOrder(root->left, lines);
    postOrder(root->right, lines);

    if(root->left && root->right) {
        string line(1,root->cutLine);
        lines.push_back(line);
        // sprintf(line, "%s\n", root->nodeId.c_str());
        // cout<<root->nodeId<<":("<<root->rectangles[index].first<<","<<root->rectangles[index].second<<"):"<<root->X<<":"<<root->Y<<endl;
        // for(auto a: root->rectangles){
        //     cout<<"("<<a.first<<","<<a.second<<")";
        // }
        // cout<<endl;
    } else {
        int len = root->rectangles.size();
        string line = root->nodeId+"(";
        for(int i=0;i<len;i++){
            int w = root->rectangles[i].first;
            int h = root->rectangles[i].second;
            line+="("+to_string(w)+","+to_string(h)+")";
        }
        line+=")";
        lines.push_back(line);
    }
}

TreeNode* binaryTreeGenerator(int n, int &leafId){
    if (n == 0) {
        return NULL;
    }

    if (n == 1) {
        int num_rect = 1+rand()%100;
        vector<pair<int, int>> rect;
        vector<pair<int, int>> ind;
        for(int i =0;i<num_rect;i++){
            int w = rand()%1000;
            int h = rand()%1000;
            rect.push_back({w,h});
            ind.push_back({-1,-1});
        }
        TreeNode* node = new TreeNode(to_string(leafId), rect, ind);
        leafId++;
        return node;
    }

    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    char cut = 'V';
    if(r > 0.5) {
        cut = 'H';
    }
    
    int leftN = rand()%(n-1);
    if(leftN %2 == 0) {
        leftN++;
    }

    // cout<<"Left"<<leftN<<":"<<"Right"<<n-leftN-1<<endl;
    TreeNode* root = new TreeNode("", cut, NULL, NULL);
    
    // Recursively build each subtree
    root ->left = binaryTreeGenerator(leftN, leafId);
    root ->right = binaryTreeGenerator(n - leftN - 1, leafId);

    return root;
}


int main(int argc, char *argv[])
{
    int leafId = 1;
    int n = stoi(argv[1]);
    if(!(n&1)) {
        cout<<"Provide an odd number of nodes."<<endl;
        return EXIT_FAILURE;
    }

    srand(42);
    TreeNode* root = binaryTreeGenerator(n, leafId);
    cout<< leafId<<endl;
    vector<string> lines;
    postOrder(root, lines);
    cout<<lines.size()<<endl;
    writeFile(argv[2], lines);

    // Release resource
    deleteTree(root);

    return EXIT_SUCCESS;
}