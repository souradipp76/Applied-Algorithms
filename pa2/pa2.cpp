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

bool compareWidth(pair<pair<int,int>,int> A, pair<pair<int,int>,int> B){
    pair<int,int> a = A.first;
    pair<int,int> b = B.first;
    if(a.first < b.first) {
        return true;
    } else if (a.first == b.first) {
        return a.second > b.second;
    }
    return false;
}

bool compareHeight(pair<pair<int,int>,int> A, pair<pair<int,int>,int> B){
    pair<int,int> a = A.first;
    pair<int,int> b = B.first;
    if(a.second < b.second) {
        return true;
    } else if (a.second == b.second) {
        return a.first > b.first;
    }
    return false;
}

vector<pair<int,int>> getRectangles(string rectangleStr) {
    int s = rectangleStr.size();
    string oneRect = "";
    int done = 0;
    int w,h;
    vector<pair<int,int>> v;
    for(int i=0;i<s;i++){
        if(rectangleStr[i] == ')'){
            oneRect+=rectangleStr[i];

            done = sscanf(oneRect.c_str(),"(%d,%d)", &w, &h);
            if(done > 0) {
                pair<int,int> p = {w, h};
                v.push_back(p);
                oneRect = "";
            }

        } else {
            oneRect+=rectangleStr[i];
        }
    }

    done = sscanf(oneRect.c_str(),"(%d,%d)", &w, &h);
    if(done > 0) {
        pair<int,int> p = {w, h};
        v.push_back(p);
    }

    return v;
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
    
        TreeNode(string id, vector<pair<int, int>> rect) : nodeId(id), rectangles(rect) {}
        TreeNode(string id, char cut, TreeNode* l, TreeNode* r) : nodeId(id), cutLine(cut), left(l), right(r) {}
};

TreeNode* constructTree(vector<string> nodes) {
    TreeNode *head = NULL;  

    stack<TreeNode*> st;
    for(auto &node: nodes) {
        bool isLeaf = true;
        char c;

        int done = node.size() == 1;
        if (done > 0) {
            isLeaf = false;
            c = node[0];
        }

        if (!isLeaf) {
            TreeNode* r = st.top();
            st.pop();
            TreeNode* l = st.top();
            st.pop();

            string newId = l->nodeId +"-"+ c + "-"+r->nodeId;
            TreeNode *treeNode = new TreeNode("", c, l, r);
            st.push(treeNode);
        } else {
            int id;
            char rect[200];
            done = sscanf(node.c_str(),"%d%s", &id, rect);
            string rectStr(rect);
            rectStr.erase(rectStr.begin());
            rectStr.pop_back();
            vector<pair<int, int>> rectangles = getRectangles(rectStr);
            TreeNode *treeNode = new TreeNode(to_string(id), rectangles);
            st.push(treeNode);
        }
    }

    if(st.empty()) {
        return NULL;
    }

    head = st.top();
    return head;
}

void deleteTree(TreeNode* root) {
    if(root == NULL) {
        return;
    }

    deleteTree(root->left);
    deleteTree(root->right);

    delete root;
}

void postOrder(TreeNode* root, int index, bool single, vector<string> &lines) {
    if(root == NULL) {
        return;
    }

    int lIdx = single ? 0 : root->indices[index].first;
    int rIdx = single ? 0 : root->indices[index].second;
    postOrder(root->left, lIdx, single, lines);
    postOrder(root->right, rIdx, single, lines);

    char line[200];
    if(root->left && root->right) {
        sprintf(line, "%s\n", root->nodeId.c_str());
        // cout<<root->nodeId<<":("<<root->rectangles[index].first<<","<<root->rectangles[index].second<<"):"<<root->X<<":"<<root->Y<<endl;
        // for(auto a: root->rectangles){
        //     cout<<"("<<a.first<<","<<a.second<<")";
        // }
        // cout<<endl;
    } else {
        int w = root->rectangles[index].first;
        int h = root->rectangles[index].second;
        sprintf(line, "%s((%d,%d)(%d,%d))\n", root->nodeId.c_str(), w, h, root->X, root->Y);
        // cout<<root->nodeId<<":("<<root->rectangles[index].first<<","<<root->rectangles[index].second<<"):"<<root->X<<":"<<root->Y<<endl;
        // for(auto a: root->rectangles){
        //     cout<<"("<<a.first<<","<<a.second<<")";
        // }
        // cout<<endl;
        
        lines.push_back(line);
    }
}

void preOrderXY(TreeNode* root, int index, bool single) {
    if(!root->left && !root->right) {
        return;
    }

    int lIdx = single ? 0 : root->indices[index].first;
    int rIdx = single ? 0 : root->indices[index].second;

    if(root->cutLine == 'V') {
        root->right->X = root->X + root->left->rectangles[lIdx].first;
        root->right->Y = root->Y;

        root->left->X = root->X;
        root->left->Y = root->Y;
    } else {
        root->left->Y = root->Y + root->right->rectangles[rIdx].second;
        root->left->X = root->X;

        root->right->X = root->X;
        root->right->Y = root->Y;
    }

    preOrderXY(root->left, lIdx, single);
    preOrderXY(root->right, rIdx, single);
}

pair<int, int> singlePacking(TreeNode* root){
    if(!root ->left && !root->right) {
        return root->rectangles[0];
    }

    pair<int, int> pl = singlePacking(root->left);
    pair<int, int> pr = singlePacking(root->right);

    int W, H;
    if(root-> cutLine == 'V') {
        W = pl.first+pr.first;
        H = max(pl.second, pr.second);
        root->rectangles.push_back({W, H});
    } else {
        W = max(pl.first, pr.first);
        H = pl.second+pr.second;
        root->rectangles.push_back({W, H});
    }

    return {W, H};
}

void doOptimalPacking(TreeNode* root){
    if(!root ->left && !root->right) {
        root->indices.push_back({-1,-1});
        return;
    }

    if(root->indices.size() > 0) {
        return;
    }

    doOptimalPacking(root->left);
    doOptimalPacking(root->right);

    if(root-> cutLine == 'V') {

        vector<pair<int, int>> vpl = root->left->rectangles;
        vector<pair<int, int>> vpr = root->right->rectangles;
        
        int llen = vpl.size();
        int rlen = vpr.size();

        vector<pair<pair<int,int>, int>> vl(llen);
        vector<pair<pair<int,int>, int>> vr(rlen);
        for(int i=0;i<llen;i++) {
            vl[i] = {vpl[i], i};
        }

        
        for(int i=0;i<rlen;i++) {
            vr[i] = {vpr[i], i};
        }

        sort(vl.begin(), vl.end(), compareWidth);
        sort(vr.begin(), vr.end(), compareWidth);

        vector<pair<int, int>> leftInds(llen);
        for(int i=0;i<llen;i++) {
            vpl[i] = vl[i].first;
            leftInds[i] = root->left->indices[vl[i].second];
        }

        root->left->rectangles = vpl;
        root->left->indices = leftInds;

        vector<pair<int, int>> rightInds(rlen);
        for(int i=0;i<rlen;i++) {
            vpr[i] = vr[i].first;
            rightInds[i] = root->right->indices[vr[i].second];
        }

        root->right->rectangles = vpr;
        root->right->indices = rightInds;

        int i=0,j=0;
        root->rectangles.clear();
        while(i<llen && j< rlen) {
            pair<int,int> p = {vpl[i].first+ vpr[j].first, max(vpl[i].second, vpr[j].second)};
            root->rectangles.push_back(p);
            root->indices.push_back({i,j});
            if(vpl[i].second > vpr[j].second) {
                i++;
            } else if(vpl[i].second < vpr[j].second) {
                j++;
            } else {
                i++;
                j++;
            }
        }
    } else {

        vector<pair<int, int>> vpl = root->left->rectangles;
        vector<pair<int, int>> vpr = root->right->rectangles;
        int llen = vpl.size();
        int rlen = vpr.size();

        vector<pair<pair<int,int>, int>> vl(llen);
        vector<pair<pair<int,int>, int>> vr(rlen);
        for(int i=0;i<llen;i++) {
            vl[i] = {vpl[i], i};
        }

        for(int i=0;i<rlen;i++) {
            vr[i] = {vpr[i], i};
        }

        sort(vl.begin(), vl.end(), compareHeight);
        sort(vr.begin(), vr.end(), compareHeight);

        vector<pair<int, int>> leftInds(llen);
        for(int i=0;i<llen;i++) {
            vpl[i] = vl[i].first;
            leftInds[i] = root->left->indices[vl[i].second];
        }

        root->left->rectangles = vpl;
        root->left->indices = leftInds;

        vector<pair<int, int>> rightInds(rlen);
        for(int i=0;i<rlen;i++) {
            vpr[i] = vr[i].first;
            rightInds[i] = root->right->indices[vr[i].second];
        }

        root->right->rectangles = vpr;
        root->right->indices = rightInds;

        int i=0,j=0;
        root->rectangles.clear();
        while(i < llen && j < rlen) {
            pair<int,int> p = {max(vpl[i].first, vpr[j].first), vpl[i].second+ vpr[j].second};
            root->rectangles.push_back(p);
            root->indices.push_back({i,j});
            if(vpl[i].first > vpr[j].first) {
                i++;
            } else if(vpl[i].first < vpr[j].first) {
                j++;
            } else {
                i++;
                j++;
            }
        }
    }
}

pair<int,int> getOptimalPacking(TreeNode* root, vector<string> &lines) {
    int n = root-> rectangles.size();
    long area = INT_MAX, index = -1;
    for(int i=0;i<n;i++) {
        pair<int,int> p = root-> rectangles[i];
        if(p.first*p.second < area) {
            area = p.first*p.second;
            index = i;
        }
    }

    preOrderXY(root, index, false);
    postOrder(root, index, false, lines);

    return root-> rectangles[index];
}


int main(int argc, char *argv[])
{
    vector<string>nodes = readFile(argv[1]);
    TreeNode* root = constructTree(nodes);
    
    if(root == NULL) {
        return EXIT_FAILURE;
    }

    vector<string> lines;

    pair<int, int> p = singlePacking(root);
    cout<< p.first <<":" << p.second <<endl;
    char pstr[100];
    sprintf(pstr, "(%d,%d)\n", p.first, p.second);
    lines.push_back(pstr);
    writeFile(argv[2], lines);
    lines.clear();
    
    preOrderXY(root, 0, true);
    postOrder(root, 0, true, lines);
    writeFile(argv[3], lines);
    lines.clear();

    vector<string> newLines;
    doOptimalPacking(root);
    p = getOptimalPacking(root, newLines);

    cout<< p.first <<":" << p.second <<endl;
    sprintf(pstr, "(%d,%d)\n", p.first, p.second);
    lines.push_back(pstr);
    writeFile(argv[4], lines);
    lines.clear();

    writeFile(argv[5], newLines);
    newLines.clear();

    // Release resource
    deleteTree(root);

    return EXIT_SUCCESS;
}