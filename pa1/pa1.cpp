#include <iostream>
#include <vector>
#include <string>
#include <stack>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <unordered_map>
#include <cmath>

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

void readBinaryFile(const string filename, int n) {
    FILE* f = fopen(filename.c_str(), "rb");
    for (int i = 0; i < n; i++) {
        int id;
        double t;
        fread(&id, sizeof(int), 1, f);
        fread(&t, sizeof(double), 1, f);
        printf("%d: %.10le\n", id, t);
    }
    fclose(f);
}

void writeBinaryFile(const string filename, vector<pair<int, double>> vals) {
    FILE* f = fopen(filename.c_str(), "wb");
    int n = vals.size();
    for (int i = 0; i < n; i++) {
        int id[] = {vals[i].first};
        double t[] = {vals[i].second};
        fwrite(&id, sizeof(int), 1, f);
        fwrite(&t, sizeof(double), 1, f);
    }
    fclose(f);
}

void writeBinaryTopoFile(const string filename, vector<string> vals) {
    FILE* f = fopen(filename.c_str(), "wb");
    int id[] = {0}, count[] = {0};
    double cap[] = {0}, llen[] = {0}, rlen[] = {0};
    for (auto val:vals) {
        const char* v = val.c_str();
        int len = sscanf(v, "%d(%le)", &id[0], &cap[0]);
        if(len> 0) {
            fwrite(&id, sizeof(int), 1, f);
            fwrite(&cap, sizeof(double), 1, f);
        } else {
            len = sscanf(v, "(%le %le %d)", &llen[0], &rlen[0], &count[0]);
            id[0] = -1;
            fwrite(&id, sizeof(int), 1, f);
            fwrite(&llen, sizeof(double), 1, f);
            fwrite(&rlen, sizeof(double), 1, f);
            fwrite(&count, sizeof(int), 1, f);
        }
    }
    fclose(f);
}

class TreeNode {
    public:
        string nodeId;
        double resistance = 0.0;
        double capacitance = 0.0;
        double sinkCapacitance = 0.0;
        double downstreamCapacitance = 0.0;
        double llen = 0.0;
        double rlen = 0.0;
        double delay = 0.0;
        int invCount = 0;
        TreeNode* left = NULL;
        TreeNode* right = NULL;
    
    TreeNode(string id, double cap) : nodeId(id), sinkCapacitance(cap) {}
    TreeNode(string id, double res, double cap, TreeNode *l, TreeNode *r) : nodeId(id), resistance(res), capacitance(cap), left(l), right(r) {}
};

TreeNode* constructTree(vector<string> nodes, double wireR, double wireC, double Co, double Ro) {
    TreeNode *head = NULL;  
    bool success = true;

    stack<TreeNode*> st;
    for(auto &node: nodes) {
        bool isLeaf = false;
        int id = 0;
        double c = 0;
        int done = sscanf(node.c_str(),"%d(%le)", &id, &c);
        if (done > 0) {
            isLeaf = true;
        }

        if (isLeaf) {
            TreeNode* treeNode = new TreeNode(to_string(id), c);
            treeNode->capacitance+=c;
            st.push(treeNode);
        } else {
            double llen = 0, rlen = 0;
            done = sscanf(node.c_str(),"(%le %le)", &llen, &rlen);
            if (done <0 || st.size() < 2) {
                success = false;
                break;
            };

            TreeNode* r = st.top();
            st.pop();
            TreeNode* l = st.top();
            st.pop();

            l->capacitance += wireC*llen*0.5;
            l->resistance += wireR*llen;
            r->capacitance += wireC*rlen*0.5;
            r->resistance += wireR*rlen;

            string newId = l->nodeId +"-"+r->nodeId;
            TreeNode *treeNode = new TreeNode(newId, 0, 0.5*wireC*(llen+rlen), l, r);
            treeNode -> llen = llen;
            treeNode -> rlen = rlen;
            st.push(treeNode);
        }
    }

    if (!success) {
        return NULL;
    }

    if(st.empty()) {
        return NULL;
    }

    head = st.top();
    head->capacitance += Co;
    head->resistance = Ro;
    head->invCount = 1;
    return head;
}

void deleteTree(TreeNode* root) {
    if(root == NULL) {
        return;
    }

    deleteTree(root->left);
    deleteTree(root->right);

    delete(root);
}

void preOrder(TreeNode* root, vector<string> &lines) {
    if(root == NULL || root->nodeId == "dummy") {
        return;
    }

    char line[100];
    if(root->left && root->right) {
        sprintf(line, "(%.10le %.10le)\n", root->llen, root-> rlen);
        // cout<<line<<endl;
        lines.push_back(line);
    } else {
        sprintf(line, "%d(%.10le)\n", stoi(root->nodeId), root-> sinkCapacitance);
        // cout<<line<<endl;
        lines.push_back(line);
    }
    
    preOrder(root->left, lines);
    preOrder(root->right, lines);
}

void postOrder(TreeNode* root, vector<string> &lines) {
    if(root == NULL || root->nodeId == "dummy") {
        return;
    }

    postOrder(root->left, lines);
    postOrder(root->right, lines);

    char line[100];
    if(root->left && root->right) {
        sprintf(line, "(%.10le %.10le %d)\n", root->llen, root-> rlen, root->invCount);
        //cout<<line<<endl;
        lines.push_back(line);
    } else {
        sprintf(line, "%d(%.10le)\n", stoi(root->nodeId), root-> sinkCapacitance);
        //cout<<line<<endl;
        lines.push_back(line);
    }

}

void calculateDownstreamCapacitance(TreeNode* root) {
    if(!root) {
        return;
    }

    TreeNode* l = root->left;
    TreeNode* r = root->right;

    root->downstreamCapacitance = root->capacitance;

    if(root->sinkCapacitance > 0) {
        return;
    }

    if(l && root->llen >=0) {
        calculateDownstreamCapacitance(l);
        double c_l = l->downstreamCapacitance;
        root->downstreamCapacitance += c_l;
    }
    
    if(r && root->rlen >=0){
        calculateDownstreamCapacitance(r);
        double c_r = r->downstreamCapacitance;
        root->downstreamCapacitance +=c_r;
    } 

    //cout<<root->nodeId<<":"<< root->downstreamCapacitance<<endl;
} 

void calculateElmoreDelay(TreeNode* root) {
    if(!root) {
        return;
    }

    TreeNode* l = root->left;
    TreeNode* r = root->right;

    if(l && root->llen >=0) {
        l-> delay = root->delay + l->resistance*l->downstreamCapacitance;
        // cout<<l->nodeId<<":"<<l->resistance<<":"<<l->downstreamCapacitance<<endl;
        // cout<<l->delay<<endl;
        calculateElmoreDelay(l);
    }
    
    if(r && root->rlen >=0){
        r-> delay = root->delay + r->resistance*r->downstreamCapacitance;
        // cout<<r->nodeId<<":"<<r->resistance<<":"<<r->downstreamCapacitance<<endl;
        // cout<<r->delay<<endl;
        calculateElmoreDelay(r);
    }   
} 

void getElmoreDelay(TreeNode* root, vector<pair<int, double>> &delays){
    if(!root->left && !root->right) {
        pair<int, double> p = {stoi(root->nodeId), root->delay};
        delays.push_back(p);
    }

    if(root->left) {
        getElmoreDelay(root->left, delays);
    }

    if (root-> right) {
        getElmoreDelay(root->right, delays);
    }
}

void displayElmoreDelay(TreeNode* root){
    cout<<root->nodeId<<":"<< root->delay<<endl;
    if(!root->left && !root->right) {
        return;
    }

    if(root->left) {
        displayElmoreDelay(root->left);
    }

    if (root-> right) {
        displayElmoreDelay(root->right);
    }
}

double calculateInsertionLength(vector<double> invParams, vector<double> wireParams, double L, double C, double t){
    //double Ci = invParams[0];
    //double Co = invParams[1];
    double Ro = invParams[2];

    double re = wireParams[0];
    double ce = wireParams[1];

    double a = 0.5*re*ce*L*L;
    double b = re*L*(C - 0.5*ce*L) + 0.5*ce*Ro*L;
    double c = Ro*(C - 0.5*ce*L) - t;

    double D = b*b - 4*a*c;
    if(D < 0) {
        return -1;
    }

    double x = (-b + sqrt(D))/(2.*a);
    double y = (-b - sqrt(D))/(2.*a);
    //cout<<x<<endl;
    if(x>=0 && x<=1) {
        return x;
    } 

    if(y>=0 && y<=1) {
        return y;
    }

    if(x>1 || y>1) {
        return 1.1;
    }

    return -1;
}


TreeNode* insertInvNode(TreeNode* head, 
                        TreeNode* node, 
                        TreeNode* parent, 
                        vector<double> &invParams, 
                        vector<double> &wireParams, 
                        double x) {
    string id = node->nodeId;
    string pid = parent -> nodeId;

    double Ci = invParams[0];
    double Co = invParams[1];
    double Ro = invParams[2];

    double re = wireParams[0];
    double ce = wireParams[1];

    bool isLeft = parent->left == node;
    double L = isLeft ? parent -> llen : parent ->rlen;

    if(x > 1) {
        // Place inverter node at parent node
        parent->invCount++;
        parent->sinkCapacitance = parent->invCount*Ci;
        parent->capacitance = parent->sinkCapacitance;
        parent->downstreamCapacitance = Co + 0.5*ce*L + node->downstreamCapacitance;

        calculateDownstreamCapacitance(head);
        calculateElmoreDelay(head);
        return parent;

    } else {
        // Create inverter node
        TreeNode* treeNode = new TreeNode(pid+"-inv-"+id, Ci);
        treeNode->llen = x*L;
        treeNode->rlen = -1;
        treeNode->resistance = (1-x)*L*re;
        treeNode->capacitance = Ci + 0.5*ce*L*(1-x);
        treeNode->downstreamCapacitance = Co + 0.5*x*ce*L + node->downstreamCapacitance;

        node->resistance = Ro+x*L*re;
        node->capacitance -= 0.5*(1-x)*ce*L;
        node->downstreamCapacitance -= 0.5*(1-x)*ce*L;

        treeNode->invCount++;
        treeNode->left = node;
        treeNode->right = new TreeNode("dummy", 0);
        
        // Insert Inverter node
        if(isLeft) {
            parent -> left = treeNode;
            parent -> llen = (1-x)*L;
        } else {
            parent -> right = treeNode;
            parent -> rlen = (1-x)*L;
        }

        calculateDownstreamCapacitance(head);
        calculateElmoreDelay(head);

        return treeNode;
    }
}

double insertInverters(TreeNode* head, 
                        TreeNode* root, 
                        vector<double> &invParams, 
                        vector<double> &wireParams, 
                        double time_constraint, 
                        bool &valid) {
    if(!root->left && !root->right) {
        return 0;
    }

    if(!valid) {
        return 0;
    }

    TreeNode* l = root->left;
    TreeNode* r = root->right;

    double ldelay = 0;
    double rdelay = 0;
    double left_delay = 0;
    double right_delay = 0;

    if(l) {
        left_delay = insertInverters(head, l, invParams, wireParams, time_constraint, valid);
        if(!valid) {
            return 0;
        }
        ldelay = l->delay - root->delay;
        // cout<<l->nodeId<<endl;
        // cout<<"Left Delay:"<<left_delay<<endl;
        // cout<<"LDelay:"<<ldelay<<endl;
        while(ldelay + left_delay > time_constraint) {
            double remain = time_constraint - left_delay;
            double x = calculateInsertionLength(invParams, wireParams, root ->llen, l->downstreamCapacitance, remain);

            if(x == -1) {
                valid = false;
                break;
            }

            TreeNode* temp = insertInvNode(head, l, root, invParams, wireParams, x);
            l = temp; 
            ldelay = l->delay - root->delay;
            left_delay = 0;
        }
    }

    if(!valid) {
        return 0;
    }

    if(r) {
        right_delay = insertInverters(head, r, invParams, wireParams, time_constraint, valid);
        if(!valid) {
            return 0;
        }
        rdelay = r->delay - root->delay;
        // cout<<r->nodeId<<endl;
        // cout<<"Right Delay:"<<right_delay<<endl;
        // cout<<"RDelay:"<<rdelay<<endl;
        while(rdelay + right_delay > time_constraint) {
            double remain = time_constraint - right_delay;
            double x = calculateInsertionLength(invParams, wireParams, root ->rlen, r->downstreamCapacitance, remain);
            
            if(x == -1) {
                break;
            }

            TreeNode* temp = insertInvNode(head, r, root, invParams, wireParams, x);
            r = temp; 
            rdelay = r->delay - root->delay;
            right_delay = 0;
        }
    }

    if(!valid) {
        return 0;
    }

    return max(ldelay + left_delay, rdelay + right_delay);
}

bool validateInverterInsertion(TreeNode* root, int pathCount) {
    if(!root->left && !root->right) {
        if((pathCount + root->invCount)&1) {
            //cout<<root->nodeId<<endl;
            return false;
        }
    }

    if(root->left && root->llen !=-1) {
        bool leftValid = validateInverterInsertion(root->left, pathCount+root->invCount);
        if(!leftValid) {
            return false;
        }
    }

    if(root->right && root->rlen !=-1) {
        bool rightValid = validateInverterInsertion(root->right, pathCount+root->invCount);
        if(!rightValid) {
            return false;
        }
    }

    return true;
}

bool correctInverterInsertion(TreeNode* head, 
                        TreeNode* root,
                        vector<double> &invParams, 
                        vector<double> &wireParams,
                        int pathCount, 
                        bool &inserted) {
    if(!root->left && !root->right) {
        if((pathCount + root->invCount)&1) {
            //cout<<root->nodeId<<endl;
            return false;
        }
    }

    if(root->left && root->llen !=-1) {
        bool leftValid = correctInverterInsertion(head, root->left, invParams, wireParams, pathCount+root->invCount, inserted);
        if(!leftValid) {
            if(!inserted) {
                insertInvNode(head, root->left, root, invParams, wireParams, 1);
                inserted = true;
            }
            return false;
        }
    }

    if(root->right && root->rlen !=-1) {
        bool rightValid = correctInverterInsertion(head, root->right, invParams, wireParams, pathCount+root->invCount, inserted);
        if(!rightValid) {
            if(!inserted) {
                insertInvNode(head, root->right, root, invParams, wireParams, 1);
                inserted = true;
            }
            return false;
        }
    }

    return true;
}

int main(int argc, char *argv[])
{
    double time_constraint = stod(argv[1]);
    //cout<<time_constraint<<endl;

    vector<string>invParamLines = readFile(argv[2]);
    vector<string> invParams = splitLine(invParamLines[0]);
    vector<double> iParams;
    for(auto x: invParams) {
        iParams.push_back(stod(x));
    }

    vector<string> wireParamLines = readFile(argv[3]);
    vector<string> wireParams = splitLine(wireParamLines[0]);
    vector<double> wParams;
    for(auto x: wireParams) {
        wParams.push_back(stod(x));
    }

    vector<string>nodes = readFile(argv[4]);
    TreeNode* root = constructTree(nodes, wParams[0], wParams[1], iParams[1], iParams[2]);
    
    if(root == NULL) {
        return EXIT_FAILURE;
    }

    // First Output
    vector<string> topology;
    preOrder(root, topology);
    writeFile(argv[5], topology);

    // Second Output
    calculateDownstreamCapacitance(root);

    // Initialize root delay
    root->delay = root->resistance*root->downstreamCapacitance;
    calculateElmoreDelay(root);

    vector<pair<int, double>> elmoreDelays;
    getElmoreDelay(root, elmoreDelays);
    writeBinaryFile(argv[6], elmoreDelays);

    // cout<<endl;
    // readBinaryFile(argv[6], 5);
    // cout<<endl;

    // string filename = "pa1_examples\\examples\\5.elmore";
    // readBinaryFile(filename, 5);
    // cout<<endl;

    // Third Output
    bool valid = true;
    insertInverters(root, root, iParams, wParams, time_constraint, valid);
    if(!valid) {
        return EXIT_FAILURE;
    }

    // bool isNonInverting = validateInverterInsertion(root, 0);
    // cout<<isNonInverting<<endl;

    //Make tree non-inverting
    bool inserted = false;
    correctInverterInsertion(root, root, iParams, wParams, 0, inserted);

    vector<string> invTopology;
    postOrder(root, invTopology);
    writeFile(argv[7], invTopology);

    // Fourth Output
    writeBinaryTopoFile(argv[8], invTopology);

    //displayElmoreDelay(root);

    // Release resource
    deleteTree(root);

    return EXIT_SUCCESS;
}