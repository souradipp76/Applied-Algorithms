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

void readBinaryTopoFile(const string filename) {
    FILE* f = fopen(filename.c_str(), "rb");
    int id;
    double cap;
    while () {
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
    double cap[] = {}, llen[] = {0}, rlen[] = {0};
    for (auto val:vals) {
        char* v = val.c_str();
        int len = sscanf(v, "%d(%le)", &id, &cap[0]);
        if(len> 0) {
            fwrite(&id, sizeof(int), 1, f);
            fwrite(&cap, sizeof(double), 1, f);
        } else {
            len = sscanf(v, "(%le %le %d)", &llen, &rlen, &count);
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
        double sink_capacitance = 0.0;
        double downstream_capacitance = 0.0;
        double llen = 0.0;
        double rlen = 0.0;
        double delay = 0.0;
        int inv_count = 0;
        TreeNode* left = NULL;
        TreeNode* right = NULL;
    
    TreeNode(string id, double cap) : nodeId(id), sink_capacitance(cap) {}
    TreeNode(string id, double res, double cap, TreeNode *l, TreeNode *r) : nodeId(id), resistance(res), capacitance(cap), left(l), right(r) {}
};

TreeNode* construct_tree(vector<string> nodes, double wireR, double wireC, double Co, double Ro) {
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
    head->inv_count = 1;
    return head;
}

void preorder(TreeNode* root, vector<string> &lines) {
    if(root == NULL || root->nodeId == "dummy") {
        return;
    }

    char line[100];
    if(root->left && root->right) {
        sprintf(line, "(%.10le %.10le)\n", root->llen, root-> rlen);
        // cout<<line<<endl;
        lines.push_back(line);
    } else {
        sprintf(line, "%d(%.10le)\n", stoi(root->nodeId), root-> sink_capacitance);
        // cout<<line<<endl;
        lines.push_back(line);
    }
    
    preorder(root->left, lines);
    preorder(root->right, lines);
}

void postorder(TreeNode* root, vector<string> &lines) {
    if(root == NULL || root->nodeId == "dummy") {
        return;
    }

    postorder(root->left, lines);
    postorder(root->right, lines);

    char line[100];
    if(root->left && root->right) {
        sprintf(line, "(%.10le %.10le %d)\n", root->llen, root-> rlen, root->inv_count);
        //cout<<line<<endl;
        lines.push_back(line);
    } else {
        sprintf(line, "%d(%.10le)\n", stoi(root->nodeId), root-> sink_capacitance);
        //cout<<line<<endl;
        lines.push_back(line);
    }

}

void deleteTree(TreeNode* root) {
    if(root == NULL) {
        return;
    }

    deleteTree(root->left);
    deleteTree(root->right);

    delete(root);
}

void calculateDownstreamCapacitance(TreeNode* root) {
    if(!root) {
        return;
    }

    TreeNode* l = root->left;
    TreeNode* r = root->right;
    root->downstream_capacitance = root->capacitance;

    if(root->sink_capacitance >0) {
        return;
    }

    if(l && root->llen >=0) {
        calculateDownstreamCapacitance(l);
        double c_l = l->downstream_capacitance;
        root->downstream_capacitance += c_l;
    }
    
    if(r && root->rlen >=0){
        calculateDownstreamCapacitance(r);
        double c_r = r->downstream_capacitance;
        root->downstream_capacitance +=c_r;
    } 

    //cout<<root->nodeId<<":"<< root->downstream_capacitance<<endl;
} 

void calculateElmoreDelay(TreeNode* root) {
    if(!root) {
        return;
    }

    if(root->sink_capacitance >0) {
        return;
    }

    TreeNode* l = root->left;
    TreeNode* r = root->right;

    if(l && root->llen >=0) {
        l-> delay = root->delay + l->resistance*l->downstream_capacitance;
        // cout<<l->nodeId<<":"<<l->resistance<<":"<<l->downstream_capacitance<<endl;
        // cout<<l->delay<<endl;
        calculateElmoreDelay(l);
    }
    
    if(r && root->rlen >=0){
        r-> delay = root->delay + r->resistance*r->downstream_capacitance;
        // cout<<r->nodeId<<":"<<r->resistance<<":"<<r->downstream_capacitance<<endl;
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
    double Co = invParams[1];
    double Ro = invParams[2];

    double re = wireParams[0];
    double ce = wireParams[1];

    double a  = 0.5*re*ce*L*L;
    double b = re*L*(Co - 0.5*ce*L) + 0.5*ce*Ro*L;
    double c = Ro*(Co - 0.5*ce*L) - t;

    double D = b*b - 4*a*c;
    cout<<"D:"<<D<<endl;
    if(D < 0) {
        return -1;
    }

    double x = (-b + sqrt(D))/(2.*a);
    cout<<x<<endl;
    if(x>=0 && x<=1) {
        return x;
    } else if (x>1) {
        return 1;
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

    if(x == 1) {
        // Place inverter node at parent node
        parent->inv_count++;
        parent->sink_capacitance = parent->inv_count*Ci;
        parent->capacitance = parent->sink_capacitance;
        parent->downstream_capacitance = Co + 0.5*ce*L + node->downstream_capacitance;

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
        treeNode->downstream_capacitance = Co + 0.5*x*ce*L + node->downstream_capacitance;

        node->resistance = Ro+x*L*re;
        node->capacitance += Co - 0.5*(1-x)*ce*L;

        treeNode->inv_count++;
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
        cout<<l->nodeId<<endl;
        cout<<"Left Delay:"<<left_delay<<endl;
        cout<<"LDelay:"<<ldelay<<endl;
        while(ldelay + left_delay > time_constraint) {
            double remain = time_constraint - left_delay;
            double x = calculateInsertionLength(invParams, wireParams, root ->llen, root->downstream_capacitance, remain);

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
        cout<<r->nodeId<<endl;
        cout<<"Right Delay:"<<right_delay<<endl;
        cout<<"RDelay:"<<rdelay<<endl;
        while(rdelay + right_delay > time_constraint) {
            double remain = time_constraint - right_delay;
            double x = calculateInsertionLength(invParams, wireParams, root ->rlen, root->downstream_capacitance, remain);
            
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

int main(int argc, char *argv[])
{
    double time_constraint = stod(argv[1]);
    cout<<time_constraint<<endl;

    vector<string>invParamLines = readFile(argv[2]);
    cout<< invParamLines.size() <<endl;
    vector<string> invParams = splitLine(invParamLines[0]);
    double inC = stod(invParams[0]), outC = stod(invParams[1]), invR = stod(invParams[2]);
    cout<<inC<< ":" << outC << ":" << invR << endl;

    vector<string> wireParamLines = readFile(argv[3]);
    cout<< wireParamLines.size() <<endl;
    vector<string> wireParams = splitLine(wireParamLines[0]);
    double wireR = stod(wireParams[0]), wireC = stod(wireParams[1]);
    cout<<wireR<< ":" << wireC <<endl;

    vector<string>nodes = readFile(argv[4]);
    TreeNode* root = construct_tree(nodes, wireR, wireC, outC, invR);
    
    if(root == NULL) {
        return EXIT_FAILURE;
    }

    vector<string> topology;
    preorder(root, topology);
    writeFile(argv[5], topology);

    cout<<"Calculate downstream capacitance ..."<<endl;
    calculateDownstreamCapacitance(root);
    cout<<endl;

    cout<<"Calculate elmore delay ..."<<endl;
    // Initialize root delay
    root->delay = root->resistance*root->downstream_capacitance;
    calculateElmoreDelay(root);

    vector<pair<int, double>> elmoreDelays;
    getElmoreDelay(root, elmoreDelays);
    writeBinaryFile(argv[6], elmoreDelays);

    // cout<<endl;
    // readBinaryFile(argv[6], 5);
    // cout<<endl;

    cout<<endl;
    string filename = "pa1_examples\\examples\\5.elmore";
    readBinaryFile(filename, 5);
    cout<<endl;

    vector<double> wParams;
    for(auto x: wireParams) {
        wParams.push_back(stod(x));
    }

    vector<double> iParams;
    for(auto x: invParams) {
        iParams.push_back(stod(x));
    }

    bool valid = true;
    insertInverters(root, root, iParams, wParams, time_constraint, valid);
    if(!valid) {
        return EXIT_FAILURE;
    }
    vector<string> invTopology;
    postorder(root, invTopology);
    writeFile(argv[7], invTopology);
    writeBinaryFile(argv[8], invTopology);


    readBinaryFile(argv[8], invTopology.size());

    displayElmoreDelay(root);

    deleteTree(root);

    return EXIT_SUCCESS;
}