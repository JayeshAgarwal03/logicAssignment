#include <iostream>
#include <sstream>
#include <fstream>
#include <unordered_set>
#include <string>
#include <stack>
#include <memory>
#include <algorithm>
#include <map>
#include <set>
#include <vector>
using namespace std;

// ==================== TASK 1: Infix to Prefix Converter ====================
// Function to assign a priority level to an operator.

int priorityofoperator(char c){
    
    if (c == '(') return 0; 
    else if( c == '~') return 1; 
    else if( c == '*') return 2;
    else if( c == '+') return 3;
    else if (c == '>') return 4;
    return -1;
}

// Function to swap '(' and ')' characters in a string.

void bracketswap(string &s){
    for (int i = 0; i < s.length(); i++){
        if (s[i] == ')') s[i] = '(';
        else if (s[i] == '(') s[i] = ')';
    }
}

string postfix(string &s) { 
    stack<char> st;
    string output = "";
    for (int i = 0; i < s.length(); i++) {
        char c = s[i];
        
        //  if character is an operand (letter or digit), add it directly to the output.
        if (isalpha(c) || isdigit(c)) { 
            output += c;
        } 
        //  if character is an opening bracket ,we  push it onto the stack.
        else if (c == '(') {
            st.push(c);
        } 
         // If character is a closing parenthesis.
        else if (c == ')') {
            // Pop operators from the stack and add to output until '(' is encountered.
            while (!st.empty() && st.top() != '(') {
                output += st.top();
                st.pop();
            }
            // Pop the '(' from the stack (don't add it to output).
            if (!st.empty()) st.pop(); 
        } 
        // If character is an operator.
        // pop stack while top operator has lower or equal precedence.
        else { 
            
            while (!st.empty() && priorityofoperator(st.top()) <= priorityofoperator(c) && st.top() != '(') {
                output += st.top();
                st.pop();
            }
            // Push the current operator onto the stack.
            st.push(c);
        }
    }
    // to pop the remaining operators 
    while (!st.empty()) {
        output += st.top();
        st.pop();
    }
    return output;
}


// Algorithm: Reverse the infix string-> Swap all '(' with ')' and vice-versa.
// -> Apply the standard Infix-to-Postfix conversion 
//  Reverse the resulting string to get the final Prefix expression.
string infixtoprefix(string s){
    reverse(s.begin(), s.end()); 
    bracketswap(s); 
    
    
    
    //Converts the modified string to a postfix form.
    string k = postfix(s);
    
    // reversal of k will yield prefix form 
    reverse(k.begin(), k.end());
    return k;
}

// Utility function to demonstrate the conversion.
void task1_InfixToPrefix(const string& infix) {
    cout << "\n========== TASK 1: Infix to Prefix Converter ==========\n";
    cout << "Infix:  " << infix << endl;
    cout << "Prefix: " << infixtoprefix(infix) << endl;
}


// ========= Task 2 : Converting a prefix expression to a Binary rooted tree =========
struct Node1 {
    char data;
    shared_ptr<Node1> left;
    shared_ptr<Node1> right;

//initialising the Node1
    Node1(char val) {
        data = val;
        left = right = nullptr;
    }
};


// to check whether we are dealing with some operator or numbers/variables
bool Valid_Op(char c) {
    return (c == '+' || c == '~' || c == '*' || c == '>');
}

// function for converting any prefix to binary parse tree
shared_ptr<Node1> Prefix_to_Tree(const string& prefix) {
    cout << "\n========== TASK 2: Converting a prefix expression to a Binary rooted tree ==========\n";
    stack<shared_ptr<Node1>> store;

    for (int i = prefix.length() - 1; i >= 0; i--) {
        char c = prefix[i];

        // to ignore any spaces in the given expression
        if (c == ' ') continue;

        //this is to deal with all the operands
        if (!Valid_Op(c)) {
            store.push(make_shared<Node1>(c));
        } 
        // this treats the unary operator different from the others as unary deals with a single operand
        else if (c == '~') { 
            if (store.empty()) {
                cout << "This prefix expression is not valid " << endl;
                return nullptr;
            }

        //this is the code for the rest of operators    
            shared_ptr<Node1> operand = store.top(); store.pop();
            shared_ptr<Node1> node1 = make_shared<Node1>(c);
            node1->left = operand;
            node1->right = nullptr;
            store.push(node1);
        } 
        else { 
            // to deal with the binary expressions
            if (store.size() < 2) {
                cout << "Invalid prefix expression " << endl;
                return nullptr;
            }
           
            // This ensures that the first popped node1 is the left one to maintain the order of the infix expression
            shared_ptr<Node1> left = store.top(); store.pop();
            // This ensures that the first popped node1 is the right one to maintain the order of the infix expression
            shared_ptr<Node1> right = store.top(); store.pop();
            
            shared_ptr<Node1> node1 = make_shared<Node1>(c);
            node1->left = left;
            node1->right = right;
            store.push(node1);
        }
    }
    //at this point only the main root of the expression should be remaining
    if (store.size() != 1) {
        cout << "Invalid prefix expression" << endl;
        return nullptr;
    }

    return store.top();
}


// ========= Task 3 : Traverse the binary rooted tree to output the infix expression =========


//Inorder traversal of the tree to output the infix expression

void Tree_to_Infix(shared_ptr<Node1> root) {
    
    //to deal with a null expression
    if (!root) return;

    // Ensures the order in unary operator case
    if (root->data == '~') {
        // Unary operator prints ~A
        cout << "(";
        cout << root->data; // 1. Print operator
        Tree_to_Infix(root->left); // 2. Print child
        cout << ")";
    } 
    // Ensures order for rest of the operators
    else if (Valid_Op(root->data)) {

        // Binary operator prints A * B
        cout << "(";
        //prints left subnode
        Tree_to_Infix(root->left);
        //prints operator
        cout << root->data;
        //prints right subnode
        Tree_to_Infix(root->right);
        cout << ")";
    } else {
        // If no opertaor is there , it just prints the corresponding data.
        cout << root->data;
    }
}


// ==================== TASK 4: Expression Tree Height ====================

struct TreeNode {
    char value;
    shared_ptr<TreeNode> left, right;
    TreeNode(char v) : value(v), left(nullptr), right(nullptr) {}
};

// Prints the tree structure in a visual hierarchical format with connecting lines.
// Right children are printed first (top), then left children (bottom).
void printTree(shared_ptr<TreeNode> node, string indent = "", bool isRight = true) {
    if (node == nullptr) return;
    
    cout << indent;
    cout << (isRight ? "└── " : "├── ");
    cout << node->value << endl;
    
    if (node->left || node->right) {
        if (node->right) {
            printTree(node->right, indent + (isRight ? "    " : "│   "), true);
        }
        if (node->left) {
            printTree(node->left, indent + (isRight ? "    " : "│   "), false);
        }
    }
}

// function to check precedance value, higher the precedence, higher the value.
int getPrecedence(char op) {
    if (op == '~') return 4; 
    if (op == '*') return 3;
    if (op == '+') return 2;
    if (op == '>') return 1; 
    return 0;
}

// for edgecase ~~ A is evaluated as ~(~A), not (~~)A.
bool isRightAssociative(char op) {
    return op == '~';  // only NOT is right-associative
}


// function to build a parse tree from infix notation using the Shunting Yard algorithm.
shared_ptr<TreeNode> buildTree(string infix) {
    stack<shared_ptr<TreeNode>> nodes;
    stack<char> ops;
    
    for (int i = 0; i < infix.length(); i++) {
        char c = infix[i];
        
        // Skip whitespace
        if (c == ' ') continue;
        
        // If variable, create a node and push to nodes stack
        if (isalpha(c)) {
            nodes.push(make_shared<TreeNode>(c));
        }

        // opening bracket: Push to operator stack
        else if (c == '(') {
            ops.push(c);
        }

        // Closing bracket: pop and build tree nodes until matching '(' is found
        else if (c == ')') {
            while (!ops.empty() && ops.top() != '(') {
                auto node = make_shared<TreeNode>(ops.top());
                ops.pop();
                
                // Unary operator (negation): only has right child
                if (node->value == '~') {
                    node->right = nodes.top();
                    nodes.pop();
                } else {
                    // Binary operator: pop two operands (right first, then left)
                    node->right = nodes.top();
                    nodes.pop();
                    node->left = nodes.top();
                    nodes.pop();
                }
                nodes.push(node);
            }
            ops.pop(); // remove '(' from ops
        }
        // Takes care of precedence and associativity rules
        else if (c == '~' || c == '+' || c == '*' || c == '>') {
            // Process operators with higher or equal precedence
            while (!ops.empty() && ops.top() != '(' &&
                   (getPrecedence(ops.top()) > getPrecedence(c) ||
                    (getPrecedence(ops.top()) == getPrecedence(c) && !isRightAssociative(c)))) {
                
                auto node = make_shared<TreeNode>(ops.top());
                ops.pop();
                
                if (node->value == '~') {
                    node->right = nodes.top();
                    nodes.pop();
                } else {
                    node->right = nodes.top();
                    nodes.pop();
                    node->left = nodes.top();
                    nodes.pop();
                }
                nodes.push(node);
            }
            ops.push(c);
        }
    }
    
    // Process remaining operators
    while (!ops.empty()) {
        auto node = make_shared<TreeNode>(ops.top());
        ops.pop();
        
        if (node->value == '~') {
            node->right = nodes.top();
            nodes.pop();
        } else {
            node->right = nodes.top();
            nodes.pop();
            node->left = nodes.top();
            nodes.pop();
        }
        nodes.push(node);
    }
    
    return nodes.top();
}

// Calculates the height of the tree recursively.
int calculateHeight(shared_ptr<TreeNode> root) {
    if (!root) return 0;
    return 1 + max(calculateHeight(root->left), calculateHeight(root->right));
}

// The main controller function.
void task4_ExpressionTreeHeight(const string& infix) {
    cout << "\n========== TASK 4: Expression Tree Height ==========\n";
    cout << "Infix expression: " << infix << endl;
    
    auto root = buildTree(infix);
    printTree(root);
    cout << "Height: " << calculateHeight(root) << endl;
}

// ==================== TASK 5: Truth Table ====================
struct Node2{
    char data;
    Node2* left;
    Node2* right;
    //Node initialisation
    Node2(char c) : data(c), left(nullptr), right(nullptr) {}
};

// to check whether we are dealing with some operator or operands
bool isOperator(char c) {
    return (c == '+' || c == '~' || c == '*' || c == '>');
}

//prefix to tree conversion
Node2* buildTreeFromPrefix1(const string& prefix) {
    //creating a stack that holds pointer to nodes
    cout << "\n========== TASK 5: Evaluate the truth value of a propositional logic formula ==========\n";
    stack<Node2*> store;

    // Traverse from right to left in the given prefix string
    for (int i = prefix.length() - 1; i >= 0; i--) {
        char c = prefix[i];

        // If the character is an operator
        if (isOperator(c)) {
            // Unary operator: we only pop one opperand which is the last one in stack and create a new node for the operator
            if (c == '~') {
                Node2* operand = store.top(); 
                store.pop();
                // Create a new node for the unary operator
                Node2* operatorNode = new Node2(c);
                operatorNode->left = operand;
                operatorNode->right = nullptr;
                store.push(operatorNode);
            }
            // Binary operators
            else {
                //check to prevent crash if stack is empty
                if (store.size() < 2) {
                    cout << "Invalid prefix: Missing operands for operator '" << c << "'" << endl;
                    return nullptr;
                }
                Node2* lnode = store.top(); 
                store.pop();
                Node2* rnode = store.top(); 
                store.pop();
                // Create a new node for the binary operator
                Node2* operatorNode = new Node2(c);
                operatorNode->left = lnode;
                operatorNode->right = rnode;
                store.push(operatorNode);
            }
        }
        // Operand
        else {
            store.push(new Node2(c));
        }
    }
    return store.top();
}


//Truth Value Evaluation
bool evaluateTruthValue(Node2* root, const map<char, bool>& truthValues) {
    //Base case: If the node is a leaf, return the truth value from the map
    if (!isOperator(root->data)) {
        return truthValues.at(root->data);
    }

    // Recursive case: for '~' operator
    if(root->data == '~') {
        bool leftValue = evaluateTruthValue(root->left, truthValues);
        return !leftValue;
    } 

    // Recursive case: for binary operators, calculating truth values for left and the right subtrees 
    bool leftValue = evaluateTruthValue(root->left, truthValues);
    bool rightValue = evaluateTruthValue(root->right, truthValues);

    //using switch:case for calculating the final truth value of the logic formula
    switch(root->data) {
        case '*': return leftValue&&rightValue;
        case '+': return leftValue||rightValue;
        case '>': return (!leftValue)||rightValue;
        default : throw runtime_error("Unkown operator encountered inthe tree");
    }
}

// Collecting unique variables for truth table by using Recusrion, and Set to store only unique items
void getVariables(Node2* root, set<char>& variables){
    if(!root) {
        return;
    }
    // If current node is an operand, store its data into set
    if(!isOperator(root->data)) {
        variables.insert(root->data);
    }
    // Recursive step for the right and left subtree of the current node, and store its data into set
    getVariables(root->left, variables);
    getVariables(root->right, variables);
}

// Printing Truth Table
void printTruthTable(Node2* root, const string& formulaString) {
    set<char> variableSet;
    getVariables(root, variableSet);

    // Converting the set into vector to access the variables using index
    vector<char> variables(variableSet.begin(), variableSet.end());
    int numVariables = variables.size();
    long long numRows = 1 << numVariables; //no. of rows = 2 to the power of number of variables

    //print headings for the columns
    cout<< "Truth Table: "<<endl;
    for(char c: variables) { 
        cout << "| " << c << "  ";
    }
    cout<< "|| " << formulaString << " |" <<endl;
    cout << string(numVariables* 5 + 4 + formulaString.length() + 2, '-') << endl;

    //printing rows
    for(long long i=0; i<numRows; i++) {
        map<char, bool> truthValues; //creating truthValue map for the current row
            
        for(int j=0; j<numVariables; j++) {
            bool value = ((i>>(numVariables-1-j)) & 1);
            truthValues[variables[j]]= value;
            cout<< "| " <<(value ? "T": "F")<< "  ";
        }

        bool result = evaluateTruthValue(root, truthValues);
        cout<< "||" << string(formulaString.length()/2, ' ')<< (result ? "T": "F")<< string(formulaString.length()/2,' ')<< " |"<< endl;
    }
}


// ==================== TASK 6: CNF Converter ====================
struct Node {
    char data;
    Node* left;
    Node* right;

    Node(char var, Node* leftptr = nullptr, Node* rightptr = nullptr)
        : data(var), left(leftptr), right(rightptr) {}
 
    ~Node() {
        // delete left;
        // delete right;
    }
};

// Helper function to build tree from prefix notation
Node* buildTreeFromPrefix(const string& prefix) {
    if (prefix.empty()) return nullptr;
    
    stack<Node*> st;
    
    for (int i = prefix.length() - 1; i >= 0; i--) {
        char c = prefix[i];
        
        if (isalpha(c)) {
            st.push(new Node(c));
        }
        else if (c == '~') {
            if (st.empty()) return nullptr;
            Node* operand = st.top();
            st.pop();
            st.push(new Node(c, operand, nullptr));
        }
        else if (c == '+' || c == '*' || c == '>') {
            if (st.size() < 2) return nullptr;
            Node* left = st.top();
            st.pop();
            Node* right = st.top();
            st.pop();
            st.push(new Node(c, left, right));
        }
    }
    
    return st.empty() ? nullptr : st.top();
}

// cloneNode creates copy of a subtree so that our desired functions 
// can modify formulas without changing the original parse tree
Node* cloneNode(Node* root) {
    if (!root) return nullptr;
    return new Node(root->data, cloneNode(root->left), cloneNode(root->right));
}

// Step 1: IMPLICATION FREE FORM (a>b becomes ~a+b)
Node* IMPL_FREE(Node* root) {
    if (root == nullptr) { return nullptr; }

    // Recursive step for children
    Node* left_res = IMPL_FREE(root->left);
    Node* right_res = IMPL_FREE(root->right);
    root->left = left_res;
    root->right = right_res;

    // If implication operator is found
    if (root->data == '>') {
        Node* lefttemporary = root->left;                     // store old left  eg:A
        Node* negationNode = new Node('~', lefttemporary, nullptr);    // creates a negation node eg:~A
        Node* orNode = new Node('+', negationNode, root->right); // creates an ornode   eg: ~A+B

        // Detach children before deletion
        root->left = nullptr;
        root->right = nullptr;
        // delete root;

        return orNode;
    }
    return root;
}

// Helper to check if a node is a literal (variable or negated variable)
bool isLiteral(Node* node) {
    if (!node) return false;
    // Variable
    if (isalpha(node->data) && !node->left && !node->right)
        return true;
    // Negated Variable
    if (node->data == '~' && node->left && isalpha(node->left->data) && !node->left->left)
        return true;
    return false;
}

// Step 2: Negation Normal Form (NNF)
// REQUIRES INPUT TO IT TO BE IN IMPLICATION FREE FORM 
Node* NNF(Node* root) {
    if (root == nullptr) return nullptr;

    // input is a literal 
    if (isLiteral(root))
        return new Node(root->data, cloneNode(root->left), cloneNode(root->right));

    
    if (root->data == '~') {
        Node* child = root->left;

        //NNF(~~phi) is same as NNF(phi)
        if (child->data == '~') {
            Node* newroot = child->left;
            child->left = nullptr; 
            root->left = nullptr;
            // delete root;
            // delete child;
            return NNF(newroot);
        }

        //NNF(~(phi_1*phi_2)) is same as NNF(~phi_1)+NNF(~phi_2)
        if (child->data == '*') {
            Node* not_left = new Node('~', child->left, nullptr); child->left = nullptr;
            Node* not_right = new Node('~', child->right, nullptr); child->right = nullptr;

            Node* newroot = new Node('+', NNF(not_left), NNF(not_right));
            root->left = nullptr; 
            // delete root; delete child;
            return newroot;
        }

         //NNF(~(phi_1+phi_2)) is same as NNF(~phi_1)*NNF(~phi_2)
        if (child->data == '+') {
            Node* not_left = new Node('~', child->left, nullptr); child->left = nullptr;
            Node* not_right = new Node('~', child->right, nullptr); child->right = nullptr;

            Node* newroot = new Node('*', NNF(not_left), NNF(not_right));
            root->left = nullptr; 
            // delete root; delete child;
            return newroot;
        }

        // Standard negation case
        return new Node('~', NNF(child), nullptr);
    }

    // NNF(phi_1 * phi_2)is NNF(phi_1) *NNF(phi_2)
    if (root->data == '*') {
        Node* newroot = new Node('*', NNF(root->left), NNF(root->right));
        root->left = nullptr; root->right = nullptr; 
        // delete root;
        return newroot;
    }

     // NNF(phi_1 + phi_2)is NNF(phi_1)+NNF(phi_2)
    if (root->data == '+') {
        Node* newroot = new Node('+', NNF(root->left), NNF(root->right));
        root->left = nullptr; root->right = nullptr; 
        // delete root;
        return newroot;
    }

    return new Node(root->data, cloneNode(root->left), cloneNode(root->right));
}


// Step 3: Distribute OR over AND for CNF
Node* DISTR(Node* left, Node* right) {
    // example: // (A + (B * C))  will become  ((A + B) * (A + C))
    if (left && left->data == '*') { // left part contains conjunction
        Node* distLeft = DISTR(left->left, right);
        Node* distRight = DISTR(left->right, right);
        return new Node('*', distLeft, distRight);
    }
    if (right && right->data == '*') {// right part contains conjunction
        Node* distLeft = DISTR(left, right->left);
        Node* distRight = DISTR(left, right->right);
        return new Node('*', distLeft, distRight);
    }

    return new Node('+', left, right);
}

// Step 4: CNF transformation
Node* CNF(Node* root) {
    if (root == nullptr || isLiteral(root))
        return root;

    // cnf of conjunction
    if (root->data == '*') {
        root->left = CNF(root->left);
        root->right = CNF(root->right);
        return root;
    }
    // cnf of disjunction
    if (root->data == '+') {
        Node* leftCNF = CNF(root->left);
        Node* rightCNF = CNF(root->right);

        root->left = nullptr;
        root->right = nullptr;
        // delete root;
        return DISTR(leftCNF, rightCNF);
    }

    return root;
}

// converting to cnf if basically applying cnf on nnf on impl_free form of the formula ;
Node* CONVERTTOCNF(Node* root) {
    Node* implFreeRoot = IMPL_FREE(root);
    Node* nnfRoot = NNF(implFreeRoot);
    Node* cnfRoot = CNF(nnfRoot);
    return cnfRoot;
}

// Print tree in infix notation
void printTreeInfix(Node* root) {
    if (!root) return;
    
    // For infix output
    if (root->left || root->right) {
        if (root->data == '~') {
            cout << root->data;
            if (root->left && (root->left->data == '+' || root->left->data == '*' || root->left->data == '>')) {
                cout << "(";
                printTreeInfix(root->left);
                cout << ")";
            } else {
                printTreeInfix(root->left);
            }
        } else {
            cout << "(";
            printTreeInfix(root->left);
            cout << root->data;
            printTreeInfix(root->right);
            cout << ")";
        }
    } else {
        cout << root->data;
    }
}

void task6_CNFConverter(const string& prefix_formula, const string& infix_formula) {
    cout << "\n========== TASK 6: CNF Converter ==========\n";
    
    // Build tree from prefix for CNF conversion
    Node* cnf_root_input = buildTreeFromPrefix(prefix_formula);
    if (!cnf_root_input) { 
        cerr << "Error: Could not build tree for CNF conversion." << endl; 
        return; 
    }

    Node* cnfTree = CONVERTTOCNF(cnf_root_input);

    cout << "Original Formula (Infix): " << infix_formula << endl;
    cout << "CNF Infix Representation: ";
    printTreeInfix(cnfTree);
    cout << endl;

    delete cnfTree;
}

// ==================== TASK 7: CNF Clause Validator ====================

// A clause is satisfiable if it contains at least one pair of complementary literals (e.g., x and ~x).
// Returns true if the clause contains both a literal and its negation (making the clause always true).
// Returns false if no complementary pair is found (clause may be false under some assignment).
bool isClauseValid(const string& clause_line) {
    stringstream ss(clause_line);
    unordered_set<int> literals;
    int literal;
    while (ss >> literal && literal != 0) {
        if (literals.count(-literal)) {
            return true;
        }
        literals.insert(literal);
    }
    return false;
}

// The main controller function.
void task7_CNF_Validator(const string& filepath) {
    cout << "\n========== TASK 7: CNF Clause Validator ==========\n";
    
    // This function reads the .cnf file
    ifstream file(filepath);
    if (!file.is_open()) {
        cout << "Error: Could not open file " << filepath << endl;
        return;
    }
    
    // Initialize the variables that count the number of valid and invalid clauses.
    string line;
    int valid_clauses = 0;
    int invalid_clauses = 0;
    
    // Loop through each line in the file.
    while (getline(file, line)) {

        // Ignore lines starting from 'c' or 'p'.
        if (line.empty() || line[0] == 'c' || line[0] == 'p') {
            continue;
        }

        // Check if each clause is valid or not and count the number of valid and invalid clauses.
        if (isClauseValid(line)) {
            valid_clauses++;
        } else {
            invalid_clauses++;
        }
    }
    
    // Check if all clauses are valid (CNF formula is valid)
    if (invalid_clauses == 0 && valid_clauses > 0) {
        cout << "\nCNF Formula: VALID (All clauses are satisfiable)" << endl;
    } else if (valid_clauses == 0 && invalid_clauses == 0) {
        cout << "\nCNF Formula: No clauses found" << endl;
    } else {
        cout << "\nCNF Formula: INVALID (Contains unsatisfiable clauses)" << endl;
    }
    
    file.close();
    
    cout << "Valid Clauses:   " << valid_clauses << endl;
    cout << "Invalid Clauses: " << invalid_clauses << endl;
}

// ==================== MAIN FUNCTION ====================
int main() {
    // Example formula from assignment
    string infix_formula = "((p*~q)>r)>(~r>(~p+q))";
    string prefix_formula = ">>*p~qr>~r+~pq";

    
    // Task 1: Infix to Prefix Converter
    task1_InfixToPrefix("q*b+(c>d)");

    // Task 2: Converting a prefix expression to a Binary rooted tree
    shared_ptr<Node1> root = Prefix_to_Tree(prefix_formula);
    cout << "Prefix : " << prefix_formula << endl;

    // Task 3: Converting a prefix expression to a Binary rooted tree
    // to prevent error in case of nullptr
    cout << "\n========== TASK 3: Traverse the binary rooted tree to output the infix expression ==========\n";
    if (root != nullptr) {
        Tree_to_Infix(root);
        cout << endl;
    }
    
    // Task 4: Expression Tree Height
    task4_ExpressionTreeHeight("(A*B)+(~C*D)");

    //Task 5: Evaluate the truth value of a propositional logic formula
    Node2* root1 = buildTreeFromPrefix1(prefix_formula);
    map<char, bool> truthValues;
    truthValues['p'] = true;
    truthValues['q'] = false;
    truthValues['r'] = true;
    bool result = evaluateTruthValue(root1, truthValues);
    cout << "Truth Value for the propositional logic formula is: " << (result ? "True" : "False") << endl;
    printTruthTable(root1, infix_formula);
    
    // Task 6: CNF Converter
    task6_CNFConverter(prefix_formula, infix_formula);
    
    // Task 7: CNF Clause Validator
    task7_CNF_Validator("/Users/jayesh/logic_assignment/tasks/test.cnf");
    
    return 0;
}