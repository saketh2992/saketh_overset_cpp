#include "adt.h"
#include <iostream>
#include <fstream>
#include <queue>
#include <stack>
#include <algorithm>

using namespace std;

// node_adt implementation
node_adt::node_adt(const size_t elementIdx, const array<su2double, SU2_BBOX_SIZE> &bBoxCoords) 
    : elementIndex(elementIdx), bBoxCoordinates(bBoxCoords), left(nullptr), right(nullptr) {}

// ADT implementation
ADT::ADT() : root(nullptr), treeHeirarchy(0) {}

void ADT::insertNode(node_adt *node) {
    node_adt *current = root;
    node_adt *parent = nullptr;
    size_t elementHeirarchy = 0, i = 0;

    /* Traverse the tree to find the insertion point */
    while (current != nullptr) {
        parent = current;
        /* Heirarchy cycling -> (x_min, y_min, z_min, x_max, y_max, z_max)*/
        i = elementHeirarchy % 4; /* For 2d cycle with 4 values*/
        /* node with equal values are put on left branch*/
        if (current->bBoxCoordinates[i] < node->bBoxCoordinates[i]) {
            current = current->right;
        } else {
            current = current->left;
        }
        elementHeirarchy = elementHeirarchy + 1;
    }
    if (parent == nullptr) {
        root = node;
    } else if (parent->bBoxCoordinates[i] < node->bBoxCoordinates[i]) {
        parent->right = node;
    } else {
        parent->left = node;
    }

    current = node;
    treeHeirarchy = max(treeHeirarchy, elementHeirarchy);
    // cout << "Element added " << current->elementIndex << " at h = " << heirarchy << endl;
}

vector<size_t> ADT::searchADT(const array<su2double, SU2_BBOX_SIZE> &testBBox) const {
    /*return a vector of indices of elements which intersect with the test bounding box*/
    size_t currHeirarchy = 0, i = 0;
    vector<size_t> intersectingBBox;
    node_adt *current = nullptr;
    bool intersect = true;

    stack<pair<node_adt *, size_t> > searchQ;
    searchQ.push(make_pair(root, 0));

    while (searchQ.empty() == false) {
        current = searchQ.top().first;
        currHeirarchy = searchQ.top().second;
        searchQ.pop();

        while (current != nullptr) {
            // if (current->elementIndex == );
            /* check intersection of current node*/
            intersect = true;
            for (unsigned short iDim = 0; iDim < SU2_BBOX_SIZE / 2; iDim++) {
                intersect = intersect && (testBBox[iDim] <= current->bBoxCoordinates[iDim + SU2_BBOX_SIZE / 2]);
                intersect = intersect && (testBBox[iDim + SU2_BBOX_SIZE / 2] >= current->bBoxCoordinates[iDim]);
            }

            if (intersect) {
                intersectingBBox.push_back(current->elementIndex);
                // cout << "intersection found with element: " << current->elementIndex << endl;
            }
            i = currHeirarchy % 4;

            /*branching based on minimum coordinate*/
            if (i < SU2_BBOX_SIZE / 2) {
                /*curr->left is always searched as the left has min coords lower than current which gives no info on intersection */
                searchQ.push(make_pair(current->left, currHeirarchy + 1));
                /*Test _max < Current _min*/
                if (testBBox[i + SU2_BBOX_SIZE / 2] < current->bBoxCoordinates[i]) {
                    current = nullptr;
                } else {
                    current = current->right;
                }
            }
            /*branching based on maximum coordinate*/
            else {
                /*curr->right is always searched as the right has max coords higher than current which gives no info on intersection */
                searchQ.push(make_pair(current->right, currHeirarchy + 1));
                /*Test _min > Current _max*/
                if (testBBox[i - SU2_BBOX_SIZE / 2] > current->bBoxCoordinates[i]) {
                    current = nullptr;
                } else {
                    current = current->left;
                }
            }
            // current = nullptr;
            currHeirarchy = currHeirarchy + 1;
        }
    }
    if (intersectingBBox.size() == 0) {
        // cout << "No intersecting element BBox found." << endl;
    }
    return intersectingBBox;
}

/* Level order output ADT */
void ADT::printLevelOrder() const {
    if (root == nullptr) {
        cerr << "Empty ADT" << endl;
        return;
    }
    queue<node_adt *> q;
    q.push(root);

    while (q.empty() == false) {
        int count = q.size();
        while (count > 0) {
            node_adt *node = q.front();
            cout << node->elementIndex << " ";
            q.pop();
            if (node->left != NULL)
                q.push(node->left);
            if (node->right != NULL)
                q.push(node->right);
            count--;
        }
        cout << endl;
    }
}

// Function to generate DOT syntax for the BST (helper function)
string ADT::generateDot(node_adt *node) const {
    if (node == nullptr) {
        return "";
    }
    string dot;
    string nodeName = to_string(node->elementIndex);
    dot += "\"" + nodeName + "\"";
    dot += "[label=\"" + nodeName + "\"]\n";

    if (node->left != nullptr) {
        dot += "\"" + nodeName + "\"" + " -> \"" + to_string(node->left->elementIndex) + "\"";
        dot += "[color=red]\n";
        dot += generateDot(node->left);
    } else if (node->right != nullptr) {
        dot += "\"" + nodeName + "_NULL\" ";
        dot += "[shape=point]\n";
        dot += "\"" + nodeName + "\"" + " -> \"" + nodeName + "_NULL\"";
        dot += "[color=red]\n";
    }
    if (node->right != nullptr) {
        dot += "\"" + nodeName + "\"" + " -> \"" + to_string(node->right->elementIndex) + "\"";
        dot += "[color=blue]\n";
        dot += generateDot(node->right);
    } else if (node->left != nullptr) {
        dot += "\"" + nodeName + "_NULL\" ";
        dot += "[shape=point]\n";
        dot += "\"" + nodeName + "\"" + " -> \"" + nodeName + "_NULL\"";
        dot += "[color=blue]\n";
    }
    return dot;
}

// Generate DOT syntax and write to a file
void ADT::writeDotToFile(const string &filename) const {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }
    file << "digraph BST {" << endl;
    file << generateDot(root);
    file << "}" << endl;
    file.close();
}