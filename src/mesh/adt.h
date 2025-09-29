#ifndef ADT_H
#define ADT_H

#include <array>
#include <vector>
#include <string>
#include "util/constants.h"

class node_adt {
   public:
    size_t elementIndex;
    std::array<su2double, SU2_BBOX_SIZE> bBoxCoordinates;
    node_adt *left;
    node_adt *right;

    node_adt(const size_t elementIdx, const std::array<su2double, SU2_BBOX_SIZE> &bBoxCoords);
    virtual ~node_adt() = default;
};

class ADT {
   public:
    node_adt *root;
    size_t treeHeirarchy;
    
    ADT();
    virtual ~ADT() = default;

    void insertNode(node_adt *node);
    std::vector<size_t> searchADT(const std::array<su2double, SU2_BBOX_SIZE> &testBBox) const;
    void printLevelOrder() const;
    void writeDotToFile(const std::string &filename) const;

private:
    std::string generateDot(node_adt *node) const;
};

#endif // ADT_H