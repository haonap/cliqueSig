//
// Created by Hao Pan on 9/2/20.
//

#ifndef CLIQUESIG_CLASSES_H
#define CLIQUESIG_CLASSES_H


#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
using namespace std;

class node{
public:
    int name;
    int degree = 0;
    vector<int> neighbors;
    vector<int> adjacency; // used in CommunityPeel
    vector<int> F; // used in CommunityPeel
    vector<int> numCommonNeighbors; // used in CommunityPeel
    vector<vector<int>> commonNeighbors; // used in CommunityPeel

    //constructor
    node(int);
};

class graph {
public:

    string name;
    int n, m;
    vector<node> nodeList;
    int maxDeg;
    int maxDegNode;
    vector<int> kCore;

    // constructor
    graph(string);
    graph(){};

    // functions
    void GetKCore(int); // this function gets kCore of a single static graph without reordering vertices
    bool IsAdj(int, int); // this function checks if two nodes are adjacent
    vector<int> FindCommonNeighbors(int, int);
    vector<int> BFS(int);
};


#endif //CLIQUESIG_CLASSES_H
