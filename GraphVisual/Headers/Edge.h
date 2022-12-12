#ifndef EDGE_H
#define EDGE_H

#include<unordered_map>
#include<map>
#include "Headers/Node.h"

typedef std::pair<Node*, Node*> pair;
typedef std::unordered_map<pair, int> edges;

class Edge{

public:
    Edge(pair&, int&);

    explicit Edge(edges::iterator&);
    Node* first() const;
    Node* second() const;
    int weight() const;

    friend class Graph;

private:
    const int &m_weight{};
    const pair &m_node_pair{};
};

#endif // EDGE_H