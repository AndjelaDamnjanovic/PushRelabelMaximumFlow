#include <utility>

#include "Headers/edge.h"
#include "Headers/node.h"

Edge::Edge(std::pair<Node*, Node*> nodePair, int weight, int flow, bool belongsToResidualGraph):
    m_weight(weight), m_nodePair(std::move(nodePair)), m_flow(0), m_belongsToResidualGraph(false){}

Edge::Edge(const Edge *other) {
    m_weight = other->m_weight;
    m_flow = other->m_flow;
    m_belongsToResidualGraph = other->m_belongsToResidualGraph;
    auto *n1 = new Node(other->m_nodePair.first);
    auto *n2 = new Node(other->m_nodePair.second);
    m_nodePair = std::pair<Node*, Node*>(n1, n2);
}

Edge::Edge(edges::iterator& iter): m_weight(iter->second), m_nodePair(iter->first){}

Edge::~Edge() {
    delete m_nodePair.first;
    delete m_nodePair.second;
}

Node* Edge::first()const{
    return m_nodePair.first;
}

Node* Edge::second() const{
    return m_nodePair.second;
}

int Edge::weight() const{
    return m_weight;
}

int Edge::flow() const
{
    return m_flow;
}

bool Edge::belongsToResidualGraph() const
{
    return m_belongsToResidualGraph;
}

QVariant Edge::toVariant()
{
    QVariantMap map;
    QString nodeName=QString::fromStdString(this->m_nodePair.first->name());
    map.insert("node1", nodeName);

    nodeName=QString::fromStdString(this->m_nodePair.second->name());
    map.insert("node2", nodeName);

    int nodeWeight= this->weight();
    map.insert("weight", nodeWeight);

    return map;
}

std::ostream &operator<<(std::ostream &os, const Edge &e){
    os << e.first()->name() << " " << e.second()->name() << e.weight() << std::endl;
    return os;
}
