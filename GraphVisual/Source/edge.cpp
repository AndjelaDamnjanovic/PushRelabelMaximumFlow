#include "Headers/edge.h"
#include "Headers/node.h"

Edge::Edge(std::pair<Node*, Node*> nodePair, int weight):m_weight(weight), m_nodePair(nodePair){}

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

void Edge::fromVariant(const QVariant &variant)
{
    const auto map = variant.toMap();
    QString node1 = map.value("node1").toString();
    QString node2 = map.value("node2").toString();
    int weight=map.value("weight").toInt();


}



std::ostream &operator<<(std::ostream &os, const Edge &e){
    os << e.first()->name() << " " << e.second()->name() << e.weight() << std::endl;
    return os;
}
