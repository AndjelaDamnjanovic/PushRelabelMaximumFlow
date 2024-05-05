#ifndef EDGE_H
#define EDGE_H

#include <QVariant>

#include <map>
#include "Headers/node.h"

typedef std::map<std::pair<Node*, Node*>, int> edges;

class Edge {
 public:
  Edge(std::pair<Node*, Node*>, int weight, int flow, bool belongsToResidualGraph);
  Edge(const Edge* other);
  ~Edge();

  explicit Edge(edges::iterator& iter);
  Node* first() const;
  Node* second() const;
  int weight() const;
  int flow() const;
  bool belongsToResidualGraph() const;
  QVariant toVariant();
  friend std::ostream& operator<<(std::ostream& os, const Edge e);

  friend class Graph;

  void setWeight(int n) { m_weight = n; }
  void setBelongsToResidual(bool belongs) { m_belongsToResidualGraph = belongs; }
  void setFlow(int n) { m_flow = n; }

 private:
  int m_weight;
  int m_flow;
  bool m_belongsToResidualGraph;
  std::pair<Node*, Node*> m_nodePair;
};

#endif  // EDGE_H
