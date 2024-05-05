#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <iostream>
#include "Headers/edge.h"
#include "Headers/graph.h"
#include "Headers/node.h"

#include <QList>
#include <QPair>

class Algorithm {
 public:
  Algorithm() {}

  QList<Node*> BFS(Node* current);
  void DFS(Node* current, QHash<Node*, bool>& visited, QList<Node*>& steps);
  bool isConnected(Node& u, Node& v);
  bool isAllConnected(Graph& graph);
  int Dijkstra(Graph& graph,
               Node* start,
               Node* end,
               QList<Node*>& path,
               QList<Node*>& visit,
               QList<QPair<Node*, Node*>>& edges);
  Node* minDist(QHash<Node*, int> dist, QHash<Node*, bool> visited);
  std::map<Node*, Node*> MST(Graph& graph);
  QList<Edge*> getBridges(Graph& graph);
  void bridge(Graph& graph,
              Node* node,
              QHash<Node*, bool>& visited,
              QHash<Node*, int>& in,
              QHash<Node*, int>& low_link,
              QHash<Node*, Node*>& parent,
              int time,
              QList<Edge*>& result);
  QSet<Node*> getArticulationNodes(Graph& graph);
  void articulationNodes(Node* node,
                         QHash<Node*, bool>& visited,
                         QHash<Node*, int>& in,
                         QHash<Node*, int>& low_link,
                         QHash<Node*, Node*>& parent,
                         int time,
                         QSet<Node*>& result);
  QList<std::string> Hierholzer(Graph graph);
  bool hasEulerianCircuit(Graph& graph);
  QList<std::string> getEulerianCircuit(Graph graph);
  Node* canPush(Graph& graph, Node* node);
  void push(Graph& graph, Node* first, Node* second, QList<int>& flows);
  void relabel(Graph& graph, Node* node);
  void relabelCurrentArc(Graph& graph, Node* node, QList<Node*> activeNodes);
  void relabelFIFO(Graph& graph, Node* node, QList<Node*> activeNodes);
  int genericPushRelabel(Graph& graph,
                         Node* source,
                         Node* sink,
                         QList<QPair<Node*, Node*>>& edges,
                         QList<int>& flows);
  int currentArcPushRelabel(Graph& graph,
                            Node* source,
                            Node* sink,
                            QList<QPair<Node*, Node*>>& edges,
                            QList<int>& flows);
  int FIFOPushRelabel(Graph& graph,
                            Node* source,
                            Node* sink,
                            QList<QPair<Node*, Node*>>& edges,
                            QList<int>& flows);
};

#endif  // ALGORITHM_H
