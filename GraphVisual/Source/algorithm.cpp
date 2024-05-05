#include"Headers/algorithm.h"

#include<queue>
#include<QHash>
#include<QQueue>
#include<QStack>
#include<algorithm>
#include<limits>
#include<cassert>
#include<limits>
#include<string>

QList<Node*> Algorithm::BFS (Node* current){

    QList<Node*> steps;
    if(current == nullptr)
        return steps;

    QQueue<Node*> q;
    QHash<Node*, bool> visited;

    q.push_back(current);
    visited[current] = true;

    while(!q.empty()){
        current = q.front();
        q.pop_front();
        steps.push_back(current);

        for(auto neighb : current->neighbours())
            if(visited.find(neighb)==visited.end()){
                visited[neighb] = true;
                q.push_back(neighb);
            }
    }

    return steps;
}


void Algorithm::DFS (Node* current, QHash<Node*, bool> &visited, QList<Node*> &steps){

    visited[current] = true;
    steps.push_back(current);

    for(Node* neighb : current->neighbours())
        if(visited.find(neighb)==visited.end())
            DFS(neighb, visited, steps);
}


bool Algorithm::isConnected (Node &u, Node &v){

    QHash<Node*, bool> visited;
    QList<Node*> steps;

    DFS(&u, visited, steps);

    return visited[&v];
}


bool Algorithm::isAllConnected (Graph &graph){

    for (auto node : graph.nodeSet()) {

        QHash<Node*, bool> visited;
        QList<Node*> steps;

        DFS(node, visited, steps);

        for(auto v : graph.nodeSet())
            if(visited.find(v)==visited.end())
                return false;
    }

    return true;
}


int Algorithm::Dijkstra (Graph &graph, Node* start, Node* end, QList<Node*> &path, QList<Node*> &visit, QList<QPair<Node*, Node*>> &edges){

    QHash<Node*, bool> visited;
    QHash<Node*, int> dist;
    QHash<Node*, Node*> parent;

    dist[start] = 0;

    Node* current = start;
    parent[current] = nullptr;

    while(current != nullptr && !(current == end)){
        visited[current] = true;
        visit.append(current);

        for(auto neighb : current->neighbours()){
            if(visited.find(neighb)==visited.end() && (dist.find(neighb)==dist.end() || dist[current] + graph.weight(current, neighb) < dist[neighb])) {
                dist[neighb] = dist[current] + graph.weight(current, neighb);
                parent[neighb] = current;
                edges.append(qMakePair(current, neighb));
            }
        }

        current = minDist(dist, visited);
    }
    visit.append(end);

    current = end;
    while(current != nullptr){
        path.append(current);
        current = parent[current];
    }
    std::reverse(path.begin(), path.end());

    if(dist.find(end)==dist.end())
        return -1;

    return dist[end];
}


Node* Algorithm::minDist(QHash<Node*, int> dist, QHash<Node*, bool> visited){
    int min = std::numeric_limits<int>::max();
    Node* minNode = nullptr;
    for (auto iter = dist.constBegin(); iter != dist.constEnd(); ++iter) {
        if(visited.find(iter.key()) == visited.end() && iter.value()<min){
            min = iter.value();
            minNode = iter.key();
        }
    }
    return minNode;
}


std::map<Node*, Node*> Algorithm::MST (Graph &graph){

    std::map<Node*, Node*> parent;
    if(graph.edgeSet().size()==0){
        if(graph.nodeSet().size()==0)
            return parent;
        else{
            parent[graph.randomNode()] = nullptr;
            return parent;
        }
    }


    std::map<Node*, bool> visited;
    std::map<Node*, int> minEdge;

    std::priority_queue<std::pair<int, Node*>, std::vector<std::pair<int, Node*>>, std::greater<>> minDist;

    int min = std::numeric_limits<int>::max();
    Node* begin;
    Node* end;

    for(Edge *e : graph.edgeSet()){
        if(e->weight() < min){
            min = e->weight();
            begin = e->first();
            end = e->second();
        }
    }

    minDist.push(std::make_pair(0, begin));

    for(Node* n : graph.nodeSet()){
        if(!(n == begin)){
            minDist.push(std::make_pair(std::numeric_limits<int>::max(), n));
            minEdge[n] = std::numeric_limits<int>::max();
        }
    }

    minEdge[begin] = 0;
    parent[end] = begin;


    while(!minDist.empty()) {
        auto node = minDist.top();
        minDist.pop();

        Node* current = node.second;

        if(visited.find(current) != visited.end())
            continue;


        visited[current] = true;
        minEdge[current] = node.first;

        for(auto neighb : current->neighbours()){
            if(visited.find(neighb)==visited.end() && minEdge[neighb] > graph.weight(current, neighb)){
                parent[neighb] = current;
                minDist.push(std::make_pair(graph.weight(current, neighb), neighb));
                minEdge[neighb] = graph.weight(current, neighb);
            }
        }
    }

    parent[begin] = nullptr;

    return parent;
}




QList<Edge*> Algorithm::getBridges (Graph &graph){

    QHash<Node*, bool> visited;
    QHash<Node*, int> in;
    QHash<Node*, int> low_link;
    QHash<Node*, Node*> parent;

    QList<Edge*> result;
    int time = 0;

    for(auto node : graph.nodeSet())
        if(visited.find(node)==visited.end()){
            bridge(graph, node, visited, in, low_link, parent, time, result);
        }

    return result;
}


void Algorithm::bridge (Graph &graph, Node* node, QHash<Node*, bool> &visited,
                    QHash<Node*, int> &in, QHash<Node*, int> &low_link,
                        QHash<Node*, Node*> &parent, int time, QList<Edge*> &result){

    visited[node] = true;
    in[node] = low_link[node] = ++time;

    for(auto neighb : node->neighbours()){
        if(visited.find(neighb)==visited.end()){
            parent[neighb] = node;

            bridge(graph, neighb, visited, in, low_link, parent, time, result);

            low_link[node] = std::min(low_link[node], low_link[neighb]);

            if(low_link[neighb] > in[node])
                result.push_back(graph.getEdge(node, neighb));
        }
        else if(!(neighb == parent[node]))
            low_link[node] = std::min(low_link[node], in[neighb]);
    }


}

QSet<Node*> Algorithm::getArticulationNodes(Graph &graph){

    QHash<Node*, bool> visited;
    QHash<Node*, int> in;
    QHash<Node*, int> low_link;
    QHash<Node*, Node*> parent;

    QSet<Node*> result;
    int time = 0;

    for(auto node : graph.nodeSet()){
        if(visited.find(node)==visited.end()){
            articulationNodes(node, visited, in, low_link, parent, time, result);
        }
    }

    return result;
}


void Algorithm::articulationNodes (Node* node, QHash<Node*, bool> &visited,
                                 QHash<Node*, int> &in, QHash<Node*, int> &low_link,
                                   QHash<Node*, Node*> &parent, int time, QSet<Node*> &result){

    int children = 0;
    visited[node] = true;
    in[node] = low_link[node] = ++time;

    for(auto neighb : node->neighbours()){
        if(visited.find(neighb)==visited.end()){
            children++;
            parent[neighb] = node;

            articulationNodes(neighb, visited, in, low_link, parent, time, result);

            low_link[node] = std::min(low_link[node], low_link[neighb]);

            if(!(parent.find(node)==parent.end()) && low_link[neighb] >= in[node])
                result.insert(node);
        }
        else if(!(neighb == parent[node]))
            low_link[node] = std::min(low_link[node], in[neighb]);
    }

    if(parent.find(node)==parent.end() && children > 1)
        result.insert(node);
}


QList<std::string> Algorithm::Hierholzer (Graph graph){

    QList<std::string> result;

    QStack<Node*> currPath;
    std::vector<Node*> cycle;

    Node* nextNode;

    Node* currNode = graph.randomNode();
    currPath.push(currNode);

    while (!currPath.empty()) {
        if (currNode->neighbours().size()) {
            currPath.push(currNode);

            nextNode = currNode->neighbours().back();

            graph.removeEdge(currNode, nextNode);

            currNode = nextNode;
        } else {
            cycle.push_back(currNode);
            currNode = currPath.top();
            currPath.pop();
        }
    }
    for (int i = (int) cycle.size() - 1; i >= 0; --i){
        result.push_back(cycle[i]->name());
    }

    return result;
}




bool Algorithm::hasEulerianCircuit (Graph &graph){

    if(!isAllConnected(graph)){
        return false;
    }

    for(auto node : graph.nodeSet()){
        if(graph.isDirected() && node->inDeg() != node->outDeg())
            return false;

        if(graph.isUndirected() && node->deg() % 2 != 0)
            return false;
    }

    return true;
}



QList<std::string> Algorithm::getEulerianCircuit (Graph graph){

    if(!hasEulerianCircuit(graph)){
        return QList<std::string>();
    }

    QList<std::string> result = Hierholzer (graph);

    return result;
}

Node *Algorithm::canPush(Graph &graph, Node *node)
{
    QList<Node*> neighbours = node->neighbours();

    for(const auto& neighbour : neighbours){
        Edge* edge = graph.getEdge(node, neighbour);

        // if the flow is equal to the weight of an edge, additional flow cannot be pushed
        if(edge->flow() == edge->weight())
            continue;

        // the push can only be performed if the height of the neighbouring node is smaller than the height of initial node
        if(neighbour->height() < node->height() && !edge->belongsToResidualGraph()){
            return neighbour;
        }

    }
    return nullptr;
}

void Algorithm::push(Graph& graph, Node *first, Node *second, QList<int>& flows)
{
    int height1 = first->height();
    int height2 = second->height();
    int excessFlow = first->excessFlow();

    assert((( height1 == height2 + 1) && excessFlow > 0));
    Edge* edge = graph.getEdge(first, second);
    Edge* residualEdge = graph.getEdge(second, first);

    if(residualEdge == nullptr){
        graph.addEdge(second, first);
        residualEdge = graph.getEdge(second, first);
        residualEdge->setWeight(edge->weight());
        residualEdge->setBelongsToResidual(true);
    }
    int delta = 0;
    if(!edge->belongsToResidualGraph()){
        delta = std::min(first->excessFlow(), edge->weight() - edge->flow());
    }else{
        delta = std::min(first->excessFlow(), edge->flow());
    }

    flows.append(delta);

    int edgeFlow = edge->flow();
    int residualEdgeFlow = residualEdge->flow();
    int excessFlowFirst = first->excessFlow();
    int excessFlowSecond = second->excessFlow();

    edge->setFlow(edgeFlow + delta);
    residualEdge->setFlow(residualEdgeFlow - delta);
    first->setExcessFlow(excessFlowFirst - delta);
    second->setExcessFlow(excessFlowSecond + delta);
}

void Algorithm::relabel(Graph &graph, Node *node)
{
    QList<Node*> neighbours = node->neighbours();
    int minHeight = INT_MAX;
    for(const auto& neighbour : neighbours){
        Edge* edge = graph.getEdge(node, neighbour);
        if(neighbour->height() < minHeight && (edge->weight() - edge->flow() > 0)){
            minHeight = neighbour->height();
        }
    }
    //elimination of excess flow
    if(minHeight + 1 == node->height()){
        node->setExcessFlow(0);
    }

    node->setHeight(minHeight + 1);
}

void Algorithm::relabelCurrentArc(Graph &graph, Node *node, QList<Node*> activeNodes)
{
    QList<Node*> neighbours = node->neighbours();
    int minHeight = INT_MAX;
    for(const auto& neighbour : neighbours){
        Edge* edge = graph.getEdge(node, neighbour);
        if(neighbour->height() < minHeight && (edge->weight() - edge->flow() > 0)){
            minHeight = neighbour->height();
        }
    }
    //elimination of excess flow
    if(minHeight + 1 == node->height()){
        node->setExcessFlow(0);
    }

    node->setHeight(minHeight + 1);
    activeNodes.push_front(node);
}

void Algorithm::relabelFIFO(Graph &graph, Node *node, QList<Node *> activeNodes)
{
    QList<Node*> neighbours = node->neighbours();
    int minHeight = INT_MAX;
    for(const auto& neighbour : neighbours){
        Edge* edge = graph.getEdge(node, neighbour);
        if(neighbour->height() < minHeight && (edge->weight() - edge->flow() > 0)){
            minHeight = neighbour->height();
        }
    }
    //elimination of excess flow
    if(minHeight + 1 == node->height()){
        node->setExcessFlow(0);
    }

    node->setHeight(minHeight + 1);
    activeNodes.removeFirst();
}

int Algorithm::genericPushRelabel(Graph &graph, Node *source, Node *sink, QList<QPair<Node *, Node *> > &edges, QList<int>& flows)
{
    QList<Node*> activeNodes;
    QList<Node*> sourceNeighbours = source->neighbours();

    for(const auto& neighbour : sourceNeighbours){
        activeNodes.append(neighbour);
        edges.append(qMakePair(source, neighbour));
        Edge* edge = graph.getEdge(source, neighbour);
        flows.append(edge->weight());
    }

    while(!activeNodes.isEmpty()){
        Node* activeNode = activeNodes.first();

        while(activeNode->excessFlow() > 0){
            Node* destination = canPush(graph, activeNode);

            if(destination != nullptr){
                edges.append(qMakePair(activeNode, destination));
                push(graph, activeNode, destination, flows);
                if(destination != source && destination != sink){
                    activeNodes.append(destination);
                }
            }else{
                relabel(graph, activeNode);
            }
        }
        activeNodes.removeFirst();
    }
    return sink->excessFlow();
}

int Algorithm::currentArcPushRelabel(Graph &graph, Node *source, Node *sink, QList<QPair<Node *, Node *> > &edges, QList<int>& flows)
{
    QList<Node*> activeNodes;
    QList<Node*> sourceNeighbours = source->neighbours();

    for(const auto& neighbour : sourceNeighbours){
        activeNodes.append(neighbour);
        edges.append(qMakePair(source, neighbour));
        Edge* edge = graph.getEdge(source, neighbour);
        flows.append(edge->weight());
    }

    while(!activeNodes.isEmpty()){
        Node* activeNode = activeNodes.first();

        while(activeNode->excessFlow() > 0){
            Node* destination = canPush(graph, activeNode);

            if(destination != nullptr){
                edges.append(qMakePair(activeNode, destination));
                push(graph, activeNode, destination, flows);
                if(destination != source && destination != sink){
                    activeNodes.append(destination);
                }
            }else{
                relabelCurrentArc(graph, activeNode, activeNodes);
            }
        }
        activeNodes.removeFirst();
    }

    return sink->excessFlow();
}

int Algorithm::FIFOPushRelabel(Graph &graph, Node *source, Node *sink, QList<QPair<Node *, Node *> > &edges, QList<int>& flows)
{
    QList<Node*> activeNodes;
    QList<Node*> sourceNeighbours = source->neighbours();

    for(const auto& neighbour : sourceNeighbours){
        activeNodes.append(neighbour);
        edges.append(qMakePair(source, neighbour));
        Edge* edge = graph.getEdge(source, neighbour);
        flows.append(edge->weight());
    }

    while(!activeNodes.isEmpty()){
        Node* activeNode = activeNodes.first();

        while(activeNode->excessFlow() > 0){
            Node* destination = canPush(graph, activeNode);

            if(destination != nullptr){
                edges.append(qMakePair(activeNode, destination));
                push(graph, activeNode, destination, flows);
                if(destination != source && destination != sink){
                    activeNodes.append(destination);
                }
            }else{
                relabelFIFO(graph, activeNode, activeNodes);
            }
        }
        activeNodes.removeFirst();
    }

    return sink->excessFlow();
}
