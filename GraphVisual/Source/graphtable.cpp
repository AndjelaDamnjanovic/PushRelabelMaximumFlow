#include "Headers/graphtable.h"
#include "Headers/graphicnode.h"
#include "Headers/graphicedge.h"

#include "Headers/popup.h"
#include <QMessageBox>
#include <iostream>


GraphTable::GraphTable(QObject *parent)
    : QGraphicsScene(parent) {

}

void GraphTable::AddNewNodeOnTable(GraphicNode *node) {
    m_Nodes.append(node);

    connect(node, &GraphicNode::needRedraw, this, &GraphTable::Redraw);

    PlaceNodeOnTable(node);
}

void GraphTable::AddNewEdgeOnTable(GraphicEdge *edge) {
    m_Edges.append(edge);
    edge->setFlag(QGraphicsItem::ItemIgnoresTransformations, false);
    connect(edge, &GraphicEdge::weightEdited, this, &GraphTable::editWeight);
    connect(edge, &GraphicEdge::needRedraw, this, &GraphTable::Redraw);
    connect(edge, &GraphicEdge::needWarning, this, &GraphTable::Warning);

    addItem(edge);
    addWidget(edge->getLineEdit());
}

void GraphTable::DeleteGraphFromTable() {
    for(auto edge : m_Edges) {
        removeItem(edge);
        delete edge;
    }
    m_Edges.clear();

    for(auto node : m_Nodes) {
        removeItem(node);
        delete node;
    }
    m_Nodes.clear();

    update();
}

void GraphTable::Redraw(){
    this->update();
}

void GraphTable::PlaceNodeOnTable(GraphicNode *node) {
    const int tableWidth = static_cast<int>(width());
    auto newNodeIndex = m_Nodes.size() - 1;

    const auto xPos = (node->m_width * newNodeIndex) % tableWidth;
    const auto yPos = node->m_height * ((node->m_width * newNodeIndex) / tableWidth);
    node->setPos(xPos, yPos);

    node->setFlag(QGraphicsItem::ItemIgnoresTransformations, false);
}


QVector<GraphicNode *> GraphTable::getNodes(){
    return m_Nodes;
}
QVector<GraphicEdge *> GraphTable::getEdges(){
    return m_Edges;
}

void GraphTable::editWeight(GraphicEdge* edge, int w){
    emit edgeWeightChanged(edge->getStart()->getNode(), edge->getEnd()->getNode(), w);
}

void GraphTable::Warning (QString s){
    emit needWarning(s);
}

bool GraphTable::hasGraphicEdge(GraphicNode *u, GraphicNode *v) {
    auto it = m_Edges.begin();
    for(;it != m_Edges.end(); ++it){
        if ((*it)->getStart()->getNode() == u->getNode() && (*it)->getEnd()->getNode() == v->getNode()){
            return true;
        }
    }
    return false;
}

void GraphTable::mousePressEvent ( QGraphicsSceneMouseEvent * event ){
//    TODO
//    if(event->pos() < 0) {
//        setPos(0, event->pos().y());
//    } else if (x() + boundingRect().right() > scene()->width()) {
//        setPos(scene()->width() - boundingRect().width(), y());
//    }

//    if (y() < 0) {
//        setPos(x(), 0);
//    } else if (y()+ boundingRect().bottom() > scene()->height()) {
//        setPos(x(), scene()->height() - boundingRect().height());
//    }



    if(m_drawingMode && itemAt(event->scenePos(), QTransform()) == NULL){
        if(m_hasTmp){
            setHasTmp(false);
        }
        else{ //add new node
            Popup* p = new Popup();

            if(p->exec() == QDialog::Accepted){

                QString nodeName = p->getNodeName();

                bool alreadyExists = false;
                for(GraphicNode* n : m_Nodes){
                    if(n->getNode()->name() == nodeName.toStdString()){
                        alreadyExists = true;
                        emit needWarning("Node with that name already exists");
                        break;
                    }
                }


                if(!alreadyExists && !nodeName.length()==0){

                    Node* node = new Node(nodeName.toStdString());
                    GraphicNode* graphicNode = new GraphicNode(node);

                    AddNewNodeOnTable(graphicNode);
                    graphicNode->setPos(event->scenePos() - QPointF(graphicNode->m_width / 2, graphicNode->m_height / 2));
                    addItem(graphicNode);

                    emit addedNewNode(node);
                }
            }
        }
    }
    else if(m_drawingMode && (itemAt(event->scenePos(), QTransform())->type() == 1)){ //click on node
        if(!m_hasTmp){
            this->m_tmp = dynamic_cast<GraphicNode*>(itemAt(event->scenePos(), QTransform()));
            setHasTmp(true);
        }
        else{
            GraphicNode* node = dynamic_cast<GraphicNode*>(itemAt(event->scenePos(), QTransform()));

            // open window to insert node name
            if(!hasGraphicEdge(m_tmp, node)) {
                GraphicEdge* edge = new GraphicEdge(m_tmp, node, 1);
                AddNewEdgeOnTable(edge);


                emit addedNewEdge(m_tmp->getNode(), node->getNode());
            }

            setHasTmp(false);
        }
        this->update();
    }
    else if(m_deleteMode && itemAt(event->scenePos(), QTransform()) == NULL){
        setHasTmp(false);
        QGraphicsScene::mousePressEvent(event);
    }
    else if(m_deleteMode && (itemAt(event->scenePos(), QTransform())->type() == 1)){ //delete node
        const auto graphicNode = dynamic_cast<GraphicNode*>(itemAt(event->scenePos(), QTransform()));
        removeItem(graphicNode);
        m_Nodes.removeOne(graphicNode);

        Node* node = graphicNode->getNode();

        // delete edges around the node
        auto i = m_Edges.begin();
        while (i != m_Edges.end()) {
            if((*i)->getStart()->getNode()==node || (*i)->getEnd()->getNode()==node){
                removeItem(*i);
                removeItem((*i)->getLineEdit()->graphicsProxyWidget());
                i = m_Edges.erase(i);
                Redraw();
            }
            else
                ++i;
        }

        emit deletedNode(node);
    }
    else if(m_deleteMode && (itemAt(event->scenePos(), QTransform())->type() == 2)){ //delete edge
        const auto graphicEdge = dynamic_cast<GraphicEdge*>(itemAt(event->scenePos(), QTransform()));
        removeItem(graphicEdge);
        removeItem(graphicEdge->getLineEdit()->graphicsProxyWidget());
        m_Edges.removeOne(graphicEdge);

        Redraw();

        emit deletedEdge(graphicEdge->getStart()->getNode(), graphicEdge->getEnd()->getNode());
    }
    else{
        setHasTmp(false);
        QGraphicsScene::mousePressEvent(event);
    }
}

//TODO: drawing a line when clicked on a node

//void GraphTable::mouseMoveEvent (QGraphicsSceneMouseEvent * event ){


//    if(m_drawingMode && m_hasTmp){
//        for (auto veza : _veze){
//            this->removeItem(veza);
//            delete veza;
//        }
//        _veze.clear();


//        QLineF linija(m_tmp->TopCenterPosition(), event->scenePos());
//        auto veza = new QGraphicsLineItem(linija);
//        _veze.append(veza);
//        this->addItem(veza);

//        this->update();
//    }
//    QGraphicsScene::mouseMoveEvent(event);

//}

void GraphTable::setDrawingMode(bool x) {
    m_drawingMode = x;
}
void GraphTable::setHasTmp(bool x) {
    m_hasTmp = x;
}
void GraphTable::setDeleteMode(bool x) {
    m_deleteMode = x;
}

