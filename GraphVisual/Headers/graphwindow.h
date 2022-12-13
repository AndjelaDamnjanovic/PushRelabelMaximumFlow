#ifndef GRAPHWINDOW_H
#define GRAPHWINDOW_H

#include <QMainWindow>
#include<QPointF>

#include"Headers/graphicedge.h"
#include"Headers/graph.h"

class Node;
class GraphicNode;
class QGraphicsScene;

QT_BEGIN_NAMESPACE
namespace Ui { class GraphWindow; }
QT_END_NAMESPACE

class GraphWindow : public QMainWindow
{
    Q_OBJECT

public:
    GraphWindow(QWidget *parent = nullptr);
    ~GraphWindow();

signals:
    void AddedNewNode(GraphicNode *);
    void AddedNewEdge(GraphicEdge *);
    void DeletedGraph();
    void NeedRedraw();


private slots:
    void AddNewEdge();
    void DeleteGraphFromTable();
    void ChangeMode(int index);
    void AddNode(Node* node);
    void AddEdge(Node* n1, Node* n2);
    void changeWeight(Node* n1, Node* n2, int weight);


    void graphDirected();
    void graphUndirected();

    void on_actionSaveAsPng_triggered();
    void on_actionSaveAsJpg_triggered();

    void on_pbUndirected_pressed();

    void on_pbUndirected_released();

    void on_pbDirected_pressed();

    void on_pbUndirected_clicked();

    void on_actionClose_triggered();

    void on_pbAddNode_clicked();

private:
    Ui::GraphWindow *ui;
    QGraphicsScene *m_GraphTable;
    Graph *m_graph;
    bool shouldPopUpUndir = false;
    bool shouldPopUpDir = true;

};
#endif // GRAPHWINDOW_H
