#ifndef GRAPHICEDGE_H
#define GRAPHICEDGE_H

#include <QGraphicsObject>
#include"Headers/graphicnode.h"
#include<QLineEdit>

class GraphicEdge : public QGraphicsObject {
Q_OBJECT
public:
    GraphicEdge(GraphicNode* start, GraphicNode* end, int weight);

    virtual ~GraphicEdge(){}

    QRectF boundingRect() const override;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget) override;
    int getWeight(){return m_weight;}
    GraphicNode* getStart(){return m_start;}
    GraphicNode* getEnd(){return m_end;}
    QPointF getCenter();
    QLineEdit* getLineEdit(){return m_weightLineEdit;}

public slots:
    void editWeight(const QString &text);
signals:
    void weightEdited(GraphicEdge* e, int w);
    void needRedraw();

private:
    GraphicNode* m_start;
    GraphicNode* m_end;
    int m_weight;
    QLineEdit* m_weightLineEdit;
};


#endif // GRAPHICEDGE_H
