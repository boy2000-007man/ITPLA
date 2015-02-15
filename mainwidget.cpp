#include "mainwidget.h"
#include <QPainter>
#include "sdk.h"
using namespace sdk;
Points tmp;
Points res;
double edge_length = 80;
MainWidget::MainWidget(QWidget *parent) :
    QWidget(parent)
{
    tmp = test();
    res = place(tmp, edge_length);
}

void MainWidget::paintEvent(QPaintEvent *) {
    QPainter painter(this);

    painter.translate(10, 10);
    painter.scale(1,1);

    painter.setPen(Qt::black);
    for (int i = 0; i < tmp.size(); i++) {
        Point &p1 = tmp[i],
              &p2 = tmp[(i + 1) % tmp.size()];
        painter.drawLine(p1.x, p1.y, p2.x, p2.y);
    }
    painter.setPen(Qt::green);
    for (int i = 0; i < res.size(); i++)
        painter.drawEllipse(res[i].x - edge_length, res[i].y - edge_length, 2 * edge_length, 2 * edge_length);
    painter.setPen(Qt::red);
    for (int i = 0; i < res.size(); i++)
        painter.drawRect(res[i].x - 1, res[i].y - 1, 2, 2);
}
