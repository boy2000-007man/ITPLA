#include "mainwidget.h"
#include <QPainter>
#include "sdk.h"
using namespace sdk;
Points tmp;
Points res;
vector<double> ang;
double edge_length = 100;
const int attention = -1;
MainWidget::MainWidget(QWidget *parent) :
    QWidget(parent)
{
    tmp = test();
    pair<Points, vector<double> > t = place(tmp, edge_length);
    res = t.first;
    ang = t.second;
}

void MainWidget::paintEvent(QPaintEvent *) {
    QPainter painter(this);

    //painter.translate(10, 10);
    //painter.scale(1,1);
    //painter.setWindow(0, 0, 2, 2);
    //painter.setViewport(0, 0, this->width(), this->height());

    painter.setPen(Qt::black);
    for (int i = 0; i < tmp.size(); i++) {
        Point &p1 = tmp[i],
              &p2 = tmp[(i + 1) % tmp.size()];
        painter.drawLine(p1.x, p1.y, p2.x, p2.y);
    }
    painter.setPen(Qt::green);
    double r = edge_length / 2 / sqrt(3);
    for (int i = 0; i < res.size(); i++)
        painter.drawEllipse(res[i].x - r, res[i].y - r, 2 * r, 2 * r);
    painter.setPen(Qt::red);
    for (int i = 0; i < ang.size(); i++) {
        painter.setPen(i == attention || attention == -1 ? Qt::red : Qt::white);
        double a = to_rad(ang[i]);
        double ttt = to_rad(120);
        Point &p = res[i];
        painter.drawLine(p.x, p.y, p.x + edge_length * sin(a) / (2 * sqrt(3)), p.y + edge_length * cos(a) / (2 * sqrt(3)));
        painter.drawLine(p.x, p.y, p.x + edge_length * sin(a + ttt) / (2 * sqrt(3)), p.y + edge_length * cos(a + ttt) / (2 * sqrt(3)));
        painter.drawLine(p.x, p.y, p.x + edge_length * sin(a + 2 * ttt) / (2 * sqrt(3)), p.y + edge_length * cos(a + 2 * ttt) / (2 * sqrt(3)));
    }
    painter.setPen(Qt::blue);
    for (int i = 0; i < ang.size(); i++) {
        painter.setPen(i == attention || attention == -1 ? Qt::blue : Qt::white);
        double a = to_rad(ang[i] + 60);
        double ttt = to_rad(120);
        Point &p = res[i],
               p0 = make_pair(p.x + edge_length * sin(a) / sqrt(3), p.y + edge_length * cos(a) / sqrt(3)),
               p1 = make_pair(p.x + edge_length * sin(a + ttt) / sqrt(3), p.y + edge_length * cos(a + ttt) / sqrt(3)),
               p2 = make_pair(p.x + edge_length * sin(a + 2 * ttt) / sqrt(3), p.y + edge_length * cos(a + 2 * ttt) / sqrt(3));
        painter.drawLine(p0.x, p0.y, p1.x, p1.y);
        painter.drawLine(p1.x, p1.y, p2.x, p2.y);
        painter.drawLine(p2.x, p2.y, p0.x, p0.y);
    }
//    painter.setPen(Qt::yellow);
//    for (int i = 0; i < res.size(); i++) {
//        Point &p1 = res[i];
//        for (int j = 0; j < i; j++) {
//            Point &p2 = res[j];
//            double r_distance = sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2)) / edge_length;
//            if (0.9 < r_distance && r_distance < 1.1)
//                painter.drawLine(p1.x, p1.y, p2.x, p2.y);
//        }
//    }
    painter.setPen(Qt::white);
    for (int i = 0; i < res.size(); i++)
        painter.drawRect(res[i].x - 1, res[i].y - 1, 2, 2);
}
