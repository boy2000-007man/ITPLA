#include "mainwidget.h"
#include <QPainter>
#include <QTimer>
#include "ITPLA.h"
using namespace ITPLA;
Points polygon, normalized_polygon;
double edge_length;
vector<b2Body *> points;
const int attention = -1;
const int FRAMES_PER_SEC = 240;
MainWidget::MainWidget(QWidget *parent) :
    QWidget(parent)
{
    edge_length = 100;
    polygon = read("/mnt/Zero_Data/zero/workspace/placement/test.csv");//test;
    for (Point p: polygon)
        cout << p.x << "|" << p.y << endl;

    edge_length = 158.88;
    polygon = read("/mnt/Zero_Data/zero/workspace/placement/tr1_1.csv");//tr1_1;
    for (int i = 0; i < polygon.size()/2; i++)
        swap(polygon[i], polygon[polygon.size()-1-i]);

    edge_length = 63.17796;//110.204;
    polygon = read("/mnt/Zero_Data/zero/workspace/placement/tr1_2.csv");//tr1_2;
    for (int i = 1; i < polygon.size(); i++)
        polygon[i] = {polygon[i-1].x+polygon[i].x, polygon[i-1].y+polygon[i].y};
    for (int i = 0; i < polygon.size()/2; i++)
        swap(polygon[i], polygon[polygon.size()-1-i]);

    edge_length = 85.86865;//110.204;
    polygon = read("/mnt/Zero_Data/zero/workspace/placement/tr1_4.csv");//tr1_4;
    for (int i = 1; i < polygon.size(); i++)
        polygon[i] = {polygon[i-1].x+polygon[i].x, polygon[i-1].y+polygon[i].y};
    for (int i = 0; i < polygon.size()/2; i++)
        swap(polygon[i], polygon[polygon.size()-1-i]);

    edge_length = 68.983385775;//110.204;
    polygon = read("/mnt/Zero_Data/zero/workspace/placement/tr1_5.csv");//tr1_5;
    for (int i = 1; i < polygon.size(); i++)
        polygon[i] = {polygon[i-1].x+polygon[i].x, polygon[i-1].y+polygon[i].y};
    for (int i = 0; i < polygon.size()/2; i++)
        swap(polygon[i], polygon[polygon.size()-1-i]);

    normalized_polygon = normalize_polygon(polygon, edge_length / sqrt(3));
    points = place(polygon, edge_length);
    QTimer *timer = new QTimer();
    timer->start(1000.0 / FRAMES_PER_SEC);
    connect(timer, SIGNAL(timeout()), this, SLOT(repaint()));
//    pair<Points, vector<double> > t = place(tmp, edge_length);
//    res = t.first;
//    ang = t.second;
}

void paint(QPaintDevice *qpd) {
    QPainter painter(qpd);
    painter.setBackgroundMode(Qt::OpaqueMode);
    painter.setBackground(QBrush(Qt::white));
    painter.eraseRect(painter.window());

    pair<Points, vector<double> > ttt;
    Vectors vel;
    for (int i = 0; i < points.size(); i++) {
        ttt.first.push_back(normalize_point(points[i]->GetPosition(), sqrt(3) / edge_length));
//        printf("%lf,%lf\n",points[i]->GetAngle(),points[i]->GetPosition().y);
        ttt.second.push_back(to_deg(points[i]->GetAngle()));
        vel.push_back(normalize_point(points[i]->GetLinearVelocity(), sqrt(3) / edge_length));
    }
    painter.setPen(Qt::black);
    painter.drawText(
                0,
                320,
                QString().sprintf("Point:%d,Frame:%6d,Time:%8.3lfs,K:%.6lf", points.size(), frame, 1.0*frame/FRAMES_PER_SEC, ITPLA::K)
    );
    Points res = ttt.first;
    vector<double> ang = ttt.second;
    painter.setPen(Qt::black);
    for (int i = 0; i < polygon.size(); i++) {
        Point &p1 = polygon[i],
              &p2 = polygon[(i + 1) % polygon.size()];
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
        painter.drawLine(p.x, p.y, p.x + r * sin(a), p.y + r * cos(a));
        painter.drawLine(p.x, p.y, p.x + r * sin(a + ttt), p.y + r * cos(a + ttt));
        painter.drawLine(p.x, p.y, p.x + r * sin(a + 2 * ttt), p.y + r * cos(a + 2 * ttt));
    }
    painter.setPen(Qt::black);
    for (int i = 0; i < points.size(); i++) {
        const Point &p = res[i];
        const Vector &v = vel[i];
        painter.drawLine(p.x, p.y, p.x + v.x, p.y + v.y);
    }
    painter.setPen(Qt::blue);
    for (int i = 0; i < ang.size(); i++) {
        painter.setPen(i == attention || attention == -1 ? Qt::blue : Qt::white);
        double a = to_rad(ang[i] + 60);
        double ttt = to_rad(120);
        Point &p = res[i],
               p0 = b2Vec2(p.x + 2 * r * sin(a), p.y + 2 * r * cos(a)),
               p1 = b2Vec2(p.x + 2 * r * sin(a + ttt), p.y + 2 * r * cos(a + ttt)),
               p2 = b2Vec2(p.x + 2 * r * sin(a + 2 * ttt), p.y + 2 * r * cos(a + 2 * ttt));
        painter.drawLine(p0.x, p0.y, p1.x, p1.y);
        painter.drawLine(p1.x, p1.y, p2.x, p2.y);
        painter.drawLine(p2.x, p2.y, p0.x, p0.y);
    }
    painter.setPen(Qt::white);
    for (int i = 0; i < res.size(); i++)
        painter.drawRect(res[i].x - 1, res[i].y - 1, 2, 2);
}

void MainWidget::paintEvent(QPaintEvent *) {
    if (frame % 1000 == 0) {
        QPixmap qi(this->width(), this->height());
        paint(&qi);
        qi.save(QString().sprintf("/home/zero/%d.bmp", frame));
    }
    QPainter painter(this);

    pair<Points, vector<double> > ttt;
    Vectors vel;
    for (int i = 0; i < points.size(); i++) {
        ttt.first.push_back(normalize_point(points[i]->GetPosition(), sqrt(3) / edge_length));
//        printf("%lf,%lf\n",points[i]->GetAngle(),points[i]->GetPosition().y);
        ttt.second.push_back(to_deg(points[i]->GetAngle()));
        vel.push_back(normalize_point(points[i]->GetLinearVelocity(), sqrt(3) / edge_length));
    }
    calc_next_step(normalized_polygon, points);
    painter.setPen(Qt::black);
    painter.drawText(
                0,
                320,
                QString().sprintf("Point:%d,Frame:%6d,Time:%8.3lfs,K:%.6lf", points.size(), frame, 1.0*frame/FRAMES_PER_SEC, ITPLA::K)
    );
    painter.drawText(
                0,
                335,
                QString().sprintf("E:%.6lf,T:%.6lf", ITPLA::E, ITPLA::E / points.size())
    );
//    if (time++ < 120)
        points.back()->GetWorld()->Step(1.0 / 60, 6, 2);
    Points res = ttt.first;
    vector<double> ang = ttt.second;
    //painter.translate(10, 10);
    //painter.scale(1,1);
    //painter.setWindow(0, 0, 2, 2);
    //painter.setViewport(0, 0, this->width(), this->height());

    painter.setPen(Qt::black);
    for (int i = 0; i < polygon.size(); i++) {
        Point &p1 = polygon[i],
              &p2 = polygon[(i + 1) % polygon.size()];
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
        painter.drawLine(p.x, p.y, p.x + r * sin(a), p.y + r * cos(a));
        painter.drawLine(p.x, p.y, p.x + r * sin(a + ttt), p.y + r * cos(a + ttt));
        painter.drawLine(p.x, p.y, p.x + r * sin(a + 2 * ttt), p.y + r * cos(a + 2 * ttt));
    }
    painter.setPen(Qt::black);
    for (int i = 0; i < points.size(); i++) {
        const Point &p = res[i];
        const Vector &v = vel[i];
        painter.drawLine(p.x, p.y, p.x + v.x, p.y + v.y);
    }
    painter.setPen(Qt::blue);
    for (int i = 0; i < ang.size(); i++) {
        painter.setPen(i == attention || attention == -1 ? Qt::blue : Qt::white);
        double a = to_rad(ang[i] + 60);
        double ttt = to_rad(120);
        Point &p = res[i],
               p0 = b2Vec2(p.x + 2 * r * sin(a), p.y + 2 * r * cos(a)),
               p1 = b2Vec2(p.x + 2 * r * sin(a + ttt), p.y + 2 * r * cos(a + ttt)),
               p2 = b2Vec2(p.x + 2 * r * sin(a + 2 * ttt), p.y + 2 * r * cos(a + 2 * ttt));
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
