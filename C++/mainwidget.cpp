#include "mainwidget.h"
#include <QPainter>
#include <QTimer>
#include <QThread>
#include <QMutex>
#include "ITPLA.h"
using namespace ITPLA;
Points polygon, normalized_polygon;
double edge_length;
vector<b2Body *> points;
const int attention = -1;
const int FRAMES_PER_SEC = 60;

class placement_thread : public QThread {
//    Q_OBJECT
signals:
public:
    bool status = false;
    double run_time = 0;
    QMutex lock;
    struct {
        int status = 0;
        std::pair<Points, std::vector<double> > ttt;
        Vectors vel;
        double K, E;
    } buffer[3];
    void run() {
        double start_time = clock() * 1.0 / CLOCKS_PER_SEC;
        while (status && ITPLA::calc_next_step(normalized_polygon, points)) {
            points.back()->GetWorld()->Step(1.0 / FRAMES_PER_SEC, 6, 2);
            lock.lock();
            int writing = 0;
            while (writing < 2 && buffer[writing].status != 1)
                writing++;
            int idle = 0;
            while (idle < 2 && buffer[idle].status)
                idle++;
            buffer[writing].status = 0;
            buffer[idle].status = 1;
            run_time = clock() * 1.0 / CLOCKS_PER_SEC - start_time;
            lock.unlock();
            buffer[idle].ttt.first.clear();
            buffer[idle].ttt.second.clear();
            buffer[idle].vel.clear();
            for (int i = 0; i < points.size(); i++) {
                buffer[idle].ttt.first.push_back(ITPLA::normalize_point(points[i]->GetPosition(), sqrt(3) / edge_length));
                //        printf("%lf,%lf\n",points[i]->GetAngle(),points[i]->GetPosition().y);
                buffer[idle].ttt.second.push_back(to_deg(points[i]->GetAngle()));
                buffer[idle].vel.push_back(ITPLA::normalize_point(points[i]->GetLinearVelocity(), sqrt(3) / edge_length));
            }
            buffer[idle].K = ITPLA::K;
            buffer[idle].E = ITPLA::E;
        }
    }
} *pt;

MainWidget::MainWidget(QWidget *parent) :
    QWidget(parent) {
//    edge_length = 63.17796;
//    polygon = tr1_2;

//    edge_length = 85.86865;
//    polygon = tr1_4;

//    edge_length = 68.983385775;
//    polygon = tr1_5;
    for (int i = 1; i < polygon.size(); i++)
        polygon[i] = {polygon[i-1].x+polygon[i].x, polygon[i-1].y+polygon[i].y};

    edge_length = 158.88;
    polygon = tr1_1;

    for (int i = 0; i < polygon.size()/2; i++)
        swap(polygon[i], polygon[polygon.size()-1-i]);

    edge_length = 100;
    polygon = test;

    edge_length = 99.9533;
    polygon = LH;

    normalized_polygon = normalize_polygon(polygon, edge_length / sqrt(3));
    points = place(polygon, edge_length);
    QTimer *timer = new QTimer();
    timer->start(1000.0 / FRAMES_PER_SEC);
    connect(timer, SIGNAL(timeout()), this, SLOT(repaint()));
    pt = new placement_thread();
    pt->status = true;
    pt->start();
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
                QString().sprintf("Run Time:%8.3lfs", pt->run_time)
    );
    painter.drawText(
                0,
                335,
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
//    if (frame % 1000 == 0) {
//        QPixmap qi(this->width(), this->height());
//        paint(&qi);
//        qi.save(QString().sprintf("/home/zero/%d.bmp", frame));
//    }
    static bool output = false;
    if (!output && pt->isFinished()) {
        QPixmap qi(this->width(), this->height());
        paint(&qi);
        qi.save(QString().sprintf("/Users/zero/%d.png", stime));
        output = true;
    }
    QPainter painter(this);

#if 1
    pt->lock.lock();
    int reading = 0;
    while (reading < 2 && pt->buffer[reading].status != -1)
        reading++;
    int idle = 0;
    while (idle < 2 && pt->buffer[idle].status)
        idle++;
    pt->buffer[reading].status = 0;
    pt->buffer[idle].status = -1;
    pt->lock.unlock();

    pair<Points, vector<double> > &ttt = pt->buffer[idle].ttt;
    Vectors &vel = pt->buffer[idle].vel;
    double K = pt->buffer[idle].K,
            E = pt->buffer[idle].E;
#else
    pair<Points, vector<double> > ttt;
    Vectors vel;
    for (int i = 0; i < points.size(); i++) {
        ttt.first.push_back(normalize_point(points[i]->GetPosition(), sqrt(3) / edge_length));
//        printf("%lf,%lf\n",points[i]->GetAngle(),points[i]->GetPosition().y);
        ttt.second.push_back(to_deg(points[i]->GetAngle()));
        vel.push_back(normalize_point(points[i]->GetLinearVelocity(), sqrt(3) / edge_length));
    }
    calc_next_step(normalized_polygon, points);
    //    if (time++ < 120)
            points.back()->GetWorld()->Step(1.0 / 60, 6, 2);
#endif
    painter.setPen(Qt::black);
    painter.drawText(
                0,
                305,
                QString().sprintf("Run Time:%8.3lfs", pt->run_time)
    );
    painter.drawText(
                0,
                320,
                QString().sprintf("Point:%d,Frame:%6d,Time:%8.3lfs,K:%.6lf", vel.size(), frame, 1.0*frame/FRAMES_PER_SEC, K)
    );
    painter.drawText(
                0,
                335,
                QString().sprintf("E:%.6lf,avg E:%.6lf", E, E / vel.size())
    );
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
    for (int i = 0; i < res.size(); i++) {
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

MainWidget::~MainWidget() {
    pt->status = false;
    while (pt->isRunning());
    delete pt;
}
