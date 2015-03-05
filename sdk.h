#ifndef SDK_H
#define SDK_H
#include <vector>
#include <climits>
#include <cassert>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <Box2D/Box2D.h>
#define Point b2Vec2
#define Points vector<Point >
#define Vector Point
#define Vectors Points
#define ZERO 1e-9
#define CURRENT_TIME (10.0 * clock() / CLOCKS_PER_SEC)
#define show_time() printf("[%.3lf]", CURRENT_TIME)
#define pi b2_pi
#define to_rad(x) ((x) / 180.0 * pi)
#define to_arc(x) ((x) * 180.0 / pi)

namespace sdk {
using namespace std;

Points test() {
    Points polygon;
    polygon.push_back(Point(0, 50));
    polygon.push_back(Point(0, 0));
//    polygon.push_back(Point(600, -100));
    polygon.push_back(Point(400, 0));
    polygon.push_back(Point(300, 300));
//    polygon.push_back(Point(200, 175));
    return polygon;
}

Point rand_point(const Point &plb, const Point &prt) {
    assert(plb.x <= prt.x && plb.y <= prt.y);
    return Point((prt.x - plb.x) * rand() / RAND_MAX + plb.x, (prt.y - plb.y) * rand() / RAND_MAX + plb.y);
}

double normalize_angle(double angle) {  //-180~+180
    while (180 < abs(angle))
        if (0 < angle)
            angle -= 360;
        else
            angle += 360;
    return angle;
}

double angle_diff(double init, double fin) {
    return normalize_angle(fin - init);
}

double distance_to_line(const Point &p, const Point &u, const Point &v) {
    Vector e1 = u - v,
           e2 = p - v;
    return abs(e1.x * e2.y - e1.y * e2.x) / e1.Length();
}

Point normalize_point(const Point &p, const double edge_length) {
    return Point(p.x / edge_length, p.y / edge_length);
}

Points normalize_polygon(const Points &polygon, const double edge_length) {
    Points normalized_polygon;
    for (int i = 0; i < polygon.size(); i++)
        normalized_polygon.push_back(normalize_point(polygon[i], edge_length));
    return normalized_polygon;
}

bool in_polygon(const Points &polygon/*normalized_polygon*/, const Point &p/*normalized_point*/) {
    bool result = false;
    for (int i = 0; i < polygon.size(); i++) {
        const Point &p1 = polygon[i],
                    &p2 = polygon[(i + 1) % polygon.size()];
        if ((p1.y - p.y) * (p2.y - p.y) <= 0
            && ZERO < abs(p1.y - p2.y)
            && !(p.x < p1.x && abs(p.y - p1.y) < ZERO)
            && p1.x + (p.y - p1.y) * (p1.x - p2.x) / (p1.y - p2.y) <= p.x)
            if (abs(p1.x + (p.y - p1.y) * (p1.x - p2.x) / (p1.y - p2.y) - p.x) < ZERO)
                return true;
            else
                result = !result;
    }
    return result;
}

double K;
int frame;
void calc_next_step(const Points &normalized_polygon, /*const*/ vector<b2Body *> &points) {
    assert(0 < points.size());

    vector<vector<int> > nearest_point(points.size(), vector<int>(3, -1));

    for (int i = 0; i < points.size(); i++) {
        b2Body *p1 = points[i];

        for (int j = 0; j < points.size(); j++) {
            b2Body *p2 = points[j];

            if (i == j)
                continue;

            Vector v = p2->GetPosition() - p1->GetPosition();
            int k = -1;
            for (int l = 0; l < 3; l++) {
                double tmp = normalize_angle(to_arc(p1->GetAngle()) + 120 * l);
                Vector n = Point(sin(to_rad(tmp)), cos(to_rad(tmp)));
                if (cos(to_rad(60)) * v.Length() <= v.x * n.x + v.y * n.y)
                    k = l;
            }
            assert(k != -1);
            if (nearest_point[i][k] == -1)
                nearest_point[i][k] = j;
            else {
                b2Body *p3 = points[nearest_point[i][k]];
                if (v.Length() < (p1->GetPosition() - p3->GetPosition()).Length())
                    nearest_point[i][k] = j;
            }
        }
    }

    Vectors force(points.size(), Point(0, 0));
    vector<double> angle(points.size(), 0);

    K = 1;
    int K_n = -1;
    for (int i = 0; i < points.size(); i++) {
        b2Body *p1 = points[i];

        double weight_sum = 0;

        for (int k = 0; k < 3; k++) {
            const double ak = to_arc(p1->GetAngle()) + k * 120;

            if (nearest_point[i][k] != -1) {
                const int &j = nearest_point[i][k];
                b2Body *p2 = points[j];

                Vector v = p2->GetPosition() - p1->GetPosition(),
                       n = Point(sin(to_rad(ak)), cos(to_rad(ak))),
                       t = Point(sin(to_rad(ak + 90)), cos(to_rad(ak + 90)));

                Point p1_line_middle = 0.5 * n;
                int l = -1;
                for (int m = 0; m < 3 && l == -1; m++) {
                    if (abs(angle_diff(ak, to_arc(p2->GetAngle()) + m * 120 + 180)) <= 60)
                        l = m;
                }
                assert(l != -1);
//                double dii = INT_MAX;
//                for (int l = 0; l < 3; l++) {
//                    Point p2_line_middle = v + 0.5 * Point(sin(to_rad(120 * l) + p2->GetAngle()), cos(to_rad(120 * l) + p2->GetAngle()));
//                    if ((p1_line_middle - p2_line_middle).Length() < dii) {
//                        dii = (p1_line_middle - p2_line_middle).Length();
//                        min12 = l;
//                    }
//                }
                Point p2_line_middle = v + 0.5 * Point(sin(to_rad(120 * l) + p2->GetAngle()), cos(to_rad(120 * l) + p2->GetAngle()));
                double ang_diff = angle_diff(ak, to_arc(p2->GetAngle()) + l * 120 + 180);

                double min_distance = 0.5 + sin(to_rad(30 + abs(ang_diff)));
                double v_n_length = v.x * n.x + v.y * n.y,
                       v_t_length = v.x * t.x + v.y * t.y,
                       kn = 1 - pow(v_n_length / min_distance, -2),
                       kt = 0.5 * v_t_length / v_n_length / sqrt(3) * min(1.0, pow(v_n_length / min_distance, -2)) * cos(to_rad(ang_diff));
                assert(0 <= v_n_length);
                double weight = pow(v_n_length / min_distance, -2) + pow((p1_line_middle - p2_line_middle).Length() + 0.1, -2);
                force[i] += weight * (kn * n + kt * t);
                angle[i] += ang_diff / 2 * weight;
                weight_sum += weight;
                if (v_n_length / min_distance /** cos(to_rad(ang_diff))*/ < K) {
                    K = v_n_length / min_distance;
                    K_n = i;
                }
                if (cos(to_rad(ang_diff)) < 0)
                    printf("%lf\n", ang_diff);
            }
        }

        for (int j = 0; j < normalized_polygon.size(); j++) {
            const Point &u = normalized_polygon[j],
                        &v = normalized_polygon[(j + 1) % normalized_polygon.size()],
                        &w = normalized_polygon[(j + 2) % normalized_polygon.size()];
            Vector e1 = u - v,
                   e2 = v - w;
            double dis = distance_to_line(p1->GetPosition(), u, v);
            Vector n = Point(u.y - v.y, v.x - u.x);
            int k = -1;
            double min_dist = INT_MAX;
            for (int l = 0; l < 3; l++) {
                const double al = normalize_angle(to_arc(p1->GetAngle()) + l * 120);
                Point t = p1->GetPosition() + 0.5 * Point(sin(to_rad(al)), cos(to_rad(al)));
                double dist = distance_to_line(t, u, v);
                if (dist < min_dist) {
                    min_dist = dist;
                    k = l;
                }
            }
            assert(k != -1);
            double ang_diff = angle_diff(to_arc(p1->GetAngle()) + k * 120, 180 - to_arc(atan2(v.y - u.y, v.x - u.x)));
            double min_distance = sin(to_rad(30 + abs(ang_diff)));
            double kn = 1 - pow(dis / min_distance, -2);
            if (0 < kn)
                continue;
            assert(ZERO < abs(n.Normalize()));
            double weight = pow(min_dist + 0.2, -2);
            force[i] -= weight * kn * n;
            angle[i] += ang_diff * weight;
            weight_sum += weight;
            if (dis / min_distance < K) {
                K = dis / min_distance;
                K_n = i;
            }
        }
        if (ZERO < weight_sum) {
            force[i] *= 1 / weight_sum;
            angle[i] /= weight_sum;
        }
    }

    for (int i = 0; i < points.size(); i++) {
        b2Body *p = points[i];
        p->SetLinearVelocity(force[i]);
        p->SetAngularVelocity(to_rad(angle[i]));
    }

    frame++;
    if (frame > 100 && frame--)
        for (int i = 0; i < points.size(); i++) {
            b2Body *p = points[i];
            p->SetLinearVelocity(Vector(0, 0));
            p->SetAngularVelocity(0);
        }
    if (frame % 2000 == 0 && K_n != -1 && K < 0.9) {
        points.back()->GetWorld()->DestroyBody(points[K_n]);
        points[K_n] = points.back();
        points.pop_back();
        for (int i = 0; i < points.size(); i++) {
            b2Body *p = points[i];
            p->SetLinearVelocity(Vector(0, 0));
            p->SetAngularVelocity(0);
        }
    }
}

vector<b2Body *>/*pair<Points, vector<double> >*/ place(const Points &polygon, double edge_length) {
    assert(2 < polygon.size());
    const Points normalized_polygon = normalize_polygon(polygon, edge_length / sqrt(3));

    show_time();
    printf("start place ...\n");
    Point plb = normalized_polygon[0],
          prt = normalized_polygon[0];
    for (int i = 1; i < normalized_polygon.size(); i++) {
        const Point &p = normalized_polygon[i];
        plb.x = min(plb.x, p.x);
        plb.y = min(plb.y, p.y);
        prt.x = max(prt.x, p.x);
        prt.y = max(prt.y, p.y);
    }
    printf("(%lf, %lf)<->(%lf, %lf)\n", plb.x, plb.y, prt.x, prt.y);

    const int xxx = 3;-1;10;
    time_t stime = 1425560927;
    stime = (stime != -1 ? stime : time(NULL));
    srand(stime);//time(NULL));
    int total = 0xfff,
        inside = 0;
    for (int i = 0; i < total; i++)
        inside += in_polygon(normalized_polygon, rand_point(plb, prt));
    double area = (prt.x - plb.x) * (prt.y - plb.x) * inside / total;
    printf("approximate area = %.6lf\n", area);

    show_time();
    printf("create world ...\n", area);
    Vector gravity(0, 0);
    b2World *world = new b2World(gravity);

    b2BodyDef border_def;
    border_def.position.Set(0, 0);
    b2Body *border = world->CreateBody(&border_def);
    for (int i = 0; i < normalized_polygon.size(); i++) {
        const Point &u = normalized_polygon[i],
                    &v = normalized_polygon[(i + 1) % normalized_polygon.size()];
        b2EdgeShape border_edge_shape;
        border_edge_shape.Set(u, v);
        border->CreateFixture(&border_edge_shape, 0);
    }

    int point_number = xxx == -1 ? int(area / (3 * sqrt(3) / 4)) : xxx;
    vector<b2Body *> points(point_number, NULL);
    for (int i = 0; i < points.size(); i++) {
        Point p;
        while (!in_polygon(normalized_polygon, p = rand_point(plb, prt)));

        b2BodyDef point_def;
        point_def.type = b2_dynamicBody;
        point_def.position.Set(p.x, p.y);
        point_def.angle = 2 * pi * rand() / RAND_MAX;
        points[i] = world->CreateBody(&point_def);
        b2CircleShape point_shape;
        point_shape.m_p.Set(0, 0);
        point_shape.m_radius = 0.1;
        points[i]->CreateFixture(&point_shape, 1);
    }

    show_time();
    printf("start evolove ...\n");

//    float32 timeStep = 1.0f / 60.0f;
//    int32 velocityIterations = 6;
//    int32 positionIterations = 2;

//    for (int i = 0; i < 600; i++) {
//        calc_next_step(points);
//        world->Step(timeStep, velocityIterations, positionIterations);
//    }

    show_time();
    printf("end evolove\n");
    printf("contain %d points\n", points.size());
    printf("stime = %d\n", stime);
    return points;
}

}

#endif // SDK_H
