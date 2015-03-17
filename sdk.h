#ifndef SDK_H
#define SDK_H
#include <vector>
#include <climits>
#include <cassert>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <Box2D/Box2D.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#define Point b2Vec2
#define Points vector<Point >
#define Vector Point
#define Vectors Points
#define Point_2s vector<Point_2 >
#define ZERO 1e-9
#define CURRENT_TIME (10.0 * clock() / CLOCKS_PER_SEC)
#define show_time() printf("[%.3lf]", CURRENT_TIME)
#define pi b2_pi
#define to_rad(x) ((x) / 180.0 * pi)
#define to_arc(x) ((x) * 180.0 / pi)
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Triangle_2<K> Triangle_2;
namespace sdk {
using namespace std;

Points test() {
    Points polygon;
    polygon.push_back(Point(0, 50));
    polygon.push_back(Point(0, 0));
//    polygon.push_back(Point(600, -100));
    polygon.push_back(Point(400, 0));
    polygon.push_back(Point(300, 300));
    polygon.push_back(Point(200, 175));
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

Point_2 convert_to_p2(const Point &p) {
    return Point_2(p.x, p.y);
}

Point_2s convert_to_p2s(const Points &polygon) {
    Point_2s polygon_2;
    for (int i = 0; i < polygon.size(); i++)
        polygon_2.push_back(convert_to_p2(polygon[i]));
    return polygon_2;
}

double area_polygon(const Points &polygon) {
    Point_2s polygon_2 = convert_to_p2s(polygon);
    return Polygon_2(&*polygon_2.begin(), &*polygon_2.end()).area();
}

bool in_polygon(const Points &polygon, const Point &p) {
    Point_2s polygon_2 = convert_to_p2s(polygon);
    Point_2 p_2 = convert_to_p2(p);
    switch (CGAL::bounded_side_2(&*polygon_2.begin(), &*polygon_2.end(), p_2, K())) {
        case CGAL::ON_BOUNDED_SIDE :
          return true;
        case CGAL::ON_BOUNDARY:
          return false;
        case CGAL::ON_UNBOUNDED_SIDE:
          return false;
    }
}

Triangle_2 create_triangle(const Point &c, const double arc) {
    Point_2 p0 = convert_to_p2(c + Point(sin(to_rad(arc - 60)), cos(to_rad(arc - 60)))),
            p1 = convert_to_p2(c + Point(sin(to_rad(arc - 180)), cos(to_rad(arc - 180)))),
            p2 = convert_to_p2(c + Point(sin(to_rad(arc + 60)), cos(to_rad(arc + 60))));
    return Triangle_2(p0, p1, p2);
}

bool intersect_each(const Point &c1, const double a1, const Point &c2, const double a2) {
    return CGAL::do_intersect(create_triangle(c1, a1), create_triangle(c2, a2));
}

int calc_direction(const Point &c, const double arc, const Point &p) {
    Vector v = p - c;
    for (int i = 0; i < 3; i++) {
        double ai = arc + 120 * i;
        Vector n = Point(sin(to_rad(ai)), cos(to_rad(ai)));
        if (cos(to_rad(60)) * v.Length() <= v.x * n.x + v.y * n.y)
            return i;
    }
    for (int i = 0; i < 3; i++) {
        double ai = arc + 120 * i;
        Vector n = Point(sin(to_rad(ai)), cos(to_rad(ai)));
        if (cos(to_rad(61)) * v.Length() <= v.x * n.x + v.y * n.y)
            return i;
    }
    return -1;
}

double K, E = INT_MAX, pre_E, min_E = INT_MAX;
int frame;
vector<Point> min_p;
vector<double> min_a;
void calc_next_step(const Points &normalized_polygon, /*const*/ vector<b2Body *> &points) {
    assert(0 < points.size());

    vector<vector<int> > nearest_point(points.size(), vector<int>(3, -1));
    for (int i = 0; i < points.size(); i++) {
        b2Body *p1 = points[i];

        for (int j = 0; j < points.size(); j++) {
            b2Body *p2 = points[j];

            if (i == j)
                continue;

            int k = calc_direction(p1->GetPosition(), to_arc(p1->GetAngle()), p2->GetPosition());
            assert(k != -1);
            if (nearest_point[i][k] == -1)
                nearest_point[i][k] = j;
            else {
                b2Body *p3 = points[nearest_point[i][k]];
                if ((p2->GetPosition() - p1->GetPosition()).Length() < (p1->GetPosition() - p3->GetPosition()).Length())
                    nearest_point[i][k] = j;
            }
        }
    }

    vector<vector<int> > overlap(points.size());
    for (int i = 0; i < points.size(); i++) {
        vector<bool> intersect(points.size(), false);
        b2Body *p1 = points[i];
        for (int j = 0; j < points.size(); j++) {
            b2Body *p2 = points[j];
            if (i == j)
                continue;
            intersect[j] = intersect_each(p1->GetPosition(), to_arc(p1->GetAngle()), p2->GetPosition(), to_arc(p2->GetAngle()));
        }
//        for (int j = 0; j < 3; j++)
//            if (nearest_point[i][j] != -1)
//                intersect[nearest_point[i][j]] = false;
        for (int j = 0; j < points.size(); j++)
            if (intersect[j])
                overlap[i].push_back(j);
    }

    Vectors force(points.size(), Point(0, 0));
    vector<double> angle(points.size(), 0);
    K = 1;
    pre_E = E;
    E = 0;
    vector<pair<pair<int, double>, int> > del_rank;
    for (int i = 0; i < points.size(); i++) {
        del_rank.push_back(make_pair(make_pair(0.0, 0), i));
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
                int l = calc_direction(p2->GetPosition(), to_arc(p2->GetAngle()), p1->GetPosition());
                assert(l != -1);

                const double al = to_arc(p2->GetAngle()) + l * 120;
                Vector n2 = Point(sin(to_rad(al)), cos(to_rad(al)));
                Point p2_line_middle = v + 0.5 * n2;
                double ang_diff = angle_diff(ak, al + 180);

                double v_n_length = v.x * n.x + v.y * n.y,
                       v_n2_length = - (v.x * n2.x + v.y * n2.y),
                       p2_line_middle_t_length = p2_line_middle.x * t.x + p2_line_middle.y * t.y,
                       min_distance = (0.5 + sin(to_rad(30 + abs(ang_diff)))) / (max(v_n_length, v_n2_length) / v.Length()),
                       kr = 1 - pow(v.Length() / min_distance, -2),
                       kt = 0.5 * p2_line_middle_t_length;
                assert(0 <= v_n_length);
                Vector r = 1 / v.Length() * v;
                double weight = pow(v.Length() / min_distance, -2) + pow((p1_line_middle - p2_line_middle).Length() + 0.1, -2);
                force[i] += weight * (kr * r + kt * t);
                angle[i] += 0.5 * ang_diff * weight;
                weight_sum += weight;
                K = min(K, v.Length() / min_distance);
                if ((p1_line_middle - p2_line_middle).Length() < 0.15)
                    del_rank.back().first.first += 10;
                del_rank.back().first.second -= max(0.0, 1 - v.Length() / min_distance);
            }
        }

        for (int k = 0; k < overlap[i].size(); k++) {
            const int &j = overlap[i][k];
            b2Body *p2 = points[j];
            double ak = calc_direction(p1->GetPosition(), to_arc(p1->GetAngle()), p2->GetPosition()) * 120 + to_arc(p1->GetAngle());

            Vector v = p2->GetPosition() - p1->GetPosition(),
                   n = Point(sin(to_rad(ak)), cos(to_rad(ak))),
                   t = Point(sin(to_rad(ak + 90)), cos(to_rad(ak + 90)));

            Point p1_line_middle = 0.5 * n;
            int l = -1;
            for (int m = 0; m < 3 && l == -1; m++) {
                double tmp = to_arc(p2->GetAngle()) + 120 * m;
                Vector n = Point(sin(to_rad(tmp)), cos(to_rad(tmp)));
                if (cos(to_rad(60)) * v.Length() <= -v.x * n.x + -v.y * n.y)
                    l = m;
            }
            assert(l != -1);

            const double al = to_arc(p2->GetAngle()) + l * 120;
            Vector n2 = Point(sin(to_rad(al)), cos(to_rad(al)));
            Point p2_line_middle = v + 0.5 * n2;
            double ang_diff = angle_diff(ak, al + 180);

            double v_n_length = v.x * n.x + v.y * n.y,
                   v_n2_length = - (v.x * n2.x + v.y * n2.y),
                   p2_line_middle_t_length = p2_line_middle.x * t.x + p2_line_middle.y * t.y,
                   min_distance = (0.5 + sin(to_rad(30 + abs(ang_diff)))) / (max(v_n_length, v_n2_length) / v.Length()),
                   kr = 1 - pow(v.Length() / min_distance, -2),
                   kt = 0.5 * p2_line_middle_t_length;
            assert(0 <= v_n_length);
            Vector r = 1 / v.Length() * v;
            double weight = pow(v.Length() / min_distance, -2) + pow((p1_line_middle - p2_line_middle).Length() + 0.1, -2);
//            force[i] += weight * (kr * r);
//            weight_sum += weight;
//            K = min(K, v.Length() / min_distance);
//            del_rank.back().first.second -= max(0.0, 1 - v.Length() / min_distance);
            E += max(0.0, 1 / (v.Length() / min_distance) - 1);
        }

        for (int j = 0; j < normalized_polygon.size(); j++) {
            const Point &u = normalized_polygon[j],
                        &v = normalized_polygon[(j + 1) % normalized_polygon.size()];
            double dis = distance_to_line(p1->GetPosition(), u, v);
            Vector n = Point(u.y - v.y, v.x - u.x);
            int k = -1;
            double min_dist = INT_MAX;
            for (int l = 0; l < 3; l++) {
                const double al = to_arc(p1->GetAngle()) + l * 120;
                Point t = p1->GetPosition() + 0.5 * Point(sin(to_rad(al)), cos(to_rad(al)));
                double dist = distance_to_line(t, u, v);
                if (dist < min_dist) {
                    min_dist = dist;
                    k = l;
                }
            }
            assert(k != -1);
            const double ak = to_arc(p1->GetAngle()) + k * 120;
            double ang_diff = angle_diff(ak, 180 - to_arc(atan2(v.y - u.y, v.x - u.x)));
            double min_distance = sin(to_rad(30 + abs(ang_diff)));

            if (dis < min_distance) {
//                printf("before mis_distance: %lf\n", min_distance);
                const double a0 = to_arc(p1->GetAngle());
                Point t0 = p1->GetPosition() + Point(sin(to_rad(a0 + 60)), cos(to_rad(a0 + 60))),
                      t1 = p1->GetPosition() + Point(sin(to_rad(a0 - 60)), cos(to_rad(a0 - 60))),
                      t2 = p1->GetPosition() + Point(sin(to_rad(a0 - 180)), cos(to_rad(a0 - 180)));
                Points triangle;
                triangle.push_back(t0);
                triangle.push_back(t1);
                triangle.push_back(t2);
                if (in_polygon(triangle, u))
                    for (int l = 0; l < 3; l++) {
                        const Point &tu = triangle[l],
                                    &tv = triangle[(l + 1) % triangle.size()];
                        Vector tn = Point(tu.y - tv.y, tv.x - tu.x);
                        if (0 < tn.x * n.x + tn.y + n.y) {
                            Vector uv = v - u,
                                   utu = tu - u,
                                   utv = tv - u,
                                   tutv = tv - tu;
                            if ((uv.x * utu.x + uv.y * utu.y) * (uv.x * utv.x + uv.y * utv.y) < 0)
                                min_distance = min(min_distance, dis + distance_to_line(u, tu, tv) / abs(cos(atan2(uv.y, uv.x) - atan2(tutv.y, tutv.x))));
                        }
                    }
                if (in_polygon(triangle, v))
                    for (int l = 0; l < 3; l++) {
                        const Point &tu = triangle[l],
                                    &tv = triangle[(l + 1) % triangle.size()];
                        Vector tn = Point(tu.y - tv.y, tv.x - tu.x);
                        if (0 < tn.x * n.x + tn.y + n.y) {
                            Vector uv = v - u,
                                   vtu = tu - v,
                                   vtv = tv - v,
                                   tutv = tv - tu;
                            if ((uv.x * vtu.x + uv.y * vtu.y) * (uv.x * vtv.x + uv.y * vtv.y) < 0)
                                min_distance = min(min_distance, dis + distance_to_line(v, tu, tv) / abs(cos(atan2(uv.y, uv.x) - atan2(tutv.y, tutv.x))));
                        }
                    }
//                printf("after mis_distance: %lf\n", min_distance);
                bool intersect = in_polygon(triangle, u) || in_polygon(triangle, v);
                for (int l = 0; l < 3 && !intersect; l++) {
                    const Point &tu = triangle[l],
                                &tv = triangle[(l + 1) % triangle.size()];
                    Vector uv = v - u,
                           vtu = tu - v,
                           vtv = tv - v,
                           tutv = tv - tu,
                           tvu = u - tv,
                           tvv = v - tv;
                    if ((uv.x * vtu.y - uv.y * vtu.x) * (uv.x * vtv.y - uv.y * vtv.x) < 0
                     && (tutv.x * tvu.y - tutv.y * tvu.x) * (tutv.x * tvv.y - tutv.y * tvv.x) < 0)
                        intersect = true;
                }
                if (!intersect)
                    continue;
            }

            double kn = 1 - pow(dis / min_distance, -2);
            if (0 < kn)
                continue;
            assert(ZERO < abs(n.Normalize()));
            double weight = pow(min_dist + 0.2, -2);
            force[i] -= weight * kn * n;
            angle[i] += ang_diff * weight;
            weight_sum += weight;
            K = min(K, dis / min_distance);
            if (min_dist < 0.1)
                del_rank.back().first.first++;
            del_rank.back().first.second -= max(0.0, 1 - dis / min_distance);
            E += max(0.0, 1 / (dis / min_distance) - 1);
        }

        if (ZERO < weight_sum) {
            force[i] *= 1 / weight_sum;
            angle[i] /= weight_sum;
        }
    }
    E /= 2;

    for (int i = 0; i < points.size(); i++) {
        b2Body *p = points[i];
        p->SetLinearVelocity(force[i]);
        p->SetAngularVelocity(to_rad(angle[i]));
    }

    sort(del_rank.begin(), del_rank.end());
    int del = -1;
    for (int i = 0; i < del_rank.size() && del == -1; i++)
        if (ZERO < abs(del_rank[i].first.second))
            del = del_rank[i].second;
    static pair<int, int> pre_del(-1, 0);
    if (pre_del.first == del)
        pre_del.second++;
    else
        pre_del = pair<int, int>(del, 1);
    frame++;
    static int pause_time = 0;
//    if (E < min_E) {
//        min_p.clear();
//        min_a.clear();
//        for (int i = 0; i < points.size(); i++) {
//            b2Body *p = points[i];
//            min_p.push_back(p->GetPosition());
//            min_a.push_back(p->GetAngle());
//        }
//    }
//    min_E = min(min_E, E);
    pause_time += exp(1 - E / pre_E) < (double)rand() / RAND_MAX;
    if (frame > 30000 && frame--)
        for (int i = 0; i < points.size(); i++) {
            b2Body *p = points[i];
            p->SetLinearVelocity(Vector(0, 0));
            p->SetAngularVelocity(0);
        }
    if (pre_del.first != -1 && pow(points.size(), 2) < pause_time && K < 0.85) {
        pre_del.second = 0;
        pause_time = 0;
        pre_E = INT_MAX;
        points.back()->GetWorld()->DestroyBody(points[del]);
        points[del] = points.back();
        points.pop_back();
        for (int i = 0; i < points.size(); i++) {
            b2Body *p = points[i];
            p->SetLinearVelocity(Vector(0, 0));
            p->SetAngularVelocity(0);
        }
    }
}

    const int xxx = -1;//10;
    time_t stime = 1425904342;//-1;//1425813081;//-1;//1425746144;//-1;//1425641876;
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

    stime = (stime != -1 ? stime : time(NULL));
    srand(stime);//time(NULL));

    double area = area_polygon(normalized_polygon);
    printf("accurate area = %.6lf\n", area);

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

/*    float32 timeStep = 1.0f / 60.0f;
    int32 velocityIterations = 6;
    int32 positionIterations = 2;

    for (int i = 0; i < 600; i++) {
        calc_next_step(points);
        world->Step(timeStep, velocityIterations, positionIterations);
    }
*/
    show_time();
    printf("end evolove\n");
    printf("contain %d points\n", points.size());
    printf("stime = %d\n", stime);
    return points;
}

}

#endif // SDK_H
