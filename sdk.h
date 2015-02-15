#ifndef SDK_H
#define SDK_H
#include <vector>
#include <climits>
#include <cassert>
#include <cmath>
#include <ctime>
#define Point pair<double, double>
#define Points vector<Point >
#define Vector Point
#define x first
#define y second
#define ZERO 1e-6

namespace sdk {
using namespace std;
bool inPolygon(const Points &polygon, const Point &p) {
    bool result = false;
    for (int i = 0; i < polygon.size(); i++) {
        //printf("result = %d, check edge[%d] start...\n", result, i);
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
        /*printf("p = (%lf, %lf), p1 = (%lf, %lf), p2 = (%lf, %lf)\n", p.x,p.y,p1.x,p1.y,p2.x,p2.y);
        printf("%lf\n", (p1.y - p.y) * (p2.y - p.y));
        printf("%lf\n", abs(p1.y - p2.y));
        printf("%lf\n", abs(p.y - p1.y));
        printf("tmp = %.24lf\n", p1.x + (p.y - p1.y) * (p1.x - p2.x) / (p1.y - p2.y));*/
    }
    /*printf("final result = %d\n", result);
    if (result != (0 <= p.x && p.x <= 100 && 0 <= p.y && p.y <= 100)) {
        printf("(%.24lf, %.24lf) ERROR\n", p.x, p.y);
        exit(0);
    }*/
    return result;
}

Point moveInPolygon(const Points &polygon, const Point &p, const Vector &v) {
    //printf("P(%lf, %lf) V(%lf, %lf)\n", p.x, p.y, v.x, v.y);
    if (!inPolygon(polygon, p)) {
        printf("(%.24lf, %.24lf) not in polygon.", p.x, p.y);
        exit(0);
    } //else
        //printf("(%.24lf, %.24lf) in polygon.", p.x, p.y);
    assert(inPolygon(polygon, p));
    double t = 1;
    //const double vrl = sqrt(pow(v.x, 2) + pow(v.y, 2)) / edge_length;// v relative length
    Vector vn = make_pair(0, 0);
    Point pn = make_pair(p.x + v.x, p.y + v.y);
    for (int i = 0; i < polygon.size(); i++) {
        const Point &p1 = polygon[i],
                    &p2 = polygon[(i + 1) % polygon.size()];
        Vector e12 = make_pair(p2.x - p1.x, p2.y - p1.y),
               e = make_pair(p.x - p1.x, p.y - p1.y),
               en = make_pair(pn.x - p1.x, pn.y - p1.y);
        double n = (e12.x * e.y - e12.y * e.x) * (e12.x * en.y - e12.y * en.x);
        if (abs(v.x * e12.y - v.y * e12.x) < ZERO) {
        } else if (n < 0) {
            double t_ = ((p.y - p1.y) * e12.x - (p.x - p1.x) * e12.y) / (v.x * e12.y - v.y * e12.x);
            assert(0 <= t_);
            if (t_ < t) {
                t = t_;
                double vnl = (e12.x * v.x + e12.y * v.y) * (1 - t) / (pow(e12.x, 2) + pow(e12.y, 2));
                vn = make_pair(vnl * e12.x, vnl * e12.y);
            }
        } else if (n < ZERO) {
            double v_length = sqrt(pow(v.x, 2) + pow(v.y, 2));
            Point pt = make_pair(p.x + v.x / v_length * ZERO, p.y + v.y / v_length * ZERO);
            Vector e12 = make_pair(p2.x - p1.x, p2.y - p1.y),
                   e = make_pair(p.x - p1.x, p.y - p1.y);
            if (0 < e12.x * e.y - e12.y * e.x)
                continue;
            t = 0;
            double vnl = (e12.x * v.x + e12.y * v.y) * (1 - t) / (pow(e12.x, 2) + pow(e12.y, 2));
            vn = make_pair(vnl * e12.x, vnl * e12.y);
        }
    }
    if (t < ZERO)
        t = 0;
    else
        t *= (1 - ZERO);
    pn = make_pair(p.x + v.x * t, p.y + v.y * t);
    if (sqrt(pow(vn.x, 2) + pow(vn.y, 2)) < ZERO)
        return pn;
    else
        return moveInPolygon(polygon, pn, vn);
}

Points test() {
    Points polygon;
    polygon.push_back(make_pair(0, 0));
    polygon.push_back(make_pair(500, 0));
    /*polygon.push_back(make_pair(250, 0));
    polygon.push_back(make_pair(200, 300));
    polygon.push_back(make_pair(500, 250));
    */
    polygon.push_back(make_pair(500, 500));
    //polygon.push_back(make_pair(100, 100));
    polygon.push_back(make_pair(0, 500));
    return polygon;
}

Point randPoint(const Point &p1, const Point &p2) {
    assert(p1.x < p2.x && p1.y < p2.y);
    return make_pair((p2.x - p1.x) * rand() / RAND_MAX + p1.x, (p2.y - p1.y) * rand() / RAND_MAX + p1.y);
}

Points place(const Points &polygon, double edge_length) {
    assert(2 < polygon.size());
    printf("start place\n");
    double x_min = polygon[0].x, x_max = polygon[0].x,
           y_min = polygon[0].y, y_max = polygon[0].y;
    for (int i = 1; i < polygon.size(); i++) {
        const Point &p = polygon[i];
        x_min = min(x_min, p.x);
        x_max = max(x_max, p.x);
        y_min = min(y_min, p.y);
        y_max = max(y_max, p.y);
    }

    //srand(time(NULL));
    int total = 0xfff,
        inside = 0;
    for (int i = 0; i < total; i++)
        inside += inPolygon(polygon, randPoint(make_pair(x_min, y_min), make_pair(x_max, y_max)));
    double area = (x_max - x_min) * (y_max - y_min) * inside / total;
    printf("approximate area = %.6lf\n", area);

    int point_number = int(0.8 * area / (sqrt(3) / 4 * pow(edge_length, 2)));
    Points point_set;
    for (int i = 0; i < point_number; i++) {
        Point p;
        while (!inPolygon(polygon, p = randPoint(make_pair(x_min, y_min), make_pair(x_max, y_max))));
        point_set.push_back(p);
    }
    printf("start evolove\n");
    double tmp = 0.1;
        double max_f = 0, max_num = 0;
    for (double max_distance = edge_length; tmp < max_distance / edge_length; (point_set[max_num] = point_set.back()),point_set.pop_back()) {
        max_f = 0;
        max_num = 0;
    for (int n = 0; n < 500; n++)
        for (int i = 0; i < point_set.size(); i++) {
            Point &p1 = point_set[i];
            Vector f = make_pair(0, 0);
            for (int j = 0; j < point_set.size(); j++) {
                Point &p2 = point_set[j];
                if (i != j) {
                    Vector v = make_pair(p1.x - p2.x, p1.y - p2.y);
                    double v_length = sqrt(pow(v.x, 2) + pow(v.y, 2)),
                           k = pow(edge_length / v_length, 3);
                    if (v_length <= edge_length)
                        f = make_pair(f.x + k * v.x, f.y + k * v.y);
                    else if (ZERO < abs(v_length - edge_length))
                        f = make_pair(f.x - k * v.x / 8, f.y - k * v.y / 8);
                }
            }
            double f_length = sqrt(pow(f.x, 2) + pow(f.y, 2)),
                   k = min(0.01 * edge_length / f_length, 1.0 / 27);
            if (max_f < f_length) {
                max_f = f_length;
                max_num = i;
            }
            f = make_pair(f.x * k, f.y * k);
            //printf("Round: %d, point[%d](%lf, %lf) move(%lf, %lf)", n, i, p1.x, p1.y, f.x, f.y);
            p1 = moveInPolygon(polygon, p1, f);
            /*double limit = 0.005 * edge_length;
            Point plb = make_pair(-limit, -limit),
                  prt = make_pair(limit, limit);
            p1 = moveInPolygon(polygon, p1, randPoint(plb, prt));*/
            //printf(" now at(%.6lf, %.6lf)\n", p1.x, p1.y);
        }
        max_distance = 0;
        for (int i = 0; i < point_set.size(); i++) {
            Point &p1 = point_set[i];
            for (int j = 0; j < i; j++) {
                Point &p2 = point_set[j];
                max_distance = max(max_distance, max(0.0, edge_length - sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2))));
            }
        }
    }
    printf("end evolove\n");
    return point_set;
}
}

#endif // SDK_H
