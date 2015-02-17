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
#define Vectors Points
#define x first
#define y second
#define ZERO 1e-9
#define show_time() printf("[%.3lf]", (double)clock() / CLK_TCK)

namespace sdk {
using namespace std;

Points test() {
    Points polygon;
    polygon.push_back(make_pair(0, 0));
    polygon.push_back(make_pair(500, 0));
    /*polygon.push_back(make_pair(250, 0));
    polygon.push_back(make_pair(200, 300));
    polygon.push_back(make_pair(500, 250));
    */
    //polygon.push_back(make_pair(500, 500));
    //polygon.push_back(make_pair(100, 100));
    polygon.push_back(make_pair(0, 500));
    return polygon;
}

Point rand_point(const Point &plb, const Point &prt) {
    assert(plb.x <= prt.x && plb.y <= prt.y);
    return make_pair((prt.x - plb.x) * rand() / RAND_MAX + plb.x, (prt.y - plb.y) * rand() / RAND_MAX + plb.y);
}

Point normalize_point(const Point &p, const double edge_length) {
    return make_pair(p.x / edge_length, p.y / edge_length);
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

Point move_in_polygon(const Points &polygon/*normalized_polygon*/, const Point &p, const Vector &v) {
    const double v_length = sqrt(pow(v.x, 2) + pow(v.y, 2));
    if (v_length < ZERO)
        return p;
    //printf("P(%lf, %lf) V(%lf, %lf)\n", p.x, p.y, v.x, v.y);
    if (!in_polygon(polygon, p)) {
        printf("(%.16lf, %.16lf) not in polygon.\n", p.x, p.y);
        exit(0);
    } //else
        //printf("(%.24lf, %.24lf) in polygon.", p.x, p.y);
    assert(in_polygon(polygon, p));
    double t = 1;
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
            //printf("case 1\n");
        } else if (n < 0) {
            //printf("case 2\n");
            double t_ = ((p.y - p1.y) * e12.x - (p.x - p1.x) * e12.y) / (v.x * e12.y - v.y * e12.x);
            assert(0 <= t_);
            if (t_ < t) {
                t = t_;
                double vnl = (e12.x * v.x + e12.y * v.y) * (1 - t) / (pow(e12.x, 2) + pow(e12.y, 2));
                vn = make_pair(vnl * e12.x, vnl * e12.y);
            }
        } else if (n < ZERO) {
            //printf("case 3\n");
            Vector e12 = make_pair(p2.x - p1.x, p2.y - p1.y);
            //printf("e12:(%lf, %lf), e:(%lf, %lf)\n", e12.x, e12.y, e.x, e.y);
            if (0 < e12.x * v.y - e12.y * v.x)
                continue;
            t = 0;
            double vnl = (e12.x * v.x + e12.y * v.y) * (1 - t) / (pow(e12.x, 2) + pow(e12.y, 2));
            vn = make_pair(vnl * e12.x, vnl * e12.y);
        }
    }
    if (ZERO < abs(t - 1))
        t *= (1 - ZERO);
    pn = make_pair(p.x + v.x * t, p.y + v.y * t);
    if (sqrt(pow(vn.x, 2) + pow(vn.y, 2)) < ZERO)
        return pn;
    else
        return move_in_polygon(polygon, pn, vn);
}

Points place(const Points &polygon, double edge_length) {
    assert(2 < polygon.size());
    const Points normalized_polygon = normalize_polygon(polygon, edge_length);

    printf("start place\n");
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

    //srand(time(NULL));
    int total = 0xfff,
        inside = 0;
    for (int i = 0; i < total; i++)
        inside += in_polygon(normalized_polygon, rand_point(plb, prt));
    double area = (prt.x - plb.x) * (prt.y - plb.x) * inside / total;
    printf("approximate area = %.6lf\n", area);

    int point_number = int(0.8 * area / (sqrt(3) / 4));
    Points point_set;
    for (int i = 0; i < point_number; i++) {
        Point p;
        while (!in_polygon(normalized_polygon, p = rand_point(plb, prt)));
        point_set.push_back(p);
    }

    show_time();
    printf("start evolove\n");
    double tmp = 0.1;
    for (;;) {
        int number = 0;
        for (int n = 0; n < 500; n++) {
            Vectors f(point_set.size(), make_pair(0, 0));
            double max_repulsion_force = 0;
            for (int i = 0; i < point_set.size(); i++) {
                const Point &p1 = point_set[i];
                Vector &f1 = f[i];
                for (int j = 0; j < i; j++) {
                    const Point &p2 = point_set[j];
                    Vector &f2 = f[j];

                    Vector v = make_pair(p1.x - p2.x, p1.y - p2.y);
                    double v_length = sqrt(pow(v.x, 2) + pow(v.y, 2)),
                           k = pow(v_length, -3) - 1 / v_length;
                    /*if (v_length <= 1) {
                        f1 = make_pair(f1.x + k * v.x, f1.y + k * v.y);
                        f2 = make_pair(f2.x - k * v.x, f2.y - k * v.y);
                        if (max_repulsion_force < k) {
                            max_repulsion_force = k;
                            number = i;
                        }
                    } else {
                        f1 = make_pair(f1.x - k * v.x / 8, f1.y - k * v.y / 8);
                        f2 = make_pair(f2.x + k * v.x / 8, f2.y + k * v.y / 8);
                    }*/
                    if (max_repulsion_force < k) {
                        max_repulsion_force = k;
                        number = i;
                    }
                    if (k < 0)
                        k *= 1.0 / (n + 1);
                    f1 = make_pair(f1.x + k * v.x, f1.y + k * v.y);
                    f2 = make_pair(f2.x - k * v.x, f2.y - k * v.y);
                }
            }
            double max_f_length = 0;
            for (int i = 0; i < f.size(); i++) {
                double f_length = sqrt(pow(f[i].x, 2) + pow(f[i].y, 2));
                max_f_length = max(max_f_length, f_length);
            }
            if (max_f_length < 0.001)
                break;
                //return normalize_polygon(point_set, 1 / edge_length);
            //double k = 0.01 * min(1.0, 1 / max_f_length);
            double k = 0.01 < max_f_length ? 0.01 / max_f_length : 1;
            for (int i = 0; i < f.size(); i++) {
                Point &p = point_set[i];
                Vector &v = f[i];
                p = move_in_polygon(normalized_polygon, p, make_pair(k * v.x, k * v.y));
            }
            //printf("Round: %d, point[%d](%lf, %lf) move(%lf, %lf)", n, i, p1.x, p1.y, f.x, f.y);
            /*double limit = 0.005 * edge_length;
            Point plb = make_pair(-limit, -limit),
                  prt = make_pair(limit, limit);
            p1 = moveInnormalized_Polygon(normalized_polygon, p1, randPoint(plb, prt));*/
            //printf(" now at(%.6lf, %.6lf)\n", p1.x, p1.y);
        }
        double max_distance = 0;
        for (int i = 0; i < point_set.size(); i++) {
            Point &p1 = point_set[i];
            for (int j = 0; j < i; j++) {
                Point &p2 = point_set[j];
                max_distance = max(max_distance, max(0.0, 1 - sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2))));
            }
        }
        if (max_distance < tmp) {
            //break;
            bool go_on = false;
            for (int i = 0; i < normalized_polygon.size(); i++) {
                const Point &p1 = normalized_polygon[i];
                bool k = true;
                for (int j = 0; j < point_set.size() && k; j++) {
                    Point &p2 = point_set[j];
                    k = (1 <= sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2)));
                }
                go_on |= k;
                if (k) {
                    printf("add point[%d](%lf, %lf)\n", point_set.size(), p1.x, p1.y);
                    point_set.push_back(p1);
                }
            }
            if (!go_on)
                break;
        } else {
            printf("delete point[%d](%lf, %lf)\n", number, point_set[number].x, point_set[number].y);
            point_set[number] = point_set.back();
            point_set.pop_back();
        }
    }
    show_time();
    printf("end evolove\n");
    return normalize_polygon(point_set, 1 / edge_length);
}

}

#endif // SDK_H
