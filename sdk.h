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
#define to_rad(x) ((x) / 180.0 * 3.1415926535)
#define to_arc(x) ((x) * 180.0 / 3.1415926535)
#define length(v) sqrt(pow(v.x, 2) + pow(v.y, 2))

namespace sdk {
using namespace std;

Points test() {
    Points polygon;
    polygon.push_back(make_pair(0, 0));
    polygon.push_back(make_pair(400, 0));
    polygon.push_back(make_pair(300, 300));
//    polygon.push_back(make_pair(200, 175));
    polygon.push_back(make_pair(0, 50));
    return polygon;
}

Point rand_point(const Point &plb, const Point &prt) {
    assert(plb.x <= prt.x && plb.y <= prt.y);
    return make_pair((prt.x - plb.x) * rand() / RAND_MAX + plb.x, (prt.y - plb.y) * rand() / RAND_MAX + plb.y);
}

double normalize_angle(double angle) {  //-180~+180
    while (180 < abs(angle))
        if (0 < angle)
            angle -= 360;
        else
            angle += 360;
    return angle;
}

double distance_to_line(const Point &p, const Point &u, const Point &v) {
    Vector e1 = make_pair(u.x - v.x, u.y - v.y),
           e2 = make_pair(p.x - v.x, p.y - v.y);
    return abs(e1.x * e2.y - e1.y * e2.x) / length(e1);
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

Point move_in_polygon(const Points &normalized_polygon, const Point &p, const Vector &v) {
    const double v_length = sqrt(pow(v.x, 2) + pow(v.y, 2));
    if (v_length < ZERO)
        return p;
    //printf("P(%lf, %lf) V(%lf, %lf)\n", p.x, p.y, v.x, v.y);
    if (!in_polygon(normalized_polygon, p)) {
        printf("(%.16lf, %.16lf) not in polygon.\n", p.x, p.y);
        exit(0);
    } //else
        //printf("(%.24lf, %.24lf) in polygon.", p.x, p.y);
    assert(in_polygon(normalized_polygon, p));
    double t = 1;
    Vector vn = make_pair(0, 0);
    Point pn = make_pair(p.x + v.x, p.y + v.y);
    for (int i = 0; i < normalized_polygon.size(); i++) {
        const Point &p1 = normalized_polygon[i],
                    &p2 = normalized_polygon[(i + 1) % normalized_polygon.size()];
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
    //if (ZERO < abs(t - 1))
        //t *= (1 - ZERO);
        //t -= ZERO;
    if (ZERO < t)
        t -= ZERO;
    else
        t = 0;
    pn = make_pair(p.x + v.x * t, p.y + v.y * t);
    return move_in_polygon(normalized_polygon, pn, vn);
}

pair<Points, vector<double> > place(const Points &polygon, double edge_length) {
    assert(2 < polygon.size());
    const Points normalized_polygon = normalize_polygon(polygon, edge_length / sqrt(3));

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

    const int xxx = -1;
    time_t stime = 1425350902;-1;1425207538;1425215711;//-1;1425214647;-1;1425210080;//1425208981;//1425207950;//1425207538;
    //1425196251;
    stime = (stime != -1 ? stime : time(NULL));
    srand(stime);//time(NULL));
    int total = 0xfff,
        inside = 0;
    for (int i = 0; i < total; i++)
        inside += in_polygon(normalized_polygon, rand_point(plb, prt));
    double area = (prt.x - plb.x) * (prt.y - plb.x) * inside / total;
    printf("approximate area = %.6lf\n", area);

    int point_number = xxx == -1 ? int(area / (3 * sqrt(3) / 4)) : xxx;
    Points point_set;
    vector<double> angle_set;
    for (int i = 0; i < point_number; i++) {
        Point p;
        while (!in_polygon(normalized_polygon, p = rand_point(plb, prt)));
        point_set.push_back(p);
        angle_set.push_back(normalize_angle(360.0 * rand() / RAND_MAX));
    }

    show_time();
    printf("start evolove\n");
    double tmp = 0.1;
    for (;;) {
        int number = 0;
        double max_distance = 0;
    double max_f_length_global = INT_MAX;
//    point_set[1] = point_set.back();
//    angle_set[1] = angle_set.back();
//    point_set.pop_back();
//    angle_set.pop_back();
        for (int n = 0; n < 1000; ) {
//            if (n == -1 && point_set.size() == 3) {
//                point_set.pop_back();
//                angle_set.pop_back();
//                n -= 63;
//            }
            vector<vector<int> > nearest_point(point_set.size(), vector<int>(3, -1));
            for (int i = 0; i < point_set.size(); i++) {
                const Point &p1 = point_set[i];
                const double &a = angle_set[i];

                for (int j = 0; j < point_set.size(); j++) {
                    const Point &p2 = point_set[j];

                    if (i == j)
                        continue;

                    Vector v = make_pair(p2.x - p1.x, p2.y - p1.y);
                    int k = -1;
                    for (int l = 0; l < 3; l++) {
                        double tmp = normalize_angle(a + 120 * l);
                        Vector n = make_pair(sin(to_rad(tmp)), cos(to_rad(tmp)));
                        if (cos(to_rad(60)) * length(v) <= v.x * n.x + v.y * n.y)
                            k = l;
                    }
                    assert(k != -1);
                    if (nearest_point[i][k] == -1)
                        nearest_point[i][k] = j;
                    else {
                        const Point &p3 = point_set[nearest_point[i][k]];
                        /*double a12 = abs(angle_set[i] - angle_set[j]),
                               a13 = abs(angle_set[i] - angle_set[nearest_point[i][k]]);
                        while (60 < a12)
                            a12 = abs(a12 - 120);
                        while (60 < a13)
                            a13 = abs(a13 - 120);*/
                        if (length(v) < length(make_pair(p1.x - p3.x, p1.y - p3.y)))
                            nearest_point[i][k] = j;
                    }
                }
            }
            Vectors f(point_set.size(), make_pair(0, 0));
            vector<double> da(point_set.size(), 0);
            double max_repulsion_force = 0;
            max_distance = 0;
            for (int i = 0; i < point_set.size(); i++) {
                const Point &p1 = point_set[i];
                const double &a1 = angle_set[i];
                Vector &fi = f[i];

                for (int j = 0; j < normalized_polygon.size(); j++) {
                    const Point &u = normalized_polygon[j],
                                &v = normalized_polygon[(j + 1) % normalized_polygon.size()],
                                &w = normalized_polygon[(j + 2) % normalized_polygon.size()];
                    Vector e1 = make_pair(u.x - v.x, u.y - v.y),
                           e2 = make_pair(v.x - w.x, v.y - w.y);
                    double dis = distance_to_line(p1, u, v);
                    if (1 < dis || (e1.x * e2.y - e1.y * e2.x < 0 && 1 < distance_to_line(p1, v, w)))
                        continue;
                    if (i == -1) {
                        printf("dis1 = %lf\n", dis);
                        printf("u(%lf, %lf)\n", u.x, u.y);
                        printf("v(%lf, %lf)\n", v.x, v.y);
                        printf("p[%d](%lf, %lf)\n", i, p1.x, p1.y);
                    }
//                    if (abs(dis) < ZERO)
//                        return make_pair(normalize_polygon(point_set, sqrt(3) / edge_length), angle_set);
                    dis = max(dis, 0.1);
                    Vector n = make_pair(u.y - v.y, v.x - u.x);
                    double kn = 1 - pow(dis / 1, -2);
                    if (0 < kn) {
                        kn = 0;
                        int k = -1;
                        double min_dist = INT_MAX;
                        for (int l = 0; l < 3; l++) {
                            const double al = normalize_angle(a1 + l * 120);
                            Point t = make_pair(0.5 * sin(to_rad(al)) + p1.x, 0.5 * cos(to_rad(al)) + p1.y);
//                            Vector el = make_pair(t.x - v.x, t.y - v.y);
                            double dist = distance_to_line(t, u, v);//abs(e1.x * el.y - e1.y * el.x) / length(e1);
                            if (dist < min_dist) {
                                min_dist = dist;
                                k = l;
                            }
                        }
                        assert(k != -1);
                        double aaa = normalize_angle(a1 + k * 120);
                        double tmp = (- to_arc(atan2(v.y - u.y, v.x - u.x)) - aaa);
                        while (90 <= abs(tmp))
                            if (0 < tmp)
                                tmp -= 90;
                            else
                                tmp += 90;
                        double min_distance = sin(to_rad(30 + abs(tmp)));
                        kn = 1 - pow(dis / min_distance, -2);
                        if (0 < kn)
                            kn = 0;
                        else {
                            da[i] += 0.5*tmp / pow(2 * dis, 2);
                        }
                        max_distance = max(max_distance, - dis + min_distance);
                    }
//                    kn *= 10;
                    if (max_repulsion_force < kn) {
                        max_repulsion_force = kn;
                        number = i;
                    }
                    max_distance = max(max_distance, - dis + 0.5);
                    assert(ZERO < abs(length(n)));
                    fi = make_pair(fi.x - kn * n.x / length(n), fi.y - kn * n.y / length(n));
                }

                for (int k = 0; k < 3; k++) {
                    const double ak = normalize_angle(a1 + k * 120);

                    if (nearest_point[i][k] != -1) {
                        const int &j = nearest_point[i][k];
                        const Point &p2 = point_set[j];
                        const double &a2 = angle_set[j];

                        double doa = normalize_angle(a2 - a1 - 180);
                        for (int l = 0; l < 3; l++)
                            if (abs(doa = normalize_angle(ak - (a2 + l * 120 + 180))) <= 60)
                                break;

                        Vector v = make_pair(p2.x - p1.x, p2.y - p1.y),
                               n = make_pair(sin(to_rad(ak)), cos(to_rad(ak))),
                               t = make_pair(sin(to_rad(ak + 90)), cos(to_rad(ak + 90)));

                        double min_distance = 0.5 + sin(to_rad(30 + abs(doa)));
                        double v_n_length = v.x * n.x + v.y * n.y,
                               v_t_length = v.x * t.x + v.y * t.y,
                               kn = 1 - pow(v_n_length / min_distance, -2),
                               kt = v_t_length * min(0.5, pow(v_n_length / min_distance, -3));
                        assert(0 <= v_n_length);
//                        kn *= 10;
                        if (v_n_length < 1) {
//                            kt = 0;
//                            kt = kn / v_n_length * sqrt(pow(length(v), 2) - pow(v_n_length, 2)) * (v_t_length / abs(v_t_length));
//                            kt = kn * tan(to_rad(doa));
                        }
                        if (max_repulsion_force < sqrt(pow(kn, 2) + pow(kt, 2))) {
                            max_repulsion_force = sqrt(pow(kn, 2) + pow(kt, 2));
                            number = i;
                        }
//                        if (k < 0)
//                            k *= 0.1;//1.0 / (n + 1);
                        fi = make_pair(fi.x + kn * n.x + kt * t.x, fi.y + kn * n.y + kt * t.y);
                        da[i] -= doa / 2 / max(1.0, pow(length(v), 2));
                        //f2 = make_pair(f2.x - k * v.x, f2.y - k * v.y);

                        if (i == -1)
                            printf("%d(%lf)exist %d(%lf), da[i] = %lf, n=%lf, t=%lf\n", i, ak, j, doa, da[i], kn, kt);
//                        printf("p[%d] <-> p[%d]\n", i, j);
                        max_distance = max(max_distance, - length(v) + min_distance);
                    }
                }
            }
            double max_f_length = 0;
            for (int i = 0; i < f.size(); i++) {
                max_f_length = max(max_f_length, length(f[i]));
                if (i == -1)
                    printf("f[%d](%lf, %lf)\n", i, f[i].x, f[i].y);
            }
            if (max_f_length < max_f_length_global) {
                max_f_length_global = max_f_length;
                n = 0;
            } else {
                n++;
            }
            if (max_f_length < 0.01)
                break;
            if (n > 995)
                printf("[%d]max_f_length = %lf\n", n, max_f_length);
                //return normalize_polygon(point_set, 1 / edge_length);
//            double k = 0.01 * min(1.0, 1 / max_f_length);
            double k = 0.02 / max_f_length;
            double max_a = 0;
            for (int i = 0; i < angle_set.size(); i++)
                max_a = max(max_a, abs(angle_set[i]));
            double ka = 1;
            for (int i = 0; i < f.size(); i++) {
                Point &p = point_set[i];
                Vector &v = f[i];
                p = move_in_polygon(normalized_polygon, p, make_pair(k * v.x, k * v.y));
                angle_set[i] = normalize_angle(angle_set[i] + ka * da[i]);
//                printf("%lf\n", angle_set[i]);
            }
            //printf("Round: %d, point[%d](%lf, %lf) move(%lf, %lf)", n, i, p1.x, p1.y, f.x, f.y);
            /*double limit = 0.005 * edge_length;
            Point plb = make_pair(-limit, -limit),
                  prt = make_pair(limit, limit);
            p1 = moveInnormalized_Polygon(normalized_polygon, p1, randPoint(plb, prt));*/
            //printf(" now at(%.6lf, %.6lf)\n", p1.x, p1.y);
        }
        if (max_distance < tmp) {
            break;
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
                    angle_set.push_back(normalize_angle(360.0 * rand() / RAND_MAX));
                }
            }
            if (!go_on)
                break;
        } else {
//    return make_pair(normalize_polygon(point_set, sqrt(3) / edge_length), angle_set);
            //continue;
//            break;
            printf("delete point[%d](%lf, %lf)\n", number, point_set[number].x, point_set[number].y);
            point_set[number] = point_set.back();
            angle_set[number] = angle_set.back();
            point_set.pop_back();
            angle_set.pop_back();
        }
    }
    show_time();
    printf("end evolove\n");
    printf("contain %d points\n", point_set.size());
    printf("stime = %d\n", stime);
    return make_pair(normalize_polygon(point_set, sqrt(3) / edge_length), angle_set);
}

}

#endif // SDK_H
