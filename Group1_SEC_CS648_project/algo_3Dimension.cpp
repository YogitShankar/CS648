#include <bits/stdc++.h>
#include <vector>
#include <fstream>
#include <chrono>
#include <random>
using namespace std;

// Helper functions and structs
struct Sphere {
    double x;
    double y;
    double z;
    double r;
};

struct Point3d {
    double x;
    double y;
    double z;
};

std::vector<Point3d> generateRandom3DPoints(int n) {
    static std::default_random_engine generator;
    std::normal_distribution<double> distributionX(0, 10000);
    std::normal_distribution<double> distributionY(0, 10000);
    std::normal_distribution<double> distributionZ(0, 10000);

    std::vector<Point3d> points;
    for (int i = 0; i < n; ++i) {
        double x = distributionX(generator);
        double y = distributionY(generator);
        double z = distributionZ(generator);
        points.push_back({x, y, z});
    }
    return points;
}

double distance(Point3d A, Point3d B) {
    return sqrt(pow(A.x - B.x, 2) + pow(A.y - B.y, 2) + pow(A.z - B.z, 2));
}

Point3d subtract(Point3d A, Point3d B) {
    return {A.x - B.x, A.y - B.y, A.z - B.z};
}

// Function to calculate the cross product of two vectors
Point3d crossProduct(Point3d A, Point3d B) {
    return {A.y * B.z - A.z * B.y,
            A.z * B.x - A.x * B.z,
            A.x * B.y - A.y * B.x};
}

// Function to calculate the magnitude of a vector
double magnitude(Point3d A) {
    return sqrt(A.x * A.x + A.y * A.y + A.z * A.z);
}

// Functions to check if points are inside the specified shape
bool is_inside(const Sphere& sphere, const Point3d& point) {
    return (point.x - sphere.x) * (point.x - sphere.x) + (point.y - sphere.y) * (point.y - sphere.y) + (point.z - sphere.z) * (point.z - sphere.z) <= sphere.r * sphere.r + 1e-4;
}

// make_sphere functions return smallest sphere passing through given set of points
Sphere make_sphere(const Point3d& p1, const Point3d& p2) {
    return { (p1.x + p2.x) / 2, (p1.y + p2.y) / 2, (p1.z + p2.z) / 2, sqrt((p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y) + (p2.z - p1.z) * (p2.z - p1.z)) / 2 };
}

Sphere make_sphere(const Point3d& a, const Point3d& b, const Point3d& c) {
    Point3d ac = subtract(c, a);
    Point3d ab = subtract(b, a);

    // Calculate cross product of ab and ac
    Point3d abXac = crossProduct(ab, ac);

    // Calculate vector from a to the circumsphere center
    Point3d toCircumsphereCenter = crossProduct(abXac, ab);
    toCircumsphereCenter.x *= magnitude(ac) * magnitude(ac);
    toCircumsphereCenter.y *= magnitude(ac) * magnitude(ac);
    toCircumsphereCenter.z *= magnitude(ac) * magnitude(ac);
    
    Point3d cross_abXac = crossProduct(abXac, ac);
    cross_abXac.x *= magnitude(ab) * magnitude(ab);
    cross_abXac.y *= magnitude(ab) * magnitude(ab);
    cross_abXac.z *= magnitude(ab) * magnitude(ab);
    
    toCircumsphereCenter = subtract(toCircumsphereCenter, cross_abXac);
    
    toCircumsphereCenter.x /= 2 * pow(magnitude(abXac), 2);
    toCircumsphereCenter.y /= 2 * pow(magnitude(abXac), 2);
    toCircumsphereCenter.z /= 2 * pow(magnitude(abXac), 2);

    // Calculate the radius of the circumsphere
    double circumsphereRadius = magnitude(toCircumsphereCenter);

    // Calculate the 3D coordinates of the circumsphere center
    Point3d c1 = {a.x + toCircumsphereCenter.x, a.y + toCircumsphereCenter.y, a.z + toCircumsphereCenter.z};

    double cx = c1.x;
    double cy = c1.y;
    double cz = c1.z;

  // Calculate the radius as the distance from the center to any of the points
  double radius = sqrt(pow(cx - a.x, 2) + pow(cy - a.y, 2) + pow(cz - a.z, 2));

  // Create and return the sphere
  return Sphere{cx, cy, cz, radius};
}

Sphere make_sphere(const Point3d& a, const Point3d& b, const Point3d& c, const Point3d& d) {
    #define U(a,b,c,d,e,f,g,h) (a.z - b.z)*(c.x*d.y - d.x*c.y) - (e.z - f.z)*(g.x*h.y - h.x*g.y)
    #define D(x,y,a,b,c) (a.x*(b.y-c.y) + b.x*(c.y-a.y) + c.x*(a.y-b.y))
    #define E(x,y) ((ra*D(x,y,b,c,d) - rb*D(x,y,c,d,a) + rc*D(x,y,d,a,b) - rd*D(x,y,a,b,c)) / uvw)
    double u = U(a,b,c,d,b,c,d,a);
    double v = U(c,d,a,b,d,a,b,c);
    double w = U(a,c,d,b,b,d,a,c);
    double uvw = 2 * (u + v + w);
    auto sq = [] (Point3d p) { return p.x*p.x + p.y*p.y + p.z*p.z; };
    double ra = sq(a);
    double rb = sq(b);
    double rc = sq(c);
    double rd = sq(d);
    double x0 = E(y,z);
    double y0 = E(z,x);
    double z0 = E(x,y);
    double radius = sqrt(sq({a.x - x0, a.y - y0, a.z - z0}));

    return{x0,y0,z0,radius};
}


// Randomised algorithm to find smallest enclosing sphere from given points
Sphere algo_3d(vector<Point3d> & points, int n){
    random_device rd; 
    mt19937 g(rd()); 
    shuffle(points.begin(), points.end(), g); // Initial random shuffle

    Sphere S = make_sphere(points[0], points[1]);
    
    for(int i = 2; i < n; i++){
        if(!is_inside(S, points[i])){
            S = make_sphere(points[i], points[0]); // if ith point lies outside, we restart with ith point as defining point
            for(int j = 1; j < i; j++){
                if(!is_inside(S, points[j])){
                    S = make_sphere(points[i], points[j]); // jth point lies outside, so we restart with 2 defining points, i and j
                    for(int k = 0; k < j; k++){
                        if(!is_inside(S, points[k])){
                            S = make_sphere(points[i], points[j], points[k]); // kth point lies outside current sphere, so we restart with 3 defining points i, j and k
                            for(int l = 0; l < k; l++){
                                if(!is_inside(S, points[l])){
                                    S = make_sphere(points[i], points[j], points[k], points[l]); // 4 points uniquely define a sphere, so we get a sphere that contains all points till l and has i, j and k on it's boundary 
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return S;

}

int main() {

    cout << "Number of Points = ";

    int n;
    cin >> n;

    cout << "Generate Random Points(y/n): ";
    char c;
    cin >> c;
    Sphere output;
    if(c == 'y'){
        vector<Point3d> points = generateRandom3DPoints(n);
        cout << "Do you want to see the Randomly generated points?(y/n): ";
        cin >> c;
        if(c == 'y'){
            for(auto i : points){
                cout << "(" << i.x << ", " << i.y << ", " << i.z << ")" << endl;
            }
        }
        cout << endl;
        output = algo_3d(points, n);
    }else{
        vector<Point3d> points(n);
        cout << "Please give coordinates of the " << n << " points(in the form of x1 y1 z1 x2 y2 z2 ...)" << endl;
        for(int i = 0; i < n ;i++){
            cin >> points[i].x >> points[i].y >> points[i].z;
        }
        cout << endl;
        output = algo_3d(points, n);
    }

    cout << "Smallest Enclosing Sphere:" << endl;
    cout << "Center = (" << output.x << ", " << output.y << ", " << output.z << ")" << endl;
    cout << "Radius = " << output.r << endl;

    return 0;
}