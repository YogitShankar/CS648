#include <bits/stdc++.h>
#include <vector>
#include <fstream>
#include <chrono>
#include <random>
using namespace std;

// Helper functions and structs
struct Point {
    double x;
    double y;
};

struct Circle {
    double x;
    double y;
    double r;
};

// Function to generate random points following a Gaussian distribution
std::vector<Point> generateRandomPoints(int n) {
    static std::default_random_engine generator;
    std::normal_distribution<double> distributionX(0, 10000);
    std::normal_distribution<double> distributionY(0, 10000);

    std::vector<Point> points;
    for (int i = 0; i < n; ++i) {
        double x = distributionX(generator);
        double y = distributionY(generator);
        points.push_back({x, y});
    }
    return points;
}

bool is_inside(const Point& point, const Circle& circle){
    return (point.x - circle.x) * (point.x - circle.x) + (point.y - circle.y) * (point.y - circle.y) <= circle.r * circle.r + 1e-8;
}

// make_circle functions return smalles Circle that passes through given points
Circle make_circle(const Point& a,const Point& b){
    Circle out;
    out.x = (a.x+b.x)/2;
    out.y = (a.y+b.y)/2;
    out.r = sqrtl((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y))/2;
    return out;
}

Circle make_circle(const Point& p1, const Point& p2, const Point& p3) {
    double offset = pow(p2.x, 2) + pow(p2.y, 2);
    double bc = (pow(p1.x, 2) + pow(p1.y, 2) - offset) / 2.0;
    double cd = (offset - pow(p3.x, 2) - pow(p3.y, 2)) / 2.0;
    double det = (p1.x - p2.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p2.y);

    double idet = 1 / det;

    double cx = (bc * (p2.y - p3.y) - cd * (p1.y - p2.y)) * idet;
    double cy = (cd * (p1.x - p2.x) - bc * (p2.x - p3.x)) * idet;
    double r = sqrt(pow(p2.x - cx, 2) + pow(p2.y - cy, 2));

    return {cx, cy, r};
}


Circle algo_n4(vector<Point> & points, int n){
    Circle C = {0,0,1e15};
    for(int i = 0; i < n; i++){ // Iterates over all nC2 points and checks if the circle with them on diameter is SEC or not
        for(int j = i+1; j < n; j++){
            Circle C1 = make_circle(points[i], points[j]);
            int flag = 1; // Check if all points lie inside or not
            for(auto p : points){
                if(!is_inside(p, C1)) flag = 0;
            }
            if(flag && C1.r < C.r){
                C = C1;
            }
        }
    }
    for(int i = 0; i < n; i++){ // Iterates over all nC3 points and checks if circle passing through them is SEC or not
        for(int j = i+1; j < n; j++){
            for(int k = j+1; k < n; k++){
                Circle C1 = make_circle(points[i], points[j], points[k]);
                int flag = 1; // Check if all points lie inside or not
                for(auto p : points){
                    if(!is_inside(p, C1)) flag = 0;
                }
                if(flag && C1.r < C.r){
                    C = C1;
                }
            }
        }
    }
    return C;
}

int main(){
    cout << "Number of Points = ";

    int n;
    cin >> n;

    cout << "Generate Random Points(y/n): ";
    char c;
    cin >> c;
    Circle output;
    if(c == 'y'){
        vector<Point> points = generateRandomPoints(n);
        cout << "Do you want to see the Randomly generated points?(y/n): ";
        cin >> c;
        if(c == 'y'){
            for(auto i : points){
                cout << "(" << i.x << ", " << i.y << ")" << endl;
            }
        }
        cout << endl;
        output = algo_n4(points, n);
    }else{
        vector<Point> points(n);
        cout << "Please give coordinates of the " << n << " points(in the form of x1 y1 x2 y2 ...)" << endl;
        for(int i = 0; i < n ;i++){
            cin >> points[i].x >> points[i].y;
        }
        cout << endl;
        output = algo_n4(points, n);
    }

    cout << "Smallest Enclosing Cirle:" << endl;
    cout << "Center = (" << output.x << ", " << output.y << ")" << endl;
    cout << "Radius = " << output.r << endl;

    return 0;
}