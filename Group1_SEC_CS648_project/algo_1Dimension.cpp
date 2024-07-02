#include <bits/stdc++.h>
#include <vector>
#include <fstream>
#include <chrono>
#include <random>
using namespace std;

// Helper Functions and Structs
struct Line { // This line segment is defined as a center x and distance r, the line subtends from x-r to x+r
    double x;
    double r;
};

struct Point1d{
    double x;
};

std::vector<Point1d> generateRandom1DPoints(int n) {
    static std::default_random_engine generator;
    std::normal_distribution<double> distributionX(0, 10000);

    std::vector<Point1d> points;
    for (int i = 0; i < n; ++i) {
        double x = distributionX(generator);
        points.push_back({x});
    }
    return points;
}


// Tells if point lies inside the given line segment
bool is_inside(const Line& L, const Point1d& point) {
    return abs(point.x - L.x) <= L.r + 1e-8;
}

// Return the line segment with p1 and p2 as end points
Line make_line(const Point1d& p1, const Point1d& p2){
    return {(p1.x + p2.x)/2, abs(p1.x - p2.x)/2};
}


Line algo_1d(vector<Point1d> & points, int n){
    random_device rd; 
    mt19937 g(rd()); 
    shuffle(points.begin(), points.end(), g); // Initial Random Shuffle

    Line L = make_line(points[0], points[1]);
    for(int i = 2; i < n; i++){
        if(!is_inside(L, points[i])){
            L = make_line(points[0], points[i]); // ith point is not a part of current line segment hence we restart with ith point as defining point
            for(int j = 1; j < i; j++){
                if(!is_inside(L, points[j])){
                    L = make_line(points[i], points[j]); // 2 points define the line uniquely
                }
            }
        }
    }

    return L;

}

int main() {

    cout << "Number of Points = ";

    int n;
    cin >> n;

    cout << "Generate Random Points(y/n): ";
    char c;
    cin >> c;
    Line output;
    if(c == 'y'){
        vector<Point1d> points = generateRandom1DPoints(n);
        cout << "Do you want to see the Randomly generated points?(y/n): ";
        cin >> c;
        if(c == 'y'){
            for(auto i : points){
                cout << i.x << endl;
            }
        }
        cout << endl;
        output = algo_1d(points, n);
    }else{
        vector<Point1d> points(n);
        cout << "Please give coordinates of the " << n << " points(in the form of x1 x2 x3 ...)" << endl;
        for(int i = 0; i < n ;i++){
            cin >> points[i].x;
        }
        cout << endl;
        output = algo_1d(points, n);
    }

    cout << "Smallest Enclosing Line Segment:" << endl;
    cout << "Center = " << output.x << endl;
    cout << "Radius = " << output.r << endl;

    return 0;
}