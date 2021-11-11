#include <iostream>  
#include <stdlib.h>
#include <cmath>
#include <iomanip>

using namespace std;

// Euclidean distance rounded to nearest integer
double dist(double x1, double y1, double x2, double y2) {
    return round(sqrt(pow(x2-x1, 2) + pow(y2-y1, 2)));
}

int main() {
    // Read input
    std::cout << std::fixed;
    std::cout << std::setprecision(10);
    
    while(true) {
        int numCities;
        if (cin >> numCities) {

            double * cities_x = new double[numCities];
            double * cities_y = new double[numCities];
            
            for (int i = 0; i < numCities; i++) {
                cin >> cities_x[i];
                cin >> cities_y[i];
            }

            cout << numCities << endl; 
            for (int i = 0; i < numCities; i++){
                cout << cities_x[i] << " " << cities_y[i] << endl;
            } 

            // Find naive tour
            int * tour = new int[numCities];
            bool * used = new bool[numCities];
            int best;

            tour[0] = 0;
            used[0] = true;
            for (int i = 1; i < numCities; i++) {
                best = -1;
                for (int j = 0; j < numCities; j++) { 
                    if (!used[j] && (best == -1 || dist(cities_x[tour[i-1]], cities_y[tour[i-1]], cities_x[j], cities_y[j]) < 
                                                dist(cities_x[tour[i-1]], cities_y[tour[i-1]], cities_x[best], cities_y[best]))) {
                        best = j;
                    }
                }
                tour[i] = best;
                used[best] = true;
            }

            for (int i = 0; i < numCities; i++) {
                cout << tour[i] << endl;
            }

            delete[] cities_x;
            delete[] cities_y;
            delete[] tour;
            delete[] used;
        }
    }
}

/* 
tour[0] = 0
used[0] = true
for i = 1 to n-1
    best = -1 // index of best city to go to next?
    for j = 0 to n-1
        if not used[j] and (best = -1 or dist(tour[i-1], j) < dist(tour[i-1], best))
        best = j
    tour[i] = best
    used[best] = true
return tour 
*/