#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <iomanip>
#include <vector>

using namespace std;

// Euclidean distance rounded to nearest integer
double dist(pair<double, double> city_a, pair<double, double> city_b)
{
    return round(sqrt(pow(city_a.first - city_b.first, 2) + pow(city_a.second - city_b.second, 2)));
}

int main()
{
    // Read input
    std::cout << std::fixed;
    std::cout << std::setprecision(10);
    
    int num_cities;
    if (cin >> num_cities)
    {
        pair<double, double> cities[num_cities];

        for (int i = 0; i < num_cities; i++)
        {
            double xcord, ycord;
            cin >> xcord >> ycord;
            cities[i] = make_pair(xcord, ycord);
        }

        // Debug print what has been input from stdin
        cout << num_cities << endl;
        for (int i = 0; i < num_cities; i++)
        {
            cout << to_string(cities[i].first) + " " + to_string(cities[i].second) + "\n";
        }

        // Find naive tour
        vector<int> tour;
        vector<bool> used;
        tour.resize(num_cities, -1);
        used.resize(num_cities, false);
        int best;

        tour[0] = 0;
        used[0] = true;
        for (int i = 1; i < num_cities; i++)
        {
            best = -1;
            for (int j = 0; j < num_cities; j++)
            {
                if (!used[j] && (best == -1 || dist(cities[tour[i - 1]], cities[j]) <
                                                   dist(cities[tour[i - 1]], cities[best])))
                {
                    best = j;
                }
            }
            tour[i] = best;
            used[best] = true;
        }

        for (int i = 0; i < num_cities; i++)
        {
            cout << tour[i] << endl;
        }
    }
}