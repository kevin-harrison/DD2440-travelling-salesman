#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <iomanip>
#include <vector>
#include <string>

using namespace std;

void print_vector(vector<int>& vect){
    for(auto i: vect){
        cout << i << ", ";
    }
}

// Euclidean distance rounded to nearest integer
/*double dist(pair<double, double>& city_a, pair<double, double>& city_b)
{
    return round(sqrt(pow(city_a.first - city_b.first, 2) + pow(city_a.second - city_b.second, 2)));
}*/

double dist(int city_a, int city_b, vector<pair<double, double> >& cities)
{
    return round(sqrt(pow(cities[city_a].first  - cities[city_b].first, 2) +
                      pow(cities[city_a].second - cities[city_b].second, 2)));
}


double total_distance(vector<int>& tour, vector<pair<double, double> >& cities) {
    double total_dist = 0;
    int num_tour = tour.size();
    //cout << "num tour: " << num_tour << endl;
    for(int i = 0; i<num_tour; i++){
        //cout << "total dist: " << total_dist << endl;
        total_dist += dist(tour[i], tour[(i+1)%num_tour], cities);
    }
    //cout << "total dist: " << total_dist << endl;

    return total_dist;
}


// Find naive tour
void getNaiveTour(int num_cities, vector<pair<double, double> >& cities, vector<int>& tour) {
    vector<bool> used;
    used.resize(num_cities, false);
    int best;

    tour[0] = 0;
    used[0] = true;
    for (int i = 1; i < num_cities; i++)
    {
        best = -1;
        for (int j = 0; j < num_cities; j++)
        {
            if (!used[j] && (best == -1 || dist(tour[i - 1], j, cities) < dist(tour[i - 1], best, cities)))
            {
                best = j;
            }
        }
        tour[i] = best;
        used[best] = true;
    }
}


void almostTwoOpt(vector<int>& tour, int num_cities, vector<pair<double, double> >& cities) {
    double new_dist;
    double best_dist = total_distance(tour, cities);
    //cout << "Initial distance: " << best_dist << endl;
    int temp;

    // try to find an improvement
    label: for (int i = 0; i < num_cities-1; i++) {
        for (int j = i+1; j < num_cities; j++) {
            //cout << "-------------------------------------" << endl;
            //cout << i << ", " << j << endl;

            // Check if new tour by swapping 2 edge endpoints is an improvement
            //(a%b+b)%b turning remainder operator into mod operator

            int wrapped_index_below = (((i-1) % num_cities) + num_cities) % num_cities;
            int wrapped_index_above = (j+1) % num_cities;
        
            //cout << "i-1 = " << wrapped_index_below << endl;
            //cout << "j+1 = " << wrapped_index_above << endl;

            if ((j-i) < 2) { // Case 1: i and j are neighbors or neighbors of neighbors ==> two distances to update
                new_dist = best_dist - dist(tour[i], tour[wrapped_index_below], cities)  // -d(i, i-1)
                                        - dist(tour[j], tour[wrapped_index_above], cities)  // -d(j, j+1)
                                        + dist(tour[j], tour[wrapped_index_below], cities)  // +d(j, i-1)
                                        + dist(tour[i], tour[wrapped_index_above], cities); // +d(i, j+1)
                /* cout << "i and j 2-neighbors"<< endl;
                cout << "-dist(i, i-1) " << dist(tour[i], tour[wrapped_index_below], cities) << endl; 
                cout << "-dist(j, j+1) " << dist(tour[j], tour[wrapped_index_above], cities) << endl; 
                cout << "+dist(j, i-1) " << dist(tour[j], tour[wrapped_index_below], cities) << endl; 
                cout << "+dist(i, j+1) " << dist(tour[i], tour[wrapped_index_above], cities) << endl; */
            }
            else if ((i - (j - num_cities)) < 2) { // Case 2: i and j are neighbors or neighbors of neighbors because of wrapping indices ==> 2 distances to update
                new_dist = best_dist - dist(tour[i], tour[i+1], cities)  // -d(i, i+1)
                                        - dist(tour[j], tour[j-1], cities)  // -d(j, j-1)
                                        + dist(tour[j], tour[i+1], cities)  // +d(j, i+1)
                                        + dist(tour[i], tour[j-1], cities); // +d(i, j-1)
                /* cout << "i and j wrapped-2-neighbors"<< endl;
                cout << "-dist(i, i+1) " << dist(tour[i], tour[i+1], cities) << endl; 
                cout << "-dist(j, j-1) " << dist(tour[j], tour[j-1], cities) << endl; 
                cout << "+dist(j, i+1) " << dist(tour[j], tour[i+1], cities) << endl; 
                cout << "+dist(i, j-1) " << dist(tour[i], tour[j-1], cities) << endl; */
            }
            else { // Case 3: four distances to update
                new_dist = best_dist - dist(tour[wrapped_index_below], tour[i], cities) // i-1 and i
                                        - dist(tour[i], tour[i+1], cities)
                                        - dist(tour[j-1], tour[j], cities)
                                        - dist(tour[j], tour[wrapped_index_above], cities)
                                        + dist(tour[wrapped_index_below], tour[j], cities) // i-1 and j
                                        + dist(tour[j], tour[i+1], cities) // j and i+1
                                        + dist(tour[j-1], tour[i], cities) // j-1 and i
                                        + dist(tour[i], tour[wrapped_index_above], cities); // i and j+1
                //cout << "i and j aren't neighbors" << endl;
                //cout << "dist(i-1, i) "  << dist(tour[wrapped_index_below], tour[i], cities) << endl; // i-1 and i 
                //cout << "dist(i, i+1) "  << dist(tour[i], tour[i+1], cities) << endl;
                //cout << "dist(j-1, j) "  << dist(tour[j-1], tour[j], cities) << endl;
                //cout << "dist(j, j+1) "  << dist(tour[j], tour[wrapped_index_above], cities) << endl;
                //cout << "dist(i-1, j) "  << dist(tour[wrapped_index_below], tour[j], cities) << endl; // i-1 and j
                //cout << "dist(j, i+1) "  << dist(tour[j], tour[i+1], cities) << endl; // j and i+1
                //cout << "dist(j-1, i) "  << dist(tour[j-1], tour[i], cities) << endl; // j-1 and i
                //cout << "dist(i, j+1) "  << dist(tour[i], tour[wrapped_index_above], cities) << endl; // i and j+1
            }

            //cout << "tour: ";
            //print_vector(tour);
            //cout << "dist=" << best_dist << endl;
            //cout << "new_dist=" << new_dist << endl << endl;

            //return;

            if (new_dist < best_dist) {
                // Create new tour by swapping 2 edge endpoints
                //cout << "FOUND AN IMPROVEMENT" << endl;
                temp = tour[i];
                tour[i] = tour[j];
                tour[j] = temp;  
                best_dist = new_dist;
                //goto label;
            }
        }
    }
}

void swap_tour(vector<int>& tour, int start, int end){
    int temp;

    if (start > end) {
        temp = start;
        start = end;
        end = temp;
    }

    for(int i= 0; i < floor((end - start) / 2); i++) {
        // swap tour[start+i] and tour[end-i]
        temp = tour[start+i];
        tour[start+i] = tour[end-i];
        tour[end-i] = temp;
    }
}

void twoOpt(vector<int>& tour, int num_cities, vector<pair<double, double> >& cities) {
    double new_dist;
    double best_dist = total_distance(tour, cities);
    //cout << "Initial distance: " << best_dist << endl;
    int temp;

    // try to find an improvement
    for (int i = 0; i < num_cities-1; i++) {
        for (int j = i+1; j < num_cities; j++) {

            /* c1 = i   c2 = i+1           c3=j     c4=j+1 */
            // d(c1,c2) + d(c3,c4) > d(c2, c3) + d(c1, c4)
            double edge1 = dist(tour[i], tour[i+1], cities);
            double edge2 = dist(tour[j], tour[(j+1) % num_cities], cities);
            double edge3 = dist(tour[i+1], tour[j], cities);
            double edge4 = dist(tour[i], tour[(j+1) % num_cities], cities);
            if (edge1 + edge2 >  edge3 + edge4) {
                cout << "FOUND AN IMPROVEMENT" << endl;
                swap_tour(tour, i+1, (j+1) % num_cities);
            }
        }
    }
}


int main()
{
    
    std::cout << std::fixed;
    std::cout << std::setprecision(10);
    
    int num_cities;
    if (cin >> num_cities)
    {
        // Get input
        vector<pair<double, double> > cities;
        cities.reserve(num_cities);
        for (int i = 0; i < num_cities; i++)
        {
            double xcord, ycord;
            cin >> xcord >> ycord;
            cities[i] = make_pair(xcord, ycord);
        }

        // Debug print what has been input from stdin
        /* cout << num_cities << endl;
        for (int i = 0; i < num_cities; i++)
        {
            cout << to_string(cities[i].first) + " " + to_string(cities[i].second) + "\n";
        } */

        // Get initial tour
        vector<int> initial_tour;
        initial_tour.resize(num_cities);
        getNaiveTour(num_cities, cities, initial_tour);
        double naive_dist = total_distance(initial_tour, cities);

        // Improve tour with heuristic
        twoOpt(initial_tour, num_cities, cities);
        double new_dist = total_distance(initial_tour, cities);

        // Print answer
        for (int i = 0; i < num_cities; i++)
        {
            cout << initial_tour[i] << endl;
        }

        cout << naive_dist << endl;
        cout << new_dist << endl;
    }
}