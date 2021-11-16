#include <iostream>
#include <fstream>
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
    cout << endl;
}


// Rounded Euclidean distance
double dist(int city_a, int city_b, vector<pair<double, double> >& cities)
{
    return round(sqrt(pow(cities[city_a].first  - cities[city_b].first, 2) +
                      pow(cities[city_a].second - cities[city_b].second, 2)));
}


// Calculate the distance of an entire tour
double total_distance(vector<int>& tour, vector<pair<double, double> >& cities) {
    double total_dist = 0;
    int num_tour = tour.size();

    for(int i = 0; i<num_tour; i++) {
        total_dist += dist(tour[i], tour[(i+1)%num_tour], cities);
    }
    return total_dist;
}


// Save tour to file filename in a .dot file format
void save_tour(vector<int>& tour, vector<pair<double, double> >& cities, string filename) {
    ofstream myfile;
    myfile.open (filename);
    myfile << "strict digraph {\nlabel=\"" << filename << "\";" << endl;
    myfile << "node [shape=plaintext, fontcolor=red, fixedsize=true, width=0.05];" << endl;

    // Add nodes
    for (int i = 0; i < tour.size(); i++) {
        myfile << i << " [pos=\"" << cities[i].first << "," << cities[i].second << "!\"];" << endl;
    }

    // Add edges
    for (int i = 0; i < tour.size(); i++) {
        int from = tour[i];
        int to = tour[(i+1) % tour.size()];
        myfile << from << " -> " << to << " [label=\"" << dist(from, to, cities) << "\"];" << endl;
    }

    myfile << "}\n";
    myfile.close();
    return;
}


// Find naive tour
void getNaiveTour(int num_cities, vector<pair<double, double> >& cities, vector<int>& tour) {
    vector<bool> used;
    used.resize(num_cities, false);
    int best;

    tour[0] = 0;
    used[0] = true;
    for (int i = 1; i < num_cities; i++) {
        best = -1;
        for (int j = 0; j < num_cities; j++) {
            if (!used[j] && (best == -1 || dist(tour[i - 1], j, cities) < dist(tour[i - 1], best, cities))) {
                best = j;
            }
        }
        tour[i] = best;
        used[best] = true;
    }
}


// Reverse elements in a tour from index "start" to index "end" (inclusive)
void swap_tour(vector<int>& tour, int start, int end){
    int temp;

    // Make sure bounds are correct
    if (start > end) {
        temp = start;
        start = end;
        end = temp;
    }

    // Swap elements pair by pair between [start,...,end]
    for(int i= 0; i < floor((end - start + 1) / 2); i++) {
        // swap tour[start+i] and tour[end-i]
        temp = tour[start+i];
        tour[start+i] = tour[end-i];
        tour[end-i] = temp;
    }
}



void twoOpt(vector<int>& tour, int num_cities, vector<pair<double, double> >& cities) {   
    int improvements = 1;
    int counter = 0;

    
    label: if (counter > 30 || improvements == 0) return;
    improvements = 0;

    // Iterate over all cities to find swapping improvements
    for (int i = 0; i < num_cities-1; i++) {
        for (int j = i+1; j < num_cities; j++) {
            if (j == i+1 || i == ((j+1) % num_cities)) continue; // Nodes are neighbors so swap is meaningless

            double edge1 = dist(tour[i],   tour[i+1], cities);                  // dist(c1,c2)
            double edge2 = dist(tour[j],   tour[(j+1) % num_cities], cities);   // dist(c3,c4)
            double edge3 = dist(tour[i],   tour[j], cities);                    // dist(c1,c3)
            double edge4 = dist(tour[i+1], tour[(j+1) % num_cities], cities);   // dist(c2,c4)
            if (edge1 + edge2 >  edge3 + edge4) { // Only swap if it will result in a shorter tour
                swap_tour(tour, i+1, j);
                improvements += 1;
                //save_tour(tour, cities, "./graphs/tour" + to_string(counter) + ".dot");
            } 
        }
    }
    counter += 1;
    goto label;
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
        //double naive_dist = total_distance(initial_tour, cities);
        //save_tour(initial_tour, cities, "./graphs/tour_initial.dot");

        // Improve tour with heuristic
        twoOpt(initial_tour, num_cities, cities);
        //double new_dist = total_distance(initial_tour, cities);
        //save_tour(initial_tour, cities, "./graphs/tour_final.dot");

        // Print answer
        for (int i = 0; i < num_cities; i++)
        {
            cout << initial_tour[i] << endl;
        }

        //cout << naive_dist << endl;
        //cout << new_dist << endl;
    }
}