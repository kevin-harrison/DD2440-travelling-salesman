#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <ctime>
//#include <iomanip>
#include <vector>
#include <tuple>
#include <algorithm>
#include <set>
#include <stack>

using std::tuple;
using namespace std;

void print_vector(vector<int> &vect) {
    for (auto i : vect) {
        cout << i << ", ";
    }
    cout << endl;
}

// Rounded Euclidean distance
int dist(int city_a, int city_b, vector<pair<double, double>> &cities) {
    return round(sqrt(pow(cities[city_a].first - cities[city_b].first, 2) +
                      pow(cities[city_a].second - cities[city_b].second, 2)));
}


// Calculate the distance of an entire tour
int total_distance(vector<int> &tour, int **&dist_matrix) {
    int total_dist = 0;
    for (int i = 0; i < tour.size() - 1; i++) 
        total_dist += dist_matrix[tour[i]][tour[i + 1]];
    return total_dist + dist_matrix[tour[tour.size() - 1]][tour[0]]; // Add distance linking back to cycle start
}


// Save tour to file filename in a .dot file format
void save_tour(vector<int> &tour, vector<pair<double, double>> &cities, string filename) {
    ofstream myfile;
    myfile.open(filename);

    // Start graph and define node style
    myfile << "strict digraph {\nlabel=\"" << filename << "\";" << endl;
    myfile << "node [shape=plaintext, fontcolor=red, fixedsize=true, width=0.05];" << endl;

    // Add nodes
    for (int i = 0; i < tour.size(); i++)
        myfile << i << " [pos=\"" << cities[i].first << "," << cities[i].second << "!\"];" << endl;

    // Add edges
    for (int i = 0; i < tour.size(); i++) {
        int from = tour[i];
        int to = tour[(i + 1) % tour.size()];
        myfile << from << " -> " << to << " [label=\"" << dist(from, to, cities) << "\"];" << endl;
    }

    myfile << "}\n";
    myfile.close();
    return;
}


// Find naive tour
void getNaiveTour(vector<int> &tour, int num_cities, int **&dist_matrix, int start) {
    vector<bool> used;
    used.resize(num_cities, false);
    int best;

    tour[0] = start;
    used[start] = true;
    for (int i = 1; i < num_cities; i++) {
        best = -1;
        for (int j = 0; j < num_cities; j++) {
            if (!used[j] && (best == -1 || dist_matrix[tour[i-1]][j] < dist_matrix[tour[i-1]][best]))
                best = j;
        }
        tour[i] = best;
        used[best] = true;
    }
}

// Reverse elements in a tour from index "start" to index "end" (inclusive)
void twoOptSwapA(vector<int> &tour, int start, int end) {
    int temp;
    // Make sure bounds are correct
    if (start > end) { // TODO: Can remove this because this case never gets called i think?
        temp = start;
        start = end;
        end = temp;
    }

    // Swap elements pair by pair between [start,...,end]
    for (int i = 0; i < floor((end - start + 1) / 2); i++) {
        // swap tour[start+i] and tour[end-i]
        temp = tour[start + i];
        tour[start + i] = tour[end - i];
        tour[end - i] = temp;
    }
}

void twoOptSwapB(vector<int> &tour, int start, int end) {
    int temp = tour[start+1];
    for(int m = start+1; m < end; m++)
        tour[m] = tour[m+1];
    tour[end] = temp;
}

void twoOptSwapC(vector<int> &tour, int start, int end) {
    int temp = tour[end];
    //for (int m = start+2; m <= end; m++)
    for (int m = end; m >= start+2; m--)
        tour[m] = tour[m-1];
    tour[start + 1] = temp;
}

int getCostDiff(vector<int> &tour, int **& distMatrix, int i, int j, int num_cities) {
    int newEdges = distMatrix[tour[i]][tour[j]] + distMatrix[tour[i + 1]][tour[(j + 1) % num_cities]];
    int prevEdges = distMatrix[tour[i]][tour[i + 1]] + distMatrix[tour[j]][tour[(j + 1) % num_cities]];
    return newEdges - prevEdges;
}

int getCostDiffB(vector<int> &tour, int **& distMatrix, int i, int j, int num_cities) {
    int newEdges = distMatrix[tour[i]][tour[i+2]] + distMatrix[tour[j]][tour[i+1]] + distMatrix[tour[i + 1]][tour[(j + 1) % num_cities]];
    int prevEdges = distMatrix[tour[i]][tour[i + 1]] + distMatrix[tour[i+1]][tour[i + 2]] + distMatrix[tour[j]][tour[(j + 1) % num_cities]];
    return newEdges - prevEdges;
}

int getCostDiffC(vector<int> &tour, int **& distMatrix, int i, int j, int num_cities) {
    int newEdges = distMatrix[tour[j-1]][tour[(j + 1) % num_cities]] + distMatrix[tour[i]][tour[j]] + distMatrix[tour[j]][tour[i + 1]];
    int prevEdges = distMatrix[tour[i]][tour[i+1]]+ distMatrix[tour[j]][tour[(j + 1) % num_cities]] +  distMatrix[tour[j]][tour[j - 1]] ;
    return newEdges - prevEdges;
}

void twoOpt(vector<int> &tour, int num_cities, int **&distMatrix) {
    int counter = 0, cost = 0;
    int best_cost = 0, best_swap_i = 0, best_swap_j = 0;

    while (counter < 600) {
        best_cost = 0;
        // Iterate over all cities to find swapping improvements
        for (int i = 0; i < num_cities - 1; i++) {
            for (int j = i + 2; j < num_cities; j++) {
                if (i == ((j + 1) % num_cities)) continue; // Nodes are neighbors so swap is meaningless

                // Save swap arguments if it will result in the shortest tour we've seen so far
                cost = getCostDiff(tour, distMatrix, i, j, num_cities);
                if (cost < best_cost) { 
                    best_cost = cost;
                    best_swap_i = i + 1;
                    best_swap_j = j;
                }
            }
        }
        if (best_cost == 0) break; // Didn't find improving 2opt move
        twoOptSwapA(tour, best_swap_i, best_swap_j);
        counter += 1;
    }
}

void twoPointFiveOpt(vector<int> &tour, int num_cities, int **&distMatrix) {
    int improvements = 1;
    int counter = 0;
    int minCost = 1;
    int iMin = 0;
    int jMin = 0;
    int B = 0;
    int C = 0;
    bool Bbest = false;
    bool Cbest = false;


    label:
    if (counter > 600 || minCost == 0) return;
    minCost = 0;
    Bbest = false;
    Cbest = false;
    // Iterate over all cities to find swapping improvements

    for (int i = 0; i < num_cities - 3; i++) {
        for (int j = i + 3; j < num_cities; j++) { // j = i+2 is a 2-opt move
            if (i == ((j + 1) % num_cities)) continue; // We think this is just a 2opt move.
            
            B = getCostDiffB(tour, distMatrix, i, j, num_cities);
            C = getCostDiffC(tour, distMatrix, i, j, num_cities);

            if (minCost > min(B,C)) { // Only swap if it will result in a shorter tour
                minCost = min(B,C);
                iMin = i;
                jMin = j;
                if(B < C) {
                    Bbest = true;
                    Cbest = false;
                } else {
                    Bbest = false;
                    Cbest = true;
                }
                //save_tour(tour, cities, "./graphs/tour" + to_string(counter) + ".dot");
            }
        }
    }
    if (Bbest) {
        //print_vector(tour);
        twoOptSwapB(tour, iMin, jMin);
        //print_vector(tour);
        //cout << "Switched " << iMin << " to " << jMin << " with swap B" << endl;
        //cout << endl;
    }else if(Cbest) {
        //print_vector(tour);
        twoOptSwapC(tour, iMin, jMin);
        //print_vector(tour);
        //cout << "Switched " << iMin << " to " << jMin << " with swap C" << endl;
        //cout << endl;
    }else 
    counter += 1;
    goto label;
}


bool sortByDist(const tuple<int, int, int> &a, const tuple<int, int, int> &b) {
    return (get<2>(a) < get<2>(b));
}

int find_cluster(int i, int *vertexClusters) {
    if (i == vertexClusters[i]) {
        return i;
    }
    return find_cluster(vertexClusters[i], vertexClusters);
}

bool contains(vector<int> v, int x) {
    if (std::find(v.begin(), v.end(), x) != v.end()) {
        return true;
    }
    return false;
}

void printGraph(vector<vector<int>> v) {
    int i = 0;
    vector<vector<int>>::iterator row;
    vector<int>::iterator col;
    for (row = v.begin(); row != v.end(); row++) {
        cout << "VERTEX " << i << ": ";
        for (col = row->begin(); col != row->end(); col++) {
            cout << (*col) << " ";
        }
        cout << endl;
        i += 1;
    }
}

void printGraph(vector<int> *&v) {
    for (size_t i = 0; i < v->size(); i++) {
        cout << "VERTEX " << i << ": ";
        for (int j = 0; j < v[i].size(); j++) {
            cout << v[i][j] << " ";
        }
        cout << endl;
    }
}

/* returns a vector of vectors of ints. This represents a list of vertices, where each of them has a list of edges that attach them to the MST */
vector<int> christofides(vector<pair<double, double>> &cities, int num_cities, vector<tuple<int, int, int>> edges)
{
    // -------------------------------Kruskals-----------------------------------
    
    vector<vector<int>> mst{};
    mst.reserve(num_cities);
    int *vertexClusters = new int[num_cities];
    for (int index = 0; index < num_cities; index++)
    {
        vertexClusters[index] = index;

        vector<int> innerVector;
        mst.push_back(innerVector);
        //cout << "insert " << index << endl;
        //set<int> cluster;
        //cluster.insert(index);
        //vertexClusters.insert(cluster);
    }

    /*
    // Print vertexClusters
    for (set<set<int> >::iterator clusterIt = vertexClusters.begin(); clusterIt != vertexClusters.end(); ++clusterIt) {
        for (set<int>::iterator vertexIt = (*clusterIt).begin(); vertexIt != (*clusterIt).end(); ++vertexIt) {
            cout << *vertexIt << " ";
        }
        cout << endl;
    }*/

    /*     vector<vector<int> >::iterator row = mst.begin();
    vector<vector<int> >::iterator col = row.begin(); */
    bool merge = false;
    for (tuple<int, int, int> edge : edges)
    {

        int cluster1 = find_cluster(get<0>(edge), vertexClusters);
        int cluster2 = find_cluster(get<1>(edge), vertexClusters);

        if (cluster1 != cluster2)
        {
            // Insert into MST
            //cout << "type: " << typeid(get<0>(edge)).name() << endl;
            //cout << get<0>(mst);

            mst[get<0>(edge)].push_back(get<1>(edge));
            mst[get<1>(edge)].push_back(get<0>(edge));

            vertexClusters[cluster1] = vertexClusters[cluster2]; // Merge

            /* for (int i = 0; i < num_cities; i++)
            {
                cout << "VERTEX " << i << ": " << vertexClusters[i] << endl;
            } */
        }

        /* merge = checkAndMergeClusters(vertexClusters, get<0>(edge), get<1>(edge));
        cout << "edge: (" << get<0>(edge) << ", " << get<1>(edge) << ", " << get<2>(edge) << ")" << endl;
        if (merge) {
            
             it = get<get<0>(edge)>(mst).begin();
            get<(get<0>(edge))>(mst).insert(it, get<1>(edge));
            it = get<get<1>(edge)>(mst).begin();
            get<get<1>(edge)>(mst).insert(it, get<0>(edge));  
        } */
    }

    // ------------------------------ MST found ---------------------------------------------

    /* cout << "--------------------------" << endl;
    cout << "MST: " << endl;
    printGraph(mst); */

    // -------------------------------Find minimal matching-----------------------------------

    // ha lista med odd degree vertices
    vector<int> oddVertices;
    oddVertices.reserve(num_cities); // Could be improved
    for (int i = 0; i < mst.size(); i++)
    {
        if (mst[i].size() % 2 != 0)
        { // Odd degree
            oddVertices.push_back(i);
        }
    }

    // gå igenom edges kortast -> längst
    for (size_t i = 0; i < edges.size(); i++)
    {
        if (contains(oddVertices, get<0>(edges[i])) && contains(oddVertices, get<1>(edges[i])))
        {
            // lägg till den i mst
            mst[get<0>(edges[i])].push_back(get<1>(edges[i]));
            mst[get<1>(edges[i])].push_back(get<0>(edges[i]));

            // ta bort elementen ur listan med odd degree vertices
            oddVertices.erase(find(oddVertices.begin(), oddVertices.end(), get<0>(edges[i])));
            oddVertices.erase(find(oddVertices.begin(), oddVertices.end(), get<1>(edges[i])));
        }
    }

    /* cout << "--------------------------" << endl;
    cout << "Minimal matching: " << endl;
    printGraph(mst); */

    // ------------------------Find Euler circuit---------------------

    // Copy mst
    vector<int> *mstCopy = new vector<int>[num_cities];
    //mstCopy.reserve(num_cities);
    for (int i = 0; i < num_cities; i++)
    {
        //mstCopy[i].reserve(mst[i].size());
        mstCopy[i].resize(mst[i].size());
        mstCopy[i] = mst[i];
        /*for (int j = 0; j < mst[i].size(); j++) {
            mstCopy[i].push_back(mst[i][j]);
        } */
    }

    int start = 0;
    stack<int> stack;
    int pos = start;
    vector<int> path;
    path.reserve(2 * num_cities);
    path.push_back(start);
    while (!stack.empty() || mstCopy[pos].size() > 0)
    {
        //Current node has no neighbors
        if (mstCopy[pos].empty())
        {
            //add to circuit
            path.push_back(pos);
            //remove last vertex from stack and set it to current
            pos = stack.top();
            stack.pop();
        }
        //If current node has neighbors
        else
        {
            //Add vertex to stack
            stack.push(pos);
            //Take a neighbor
            int neighbor = mstCopy[pos].back();
            //Remove edge between neighbor and current vertex
            mstCopy[pos].pop_back();
            for (int i = 0; i < mstCopy[neighbor].size(); i++)
            {
                if (mstCopy[neighbor][i] == pos)
                {
                    mstCopy[neighbor].erase(mstCopy[neighbor].begin() + i);
                }
            }
            //Set neighbor as current vertex
            pos = neighbor;
        }
    }
    path.push_back(pos);

    /* cout << "--------------------------" << endl;
    cout << "Euler cycle: " << endl;
    for (int i = 0; i < path.size(); i++) {
        cout << path[i] << " ";
    }
    cout << endl; */

    //---------------------- Euler circuit without duplicates (Hamiltonian) -----------------------

    bool *visited = new bool[num_cities];
    for (int i = 0; i < num_cities; i++)
    {
        visited[i] = 0;
    }

    vector<int> tour;
    tour.reserve(num_cities);

    for (int i = 0; i < path.size(); i++)
    {
        if (!visited[path[i]])
        {
            tour.push_back(path[i]);
            visited[path[i]] = true;
        }
    }

    delete[] vertexClusters;
    delete[] mstCopy;
    delete[] visited;
    return tour;
}

void prepareDistMarixAndEdges(int **&distMatrix, vector<tuple<int, int, int> > &edges, int num_cities, vector<pair<double, double> > &cities)
{
    // Reserve edges
    int to_reserve = 0;
    // Calculate distances, populate distance matrix and edges
    for (int i = 0; i < num_cities; i++) {
        distMatrix[i] = new int[num_cities];
        to_reserve += i;
    }
    edges.reserve(to_reserve);


    tuple<int, int, int> edge;
    for (int i = 0; i < num_cities; i++) {
        for (int j = i+1; j < num_cities; j++) {
            int distance = dist(i, j, cities);
            distMatrix[i][j] = distance;
            distMatrix[j][i] = distance;
            edge = make_tuple(i, j, distance);
            edges.push_back(edge);
        }
    }
    sort(edges.begin(), edges.end(), sortByDist);
}

int main() {
    //std::cout << std::fixed;
    //std::cout << std::setprecision(10);

    clock_t program_start = clock();

    int num_cities;
    if (cin >> num_cities) {
        // Get input
        vector<pair<double, double>> cities;
        cities.reserve(num_cities);
        for (int i = 0; i < num_cities; i++) {
            double xcord, ycord;
            cin >> xcord >> ycord;
            cities[i] = make_pair(xcord, ycord);
        }

        // Create datastructures
        int **distMatrix;
        distMatrix = new int *[num_cities];
        vector<tuple<int, int, int>> edges;
        prepareDistMarixAndEdges(distMatrix, edges, num_cities, cities);

        // try out christofides
        //vector<int> christ_tour = christofides


        /*clock_t iteration_start = clock();
        vector<int> naive_tour;
        naive_tour.resize(num_cities);
        getNaiveTour(naive_tour, num_cities, distMatrix, 0);
        twoOpt(naive_tour, num_cities, distMatrix);
        twoPointFiveOpt(naive_tour, num_cities, distMatrix);
        int best_dist = total_distance(naive_tour, distMatrix);
        vector<int> best_tour = naive_tour;
        clock_t iteration_end = clock();*/



        // Initial run
        clock_t iteration_start = clock();
        vector<int> naive_tour;
        naive_tour.resize(num_cities);
        getNaiveTour(naive_tour, num_cities, distMatrix, 0);
        twoOpt(naive_tour, num_cities, distMatrix);
        twoPointFiveOpt(naive_tour, num_cities, distMatrix);
        int best_dist = total_distance(naive_tour, distMatrix);
        vector<int> best_tour = naive_tour;
        clock_t iteration_end = clock();
        
        // Gather time info
        double approx_iteration_time = 1.3 * 1000.0 * (iteration_end-iteration_start) / CLOCKS_PER_SEC;
        double time_left = 2000 - (1000.0 * (iteration_start-program_start) / CLOCKS_PER_SEC);
        int iterations_left = (int) (time_left / approx_iteration_time);

        /*cout << "init dist: " << best_dist << endl;
        cout << "time left: " << time_left << "ms" << endl;
        cout << "iter time: " << approx_iteration_time << "ms" << endl;
        cout << "iterations: " << iterations_left << endl;*/


        for (int i = 1; i < min(iterations_left, num_cities); i ++) {
            getNaiveTour(naive_tour, num_cities, distMatrix, i);
            twoOpt(naive_tour, num_cities, distMatrix);
            twoPointFiveOpt(naive_tour, num_cities, distMatrix);
            int tour_dist = total_distance(naive_tour, distMatrix);

            if (tour_dist < best_dist) {
                //cout << i << ": " << tour_dist << " < " << best_dist << endl;
                best_tour = naive_tour;
                best_dist = tour_dist;
            }
        }



        // Get initial tour

        //vector<int> naive_tour;
        //naive_tour.resize(num_cities);
        //getNaiveTour(naive_tour, num_cities, distMatrix);
        //int naive_dist = total_distance(naive_tour, distMatrix);

        //vector<int> chris_tour = christofides(cities, num_cities, edges);
        //int chris_dist = total_distance(chris_tour, distMatrix);
        //save_tour(chris_tour, cities, "./graphs/christofides_tour.dot");

        // Improve tour with heuristic
        //twoOpt(chris_tour, num_cities, distMatrix);
        //int two_opt_chris_dist = total_distance(chris_tour, distMatrix);
        //save_tour(chris_tour, cities, "./graphs/2opt_chris_tour.dot");

        // 2.5-Opt heuristic
        //twoPointFiveOpt(chris_tour, num_cities, distMatrix);
        //int twopointfive_opt_chris_dist = total_distance(chris_tour, distMatrix);
        //save_tour(chris_tour, cities, "./graphs/2-5opt_chris_tour.dot");

        //cout << "---------------------------------ANSWER---------------------------------" << endl;
        // Print answer
        for (int i = 0; i < num_cities; i++) {
            cout << best_tour[i] << endl;
        }

        

    }
}