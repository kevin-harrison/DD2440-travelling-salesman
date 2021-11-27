#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <iomanip>
#include <vector>
#include <string>
#include <tuple>
#include <algorithm>
#include <set>
#include <stack>


using std::tuple;
using namespace std;

void print_vector(vector<int> &vect)
{
    for (auto i : vect)
    {
        cout << i << ", ";
    }
    cout << endl;
}

// Rounded Euclidean distance
double dist(int city_a, int city_b, vector<pair<double, double>> &cities)
{
    return round(sqrt(pow(cities[city_a].first - cities[city_b].first, 2) +
                      pow(cities[city_a].second - cities[city_b].second, 2)));
}

// Calculate the distance of an entire tour
double total_distance(vector<int> &tour, float **&dist_matrix)
{
    double total_dist = 0;
    int num_tour = tour.size();

    for (int i = 0; i < num_tour - 1; i++)
    {
        total_dist += dist_matrix[tour[i]][tour[i + 1]];
    }
    return total_dist + dist_matrix[tour[num_tour - 1]][tour[0]]; // Add distance linking back to cycle start
}

// Save tour to file filename in a .dot file format
void save_tour(vector<int> &tour, vector<pair<double, double>> &cities, string filename)
{
    ofstream myfile;
    myfile.open(filename);
    myfile << "strict digraph {\nlabel=\"" << filename << "\";" << endl;
    myfile << "node [shape=plaintext, fontcolor=red, fixedsize=true, width=0.05];" << endl;

    // Add nodes
    for (int i = 0; i < tour.size(); i++)
    {
        myfile << i << " [pos=\"" << cities[i].first << "," << cities[i].second << "!\"];" << endl;
    }

    // Add edges
    for (int i = 0; i < tour.size(); i++)
    {
        int from = tour[i];
        int to = tour[(i + 1) % tour.size()];
        myfile << from << " -> " << to << " [label=\"" << dist(from, to, cities) << "\"];" << endl;
    }

    myfile << "}\n";
    myfile.close();
    return;
}

// Find naive tour
void getNaiveTour(int num_cities, vector<pair<double, double>> &cities, vector<int> &tour)
{
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

// Reverse elements in a tour from index "start" to index "end" (inclusive)
void swap_tour(vector<int> &tour, int start, int end)
{
    int temp;

    // Make sure bounds are correct
    if (start > end)
    {
        temp = start;
        start = end;
        end = temp;
    }

    // Swap elements pair by pair between [start,...,end]
    for (int i = 0; i < floor((end - start + 1) / 2); i++)
    {
        // swap tour[start+i] and tour[end-i]
        temp = tour[start + i];
        tour[start + i] = tour[end - i];
        tour[end - i] = temp;
    }
}

float getCostDiff(vector<int> &tour, float **& distMatrix, int i, int j, int num_cities)
{
    float newEdges = distMatrix[tour[i]][tour[j]] + distMatrix[tour[i + 1]][tour[(j + 1) % num_cities]];
    float prevEdges = distMatrix[tour[i]][tour[i + 1]] + distMatrix[tour[j]][tour[(j + 1) % num_cities]];
    return newEdges - prevEdges;
}

void twoOpt(vector<int> &tour, int num_cities, vector<pair<double, double>> &cities, float **&distMatrix)
{
    int improvements = 1;
    int counter = 0;
    int minCost = 1;
    int iMin = 0;
    int jMin = 0;
    int cost = 0;

label:
    if (counter > 600 || minCost == 0)
        return;
    minCost = 0;
    // Iterate over all cities to find swapping improvements
    for (int i = 0; i < num_cities - 1; i++)
    {
        for (int j = i + 2; j < num_cities; j++)
        {
            if (i == ((j + 1) % num_cities))
                continue; // Nodes are neighbors so swap is meaningless

            //double edge1 = dist(tour[i],   tour[i+1], cities);                  // dist(c1,c2)
            //double edge2 = dist(tour[j],   tour[(j+1) % num_cities], cities);   // dist(c3,c4)
            //double edge3 = dist(tour[i],   tour[j], cities);                    // dist(c1,c3)
            //double edge4 = dist(tour[i+1], tour[(j+1) % num_cities], cities);   // dist(c2,c4)
            cost = getCostDiff(tour, distMatrix, i, j, num_cities);

            if (minCost > cost)
            { // Only swap if it will result in a shorter tour
                minCost = cost;
                iMin = i + 1;
                jMin = j;
                //save_tour(tour, cities, "./graphs/tour" + to_string(counter) + ".dot");
            }
        }
    }
    if (minCost < 0)
    {
        swap_tour(tour, iMin, jMin);
    }
    counter += 1;
    goto label;
}

bool sortByDist(const tuple<int, int, float> &a,
                const tuple<int, int, float> &b)
{
    return (get<2>(a) < get<2>(b));
}

int find_cluster(int i, int *vertexClusters)
{
    if (i == vertexClusters[i])
    {
        return i;
    }
    return find_cluster(vertexClusters[i], vertexClusters);
}

bool contains(vector<int> v, int x)
{
    if (std::find(v.begin(), v.end(), x) != v.end())
    {
        return true;
    }
    return false;
}

void printGraph(vector<vector<int>> v)
{
    int i = 0;
    vector<vector<int>>::iterator row;
    vector<int>::iterator col;
    for (row = v.begin(); row != v.end(); row++)
    {
        cout << "VERTEX " << i << ": ";
        for (col = row->begin(); col != row->end(); col++)
        {
            cout << (*col) << " ";
        }
        cout << endl;
        i += 1;
    }
}

void printGraph(vector<int> *&v)
{
    for (size_t i = 0; i < v->size(); i++)
    {
        cout << "VERTEX " << i << ": ";
        for (int j = 0; j < v[i].size(); j++)
        {
            cout << v[i][j] << " ";
        }
        cout << endl;
    }
}

/* returns a vector of vectors of ints. This represents a list of vertices, where each of them has a list of edges that attach them to the MST */
vector<int> christofides(vector<pair<double, double>> &cities, int num_cities, vector<tuple<int, int, float>> edges)
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
    for (tuple<int, int, float> edge : edges)
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

void prepareDistMarixAndEdges(float **&distMatrix, vector<tuple<int, int, float> > &edges, int num_cities, vector<pair<double, double> > &cities)
{
    // Reserve edges
    int to_reserve = 0;
    // Calculate distances, populate distance matrix and edges
    for (int i = 0; i < num_cities; i++)
    {
        distMatrix[i] = new float[num_cities];
        to_reserve += i;
    }
    edges.reserve(to_reserve);


    tuple<int, int, float> edge;
    for (int i = 0; i < num_cities; i++)
    {
        for (int j = i+1; j < num_cities; j++)
        {
            int distance = dist(i, j, cities);
            distMatrix[i][j] = distance;
            distMatrix[j][i] = distance;
            edge = make_tuple(i, j, distance);
            edges.push_back(edge);
        }
    }
    sort(edges.begin(), edges.end(), sortByDist);
}

int main()
{

    std::cout << std::fixed;
    std::cout << std::setprecision(10);

    int num_cities;
    if (cin >> num_cities)
    {
        // Get input
        vector<pair<double, double>> cities;
        cities.reserve(num_cities);
        for (int i = 0; i < num_cities; i++)
        {
            double xcord, ycord;
            cin >> xcord >> ycord;
            cities[i] = make_pair(xcord, ycord);
        }

        float **distMatrix;
        distMatrix = new float *[num_cities];
        vector<tuple<int, int, float>> edges;
        prepareDistMarixAndEdges(distMatrix, edges, num_cities, cities);

        // Print edges
        /* cout << "EDGES: " << endl;
        cout << to_string(edges.size()) << endl;
        for (int i = 0; i < edges.size(); i++) {
            cout << "(" << to_string(get<0>(edges[i])) << ", " << to_string(get<1>(edges[i])) << ", " << to_string(get<2>(edges[i])) << ")" << endl;
        }
        cout << endl; */

        // Get initial tour

        vector<int> initial_tour;
        initial_tour.resize(num_cities);
        //getNaiveTour(num_cities, cities, initial_tour);
        //double naive_dist = total_distance(initial_tour, distMatrix);
        //save_tour(initial_tour, cities, "./graphs/naive_tour.dot");
        //twoOpt(initial_tour, num_cities, cities, distMatrix);
        //double two_opt_naive_dist = total_distance(initial_tour, distMatrix);
        //save_tour(initial_tour, cities, "./graphs/2opt_naive_tour.dot");

        vector<int> chris_tour = christofides(cities, num_cities, edges);
        //double chris_dist = total_distance(chris_tour, distMatrix);
        //save_tour(chris_tour, cities, "./graphs/christofides_tour.dot");

        // Improve tour with heuristic
        twoOpt(chris_tour, num_cities, cities, distMatrix);
        //double two_opt_chris_dist = total_distance(chris_tour, distMatrix);
        //save_tour(chris_tour, cities, "./graphs/2opt_chris_tour.dot");

        //cout << "---------------------------------ANSWER---------------------------------" << endl;
        // Print answer
        for (int i = 0; i < num_cities; i++)
        {
            cout << chris_tour[i] << endl;
        }

        /* cout << "Naive dist: " << naive_dist << endl;
        cout << "Chris dist: " << chris_dist << endl;
        cout << "2-opt naive dist: " << two_opt_naive_dist << endl;
        cout << "2-opt chris dist: " << two_opt_chris_dist << endl; */
    }
}