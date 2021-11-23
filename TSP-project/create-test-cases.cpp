#include <iostream>
#include <sstream>
#include <stdlib.h>     /* srand, rand */
#include <string.h>
#include <sys/time.h>

using namespace std;
using std::to_string;

string makeTest() {
    string test = "";
    float HI = 180;
    float LO = -HI;
    int numCities = 10;//rand() % 150 + 50;
    test = test + to_string(numCities) + "\n";
    for (size_t i = 0; i < numCities; i++)
    {
        float x = LO + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(HI-LO)));
        float y = LO + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(HI-LO)));
        test = test + to_string(x) + " " + to_string(y) + "\n";
    }
    
    return test;
}

int main(int argc, char * argv[])
{
    // Get number of tests in program argument
    istringstream ss(argv[1]);
    int numOfTests;
    if (!(ss >> numOfTests)) {
        cerr << "Invalid number: " << argv[1] << '\n';
    } else if (!ss.eof()) {
        cerr << "Trailing characters after number: " << argv[1] << '\n';
    }
    
    // Seed random number generator with microseconds
    timeval t1;
    gettimeofday(&t1, NULL);
    srand(t1.tv_usec * t1.tv_sec);
 
    // Send test cases to stdout
    for(int i=0; i<numOfTests; i++) {
        string test = makeTest();
        cout << test;
    }
    return 0;
}