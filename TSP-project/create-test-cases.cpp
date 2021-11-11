#include <iostream>
#include <sstream>
#include <stdlib.h>     /* srand, rand */
#include <fstream>
#include <string.h>
#include <ctime>

using namespace std;

string makeTest() {
    string test = "";
    float HI = 180;
    float LO = -HI;
    int numCities = rand() % 10 + 1; //Random integer from 1 to 100
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

    std::istringstream ss(argv[1]);
    int x;
    if (!(ss >> x)) {
        std::cerr << "Invalid number: " << argv[1] << '\n';
    } else if (!ss.eof()) {
        std::cerr << "Trailing characters after number: " << argv[1] << '\n';
    }


    srand((uint) * argv[1]);
    int numOfTests = 100;
    for(int i=0; i<numOfTests; i++) {
        string test = makeTest();
        cout << test;
        // TODO: Check if output is malformed somehow
    }
    return 0;
}