#include <iostream>
#include <vector>
#include "remus.h"
#include "controller.h"

using namespace std;

int main() {
    size_t step_number;
    double heading_ref;
    Remus a({1.5, 0.0, 0, 0 ,0 ,0});

    cout << "Please input the heading reference: ";
    cin >> heading_ref;
    cout << "Please input the step number: ";
    cin >> step_number;
    RunRemus(a, heading_ref, step_number);
    std::cout << "The velocity of the vehilce is: ";
    std::cout << a.velocity_.transpose() << std::endl;

    return 0;
}
