#include <iostream>
#include <vector>
#include "remus.h"
#include "controller.h"

using namespace std;

int main() {
    size_t step_number;
    double heading_ref;
    Remus a({1.0, 0, 0, 0 ,0 ,0}, {0, 0, 0, 0, 0, 0}, {0, 0.2, 0}, 4);
    MovingMassController controller(a, 0.3, 0.8, 1, 1, 1);

    cout << "Please input the heading reference: ";
    cin >> heading_ref;
    cout << "Please input the step number: ";
    cin >> step_number;

    controller.ControlRemus(a, heading_ref, step_number);

    return 0;
}
