#include <iostream>
#include <vector>
#include "remus.h"

using namespace std;

int main() {
    size_t step_number;
    Remus a({1.5, 0.1, 0, 0 ,0 ,0});

//    std::cout << "The velocity of the vehilce is: ";
//    std::cout << a.velocity_.transpose() << std::endl;
//    std::cout << "The transformation matrix is given by " << std::endl;
//    std::cout << Remus::TransformationMatrix({pi/6, 0, 0}) << std::endl;
//    std::cout << "The inertia matrix of the vehicle is: " << std::endl;
//    std::cout << a.RigidBodyInertiaMatrix() << std::endl;
//    std::cout << "The centripetal-coriolis matrix of the rigid body is" << std::endl;
//    std::cout << a.RigidBodyCCMatrix(a.velocity_) << std::endl;
//    std::cout << "The centripetal-coriolis matrix of the added mass is" << std::endl;
//    std::cout << a.AddedMassCCMatrix(a.velocity_) << std::endl;
//    std::cout << "The centripetal-coriolis matrix of the displaced mass is" << std::endl;
//    std::cout << a.DisplacedMassCCMatrix(a.velocity_) << std::endl;
//    std::cout << "The viscous damping vector is: ";
//    std::cout << a.ViscousDampingVector(a.velocity_).transpose() << std::endl;
//    std::cout << "The restoring force vector is: ";
//    std::cout << a.RestoringForceVector(a.position_).transpose() << std::endl;
//    std::cout << "The state derivative vector is: " << std::endl;
//    std::cout << a.StateDerivative(a.velocity_, a.position_, 0).transpose() << std::endl;

    cout << "Please input the step number: ";
    cin >> step_number;
    RunRemus(a, step_number);
    std::cout << "The velocity of the vehilce is: ";
    std::cout << a.velocity_.transpose() << std::endl;

    return 0;
}
