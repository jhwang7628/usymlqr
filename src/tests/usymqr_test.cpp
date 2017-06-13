#include <iostream>
#include "Widget.h"
#include "Eigen/Dense"
//##############################################################################
int main() {
    Widget w; 
    w.print(); 

    Eigen::Matrix3d m;
    m << 1, 2, 3,
         4, 5, 6,
         7, 8, 9;
    std::cout << m << std::endl;
}
//##############################################################################
