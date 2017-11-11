#include <iostream>
#include "gtest/gtest.h"
#include <vector>

#include "helper_functions.h"

#define MAX_ABSOLUTE_ERROR 1e-4

using Eigen::MatrixXd;
using Eigen::VectorXd;

TEST(HomogenousTransformation, NoTransformationTest){

    // Case 1
    float theta = 0.0F;
    float x = 0.0F;
    float y = 0.0F;

    MatrixXd rotMatrix = getHomogenousTransformationMatrix(theta, x, y);

    Eigen::VectorXd someVector(3);
    someVector << 1.6F, 2.3F, 1.0F;

    Eigen::VectorXd result = rotMatrix * someVector;

    // Result and someVector should be equal
    ASSERT_TRUE((result - someVector).norm() < MAX_ABSOLUTE_ERROR);
}

TEST(HomogenousTransformation, RotationTest){

    // Case 2
    float theta = -M_PI/2.0F;
    float x = 4.0F;
    float y = 5.0F;

    MatrixXd rotMatrix = getHomogenousTransformationMatrix(theta, x, y);

    Eigen::VectorXd someVector(3);
    someVector << 2.0F, 2.0F, 1.0F;

    Eigen::VectorXd result = rotMatrix * someVector;

    Eigen::VectorXd expectedOutput(3);
    expectedOutput << 6.0F, 3.0F, 1.0F;

    ASSERT_TRUE((result - expectedOutput).norm() < MAX_ABSOLUTE_ERROR);
}
