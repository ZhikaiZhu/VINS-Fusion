/*******************************************************
 * Copyright (C) 2019, Aerial Robotics Group, Hong Kong University of Science and Technology
 * 
 * This file is part of VINS.
 * 
 * Licensed under the GNU General Public License v3.0;
 * you may not use this file except in compliance with the License.
 *
 * Author: Qin Tong (qintonguav@gmail.com)
 *******************************************************/

#include <stdio.h>
#include <math.h>
#include "utility/visualization.h"
#include <iostream>
#include <memory>
#include <queue>
#include <map>
#include <thread>
#include <mutex>
#include <random>
#include <ros/ros.h>
#include <cv_bridge/cv_bridge.h>
#include <opencv2/opencv.hpp>
#include "estimator/estimator.h"
#include "estimator/parameters.h"

std::random_device rd{};
std::mt19937 gen{rd()};

std::normal_distribution<> position_dist{0, 10.0};
std::normal_distribution<> orientation_dist{0, 10.0};

int main(int argc, char **argv)
{
    double** para_Pose = new double *[2];
    para_Pose[0] = new double[7];
    para_Pose[1] = new double[7];
    Eigen::Map<Eigen::Matrix<double, 7, 1>> p1(para_Pose[0]);
    Eigen::Map<Eigen::Matrix<double, 7, 1>> p2(para_Pose[1]);
    Eigen::Matrix<double, 6, 6> sqrt_info = Eigen::MatrixXd::Identity(6, 6);

    //generate different rel_P, rel_Q
    Eigen::Vector3d rel_P(4.0, 4.0, 5.0);
    Eigen::Quaterniond rel_Q(1.0, 0.0, 0.0, 0.0);
    std::shared_ptr<RelativePoseFactor> rp = std::make_shared<RelativePoseFactor>(rel_P, rel_Q, sqrt_info);
    std::shared_ptr<AbsPositionFactor> ap = std::make_shared<AbsPositionFactor>(rel_P, Eigen::MatrixXd::Identity(3, 3));
    std::shared_ptr<RollPitchFactor> rpf = std::make_shared<RollPitchFactor>(rel_Q, Eigen::MatrixXd::Identity(2, 2));
    std::shared_ptr<YawFactor> yf = std::make_shared<YawFactor>(rel_P, Eigen::MatrixXd::Identity(1, 1));
    std::cout << "starting test\n"; 

    for (int i = 0; i < 10; ++i) {
        //generate different para_Pose[0] && para_Pose[1]
        p1.block<3, 1>(0, 0) = Eigen::MatrixXd::Zero(3, 1);
        p1.block<4, 1>(3, 0) = Eigen::MatrixXd::Random(4, 1);
        p1.normalize();
        p1.block<3, 1>(0, 0) = Eigen::MatrixXd::Random(3, 1);
        p2.block<3, 1>(0, 0) = Eigen::MatrixXd::Zero(3, 1);
        p2.block<4, 1>(3, 0) = Eigen::MatrixXd::Random(4, 1);
        p2.normalize();
        p2.block<3, 1>(0, 0) = Eigen::MatrixXd::Random(3, 1);
        std::cout << "Checking relative pose\n" << std::endl;
        rp->check(para_Pose);
        std::cout << "Checking abs position\n" << std::endl;
        ap->check(para_Pose);
        std::cout << "Checking rollpitch position\n" << std::endl;
        rpf->check(para_Pose);
        std::cout << "Checking yaw position\n" << std::endl;
        yf->check(para_Pose);
    }


    delete[] para_Pose[0];
    delete[] para_Pose[1];
    delete[] para_Pose;
    return 0;
}
