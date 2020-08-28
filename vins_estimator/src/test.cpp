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
#include "vins/NonlinearFactor.h"


int main(int argc, char** argv)
{
    ros::init(argc, argv, "vins_estimator");
    /* ros::NodeHandle nh_("~");
    ros::Publisher nf_pub =  nh_.advertise<vins::NonlinearFactor>("NonlinearFactor", 1000);
    ros::Rate loop_rate(1);

    while (ros::ok()) {
        
        vins::NonlinearFactor nf;

        nf_pub.publish(nf);
        ROS_INFO("Publishing NonlinearFactor");

        loop_rate.sleep();

    } */
    
    double** para_Pose = new double* [2];
    para_Pose[0] = new double[7];
    para_Pose[1] = new double[7];
    Eigen::Map<Eigen::Matrix<double, 7, 1>> p1(para_Pose[0]);
    Eigen::Map<Eigen::Matrix<double, 7, 1>> p2(para_Pose[1]);
    //Eigen::Matrix<double, 6, 6> sqrt_info = Eigen::MatrixXd::Identity(6, 6);


    //std::shared_ptr<RelativePoseFactor> rp = std::make_shared<RelativePoseFactor>(z_rel_P, z_rel_Q, sqrt_info);
    //std::shared_ptr<AbsPositionFactor> ap = std::make_shared<AbsPositionFactor>(z_rel_P, Eigen::MatrixXd::Identity(3, 3));
    std::cout << "starting test\n" << std::endl; 

    const double delta_theta = 0.5;
    for (int i = 0; i < 10; ++i) 
    {
        //generate different para_Pose[0] && para_Pose[1]
        Eigen::Vector4d v0 = Eigen::MatrixXd::Random(4, 1);
        Eigen::Quaterniond q0(v0);
        q0.normalize();
        p1.block<3, 1>(0, 0) = Eigen::MatrixXd::Zero(3, 1);
        p1.block<4, 1>(3, 0) = q0.coeffs();
        p2.block<7, 1>(0, 0) = Eigen::MatrixXd::Zero(7, 1);

        // generate z_rel_P and z_rel_Q for test
        Eigen::Vector3d z_rel_P = q0.inverse() * Eigen::Vector3d::UnitX();
        std::shared_ptr<YawFactor> yf = std::make_shared<YawFactor>(z_rel_P, Eigen::MatrixXd::Identity(1, 1));

        std::cout << "Checking yaw position\n";
        for (int k = 0; k < 4; k++)
        {
            if (k == 0)
            {
                puts("For yawfactor, original residual is: ");
                yf->check_consistency(para_Pose);
                continue;
            }

            std::cout << "changing ypr in order, the residuals are as follow: " << std::endl;
            Eigen::Vector3d ypr = Utility::R2ypr(q0.toRotationMatrix());
            Eigen::Vector3d delta = Eigen::Vector3d(k == 1, k == 2, k == 3) * delta_theta;
            ypr += delta;
            Eigen::Quaterniond q0_new(Utility::ypr2R(ypr));
            p1.block<4, 1>(3, 0) = q0_new.coeffs();
            p1.normalize();
            yf->check_consistency(para_Pose);
        }
        
        std::cout << std::endl;
        Eigen::Vector4d v1 = Eigen::MatrixXd::Random(4, 1);
        Eigen::Quaterniond q1(v1);
        q1.normalize();
        p1.block<3, 1>(0, 0) = Eigen::MatrixXd::Zero(3, 1);
        p1.block<4, 1>(3, 0) = q1.coeffs();

        Eigen::Quaterniond z_rel_Q = q1;
        std::shared_ptr<RollPitchFactor> rpf = std::make_shared<RollPitchFactor>(z_rel_Q, Eigen::MatrixXd::Identity(2, 2));

        std::cout << "Checking rollpitch position\n";
        for (int k = 0; k < 4; k++)
        {
            if (k == 0)
            {
                puts("For rollpitchfactor, original residual is: ");
                rpf->check_consistency(para_Pose);
                continue;
            }

            std::cout << "changing ypr in order, the residuals are as follow: " << std::endl;
            Eigen::Vector3d ypr = Utility::R2ypr(q1.toRotationMatrix());
            Eigen::Vector3d delta = Eigen::Vector3d(k == 1, k == 2, k == 3) * delta_theta;
            ypr += delta;
            Eigen::Quaterniond q1_new(Utility::ypr2R(ypr));
            p1.block<4, 1>(3, 0) = q1_new.coeffs();
            p1.normalize();
            rpf->check_consistency(para_Pose);
        }
        std::cout << std::endl;
        
    }


    delete[] para_Pose[0];
    delete[] para_Pose[1];
    delete[] para_Pose;
    return 0;
}