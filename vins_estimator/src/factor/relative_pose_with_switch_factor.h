#pragma once
#include <ros/assert.h>
#include <iostream>
#include <eigen3/Eigen/Dense>

#include "../utility/utility.h"
#include "../estimator/parameters.h"
#include "../utility/tic_toc.h"

#include <ceres/ceres.h>


class RelPoseWithSwitchFactor: public ceres::SizedCostFunction<6, 7, 7, 1>
{
public:
    RelPoseWithSwitchFactor(const Vector3d& _z_rel_P, const Quaterniond& _z_rel_Q, const Eigen::Matrix<double, 6, 6>& _sqrt_info): z_rel_P(_z_rel_P), z_rel_Q(_z_rel_Q), sqrt_info(_sqrt_info) {}
    
    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const
    {
        Eigen::Vector3d P_i(parameters[0][0], parameters[0][1], parameters[0][2]);
    	Eigen::Quaterniond Q_i(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

		Eigen::Vector3d P_j(parameters[1][0], parameters[1][1], parameters[1][2]);
    	Eigen::Quaterniond Q_j(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

        double S_ij{parameters[2][0]};

        if (S_ij > 1.0)
        {
            S_ij = 1.0;
        }
        else if (S_ij < 0.0)
        {
            S_ij = 0.0;
        }

		Eigen::Vector3d P_w_ij = P_j - P_i;

		Eigen::Quaterniond Q_i_inverse = Q_i.inverse();

		Eigen::Vector3d P_i_ij = Q_i_inverse * P_w_ij;

		Eigen::Quaterniond Q_ij = Q_i_inverse * Q_j;

    	Eigen::Map<Eigen::Matrix<double, 6, 1>> residual(residuals);
    	residual.block<3, 1>(0, 0) = P_i_ij - z_rel_P;
    	residual.block<3, 1>(3, 0) = 2 * (z_rel_Q.inverse() * Q_ij).vec();
    	residual = S_ij * sqrt_info * residual;

    	if (jacobians)
    	{
    		if (jacobians[0])
    		{
    		    Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);
    		    jacobian_pose_i.setZero();
    		    jacobian_pose_i.block<3, 3>(0, 0) = -Q_i_inverse.toRotationMatrix();
				jacobian_pose_i.block<3, 3>(0, 3) = Utility::skewSymmetric(P_i_ij);
    		    jacobian_pose_i.block<3, 3>(3, 3) = -(Utility::Qright(Q_ij) * Utility::Qleft(z_rel_Q.inverse())).bottomRightCorner<3, 3>();
    		    jacobian_pose_i = S_ij * sqrt_info * jacobian_pose_i;
    		}
			if (jacobians[1]) 
			{
				Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::RowMajor>> jacobian_pose_j(jacobians[1]);
    		    jacobian_pose_j.setZero();
    		    jacobian_pose_j.block<3, 3>(0, 0) = Q_i_inverse.toRotationMatrix();
    		    jacobian_pose_j.block<3, 3>(3, 3) = Utility::Qleft(z_rel_Q.inverse() * Q_ij).bottomRightCorner<3, 3>();
    		    jacobian_pose_j = S_ij * sqrt_info * jacobian_pose_j;
			}
            if (jacobians[2])
            {
                Eigen::Map<Eigen::Matrix<double, 6, 1>> jacobian_switch(jacobians[2]);
                jacobian_switch.setZero();
                jacobian_switch.block<3, 1>(0, 0) = P_i_ij - z_rel_P;
                jacobian_switch.block<3, 1>(3, 0) = 2 * (z_rel_Q.inverse() * Q_ij).vec();
                jacobian_switch = sqrt_info * jacobian_switch;
            }
    	}
    	return true;
    }

    void check(double** parameters)
    {
	    double* res = new double[6];
	    double** jaco = new double* [3];
	    jaco[0] = new double[6 * 7];
		jaco[1] = new double[6 * 7];
        jaco[2] = new double[6 * 1];
	    Evaluate(parameters, res, jaco);
	    puts("check begins");

	    puts("my");

	    std::cout << Eigen::Map<Eigen::Matrix<double, 6, 1>>(res).transpose() << std::endl
	              << std::endl;
	    std::cout << Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::RowMajor>>(jaco[0]) << std::endl
	              << std::endl;
		std::cout << Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::RowMajor>>(jaco[1]) << std::endl
	              << std::endl;
        std::cout << Eigen::Map<Eigen::Matrix<double, 6, 1>>(jaco[2]) << std::endl
	              << std::endl;

		Eigen::Vector3d P_i(parameters[0][0], parameters[0][1], parameters[0][2]);
    	Eigen::Quaterniond Q_i(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

		Eigen::Vector3d P_j(parameters[1][0], parameters[1][1], parameters[1][2]);
    	Eigen::Quaterniond Q_j(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

        double S_ij{parameters[2][0]};

		Eigen::Vector3d P_w_ij = P_j - P_i;

		Eigen::Quaterniond Q_i_inverse = Q_i.inverse();

		Eigen::Vector3d P_i_ij = Q_i_inverse * P_w_ij;

		Eigen::Quaterniond Q_ij = Q_i_inverse * Q_j;

    	Eigen::Matrix<double, 6, 1> residual;
    	residual.block<3, 1>(0, 0) = P_i_ij - z_rel_P;
    	residual.block<3, 1>(3, 0) = 2 * (z_rel_Q.inverse() * Q_ij).vec();
    	residual = S_ij * sqrt_info * residual;

	    puts("num");
	    std::cout << residual.transpose() << std::endl << std::endl;

	    const double eps = 1e-6;
	    Eigen::Matrix<double, 6, 6> num_jacobian_i;
	    for (int k = 0; k < 6; k++)
	    {
	    	Eigen::Vector3d P_i_new(parameters[0][0], parameters[0][1], parameters[0][2]);
	    	Eigen::Quaterniond Q_i_new(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

	        int a = k / 3, b = k % 3;
	        Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;

	        if (a == 0)
	            P_i_new += delta;
	        else if (a == 1)
	            Q_i_new = Q_i_new * Utility::deltaQ(delta);

			P_w_ij = P_j - P_i_new;

			Q_i_inverse = Q_i_new.inverse();

			P_i_ij = Q_i_inverse * P_w_ij;

			Q_ij = Q_i_inverse * Q_j;

    		Eigen::Matrix<double, 6, 1> tmp_residual;
    		tmp_residual.block<3, 1>(0, 0) = P_i_ij - z_rel_P;
    		tmp_residual.block<3, 1>(3, 0) = 2 * (z_rel_Q.inverse() * Q_ij).vec();
    		tmp_residual = S_ij * sqrt_info * tmp_residual;

	        num_jacobian_i.col(k) = (tmp_residual - residual) / eps;
	    }
	    std::cout << num_jacobian_i << std::endl << std::endl;
		Eigen::Matrix<double, 6, 6> num_jacobian_j;
	    for (int k = 0; k < 6; k++)
	    {
	    	Eigen::Vector3d P_j_new(parameters[1][0], parameters[1][1], parameters[1][2]);
	    	Eigen::Quaterniond Q_j_new(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

	        int a = k / 3, b = k % 3;
	        Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;

	        if (a == 0)
	            P_j_new += delta;
	        else if (a == 1)
	            Q_j_new = Q_j_new * Utility::deltaQ(delta);

			P_w_ij = P_j_new - P_i;

			Q_i_inverse = Q_i.inverse();

			P_i_ij = Q_i_inverse * P_w_ij;

			Q_ij = Q_i_inverse * Q_j_new;

    		Eigen::Matrix<double, 6, 1> tmp_residual;
    		tmp_residual.block<3, 1>(0, 0) = P_i_ij - z_rel_P;
    		tmp_residual.block<3, 1>(3, 0) = 2 * (z_rel_Q.inverse() * Q_ij).vec();
    		tmp_residual = S_ij * sqrt_info * tmp_residual;

	        num_jacobian_j.col(k) = (tmp_residual - residual) / eps;
	    }
	    std::cout << num_jacobian_j << std::endl;
        Eigen::Matrix<double, 6, 1> num_jacobian_switch;
	    for (int k = 0; k < 1; k++)
	    {
	    	Eigen::Vector3d P_i(parameters[0][0], parameters[0][1], parameters[0][2]);
    	    Eigen::Quaterniond Q_i(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

		    Eigen::Vector3d P_j(parameters[1][0], parameters[1][1], parameters[1][2]);
    	    Eigen::Quaterniond Q_j(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

            double S_ij_new{parameters[2][0]};

            Eigen::Vector3d P_w_ij = P_j - P_i;

		    Eigen::Quaterniond Q_i_inverse = Q_i.inverse();

		    Eigen::Vector3d P_i_ij = Q_i_inverse * P_w_ij;

		    Eigen::Quaterniond Q_ij = Q_i_inverse * Q_j;

            S_ij_new += eps;

    		Eigen::Matrix<double, 6, 1> tmp_residual;
    		tmp_residual.block<3, 1>(0, 0) = P_i_ij - z_rel_P;
    		tmp_residual.block<3, 1>(3, 0) = 2 * (z_rel_Q.inverse() * Q_ij).vec();
    		tmp_residual = S_ij_new * sqrt_info * tmp_residual;

	        num_jacobian_switch.col(k) = (tmp_residual - residual) / eps;
	    }
	    std::cout << num_jacobian_switch << std::endl << std::endl;
    }

    Eigen::Vector3d z_rel_P;
    Eigen::Quaterniond z_rel_Q;
    Eigen::Matrix<double, 6, 6> sqrt_info;
};
