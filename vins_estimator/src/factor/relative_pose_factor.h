#pragma once

#include <ros/assert.h>
#include <ceres/ceres.h>
#include <Eigen/Dense>
#include "../utility/utility.h"
#include "../utility/tic_toc.h"
#include "../estimator/parameters.h"

struct RelPoseFactor {
  double Header_i, Header_j; // or long

  Eigen::Vector3d rel_P;
  Eigen::Quaterniond rel_Q;
  Eigen::Matrix<double, 6, 6> cov_inv;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class RelativePoseFactor : public ceres::SizedCostFunction<6, 7, 7>
{
  public:
    RelativePoseFactor(const Eigen::Vector3d &_rel_P, const Eigen::Quaterniond &_rel_Q, const Eigen::Matrix<double, 6, 6> _sqrt_info) : rel_P(_rel_P), rel_Q(_rel_Q), sqrt_info(_sqrt_info) {}

    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
    {
    	Eigen::Vector3d P_i(parameters[0][0], parameters[0][1], parameters[0][2]);
    	Eigen::Quaterniond Q_i(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

		Eigen::Vector3d P_j(parameters[1][0], parameters[1][1], parameters[1][2]);
    	Eigen::Quaterniond Q_j(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

		Eigen::Vector3d P_w_ij = P_j - P_i;

		Eigen::Quaterniond Q_i_inverse = Q_i.inverse();

		Eigen::Vector3d P_i_ij = Q_i_inverse * P_w_ij;

		Eigen::Quaterniond Q_ij = Q_i_inverse * Q_j;

    	Eigen::Map<Eigen::Matrix<double, 6, 1>> residual(residuals);
    	residual.block<3, 1>(0, 0) = P_i_ij - rel_P;
    	residual.block<3, 1>(3, 0) = 2 * (rel_Q.inverse() * Q_ij).vec();
    	residual = sqrt_info * residual;

    	if (jacobians)
    	{
    		if (jacobians[0])
    		{
    		    Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);
    		    jacobian_pose_i.setZero();
    		    jacobian_pose_i.block<3, 3>(0, 0) = -Q_i_inverse.toRotationMatrix();
				jacobian_pose_i.block<3, 3>(0, 3) = Utility::skewSymmetric(P_i_ij);
    		    jacobian_pose_i.block<3, 3>(3, 3) = -(Utility::Qright(Q_ij) * Utility::Qleft(rel_Q.inverse())).bottomRightCorner<3, 3>();
    		    jacobian_pose_i = sqrt_info * jacobian_pose_i;
    		}
			if (jacobians[1]) {
				Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::RowMajor>> jacobian_pose_j(jacobians[1]);
    		    jacobian_pose_j.setZero();
    		    jacobian_pose_j.block<3, 3>(0, 0) = Q_i_inverse.toRotationMatrix();
    		    jacobian_pose_j.block<3, 3>(3, 3) = Utility::Qleft(rel_Q.inverse() * Q_ij).bottomRightCorner<3, 3>();
    		    jacobian_pose_j = sqrt_info * jacobian_pose_j;
			}
    	}
    	return true;
    }

    void check(double **parameters)
    {
	    double *res = new double[6];
	    double **jaco = new double *[2];
	    jaco[0] = new double[6 * 7];
		jaco[1] = new double[6 * 7];
	    Evaluate(parameters, res, jaco);
	    puts("check begins");

	    puts("my");

	    std::cout << Eigen::Map<Eigen::Matrix<double, 6, 1>>(res).transpose() << std::endl
	              << std::endl;
	    std::cout << Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::RowMajor>>(jaco[0]) << std::endl
	              << std::endl;
		std::cout << Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::RowMajor>>(jaco[1]) << std::endl
	              << std::endl;

		Eigen::Vector3d P_i(parameters[0][0], parameters[0][1], parameters[0][2]);
    	Eigen::Quaterniond Q_i(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

		Eigen::Vector3d P_j(parameters[1][0], parameters[1][1], parameters[1][2]);
    	Eigen::Quaterniond Q_j(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

		Eigen::Vector3d P_w_ij = P_j - P_i;

		Eigen::Quaterniond Q_i_inverse = Q_i.inverse();

		Eigen::Vector3d P_i_ij = Q_i_inverse * P_w_ij;

		Eigen::Quaterniond Q_ij = Q_i_inverse * Q_j;

    	Eigen::Matrix<double, 6, 1> residual;
    	residual.block<3, 1>(0, 0) = P_i_ij - rel_P;
    	residual.block<3, 1>(3, 0) = 2 * (rel_Q.inverse() * Q_ij).vec();
    	residual = sqrt_info * residual;

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
    		tmp_residual.block<3, 1>(0, 0) = P_i_ij - rel_P;
    		tmp_residual.block<3, 1>(3, 0) = 2 * (rel_Q.inverse() * Q_ij).vec();
    		tmp_residual = sqrt_info * tmp_residual;

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
    		tmp_residual.block<3, 1>(0, 0) = P_i_ij - rel_P;
    		tmp_residual.block<3, 1>(3, 0) = 2 * (rel_Q.inverse() * Q_ij).vec();
    		tmp_residual = sqrt_info * tmp_residual;

	        num_jacobian_j.col(k) = (tmp_residual - residual) / eps;
	    }
	    std::cout << num_jacobian_j << std::endl;
    }

    Eigen::Vector3d rel_P;
    Eigen::Quaterniond rel_Q;
    Eigen::Matrix<double, 6, 6> sqrt_info;
};
