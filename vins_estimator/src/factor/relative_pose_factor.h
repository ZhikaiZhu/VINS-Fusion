#pragma once
#include <ros/assert.h>
#include <iostream>
#include <eigen3/Eigen/Dense>

#include "../utility/utility.h"
#include "../estimator/parameters.h"
#include "../utility/tic_toc.h"

#include <ceres/ceres.h>

struct RelPoseFactor
{
    double Header_i, Header_j;
    Eigen::Vector3d z_rel_P;
    Eigen::Quaterniond z_rel_Q;
	//Eigen::Vector3d z_rel_Yaw;
    Eigen::Matrix<double, 6, 6> cov_inv;

	/*Eigen::Matrix<double, 3, 3> relP_cov_inv;
	Eigen::Matrix<double, 2, 2> relRP_cov_inv;
	Eigen::Matrix<double, 1, 1> relYaw_cov_inv; */

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class RelativePoseFactor: public ceres::SizedCostFunction<6, 7, 7>
{
public:
    RelativePoseFactor(const Vector3d& _z_rel_P, const Quaterniond& _z_rel_Q, const Eigen::Matrix<double, 6, 6>& _sqrt_info): z_rel_P(_z_rel_P), z_rel_Q(_z_rel_Q), sqrt_info(_sqrt_info) {}
    
    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const
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
    	residual.block<3, 1>(0, 0) = P_i_ij - z_rel_P;
    	residual.block<3, 1>(3, 0) = 2 * (z_rel_Q.inverse() * Q_ij).vec();
    	residual = sqrt_info * residual;

    	if (jacobians)
    	{
    		if (jacobians[0])
    		{
    		    Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);
    		    jacobian_pose_i.setZero();
    		    jacobian_pose_i.block<3, 3>(0, 0) = -Q_i_inverse.toRotationMatrix();
				jacobian_pose_i.block<3, 3>(0, 3) = Utility::skewSymmetric(P_i_ij);
    		    jacobian_pose_i.block<3, 3>(3, 3) = -(Utility::Qright(Q_ij) * Utility::Qleft(z_rel_Q.inverse())).bottomRightCorner<3, 3>();
    		    jacobian_pose_i = sqrt_info * jacobian_pose_i;
    		}
			if (jacobians[1]) 
			{
				Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::RowMajor>> jacobian_pose_j(jacobians[1]);
    		    jacobian_pose_j.setZero();
    		    jacobian_pose_j.block<3, 3>(0, 0) = Q_i_inverse.toRotationMatrix();
    		    jacobian_pose_j.block<3, 3>(3, 3) = Utility::Qleft(z_rel_Q.inverse() * Q_ij).bottomRightCorner<3, 3>();
    		    jacobian_pose_j = sqrt_info * jacobian_pose_j;
			}
    	}
    	return true;
    }

    void check(double** parameters)
    {
	    double* res = new double[6];
	    double** jaco = new double* [2];
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
    	residual.block<3, 1>(0, 0) = P_i_ij - z_rel_P;
    	residual.block<3, 1>(3, 0) = 2 * (z_rel_Q.inverse() * Q_ij).vec();
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
    		tmp_residual.block<3, 1>(0, 0) = P_i_ij - z_rel_P;
    		tmp_residual.block<3, 1>(3, 0) = 2 * (z_rel_Q.inverse() * Q_ij).vec();
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
    		tmp_residual.block<3, 1>(0, 0) = P_i_ij - z_rel_P;
    		tmp_residual.block<3, 1>(3, 0) = 2 * (z_rel_Q.inverse() * Q_ij).vec();
    		tmp_residual = sqrt_info * tmp_residual;

	        num_jacobian_j.col(k) = (tmp_residual - residual) / eps;
	    }
	    std::cout << num_jacobian_j << std::endl;
    }

    Eigen::Vector3d z_rel_P;
    Eigen::Quaterniond z_rel_Q;
    Eigen::Matrix<double, 6, 6> sqrt_info;
};

class RelPositionFactor: public ceres::SizedCostFunction<3, 7, 7>
{
	public:
    RelPositionFactor(const Eigen::Vector3d &_rel_P, const Eigen::Matrix<double, 3, 3> _sqrt_info) : rel_P(_rel_P), sqrt_info(_sqrt_info) {}

    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
    {
    	Eigen::Vector3d P_i(parameters[0][0], parameters[0][1], parameters[0][2]);
		Eigen::Vector3d P_j(parameters[1][0], parameters[1][1], parameters[1][2]);

    	Eigen::Map<Eigen::Matrix<double, 3, 1>> residual(residuals);
    	residual.block<3, 1>(0, 0) = P_j - P_i - rel_P;
    	residual = sqrt_info * residual;

    	if (jacobians)
    	{
    		if (jacobians[0])
    		{
    		    Eigen::Map<Eigen::Matrix<double, 3, 7, Eigen::RowMajor>> jacobian_position_i(jacobians[0]);
    		    jacobian_position_i.setZero();
    		    jacobian_position_i.block<3, 3>(0, 0) = -Eigen::Matrix<double, 3, 3>::Identity();
    		    jacobian_position_i = sqrt_info * jacobian_position_i;
    		}
			if (jacobians[1])
			{
				Eigen::Map<Eigen::Matrix<double, 3, 7, Eigen::RowMajor>> jacobian_position_j(jacobians[1]);
    		    jacobian_position_j.setZero();
    		    jacobian_position_j.block<3, 3>(0, 0) = Eigen::Matrix<double, 3, 3>::Identity();
    		    jacobian_position_j = sqrt_info * jacobian_position_j;
			}
    	}
    	return true;
    }

    void check(double **parameters)
    {
	    double *res = new double[3];
	    double **jaco = new double *[2];
	    jaco[0] = new double[3 * 7];
		jaco[1] = new double[3 * 7];
	    Evaluate(parameters, res, jaco);
	    puts("check begins");

	    puts("my");

	    std::cout << Eigen::Map<Eigen::Matrix<double, 3, 1>>(res).transpose() << std::endl
	              << std::endl;
	    std::cout << Eigen::Map<Eigen::Matrix<double, 3, 7, Eigen::RowMajor>>(jaco[0]) << std::endl
	              << std::endl;
		std::cout << Eigen::Map<Eigen::Matrix<double, 3, 7, Eigen::RowMajor>>(jaco[1]) << std::endl
	              << std::endl;

		Eigen::Vector3d P_i(parameters[0][0], parameters[0][1], parameters[0][2]);
		Eigen::Vector3d P_j(parameters[1][0], parameters[1][1], parameters[1][2]);

    	Eigen::Matrix<double, 3, 1> residual;
    	residual.block<3, 1>(0, 0) = P_j - P_i - rel_P;
    	residual = sqrt_info * residual;

	    puts("num");
	    std::cout << residual.transpose() << std::endl;

	    const double eps = 1e-6;
	    Eigen::Matrix<double, 3, 3> num_jacobian_i;
	    for (int k = 0; k < 3; k++)
	    {
	    	Eigen::Vector3d P_i_new(parameters[0][0], parameters[0][1], parameters[0][2]);

	        int b = k % 3;
	        Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;

	        P_i_new += delta;

    		Eigen::Matrix<double, 3, 1> tmp_residual;
    		tmp_residual.block<3, 1>(0, 0) = P_j - P_i_new - rel_P;
    		tmp_residual = sqrt_info * tmp_residual;

	        num_jacobian_i.col(k) = (tmp_residual - residual) / eps;
	    }
	    std::cout << num_jacobian_i << std::endl;
		Eigen::Matrix<double, 3, 3> num_jacobian_j;
	    for (int k = 0; k < 3; k++)
	    {
	    	Eigen::Vector3d P_j_new(parameters[1][0], parameters[1][1], parameters[1][2]);

	        int b = k % 3;
	        Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;

	        P_j_new += delta;

    		Eigen::Matrix<double, 3, 1> tmp_residual;
    		tmp_residual.block<3, 1>(0, 0) = P_j_new - P_i - rel_P;
    		tmp_residual = sqrt_info * tmp_residual;

	        num_jacobian_j.col(k) = (tmp_residual - residual) / eps;
	    }
	    std::cout << num_jacobian_j << std::endl;
    }

    Eigen::Vector3d rel_P;
    Eigen::Matrix<double, 3, 3> sqrt_info;
};

class RelRollPitchFactor : public ceres::SizedCostFunction<2, 7, 7>
{
  public:
    RelRollPitchFactor(const Eigen::Quaterniond &_rel_rp, const Eigen::Matrix<double, 2, 2> _sqrt_info) : rel_rp(_rel_rp), sqrt_info(_sqrt_info) {}

    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
    {
    	Eigen::Quaterniond Q_i(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
		Eigen::Quaterniond Q_j(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
		Eigen::Quaterniond Q_ij = Q_i.inverse() * Q_j;
		Eigen::Vector3d g_ij = Q_ij.inverse() * (-Eigen::Vector3d::UnitZ());
		Eigen::Vector3d g_w = rel_rp * g_ij;
    	Eigen::Map<Eigen::Matrix<double, 2, 1>> residual(residuals);
    	residual.block<2, 1>(0, 0) = g_w.segment(0, 2);
    	residual = sqrt_info * residual;

    	if (jacobians)
    	{
    		if (jacobians[0])
    		{
    		    Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_orientation_i(jacobians[0]);
    		    jacobian_orientation_i.setZero();
				Eigen::Matrix<double, 3, 3> jacobian_orientation = -rel_rp.toRotationMatrix() * Utility::skewSymmetric(g_ij) * Q_ij.toRotationMatrix().transpose();
        		Eigen::Matrix<double, 2, 3> reduce(2, 3);
				reduce << 1.0, 0.0, 0.0,
				          0.0, 1.0, 0.0;
				jacobian_orientation_i.block<2, 3>(0, 3) = sqrt_info * reduce * jacobian_orientation;
    		}
			if (jacobians[1])
    		{
    		    Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_orientation_j(jacobians[1]);
    		    jacobian_orientation_j.setZero();
				Eigen::Matrix<double, 3, 3> jacobian_orientation = rel_rp.toRotationMatrix() * Utility::skewSymmetric(g_ij);
        		Eigen::Matrix<double, 2, 3> reduce(2, 3);
				reduce << 1.0, 0.0, 0.0,
				          0.0, 1.0, 0.0;
				jacobian_orientation_j.block<2, 3>(0, 3) = sqrt_info * reduce * jacobian_orientation;
    		}

    	}
    	return true;
    }

    void check(double **parameters)
    {
	    double *res = new double[2];
	    double **jaco = new double *[2];
	    jaco[0] = new double[2 * 7];
		jaco[1] = new double[2 * 7];
	    Evaluate(parameters, res, jaco);
	    puts("check begins");

	    puts("my");

	    std::cout << Eigen::Map<Eigen::Matrix<double, 2, 1>>(res).transpose() << std::endl
	              << std::endl;
	    std::cout << Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>>(jaco[0]) << std::endl
	              << std::endl;
		std::cout << Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>>(jaco[1]) << std::endl
	              << std::endl;

		Eigen::Quaterniond Q_i(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
		Eigen::Quaterniond Q_j(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
		Eigen::Quaterniond Q_ij = Q_i.inverse() * Q_j;		
		Eigen::Vector3d g_ij = Q_ij.inverse() * (-Eigen::Vector3d::UnitZ());
		Eigen::Vector3d g_w = rel_rp * g_ij;
    	Eigen::Matrix<double, 2, 1> residual;
    	residual.block<2, 1>(0, 0) = g_w.segment(0, 2);
    	residual = sqrt_info * residual;

	    puts("num");
	    std::cout << residual.transpose() << std::endl;

	    const double eps = 1e-6;
	    Eigen::Matrix<double, 2, 3> num_jacobian_i;
	    for (int k = 0; k < 3; k++)
	    {
	    	Eigen::Quaterniond Q_i_new(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

	        int b = k % 3;
	        Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;

			Q_i_new = Q_i_new * Utility::deltaQ(delta);
			Eigen::Quaterniond Q_ij_new = Q_i_new.inverse() * Q_j;
			g_ij = Q_ij_new.inverse() * (-Eigen::Vector3d::UnitZ());
			g_w = rel_rp * g_ij;
    		Eigen::Matrix<double, 2, 1> tmp_residual;
    		tmp_residual.block<2, 1>(0, 0) = g_w.segment(0, 2);
    		tmp_residual = sqrt_info * tmp_residual;

	        num_jacobian_i.col(k) = (tmp_residual - residual) / eps;
	    }
	    std::cout << num_jacobian_i << std::endl;
		Eigen::Matrix<double, 2, 3> num_jacobian_j;
	    for (int k = 0; k < 3; k++)
	    {
	    	Eigen::Quaterniond Q_j_new(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

	        int b = k % 3;
	        Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;

			Q_j_new = Q_j_new * Utility::deltaQ(delta);
			Eigen::Quaterniond Q_ij_new = Q_i.inverse() * Q_j_new;
			g_ij = Q_ij_new.inverse() * (-Eigen::Vector3d::UnitZ());
			g_w = rel_rp * g_ij;
    		Eigen::Matrix<double, 2, 1> tmp_residual;
    		tmp_residual.block<2, 1>(0, 0) = g_w.segment(0, 2);
    		tmp_residual = sqrt_info * tmp_residual;

	        num_jacobian_j.col(k) = (tmp_residual - residual) / eps;
	    }
	    std::cout << num_jacobian_j << std::endl;
    }

    Eigen::Quaterniond rel_rp;
    Eigen::Matrix<double, 2, 2> sqrt_info;
};

class RelYawFactor : public ceres::SizedCostFunction<1, 7, 7>
{
  public:
    RelYawFactor(const Eigen::Vector3d &_rel_yaw, const Eigen::Matrix<double, 1, 1> _sqrt_info) : rel_yaw(_rel_yaw), sqrt_info(_sqrt_info) {}

    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
    {
    	Eigen::Quaterniond Q_i(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
		Eigen::Quaterniond Q_j(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
		Eigen::Quaterniond Q_ij = Q_i.inverse() * Q_j;
		Eigen::Vector3d yaw_ij = Q_ij * rel_yaw;
    	Eigen::Map<Eigen::Matrix<double, 1, 1>> residual(residuals);
    	residual.block<1, 1>(0, 0) = yaw_ij.segment(1, 1);
    	residual = sqrt_info * residual;

    	if (jacobians)
    	{
    		if (jacobians[0])
    		{
    		    Eigen::Map<Eigen::Matrix<double, 1, 7, Eigen::RowMajor>> jacobian_orientation_i(jacobians[0]);
    		    jacobian_orientation_i.setZero();
				Eigen::Matrix<double, 3, 3> jacobian_orientation = Q_ij.toRotationMatrix() * Utility::skewSymmetric(rel_yaw) * Q_ij.toRotationMatrix().transpose();
        		Eigen::Matrix<double, 1, 3> reduce(1, 3);
				reduce << 0.0, 1.0, 0.0;
				jacobian_orientation_i.block<1, 3>(0, 3) = sqrt_info * reduce * jacobian_orientation;
    		}
			if (jacobians[1])
    		{
    		    Eigen::Map<Eigen::Matrix<double, 1, 7, Eigen::RowMajor>> jacobian_orientation_j(jacobians[1]);
    		    jacobian_orientation_j.setZero();
				Eigen::Matrix<double, 3, 3> jacobian_orientation = -Q_ij.toRotationMatrix() * Utility::skewSymmetric(rel_yaw);
        		Eigen::Matrix<double, 1, 3> reduce(1, 3);
				reduce << 0.0, 1.0, 0.0;
				jacobian_orientation_j.block<1, 3>(0, 3) = sqrt_info * reduce * jacobian_orientation;
    		}
    	}
    	return true;
    }

    void check(double **parameters)
    {
	    double *res = new double[1];
	    double **jaco = new double *[2];
	    jaco[0] = new double[1 * 7];
		jaco[1] = new double[1 * 7];
	    Evaluate(parameters, res, jaco);
	    puts("check begins");

	    puts("my");

	    std::cout << Eigen::Map<Eigen::Matrix<double, 1, 1>>(res).transpose() << std::endl
	              << std::endl;
	    std::cout << Eigen::Map<Eigen::Matrix<double, 1, 7, Eigen::RowMajor>>(jaco[0]) << std::endl
	              << std::endl;
		std::cout << Eigen::Map<Eigen::Matrix<double, 1, 7, Eigen::RowMajor>>(jaco[1]) << std::endl
	              << std::endl;

		Eigen::Quaterniond Q_i(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
		Eigen::Quaterniond Q_j(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
		Eigen::Quaterniond Q_ij = Q_i.inverse() * Q_j;
		Eigen::Vector3d yaw_ij = Q_ij * rel_yaw;
    	Eigen::Matrix<double, 1, 1> residual;
    	residual.block<1, 1>(0, 0) = yaw_ij.segment(1, 1);
    	residual = sqrt_info * residual;

	    puts("num");
	    std::cout << residual.transpose() << std::endl;

	    const double eps = 1e-6;
	    Eigen::Matrix<double, 1, 3> num_jacobian_i;
	    for (int k = 0; k < 3; k++)
	    {
	    	Eigen::Quaterniond Q_i_new(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

	        int b = k % 3;
	        Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;

			Q_i_new = Q_i_new * Utility::deltaQ(delta);
			Eigen::Quaterniond Q_ij_new = Q_i_new.inverse() * Q_j;
			yaw_ij = Q_ij_new * rel_yaw;
    		Eigen::Matrix<double, 1, 1> tmp_residual;
    		tmp_residual.block<1, 1>(0, 0) = yaw_ij.segment(1, 1);
    		tmp_residual = sqrt_info * tmp_residual;

	        num_jacobian_i.col(k) = (tmp_residual - residual) / eps;
	    }
	    std::cout << num_jacobian_i << std::endl;
		Eigen::Matrix<double, 1, 3> num_jacobian_j;
	    for (int k = 0; k < 3; k++)
	    {
	    	Eigen::Quaterniond Q_j_new(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

	        int b = k % 3;
	        Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;

			Q_j_new = Q_j_new * Utility::deltaQ(delta);
			Eigen::Quaterniond Q_ij_new = Q_i.inverse() * Q_j_new;
			yaw_ij = Q_ij_new * rel_yaw;
    		Eigen::Matrix<double, 1, 1> tmp_residual;
    		tmp_residual.block<1, 1>(0, 0) = yaw_ij.segment(1, 1);
    		tmp_residual = sqrt_info * tmp_residual;

	        num_jacobian_j.col(k) = (tmp_residual - residual) / eps;
	    }
	    std::cout << num_jacobian_j << std::endl;
    }

    Eigen::Vector3d rel_yaw;
    Eigen::Matrix<double, 1, 1> sqrt_info;
};
