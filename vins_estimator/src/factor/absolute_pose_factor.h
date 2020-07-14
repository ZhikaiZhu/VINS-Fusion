#pragma once

#include <ros/assert.h>
#include <ceres/ceres.h>
#include <Eigen/Dense>
#include "../utility/utility.h"
#include "../utility/tic_toc.h"
#include "../estimator/parameters.h"

struct RPFactor {
  double Header_i;

  Eigen::Quaterniond zrp;
  Eigen::Matrix2d cov_inv;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class AbsPositionFactor : public ceres::SizedCostFunction<3, 7>
{
  public:
    AbsPositionFactor(const Eigen::Vector3d &_abs_P, const Eigen::Matrix<double, 3, 3> _sqrt_info) : abs_P(_abs_P), sqrt_info(_sqrt_info) {}

    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
    {
    	Eigen::Vector3d P_i(parameters[0][0], parameters[0][1], parameters[0][2]);

    	Eigen::Map<Eigen::Matrix<double, 3, 1>> residual(residuals);
    	residual.block<3, 1>(0, 0) = P_i - abs_P;
    	residual = sqrt_info * residual;

    	if (jacobians)
    	{
    		if (jacobians[0])
    		{
    		    Eigen::Map<Eigen::Matrix<double, 3, 7, Eigen::RowMajor>> jacobian_position_i(jacobians[0]);
    		    jacobian_position_i.setZero();
    		    jacobian_position_i.block<3, 3>(0, 0) = Eigen::Matrix<double, 3, 3>::Identity();
    		    jacobian_position_i = sqrt_info * jacobian_position_i;
    		}
    	}
    	return true;
    }

    void check(double **parameters)
    {
	    double *res = new double[3];
	    double **jaco = new double *[1];
	    jaco[0] = new double[3 * 7];
	    Evaluate(parameters, res, jaco);
	    puts("check begins");

	    puts("my");

	    std::cout << Eigen::Map<Eigen::Matrix<double, 3, 1>>(res).transpose() << std::endl
	              << std::endl;
	    std::cout << Eigen::Map<Eigen::Matrix<double, 3, 7, Eigen::RowMajor>>(jaco[0]) << std::endl
	              << std::endl;

		Eigen::Vector3d P_i(parameters[0][0], parameters[0][1], parameters[0][2]);

    	Eigen::Matrix<double, 3, 1> residual;
    	residual.block<3, 1>(0, 0) = P_i - abs_P;
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
    		tmp_residual.block<3, 1>(0, 0) = P_i_new - abs_P;
    		tmp_residual = sqrt_info * tmp_residual;

	        num_jacobian_i.col(k) = (tmp_residual - residual) / eps;
	    }
	    std::cout << num_jacobian_i << std::endl;
    }

    Eigen::Vector3d abs_P;
    Eigen::Matrix<double, 3, 3> sqrt_info;
};

class RollPitchFactor : public ceres::SizedCostFunction<2, 7>
{
  public:
    RollPitchFactor(const Eigen::Quaterniond &_zrp, const Eigen::Matrix<double, 2, 2> _sqrt_info) : zrp(_zrp), sqrt_info(_sqrt_info) {}

    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
    {
    	Eigen::Quaterniond Q_i(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
		Eigen::Vector3d g_i = Q_i.inverse() * (-Eigen::Vector3d::UnitZ());
		Eigen::Vector3d g_w = zrp * g_i;
    	Eigen::Map<Eigen::Matrix<double, 2, 1>> residual(residuals);
    	residual.block<2, 1>(0, 0) = g_w.segment(0, 2);
    	residual = sqrt_info * residual;

    	if (jacobians)
    	{
    		if (jacobians[0])
    		{
    		    Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_orientation_i(jacobians[0]);
    		    jacobian_orientation_i.setZero();
				Eigen::Matrix<double, 3, 3> jacobian_orientation = zrp.toRotationMatrix() * Utility::skewSymmetric(g_i);
        		Eigen::Matrix<double, 2, 3> reduce(2, 3);
				reduce << 1.0, 0.0, 0.0,
				          0.0, 1.0, 0.0;
				jacobian_orientation_i.block<2, 3>(0, 3) = sqrt_info * reduce * jacobian_orientation;
    		}
    	}
    	return true;
    }

    void check(double **parameters)
    {
	    double *res = new double[2];
	    double **jaco = new double *[1];
	    jaco[0] = new double[2 * 7];
	    Evaluate(parameters, res, jaco);
	    puts("check begins");

	    puts("my");

	    std::cout << Eigen::Map<Eigen::Matrix<double, 2, 1>>(res).transpose() << std::endl
	              << std::endl;
	    std::cout << Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>>(jaco[0]) << std::endl
	              << std::endl;

		Eigen::Quaterniond Q_i(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
		Eigen::Vector3d g_i = Q_i.inverse() * (-Eigen::Vector3d::UnitZ());
		Eigen::Vector3d g_w = zrp * g_i;
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
			g_i = Q_i_new.inverse() * (-Eigen::Vector3d::UnitZ());
			g_w = zrp * g_i;
    		Eigen::Matrix<double, 2, 1> tmp_residual;
    		tmp_residual.block<2, 1>(0, 0) = g_w.segment(0, 2);
    		tmp_residual = sqrt_info * tmp_residual;

	        num_jacobian_i.col(k) = (tmp_residual - residual) / eps;
	    }
	    std::cout << num_jacobian_i << std::endl;
    }

    Eigen::Quaterniond zrp;
    Eigen::Matrix<double, 2, 2> sqrt_info;
};

class YawFactor : public ceres::SizedCostFunction<1, 7>
{
  public:
    YawFactor(const Eigen::Vector3d &_z_yaw, const Eigen::Matrix<double, 1, 1> _sqrt_info) : z_yaw(_z_yaw), sqrt_info(_sqrt_info) {}

    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
    {
    	Eigen::Quaterniond Q_i(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
		Eigen::Vector3d yaw_i = Q_i * z_yaw;
    	Eigen::Map<Eigen::Matrix<double, 1, 1>> residual(residuals);
    	residual.block<1, 1>(0, 0) = yaw_i.segment(1, 1);
    	residual = sqrt_info * residual;

    	if (jacobians)
    	{
    		if (jacobians[0])
    		{
    		    Eigen::Map<Eigen::Matrix<double, 1, 7, Eigen::RowMajor>> jacobian_orientation_i(jacobians[0]);
    		    jacobian_orientation_i.setZero();
				Eigen::Matrix<double, 3, 3> jacobian_orientation = -Q_i.toRotationMatrix() * Utility::skewSymmetric(z_yaw);
        		Eigen::Matrix<double, 1, 3> reduce(1, 3);
				reduce << 0.0, 1.0, 0.0;
				jacobian_orientation_i.block<1, 3>(0, 3) = sqrt_info * reduce * jacobian_orientation;
    		}
    	}
    	return true;
    }

    void check(double **parameters)
    {
	    double *res = new double[1];
	    double **jaco = new double *[1];
	    jaco[0] = new double[1 * 7];
	    Evaluate(parameters, res, jaco);
	    puts("check begins");

	    puts("my");

	    std::cout << Eigen::Map<Eigen::Matrix<double, 1, 1>>(res).transpose() << std::endl
	              << std::endl;
	    std::cout << Eigen::Map<Eigen::Matrix<double, 1, 7, Eigen::RowMajor>>(jaco[0]) << std::endl
	              << std::endl;

		Eigen::Quaterniond Q_i(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
		Eigen::Vector3d yaw_i = Q_i * z_yaw;
    	Eigen::Matrix<double, 1, 1> residual;
    	residual.block<1, 1>(0, 0) = yaw_i.segment(1, 1);
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
			yaw_i = Q_i_new * z_yaw;
    		Eigen::Matrix<double, 1, 1> tmp_residual;
    		tmp_residual.block<1, 1>(0, 0) = yaw_i.segment(1, 1);
    		tmp_residual = sqrt_info * tmp_residual;

	        num_jacobian_i.col(k) = (tmp_residual - residual) / eps;
	    }
	    std::cout << num_jacobian_i << std::endl;
    }

    Eigen::Vector3d z_yaw;
    Eigen::Matrix<double, 1, 1> sqrt_info;
};