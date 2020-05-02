#ifndef FASTCONTROLLER_V12020_H
#define FASTCONTROLLER_V12020_H

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <linear_algebra_utils/utils.h>

class FastController
{
	public:
		FastController(double mass, double max_thrust, double min_thrust,
			const Eigen::Matrix3d & inertia_tensor,
			const Eigen::Vector3d & Kp,
			const Eigen::Vector3d & Kd,
			const Eigen::Vector3d & Ki,
			double Kr, double gravity, double cd1, int rotor_count);
		Eigen::Vector3d computeDesiredAcceleration(const Eigen::Vector3d p, const  Eigen::Vector3d p_ref, 
			const Eigen::Vector3d v, const	Eigen::Vector3d v_ref, const Eigen::Vector3d a_ref);
	 	double computeCollectiveThrust(const Eigen::Quaterniond q, const Eigen::Vector3d a_des, 
	 		const Eigen::Vector3d v, double thrust_ref);
		double computeCommandThrust(const Eigen::Quaterniond q, const Eigen::Vector3d v, 
			double collective_thrust);
		Eigen::Vector3d computeCollectiveThrustVector(const Eigen::Vector3d a_des, const Eigen::Vector3d v, 
			double collective_thrust);
		Eigen::Matrix3d computeDesiredOrientation(const Eigen::Quaterniond q, 
			const Eigen::Vector3d thrust_vector, double yaw_ref); 
		Eigen::Vector3d computeCommandAcceleration(const Eigen::Quaterniond q, double thrust);
		Eigen::Vector3d computeCollectiveThrustVectorDot(const Eigen::Vector3d v, const	Eigen::Vector3d v_ref,
			const Eigen::Vector3d a, const Eigen::Vector3d a_ref, const Eigen::Vector3d j_ref);
		Eigen::Vector3d computeDesiredAngularVelocity(const Eigen::Vector3d a_cmd,
		 	const Eigen::Vector3d thrust_vector, const Eigen::Vector3d thrust_vector_dot, double yaw_ref,
		 	double yaw_dot_ref);

	private:
		// vehicle parameters
		double mass_;
		Eigen::Matrix3d inertia_tensor_;
		double max_thrust_;
		double min_thrust_;
		double gravity_;
		int rotor_count_;

		// control parameters
		double k_; // thrust linearization slop
		double b_; // thrust linearization intercept
		Eigen::Vector3d Kp_ = Eigen::Vector3d::Zero();
		Eigen::Vector3d Kd_ = Eigen::Vector3d::Zero();
		Eigen::Vector3d Ki_ = Eigen::Vector3d::Zero();
 		double Kr_; 
 		double cd1_;

 		// basis vector generating R^3
 		Eigen::Vector3d e1_;
 		Eigen::Vector3d e2_;
 		Eigen::Vector3d e3_;

 		Eigen::Vector3d computeVelocityDrag(const Eigen::Vector3d v, double thrust);

 		inline double getThrustLinearizationSlope(double thrust)
 		{
 			double k = 1./(2. * std::sqrt(thrust));
 			return k;
 		}
 		inline double getThrustLinearizationIntercept(double thrust)
 		{
 			double b = thrust / (2. * std::sqrt(thrust));
 			return b;
 		}





};

#endif
