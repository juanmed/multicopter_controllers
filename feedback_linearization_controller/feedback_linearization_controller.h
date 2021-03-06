#ifndef FEEDBACKLINEARIZATIONCONTROLLER_V12020_H
#define FEEDBACKLINEARIZATIONCONTROLLER_V12020_H

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <linear_algebra_utils/utils.h>

class FeedbackLinearizationController
{
	public:

		FeedbackLinearizationController(double mass, double max_thrust, double min_thrust,
																		const Eigen::Matrix3d & inertia_tensor,
																		const Eigen::Vector3d & Kp,
																		const Eigen::Vector3d & Kd,
																		const Eigen::Vector3d & Ki,
																		double Kr, double gravity = 9.8);

		Eigen::Vector3d computeDesiredAcceleration(const Eigen::Vector3d p, const  Eigen::Vector3d p_ref, 
																							 const Eigen::Vector3d v, const Eigen::Vector3d v_ref,
																							 const Eigen::Vector3d a_ref);
		Eigen::Vector3d computeDesiredThrustVector(const Eigen::Quaterniond q, const Eigen::Vector3d a_des);
		Eigen::Matrix3d computeDesiredOrientation(const Eigen::Vector3d thrust_vector, double yaw_ref);
		Eigen::Vector3d computeDesiredAngularVelocity(const Eigen::Matrix3d Rbw, 
																									const Eigen::Matrix3d Rbw_des, 
																									const Eigen::Vector3d euler_dot_ref);
		Eigen::Vector3d computeDesiredTorque(const Eigen::Vector3d angular_velocity, const Eigen::Vector3d angular_velocity_ref,
																				 const Eigen::Vector3d torque_ref);
		Eigen::Vector3d computeDesiredTorque2(const Eigen::Vector3d angular_velocity, const Eigen::Vector3d angular_velocity_ref,
																				 const Eigen::Vector3d angular_velocity_dot_ref);
		Eigen::Vector4d computeRotorRPM(double thrust, const Eigen::Vector3d torque, const Eigen::Matrix4d mixer_matrix_inv); 	

	private:

		// vehicle parameters
		double mass_;
		Eigen::Matrix3d inertia_tensor_;
		double max_thrust_;
		double min_thrust_;		
		double gravity_;

		// PID control diagonal gain matrices
		Eigen::Vector3d Kp_;
		Eigen::Vector3d Kd_;
		Eigen::Vector3d Ki_;
 		double Kr_; 

 		// basis vector generating R^3
 		Eigen::Vector3d e1_;
 		Eigen::Vector3d e2_;
 		Eigen::Vector3d e3_;


 		//Eigen::Vector3d	matrixToEulerZYX(const Eigen::Matrix3d R);
};

#endif
	