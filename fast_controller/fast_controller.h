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
			double Kr, double Ko, double gravity, double cd1, int rotor_count);
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
		Eigen::Vector3d computeCommandAcceleration(const Eigen::Quaterniond q, const Eigen::Vector3d v, 
			double thrust);
		Eigen::Vector3d computeCollectiveThrustVectorDot(const Eigen::Vector3d v, const	Eigen::Vector3d v_ref,
			const Eigen::Vector3d a, const Eigen::Vector3d a_ref, const Eigen::Vector3d j_ref);
		double computeCollectiveThrustDot(const Eigen::Vector3d thrust_vector_dot, const Eigen::Vector3d wzb_dot,
			const Eigen::Vector3d wzb);
		Eigen::Matrix3d computeDesiredOrientationDot(const Eigen::Vector3d thrust_vector, 
			const Eigen::Vector3d thrust_vector_dot, double yaw_ref, double yaw_dot_ref);
		Eigen::Matrix3d computeDesiredOrientationDot2(const Eigen::Vector3d p, const  Eigen::Vector3d p_ref, 
			const Eigen::Vector3d v, const	Eigen::Vector3d v_ref, const Eigen::Vector3d a, 
			const Eigen::Vector3d a_ref, const Eigen::Vector3d j, const Eigen::Vector3d j_ref, 
			const Eigen::Vector3d s_ref, double thrust_ref, double yaw_ref, double yaw_dot_ref);
		Eigen::Vector3d computeDesiredAngularVelocity(const Eigen::Matrix3d desired_orientation,
			const Eigen::Matrix3d desired_orientation_dot);
		Eigen::Vector3d computeDesiredAngularVelocity2(const Eigen::Matrix3d Rbw, const Eigen::Matrix3d Rbw_des, 
			const Eigen::Vector3d euler_dot_ref);
		Eigen::Vector3d computeDesiredAngularVelocityDot(const Eigen::Matrix3d desired_orientation,
			const Eigen::Matrix3d desired_orientation_ddot, const Eigen::Vector3d desired_angular_velocity);
		Eigen::Vector3d computeCollectiveThrustVectorDDot(const Eigen::Vector3d a, const Eigen::Vector3d a_ref,
			const Eigen::Vector3d j, const Eigen::Vector3d j_ref, const Eigen::Vector3d s_ref);
		Eigen::Matrix3d computeDesiredOrientationDDot(const Eigen::Vector3d thrust_vector, 
			const Eigen::Vector3d thrust_vector_dot, const Eigen::Vector3d thrust_vector_ddot,
			const Eigen::Matrix3d desired_orientation, const Eigen::Matrix3d desired_orientation_dot, double yaw_ref, 
			double yaw_dot_ref, double yaw_ddot_ref);
		Eigen::Matrix3d computeDesiredOrientationDDot2(const Eigen::Vector3d p, const  Eigen::Vector3d p_ref, 
			const Eigen::Vector3d v, const	Eigen::Vector3d v_ref, const Eigen::Vector3d a, 
			const Eigen::Vector3d a_ref, const Eigen::Vector3d j, const Eigen::Vector3d j_ref, 
			const Eigen::Vector3d s_ref, const Eigen::Matrix3d desired_orientation, const Eigen::Matrix3d desired_orientation_dot,
			double yaw_ref, double yaw_dot_ref, double yaw_ddot_ref, double thrust_ref);
		Eigen::Vector3d computeDesiredTorque(const Eigen::Vector3d angular_velocity, 
			const Eigen::Vector3d angular_velocity_ref, const Eigen::Vector3d angular_velocity_dot_ref);
		Eigen::Vector3d computeDesiredTorque2(const Eigen::Matrix3d orientation, 
			const Eigen::Matrix3d desired_orientation, const Eigen::Vector3d desired_angular_acceleration,
			const Eigen::Vector3d angular_velocity, const Eigen::Vector3d desired_angular_velocity);
		Eigen::Vector3d computeDesiredTorque3(const Eigen::Matrix3d orientation, 
			const Eigen::Matrix3d desired_orientation, const Eigen::Vector3d angular_velocity, 
			const Eigen::Vector3d desired_angular_velocity);
		Eigen::Vector4d computeRotorRPM(double thrust, const Eigen::Vector3d torque, 
			const Eigen::Matrix4d mixer_matrix_inv); 

		// controller parameters
		double k_; // thrust linearization slop
		double b_; // thrust linearization intercept

	private:
		// vehicle parameters
		double mass_;
		Eigen::Matrix3d inertia_tensor_;
		double max_thrust_;
		double min_thrust_;
		double gravity_;
		int rotor_count_;


		Eigen::Vector3d Kp_ = Eigen::Vector3d::Zero();
		Eigen::Vector3d Kd_ = Eigen::Vector3d::Zero();
		Eigen::Vector3d Ki_ = Eigen::Vector3d::Zero();
 		double Kr_, Ko_; 
 		double cd1_;

 		// controller values
 		Eigen::Vector3d A_, A_dot_, A_ddot_ = Eigen::Vector3d::Zero();
 		double A_norm_, A_norm_dot_ = 0.0;
 		Eigen::Vector3d c_, c_dot_, c_ddot_, wxc_des_, wyc_des_ = Eigen::Vector3d::Zero();
 		double c_norm_, c_norm_dot_ = 0.0;


 		// basis vector generating R^3
 		Eigen::Vector3d e1_;
 		Eigen::Vector3d e2_;
 		Eigen::Vector3d e3_;

 		Eigen::Vector3d computeVelocityDrag(const Eigen::Vector3d v, double thrust);

		void computeCvalues(const Eigen::Vector3d wzb, const Eigen::Vector3d wzb_dot, 
			double yaw, double yaw_dot);

		void computeAvalues(const Eigen::Vector3d p, const  Eigen::Vector3d p_ref, 
			const Eigen::Vector3d v, const	Eigen::Vector3d v_ref, const Eigen::Vector3d a, 
			const Eigen::Vector3d a_ref, const Eigen::Vector3d j, const Eigen::Vector3d j_ref, 
			const Eigen::Vector3d s_ref, double thrust_ref);

 		inline double getThrustLinearizationSlope(double thrust)
 		{
 			double k = 1./(2. * std::sqrt(thrust));
 			return 0.2;
 		}
 		inline double getThrustLinearizationIntercept(double thrust)
 		{
 			double b = thrust / (2. * std::sqrt(thrust));
 			return 1.1;
 		}

		inline Eigen::Vector3d computeOrientationError(const Eigen::Matrix3d orientation, 
			const Eigen::Matrix3d desired_orientation)
		{
			Eigen::Vector3d orientation_error = vex2( desired_orientation.transpose() * orientation - orientation.transpose() * desired_orientation ) / 2.;
			return orientation_error;
		}
		inline Eigen::Vector3d computeAngularVelocityError(const Eigen::Matrix3d orientation, 
			const Eigen::Matrix3d desired_orientation, const Eigen::Vector3d angular_velocity,
			const Eigen::Vector3d desired_angular_velocity)
		{
			Eigen::Vector3d angular_velocity_error = angular_velocity - (orientation.transpose() * desired_orientation * desired_angular_velocity);
			return angular_velocity_error;
		}

};

#endif /* FASTCONTROLLER_V12020_H */ 
