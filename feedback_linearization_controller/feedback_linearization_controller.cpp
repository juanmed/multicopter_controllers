#include "feedback_linearization_controller.h"

FeedbackLinearizationController::FeedbackLinearizationController(double mass, double max_thrust, double min_thrust, 
																																	const Eigen::Matrix3d & inertia_tensor,
																																	const Eigen::Vector3d & Kp,
																																	const Eigen::Vector3d & Kd,
																																	const Eigen::Vector3d & Ki,
																																	double Kr, double gravity)
	{
		mass_ = mass;
		max_thrust_ = max_thrust;
		min_thrust_ = min_thrust;
		inertia_tensor_ = inertia_tensor;
		gravity_ = gravity;

		Kp_ = Kp;
		Kd_ = Kd;
		Ki_ = Ki;
		Kr_ = Kr;

		
		std::cout << "Initializing Feedback Linearization Controller" << std::endl;
		std::cout << "mass: " << mass_ << std::endl;
		std::cout << "Inertia Tensor: \n" << inertia_tensor_ << std::endl;
		std::cout << "Proportional Gains: \n" << Kp_ << std::endl;
		std::cout << "Derivative Gains: \n" << Kd_ << std::endl;
		std::cout << "Integral Gains: \n" << Ki_ << std::endl;
		std::cout << "Rotation Gain: \n" << Kr_ << std::endl;
		

		e1_ << 1., 0., 0.;
		e2_ << 0., 1., 0.;
		e3_ << 0., 0., 1.;
	}

Eigen::Vector3d FeedbackLinearizationController::computeDesiredAcceleration(const Eigen::Vector3d p, const  Eigen::Vector3d p_ref, 
																																						const Eigen::Vector3d v, const Eigen::Vector3d v_ref,
																																						const Eigen::Vector3d a_ref)
{
	Eigen::Vector3d a_des = Eigen::Vector3d::Zero();
	
	// compute feedback terms
	a_des = -1.0*Kp_.asDiagonal()*(p - p_ref) -1.0*Kd_.asDiagonal()*(v - v_ref); 

	//@TODO implement anti-windup integral control

	// add feedforward terms and grativy
	a_des += a_ref + gravity_*e3_;

	return a_des;
}

Eigen::Vector3d FeedbackLinearizationController::computeDesiredThrustVector(const Eigen::Quaterniond q, const Eigen::Vector3d a_des)
{
	Eigen::Vector3d wzb, wzb_des;
	double thrust;

	// wzb is body z-axis in world frame
	wzb = q.toRotationMatrix()*e3_;
	thrust = mass_ * wzb.dot(a_des);
	wzb_des = a_des.normalized(); /// a_des.norm();

	//@TODO implement thrust value saturation
	return thrust*wzb_des;
}

Eigen::Matrix3d FeedbackLinearizationController::computeDesiredOrientation(const Eigen::Vector3d thrust_vector, double yaw_ref)
{
	Eigen::Vector3d wzb_des, wyc_des, wxb_des, wyb_des;
	Eigen::Matrix3d Rbw_des;

	wzb_des = thrust_vector.normalized();

	// wyc_des is intermediate frame projection in world frame
	wyc_des << -1.0*std::sin(yaw_ref), std::cos(yaw_ref), 0.0;
	wxb_des = wyc_des.cross(wzb_des); //.normalized();
	wxb_des = wxb_des / wxb_des.norm();
	wyb_des = wzb_des.cross(wxb_des); //.normalized();
	wyb_des = wyb_des / wyb_des.norm();

  Rbw_des << wxb_des(0), wyb_des(0), wzb_des(0),
           	 wxb_des(1), wyb_des(1), wzb_des(1),
             wxb_des(2), wyb_des(2), wzb_des(2);

	return Rbw_des;
}

Eigen::Vector3d FeedbackLinearizationController::computeDesiredAngularVelocity(const Eigen::Matrix3d Rbw, 
																																							 const Eigen::Matrix3d Rbw_des, 
																																							 const Eigen::Vector3d euler_dot_ref)
{
	Eigen::Vector3d gains, euler, euler_des, euler_dot, angular_velocity, u;
	double roll, pitch, yaw;
	// The Q matrix maps from body frame angular velocities to world frame eugler angle velocities
	Eigen::Matrix3d Q; 

  gains << Kr_, Kr_, Kr_;
  euler = matrixToEulerZYX(Rbw); //Rbw.eulerAngles(2,1,0);
  euler_des = matrixToEulerZYX(Rbw_des); //Rbw_des.eulerAngles(2,1,0);

  u = -1.0 * gains.asDiagonal() * (euler - euler_des);
  euler_dot = u + euler_dot_ref;

  /* compute w_b angular velocity commands as
     w_b = Q.inv * uc
     where  (euler dot) = K*(angular_velocity)
     Q is -not- a gain matrix, see definition below */
  roll = euler(0);
  pitch = euler(1);
  yaw = euler(2);
  Q << 1.0, std::sin(roll) * std::tan(pitch), std::cos(roll) * std::tan(pitch),
       0.0, std::cos(roll), -1.0*std::sin(roll), 
       0.0, std::sin(roll)/std::cos(pitch), std::cos(roll)/std::cos(pitch);     
  angular_velocity = Q.inverse() * euler_dot;

  return angular_velocity;	
}

Eigen::Vector3d FeedbackLinearizationController::computeDesiredTorque(const Eigen::Vector3d angular_velocity, 
	 																	 const Eigen::Vector3d angular_velocity_ref,
																		 const Eigen::Vector3d torque_ref)
{
	Eigen::Vector3d torque, gains, angular_velocity_error;
	gains << Kr_, Kr_, Kr_;
	angular_velocity_error = angular_velocity - angular_velocity_ref;
	torque = - inertia_tensor_ * gains.asDiagonal() * angular_velocity_error 
					 - angular_velocity_ref.cross(inertia_tensor_ * angular_velocity_ref)
					 + torque_ref
					 + angular_velocity.cross(inertia_tensor_ * angular_velocity);

	//@TODO Torque saturation
	return torque;
}

Eigen::Vector3d FeedbackLinearizationController::computeDesiredTorque2(const Eigen::Vector3d angular_velocity, 
	 																	 const Eigen::Vector3d angular_velocity_ref,
																		 const Eigen::Vector3d angular_velocity_dot_ref)
{
	Eigen::Vector3d torque, gains, angular_velocity_error, input;
	gains << Kr_, Kr_, Kr_;
	angular_velocity_error = angular_velocity - angular_velocity_ref;

	input = -1.0 * gains.asDiagonal() * angular_velocity_error + angular_velocity_dot_ref;
	torque = inertia_tensor_ * input + angular_velocity.cross( inertia_tensor_ * angular_velocity);
	//@TODO Torque saturation
	return torque;
}

Eigen::Vector4d FeedbackLinearizationController::computeRotorRPM(double thrust, const Eigen::Vector3d torque,
																																 const Eigen::Matrix4d mixer_matrix_inv)
{
	Eigen::Vector4d general_input, rotors_rpm;
	general_input << thrust, torque(0), torque(1), torque(2);
	rotors_rpm = mixer_matrix_inv * general_input;
	// the previous mapping returns rpm^2, then take sqrt of each element
	for (uint i = 0; i < 4; i++){
		rotors_rpm(i) = std::sqrt(rotors_rpm(i));
	}
	return rotors_rpm;
}

Eigen::Vector3d	FeedbackLinearizationController::matrixToEulerZYX(const Eigen::Matrix3d R)
{
	double roll, pitch, yaw;
	Eigen::Vector3d euler_angles;

	pitch = std::asin(-1.0 * R(2, 0));
	roll = std::atan2(R(2, 1) / std::cos(pitch), R(2, 2) / std::cos(pitch));
  yaw = std::atan2(R(1, 0) / std::cos(pitch), R(0, 0) / std::cos(pitch));

  euler_angles << roll, pitch, yaw;
  return euler_angles;
}