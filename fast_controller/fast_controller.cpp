#include "fast_controller.h"

FastController::FastController(double mass, double max_thrust, double min_thrust,
	const Eigen::Matrix3d & inertia_tensor, const Eigen::Vector3d & Kp, 
	const Eigen::Vector3d & Kd, const Eigen::Vector3d & Ki, double Kr, 
	double gravity, double cd1, int rotor_count )
{
	mass_ = mass;
	max_thrust_ = max_thrust;
	min_thrust_ = min_thrust;
	inertia_tensor_ = inertia_tensor;
	gravity_ = gravity;
	cd1_ = cd1;
	rotor_count_ = rotor_count;

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

Eigen::Vector3d FastController::computeDesiredAcceleration(const Eigen::Vector3d p, const  Eigen::Vector3d p_ref, 
	const Eigen::Vector3d v, const	Eigen::Vector3d v_ref, const Eigen::Vector3d a_ref)
{
	Eigen::Vector3d a_des = Eigen::Vector3d::Zero();
	
	// compute feedback terms
	a_des = -1.0*Kp_.asDiagonal()*(p - p_ref) -1.0*Kd_.asDiagonal()*(v - v_ref) + a_ref; 

	//@TODO implement anti-windup integral control

	return a_des;
}

double FastController::computeCollectiveThrust(const Eigen::Quaterniond q, const Eigen::Vector3d a_des, 
	const Eigen::Vector3d v, double thrust_ref)
{
	Eigen::Vector3d wzb, wzb_des, v_drag = Eigen::Vector3d::Zero();
	double collective_thrust = 0;

	// calculate velocity dependent drag
	v_drag = computeVelocityDrag(v, thrust_ref);
	// wzb is body z-axis in world frame
	wzb = q.toRotationMatrix()*e3_;

	collective_thrust = mass_ * wzb.dot(a_des + gravity_*e3_ + (v_drag / mass_));

	return collective_thrust;
}

Eigen::Vector3d FastController::computeVelocityDrag(const Eigen::Vector3d v, double thrust)
{

	Eigen::Vector3d v_drag = Eigen::Vector3d::Zero();
	double gamma = 0;
	// calculate velocity dependent drag
	k_ = getThrustLinearizationSlope(thrust);
	b_ = getThrustLinearizationIntercept(thrust);
	gamma = cd1_ * ( thrust * k_ * + rotor_count_ * b_ );
	v_drag = gamma * v ;

	return v_drag;
}

double FastController::computeCommandThrust(const Eigen::Quaterniond q, const Eigen::Vector3d v, 
	double collective_thrust)
{
	// wzb is body z-axis in world frame
	Eigen::Vector3d wzb = q.toRotationMatrix()*e3_;

	double command_thrust = (collective_thrust - (rotor_count_ * b_ * cd1_ * v.dot(wzb)) ) / ( 1. + (k_ * cd1_ * v.dot(wzb)) ); 
	return command_thrust;
}

Eigen::Vector3d FastController::computeCollectiveThrustVector(const Eigen::Vector3d a_des, const Eigen::Vector3d v, 
	double collective_thrust)
{
	Eigen::Vector3d thrust_vector, v_drag = Eigen::Vector3d::Zero();
	// recompute velocity drag with new thrust
	v_drag = computeVelocityDrag(v, collective_thrust); 
	thrust_vector = mass_ * (a_des + gravity_*e3_ + (v_drag / mass_));
	thrust_vector = thrust_vector.normalized();
	return collective_thrust * thrust_vector;
}

Eigen::Matrix3d FastController::computeDesiredOrientation(const Eigen::Quaterniond q, 
	const Eigen::Vector3d thrust_vector, double yaw_ref)
{
	Eigen::Matrix3d Rbw_des1, Rbw_des2 = Eigen::Matrix3d::Zero();
	Eigen::Vector3d wzb_des = thrust_vector.normalized();
	Eigen::Vector3d wyc_des, wxb_des, wyb_des = Eigen::Vector3d::Zero();
	double angle1, angle2 = 0;

	// wyc_des is intermediate frame projection of y_b in world frame
	wyc_des << -1.0*std::sin(yaw_ref), std::cos(yaw_ref), 0.0;
	wxb_des = wyc_des.cross(wzb_des); 
	wxb_des = wxb_des.normalized();
	wyb_des = wzb_des.cross(wxb_des); 
	wyb_des = wyb_des.normalized();	

  Rbw_des1 << wxb_des(0), wyb_des(0), wzb_des(0),
           	  wxb_des(1), wyb_des(1), wzb_des(1),
              wxb_des(2), wyb_des(2), wzb_des(2);	

  Rbw_des2 << -wxb_des(0), -wyb_des(0), wzb_des(0),
           	  -wxb_des(1), -wyb_des(1), wzb_des(1),
              -wxb_des(2), -wyb_des(2), wzb_des(2);

  angle1 = std::fabs(rotation_distance(q.toRotationMatrix(), Rbw_des1));
  angle2 = std::fabs(rotation_distance(q.toRotationMatrix(), Rbw_des2));

  if (angle1 > angle2)
  {
  	return Rbw_des2;
  }
  else
  {
  	return Rbw_des1;
  }
}

Eigen::Vector3d FastController::computeCommandAcceleration(const Eigen::Quaterniond q, const Eigen::Vector3d v,
 double thrust)
{
	// wzb is body z-axis in world frame
	Eigen::Vector3d wzb = q.toRotationMatrix()*e3_;
	Eigen::Vector3d a_cmd, v_drag = Eigen::Vector3d::Zero();
	//@TODO evalute if thrust passed to get velocity drag should be collective or command thrust
	v_drag = computeVelocityDrag(v, thrust); 
	a_cmd = - gravity_ * e3_ 	+ (thrust * wzb / mass_) - (v_drag / mass_);
	
	return a_cmd;
}

Eigen::Vector3d FastController::computeCollectiveThrustVectorDot(const Eigen::Vector3d v, const	Eigen::Vector3d v_ref,
	const Eigen::Vector3d a, const Eigen::Vector3d a_ref, const Eigen::Vector3d j_ref)
{
	Eigen::Vector3d thrust_vector_dot = Eigen::Vector3d::Zero();
	//@TODO include effect of velocity drag
	thrust_vector_dot = -1.0 * Kp_.asDiagonal() * (v - v_ref) -1.0 * Kd_.asDiagonal() * (a - a_ref) + j_ref;
	
	return thrust_vector_dot;
}

Eigen::Vector3d FastController::computeDesiredAngularVelocity(const Eigen::Vector3d thrust_vector, 
	const Eigen::Vector3d thrust_vector_dot, double yaw_ref, double yaw_dot_ref)
{
	//@TODO include effect of velocity drag
	Eigen::Vector3d wyc_des, wxc_des, wxb_des, wzb_des_dot, wxb_des_dot, wyb_des_dot, c, c_dot = Eigen::Vector3d::Zero();
	Eigen::Vector3d angular_velocity = Eigen::Vector3d::Zero();
	Eigen::Vector3d wzb_des = thrust_vector.normalized();
	double thrust_vector_norm, thrust_vector_norm_dot, c_norm, c_norm_dot = 0.0;
	Eigen::Matrix3d Rbw_des_dot = Eigen::Matrix3d::Zero();

	// wyc_des is intermediate frame projection of y_b in world frame
	wyc_des << -1.0*std::sin(yaw_ref), std::cos(yaw_ref), 0.0;
	wxc_des << std::cos(yaw_ref), std::sin(yaw_ref), 0.0;
	wxb_des = wyc_des.cross(wzb_des); 
	wxb_des = wxb_des.normalized();

	// wzb_des_dot
	thrust_vector_norm = thrust_vector.norm();
	thrust_vector_norm_dot = thrust_vector.dot(thrust_vector_dot)/thrust_vector_norm;
	wzb_des_dot = (thrust_vector_norm * thrust_vector_dot - thrust_vector * thrust_vector_norm_dot) / (std::pow(thrust_vector_norm, 2));

	c = wyc_des.cross(wzb_des);
	c_norm = c.norm();
	c_dot = -1.0 * yaw_dot_ref * wxc_des.cross(wzb_des) + wyc_des.cross(wzb_des_dot);
	c_norm_dot = c.dot(c_dot)/c_norm;

	wxb_des_dot = (c_norm * c_dot - c * c_norm_dot) / std::pow(c_norm, 2);

	wyb_des_dot = wzb_des_dot.cross(wxb_des) + wzb_des.cross(wxb_des_dot);

	Rbw_des_dot << wxb_des_dot(0), wyb_des_dot(0), wzb_des_dot(0),
           	  	 wxb_des_dot(1), wyb_des_dot(1), wzb_des_dot(1),
              	 wxb_des_dot(2), wyb_des_dot(2), wzb_des_dot(2);

  angular_velocity = vex2(Rbw_des_dot);

  return angular_velocity;
}

Eigen::Vector3d FastController::computeDesiredTorque(const Eigen::Vector3d angular_velocity, 
	const Eigen::Vector3d angular_velocity_ref, const Eigen::Vector3d angular_velocity_dot_ref)
{
	Eigen::Vector3d torque, gains, angular_velocity_error, input = Eigen::Vector3d::Zero();
	gains << Kr_, Kr_, Kr_;
	angular_velocity_error = angular_velocity - angular_velocity_ref;

	input = -1.0 * gains.asDiagonal() * angular_velocity_error + angular_velocity_dot_ref;
	torque = inertia_tensor_ * input + angular_velocity.cross( inertia_tensor_ * angular_velocity);
	//@TODO Torque saturation
	return torque;	
}

Eigen::Vector4d FastController::computeRotorRPM(double thrust, const Eigen::Vector3d torque, 
	const Eigen::Matrix4d mixer_matrix_inv)
{
	Eigen::Vector4d general_input, rotors_rpm = Eigen::Vector4d::Zero();
	general_input << thrust, torque(0), torque(1), torque(2);
	rotors_rpm = mixer_matrix_inv * general_input;
	// the previous mapping returns rpm^2, then take sqrt of each element
	for (uint i = 0; i < 4; i++){
		rotors_rpm(i) = std::sqrt(rotors_rpm(i));
	}
	return rotors_rpm;
}
