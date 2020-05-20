#include "fast_controller.h"

FastController::FastController(double mass, double max_thrust, double min_thrust,
	const Eigen::Matrix3d & inertia_tensor, const Eigen::Vector3d & Kp, 
	const Eigen::Vector3d & Kd, const Eigen::Vector3d & Ki, double Kr, 
	double Ko, double gravity, double cd1, int rotor_count )
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
	Ko_ = Ko;

	std::cout << "Initializing Fast Controller" << std::endl;
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

void FastController::computeAvalues(const Eigen::Vector3d p, const  Eigen::Vector3d p_ref, 
	const Eigen::Vector3d v, const	Eigen::Vector3d v_ref, const Eigen::Vector3d a, 
	const Eigen::Vector3d a_ref, const Eigen::Vector3d j, const Eigen::Vector3d j_ref, 
	const Eigen::Vector3d s_ref, double thrust_ref)
{
	Eigen::Vector3d v_drag = Eigen::Vector3d::Zero();
	v_drag = computeVelocityDrag(v, thrust_ref);
	A_ = mass_ *  (-1.0*Kp_.asDiagonal()*(p - p_ref) -1.0*Kd_.asDiagonal()*(v - v_ref) + a_ref + gravity_*e3_  + v_drag/mass_ );
	A_norm_ = A_.norm();
	A_dot_ = mass_ * (-1.0 * Kp_.asDiagonal() * (v - v_ref) -1.0 * Kd_.asDiagonal() * (a - a_ref) + j_ref);
	A_norm_dot_ = A_.dot(A_dot_) / A_norm_;
	A_ddot_ = mass_ * (-1.0 * Kp_.asDiagonal() * (a - a_ref) - 1.0 * Kd_.asDiagonal() * (j - j_ref) + s_ref);
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
	double gamma = 0.;
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
	thrust_vector = mass_ * (a_des + gravity_*e3_ ) + v_drag;
	thrust_vector = thrust_vector.normalized();
	return collective_thrust * thrust_vector;
}

double FastController::computeCollectiveThrustDot(const Eigen::Vector3d thrust_vector_dot, const Eigen::Vector3d wzb_dot, 
	const Eigen::Vector3d wzb)
{
	double collective_thrust_dot = wzb_dot.dot(A_) + mass_ * wzb.dot(thrust_vector_dot);
	return collective_thrust_dot;
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

Eigen::Matrix3d FastController::computeDesiredOrientationDot(const Eigen::Vector3d thrust_vector, 
	const Eigen::Vector3d thrust_vector_dot, double yaw_ref, double yaw_dot_ref)
{
	//@TODO include effect of velocity drag
	Eigen::Vector3d wyc_des, wxc_des, wxb_des, wzb_des_dot, wxb_des_dot, wyb_des_dot, wyb_des = Eigen::Vector3d::Zero();
	Eigen::Vector3d angular_velocity = Eigen::Vector3d::Zero();
	Eigen::Vector3d wzb_des = thrust_vector.normalized();
	double thrust_vector_norm, thrust_vector_norm_dot= 0.0;
	Eigen::Matrix3d Rbw_des_dot = Eigen::Matrix3d::Zero();

	// wzb_des_dot
	thrust_vector_norm = thrust_vector.norm();
	thrust_vector_norm_dot = thrust_vector.dot(thrust_vector_dot)/thrust_vector_norm;
	wzb_des_dot = (thrust_vector_norm * thrust_vector_dot - thrust_vector * thrust_vector_norm_dot) / (std::pow(thrust_vector_norm, 2));

	/* IMPLEMENTATION WITH Yc
	// wyc_des is intermediate frame projection of y_b in world frame
	wyc_des << -1.0*std::sin(yaw_ref), std::cos(yaw_ref), 0.0;
	wxc_des << std::cos(yaw_ref), std::sin(yaw_ref), 0.0;
	wxb_des = wyc_des.cross(wzb_des); 
	wxb_des = wxb_des.normalized();

	c = wyc_des.cross(wzb_des);
	c_norm = c.norm();
	c_dot = -1.0 * yaw_dot_ref * wxc_des.cross(wzb_des) + wyc_des.cross(wzb_des_dot);
	c_norm_dot = c.dot(c_dot)/c_norm;

	wxb_des_dot = (c_norm * c_dot - c * c_norm_dot) / std::pow(c_norm, 2);

	wyb_des_dot = wzb_des_dot.cross(wxb_des) + wzb_des.cross(wxb_des_dot);
  */

	/* IMPLEMENTATION WITH Xc 
  wxc_des << std::cos(yaw_ref), std::sin(yaw_ref), 0.0;
  wyc_des << -1.0*std::sin(yaw_ref), std::cos(yaw_ref), 0.0;
  c_ = wzb_des.cross(wxc_des);
  c_norm_ = c_.norm();
  c_dot_ = yaw_dot_ref * wzb_des.cross(wyc_des) + wzb_des_dot.cross(wxc_des);
  c_norm_dot_ = c_.dot(c_dot_) / c_norm_;
	*/

	computeCvalues(wzb_des, wzb_des_dot, yaw_ref, yaw_dot_ref);
  wyb_des = c_ / c_norm_;
  wyb_des_dot = (c_norm_ * c_dot_ - c_ * c_norm_dot_) / std::pow(c_norm_, 2);
  
  wxb_des_dot = wyb_des.cross(wzb_des_dot) + wyb_des_dot.cross(wzb_des);

	Rbw_des_dot << wxb_des_dot(0), wyb_des_dot(0), wzb_des_dot(0),
           	  	 wxb_des_dot(1), wyb_des_dot(1), wzb_des_dot(1),
              	 wxb_des_dot(2), wyb_des_dot(2), wzb_des_dot(2);

  return Rbw_des_dot;
}

Eigen::Matrix3d FastController::computeDesiredOrientationDot2(const Eigen::Vector3d p, const  Eigen::Vector3d p_ref, 
	const Eigen::Vector3d v, const	Eigen::Vector3d v_ref, const Eigen::Vector3d a, 
	const Eigen::Vector3d a_ref, const Eigen::Vector3d j, const Eigen::Vector3d j_ref, 
	const Eigen::Vector3d s_ref, double thrust_ref, double yaw_ref, double yaw_dot_ref)
{
	Eigen::Vector3d wzb_des, wzb_des_dot, wyb_des, wyb_des_dot, wxb_des_dot =Eigen::Vector3d::Zero();
	Eigen::Matrix3d Rbw_des_dot = Eigen::Matrix3d::Zero();

	computeAvalues(p, p_ref, v, v_ref, a, a_ref, j, j_ref, s_ref, thrust_ref);
	wzb_des_dot = (A_norm_ * A_dot_ - A_ * A_norm_dot_) / std::pow(A_norm_, 2);

	wzb_des = A_.normalized();
	computeCvalues(wzb_des, wzb_des_dot, yaw_ref, yaw_dot_ref);
  wyb_des = c_ / c_norm_;
  wyb_des_dot = (c_norm_ * c_dot_ - c_ * c_norm_dot_) / std::pow(c_norm_, 2);

  wxb_des_dot = wyb_des.cross(wzb_des_dot) + wyb_des_dot.cross(wzb_des);

	Rbw_des_dot << wxb_des_dot(0), wyb_des_dot(0), wzb_des_dot(0),
           	  	 wxb_des_dot(1), wyb_des_dot(1), wzb_des_dot(1),
              	 wxb_des_dot(2), wyb_des_dot(2), wzb_des_dot(2);

  return Rbw_des_dot;
}

void FastController::computeCvalues(const Eigen::Vector3d wzb, const Eigen::Vector3d wzb_dot, 
	double yaw, double yaw_dot)
{
  wxc_des_ << std::cos(yaw), std::sin(yaw), 0.0;
  wyc_des_ << -1.0*std::sin(yaw), std::cos(yaw), 0.0;
  c_ = wzb.cross(wxc_des_);
  c_norm_ = c_.norm();
  c_dot_ = yaw_dot * wzb.cross(wyc_des_) + wzb_dot.cross(wxc_des_);
  c_norm_dot_ = c_.dot(c_dot_) / c_norm_;
}

Eigen::Vector3d FastController::computeDesiredAngularVelocity(const Eigen::Matrix3d desired_orientation,
	const Eigen::Matrix3d desired_orientation_dot)
{
	Eigen::Vector3d angular_velocity = Eigen::Vector3d::Zero();
	angular_velocity = vex2(desired_orientation.transpose() * desired_orientation_dot);
	return angular_velocity;
}

Eigen::Vector3d FastController::computeDesiredAngularVelocity2(const Eigen::Matrix3d Rbw, const Eigen::Matrix3d Rbw_des, 
	const Eigen::Vector3d euler_dot_ref)
{
	Eigen::Vector3d gains, euler, euler_des, euler_dot, angular_velocity, u = Eigen::Vector3d::Zero();
	double roll, pitch, yaw = 0.0;
	// The Q matrix maps from body frame angular velocities to world frame eugler angle velocities
	Eigen::Matrix3d Q = Eigen::Matrix3d::Zero(); 

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

Eigen::Vector3d FastController::computeCollectiveThrustVectorDDot(const Eigen::Vector3d a, 
	const Eigen::Vector3d a_ref, const Eigen::Vector3d j, const Eigen::Vector3d j_ref, 
	const Eigen::Vector3d s_ref)
{
	Eigen::Vector3d thrust_vector_ddot = Eigen::Vector3d::Zero();
	thrust_vector_ddot = -1.0 * Kp_.asDiagonal() * (a - a_ref) - 1.0 * Kd_.asDiagonal() * (j - j_ref) + s_ref;
	return thrust_vector_ddot;
}

Eigen::Matrix3d FastController::computeDesiredOrientationDDot(const Eigen::Vector3d thrust_vector, 
	const Eigen::Vector3d thrust_vector_dot, const Eigen::Vector3d thrust_vector_ddot,
	const Eigen::Matrix3d desired_orientation, const Eigen::Matrix3d desired_orientation_dot, double yaw_ref, 
	double yaw_dot_ref, double yaw_ddot_ref)
{
	Eigen::Vector3d wzb_des_ddot, wyb_des_ddot, wxb_des_ddot, wzb_des_dot, wyb_des_dot, wyb_des, wzb_des = Eigen::Vector3d::Zero();
	double thrust_vector_norm, thrust_vector_norm_dot= 0.0;
	Eigen::Matrix3d Rbw_des_ddot = Eigen::Matrix3d::Zero();

	thrust_vector_norm = thrust_vector.norm();
	thrust_vector_norm_dot = thrust_vector.dot(thrust_vector_dot)/thrust_vector_norm;	
	wzb_des_dot = desired_orientation_dot * e3_;
	wzb_des = desired_orientation * e3_;
	wzb_des_ddot =  ( thrust_vector_norm * thrust_vector_ddot - thrust_vector * ( thrust_vector_ddot.dot(wzb_des) + thrust_vector_dot.dot(wzb_des_dot))  - 2.0 * wzb_des_dot * thrust_vector_dot.dot( thrust_vector ) ) / std::pow(thrust_vector_norm, 2);

	computeCvalues(wzb_des, wzb_des_dot, yaw_ref, yaw_dot_ref);

	wyb_des = desired_orientation * e2_;
	wyb_des_dot = desired_orientation_dot * e2_;
	c_ddot_ = yaw_dot_ref * (wzb_des_dot.cross(wyc_des_) - yaw_dot_ref * wzb_des.cross(wxc_des_) ) + yaw_ddot_ref * wzb_des.cross(wyc_des_) + wzb_des_ddot.cross(wxc_des_) + yaw_dot_ref * wzb_des_dot.cross(wyc_des_);
	wyb_des_ddot = ( (c_norm_ * c_ddot_) - c_ * (c_ddot_.dot(wyb_des) + c_dot_.dot(wyb_des_dot) ) - 2.0 * wyb_des_dot * c_dot_.dot(c_) ) / std::pow(c_norm_, 2);

	wxb_des_ddot = 2.0 * wyb_des_dot.cross(wzb_des_dot) + wyb_des.cross(wzb_des_ddot) + wyb_des_ddot.cross(wzb_des);

	Rbw_des_ddot << wxb_des_ddot(0), wyb_des_ddot(0), wzb_des_ddot(0),
           	  	 	wxb_des_ddot(1), wyb_des_ddot(1), wzb_des_ddot(1),
              	 	wxb_des_ddot(2), wyb_des_ddot(2), wzb_des_ddot(2);

  return Rbw_des_ddot;
}

Eigen::Matrix3d FastController::computeDesiredOrientationDDot2(const Eigen::Vector3d p, const  Eigen::Vector3d p_ref, 
	const Eigen::Vector3d v, const	Eigen::Vector3d v_ref, const Eigen::Vector3d a, 
	const Eigen::Vector3d a_ref, const Eigen::Vector3d j, const Eigen::Vector3d j_ref, 
	const Eigen::Vector3d s_ref, const Eigen::Matrix3d desired_orientation, const Eigen::Matrix3d desired_orientation_dot,
	double yaw_ref, double yaw_dot_ref, double yaw_ddot_ref, double thrust_ref)
{
	Eigen::Vector3d wzb_des, wzb_des_dot, wzb_des_ddot, wyb_des, wyb_des_dot, wyb_des_ddot, wxb_des_ddot = Eigen::Vector3d::Zero();
	Eigen::Matrix3d Rbw_des_ddot = Eigen::Matrix3d::Zero();

	computeAvalues(p, p_ref, v, v_ref, a, a_ref, j, j_ref, s_ref, thrust_ref);
	wzb_des = desired_orientation * e3_;
	wzb_des_dot = desired_orientation_dot * e3_;

	wzb_des_ddot = (A_norm_ * A_ddot_ - A_ * (A_ddot_.dot(wzb_des) + A_dot_.dot(wzb_des_dot) ) - 2.0 * wzb_des_dot * A_dot_.dot(A_)  ) / std::pow(A_norm_, 2);

	computeCvalues(wzb_des, wzb_des_dot, yaw_ref, yaw_dot_ref);
	wyb_des = desired_orientation * e2_;
	wyb_des_dot = desired_orientation_dot * e2_;

	wyb_des_ddot = (c_norm_ * c_ddot_ - c_ * (c_ddot_.dot(wyb_des) + c_dot_.dot(wyb_des_dot) ) - 2.0 * wyb_des_dot * c_dot_.dot(c_)  ) / std::pow(c_norm_, 2);

	wxb_des_ddot = 2.0 * wyb_des_dot.cross(wzb_des_dot) + wyb_des.cross(wzb_des_ddot) + wyb_des_ddot.cross(wzb_des);

	Rbw_des_ddot << wxb_des_ddot(0), wyb_des_ddot(0), wzb_des_ddot(0),
           	  	 	wxb_des_ddot(1), wyb_des_ddot(1), wzb_des_ddot(1),
              	 	wxb_des_ddot(2), wyb_des_ddot(2), wzb_des_ddot(2);

  return Rbw_des_ddot;

}

Eigen::Vector3d FastController::computeDesiredAngularVelocityDot(const Eigen::Matrix3d desired_orientation,
	const Eigen::Matrix3d desired_orientation_ddot, const Eigen::Vector3d desired_angular_velocity)
{
	Eigen::Vector3d angular_velocity_dot = Eigen::Vector3d::Zero();
	angular_velocity_dot = vex2( desired_orientation.transpose() * desired_orientation_ddot - hat(desired_angular_velocity) * hat(desired_angular_velocity) );
	return angular_velocity_dot;
}

Eigen::Vector3d FastController::computeDesiredTorque(const Eigen::Vector3d angular_velocity, 
	const Eigen::Vector3d angular_velocity_ref, const Eigen::Vector3d angular_velocity_dot_ref)
{
	Eigen::Vector3d torque, gains, angular_velocity_error, input = Eigen::Vector3d::Zero();
	angular_velocity_error = angular_velocity - angular_velocity_ref;

	input = - Kr_ * angular_velocity_error + angular_velocity_dot_ref;
	torque = inertia_tensor_ * input + angular_velocity.cross( inertia_tensor_ * angular_velocity);
	//@TODO Torque saturation
	return torque;	
}

Eigen::Vector3d FastController::computeDesiredTorque2(const Eigen::Matrix3d orientation, 
	const Eigen::Matrix3d desired_orientation, const Eigen::Vector3d desired_angular_acceleration,
	const Eigen::Vector3d angular_velocity, const Eigen::Vector3d desired_angular_velocity)
{
	Eigen::Vector3d torque, orientation_error, angular_velocity_error, input, coriolis_torque, angular_inertia = Eigen::Vector3d::Zero();
	orientation_error = computeOrientationError(orientation, desired_orientation);
	angular_velocity_error = computeAngularVelocityError(orientation, desired_orientation, angular_velocity, desired_angular_velocity);

	input = -Kr_ * orientation_error - Ko_ * angular_velocity_error;
	coriolis_torque = angular_velocity.cross( inertia_tensor_ * angular_velocity );
	angular_inertia = - inertia_tensor_ * ( hat(angular_velocity) * orientation.transpose() * desired_orientation * desired_angular_velocity - orientation.transpose() * desired_orientation *  desired_angular_acceleration);
	torque = input + coriolis_torque + angular_inertia;
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
		rotors_rpm(i) = clip_scalar(std::sqrt(rotors_rpm(i)), max_thrust_, min_thrust_);
	}
	return rotors_rpm;
}
