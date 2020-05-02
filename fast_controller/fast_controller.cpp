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

	collective_thrust = mass_ * wzb.dot(a_des + gravity_*e3_ + v_drag);

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
	v_drag = gamma * v / mass_;

	return v_drag;
}

double FastController::computeCommandThrust(const Eigen::Quaterniond q, const Eigen::Vector3d v, 
	double collective_thrust)
{
	// wzb is body z-axis in world frame
	Eigen::Vector3d wzb = q.toRotationMatrix()*e3_;

	double command_thrust = (collective_thrust - rotor_count_ * b_ * cd1_ * v.dot(wzb)) / ( 1. + k_ * cd1_ * v.dot(wzb)); 
	return command_thrust;
}

Eigen::Vector3d FastController::computeCollectiveThrustVector(const Eigen::Vector3d a_des, const Eigen::Vector3d v, 
	double collective_thrust)
{
	Eigen::Vector3d thrust_vector, v_drag = Eigen::Vector3d::Zero();
	v_drag = computeVelocityDrag(v, collective_thrust);
	thrust_vector = mass_ * (a_des + gravity_*e3_ + v_drag);
	return thrust_vector;
}


