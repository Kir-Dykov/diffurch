#pragma once

class DDE1 {
public:
	double b,c,d,tau;
	double h, t_finish;

	inline void eval_f(const double& x, const double& dx, const double& s, double& to_x, double& to_dx ){
		to_dx = -b * dx - c * x + d * s;
		to_x  = dx;
	}

	inline void RKstep(double x, double dx, double x_tau, double h, double & xNext, double & dxNext) {
		double s = sign(x_tau); 

		double dK1 = (-b * dx - c * x + d*s); 
		double  K1 = dx; 

		double dK2 = (-b * (dx + dK1*h*0.5) - c * (x + K1*h*0.5) + d*s);
		double  K2 = (dx +dK1*h*0.5);

		double dK3 = (-b * (dx + dK2*h*0.5) - c * (x + K2*h*0.5) + d*s);
		double  K3 = (dx +dK2*h*0.5);

		double dK4 = (-b * (dx + dK3*h) - c * (x + K3*h) + d*s);
		double  K4 = (dx + dK3*h); 

		dxNext = dx + h * (1./6.) * (dK1 + 2. * dK2 + 2. * dK3 + dK4);
		xNext  =  x + h * (1./6.) * ( K1 + 2. *  K2 + 2. *  K3 +  K4);
	}

	inline void DOPRI5step(double x, double dx, double x_tau, double h, double & xNext, double & dxNext) {
		double s = sign(x_tau); 

		double  k1, k2, k3, k4, k5, k6, k7;
		double dk1,dk2,dk3,dk4,dk5,dk6,dk7;

		eval_f(x, dx, s, k1, dk1);
		eval_f( x + h*(1./5.)* k1,
				dx + h*(1./5.)*dk1,
				s, k2, dk2);
		eval_f( x + h*((3./40.)* k1 + (9./40.)* k2),
				dx + h*((3./40.)*dk1 + (9./40.)*dk2),
				s, k3, dk3);
		eval_f( x + h*((44./45.)* k1 + (-56./15.)* k2 + (32./9.)* k3),
				dx + h*((44./45.)*dk1 + (-56./15.)*dk2 + (32./9.)*dk3),
				s, k4, dk4);

		eval_f( 
				x + h*((19372./6561.)* k1 + (-25360./2187.)* k2 + (64448./6561.)* k3 + (-212./729.)* k4),
				dx + h*((19372./6561.)*dk1 + (-25360./2187.)*dk2 + (64448./6561.)*dk3 + (-212./729.)*dk4),
				s, k5, dk5);

		eval_f( 
				x + h*((9017./3168.)* k1 + (-355./33.)* k2 + (46732./5247.)* k3 + (49./176.)* k4 + (-5103./18656)* k5),
				dx + h*((9017./3168.)*dk1 + (-355./33.)*dk2 + (46732./5247.)*dk3 + (49./176.)*dk4 + (-5103./18656)*dk5),
				s, k6, dk6);


		xNext =  x + h*((35./384.)* k1 + (500./1113.)* k3 + (125./192.)* k4 + (-2187./6784)* k5 + (11./84.)* k6);
		dxNext = dx + h*((35./384.)*dk1 + (500./1113.)*dk3 + (125./192.)*dk4 + (-2187./6784)*dk5 + (11./84.)*dk6);

		eval_f( 
				xNext,
				dxNext,
				s, k7, dk7);

	}

	void step(double x, double dx, double x_tau, double h, double & xNext, double & dxNext) {
		return DOPRI5step(x,dx,x_tau,h,xNext,dxNext);
	}
};
