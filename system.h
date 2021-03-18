#ifndef GUARD_SYSTEM_H
#define GUARD_SYSTEM_H



#include "ConfigFile.h"

#include "xy.h"
#include "vfield.h"

#include <iostream>
#include <vector>
#include <math.h>
#include <boost/random.hpp>

namespace system_func {
template<class RNDIST> 
inline void xy_random_normal(XY &r, RNDIST &rndist) {
	r.x = rndist();
	r.y = rndist();
}

};

struct System {
public:
	// initialize from ConfigFile object
	System(ConfigFile config);

	// random number generator
	const boost::normal_distribution<double> ndist;
	const boost::uniform_real<double> udist;

	int seed;
	boost::mt19937 rng;		
	boost::variate_generator<boost::mt19937&,
		boost::normal_distribution<double> > rndist;
	boost::variate_generator<boost::mt19937&,
		boost::uniform_real<double> > rudist;

	// fixed system parameters
	unsigned int N;
	double L;
	double v0;

	double T;
	double gamma;
	double gamma_theta;
    double A1xx, A1xy, A1yx, A1yy;
    double A2xx, A2xy, A2yx, A2yy;
    double l;
	double dt;


	double Dt;
	double Dtheta;
	double sqrt_2dt;
	double sqrt_2dt_Dt;
	double sqrt_2dt_Dtheta;

	double pi2;

	// state of the system
	double t;
	std::vector<XY> r;
    std::vector<XY> dr;
    std::vector<XY> v;
	std::vector<double> theta;
	
	Vfield vfield;



	// initialize with random coordinates.
	void init_random();
	void init_random(double);

	// increment time
	void step();

	bool check_x_in_box();
	bool check_y_in_box();

	void write(const char* outname);
	
	// temporary containers
	double Vri;
	XY xi;	
	double eta;


};

void System::step()
{
	double Fx, Fy, Ftheta;
    double px,py, p1x, p1y, p2x, p2y;
    double drx, dry;
    double f1, f2;
	for(unsigned int i=0;i<N;++i) {
        // FRitction!!!

        px = cos(theta[i]);
        py = sin(theta[i]);

        XY p(px,py);

        f1 = v0*vfield.get_field( r[i] + l*p/2);
        f2 = v0*vfield.get_field( r[i] - l*p/2);

        p1x = A1xx*px + A1xy*py;
        p1y = A1yx*px + A1yy*py;
        p2x = A2xx*px + A2xy*py;
        p2y = A2yx*px + A2yy*py;

        Fx  = f1*p1x + f2*p2x;
        Fy  = f1*p1y + f2*p2y;


		system_func::xy_random_normal(xi,rndist);
		xi *= std::sqrt(T*dt/gamma);

        drx = dt*Fx/(2*gamma) + xi.x;
        dry = dt*Fy/(2*gamma) + xi.y;

        v[i].x = ( dr[i].x + drx ) / (2*dt);
        v[i].y = ( dr[i].y + dry ) / (2*dt);
    
        dr[i].x = drx;
        dr[i].y = dry;

        r[i].x += drx;
        r[i].y += dry;
	    //r[i].pbc(L);	


		Ftheta = ( f1*A1yx - f2*A2yx );

		theta[i] += dt*Ftheta/(gamma*l)+ std::sqrt(4*T*dt/(gamma*l*l))*rndist();
		theta[i] = std::fmod( theta[i] , pi2 );

		

	}

	t += dt;
}



System::System(ConfigFile config)
:
	ndist(0.,1.),udist(0,1),
	seed(config.read<unsigned int>("seed")),
	rng(seed), rndist(rng,ndist), rudist(rng,udist),
	vfield( config.read<double>("w"),
			config.read<double>("L"),
			config.read<std::string>("VType") )
{
	XY rr;
    

	// init parameters
	N = config.read<unsigned int>("N");
	L = config.read<double>("L");
	v0 = config.read<double>("v0");
    l = config.read<double>("l");
	T = config.read<double>("T");
	gamma = config.read<double>("gamma");
	gamma_theta = config.read<double>("gamma_theta");

    double phi = config.read<double>("phi1");
    A1xx = std::cos(phi);
    A1yy = std::cos(phi);
    A1xy = -std::sin(phi);
    A1yx = std::sin(phi);
    phi = config.read<double>("phi2");
    A2xx = std::cos(phi);
    A2yy = std::cos(phi);
    A2xy = -std::sin(phi);
    A2yx = std::sin(phi);



	dt = config.read<double>("dt");

	Dt = T/gamma;
	Dtheta = T/gamma_theta;
	sqrt_2dt = std::sqrt(2*dt);
	sqrt_2dt_Dt = std::sqrt(2*dt*Dt);
	sqrt_2dt_Dtheta = std::sqrt(2*dt*Dtheta);

	pi2 = 2*std::acos(-1);
	
	// init state
	t = 0.0;
	r = std::vector<XY>(N);
	dr = std::vector<XY>(N, XY(0,0) );
	v = std::vector<XY>(N, XY(0,0) );
	theta = std::vector<double>(N);
    for(unsigned int i=0;i<N;++i) {
       r[i].x = rudist()*L; 
       r[i].y = rudist()*L; 

       theta[i] = rudist()*2*std::acos(-1);
    }


}

void System::init_random()
{

	for(unsigned int i=0;i<N; ++i) {
        r[i].x = rudist()*L;
        r[i].y = rudist()*L;
		theta[i] = rudist()*pi2;

	}

}

bool System::check_x_in_box()
{
	XY temp;
	bool inside = true;
	for(unsigned int i=0;i<N;++i){
		temp = r[i];
		if(temp.x < 0 or temp.x >L)
			inside = false;

		if(!inside) break;
	}
	return inside;
}


bool System::check_y_in_box()
{
	XY temp;
	bool inside = true;
	for(unsigned int i=0;i<N;++i){
		temp = r[i];
		if(temp.y < 0 or temp.y>L)
			inside = false;

		if(!inside) break;
	}
	return inside;
}



void System::write(const char* outname)
{
	std::ofstream out;
	out.open(outname);
	XY temp;
	for(unsigned int i=0;i<N;++i) {
		temp = r[i];
		temp.pbc(L);	
		out << temp.x << '\t';
		out << temp.y;
		if(i<(N-1)) out << '\n';
	}

	out.close();

}

class Integration {
public:
	Integration( ConfigFile config)
	{
		Nt_init = config.read<unsigned int>("Nt_init");
		Nt = config.read<unsigned int>("Nt");
		sample_freq = config.read<unsigned int>("sample_freq");
		print_freq = config.read<unsigned int>("print_freq");
		t_unit = config.read<unsigned int>("t_unit");
		dt = config.read<double>("dt");
		bs = config.read<double>("bs");
	}

	unsigned int Nt_init;
	unsigned int Nt;
	unsigned int sample_freq;
	unsigned int print_freq;
	unsigned int t_unit;
	double dt;
	double bs;
};


#endif
