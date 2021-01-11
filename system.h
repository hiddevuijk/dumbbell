#ifndef GUARD_SYSTEM_H
#define GUARD_SYSTEM_H



#include "ConfigFile.h"

#include "xy.h"
#include "vfield.h"

#include <iostream>
#include <vector>
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
	double Dr;
	double dt;

	double sqrt_2dt;
	double sqrt_2dt_Dr;

	// state of the system
	double t;
	std::vector<XY> r;
	std::vector<XY> p;
	
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
	XY xi,eta,dp,dv;	


};

void System::step()
{

	for(unsigned int i=0;i<N;++i) {

	    r[i].pbc(L);	
		

		system_func::xy_random_normal(xi,rndist);
		xi *= sqrt_2dt;


	
		if( v0 > 0) {
			
			system_func::xy_random_normal(eta,rndist);
			eta *= sqrt_2dt_Dr;


			p[i] += dp;
			p[i].normalize();
		}
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
	Dr = config.read<double>("Dr");
	dt = config.read<double>("dt");
	sqrt_2dt = std::sqrt(2*dt);
	sqrt_2dt_Dr = std::sqrt(2*dt*Dr);

	// init state
	t = 0.0;
	r = std::vector<XY>(N);
	p = std::vector<XY>(N);


}

void System::init_random()
{

	for(unsigned int i=0;i<N; ++i) {
        r[i].x = rudist()*L;
        r[i].y = rudist()*L;


        double d;
        XY zeta;
		do {
			zeta.x = 2*rudist() - 1.;
			zeta.y = 2*rudist() - 1.;
			d = zeta.length_sq();	
		} while (d > 1.);
		zeta.normalize();
		p[i] = zeta;
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
