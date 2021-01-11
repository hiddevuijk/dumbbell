#ifndef GUARD_VFIELD_H
#define GUARD_VFIELD_H

#include "xy.h"

#include <iostream>
#include <vector>

class VfieldAll {
public:
	virtual double f(const XY& r) = 0; 

protected:
	double w;
	double L;

};

class VNone: public VfieldAll {
public:
	VNone() {w = 0; L = 0;}
	double f(const XY& r)
			{return 1;}
};


class VSine: public VfieldAll {
public:
	VSine() {L = 0; }
	VSine(double ww) { w=ww; }

	double f(const XY& r) {
			return std::sin(w*r.x); }

};

class VSine2: public VfieldAll {
public:
	VSine2() {L = 0; }
	VSine2(double ww) { w=ww; }

	double f(const XY& r) {
			return 1+std::sin(w*r.x); }

};

class VGauss: public VfieldAll {
public:
	VGauss() {L = 0; }
	VGauss(double ww, double LL) {L = LL; w=ww; }

	double f(const XY& r)
    {
		return std::exp(-w*(r.x-0.5*L)*(r.x-0.5*L) ) - std::exp(-w*L*L/4);
    }

};


class VStep: public VfieldAll {
public:
	VStep() {L = 0; }
	VStep(double ww, double LL) {L = LL; w=ww; }

	double f(const XY& r)
    {
		return 0.5*(std::tanh( w*( r.x - 0.25*L ) )*std::tanh( w*(0.75*L - r.x ) ) - std::tanh( -0.25*w*L)*std::tanh(0.75*w*L) );
    }

};


class Vfield {
public:
	Vfield();

	Vfield(double w, double L,
			std::string VType);

	double get_field(const XY& r) const
		{return vfield_ptr->f(r);}

private:
	VNone vnone;
	VSine vsine;
	VSine2 vsine2;
    VGauss vgauss;
    VStep vstep;
	
	VfieldAll *vfield_ptr;
};



Vfield::Vfield()
: vnone(), vsine(), vsine2(), vgauss(), vstep()
{
	vfield_ptr = &vnone;
}


Vfield::Vfield(double w, double L,
			std::string VType)
:
	vnone(), vsine(w), vsine2(w), vgauss(w,L), vstep(w,L)
{
	if( VType == "none") {
		vfield_ptr = &vnone;	
	} else if( VType == "sine" ) {
		vfield_ptr = &vsine;
	} else if( VType == "sine_shift" ) {
		vfield_ptr = &vsine2;
	} else if( VType == "gauss" ) {
		vfield_ptr = &vgauss;
	} else if( VType == "step" ) {
		vfield_ptr = &vstep;
	} else {
		std::cerr << "ERROR: " << VType
			<< " is not a valid option.\n";
	}


}

#endif

