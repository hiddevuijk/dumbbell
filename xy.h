#ifndef	GUARD_XY_H
#define	GUARD_XY_H

#include <cmath>
#include <ostream>

class XY {
public:
	XY(): x(0.), y(0.) {}
	XY(double xx, double yy)
		: x(xx), y(yy){}

	double x,y;

	double length() const
		{ return std::sqrt(x*x + y*y ); }

	double length_sq() const
		{ return x*x + y*y; } 

	void normalize( double d = 1.)
	{	double len = std::sqrt(x*x+y*y);
		x *= d/len;
		y *= d/len;
	}

	void pbc(double L) {
		x -= L*std::floor(x/L);
		y -= L*std::floor(y/L);
	}

	// addition, subtraction, multiplication and division	
	XY operator+=(const XY &r) {
		x += r.x;
		y += r.y;
		return *this;
	}

	XY operator+=(const double & add) {
		x += add;
		y += add;
		return *this;
	}

	XY operator-=(const XY &r) {
		x -= r.x;
		y -= r.y;
		return *this;
	}

	XY operator-=(const double & sub) {
		x -= sub;
		y -= sub;
		return *this;
	}

	XY operator*=(const XY &r) {
		x *= r.x;
		y *= r.y;
		return *this;
	}

	XY operator*=(const double & mult) {
		x *= mult;
		y *= mult;
		return *this;
	}

	XY operator/=(const XY &r) {
		x /= r.x;
		y /= r.y;
		return *this;
	}

	XY operator/=(const double & div) {
		x /= div;
		y /= div;
		return *this;
	}
};


std::ostream& operator<< (std::ostream &out, XY const& r)
{
	out << r.x << '\t' << r.y;
	return out;
}

// arithmatic operators
XY operator+ (const XY &r1,const XY &r2)
{ return XY(r1.x+r2.x,r1.y+r2.y); }

XY operator+ (const XY &r1,double add)
{ return XY(r1.x+add,r1.y+add); }

XY operator+ (double add, const XY &r1)
{ return XY(r1.x+add,r1.y+add); }

XY operator- (const XY &r1,const XY &r2)
{ return XY(r1.x-r2.x,r1.y-r2.y); }

XY operator- (const XY &r1,double sub)
{ return XY(r1.x-sub,r1.y-sub); }

XY operator* (const XY &r1,const XY &r2)
{ return XY(r1.x*r2.x,r1.y*r2.y); }

XY operator* (const XY &r1, double mult)
{ return XY(r1.x*mult,r1.y*mult); }

XY operator* ( double mult,const XY &r1)
{ return XY(r1.x*mult,r1.y*mult); }

XY operator/ (const XY &r1,const XY &r2)
{ return XY(r1.x/r2.x,r1.y/r2.y); }

XY operator/ (const XY &r1, double div)
{ return XY(r1.x/div,r1.y/div); }



// implement
namespace xy {

	double dist(const XY &r1, const XY &r2);
	double dist_sq(const XY &r1, const XY &r2);
	double dist_pbc(const XY &r1, const XY &r2,double L);
	double dist_sq_pbc(const XY &r1, const XY &r2, double L)
		{
			XY d = r1 - r2;
			d.x -= L*round(d.x/L);
			d.y -= L*round(d.y/L);
			return d.x*d.x+d.y*d.y;
		}

	// dot product cross product
	double dot(const XY &r1, const XY &r2) {
		return r1.x*r2.x + r1.y*r2.y ;
	}

	double cross(const XY &r1, const XY &r2) {
		return r1.x*r2.y - r1.y*r2.x;
	}
};




#endif
