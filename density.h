#ifndef GUARD_DENSITY_H
#define GUARD_DENSITY_H

// density of the system projected on the xy plane

#include <iostream>

#include "system.h"
#include "xy.h"

#include <fstream>
#include <string> 

class Density_xy {
public:
	Density_xy( double bs, double max, int Nparticles);	

	void sample(const System &system);

    void reset();

	// write to out stream
	void write(std::ostream &out);
	void writeX(std::ostream &out);
	void write_bins(std::ostream &out);
		
	// write to file named outname
	void write(const char* outname);
	void writeX(const char* outname);
	void write_bins(const char* outname);
	//void write(const char* outname);

	unsigned int get_Nsample() {return Nsample;}
private:

	double bs;	// bin size
	unsigned int Nbin;	// number of bins

	std::vector<std::vector<double> > rho;
	std::vector<double> bins;

	unsigned int Nsample;
	unsigned int N; // number of particles
};

Density_xy::Density_xy(double bss,double max,int Nparticles)
{
	bs = bss;
	Nbin = (unsigned int) std::ceil(max/bs);

	rho = std::vector<std::vector<double> >(Nbin,std::vector<double>(Nbin,0.));

	bins = std::vector<double>(Nbin,0.);
	for(unsigned int i=0;i<Nbin;++i)
		bins[i] = (i+0.5)*bs;
	Nsample = 0;
	N = Nparticles;
}	



void Density_xy::sample(const System &system)
{
	++Nsample;
	//double x,y;
	XY r;
	unsigned int jx,jy;
	for(unsigned int i=0;i<system.N;++i ) {
		r = system.r[i];
		r.pbc(system.L);
		jx = std::floor(r.x/bs);
		jy = std::floor(r.y/bs);

		//x = system.r[i].x;
		//x -= system.L*std::floor(x/system.L);
		//jx = std::floor(x/bs);

		//y = system.r[i].y;
		//y -= system.L*std::floor(y/system.L);	
		//jy = std::floor(y/bs);

		if((jx<Nbin) && (jy<Nbin) ) {
			rho[jx][jy] += 1.;
		}
	}
		
}




void Density_xy::reset()
{

    Nsample = 0;
    for(unsigned int ix=0; ix< Nbin; ++ix) {
    for(unsigned int iy=0; iy< Nbin; ++iy) {
        rho[ix][iy] = 0;
    }}

}

void Density_xy::write(std::ostream &out)
{

	double norm = 1./(bs*bs*N*Nsample);
	for(unsigned int jy=0;jy<Nbin;++jy) {
		for(unsigned int jx=0;jx<Nbin;++jx) {
			out << rho[jx][jy]*norm;
			if(jx<(Nbin-1)) out << ' '; 
		}
		if(jy<(Nbin-1)) out << '\n';
	}
}

void Density_xy::write(const char* outname)
{
	std::ofstream out;
	out.open(outname);

    write(out);	

	out.close();

}

void Density_xy::writeX(std::ostream &out)
{

	double norm = 1./(bs*N*Nsample);
    double rhox;
    for(unsigned int jx=0;jx<Nbin;++jx) {
        rhox = 0;
        for(unsigned int jy=0;jy<Nbin;++jy) {
			rhox += rho[jx][jy]*norm;
		}
        out << rhox;
		if(jx<(Nbin-1)) out << '\t';
	}
    out << '\n';
}

void Density_xy::writeX(const char* outname)
{
	std::ofstream out;
	out.open(outname);

    writeX(out);	

	out.close();

}



void Density_xy::write_bins(std::ostream &out)
{

	for(unsigned int j=0;j<Nbin;++j) {
		out << bins[j];
		if(j<(Nbin-1)) out << ' '; 
	}

}

void Density_xy::write_bins(const char* outname)
{
	std::ofstream out;
	out.open(outname);

	for(unsigned int j=0;j<Nbin;++j) {
		out << bins[j];
		if(j<(Nbin-1)) out << ' '; 
	}

	out.close();

}
//
//void Density_xy::write(const char* outname)
//{
//	std::ofstream out;
//	out.open(outname);
//
//	for(unsigned int jy=0;jy<Nbin;++jy) {
//		for(unsigned int jx=0;jx<Nbin;++jx) {
//			out << bins[jx] << '\t'
//				<< bins[jy] << '\t'
//				<< rho[jx][jy] << '\n';
//		}
//	}
//	out.close();
//}
#endif

