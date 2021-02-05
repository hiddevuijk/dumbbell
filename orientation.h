#ifndef GUARD_ORIENTATION_H
#define GUARD_ORIENTATION_H


class Orientation {
public:
	Orientation(double bs, double max);

	void sample(const System &system);

	void normalize();
 
	// write to out stream
	void write(std::ostream &out);
	void write_bins(std::ostream &out);

	// write to file named outname
	void write(const char* outname);
	void write_bins(const char* outname);

	unsigned int get_Nsample() {return Nsample;}

private:
	double bs;
	unsigned int Nbin;

	std::vector<std::vector<double> > theta;
	std::vector<std::vector<unsigned int> > n_in_bin;
	std::vector<double> bins;

	unsigned int Nsample;
};




Orientation::Orientation(double bss,double max)
{
	bs = bss;
	Nbin = (unsigned int) std::ceil(max/bs);

	theta    = std::vector<std::vector<double> >(Nbin, std::vector<double>(Nbin,0) );
	n_in_bin = std::vector<std::vector<unsigned int> >(Nbin, std::vector<unsigned int>(Nbin,0) );

	bins = std::vector<double>(Nbin,0.);
	for(unsigned int i=0;i<Nbin;++i)
		bins[i] = (i+0.5)*bs;
	Nsample = 0;
}	



void Orientation::sample(const System &system)
{
	XY r;
	++Nsample;
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

		if((jx<Nbin) && (jy<Nbin) )  {
			theta[jx][jy] += system.theta[i];
			n_in_bin[jx][jy] += 1;
		}
	}
		
}

void Orientation::normalize() 
{
	for(unsigned int jx = 0;jx < Nbin; ++jx ) {
		for(unsigned int jy = 0;jy < Nbin; ++jy ) {
			if(n_in_bin[jx][jy] > 0 ) 
				theta[jx][jy] /= n_in_bin[jx][jy];
		}
	}
}

void Orientation::write(std::ostream &out)
{
	for(unsigned int jy=0;jy<Nbin;++jy) {
		for(unsigned int jx=0;jx<Nbin;++jx) {
			out << theta[jx][jy];
			if(jx<(Nbin-1)) out << ' '; 
		}
		if(jy<(Nbin-1)) out << '\n';
	}
}

void Orientation::write(const char* outname)
{
	std::ofstream out;
	out.open(outname);

	for(unsigned int jy=0;jy<Nbin;++jy) {
		for(unsigned int jx=0;jx<Nbin;++jx) {
			out << theta[jx][jy];
			if(jx<(Nbin-1)) out << ' '; 
		}
		if(jy<(Nbin-1)) out << '\n';
	}

	out.close();

}

void Orientation::write_bins(std::ostream &out)
{

	for(unsigned int j=0;j<Nbin;++j) {
		out << bins[j];
		if(j<(Nbin-1)) out << ' '; 
	}

}

void Orientation::write_bins(const char* outname)
{
	std::ofstream out;
	out.open(outname);

	for(unsigned int j=0;j<Nbin;++j) {
		out << bins[j];
		if(j<(Nbin-1)) out << ' '; 
	}

	out.close();

}




#endif
