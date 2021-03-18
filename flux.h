#ifndef GUARD_FLUX_H
#define GUARD_FLUX_H

class Flux {
public:
    Flux(double bs, double max);
    
    void sample(const System &system);


    //write to out stream
    void writeX(std::ostream &out);
    void writeX(const char* outname);
    void writeY(std::ostream &out);
    void writeY(const char* outname);
    void write_bins(std::ostream &out);
    void write_bins( const char* outname);

    unsigned int get_Nsample() {return Nsample;}

private:
    double bs;
    unsigned int Nbin;

    std::vector<std::vector<XY> > J;
    std::vector<double> bins;

    unsigned int Nsample;

};


Flux::Flux(double bs, double max)
: bs(bs), Nbin( std::ceil(max/bs) ),
    J(  Nbin, std::vector<XY>(Nbin, XY(0,0) ) ),
    bins(  Nbin, .0 ), Nsample(0) 
{
    for(unsigned int i=0; i<Nbin; ++i) {
        bins[i] = (i+0.5)*bs;
    }
}

void Flux::sample(const System &system)
{

    XY r;
    ++Nsample;
    unsigned int jx,jy;

    for(unsigned int i=0; i<system.N; ++i) {

        r = system.r[i] - system.dr[i];
        r.pbc(system.L);

        jx = std::floor(r.x/bs);
        jy = std::floor(r.y/bs);

        if( ( jx<Nbin ) && ( jy<Nbin ) ){
            J[jx][jy] += system.v[i];
        }
    }
}

void Flux::writeX( std::ostream &out )
{

    for(unsigned int iy=0; iy < Nbin; ++iy ) {
        for(unsigned int ix=0; ix < Nbin; ++ix ) {
            out << J[ix][iy].x/Nsample;
            if( ix < Nbin - 1 ) out << '\t';
        }
        if(iy < Nbin-1) out << "\n";
    }
}

void Flux::writeY( std::ostream &out )
{

    for(unsigned int iy=0; iy < Nbin; ++iy ) {
        for(unsigned int ix=0; ix < Nbin; ++ix ) {
            out << J[ix][iy].y/Nsample;
            if( ix < Nbin - 1 ) out << '\t';
        }
        if(iy < Nbin-1) out << "\n";
    }
}

void Flux::writeX(const char* outname)
{
    std::ofstream out;
    out.open(outname);
    writeX(out);
    out.close();
}

void Flux::writeY(const char* outname)
{
    std::ofstream out;
    out.open(outname);
    writeY(out);
    out.close();
}


void Flux::write_bins( std::ostream &out )
{
    for(unsigned int i=0; i<Nbin; ++i) {
        out << bins[i];
        if( i < Nbin-1 ) out << '\t';
    }

}
void Flux::write_bins( const char* outname)
{
    std::ofstream out;
    out.open(outname);
    write_bins(out);
    out.close();
}


#endif
