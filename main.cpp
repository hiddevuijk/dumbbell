
#include "ConfigFile.h"


#include "xy.h"
#include "vfield.h"
#include "system.h"
#include "density.h"
#include "orientation.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <vector>
#include <string>


using namespace std;



int main()
{


	// read input into config
	ConfigFile config("input.txt");

	// read integration parameters
	Integration int_params(config);

	// read system parameters
	System system(config);
	// start with random config. 
	system.init_random();


	system.write("initial_config.dat");

	// objects to sample density, orientation
	Density_xy density(int_params.bs,system.L,system.N);
	Orientation orientation(int_params.bs,system.L);

	// integrate Nt_init time steps
	unsigned int ti;
	cout << "Starting with equilibration ...\n";
	for( ti = 0; ti < int_params.Nt_init; ++ ti) {
		// print progress
		if( (ti%int_params.print_freq) == 0 ) {
			cout << (int_params.Nt_init + int_params.Nt) << '\t';
			cout << ti << endl;
		}

		// make t_unit time steps
		for(unsigned int tti=0; tti < int_params.t_unit; ++tti)	{
			system.step();
		}
	}

	cout << "Ended equilibration. Starting sampling ... \n";

	for(; ti < (int_params.Nt+int_params.Nt_init); ++ti) {
		// print progress 
		if( (ti%int_params.print_freq) == 0 ) {
			cout << (int_params.Nt_init + int_params.Nt) << '\t';
			cout << ti << endl;
		}

		// make t_unit time steps
		for(unsigned int tti = 0;tti<int_params.t_unit;++tti) 
			system.step();
	
		if( (ti%int_params.sample_freq) == 0 ) {
			density.sample(system);
            orientation.sample(system);
		}

	}


	cout << "Simulation finished.\nNormalizing and writing results ..." << endl;
	// normalize and save density
	density.normalize();
	density.write("rho.dat");
	density.write_bins("rho_bins.dat");

	// normalize and save orientation
	orientation.normalize();
	orientation.write("theta.dat");
	orientation.write_bins("theta_bins.dat");


	// write final configuration
	system.write("final_config.dat");
	return 0;
}



