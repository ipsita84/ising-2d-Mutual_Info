// Considering 3d ising model using a 3d array in zero magnetic field
//Metropolis algorithm employed

#include <iostream>
#include <fstream>
#include <ctime>
#include <vector>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

// gen is a variable name
// Its data-type is boost::random::mt19937
boost::random::mt19937 gen;
using namespace std;


// Define a global scope 'J' to use it across all functions
double J = 0;



//function to generate random integer
// between 2 integers a & b, including a & b
int roll_coin(int a, int b)
{
	boost::random::uniform_int_distribution <> dist(a, b);

	return dist(gen);

}


//function to calculate total energy
//for a given spin configuration
//with periodic boundary conditions

double energy_tot(unsigned int size, int ***sitespin)
{
	double energy = 0;

	for (unsigned int i = 0; i < size-1; ++i)
	{
		for (unsigned int j = 0; j < size-1; ++j)
		{
				
			for (unsigned int k = 0; k < size-1; ++k)
			{
				energy -= J*sitespin[i][j][k]*sitespin[i+1][j][k];
			    
				energy -= J*sitespin[i][j][k]*sitespin[i][j+1][k];
				
				energy -= J*sitespin[i][j][k]*sitespin[i][j][k+1];
			 }   
		}
	}
	
	//periodic boundary conditions
	for (unsigned int i=0 ; i < size ; ++i)
	{
		for (unsigned int j=0 ; j < size ; ++j)
		{
			energy -= J * sitespin[i][j][size-1] *  sitespin[i][j][0];
			energy -= J * sitespin[i][size-1][j] *  sitespin[i][0][j];
			energy -= J * sitespin[size-1][i][j] *  sitespin[0][i][j];
		}
	}
		
	
	
	
			
	return energy ;
}



//function to calculate magnetization
//for a given spin configuration
	
double mag_tot(unsigned int size, int ***sitespin)
{
	int mag = 0;

	for (unsigned int i = 0; i < size ; ++i)
	{
		for (unsigned int j = 0; j < size ; ++j)
		{
			for (unsigned int k = 0; k < size ; ++k) mag +=sitespin[i][j][k];
	
		}
	}
	

			
	return mag ;
}



//function to generate random real no.
// between 2 integers a & b, including a & excluding b
double random_real(int a, int b)
{
	boost::random::uniform_real_distribution <> dist(a, b);
	// uniform_real_distribution: continuous uniform distribution on some range [min, max) of real nos

	return dist(gen);

}










int main()
{	unsigned int size = 0;

	cout << "Enter system size as cube length" << endl;
	cin >> size;
	arrcol = size;

	int ***sitespin;
	// stores the spin configuration of the system
	//initial state chosen by random no. generator above
	//complicated allocation below required as system size is not fixed
	//before user enters  from keyboard while running code
	
	sitespin = new int ** [size];
	for (unsigned int i = 0; i < size; ++i)
	{
		sitespin[i] = new int * [size];
		for (unsigned int j = 0; j < size; ++j)
		{
			sitespin[i][j] = new int [size];
		}
	}	
	
	
	for (unsigned int i = 0; i < size; ++i)
	{
		for (unsigned int j = 0; j < size; ++j)
		{
			for (unsigned int k = 0; k < size; ++k)
				sitespin[i][j][k] = pow(-1,roll_coin(0,1));
		}
	}
	

	double kT(0), beta(0);

	cout << "Enter T (in units of k)" << endl;
	cin >> kT;
	
	beta = 1.0/kT;
	
	
	//Calculating initial energy for this configuration above
	cout << "Enter value of coupling constant J > 0 " << endl;
	cin >> J;
	
	double energy = energy_tot(size, sitespin);
	
	//Input from terminal number of equilibration steps you wnat

	unsigned int maxstep(0);

	cout << "Enter no. of equilibration steps" << endl;
	cin >> maxstep;
	
	double en_sum(0), mag_sum(0);
	
	
	ofstream f3out("3d.dat"); // Opens a file for output

	for (unsigned int i = 1; i <= maxstep; i++)
	{

		//Now choose a random spin site, say, i

		unsigned int rnd_row , rnd_col, rnd_dep;
		rnd_row = roll_coin(1, size)-1;
		rnd_col = roll_coin(1, size)-1;
		rnd_dep = roll_coin(1, size)-1;

		//Calculating new energy for flipping spin 
		//at random site->(rnd_row,rnd_col, rnd_dep)
		
		sitespin[rnd_row] [rnd_col] [rnd_dep]*= -1;

		double new_energy = energy_tot(size , sitespin);


		//Calculating change in energy for the above spin flip

		double energy_diff = new_energy - energy;
		
		

		//Generate a random no. r such that 0 < r < 1

		double r = random_real(0, 1);
		double acc_ratio = exp(-1.0 * energy_diff * beta);

		if (r > acc_ratio) sitespin[rnd_row] [rnd_col] [rnd_dep]*= -1;
		//Spin not flipped if r > acceptance ratio
			
		
		//Given the energy of ising system at a selection of times 
		// during the simulation, we can average them to find 
		//the estimators of internal energy
		//dividing it by no. of sites gives internal energy per site
		
		en_sum += energy_tot( size , sitespin ) / (size*size*size) ;
		mag_sum += mag_tot ( size  , sitespin ) * 1.0 / (size*size*size) ;
		
		if ((i % 1000) == 0)
			f3out << i * 1.0 / (size*size*size) << '\t' << en_sum / i << '\t' << mag_sum / i << endl;

	}
	
	for (unsigned int i = 0; i < size; ++i)
	{
		for (unsigned int j = 0; j < size; ++j)
			delete [] sitespin[i][j];
		delete [] sitespin[i];
	}
	delete [] sitespin;
	
	f3out.close();
	
	
	return 0;
	
	
	}
