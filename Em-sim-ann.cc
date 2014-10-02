//  g++ -Wall -O3 tot-en-vs-beta-Replica.cc -o testo
// 2nd Renyi entropy for classical 2d Ising model in zero magnetic field 
//Metropolis algorithm employed
//warming up system for N_mc updates
//taking avg in the next N_mc updates
//Parameters that can be changed for different runs:
//J, axis1, axis2, N_mc

#include <iostream>
#include <fstream>
#include <ctime>
#include <boost/multi_array.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/lexical_cast.hpp>


// gen is a variable name
// Its data-type is boost::random::mt19937
boost::random::mt19937 gen;
using boost::lexical_cast;
using boost::bad_lexical_cast;
using namespace std;

typedef boost::multi_array < int, 2 > array_2d;
// typedef keyword allows you to create an alias fo a data type

// Define global scopes to use them across all functions
double J = 1.0;
const unsigned int axis1 = 8, axis2 = 8;
//axis1 should be of even no. of sites
// above assigns length along each dimension of the 2d configuration
//No.of Monte Carlo updates we want
unsigned int N_mc = 1e6;

//Function templates
int roll_coin(int a, int b);
double random_real(int a, int b);
double energy_tot(array_2d sitespin);
double nn_energy(array_2d sitespin, unsigned int row, unsigned int col);

int main(int argc, char * argv[])
{
	if (argc != 4)
	{
		cout << "Expecting three inputs: beta_min, beta_max, del_beta."
		     << endl << "Got " << argc - 1 << endl;
		return 1;
	}

	double beta_min(0), beta_max(0), del_beta(0);

	try
	{
		beta_min = lexical_cast<double>(argv[1]);
		beta_max = lexical_cast<double>(argv[2]);
		del_beta = lexical_cast<double>(argv[3]);
	}
	catch (const bad_lexical_cast & x)
	{
		cout << "Cannot convert input to double" << endl;
		return 2;
	}

//	cout << "Enter minimum beta" << endl;
//	cin >> beta_min;
//	cout << "Enter maximum beta" << endl;
//	cin >> beta_max;
//	cout << "Enter increment of beta" << endl;
//	cin >> del_beta;
	ofstream fout("Em-old.dat"); // Opens a file for output
	
	fout << 0 << '\t' << 0 << endl;
	
	//define replica 1 spin configuration array
	array_2d sitespin1(boost::extents[axis1][axis2]);
	//define replica 2 spin configuration array
	array_2d sitespin2(boost::extents[axis1][axis2]);

	//Both replicas have same initial spin configuration
	for (unsigned int i = 0; i < axis1; ++i)
		for (unsigned int j = 0; j < axis2; ++j)
		{
			sitespin1[i][j] = 2 * roll_coin(0, 1) - 1;
			sitespin2[i][j] = sitespin1[i][j];
		}


	double energy = 2*energy_tot(sitespin1);
	

	for (double beta = beta_min + del_beta; beta < beta_max + del_beta; beta += del_beta)
	 {	
	  unsigned int sys_size = axis1 * axis2;
	  unsigned int row, col, label;
	  double r(0), acc_ratio(0) ;
	  double en_sum(0);

	  for (unsigned int i = 1; i <=1e5 + N_mc; ++i)
	  {
		for (unsigned int j = 1; j <= 3*sys_size/2; ++j)
		{	//Choose a random spin site for the entire 2 replica system
			double energy_diff(0);
			label = roll_coin(1,2*sys_size);
			
			//if the random spin site is located in layer 1
			if (label <= sys_size)
			{	if (label % axis2 == 0)
				{	row = (label / axis2) - 1;
					col = axis2 -1 ; 
				}
				else
				{	col = label % axis2 - 1;
					row = (label-col-1)/axis2;
				} 
			
				energy_diff=-2.0*nn_energy(sitespin1, row, col);
				if (row < axis1/2)
					energy_diff +=-2.0*nn_energy(sitespin2, row, col);

				//Generate a random no. r such that 0 < r < 1

				r = random_real(0, 1);
				acc_ratio = exp(-1.0 * energy_diff *beta);

				//Spin flipped if r <= acceptance ratio
				if (r <= acc_ratio)
				{
					sitespin1[row][col] *= -1;
					if (row < axis1/2) 
						sitespin2[row][col] *=-1;
							//this line on addition creates trouble
					energy += energy_diff;
				}
			}
			
			//if the random spin site is located in layer 2
			if (label > sys_size)
			{	label -= sys_size;
				if (label % axis2 == 0)
				{	row = (label / axis2) - 1;
					col = axis2 -1 ; 
				}
				else
				{	col = label % axis2 - 1;
					row = (label-col-1)/axis2;
				} 
			
			
				energy_diff=-2.0*nn_energy(sitespin2, row, col);
				if (row < axis1/2)
					energy_diff +=-2.0*nn_energy(sitespin1, row, col);

					 
				//Generate a random no. r such that 0 < r < 1

				r = random_real(0, 1);
				acc_ratio = exp(-1.0 * energy_diff *beta);

				//Spin flipped if r <= acceptance ratio
				if (r <= acc_ratio)
				{
					sitespin2[row][col] *= -1;
					if (row < axis1/2) 
						sitespin1[row][col] *=-1;
						
					energy += energy_diff;
				}
			}

		}

		if (i> 1e5) en_sum += energy;
		}
	
	
	fout << beta << '\t' << en_sum / N_mc << endl;	
	}
	
	fout.close();
	return 0;

}






//function to generate random integer
// between 2 integers a & b, including a & b
int roll_coin(int a, int b)
{
	boost::random::uniform_int_distribution <> dist(a, b);

	return dist(gen);

}

//function to generate random real no.
// between 2 integers a & b, including a & excluding b

double random_real(int a, int b)
{
	boost::random::uniform_real_distribution <> dist(a, b);
	// uniform_real_distribution: continuous uniform distribution 
	//on some range [min, max) of real number
	return dist(gen);

}

//function to calculate total energy
//for a given spin configuration
//with periodic boundary conditions

double energy_tot(array_2d sitespin)
{
	double energy(0);

	for (unsigned int i = 0; i < axis1 - 1; ++i)
	{
		for (unsigned int j = 0; j < axis2 - 1; ++j)
		{
			energy -= J * sitespin[i][j] * sitespin[i + 1][j];

			energy -= J * sitespin[i][j] * sitespin[i][j + 1];

		}
	}

	//periodic boundary conditions
	for (unsigned int j = 0; j < axis2; ++j)
		energy -= J * sitespin[axis1-1][j] * sitespin[0][j];

	for (unsigned int i = 0; i < axis1; ++i)
		energy -= J * sitespin[i][axis2-1] * sitespin[i][0];

	return energy;
}


//Calculating interaction energy change for spin 
//at random site->(row,col) with its nearest neighbours
double nn_energy(array_2d sitespin, unsigned int row, unsigned int col)
{
	double nn_en = 0;

	if ( row > 0 && row < axis1 - 1)
	{
		nn_en -= J * sitespin[row][col] * sitespin[row-1][col];
		nn_en -= J * sitespin[row][col] * sitespin[row+1][col];
	}

	if (col > 0 && col < axis2-1)
	{
		nn_en -= J * sitespin[row][col] * sitespin[row][col-1];
		nn_en -= J * sitespin[row][col] * sitespin[row][col+1];
	}

	if (row == 0)
	{
		nn_en -= J * sitespin[0][col] * sitespin[axis1-1][col];
		nn_en -= J * sitespin[0][col] * sitespin[1][col];

	}

	if (row == axis1-1)
	{
		nn_en -= J * sitespin[axis1-1][col] * sitespin[axis1-2][col];
		nn_en -= J * sitespin[axis1-1][col] * sitespin[0][col];

	}

	if (col == 0)
	{
		nn_en -= J * sitespin[row][0] * sitespin[row][axis2-1];
		nn_en -= J * sitespin[row][0] * sitespin[row][1];

	}

	if (col == axis2-1)
	{
		nn_en -= J * sitespin[row][axis2-1] * sitespin[row][axis2-2];
		nn_en -= J * sitespin[row][axis2-1] * sitespin[row][0];

	}
	return nn_en;
}

