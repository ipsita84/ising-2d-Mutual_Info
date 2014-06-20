//  g++ -Wall -O3 `pkg-config --cflags --libs gsl tabdatrw-0.3 interp2dpp` ising-2d-Renyi-by-interpol.cc -o testo
// 2nd Renyi entropy for classical 2d Ising model in zero magnetic field 
//Metropolis algorithm employed
//Parameters that can be changed for different runs:
//J, axis1, axis2, N_mc

#include <iostream>
#include <fstream>
#include <ctime>
#include <boost/multi_array.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <gsl/gsl_integration.h>
#include <interp2d.hpp> // For interpolation
#include <tabdatrw.hpp> // For tabdatr and tabdatw
#include <vector>

// gen is a variable name
// Its data-type is boost::random::mt19937
boost::random::mt19937 gen;
using namespace std;


double J=1.0;

//Function templates
double f (double , void * params);
double g (double , void * params);

int main()
{	

	double kT(0), T_min(0), T_max(0), del_T(0);

	cout << "Enter minimum T (in units of k)" << endl;
	cin >> T_min;

	cout << "Enter maximum T (in units of k)" << endl;
	cin >> T_max;

	cout << "Enter increment of T (in units of k) at each step" << endl;
	cin >> del_T;

	
	double mut_info(0); //mutual information I_2

	ofstream fout("mutual-info.dat"); // Opens a file for output
	
	vvdouble vm = tabdatr("Em-32by32.dat", 2);//modified energy data
	interp_data idm(vm,1);
	
        vvdouble vn = tabdatr("E-32by32.dat", 2);//normal energy data
	interp_data idn(vn,1);
	
	gsl_integration_workspace * w 
          = gsl_integration_workspace_alloc (1000);

	for (kT = T_min; kT < T_max + del_T; kT += del_T)
	{	double beta = 1.0/kT ;
	 	mut_info = 0 ;
	 	
	 	double term1(0), term2(0), term3(0), abs_error(0);
	 	gsl_function F;
  		F.function = &f;
		F.params = &idn;
//Function: int gsl_integration_qags (const gsl_function * f,double a,double b, double epsabs, double epsrel,size_t limit,gsl_integration_workspace * workspace,double * result, double *abserr)
//gsl_integration_cquad (const gsl_function * f, double a, double b, double epsabs, double epsrel, gsl_integration_cquad_workspace * workspace, double * result, double * abserr, size_t * nevals)		
  		gsl_integration_qags (&F, 0, beta, 1e-6, 1e-4,1000, w, &term2, &abs_error);
  		gsl_integration_qags (&F, 0, 2.0*beta, 1e-6, 1e-4,1000, w, &term3, &abs_error);
  		
  		F.function = &g;
  		F.params = &idm;
  		gsl_integration_qags (&F, 0, beta, 1e-6, 1e-4,1000, w, &term1, &abs_error);
		
	        
	        mut_info =2.0*term1 -2.0* term2 - term3;
	  
		fout << kT / J << '\t' << mut_info /32 << endl;
	
	}

	

	fout.close();
	return 0;

}

//performing numerical integration using gsl
double f (double beta, void * params) 
{
  	interp_data p = *(interp_data *) params;
	return p.interp_akima(beta);
}

double g (double beta, void * params) 
{
  
  	interp_data p = *(interp_data *) params;
	return p.interp_akima(beta);
}



