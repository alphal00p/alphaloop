/*
  demo-c.c
  test program for the Cuba library
  last modified 13 Mar 15 th
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iomanip> //setw
#include <map>
#include <string>
#include "DCD.h"
#include "integrand/PS_points.h"

#include <boost/program_options.hpp>                                                      
                                                                                         
//using namespace boost;                                                                    
namespace po = boost::program_options;  
#if REALSIZE == 16
#include "cubaq.h"
#elif REALSIZE == 10
#include "cubal.h"
#else
#include "cuba.h"
#endif

DIdeform::R4vector p1, p2, p3, p4;
std::vector<DIdeform::R4vector> Qs(4);

extern int Integrand(const int *ndim, const cubareal xx[],
		     const int *ncomp, cubareal ff[], void *userdata) ;


/*********************************************************************/

#define NDIM 4
#define NCOMP 2
#define USERDATA NULL
#define NVEC 1
#define EPSREL 5.e-3
#define EPSABS 1.0e-15
#define VERBOSE 2 //2
#define LAST 4
#define SEED 10
#define MINEVAL 0
#define MAXEVAL 5000000

#define NSTART 1000
#define NINCREASE 500
#define NBATCH 1000
#define GRIDNO 0
#define STATEFILE NULL
#define SPIN NULL

#define NNEW 1000
#define NMIN 2
#define FLATNESS 25.

#define KEY1 47
#define KEY2 1
#define KEY3 1
#define MAXPASS 5
#define BORDER 0.
#define MAXCHISQ 10.
#define MINDEVIATION .25
#define NGIVEN 0
#define LDXGIVEN NDIM
#define NEXTRA 0

#define KEY 0

void print_help(const char * op, const char * description){
    std::cout << "\t" << std::setw(20) << std::left << op;
	std::cout << std::setw(10) << "";
	std::cout << description << std::endl;
}
int which_integrand;

int main(int argc, char **argv)
{
	//Precision
	std::cout.setf(std::ios::scientific, std::ios::floatfield);
	std::cout.precision(16);

	int comp, nregions, neval, fail;
	cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];

	//Cuba options
	int maxeval, verbose, nstart, nincrease, nbatch;
	double epsrel;
	std::string integrator_name;

	//Integrand options
	int ps_seed;
	double ps_angle;
	std::string integrand_name;

	std::map<std::string, int> integrand_function;
	integrand_function["box1L_6d"] = 0;
	integrand_function["box1L_offshell"] = 1;
	integrand_function["box1L_subtracted"] = 2;

	po::options_description desc("Allowed options");
	desc.add_options()
	("help,h", "produce help message.")
	/*Integrand options*/
	("seed", po::value<int>(&ps_seed)->default_value(0), "choose seed for phase space point.")
	("angle", po::value<double>(&ps_angle)->implicit_value(M_PI/2.0), "set scattering angle for seed=0")
	("pi_angle", po::value<double>(&ps_angle)->implicit_value(1.0)->default_value(0.5),"set scattering angle as multiple of pi")
	("integrand,i", po::value<std::string>(&integrand_name)->implicit_value("")->default_value("box1L_6d"), "choose integrand, call it with no argument to see possible arguments")
	/*Cuba options*/
	("integrator,I",po::value<std::string>(&integrator_name)->implicit_value("")->default_value("Vegas"),"choose cuba integrator among Vega, Cuhre, Divonne, Suave.")
	("maxeval", po::value<int>(&maxeval)->default_value(MAXEVAL), "set cuba MAXEVAL")
	("verbose", po::value<int>(&verbose)->default_value(VERBOSE), "set cuba VERBOSE")
	("nstart", po::value<int>(&nstart)->default_value(NSTART), "set cuba NSTART")
	("nincrease", po::value<int>(&nincrease)->default_value(NINCREASE), "set cuba NINCREASE")
	("nbatch", po::value<int>(&nbatch)->default_value(NBATCH), "set cuba NBATCH")
	("epsrel", po::value<double>(&epsrel)->default_value(EPSREL), "set cuba EPSREL")
        ;
	//Read out all the options
	po::positional_options_description p;
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("help"))
	{
		std::cout << "Usage: options_description [options]\n";
		std::cout << desc;
		return 0;
	}

	if (vm.count("pi_angle"))
	{
		ps_angle = M_PI * ps_angle;
	}
	if (integrand_function.find(integrand_name) == integrand_function.end())
	{
		printf("\033[1;31mNo integrand found with name %s\033[0m\nAllowed options are:\n", integrand_name.c_str());
		for (auto it = integrand_function.begin(); it != integrand_function.end(); it++)
			std::cout << " * " << it->first << std::endl;
		std::cout << std::endl;
		exit(1);
	}
	else
	{
		printf("\nSelect the integral %s with id %d\n", integrand_name.c_str(), integrand_function[integrand_name]);
		which_integrand = integrand_function[integrand_name];
	}

	printf("Seed = %d", ps_seed);
	if (ps_seed == 0)
		printf(", with scattering angle: %.10f", ps_angle);
	printf("\n");

	//Get Phase Space Points
	PS_points(p1, p2, p3, p4, ps_seed, ps_angle);

	std::cout << "p1:\t" << p1 << std::endl;
	std::cout << "p2:\t" << p2 << std::endl;
	std::cout << "p3:\t" << p3 << std::endl;
	std::cout << "p4:\t" << p4 << std::endl;

	Qs[0] = p1;
	Qs[1] = Qs[0] + p2;
	Qs[2] = Qs[1] + p3;
	Qs[3] = 0.0 * Qs[3];

	DIdeform::R4vector shift(Qs[0]);
	for (int i = 0; i < 4; i++)
		Qs[i] = Qs[i] - shift;

	//Integrator name is converted to all upper cases
	std::transform(integrator_name.begin(), integrator_name.end(), integrator_name.begin(), ::toupper);
	printf("----------------- Integrate using %s -----------------\n", integrator_name.c_str());

	if (integrator_name.compare("VEGAS") == 0)
	{
		Vegas(NDIM, NCOMP, Integrand, USERDATA, NVEC,
			  epsrel, EPSABS, verbose, SEED,
			  MINEVAL, maxeval, nstart, nincrease, nbatch,
			  GRIDNO, STATEFILE, SPIN,
			  &neval, &fail, integral, error, prob);
}
else if (integrator_name.compare("SUAVE") == 0)
{
	Suave(NDIM, NCOMP, Integrand, USERDATA, NVEC,
		  epsrel, EPSABS, verbose | LAST, SEED,
		  MINEVAL, maxeval, NNEW, NMIN, FLATNESS,
		  STATEFILE, SPIN,
		  &nregions, &neval, &fail, integral, error, prob);
}
else if (integrator_name.compare("DIVONNE") == 0)
{
	Divonne(NDIM, NCOMP, Integrand, USERDATA, NVEC,
			epsrel, EPSABS, verbose, SEED,
			MINEVAL, maxeval, KEY1, KEY2, KEY3, MAXPASS,
			BORDER, MAXCHISQ, MINDEVIATION,
			NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
			STATEFILE, SPIN,
			&nregions, &neval, &fail, integral, error, prob);
}
else if (integrator_name.compare("CUHRE") == 0)
{
	Cuhre(NDIM, NCOMP, Integrand, USERDATA, NVEC,
		  epsrel, EPSABS, verbose | LAST,
		  MINEVAL, maxeval, KEY,
		  STATEFILE, SPIN,
		  &nregions, &neval, &fail, integral, error, prob);
}

//Print results
printf("%s RESULT:\tneval %d\tfail %d\n",
	   integrator_name.c_str(), neval, fail);
for (comp = 0; comp < NCOMP; ++comp)
	printf("%s RESULT:\t%.8g +- %.8g\tp = %.3g\n",
		   integrator_name.c_str(),
		   (double)integral[comp], (double)error[comp], (double)prob[comp]);

printf("\nFrendly Output:\n");
printf("\t    ( %+.8e %+.8e I )\n",
	   (double)integral[0], (double)integral[1]);
printf("\t+/- (  %.8e +%.8e I )\n",
	   (double)error[0], (double)error[1]);

return 0;
}

