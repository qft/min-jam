#include <iostream>
#include "State/StaticState.h"
#include "Potentials/Potentials.h"
//#include "Potentials/RegisterPotentials.h"
#include "Boundaries/Boxes.h"
#include "Resources/Exception.h"
#include "Computers/StaticComputer.h"
#include "Computers/MatrixInterface.h"
#include "Minimization/minimizer.h"
#include "Database/Database.h"

//#include "Resources/Clusterizer.h"



using namespace LiuJamming;
#define DIM 2
using namespace std;

/*
typedef Eigen::Matrix<dbl,DIM,DIM> dmat;

void get_filename(char *filename, int N, dbl phi)
{
	sprintf(filename, "../fullrange_states/states_N%05i_f%8.6f.nc", N, phi);
}

void get_data_filename(char *filename, int N, dbl phi)
{
	sprintf(filename, "../fullrange_states/data_N%05i_f%8.6f.nc", N, phi);
}
*/

int CheckMinSuccess(CStdData<DIM> const &data, dbl p0, dbl tol, dbl rel_ptol)
{
	dbl rel_dp = (data.Pressure-p0)/p0;
	printf("checking minimization success: MaxGrad/tol = %e, rel_dp = %e, success = %i\n", data.MaxGrad/tol, rel_dp, (data.MaxGrad>tol||fabs(rel_dp)>rel_ptol)?0:1);
	if(data.MaxGrad > tol) return 0;
	if(fabs(rel_dp) > rel_ptol) return 0;
	return 1;
}

int CheckIsJammed(CStdData<DIM> const &data)
{
	if(data.NPp == 0)
		return 0;
	if(data.NcmNciso() <= 0)
		return 0;
	return 1;
}

void GetStateTargetP(CStaticState<DIM> &s, CStdData<DIM> &data, vector<std::string> &stepSummary, dbl p0, dbl tol, dbl rel_ptol)
{
	stepSummary.clear();

	//Assume that s starts off properly minimized with data calculated
	CStaticState<DIM> stemp = s;

	dbl phi1, phi2;
	dbl delta_p, B, V, delta_phi;
	dbl factor = 1.0;

	char temp[512];
	bool is_jammed = CheckIsJammed(data);
	dbl rel_dp = (data.Pressure - p0)/p0;
	sprintf(temp, "%3i %6i     %e       %f    %6i     %6i    %6i     % f     %e     %e    % e         %i", 0, data.NPp,data.MaxGrad,2.*data.Nc/((dbl)s.GetParticleNumber()),data.Nc,DIM*data.NPp-DIM,data.NcmNciso(),s.GetPackingFraction(),data.Energy,data.Pressure,rel_dp,is_jammed);
	stepSummary.push_back( std::string(temp) );

	int ii = 1;
	int ii_max = 6;
	while(!CheckMinSuccess(data, p0, tol, rel_ptol) && ii <= ii_max)
	{
		//First, set the new packing fraction 
		//Assume harmonic approximation
		//delta phi = (V^2/(B*phi)) * delta p
		delta_p = p0 - data.Pressure;
		phi1 = s.GetPackingFraction();
		B = data.cijkl.CalculateBulkModulus();
		phi2 = phi1 + factor*(phi1/B)*delta_p;

		printf("setting packing fraction to %f\n", phi2);
		s.SetPackingFraction(phi2);

		//Minimize the energy
		CStaticComputer<DIM> c(s);
		CSimpleMinimizer<DIM> miner(c);
		miner.minimizeFIRE(tol);
	
//		CBondList<DIM> tempBonds;
//		CStdData<DIM> tempData;
//		c.ComputeBondList(tempBonds);
//		if(!tempBonds.Empty())
//			tempBonds.CalculateStdData(tempData);
//		printf("tempBonds.MaxGrad = %e\n", tempData.MaxGrad);

		//Prepare system for standard calculations
		//CStaticComputer<DIM>::StdPrepareSystem() sets the internal bond list and removes rattlers
		//it returns 0 if everything is done correctly.
		if(c.StdPrepareSystem())
		{
			printf("The bonds list is empty...\n");
		}else
			c.CalculateStdData(); //Have option to not write down hessian.
		
		c.Data.Print(); //This should just be a bunch of zeros if bond list is empty.
		
		is_jammed = CheckIsJammed(c.Data);
		rel_dp = (c.Data.Pressure - p0)/p0;
		sprintf(temp, "%3i %6i     %e       %f    %6i     %6i    %6i     % f     %e     %e    % e         %i", ii, c.Data.NPp,c.Data.MaxGrad,2.*c.Data.Nc/((dbl)s.GetParticleNumber()),c.Data.Nc,DIM*c.Data.NPp-DIM,c.Data.NcmNciso(),s.GetPackingFraction(),c.Data.Energy,c.Data.Pressure,rel_dp,is_jammed);
		stepSummary.push_back( std::string(temp) );

		if(is_jammed)
		{
			stemp = s;
			data = c.Data;
			factor = 1.0;
		}
		else
		{
			s = stemp;
			
			factor *= 0.5;
//			if(factor < 0.2)
//				factor = -1.0;
//			printf("WARNING! took too big of a step. This is not handled properly yet. Exiting.\n");
//			assert(false);
		}

		++ii;
	}

	printf("  i      N         max_grad     z_norattle        Nc      Nciso  NcmNciso           phi           energy         pressure           rel_dp    jammed\n");
	for(int i=0; i<(int)stepSummary.size(); ++i)
		printf("%s\n", stepSummary[i].c_str());
	printf("----------------------------------------------------------------------------------------------------------------------------------------------------\n");
	fflush(stdout);
}

void GetStatesRangeP(int N, int seed, int dispersity, dbl Lpmax, dbl Lpmin, dbl Lpinc, dbl tol, dbl rel_ptol)
{
	//set up the output
	char dirbase[256], fn_disp_part[256], fn_pot_part[256];
	sprintf(fn_pot_part, "harmonic");
	switch(dispersity)
	{
		case 0:
			sprintf(fn_disp_part, "mono");
			break;
		case 1:	
			sprintf(fn_disp_part, "bi");
			break;
		case 2:
			sprintf(fn_disp_part, "poly_uniform");
			break;
		default:
			assert(false);
	}
	sprintf(dirbase, "/data/p_ajliu/cpgoodri/get_states2/%id_%s/%s", DIM, fn_pot_part, fn_disp_part);

	char statedb_fn[256], datadb_fn[256];
	sprintf(statedb_fn, "%s/temp/state_N%08i_r%06i.nc", dirbase, N, seed);
	sprintf(datadb_fn,  "%s/temp/data_N%08i_r%06i.nc",  dirbase, N, seed);

	CStaticDatabase<DIM>  state_db(N,statedb_fn,NcFile::Replace);
	CStdDataDatabase<DIM> data_db(datadb_fn,NcFile::Replace);

	//Create the system
	CStaticState<DIM> s(N);
	CStdData<DIM> data;
	s.RandomizePositions(seed);

	switch(dispersity)
	{
		case 0:
			s.SetRadiiMono();
			break;
		case 1:	
			s.SetRadiiBi();
			break;
		case 2:
			s.SetRadiiPolyUniform();
			break;
		default:
			assert(false);
	}

	dbl initial_phi = 1.0;
	s.SetPackingFraction(initial_phi);

	//Minimize the energy
	printf("Initial minimization...\n");
	CStaticComputer<DIM> c(s);
	CSimpleMinimizer<DIM> miner(c);
	miner.minimizeFIRE(tol);

	//Prepare system for standard calculations
	//CStaticComputer<DIM>::StdPrepareSystem() sets the internal bond list and removes rattlers
	//it returns 0 if everything is done correctly.
	if(c.StdPrepareSystem())
	{
		printf("The bonds list is empty...\n");
		assert(false);
	}
	c.CalculateStdData(); //Have option to not write down hessian.
	c.Data.Print();
	data = c.Data;

	vector<std::string> finalStepSummary;
	vector<std::string> stepSummary;
	int numMinimizations = 0;

	int numLpi = ((int)((Lpmax-Lpmin)/(-Lpinc))) + 1;
	printf("numLpi = %i\n", numLpi);
	dbl Lp;
	for(int Lpi = 0; Lpi < numLpi; ++Lpi)
	{
		Lp = Lpmax + ((dbl)Lpi)*Lpinc;
		printf("Lpi = %i, Lp = %f\n", Lpi, Lp);
		GetStateTargetP(s, data, stepSummary, pow(10.,Lp), tol, rel_ptol);

		numMinimizations += (int)stepSummary.size() - 1;
		finalStepSummary.push_back( stepSummary.back() );

		state_db.Write(s);
		data_db.Write(data);
	}
	
	printf("done getting all states. Process took %i minimizations\n", numMinimizations);

	printf("  i      N         max_grad     z_norattle        Nc      Nciso  NcmNciso           phi           energy         pressure           rel_dp    jammed\n");
	for(int i=0; i<(int)finalStepSummary.size(); ++i)
		printf("%s\n", finalStepSummary[i].c_str());
	printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");

}

//Dispersity key:
//	0 -> monodisperse [100% sigma=1]
//	1 -> bidisperse [50% sigma=1, 50% sigma=1.4]
//	2 -> poly-uniform [uniform between sigma=1 and sigma=1.4]

int main(int argc, char* argv[])
{
	int N = 2048;				//n
	int seed = 1;				//r
	int dispersity = 2;			//s


	int c;
	while((c=getopt(argc, argv, "n:r:s:")) != -1)
		switch(c)
		{
			case 'n':	N = atoi(optarg); break;
			case 'r':	seed = atoi(optarg); break;
			case 's':	dispersity = atoi(optarg); break;
			case '?':	
				if(optopt == 'c') 
					std::cerr << "Option -" << optopt << "requires an argument.\n";
				else if(isprint(optopt)) 
					std::cerr << "Unknown opton '-" << optopt << "'.\n";
				else
					std::cerr << "Unknown option character '\\x" << optopt << "'.\n";
				return 1;
			default:
				abort();
		}
	for(int index = optind; index < argc; index++)
		std::cerr << "Non-option argument " << argv[index] << "\n";


	assert(dispersity == 0 || dispersity == 1 || dispersity == 2);

	dbl Lpmax = -1.;
	dbl Lpmin = -8.;
	dbl Lpinc = -0.2;
	dbl tol = 1e-13;
	dbl rel_ptol = 1e-2;
	GetStatesRangeP(N, seed, dispersity, Lpmax, Lpmin, Lpinc, tol, rel_ptol);
	




}








