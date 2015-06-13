#include <iostream>
#include <ctime>
#include "State/StaticState.h"
#include "Potentials/Potentials.h"
#include "Boundaries/Boxes.h"
#include "Resources/Exception.h"
#include "Computers/StaticComputer.h"
#include "Computers/MatrixInterface.h"
#include "Minimization/minimizer.h"
#include "Database/Database.h"
using namespace std;
using namespace LiuJamming;

#define DIM 2
typedef CStaticState<DIM> STATE;
typedef CStaticComputer<DIM> COMP;
typedef CBondList<DIM> BLIST;
typedef CSimpleMinimizer<DIM> SMINER;
typedef MatrixInterface<dbl> MI;
typedef Eigen::Matrix<dbl, DIM, 1> dvec;

/*	Test 1: Quick Diagonalization
 *
 *	This shows how to find and print the first 100 eigenvalues in very few lines.
 */
void test1(STATE &s)
{
	COMP c(s);						//Create a Computer
	c.StdPrepareSystem();			//Prepare system for standard calculations (this sets the internal BondList and removes rattlers).
	c.CalculateStdData();			//This also creates a Matrix Interface and sets the Hessian. The Matrix Interface can be accessed in c.Data.H
	c.Data.Print();					//Print the data  (not necessary)
	MI &H = c.Data.H;				//For convenience, create a reference of the Matrix Interface (not necessary)
	H.VDiagonalize(100);			//Diagonalize the hessian to find the first 100 eigenvalues, and print the report to the screen
}

/*	Test 2: Reading the data
 *
 *	This shows how to read the eigenvalues and eigenvectors
 */
void test2(STATE &s)
{
	//Setup
	COMP c(s);						//Create a Computer
	c.StdPrepareSystem();			//Prepare system for standard calculations (this sets the internal BondList and removes rattlers).
	c.CalculateStdData();			//This also creates a Matrix Interface and sets the Hessian. The Matrix Interface can be accessed in c.Data.H
	c.Data.Print();					//Print the data
	MI &H = c.Data.H;				//For convenience, create a reference of the Matrix Interface.

	//Diagonalize the system
	H.Diagonalize(100);				//Diagonalize the hessian to find the first 100 eigenvalues

	//Print out the eigenvalues and the first 3 components of the corresponding eigenvectors
	//
	//	- The Matrix Interface has a parameter call num_converged. ALWAYS use this parameter when looping over eigenvalues/eigenvectors.
	//		Do not assume that the number of converged eigenvalues equals the number of requested eigenvalues.
	//	- The eigenvalues are stored in an array of length H.num_converged called H.Eigenvalues
	//	- The eigenvectors are stored in an array of length H.num_converged*H.A.rows() called H.Eigenvectors
	//			The first eigenvector is stored in the first H.A.rows() entries, etc.
	//
	//	- The Eigenvalues and Eigenvectors can be accessed through the functions
	//			GetVal(int m) const;			//m is the mode index
	//			GetVec(int m, int i) const;		//m is the mode index, i is the component of the eigenvector
	for(int i=0; i<H.num_converged; ++i)
	{
		printf("Eigenvalue[%5i] = % e     The first three components of the eigenvector are % e % e % e\n", i, H.GetVal(i),H.GetVec(i,0),H.GetVec(i,1),H.GetVec(i,2));
	}
}

/*	Test 3: Options
 *
 *	This shows some of the options that can be used
 */
void test3(STATE &s)
{
	//Setup
	COMP c(s);						//Create a Computer
	c.StdPrepareSystem();			//Prepare system for standard calculations (this sets the internal BondList and removes rattlers).
	c.CalculateStdData();			//This also creates a Matrix Interface and sets the Hessian. The Matrix Interface can be accessed in c.Data.H
	c.Data.Print();					//Print the data
	MI &H = c.Data.H;				//For convenience, create a reference of the Matrix Interface.

	//The default/simplest way to diagonalize the system is to call:
	//				H.Diagonalize();

	//There are currently 3 primary options
	//	1. Choose the number of eigenvalues/vectors to calculate.
	//		- By default, all eigenvalues are calculated (except the largest, for unknown reasons having to do with ARPACK being stupid)
	//		- This is easily specified by as in input argumen to H.Diagonalize. For example, to calculate the lowest 100 eigenvalues:
	//				H.Diagonalize(100);
	//		- Alternatively, one can call H.SetNumRequest(100); prior to calling the H.Diagonalize() routine
	//
	//	2. Choose to print the default diagonalization report after diagonalization is complete.
	//		- This is easily done by including a "V" (for verbose) in the function call, e.g.
	//				H.VDiagonalize(100);
	//		- After diagonalization, one can also call 
	//				H.Report(45);
	//			to print the report on the first 45 modes. If the argument is omitted, 30 modes are reported on by default.
	//
	//	3. Choose to diagonalize using the Shift and Invert functionality of ARPACK (see http://docs.scipy.org/doc/scipy/reference/tutorial/arpack.html for a simple explanation).
	//		- WARNING: While this is significantly faster (up to ~25X speedup), the accuracy of the eigenvalues is significantly decreased. Use at your own risk.
	//		- This is easily done by includeing "_SI" at the end of the function call, e.g.
	//				H.Diagonalize_SI(100);
	//		- Alternatively, one can call H.SetShiftAndInvert(); prior to calling the H.Diagonalize routine
	//
	//	In Summary, there are currently 8 primary routines that can be called:
	//				H. Diagonalize   ();
	//				H. Diagonalize   (100);
	//				H.VDiagonalize   ();
	//				H.VDiagonalize   (100);
	//				H. Diagonalize_SI();
					H. Diagonalize_SI(100);
	//				H.VDiagonalize_SI();
	//				H.VDiagonalize_SI(100);
	//
	//Additionally, one can choose to only compute the eigenvectors by calling
	//				H.SetComputeVecs(false);
	//		prior to the diagonalization routine.
	//
	//TODO:
	//	- Option for changing the region where eigenvalues are searched (either via "which" or via "sigma" in the shift and invert mode).
	//	- Easily calculate the ENTIRE spectrum... this might not be necessary bc probably only applicable to small systems when Eigen dense diagonalization would suffice
	//	- Diagonalize by first converting to an Eigen dense matrix.
	//
	//
	H.Report(100);

}

void saveCentersToSSV(Eigen::VectorXd pos, int N, char* filename, double phi)
{
	FILE* out;
	stringstream ss;
	static int number = 0;
	ss << "Minimized" << filename << "_phi" << phi << "monoRadii.txt";
	char *fn = (char *) malloc(500);
	ss >> fn;
	out = fopen(fn, "w");
	delete(fn);

	for(int i = 0; i < N; i++)
	{
		fprintf(out, "%17.16e %17.16e\n", pos[2 * i], pos[2 * i + 1]);
	}
	
	fclose(out);

}

int CheckIsJammed(CStdData<DIM>& data)
{
	if (data.NPp == 0)
		return 0;
	if (data.NcmNciso() <= 0)
		return 0;
	return 1;
}

void writePsi6Data(char filename[], vector< vector<double> > psi6p, double phi)
{
	char fn[500];
	sprintf(fn, "%sPsi6Pphi_%.6f", filename, phi);

	FILE* out = fopen(fn, "w");
	if (out == NULL) exit(-50);
	const int len = psi6p.size();
	for (int i = 0; i < len; i++)
	{
		fprintf(out, "%.16e %.16e", psi6p[i][0], psi6p[i][1]);
	}
	fclose(out);
}

void run(int* N, dbl phi, int seed, char filename[], int test_num=1)
{
	//read in the data
	//Create the system
	vector<double> temp;

	char* fn = (char *) malloc(200);
	char* xS = (char *) malloc(100);
	char* yS = (char *) malloc(100);

	char* fileNameOrig  = (char *) malloc(500);
	strcpy(fileNameOrig, filename);

	sprintf(filename, "%s.txt", filename);
	ifstream in(filename);
	double x, y;

	while (in.good())
	{
		in.getline(fn, 200);
		sscanf(fn, "%s %s", xS, yS);
		istringstream a(xS);
		a >> x;
		istringstream b(yS);
		b >> y;
		if(x >= 0 && y >= 0)
		{
			temp.push_back(x);
			temp.push_back(y);
		}
	}


	in.close();
	delete(fn);
	*N = temp.size() / 2 - 1;

if (*N > 0) {
	CStaticState<DIM> current(*N);
	CStdData<DIM> data;
	vector< vector<double> > psi6p;
	Eigen::VectorXd positions(*N * DIM);

	for (int i = 0; i < *N * 2; i++)
	{
		positions[i] = temp[i];
	}

	current.SetPositionsVirtual(positions);

	phi = 0.1;
	double phistep = 0.01;

	
	char statedb_fn[500], datadb_fn[500];
	sprintf(statedb_fn, "MinimizedState%s", filename);
	sprintf(datadb_fn, "MinimizedData%s", filename);

	CStaticDatabase<DIM> state_db(*N, statedb_fn, NcFile::Replace);
	CStdDataDatabase<DIM> data_db(datadb_fn, NcFile::Replace);

	CStaticState<DIM> last(current);
	bool gotToJammed = false;

	while (phistep >= 5e-7)
	{
	
		current.SetRadiiMono();

		current.SetPackingFraction(phi);


		//Minimize the energy:
		CStaticComputer<DIM> computer(current);
		CSimpleMinimizer<DIM> miner(computer, CSimpleMinimizer<DIM>::FIRE); //Note: the miner class should be able to create the Computer (i.e. just pass it the state)

		if (!computer.StdPrepareSystem())
		{
			computer.CalculateStdData();
			computer.CalculatePsi6Mag(psi6p);
		}

		data = computer.Data;

		if (CheckIsJammed(data))
		{
			cout << "\n\n Jammed at phi = " << phi << ", phistep = " << phistep << "\n\n";
			phistep /= 2;
			current = last;
			phi -= phistep;
			gotToJammed = true;
		}
		else 
		{
			cout << "\n\n Not jammed at phi = " << phi << ", phistep = " << phistep << "\n\n";
			if (gotToJammed) phistep /= 2;
			phi += phistep;
			last = current;
		}
	}

	
	// Best estimate for phi_critical
	phi -= phistep / 2;
	cout << "\n\n Phi_C = " << phi << endl << endl << endl;
	cout << "\n\n System is jammed!\n\n\n";
	// Save current state
	state_db.Write(current);
	data_db.Write(data);
	writePsi6Data(fileNameOrig, psi6p, phi);
	double LMinDeltaPhi = -6;
	double LMaxDeltaPhi = -0.2;
	double LStep = 0.2;

	while (LMinDeltaPhi <= LMaxDeltaPhi)
	{
		phi += pow(10, LMinDeltaPhi);
		cout << "\n\n phi = " << phi << "\n\n\n";
		LMinDeltaPhi += LStep;

		current.SetRadiiMono();
		current.SetPackingFraction(phi);
		CStaticComputer<DIM> computer(current);
		CSimpleMinimizer<DIM> miner(computer, CSimpleMinimizer<DIM>::FIRE);

		if (!computer.StdPrepareSystem())
		{
			computer.CalculateStdData();
			computer.CalculatePsi6Mag(psi6p);
		}
		data = computer.Data;


		state_db.Write(current);
		data_db.Write(data);
		writePsi6Data(fileNameOrig, psi6p, phi);
	}
	
}
	delete(fileNameOrig);

}





int main(int argc, char* argv[])
{
	int N = 256;			//n
	dbl phi = 0.9;			//f
	int r = 1;				//r
	int test_num = 1;		//t
	int iter = 0;			//i

	int c;
	while((c=getopt(argc, argv, "n:f:r:t:i:m:")) != -1)
		switch(c)
		{
			case 'n':	N = atoi(optarg); break;
			case 'f':	phi = atof(optarg); break;
			case 'r':	r = atoi(optarg); break;
			case 't':	test_num = atoi(optarg); break;
			case 'i':	iter = atoi(optarg); break;
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


	printf("Begin Diagonalization Tutorial with N=%i, phi=%f, r=%i, iter=%i\n", N, phi, r, iter);


/*#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < iter; i++)
		{
			run(N, phi, r + 9, test_num);
		}
	}
	*/

	int numVoronoiCells[] = {8, 16, 32, 64, 128, 256, 512};
	
/*#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < 20; i++)
		{
		*/
	int i = iter;

	// should be [0,7)
	for (i = iter; i < 20; i++)	
	for (int j = 0; j < 7; j++)
	{	
/*#pragma omp parallel 
		{
#pragma omp for*/
		// should be [1,10)
		for (int numToCellRatioPower = 1; numToCellRatioPower < 10; numToCellRatioPower++)
		{
			int tempN = N;
			int nVCells = numVoronoiCells[j];
			int nTarget = (int) pow(2.0, numToCellRatioPower) * nVCells;
			char* filename = (char *) malloc(500);
			stringstream ss;
			ss << "DataNt_" << nTarget << "nV_" << nVCells << "iter_" << i;
			ss >> filename;
			cout << filename << endl;
			run(&tempN, phi, r + 9, filename, test_num);
			delete(filename);
		
		}
	}
	

}








