// DetlevCBoundScProject.cpp : Defines the entry point for the console application.
//


#include "headers.hpp"

//http://www.cs.hmc.edu/~geoff/classes/hmc.cs070.200401/notes/command-args.html
//char main ( char ,char *argv[]/*Eingabe, int argc, char *argv[]*/)
int main(int argc, char* argv[])
{
	
	/* Some general Counters I use */
	int  i;
	int ii;

	// and a transfer variables for MPIR
	mpz_t Zii;
	mpz_init(Zii);
	//Counter for floating point nuBounders
	mpf_t FloatMii;
	mpf_init(FloatMii);

	mpf_set_default_prec (32); // Set default precision for mpf variables

	/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxx *//* Variables used in my preparation steps */ /* xxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */
	
	// Normal Variables
	std::string Eingabe; // input variable - a standard string, char only takes 1 character

	int o1 = 0; // Refers to the complexity - 0 for the start, might need changing
	int PrimeArray; //a variable to record the size of the prime Array
	int FactorBaseArray; //a variable to record the size of the Factor Base
	
	double Root16thN; // the 16th Root of N to work out my Logarithms, small enough for a double
	double DM; // my "raw M" as it's worked out - hence a double
	double mu; // my mu

	double Prim; // a double to get my primes for conversion into Logs, 15 digits should suffice
	double FactorBaseLogs[5000]; // my array with logarithms of primes

	/* xxxxxxxxxx *//* xxxxxxxxxx *//* xxxxxxxxxx */
	/* xxxxxxxxxx *//* MPIR Stuff *//* xxxxxxxxxx */

	/* XXX - MPZ Types - XXX */
	// initialize N
	mpz_t N;
	mpz_init(N);
	// initialize my variable for floor root N
	mpz_t WurzelN;
	mpz_init(WurzelN);
	// initialize my M sqrt. N
	mpz_t MWurzelN;
	mpz_init(MWurzelN);
	// My "M" for the sieving bounds
	mpz_t MM;
	mpz_init(MM);
	// My Bound for the Factor Base
	mpz_t Bound;
	mpz_init(Bound);

	/* XXX - MPF Types - XXX */
	mpf_t FloatN; //N as a floating point number for a reasonably accurate 16th root
	mpf_init (FloatN);
		mpf_t Root1; // steppingstones to find my 16th Root
	mpf_init (Root1);
	mpf_t Root2;
	mpf_init (Root2);

	/* XXX - My Arrays - XXX */

	// Create an array for primes - 10.000 spaces is enough
	//mpz_t MPrimes[20000];
	vector< mpz_class > MPrimes;
	MPrimes.resize(20000);
	for(i = 0; i <= 19999; i++){  //initialize every value in the array by feeding it trough a loop
		//mpz_init(MPrimes[i]);
		MPrimes[i] = 0; // if the mpz_class is usedm this initialises the variable
	}

	mpz_t MFactorBase[10000]; // this array contains my factor base (Primes that 'passed' the Legendre Symbol)
	for(i = 0; i <= 9999; i++){
		mpz_init(MFactorBase[i]);
	}

	/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxx *//* End Preparation Variables *//* xxxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */


	/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxx *//* Begin Sieving Interval Variables *//* xxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

	double ArrayOfQxLogs[10000]; // the array that contains the ln of my q(x)

	// Steppingstone Variable for working out logs
	mpf_t LogsQx;
	mpf_init(LogsQx);

	// This array contains my Sieving Interval
	mpz_t SievingInterval[10000];
	mpz_t SievingIntervalSquares[10000];
	for(i=0;i<=9999;i++){
		mpz_init(SievingInterval[i]);
		mpz_init(SievingIntervalSquares[i]);
	}

	// Some Stepping Stones
	mpz_t step1;
	mpz_init(step1);
	mpz_t step2;
	mpz_init(step2);

	/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxx *//* End Sieving Interval Variables *//* xxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */


	/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx *//* The actual Program starts here *//* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */
	/* The actual Program starts here *//* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx *//* The actual Program starts here */
	/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx *//* The actual Program starts here *//* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */
	/* The actual Program starts here *//* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx *//* The actual Program starts here */
	/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx *//* The actual Program starts here *//* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

	std::cout << "\n Input a Number: "; // prompt the user for input and accept it
	std::cin >> Eingabe;
	time_t starttime;
	starttime = time (NULL);
	//mpz_set_ui(N, 3277);	// Set my Input to something fixed for testing
	mpz_init_set_str( N , Eingabe.c_str() , 10); //Input function - works, thanks to David
	//mpz_out_str(stdout,10, N); // Testing input
	//printf ("\n")

	mpz_sqrt(WurzelN, N);
	mpz_gcd(Zii, WurzelN ,N);

	if(mpz_cmp (Zii, WurzelN)==0){
		std::cout << "The number you have entered is a perfect square - the root is: " ;
		mpz_out_str(stdout,10, Zii);
		printf("\n");
		// http://www.cplusplus.com/doc/tutorial/control/
		std::cin >> Eingabe;
		return 0; // end the program as the factorization was successful already
	}

	/* Let's work out the sieving interval M */
	mpf_set_z(FloatN,N);

	mpf_sqrt (Root1,FloatN); //sqrt once
	mpf_sqrt (Root2,Root1); // sqrt twice
	mpf_set (Root1,Root2);
	mpf_sqrt (Root2,Root1); // sqrt three times
	mpf_set (Root1,Root2);
	mpf_sqrt (Root2,Root1); // sqrt four times

	Root16thN = mpf_get_d (Root2); // move my root into a double to use logarithms
	
	DM = exp((1+o1)*sqrt((16*log(Root16thN)*log(16*log(Root16thN))))); // this is my M contained in a double

	/* Now I have my M */
	//std::cout << "This is my Rooth16thN: " << Root16thN << "\n" << "This is my DM: " << DM ;

	/* The next bound is B - my factor base bound */
	mu = sqrt((2*(log(DM)+8*log(Root16thN)))/(log((log(DM)+8*log(Root16thN))))); // -> works
	//std::cout << "DM: " << DM << "\n" << "mu: " << mu << "\n"; // Just checking my DM & mu
	
	
	long M = ceil(DM); // round my DM up to an integer
	if(DM > 4999){ // force bound smaller if excessive
		M = 4999;
	}

	mpz_set_ui(MM, M); // move my M to a MPZ variable to handle

	mpz_sqrt(WurzelN, N); //obtain my floor square root of N
	mpz_mul(MWurzelN,MM,WurzelN); // find M sqrt N

	
	long MU = ceil(mu); // rounding mu up, so that numbers don't get too large
	//std::cout << "This is MU " << MU << "\n" ;

	//mpz_nthroot (Bound, MWurzelN , MU); // nth root = power 1/n
	mpz_root (Bound, MWurzelN , MU); // nth root = power 1/n
	long B = mpz_get_si (Bound); // Save Bound as an integer for use in loops
	if(B > 19999){ // force bound smaller if excessive
		mpz_set_ui (Bound, 19999);
		B = 19999;
	}
	
	std::cout << "This is M: " << M << "\n" << "This is B: " << B <<"\n \n"; // Testing my bound

	/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx *//* Let us build an array of primes *//* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */
	
	i = 2; // start my first prime in position 2, position 1 would normally be -1
	//mpz_set_ui(MPrimes[0],1); // set position 1 to 1, just in case -> causes error
	mpz_set_ui (Zii,2);
	mpz_set (MPrimes[1].get_mpz_t (), Zii);

	//loop trough a counter to find primes up to bound B
	for(ii = 2; ii <= B/2 && i<20000 ;ii++){
		mpz_set_ui (Zii,2*ii-1); // set my mpz variable to the counter, only check odd numbers
		//check = mpz_probab_prime_p (Zii , 10);// page 43 MPIR manual - do 10 tests if Zii is prime
		// if the value is definitely prime, write into array
		if(mpz_probab_prime_p (Zii , 10) > 1)
		{
			mpz_set (MPrimes[i].get_mpz_t (),Zii);
			i++; // increase array position by 1
		}
	}
	PrimeArray = i; // retain size of Array of Primes for later use
	/*for(i = 0; i<PrimeArray; i++){ // Testing - yes, I get a list of primes:
		mpz_out_str(stdout,10, MPrimes[i]);
		printf ("\n");
	}*/

	/* Let us optimize our factor base with the Legendre Symbol - page 43 in the manual */
	mpz_set_ui(Zii,2);
	mpz_set (MFactorBase[1], Zii); //a roundabout way of including 2, but easy to do
	i=2;
	for(ii = 2; ii<=PrimeArray && ii <= 19999 ;ii++){
		//std::cout << mpz_legendre(N,MPrimes[ii]);
		if( mpz_legendre(N,MPrimes[ii].get_mpz_t ()) >= 1 ){ // mpz_legendre is 1 when true
			mpz_set (MFactorBase[i],MPrimes[ii].get_mpz_t ());
			i++;
		}
	}
	FactorBaseArray = i; // Retain size of reduced array for later use
	std::cout << "Factor Base Array Size: " << FactorBaseArray << "\n" ;

	/*for(i = 0; i<FactorBaseArray; i++){  // Testing I get a list of my factor base primes
		mpz_out_str(stdout,10, MFactorBase[i]);
		printf ("\n");
	}*/

	for(i=1;i<FactorBaseArray;i++){ // this creates an array of the logs of primes
		//mpz_out_str(stdout,10,MFactorBase[i]); // Testing
		//printf ("\n");
		Prim = mpz_get_d (MFactorBase[i]);
		FactorBaseLogs[i] = log(Prim);
		//std::cout << FactorBaseLogs[i] << "\n" ; // Testing
	}

	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* Preparations  are done now *//* xxxxxxxxxxxxxxxxxxxxxxxxx */
	/* I have found my bounds, used that to determine my factor base and reduced my factor base */
	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx */


	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxx *//* Next on the list - solutions to x^2 equiv 0 modulo p *//* xxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx */

	// MFactorBase[i]

	mpz_t ZiiSquared;
	mpz_init(ZiiSquared);
	mpz_t Nmodp;
	mpz_init(Nmodp);

	mpz_t CongruentX1[5000]; // this array contains my x1^2 equiv N mod P
	mpz_t CongruentX2[5000]; // this array contains my x2^2 equiv N mod P
	for(i = 0; i <= 4999; i++){
		mpz_init(CongruentX1[i]);
		mpz_init(CongruentX2[i]);
	}
	
	for(i=2;i<FactorBaseArray;i++){
		mpz_mod (Nmodp, N , MFactorBase[i]); // find my N mod p

		long Prim = mpz_get_ui(MFactorBase[i]); //getting my prime from the mpz - thanks NBR :)

		int Zeiger = 0; // A pointer - set to 0
		mpz_set_ui(CongruentX2[1],1); // arange for prime = 2
		mpz_set_ui(CongruentX2[1],1);

		for(ii=1;ii<=Prim;ii++){
			mpz_set_ui(Zii,ii);
			mpz_mul (ZiiSquared,Zii,Zii);
			if(mpz_congruent_p (ZiiSquared,Nmodp,MFactorBase[i])){
				if(Zeiger==1){
					mpz_set(CongruentX2[i],Zii); //the first solution
					//std::cout << "Solution 2: " << ii << " - The Prime Number is: " << Prim << "\n"; //testcode
				}
				if(Zeiger==0){
					mpz_set(CongruentX1[i],Zii); //the second solution
					Zeiger=1;
					//std::cout << "Solution 1: " << ii << " - The Prime Number is: " << Prim << "\n"; //testcode
				}
			}
		}
	}

	/*for(i = 0; i<30; i++){ //test code
		mpz_out_str(stdout,10,CongruentX1[i]); // solution 1
		printf ("\n");
		mpz_out_str(stdout,10,CongruentX2[i]); // solution 2
		printf ("\n");
	}*/

	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxx *//*End - solutions to x^2 equiv 0 modulo p *//* xxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx */


	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxx *//* Next step - work out my sieving interval *//* xxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* As well as the Logarithms *//* xxxxxxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx */

	for(i=0;i<=M;i++){ // this works out root n - m
		mpz_set_ui (Zii,M-i);
		mpz_sub(step1,WurzelN,Zii);
		mpz_set(SievingInterval[i],step1);
	}
	for(i=M+1;i<=2*M;i++){ // this works out root n + m
		mpz_set_ui (Zii,i-M);
		mpz_add(step1,WurzelN,Zii);
		mpz_set(SievingInterval[i],step1);
	}

	for(i=0;i<=M;i++){ // this does the negative values and records them positively
		mpz_set_ui (Zii,M-i);
		mpz_sub(step1,WurzelN,Zii);
		mpz_mul (step2,step1,step1);
		mpz_sub(step1,N,step2);
		mpz_set(SievingIntervalSquares[i],step1);
	}
	for(i=M+1;i<=2*M;i++){ // this works out the positive values
		mpz_set_ui (Zii,i-M);
		mpz_add(step1,WurzelN,Zii);
		mpz_mul (step2,step1,step1);
		mpz_sub(step1,step2,N);
		mpz_set(SievingIntervalSquares[i],step1);
	}

	/*for(i = M-10; i<M+10 ;i++){ // Just testing
		mpz_out_str(stdout,10,SievingIntervalSquares[i]);
		printf ("\n");
	}*/

	// Let's get the logs of q(x)
	for(i=0;i<=2*M;i++){
		mpf_set_z(LogsQx,SievingIntervalSquares[i]);

		mpf_sqrt (Root1,LogsQx); //sqrt once
		mpf_sqrt (Root2,Root1); // sqrt twice
		mpf_set (Root1,Root2);
		mpf_sqrt (Root2,Root1); // sqrt three times
		mpf_set (Root1,Root2);
		mpf_sqrt (Root2,Root1); // sqrt four times

		double step = mpf_get_d (Root2);
		//std::cout << 16*log(step) << " - " << i << "\n";
		ArrayOfQxLogs[i] = 16*log(step); // move my root into a double to use logarithms
	}

	/*for(i = M-10; i<M+10 ;i++){ // Just testing
		std::cout << ArrayOfQxLogs[i] << "\n";
	}*/

	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxx *//* End working out my sieving interval *//* xxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx */


	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* Let us start the sieving *//* xxxxxxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx */

	int iii;

	std::vector<std::vector<int>> FirstPowerSieving(5000 , std::vector<int> (10000));

	//int* FirstPowerSieving = new int[5000][10000];

	/*for(i=0;i<=4999;i++){
		for(ii=0;ii<=9999;ii++){
			mpz_init(FirstPowerSieving(i)(ii));
		}
	}*/

	//mpz_congruent_p (SievingInterval[ii],CongruentX1[i],MFactorBase[i]); // first solution
	//mpz_congruent_p (SievingInterval[ii],CongruentX2[i],MFactorBase[i]); // second solution

	for(i=1;i<FactorBaseArray;i++){ // I need to sieve through all my primes
		for(ii=0;ii<Prim && ii< 2*M;ii++){ // check a range of "prime" for congruence
			long LongPrim = mpz_get_ui (MFactorBase[i]);
			if(0<mpz_congruent_p (SievingInterval[ii],CongruentX1[i],MFactorBase[i])){ //write a 1 to array where divisible
				for(iii=ii;iii<2*M;){
					FirstPowerSieving[i][iii]=1;
					iii = iii+LongPrim;
				}
			}
			if(0<mpz_congruent_p (SievingInterval[ii],CongruentX2[i],MFactorBase[i])){
				for(iii=ii;iii<2*M;){
					FirstPowerSieving[i][iii]=1;
					iii = iii+LongPrim;
				}
			}
		}
	}
	// Let's add a marker for -1 where the value should be negative:
	for(i=0;i<=M;i++){
			FirstPowerSieving[0][i]=1;
	}

	/*for(i=0;i<FactorBaseArray;i++){ // testing my first power sieving, works.
		std::cout << " -- ";
		for(ii=M-10;ii<M+10;ii++){
			std::cout << FirstPowerSieving[i][ii];
		}
		std::cout << "\n";
	}*/

	long PrimPower;

	for(i=1;i<FactorBaseArray;i++){
		Prim = mpz_get_d (MFactorBase[i]);
		//std::cout << "\n" << Prim << " - " ;
		for(ii=0;ii<2*M;ii++){
			if(FirstPowerSieving[i][ii]>=1){ // check whether it divides by p^1 -> no need to test the ones I know aren't divisible
				iii = 0; // start with zero, as the first step raises the power to 1
				do{
					iii = iii+1;
					PrimPower = pow(Prim,iii);// work out my prime^something for the modulo
					//std::cout << PrimPower << " - " << mpz_divisible_ui_p (SievingIntervalSquares[ii],PrimPower) << " - " << mpz_get_ui(SievingIntervalSquares[ii]) <<"\n";
				//}while(0<mpz_congruent_ui_p (SievingIntervalSquares[ii],0,PrimPower)); // check if congruent
				// testing for divisibility or congruence in modulo, both work equally well,both fail at powers of 2
				}while(mpz_divisible_ui_p (SievingIntervalSquares[ii],PrimPower) && iii < 11); // check if civisible
				FirstPowerSieving[i][ii]=iii-1; // set new power, - 1, as the last evaluated value would be 1 too big
				//std::cout << FirstPowerSieving[i][ii];
			}
		}
	}
	
	/*for(i=0;i<FactorBaseArray;i++){ // testing my prime power finder 
		std::cout << "\n" << mpz_get_d (MFactorBase[i]) << " - ";
		for(ii=M-10;ii<M+10;ii++){
			std::cout << FirstPowerSieving[i][ii];
		}
	}*/

	for(i=1;i<FactorBaseArray;i++){ // subtracting from the logs, the ones close to zero will be smooth
		//std::cout << FactorBaseLogs[i] << " - ";
		for(ii=0;ii<2*M;ii++){
			//std::cout  << FirstPowerSieving[i][ii];
			ArrayOfQxLogs[ii] = ArrayOfQxLogs[ii]-(FactorBaseLogs[i]*FirstPowerSieving[i][ii]);
		}
		//std::cout << "\n";
	}

	/*for(i=0;i<2*M;i++){
		std::cout << ArrayOfQxLogs[i] << " - " << i << "\n";
	}*/
	
	std::cout << "Logartihms Calculated \n";
	
	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxxxxxxx *//* Let us end the sieving *//* xxxxxxxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx */


	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxx *//* Picking my values and using them *//* xxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx */

	std::vector<std::vector<int>> Matrix(5000 , std::vector<int> (10000)); 
	// same size as first matrix, even though this is excessive
	int mapone[10000];
	for(i=0;i<10000;i++){
		mapone[i]=0;
	}
	int ReducedMatrix;

	iii = 0;
	
	for(i=0;i<2*M;i++){
		if(abs(ArrayOfQxLogs[i])<1){
			for(ii=0;ii<FactorBaseArray;ii++){
				Matrix[ii][iii]=FirstPowerSieving[ii][i];
			}
			mapone[iii]=i; // keep track of the position mappings
			iii = iii+1;
		}
	}

	ReducedMatrix = iii; // a variable for my reduced matrix "width" count

	/*for(i=0;i<FactorBaseArray;i++){ // let's test our reduced matrix
		for(ii=0;ii<ReducedMatrix;ii++){
			std::cout << Matrix[i][ii];
		}
		std::cout << "\n";
	}*/

	// now let us change the values to whatever they are in modulo 2:
	for(i=0;i<FactorBaseArray;i++){
		for(ii=0;ii<ReducedMatrix;ii++){
			Matrix[i][ii]=Matrix[i][ii] % 2;
		}
	}

	std::cout << "First column reduced matrix modulo 2 created \n";

	/*for(i=0;i<FactorBaseArray;i++){ // let's test our modulo reduction
		for(ii=0;ii<ReducedMatrix;ii++){
			std::cout << Matrix[i][ii];
		}
		std::cout << "\n";
	}*/

	/* Gaussian Elimination... or maybe some more reducing... */
	// More reducing - vectors of zero can be removed

	
	int maptwo[10000];
	for(i=0;i<10000;i++){
		maptwo[i]=0;
	}
	std::vector<std::vector<int>> MatrixReducedTwo(5000 , std::vector<int> (10000));
	int counter;
	int ReducedMatrixTwo;
	
	iii = 0;
	for(ii=0;ii<ReducedMatrix;ii++){
		counter = 0;
		for(i=0;i<FactorBaseArray;i++){
			counter = counter + Matrix[i][ii];
		}
		if(counter > 0){
			for(i=0;i<FactorBaseArray;i++){
				MatrixReducedTwo[i][iii] = Matrix[i][ii];
			}
		maptwo[iii]=mapone[ii];
		iii = iii +1;
		}
	}

	ReducedMatrixTwo = iii;

	std::cout << "Second column reduced matrix modulo 2 created \n";

	/*for(i=0;i<FactorBaseArray;i++){ // let's test our reduced matrix
		for(ii=0;ii<ReducedMatrixTwo;ii++){
			std::cout << MatrixReducedTwo[i][ii];
		}
		std::cout << "\n";
	}*/

	// let us get rid of empty rows in the matrix

	int mapprimes[5000];
	for(i=0;i<5000;i++){
		mapprimes[i]=0;
	}
	std::vector<std::vector<int>> MatrixReducedThree(5000 , std::vector<int> (10000));
	int ReducedMatrixRows;
	iii = 0;
	
	for(i=0;i<FactorBaseArray;i++){
		counter = 0;
		for(ii=0;ii<ReducedMatrixTwo;ii++){
			counter = counter + Matrix[i][ii];
		}
		if(counter > 0){
			for(ii=0;ii<ReducedMatrixTwo;ii++){
				MatrixReducedThree[iii][ii] = Matrix[i][ii];
			}
		mapprimes[iii]=i;
		iii = iii +1;
		}
	}

	ReducedMatrixRows = iii;

	std::cout << "Created row-reduced Matrix \n";

	/*for(i=0;i<ReducedMatrixRows;i++){ // let's test our reduced matrix
		std::cout << mapprimes[i] << " - ";
		for(ii=0;ii<ReducedMatrixTwo;ii++){
			std::cout << MatrixReducedThree[i][ii];
		}
		std::cout << "\n";
	}*/

	// let us clear the vectors we no longer need

	// These contain the first two steps in the reduction and are now no longer needed
	// this should also reduce memory use
	MatrixReducedTwo.clear();
	MatrixReducedTwo.clear();

	// we still need the original sieving/factor finding matrix as it contains my powers and will save us
	// factorization efforts a bit later on
	// -> source: http://msdn.microsoft.com/en-us/library/fs5a18ce%28v=vs.80%29.aspx

	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxx *//* End picking my vlaues and using them *//* xxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx */


	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxx *//* Gaussian Elimination *//* xxxxxxxxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx */

	// MatrixReducedThree[iii][ii] <- reduced small matrix
	// ReducedMatrixRows <- Number of rows
	// ReducedMatrixTwo  <- Number of columns

	
	ii = 0;
	int check;
	int selection[10000];
	int row[5000];

	for(i=0;i<4999;i++){
		selection[2*i+1]=0;
		selection[2*i]=0;
		row[i]=0;
	}

	for(i=0;i<ReducedMatrixTwo;i++){ // over all columns
		ii = 0;
		do{
			if(MatrixReducedThree[ii][i] == 1 && row[ii] != 1){ // check which row contains a 1
				check = 1;
				selection[i] = ii; // remember which row is used in which column
				row[ii] = 1;
			}
			ii = ii +1 ;
		}while(check == 0 && ii < ReducedMatrixRows);

		if(check == 1){
		for(ii=0;ii<ReducedMatrixRows;ii++){ // over all rows
			
			if(MatrixReducedThree[ii][i] == 1 && selection[i] != ii){
				if(i==0){
				for(iii=0;iii<ReducedMatrixTwo;iii++){ // clumsy way of doing it, split it...
					MatrixReducedThree[ii][iii]=(MatrixReducedThree[ii][iii]+MatrixReducedThree[selection[i]][iii]) % 2; // I need to make this run over all columns
				}
				}

				if(i != 0 && selection[i] != 0){
					for(iii=0;iii<ReducedMatrixTwo;iii++){
					MatrixReducedThree[ii][iii]=(MatrixReducedThree[ii][iii]+MatrixReducedThree[selection[i]][iii]) % 2;
				}
				}
			}
		}
		}
		check = 0;
	}

	std::cout << "Gaussian Elimination complete - The matrix has " << ReducedMatrixRows << " rows and " << ReducedMatrixTwo << " columns \n\n";

	/*for(i=0;i<ReducedMatrixRows;i++){ // let's test our reduced matrix
		std::cout << mapprimes[i] << " - ";
		for(ii=0;ii<ReducedMatrixTwo;ii++){
			std::cout << MatrixReducedThree[i][ii];
		}
		std::cout << "\n";
	}*/

	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* End Gaussian Elimination *//* xxxxxxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx */

	// redo with another congruence if factor is 1
	int Verification;
	Verification = 0;
	int Blacklist; //this will aid in sorting out non-useful vectors
	Blacklist = 1;

	// in case my first congruence gives no factor...
	do{

	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx *//* Linear Dependency *//* xxxxxxxxxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx */


	/* Variables for this section */
	// chosing vector
	int Vektorwahl; //this contains the sum of the elements in the vector -> if it's zero I have a zero vector, which is of little use

	// searchign for the right pivots
	int NumberSelectedVectors;
	int ChosenVectors[10000];
	for(i = 0; i <= 9999; i++){
		ChosenVectors[i]=0;
	}
	iii = 0;

	/* The idea: Pick a vector from the right of the matrix, then select pivots from the right
	finally add all elements together, rows in modulo 2, they should add up to zero */
	
	i = ReducedMatrixTwo-Blacklist; // starting at the right hand side of my matrix

	do{
		Vektorwahl = 0;
		for(ii=0;ii<ReducedMatrixRows;ii++){
			Vektorwahl = Vektorwahl + MatrixReducedThree[ii][i];
			//swapped i and ii round, I need to go trhough the rows not columns as I did before by accident
		}
		i = i-1;
	}while(Vektorwahl == 0); //this do loop will find a non-zero vector, starting from the right hand side of the matrix

	Vektorwahl = i; // -> Now "Vektorwahl" contains the "ID" of the selected Vector

	std::cout << "The selected vector has position " << Vektorwahl << "\n"; // Just testing & feedback for user, works

	// Next step: Find the pivots that fit "our" vector and will add to a zero vector in mod 2

	for(i=0;i<ReducedMatrixRows;i++){
		//std::cout << i << "\n";
		if(MatrixReducedThree[i][Vektorwahl] == 1){
			check = 0;
			ii = 0;

			do{
				if(selection[ii] == i){ // search for the columns with pivots in the required rows
					check = 1;
					ChosenVectors[iii] = ii;
					iii = iii+1;
					//std::cout << iii << "\n";
				}
				ii = ii + 1;
			}while(check == 0 && ii < ReducedMatrixTwo);
		}
	}

	iii = NumberSelectedVectors;

	std::cout << "Vectors Selected \n";
	//std::cout << "The Number of Selected Vectors is " << NumberSelectedVectors << "\n" ; // just testing

	/*for(i=0;i<NumberSelectedVectors;i++){ // Just testing
		std::cout << ChosenVectors[i] << "\n";
	}*/

	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxxxxxxxx *//* End Linear Dependency *//* xxxxxxxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx */


	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxx *//* Now let us make use of our vectors *//* xxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx */

	/* Variables for this section: */

	// maptwo[] contains the mappings for my "M" which I need to select the right ones
	int choices[5000];
	int choice;

	int ChoicePrimePowers[5000];
	for(i=0;i<4999;i++){
		ChoicePrimePowers[i]=0;
	}
	int ChoicePrimePowerArray;

	mpz_t LargeX;
	mpz_init(LargeX);
	mpz_set_ui(LargeX,1); // set one for later multiplications, an empty variable with zero wouldn't work
	mpz_set_ui(MFactorBase[0],1);
	
	// my large y, 1 for multiplication
	mpz_t LargeY;
	mpz_init(LargeY);
	mpz_set_ui(LargeY,1);

	/* End variables for this section */

	ii = 1;
	choices[0]=maptwo[ChosenVectors[0]]; //this case must be handled individually because of the later zeroes...
	for(i=1;i<NumberSelectedVectors;i++){
		if(ChosenVectors[i] > 0){
			choices[ii]=maptwo[ChosenVectors[i]];
			ii = ii + 1;
		}
	}
	choices[ii]=maptwo[Vektorwahl];

	choice = ii+1;
	//std::cout << choice << "\n"; // Just checking, works

	/*for(i=0;i<choice;i++){
		std::cout << "My positions are " << i << " - " << choices[i] << "\n";
	}*/

	// calculating my "LargeX"
	for(i=0;i<choice;i++){
		mpz_mul(LargeX,LargeX,SievingInterval[choices[i]]);
	}

	// Now the prime powers
	iii = 0;
	for(i=0;i<FactorBaseArray;i++){
		for(ii = 0;ii<choice;ii++){
			ChoicePrimePowers[iii] = ChoicePrimePowers[i]+FirstPowerSieving[i][choices[ii]];
		}
		iii = iii + 1;
	}
	ChoicePrimePowerArray = iii;
	// now I have all my powers
	//I need to halve all of them

	for(i=0;i<ChoicePrimePowerArray;i++){
		ChoicePrimePowers[i] = ChoicePrimePowers[i]/2; //half the values
	}
	
	// testing prime powers
	/*for(i=0;i<ChoicePrimePowerArray;i++){
		mpz_out_str(stdout,10, MFactorBase[i]);
		std::cout<< " - " << ChoicePrimePowers[i] << "\n";
	}*/

	for(i=1;i<ChoicePrimePowerArray;i++){
		mpz_pow_ui(Zii,MFactorBase[i],ChoicePrimePowers[i]);
		mpz_mul(LargeY,LargeY,Zii);
	}

	// Let's do the mod stuff
	mpz_mod(LargeX,LargeX,N);
	mpz_mod(LargeY,LargeY,N);
	
	/*mpz_sub(Zii,LargeX,LargeY);
	if(mpz_get_ui(Zii) == 0 ){
		Verification = 1;
		std::cout << "LargeX and LargeY are equal \n";
		std::cout << "\n Large X - Large Y N is 0 ";
		Verification = 1;
	}*/
	/*std::cout << "\n Large X - Large Y N is: ";
	mpz_out_str(stdout,10, Zii);
	printf ("\n");*/

	// Output LargeX and LargeY for the user to see

	std::cout << "\n The Large X in modulo N is: ";
	mpz_out_str(stdout,10, LargeX);
	printf ("\n");

	std::cout << "\n The Large Y in modulo N is: ";
	mpz_out_str(stdout,10, LargeY);
	printf ("\n");

	// subtract Y from X
	mpz_sub(Zii,LargeX,LargeY);
	// fidn the greatest common divisor with N
	mpz_gcd(Zii,Zii,N);

	// I need to check whether Zii is one, if yes, redo with another congruence

	if(mpz_get_ui(Zii)!=1){ //if factor is nontrivial
	Verification = 0;
	std::cout << "\n A factor has been found \n";
	std::cout << " The first factor is ";
	mpz_out_str(stdout,10, Zii); // Output the result
	//printf ("\n");
	mpz_divexact(Zii,N,Zii);
	std::cout << "\n The second factor is ";
	mpz_out_str(stdout,10, Zii);
	printf ("\n");
	time_t endtime;
	endtime = time (NULL);
	std::cout << " \n The program took " << endtime - starttime<< " seconds to run." << "\n";
	}
	else{ //if factor is trivial, i.e. 1
	std::cout << "The " << Blacklist << " congruence has not resulted in a non-trivial factor. \n \n";
	Verification = 1;
	Blacklist = Blacklist + 1;
	}

	// closing do for loop in case of no factor
	}while(Verification != 0 && Blacklist < ReducedMatrixTwo-1);
	
	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* End making  use of vectors *//* xxxxxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx */


	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxx *//* Clear no longer needed variables *//* xxxxxxxxxxxxxxxxxxxxxx */
	/* xxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxxx *//* xxxxxxxxxxxxxxxxxxxxxxxxx */

	/*
	for(i = 0; i <= 19999; i++){
		mpz_clear(MPrimes[i]);
	}//*/
	
	mpz_clear(step1);
	mpz_clear(step2);
		
	mpf_clear(Root1);
	mpf_clear(Root2);

	mpz_clear(Zii);
	mpz_clear(N);
	mpz_clear(WurzelN);

	//Trying to keep the window open by requesting another input
	std::cin >> Eingabe;
	return 0;
}


