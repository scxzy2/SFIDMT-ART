#include <math.h>
#include <fstream> 
#include <iostream>
#include <random>
#include <ctime>
using namespace std;
const double EPS = 1e-6; //Epsilon
const double InputDomain = 1000.0;  //The size of the entire input domain
const double Top = InputDomain / 3, Down = 0.0; //Upper and lower bounds of the entire input domain
const double centre = (Top + Down) / 2; //Centre of the entire input domain
const int candidateNum = 10; //The number of source candidates generated in each iteration
//ofstream output1; //Outputs
uniform_real_distribution<double> uRandom(Down, Top); //Used to generate source candidates
vector<double> Executedx, ExecutedxSource, ExecutedxFollow; //Executed test sets
vector<int> CandidateSubdomainList, CandidateSubdomain, farthestCandidate; //Subdomain sets and distance sets
double range = 0, Position = 0; //Size of subdomains and position of test cases in the input domain
double Candidate[candidateNum]; //Source candidate set
double candidateSTC[candidateNum]; //Source candidate set
double x1 = 0, x2 = 0, ret1 = 0, ret2 = 0; //Source input, follow-up input, source output and follow-up output
int Random = 0, N = 0; //Represent the subdomain that is selected for new source candidate generation

//Program tanh (Origin)
double tanh1(double x) {
	double epu, emu, dum, tanh1, x2;
	epu = exp(x);
	emu = 1.0 / epu; //correct
	//emu = 10 / epu; //mutant 1
	if (abs(x) < 0.1) { //correct
	//if (abs(x) < 3.1) { //mutant 2
		x2 = x * x;
		dum = 1 + x2 / 6 * (1 + x2 / 20 * (1 + x2 / 42 * (1 + x2 / 72)));
		tanh1 = 2 * x * dum / (epu + emu);
	}
	else
		tanh1 = (epu - emu) / (epu + emu);
	return tanh1;
}

//Program tanh (mutant 1)
double tanh1Mutant1(double x) {
	double epu, emu, dum, tanh1, x2;
	epu = exp(x);
	//emu = 1.0 / epu; //correct
	emu = 10 / epu; //mutant 1
	if (abs(x) < 0.1) { //correct
	//if (abs(x) < 3.1) { //mutant 2
		x2 = x * x;
		dum = 1 + x2 / 6 * (1 + x2 / 20 * (1 + x2 / 42 * (1 + x2 / 72)));
		tanh1 = 2 * x * dum / (epu + emu);
	}
	else
		tanh1 = (epu - emu) / (epu + emu);
	return tanh1;
}

//Program tanh (mutant 2)
double tanh1Mutant2(double x) {
	double epu, emu, dum, tanh1, x2;
	epu = exp(x);
	emu = 1.0 / epu; //correct
	//emu = 10 / epu; //mutant 1
	//if (abs(x) < 0.1) { //correct
	if (abs(x) < 3.1) { //mutant 2
		x2 = x * x;
		dum = 1 + x2 / 6 * (1 + x2 / 20 * (1 + x2 / 42 * (1 + x2 / 72)));
		tanh1 = 2 * x * dum / (epu + emu);
	}
	else
		tanh1 = (epu - emu) / (epu + emu);
	return tanh1;
}

//Initialize the test sets
void initialize() {
	Executedx.clear(); ExecutedxSource.clear(); ExecutedxFollow.clear();
}

//Random Testing
int RT(int randomNumber, int mutantNum) {
	default_random_engine e1(randomNumber);
	for (int n = 1; n <= 100000; n++) {
		//generate STCs and FTCs
		x1 = uRandom(e1);
		x2 = 3 * x1;
		//execution
		if (mutantNum == 1) {
			ret1 = tanh1Mutant1(x1);
			ret2 = tanh1Mutant1(x2);
		}
		else {
			ret1 = tanh1Mutant2(x1);
			ret2 = tanh1Mutant2(x2);
		}
		if (fabs((ret1 * ret1 * ret1 + 3 * ret1) / (1 + 3 * ret1 * ret1) - ret2) > EPS) {
			//output1 << "Violate: " << x1 << " " << x2 << endl; //output the STC and FTC that violate the MR
			return n;
		}
		//output1 << x1 << " " << x2 << endl; //output the generated STCs and FTCs
	}
	return 100000;
}

int MTARTMin(int randomNumber, int mutantNum) {
	// Random generator seed
	default_random_engine e1(randomNumber);
	unsigned seed = randomNumber;
	srand(seed);
	initialize();
	for (int n = 0; n < 100000; n++) {
		CandidateSubdomainList.clear();
		if (n == 0) { //randomly generate the first test case
			x1 = uRandom(e1);
		}
		else {
			//calculate the size of subdomains
			range = (Top - Down) / (Executedx.size() + 1.0);
			//generate subdomains
			CandidateSubdomain.clear();
			CandidateSubdomain.assign(Executedx.size() + 1, 0);
			//map executed test cases into the subdomains
			for (int i = 0; i < ExecutedxSource.size(); i++) {
				Position = fabs(ExecutedxSource[i] - centre);
				CandidateSubdomain[floor(Position / range * 2)] = 1;
			}
			for (int i = 0; i < Executedx.size() + 1; i++)
				if (CandidateSubdomain[i] == 0)
					CandidateSubdomainList.push_back(i);
			//generate candidates in the target (empty) subdomain
			for (int m = 0; m < candidateNum; m++) {
				N = CandidateSubdomainList[rand() % CandidateSubdomainList.size()];
				Random = rand() % 2;
				//centre
				if (N == 0) {
					uniform_real_distribution<double> ux(centre - range * (N + 1.0) / 2, centre + range * (N + 1.0) / 2);
					Candidate[m] = ux(e1);
				}
				//left
				else if (Random == 0) {
					uniform_real_distribution<double> ux(centre - range * (N + 1.0) / 2, centre - range * (N) / 2);
					Candidate[m] = ux(e1);
				}
				//right
				else if (Random == 1) {
					uniform_real_distribution<double> ux(centre + range * (N) / 2, centre + range * (N + 1.0) / 2);
					Candidate[m] = ux(e1);
				}
			}
			double distanceNearestSTC = 0.0, distanceNearestFTC = 0.0, distanceCurrent = 0.0, distanceNearest = 0.0;
			//select the test case with farthest distance
			for (int p = 0; p < candidateNum; p++) {
				//compute the distances
				for (int q = 0; q < Executedx.size(); q++) {
					distanceCurrent = fabs(Executedx[q] - Candidate[p]);
					if (q == 0)
						distanceNearestSTC = distanceCurrent;
					else if (distanceNearestSTC > distanceCurrent)
						distanceNearestSTC = distanceCurrent;
				}
				for (int q = 0; q < Executedx.size(); q++) {
					distanceCurrent = fabs(Executedx[q] - (3 * Candidate[p]));
					if (q == 0)
						distanceNearestFTC = distanceCurrent;
					else if (distanceNearestFTC > distanceCurrent)
						distanceNearestFTC = distanceCurrent;
				}
				//select the smaller one
				if (distanceNearestSTC < distanceNearestFTC)
					distanceNearestSTC = distanceNearestSTC;
				else
					distanceNearestSTC = distanceNearestFTC;
				//compare the distances
				if (p == 0) {
					distanceNearest = distanceNearestSTC;
					farthestCandidate.clear();
					farthestCandidate.push_back(p);
				}
				else {
					if (distanceNearest < distanceNearestSTC) {
						distanceNearest = distanceNearestSTC;
						farthestCandidate.clear();
						farthestCandidate.push_back(p);
					}
					else if (fabs(distanceNearest - distanceNearestSTC) < EPS)
						farthestCandidate.push_back(p);
				}
			}
			//If more than one candidate MG is selected, consider the distances between source and follow-up candidates
			int farthestCandidateNum = 0; double distanceFurthestAll = 0.0;
			if (farthestCandidate.size() > 1) {
				distanceFurthestAll = 0.0;
				for (int p = 0; p < farthestCandidate.size(); p++) {
					distanceCurrent = fabs(2 * Candidate[farthestCandidate[p]]);
					if (distanceFurthestAll < distanceCurrent) {
						distanceFurthestAll = distanceCurrent;
						farthestCandidateNum = p;
					}
				}
			}
			//select one source candidate for execution
			x1 = Candidate[farthestCandidate[farthestCandidateNum]];
		}
		//generate FTCs
		x2 = 3 * x1;
		//execution
		if (mutantNum == 1) {
			ret1 = tanh1Mutant1(x1);
			ret2 = tanh1Mutant1(x2);
		}
		else {
			ret1 = tanh1Mutant2(x1);
			ret2 = tanh1Mutant2(x2);
		}
		if (fabs((ret1 * ret1 * ret1 + 3 * ret1) / (1 + 3 * ret1 * ret1) - ret2) > EPS) {
			//output1 << "Violate: " << x1 << " " << x2 << endl; //output the STC and FTC that violate the MR
			return n;
		}
		//output1 << x1 << " " << x2 << endl; //output the generated STCs and FTCs
		//Add the executed test cases into the corresponding test sets
		//Includes all executed STCs and FTCs
		Executedx.push_back(x1);
		Executedx.push_back(x2);
		//Includes the executed STCs and FTCs that are inside the source input domain
		ExecutedxSource.push_back(x1);
		if (3 * x1 < Top) {
			ExecutedxSource.push_back(x2);
		}
	}
	return 100000;
}

int MTARTMax(int randomNumber, int mutantNum) {
	//random generator seed
	default_random_engine e1(randomNumber);
	unsigned seed = randomNumber;
	srand(seed);
	initialize();
	for (int n = 0; n < 100000; n++) {
		CandidateSubdomainList.clear();
		if (n == 0) { //randomly generate the first test case
			x1 = uRandom(e1);
		}
		else {
			//calculate the size of subdomains
			range = (Top - Down) / (Executedx.size() + 1.0);
			//generate subdomains
			CandidateSubdomain.clear();
			CandidateSubdomain.assign(Executedx.size() + 1, 0);
			//map executed test cases into the subdomains
			for (int i = 0; i < ExecutedxSource.size(); i++) {
				Position = fabs(ExecutedxSource[i] - centre);
				CandidateSubdomain[floor(Position / range * 2)] = 1;
			}
			for (int i = 0; i < Executedx.size() + 1; i++)
				if (CandidateSubdomain[i] == 0)
					CandidateSubdomainList.push_back(i);
			//generate candidates in the target (empty) subdomain
			for (int m = 0; m < candidateNum; m++) {
				N = CandidateSubdomainList[rand() % CandidateSubdomainList.size()];
				Random = rand() % 2;
				//centre
				if (N == 0) {
					uniform_real_distribution<double> ux(centre - range * (N + 1.0) / 2, centre + range * (N + 1.0) / 2);
					Candidate[m] = ux(e1);
				}
				//left
				else if (Random == 0) {
					uniform_real_distribution<double> ux(centre - range * (N + 1.0) / 2, centre - range * (N) / 2);
					Candidate[m] = ux(e1);
				}
				//right
				else if (Random == 1) {
					uniform_real_distribution<double> ux(centre + range * (N) / 2, centre + range * (N + 1.0) / 2);
					Candidate[m] = ux(e1);
				}
			}
			double distanceNearestSTC = 0.0, distanceNearestFTC = 0.0, distanceCurrent = 0.0, distanceNearest = 0.0;
			//select the test case with farthest distance
			for (int p = 0; p < candidateNum; p++) {
				//compute the distance
				for (int q = 0; q < Executedx.size(); q++) {
					distanceCurrent = fabs(Executedx[q] - Candidate[p]);
					if (q == 0)
						distanceNearestSTC = distanceCurrent;
					else if (distanceNearestSTC < distanceCurrent)
						distanceNearestSTC = distanceCurrent;
				}
				for (int q = 0; q < Executedx.size(); q++) {
					distanceCurrent = fabs(Executedx[q] - (3 * Candidate[p]));
					if (q == 0)
						distanceNearestFTC = distanceCurrent;
					else if (distanceNearestFTC < distanceCurrent)
						distanceNearestFTC = distanceCurrent;
				}
				//select the larger one
				if (distanceNearestSTC < distanceNearestFTC)
					distanceNearestSTC = distanceNearestFTC;
				else
					distanceNearestSTC = distanceNearestSTC;
				//compare the distances
				if (p == 0) {
					distanceNearest = distanceNearestSTC;
					farthestCandidate.clear();
					farthestCandidate.push_back(p);
				}
				else {
					if (distanceNearest < distanceNearestSTC) {
						distanceNearest = distanceNearestSTC;
						farthestCandidate.clear();
						farthestCandidate.push_back(p);
					}
					else if (fabs(distanceNearest - distanceNearestSTC) < EPS)
						farthestCandidate.push_back(p);
				}
			}
			//If more than one candidate MG is selected, consider the distances between source and follow-up candidates
			int farthestCandidateNum = 0; double distanceFurthestAll = 0.0;
			if (farthestCandidate.size() > 1) {
				distanceFurthestAll = 0.0;
				for (int p = 0; p < farthestCandidate.size(); p++) {
					distanceCurrent = fabs(2 * Candidate[farthestCandidate[p]]);
					if (distanceFurthestAll < distanceCurrent) {
						distanceFurthestAll = distanceCurrent;
						farthestCandidateNum = p;
					}
				}
			}
			//select one source candidate for execution
			x1 = Candidate[farthestCandidate[farthestCandidateNum]];
		}
		//generate FTCs
		x2 = 3 * x1;
		//execution
		if (mutantNum == 1) {
			ret1 = tanh1Mutant1(x1);
			ret2 = tanh1Mutant1(x2);
		}
		else {
			ret1 = tanh1Mutant2(x1);
			ret2 = tanh1Mutant2(x2);
		}
		if (fabs((ret1 * ret1 * ret1 + 3 * ret1) / (1 + 3 * ret1 * ret1) - ret2) > EPS) {
			//output1 << "Violate: " << x1 << " " << x2 << endl; //output the STC and FTC that violate the MR
			return n;
		}
		//output1 << x1 << " " << x2 << endl; //output the generated STCs and FTCs
		//add the executed test cases into the corresponding test sets
		//include all executed STCs and FTCs
		Executedx.push_back(x1);
		Executedx.push_back(x2);
		//include the executed STCs and FTCs that are inside the source input domain
		ExecutedxSource.push_back(x1);
		if (3 * x1 < Top) {
			ExecutedxSource.push_back(x2);
		}
	}
	return 100000;
}

int MTARTAve(int randomNumber, int mutantNum) {
	//random generator seed
	default_random_engine e1(randomNumber);
	unsigned seed = randomNumber;
	srand(seed);
	initialize();
	for (int n = 0; n < 100000; n++) {
		CandidateSubdomainList.clear();
		if (n == 0) { //randomly generate the first test case
			x1 = uRandom(e1);
		}
		else {
			range = (Top - Down) / (Executedx.size() + 1.0);
			//generate subdomains
			CandidateSubdomain.clear();
			CandidateSubdomain.assign(Executedx.size() + 1, 0);
			//map executed test cases into the subdomains
			for (int i = 0; i < ExecutedxSource.size(); i++) {
				Position = fabs(ExecutedxSource[i] - centre);
				CandidateSubdomain[floor(Position / range * 2)] = 1;
			}
			for (int i = 0; i < Executedx.size() + 1; i++)
				if (CandidateSubdomain[i] == 0)
					CandidateSubdomainList.push_back(i);
			//generate candidates in the target (empty) subdomain
			for (int m = 0; m < candidateNum; m++) {
				N = CandidateSubdomainList[rand() % CandidateSubdomainList.size()];
				Random = rand() % 2;
				//centre
				if (N == 0) {
					uniform_real_distribution<double> ux(centre - range * (N + 1.0) / 2, centre + range * (N + 1.0) / 2);
					Candidate[m] = ux(e1);
				}
				else if (Random == 0) {
					uniform_real_distribution<double> ux(centre - range * (N + 1.0) / 2, centre - range * (N) / 2);
					Candidate[m] = ux(e1);
				}
				else if (Random == 1) {
					uniform_real_distribution<double> ux(centre + range * (N) / 2, centre + range * (N + 1.0) / 2);
					Candidate[m] = ux(e1);
				}
			}
			double distanceNearestSTC = 0.0, distanceNearestFTC = 0.0, distanceCurrent = 0.0, distanceNearest = 0.0;
			//select the test case with farthest distance
			for (int p = 0; p < candidateNum; p++) {
				distanceNearestSTC = 0; distanceNearestFTC = 0;
				//calculate the distance
				for (int q = 0; q < Executedx.size(); q++) {
					distanceNearestSTC += fabs(Executedx[q] - Candidate[p]);
					distanceNearestFTC += fabs(Executedx[q] - (3 * Candidate[p]));
				}	
				//compute the average
				distanceNearestSTC = (distanceNearestSTC + distanceNearestFTC) / 2 / Executedx.size();
				//compare the distances
				if (p == 0) {
					distanceNearest = distanceNearestSTC;
					farthestCandidate.clear();
					farthestCandidate.push_back(p);
				}
				else {
					if (distanceNearest < distanceNearestSTC) {
						distanceNearest = distanceNearestSTC;
						farthestCandidate.clear();
						farthestCandidate.push_back(p);
					}
					else if (fabs(distanceNearest - distanceNearestSTC) < EPS)
						farthestCandidate.push_back(p);
				}
			}
			//If more than one candidate MG is selected, consider the distances between source and follow-up candidates
			int farthestCandidateNum = 0; double distanceFurthestAll = 0.0;
			if (farthestCandidate.size() > 1) {
				distanceFurthestAll = 0.0;
				for (int p = 0; p < farthestCandidate.size(); p++) {
					distanceCurrent = fabs(2 * Candidate[farthestCandidate[p]]);
					if (distanceFurthestAll < distanceCurrent) {
						distanceFurthestAll = distanceCurrent;
						farthestCandidateNum = p;
					}
				}
			}
			//select one source candidate for execution
			x1 = Candidate[farthestCandidate[farthestCandidateNum]];
		}
		//generate FTCs
		x2 = 3 * x1;
		//execution
		if (mutantNum == 1) {
			ret1 = tanh1Mutant1(x1);
			ret2 = tanh1Mutant1(x2);
		}
		else {
			ret1 = tanh1Mutant2(x1);
			ret2 = tanh1Mutant2(x2);
		}
		if (fabs((ret1 * ret1 * ret1 + 3 * ret1) / (1 + 3 * ret1 * ret1) - ret2) > EPS) {
			//output1 << "Violate: " << x1 << " " << x2 << endl; //output the STC and FTC that violate the MR
			return n;
		}
		//output1 << x1 << " " << x2 << endl; //output the generated STCs and FTCs
		//add the executed test cases into the corresponding test sets
		//includes all executed STCs and FTCs
		Executedx.push_back(x1);
		Executedx.push_back(x2);
		//includes the executed STCs and FTCs that are inside the source input domain
		ExecutedxSource.push_back(x1);
		if (3 * x1 < Top) {
			ExecutedxSource.push_back(x2);
		}
	}
	return 100000;
}

int SFIDMTART(int randomNumber, int mutantNum) {
	//random generator seed
	default_random_engine e1(randomNumber);
	unsigned seed = randomNumber;  
	srand(seed);
	initialize();
	for (int n = 0; n < 100000; n++) {
		CandidateSubdomainList.clear();
		if (n == 0) { //randomly generate the first test case
			x1 = uRandom(e1);
		}
		else {
			//calculate the size of subdomains
			range = (Top - Down) / (ExecutedxSource.size() + 1.0);
			//generate subdomains
			CandidateSubdomain.clear();
			CandidateSubdomain.assign(ExecutedxSource.size() + 1, 0);
			//map all executed STCs into the subdomains
			for (int i = 0; i < ExecutedxSource.size(); i++) {
				Position = fabs(ExecutedxSource[i] - centre);
				CandidateSubdomain[floor(Position / range * 2)] = 1;
			}
			for (int i = 0; i < ExecutedxSource.size() + 1; i++)
				if (CandidateSubdomain[i] == 0)
					CandidateSubdomainList.push_back(i);
			//generate candidates in the target (empty) subdomain
			for (int m = 0; m < candidateNum; m++) {
				N = CandidateSubdomainList[rand() % CandidateSubdomainList.size()];
				Random = rand() % 2;
				//centre
				if (N == 0) {
					uniform_real_distribution<double> ux(centre - range * (N + 1.0) / 2, centre + range * (N + 1.0) / 2);
					Candidate[m] = ux(e1);
				}
				//left
				else if (Random == 0) {
					uniform_real_distribution<double> ux(centre - range * (N + 1.0) / 2, centre - range * (N) / 2);
					Candidate[m] = ux(e1);
				}
				//right
				else if (Random == 1) {
					uniform_real_distribution<double> ux(centre + range * (N) / 2, centre + range * (N + 1.0) / 2);
					Candidate[m] = ux(e1);
				}
			}
			double distanceNearest = 0.0, distanceCurrent = 0.0, distanceFurthest = 0.0;
			//select the test case with farthest distance
			for (int p = 0; p < candidateNum; p++) {
				//compute the distances
				for (int q = 0; q < ExecutedxSource.size(); q++) {
					distanceCurrent = fabs(ExecutedxSource[q] - Candidate[p]) + fabs(ExecutedxFollow[q] - 3 * Candidate[p]);
					if (q == 0)
						distanceNearest = distanceCurrent;
					else if (distanceNearest > distanceCurrent)
						distanceNearest = distanceCurrent;
				}
				//compare the distances
				if (p == 0) {
					distanceFurthest = distanceNearest;
					farthestCandidate.clear();
					farthestCandidate.push_back(p);
				}
				else {
					if (distanceFurthest < distanceNearest) {
						distanceFurthest = distanceNearest;
						farthestCandidate.clear();
						farthestCandidate.push_back(p);
					}
					else if (fabs(distanceFurthest - distanceNearest) < EPS)
						farthestCandidate.push_back(p);
				}
			}
			//If more than one candidate MG is selected, randomly select one of them
			int farthestCandidateNum = 0;
			if (farthestCandidate.size() > 1) {
				farthestCandidateNum = rand() % farthestCandidate.size();
			}
			//select one source candidate for execution
			x1 = Candidate[farthestCandidate[farthestCandidateNum]];
		}
		//generate FTCs
		x2 = 3 * x1;
		//execution
		if (mutantNum == 1) {
			ret1 = tanh1Mutant1(x1);
			ret2 = tanh1Mutant1(x2);
		}
		else {
			ret1 = tanh1Mutant2(x1);
			ret2 = tanh1Mutant2(x2);
		}
		if (fabs((ret1 * ret1 * ret1 + 3 * ret1) / (1 + 3 * ret1 * ret1) - ret2) > EPS) {
			//output1 << "Violate: " << x1 << " " << x2 << endl; //output the STC and FTC that violate the MR
			return n;
		}
		//output1 << x1 << " " << x2 << endl; //output the generated STCs and FTCs
		//add the executed STCs and FTCs into the corresponding test sets
		ExecutedxSource.push_back(x1); ExecutedxFollow.push_back(x2);
	}
	return 100000;
}

int main() {
	//output.open("tanh1RTTestCase.txt"); //outputs test cases. Format: STC FTC
	int TempRT1 = 0, TempMax1 = 0, TempAve1 = 0, TempMin1 = 0, TempSFID1 = 0; //store the F-measure values
	int TempRT2 = 0, TempMax2 = 0, TempAve2 = 0, TempMin2 = 0, TempSFID2 = 0; //store the F-measure values
	for (int m = 0; m < 10000; m++) {
		TempRT1 += RT(m, 1); //execute against mutant 1
		TempMax1 += MTARTMax(m, 1); //execute against mutant 1
		TempAve1 += MTARTAve(m, 1); //execute against mutant 1
		TempMin1 += MTARTMin(m, 1); //execute against mutant 1
		TempSFID1 += SFIDMTART(m, 1); //execute against mutant 1
		
		TempRT2 += RT(m, 2); //execute against mutant 2
		TempMax2 += MTARTMax(m, 2); //execute against mutant 2
		TempAve2 += MTARTAve(m, 2); //execute against mutant 2
		TempMin2 += MTARTMin(m, 2); //execute against mutant 2
		TempSFID2 += SFIDMTART(m, 2); //execute against mutant 2
	}
	cout << "Mean F-measure of RT (mutant 1): " << TempRT1 / 10000 << endl; //outputs F-RT values
	cout << "Mean F-measure of MT-ART-Max (mutant 1): " << TempMax1 / 10000 << endl; //outputs F-MTARTMax values
	cout << "Mean F-measure of MT-ART-Ave (mutant 1): " << TempAve1 / 10000 << endl; //outputs F-MTARTAve values
	cout << "Mean F-measure of MT-ART-Min (mutant 1): " << TempMin1 / 10000 << endl; //outputs F-MTARTMin values
	cout << "Mean F-measure of SFIDMT-ART (mutant 1): " << TempSFID1 / 10000 << endl; //outputs F-SFIDMTART values
	
	cout << "Mean F-measure of RT (mutant 2): " << TempRT2 / 10000 << endl; //outputs F-RT values
	cout << "Mean F-measure of MT-ART-Max (mutant 2): " << TempMax2 / 10000 << endl; //outputs F-MTARTMax values
	cout << "Mean F-measure of MT-ART-Ave (mutant 2): " << TempAve2 / 10000 << endl; //outputs F-MTARTAve values
	cout << "Mean F-measure of MT-ART-Min (mutant 2): " << TempMin2 / 10000 << endl; //outputs F-MTARTMin values
	cout << "Mean F-measure of SFIDMT-ART (mutant 2): " << TempSFID2 / 10000 << endl; //outputs F-SFIDMTART values
}
