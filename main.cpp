/*
Bilateral Exchange Model
Originally written by Rob Axtell
Extended by Stefan McCabe

12/11/2015
*/

//  Rob:
//	Additional work to do on this code:
//
//	A. Multithread the main 'Equilibrate()' method
//	B. Rationalize use of Init()s, Reset()s and constructors
//	C. Use my 'timer()' object instead of what is here...
//	D. Generalize to CES preferences?
//	E. Shrink below 1000 lines?

#include <cmath>
#include <iostream>

const double Version = 0.97;

const int NumberOfAgents = 10000000; //	10^7
const int NumberOfCommodities = 2;	 // 	2 x 10^4

//const int PairwiseInteractionsPerPeriod = NumberOfAgents / 2;
const int PairwiseInteractionsPerPeriod = 50;

const double alphaMin = 0.01;
const double alphaMax = 0.99;
const int wealthMin = 50;
const int wealthMax = 150;

const bool DefaultSerialExecution = false;	// FALSE means parallel
const int AgentsToRandomize = NumberOfAgents / 10;	// This is done only for parallel interactions	
const int RandomizationMethod = 0;	// 0 means swap pairs
									// 1 means "threaded"
const int RequestedEquilibrations = 1;
const bool SameAgentInitialCondition = false;

const double trade_eps = 0.01;
const double exp_trade_eps = exp(trade_eps);

const int termination_criterion = -1; 	// -2 means termination based on time
										// -1 means use L2 norm on MRSs
										//	0 means use L∞ norm on MRSs
										//	1 means use relative increase in V
										//	2 means use absolute increase in V
const long TerminationTime = 20000;
const double termination_eps = 0.01;	//	10^-1
const long CheckTerminationThreshhold = 100000;	//	1.5 x 10^9 periods ~ 60 x 10^9 interactions
const int CheckTerminationPeriod = 10000;	//	Typical values:
											//	A = 100, 		N = 10,		use 1	
											//	A = 100, 		N = 50,		use 25
											//	A = 100, 		N = 100,	use 100
											//	A = 100, 		N = 500,	use 2500
											//	A = 100, 		N = 1000,	use 10,000
											//	A = 100, 		N = 5000,	use 250,000
											//	A = 100, 		N = 10,000,	use 1,000,000
											//	A = 1000, 		N = 2, 		use 1												
											//	A = 10,000, 	N = 2, 		use 10											
											//	A = 100,000, 	N = 2,		use 100											
											//	A = 1,000,000, 	N = 2, 		use 1000									
											//	A = 10,000,000, N = 2, 		use 10,000

const bool ShockPreferences = false;
const int ShockPeriod = 10;
const double MinShock = 1.0;
const double MaxShock = 5.0;

const bool On = true;
const bool Off = false;

const bool debug = Off;

const bool UseRandomSeed = false;
const int NonRandomSeed = 1;

const bool PrintEndowments = Off;
const bool PrintIntermediateOutput = On;
const int IntermediateOutputPrintPeriod = 10000; //	20x10^6
const bool PrintConvergenceStats = true;
const bool PrintFinalCommodityList = Off;

typedef double	CommodityArray[NumberOfCommodities+1];	//	Size: NumberOfCommodities * 8 bytes
typedef CommodityArray *CommodityArrayPtr;

double inline	Dot (CommodityArrayPtr vector1, CommodityArrayPtr vector2);

void	WriteActualTimeAndDate (char DescriptionString[63]);

class	MemoryObject
{
	long	start;
public:
	MemoryObject();
	void WriteMemoryRequirements();
} MemoryState;

class RandomNumberGenerator
{
	int last;
public:
	RandomNumberGenerator();
	int LongInteger();								//	Return INTEGER-valued random numbers in the interval (0, MaxLongInt)
	int IntegerInRange (int min, int max);			//	Return INTEGER-valued random numbers in the interval [min, max]
	double UnitReal();								//	Return REAL-valued random numbers in the interval [0, 1]
	double RealInRange (double min, double max);	//	Return REAL-valued random numbers in the interval [min, max]
} RNG;

class Data
{					//	Size:
	int N;			//		4 bytes
	double min;		//		8 bytes
	double max;		//		8 bytes
	double sum;		//		8 bytes
	double sum2;	//		8 bytes
public:
	Data();
	void Init();
	void AddDatuum (double Datuum);
	int GetN() {return N;};
	double GetMin() {return min;};
	double GetMax() {return max;};
	double GetDelta() {return max - min;};
	double GetAverage() {if (N > 0) return sum/N; else return 0.0;};
	double GetExpAverage() {return exp(GetAverage());};
	double GetVariance();
	double GetStdDev() {return sqrt(GetVariance());};
};	//	Total:  36 bytes

typedef Data *DataPtr;

class CommodityData
{														//	Size:																								
	Data data[NumberOfCommodities+1];					//		NumberOfCommodities * 8 bytes
public:
	CommodityData();
	void Clear();
	DataPtr GetData(int index) {return &data[index];};
	double L2StdDev();
	double LinfStdDev();																		
};														//	Total:	NumberOfCommodities * 8 bytes									
														//	Example:	NumberOfCommodities = 100, size = 800 bytes

class Agent
{								//	Size:
	CommodityArray alphas;		//		NumberOfCommodities * 8 bytes
	CommodityArray endowment;	//		NumberOfCommodities * 8 bytes							
	CommodityArray initialMRSs;	//		NumberOfCommodities * 8 bytes
	CommodityArray allocation;	//		NumberOfCommodities * 8 bytes
	CommodityArray currentMRSs;	//		NumberOfCommodities * 8 bytes
	double initialUtility;		//		8 bytes
	double initialWealth;		//		8 bytes
public:
	Agent();
	void Init();
	void Reset();
	double GetAlpha (int CommodityIndex) {return alphas[CommodityIndex];};
	void SetAlpha (int CommodityIndex, double alpha) {alphas[CommodityIndex] = alpha;};
	double GetEndowment (int CommodityIndex) {return endowment[CommodityIndex];};
	double MRS (int CommodityIndex, int Numeraire);
	void ComputeMRSs();
	double GetInitialMRS(int CommodityIndex) {return initialMRSs[CommodityIndex];};
	double GetAllocation(int CommodityIndex) {return allocation[CommodityIndex];};
	void IncreaseAllocation(int CommodityIndex, double amount) {allocation[CommodityIndex] += amount;};
	double GetCurrentMRS(int CommodityIndex) {return currentMRSs[CommodityIndex];};
	CommodityArrayPtr GetCurrentMRSs() {return &currentMRSs;};
	double Utility();
	double GetInitialUtility() {return initialUtility;};
	double Wealth (CommodityArrayPtr prices) {return Dot(&allocation, prices);};
	double GetInitialWealth() {return initialWealth;};
};
								//	Total:	20 + NumberOfCommodities * 40 bytes
								//	Example:	NumberOfCommodities = 2, size = 20 + 80 = 100
								//	Example:	NumberOfCommodities = 100, size = 20 + 4000 = 4020
typedef Agent *AgentPtr;

class AgentPopulation
{													//	Size:
	AgentPtr Agents[NumberOfAgents + 1];			//	NumberOfAgents * 4 bytes
	CommodityArray Volume;	
	CommodityData AlphaData, EndowmentData, LnMRSsData;
	bool LnMRSsDataUpToDate;
	void ComputeLnMRSsDistribution();
	double LastSumOfUtilities;
	double ComputeSumOfUtilities();
	double ComputeIncreaseInSumOfUtilities() {return ComputeSumOfUtilities() - LastSumOfUtilities;};
	double ComputeRelativeIncreaseInSumOfUtilities() {return ComputeIncreaseInSumOfUtilities() / LastSumOfUtilities;};
	void GetSerialAgentPair (AgentPtr& Agent1, AgentPtr& Agent2);
	int ActiveAgentIndex;
	void RandomizeAgents (int NumberToRandomize);
	void RandomizeAgents2 (int NumberToRandomize);
	void GetParallelAgentPair (AgentPtr& Agent1, AgentPtr& Agent2);
	void (AgentPopulation::*GetAgentPair) (AgentPtr& Agent1, AgentPtr& Agent2);
public:
	AgentPopulation();								//	Constructor...
	void Init();
	void Reset();
	long Equilibrate(int NumberOfEquilibrationsSoFar);
	void ConvergenceStatistics(CommodityArray VolumeStats);
	void CompareTwoAgents(AgentPtr Agent1, AgentPtr Agent2);
	void ShockAgentPreferences();
};
													//	Example:	NumberOfAgents = 100, size = 400 bytes

typedef AgentPopulation *AgentPopulationPtr;

double inline	Dot (CommodityArrayPtr vector1, CommodityArrayPtr vector2)
{
	double sum = 0.0;
	for (int i = 1; i <= NumberOfCommodities; ++i)
		sum += (*vector1)[i] * (*vector2)[i];
	return sum;
}	//	Dot

void	WriteActualTimeAndDate (char DescriptionString[63])
{
	time_t	SysTime;
	tm	*DateAndTimePtr;
	
	SysTime = time(NULL);
	DateAndTimePtr = localtime(&SysTime);
	std::cout << DescriptionString;
	if (DateAndTimePtr->tm_hour < 10) 
		std::cout << "0" << DateAndTimePtr->tm_hour;
	else
		std::cout << DateAndTimePtr->tm_hour;
	std::cout << ":";
	if (DateAndTimePtr->tm_min < 10)
		std::cout << "0" << DateAndTimePtr->tm_min;
	else
		std::cout << DateAndTimePtr->tm_min;
	std::cout << ":";
	if (DateAndTimePtr->tm_sec < 10)
		std::cout << "0" << DateAndTimePtr->tm_sec;
	else
		std::cout << DateAndTimePtr->tm_sec;
	std::cout << " on ";
	switch (DateAndTimePtr->tm_wday + 1)
	{
		case 1:
			std::cout << "Sunday ";
			break;
		case 2:
			std::cout << "Monday ";
			break;
		case 3:
			std::cout << "Tuesday ";
			break;
		case 4:
			std::cout << "Wednesday ";
			break;
		case 5:
			std::cout << "Thursday ";
			break;
		case 6:
			std::cout << "Friday ";
			break;
		case 7:
			std::cout << "Saturday ";
			break;
	}
	std::cout << DateAndTimePtr->tm_mon + 1 << "/" << DateAndTimePtr->tm_mday << "/";
	std::cout << DateAndTimePtr->tm_year + 1900 << std::endl;
}	//	WriteActualTimeAndDate

/*===========
	==Methods==
	===========*/

MemoryObject::MemoryObject():
	start(0)
{}	//	MemoryObject::MemoryObject()

void	MemoryObject::WriteMemoryRequirements()
{
	std::cout << "Memory requirements:" << std::endl;
	std::cout << "Size of Agent object: " << sizeof(Agent) << " bytes" << std::endl;
	std::cout << "Size of Data object: " << sizeof(Data) << " bytes" << std::endl;
	std::cout << "Size of CommodityData object: " << sizeof(CommodityData) << " bytes" << std::endl;
	std::cout << "Size of Population object: " << sizeof(AgentPopulation) + sizeof(Agent) * NumberOfAgents << " bytes" << std::endl;
	
	long TotalBytesRequired = (sizeof(RandomNumberGenerator) + sizeof(MemoryObject) + sizeof(CommodityData) + sizeof(AgentPopulation) + sizeof(Agent) * NumberOfAgents);
	
	if (TotalBytesRequired < 1000000)
		std::cout << "Total memory requirements: " << double(TotalBytesRequired/1000.0) << " kilobytes" << std::endl;
	else
		std::cout << "Total memory requirements: " << double(TotalBytesRequired/1000000.0) << " megabytes" << std::endl;
	std::cout << std::endl;
}	//	MemoryObject::WriteMemoryRequirements()

RandomNumberGenerator::RandomNumberGenerator():
  last(0)
{
	if (UseRandomSeed)
		last = (int)time(NULL);
	else
		last = NonRandomSeed;
}	//	RandomNumberGenerator::Init()

int	RandomNumberGenerator::LongInteger()
{/*
 This method generates INTEGER-valued random numbers in the interval [1, INT_MAX - 1]}
 Source: Bratley, Fox and Schrage, 1987
 */
	int k = last / 127773;
	last = 16807 * (last - k	*	127773) - k	*	2836;
	if (last < 0)
		last += INT_MAX;
	return last;
}	//	RandomNumberGenerator::LongInteger()

int	RandomNumberGenerator::IntegerInRange (int min, int max)
{/*
 This method generates INTEGER-valued random numbers in the interval [min, max]
 */
	return min + LongInteger() % (max - min + 1);
}	//	RandomNumberGenerator::IntegerInRange()

double	RandomNumberGenerator::UnitReal()
{/*
 This method generates REAL-valued random numbers in the interval [0, 1]
 */
	return double(LongInteger()) / double(INT_MAX);
}	//	RandomNumberGenerator::UnitReal()

double	RandomNumberGenerator::RealInRange (double min, double max)
{/*
 This method generates REAL-valued random numbers in the interval [min, max]}
 */
	return min + (max - min) * UnitReal();
}	//	RandomNumberGenerator::RealInRange()

Data::Data():
  N(0), min(1000000.0), max (0.0), sum (0.0), sum2(0.0)
{}	//	Data::Data()

void Data::Init()
{
	N = 0;
	min = 1000000.0;	//	10^6
	max = 0.0;
	sum = 0.0;
	sum2 = 0.0;
}	//	Data::Data()

void	Data::AddDatuum (double Datuum)
{
	N = N + 1;
	if (Datuum < min)
		min = Datuum;
	if (Datuum > max)
		max = Datuum;
	sum += Datuum;
	sum2 += Datuum * Datuum;
}	//	Data::AddDatuum()

double	Data::GetVariance()
{
	double avg, arg;
	
	if (N > 1)
	{
		avg = GetAverage();
		arg = sum2 - N * avg * avg;
		return arg / (N - 1);
	}
	else
		return 0.0;
}	//	Data::GetVariance()

CommodityData::CommodityData()
{}

void	CommodityData::Clear()
{
	//	Initialize data objects
	//
	for (int i = 1; i <= NumberOfCommodities; ++i)
		data[i].Init();
}	//	CommodityData::Clear()

double	CommodityData::L2StdDev()
{
	double sum = 0.0;
	
	// Instead of computing the standard deviation for each commodity (and calling sqrt() N times), find the largest variance and from it get the std dev
	//
	for (int i = 1; i <= NumberOfCommodities; ++i)
		sum += data[i].GetVariance();
	return sqrt(sum);
}	//	CommodityData::L2StdDev

double	CommodityData::LinfStdDev()
{
	double var, max = 0.0;
	
	// Instead of computing the standard deviation for each commodity (and  calling sqrt() N times), find the largest variance and from it get the std dev
	//
	for (int i = 1; i <= NumberOfCommodities; ++i)
	{
		var = data[i].GetVariance();
		if (var > max)
			max = var;
	}
	return sqrt(max);
}	//	CommodityData::LinfStdDev

Agent::Agent():
  initialUtility(0.0), initialWealth(0.0)
{
	Init();
}	

void	Agent::Init()
{
	int CommodityIndex;
	
	//	First generate and normalize the exponents...
	//
	double sum = 0.0;
	for (CommodityIndex = 1; CommodityIndex <= NumberOfCommodities; ++CommodityIndex)
	{
		alphas[CommodityIndex] = RNG.RealInRange(alphaMin, alphaMax);
		sum += alphas[CommodityIndex];
	}
	//	Next, fill up the rest of the agent fields...
	//
	for (CommodityIndex = 1; CommodityIndex <= NumberOfCommodities; ++CommodityIndex)
	{
		alphas[CommodityIndex] = alphas[CommodityIndex] / sum;
		endowment[CommodityIndex] = RNG.RealInRange(wealthMin, wealthMax);
		allocation[CommodityIndex] = endowment[CommodityIndex];
		initialMRSs[CommodityIndex] = MRS(CommodityIndex, 1);
		currentMRSs[CommodityIndex] = initialMRSs[CommodityIndex];
	}
	initialUtility = Utility();
	initialWealth = Wealth(&initialMRSs);
}	//	Agent:Init()

void	Agent::Reset()
{
	for (int CommodityIndex = 1; CommodityIndex <= NumberOfCommodities; ++CommodityIndex)
	{
		allocation[CommodityIndex] =  endowment[CommodityIndex];
		currentMRSs[CommodityIndex] = initialMRSs[CommodityIndex];
	}
}	//	Agent::Reset

double	Agent::MRS (int CommodityIndex, int Numeraire)
{
	return (alphas[CommodityIndex] * allocation[Numeraire]) / (alphas[Numeraire] * allocation[CommodityIndex]);
}	//	Agent::MRS()

void	Agent::ComputeMRSs()
{
	int i;
	currentMRSs[1] = 1.0;
	for (i = 2; i <= NumberOfCommodities; ++i)
		currentMRSs[i] = MRS(i, 1);
}	//		Agent::ComputeMRSs()

double	Agent::Utility()
{
	double product = 1.0;
	
	for (int i = 1; i <= NumberOfCommodities; ++i)
		product *= pow(allocation[i], alphas[i]);
	return product;
}	//	Agent::Utility()

void	AgentPopulation::ComputeLnMRSsDistribution()
{
	LnMRSsData.Clear();
	
	for (int i = 1; i <= NumberOfAgents; ++i)
	{
		Agents[i]->ComputeMRSs();
		for (int j = 1; j <= NumberOfCommodities; ++j)
			LnMRSsData.GetData(j)->AddDatuum(log(Agents[i]->GetCurrentMRS(j)));		
	}		//	for i...
	LnMRSsDataUpToDate = true;
}	//	AgentPopulation::ComputeLnMRSsDistribution()

double	AgentPopulation::ComputeSumOfUtilities()
{
	double sum = 0.0;
	
	for (int i = 1; i <= NumberOfAgents; ++i)
		sum += Agents[i]->Utility();
	return sum;
}	//	AgentPopulation::ComputeSumOfUtilities()

void	inline	AgentPopulation::GetSerialAgentPair (AgentPtr& Agent1, AgentPtr& Agent2)
{
	Agent1 = Agents[RNG.IntegerInRange(1, NumberOfAgents)];
	do
		Agent2 = Agents[RNG.IntegerInRange(1, NumberOfAgents)];
	while (Agent2 == Agent1);
}	//	AgentPopulation::GetSerialAgentPair()

void	inline AgentPopulation::RandomizeAgents (int NumberToRandomize)
{
	// Note:  If NumberToRandomize < 2 then it is increased to 2, so that 1 pair is swapped each period
	//
	int AgentIndex1, AgentIndex2;
	int PairsToRandomize; 
	AgentPtr TempPtr;
	
	PairsToRandomize = NumberToRandomize / 2;
	if (PairsToRandomize < 1)
		PairsToRandomize = 1;
	
	for (int i = 1; i <= PairsToRandomize; ++i)
	{
		AgentIndex1 = RNG.IntegerInRange(1, NumberOfAgents);
		do
			AgentIndex2 = RNG.IntegerInRange(1, NumberOfAgents);
		while (AgentIndex1 == AgentIndex2);
		TempPtr = Agents[AgentIndex1];
		Agents[AgentIndex1] = Agents[AgentIndex2];
		Agents[AgentIndex2] = TempPtr;
	}
}	//	AgentPopulation::RandomizeAgents()

void	inline AgentPopulation::RandomizeAgents2 (int NumberToRandomize)
{
	// Note:  The NumberToRandomize need not be even, just > 1
	//
	int AgentIndex1, AgentIndex2;
	AgentPtr TempPtr;
	
	AgentIndex1 = RNG.IntegerInRange(1, NumberOfAgents);
	TempPtr = Agents[AgentIndex1];
	
	if (NumberToRandomize == 1)
		NumberToRandomize = 2;
	for (int i = 1; i < NumberToRandomize; ++i)
	{
		AgentIndex2 = RNG.IntegerInRange(1, NumberOfAgents);
		Agents[AgentIndex1] = Agents[AgentIndex2];
		AgentIndex1 = AgentIndex2;
	}
	Agents[AgentIndex2] = TempPtr;
	
}	//	AgentPopulation::RandomizeAgents2()

void	inline AgentPopulation::GetParallelAgentPair (AgentPtr& Agent1, AgentPtr& Agent2)
{
	if (ActiveAgentIndex > NumberOfAgents)
	{
		ActiveAgentIndex = 1;
		RandomizeAgents(AgentsToRandomize);
	}
	Agent1 = Agents[ActiveAgentIndex];
	Agent2 = Agents[++ActiveAgentIndex];
	ActiveAgentIndex++;
}	//	AgentPopulation::GetParallelAgentPair()

AgentPopulation::AgentPopulation():
	AlphaData(), EndowmentData(), LnMRSsData(), LnMRSsDataUpToDate(true), LastSumOfUtilities(0.0), ActiveAgentIndex(0), GetAgentPair(NULL)
{
	for (int i = 1; i <= NumberOfAgents; i++)
		Agents[i] = new Agent;

	if (DefaultSerialExecution)
		GetAgentPair = &AgentPopulation::GetSerialAgentPair;
	else
		GetAgentPair = &AgentPopulation::GetParallelAgentPair;
}	// Constructor...

void	AgentPopulation::Init()
{
	AgentPtr ActiveAgent;
	int CommodityIndex;
	Data	InitialOwnWealthData;
	
	ActiveAgentIndex = 1;
	AlphaData.Clear();
	EndowmentData.Clear();
	LnMRSsData.Clear();
	
	//	Cycle through the agents...
	//
	for (int AgentIndex = 1; AgentIndex <= NumberOfAgents; ++AgentIndex)
	{
		ActiveAgent = Agents[AgentIndex];
		ActiveAgent->Init();
		
		//	Next, fill up the rest of the agent fields...
		//
		for (CommodityIndex = 1; CommodityIndex <= NumberOfCommodities; ++CommodityIndex)
		{
			AlphaData.GetData(CommodityIndex)->AddDatuum(ActiveAgent->GetAlpha(CommodityIndex));
			EndowmentData.GetData(CommodityIndex)->AddDatuum(ActiveAgent->GetEndowment(CommodityIndex));
			LnMRSsData.GetData(CommodityIndex)->AddDatuum(log(ActiveAgent->GetInitialMRS(CommodityIndex)));
		}	//	for (CommodityIndex...
		InitialOwnWealthData.AddDatuum(ActiveAgent->GetInitialWealth());
		
	}	//	for (AgentIndex...
	
	LnMRSsDataUpToDate = true;
	LastSumOfUtilities = ComputeSumOfUtilities();
	
	//	Finally, display stats on the instantiated population...
	//
	std::cout << "******************************" << std::endl;
	std::cout << std::endl;			
	if (PrintEndowments)
	{
		std::cout << "Initial endowments:" << std::endl;
		for (CommodityIndex = 1; CommodityIndex <= NumberOfCommodities; ++CommodityIndex)
		{
			std::cout << "Commodity " << CommodityIndex << ": <exp.> = " << AlphaData.GetData(CommodityIndex)->GetAverage() << "; s.d. = " << AlphaData.GetData(CommodityIndex)->GetStdDev();
			std::cout << "; <endow.> = " << EndowmentData.GetData(CommodityIndex)->GetAverage() << "; s.d. = " << EndowmentData.GetData(CommodityIndex)->GetStdDev();
			std::cout << "; <MRS> = " << LnMRSsData.GetData(CommodityIndex)->GetExpAverage() << "; s.d. = " << LnMRSsData.GetData(CommodityIndex)->GetStdDev() << std::endl;
		}
		std::cout << std::endl;
	}
	std::cout << "Average initial wealth (@ own prices) = " << InitialOwnWealthData.GetAverage() << "; standard deviation = " << InitialOwnWealthData.GetStdDev() << std::endl;
	std::cout << "Initial sum of utilities = " << LastSumOfUtilities << std::endl;;
	std::cout << std::endl;
}	//	AgentPopulation::Init()

void	AgentPopulation::Reset()
{
	for (int AgentIndex = 1; AgentIndex <= NumberOfAgents; ++AgentIndex)
		Agents[AgentIndex]->Reset();
}	//	AgentPopulation::Reset

long	AgentPopulation::Equilibrate (int NumberOfEquilibrationsSoFar)
{
	char OutputStr[63];
	bool Converged = false;
	long theTime = 0;
	long TotalInteractions = 0;
	AgentPtr Agent1, Agent2, LargerMRSAgent, SmallerMRSAgent;
	int Commodity1, Commodity2;
	double MRSratio12, MRSratio;
	register double LAgentalpha1, LAgentalpha2, LAgentx1, LAgentx2, SAgentalpha1, SAgentalpha2, SAgentx1, SAgentx2;
	register double num, denom, delta1, delta2;
	double Agent1PreTradeUtility, Agent2PreTradeUtility;
	
	std::cout << "******************************" << std::endl;
	std::cout << std::endl;
	sprintf(OutputStr, "Equilibration #%d starting at time: ", NumberOfEquilibrationsSoFar);
	WriteActualTimeAndDate(OutputStr);
	std::cout << std::endl;
	
	//	Next, initialize some variables...
	//
	for (int i = 1; i <= NumberOfCommodities; ++i)
		Volume[i] = 0.0;
	
	//	Start up the exchange process here...
	//
	do {
		LnMRSsDataUpToDate = false;	//	Since these data are gonna change...
		
		++theTime;
		if ((ShockPreferences) && (theTime % ShockPeriod == 0))
			ShockAgentPreferences();
		for (int i = 1; i <= PairwiseInteractionsPerPeriod; ++i)
		{
			//	First, select the agents to be active...
			//
			(this->*GetAgentPair) (Agent1, Agent2);
			
			if (debug)
			{
				Agent1PreTradeUtility = Agent1->Utility();
				Agent2PreTradeUtility = Agent2->Utility();
			}
			//	Next, select the commodities to trade...
			//
			if (NumberOfCommodities == 2)
			{
				Commodity1 = 1;
				Commodity2 = 2;
			}
			else
			{
				Commodity1 = RNG.IntegerInRange(1, NumberOfCommodities);
				do
					Commodity2 = RNG.IntegerInRange(1, NumberOfCommodities);
				while (Commodity2 == Commodity1);
			}
			//	Compare MRSs...
			//
			MRSratio12 = Agent1->MRS(Commodity2, Commodity1) / Agent2->MRS(Commodity2, Commodity1);
			if (MRSratio12 > 1.0)
				MRSratio = MRSratio12;
			else
				MRSratio = 1.0/MRSratio12;
			if (MRSratio >= exp_trade_eps)	//	do exchange...
			{
				if (MRSratio12 > 1.0)
				{
					LargerMRSAgent = Agent1;
					SmallerMRSAgent = Agent2;
				}
				else
				{
					LargerMRSAgent = Agent2;
					SmallerMRSAgent = Agent1;
				}
				//	Here are the guts of bilateral Walrasian exchange
				//
				SAgentalpha1 = SmallerMRSAgent->GetAlpha(Commodity1);
				SAgentalpha2 = SmallerMRSAgent->GetAlpha(Commodity2);
				SAgentx1 = SmallerMRSAgent->GetAllocation(Commodity1);
				SAgentx2 = SmallerMRSAgent->GetAllocation(Commodity2);
				LAgentalpha1 = LargerMRSAgent->GetAlpha(Commodity1);
				LAgentalpha2 = LargerMRSAgent->GetAlpha(Commodity2);
				LAgentx1 = LargerMRSAgent->GetAllocation(Commodity1);
				LAgentx2 = LargerMRSAgent->GetAllocation(Commodity2);
				
				num = (SAgentalpha1 * LAgentalpha2 * LAgentx1 * SAgentx2) - (LAgentalpha1 * SAgentalpha2 * LAgentx2 * SAgentx1);
				denom = LAgentalpha1 * LAgentx2 + SAgentalpha1 * SAgentx2;
				delta1 = num / denom;
				denom = LAgentalpha2 * LAgentx1 + SAgentalpha2 * SAgentx1;
				delta2 = num / denom;
				
				SmallerMRSAgent->IncreaseAllocation(Commodity1, delta1);
				SmallerMRSAgent->IncreaseAllocation(Commodity2, -delta2);
				LargerMRSAgent->IncreaseAllocation(Commodity1, -delta1);
				LargerMRSAgent->IncreaseAllocation(Commodity2, delta2);
				
				++TotalInteractions;
				Volume[Commodity1] += delta1;
				Volume[Commodity2] += delta2;
			}	//	if (MRSratio...
			
			if (debug)
			{
				if (Agent1->Utility() < Agent1PreTradeUtility)
					std::cout << "!!!Utility decreasing trade by agent #1!!!  Actual utility change = " << Agent1->Utility() - Agent1PreTradeUtility << std::endl;
				if (Agent2->Utility() < Agent2PreTradeUtility)
					std::cout << "!!!Utility decreasing trade by agent #2!!!  Actual utility change = " << Agent2->Utility() - Agent2PreTradeUtility << std::endl;
			}
		}	//	for i...
		
		//	Check for termination...
		//
		if ((theTime > CheckTerminationThreshhold) && (theTime % CheckTerminationPeriod == 0))
			switch (termination_criterion)
		{
			case -2:
				if (theTime >= TerminationTime)
					Converged = true;
				break;
			case -1:		// 	Termination based on L2 norm of MRS distribution
				ComputeLnMRSsDistribution();
				if (LnMRSsData.L2StdDev() < termination_eps)
					Converged = true;
				break;
			case 0:		// 	Termination based on L∞ norm of MRS distribution
				ComputeLnMRSsDistribution();
				if (LnMRSsData.LinfStdDev() < termination_eps)
					Converged = true;
				break;
			case 1:		//	Termination based on relative increase in V
				if (ComputeRelativeIncreaseInSumOfUtilities() < termination_eps)
					Converged = true;
				break;
			case 2:		//	Termination based on absolute increase in V
				if (ComputeIncreaseInSumOfUtilities() < termination_eps)
					Converged = true;
				break;
		}	//	switch...
		
		//	Display stats if the time is right...
		//
		if (PrintIntermediateOutput)
			if (theTime % IntermediateOutputPrintPeriod == 0)
			{
				//printf("Through time %ld, %ld total exchanges; ", theTime, TotalInteractions);
				std::cout << "Through time " << theTime << ", " << TotalInteractions << " total exchanges; ";
				switch (termination_criterion)
				{
					case -2:
						if (!LnMRSsDataUpToDate)
							ComputeLnMRSsDistribution();
						std::cout << "current L2 s.d. in MRS = " << LnMRSsData.L2StdDev() << std::endl;
						break;
					case -1:
						if (!LnMRSsDataUpToDate)
							ComputeLnMRSsDistribution();
						std::cout << "current L2 s.d. in MRS = " << LnMRSsData.L2StdDev() << std::endl;
						break;
					case 0:
						if (!LnMRSsDataUpToDate)
							ComputeLnMRSsDistribution();
						std::cout << "current max s.d. in MRS = " << LnMRSsData.LinfStdDev() << std::endl;
						break;
					case 1:
						std::cout << "relative increase in ∑U = ";
						std::cout.width(10);
						std::cout.precision(10);
						std::cout << ComputeRelativeIncreaseInSumOfUtilities() << std::endl;
						break;
					case 2:
						std::cout << "increase in ∑U = "; 
						std::cout.width(10);
						std::cout.precision(10);
						std::cout << ComputeIncreaseInSumOfUtilities() << std::endl;
						break;
				}	// switch...
			}	//	theTime...
		
		//	Store the sum of utilities if it will be needed next period for either termincation check or printing
		//
		if (termination_criterion > 0)
			if (((theTime > CheckTerminationThreshhold) && ((theTime + 1) % CheckTerminationPeriod == 0)) || ((PrintIntermediateOutput) && ((theTime + 1) % IntermediateOutputPrintPeriod == 0)))
				LastSumOfUtilities = ComputeSumOfUtilities();
	}	//	do...
	while (!Converged);
	//
	//	Agents are either equilibrated or user has asked for termination; display stats for the former case
	//
	if (!Converged)
	{	
		std::cout << "Terminated by user!" << std::endl;
		return 0;
	}
	else	//	the economy has converged...
	{
		std::cout << "Equilibrium achieved at time " << theTime << " via " << TotalInteractions << " interactions" << std::endl;
		std::cout << std::endl;
		sprintf(OutputStr, "Equilibration #%d ended at time: ", NumberOfEquilibrationsSoFar);
		WriteActualTimeAndDate(OutputStr);
		std::cout << std::endl;
		if (PrintConvergenceStats)
			ConvergenceStatistics(Volume);
		return TotalInteractions;
	}
}	//	AgentPopulation::Equilibrate

void	AgentPopulation::ConvergenceStatistics(CommodityArray VolumeStats)
{
	Data InitialMarketWealthData, FinalMarketWealthData, DeltaMarketWealthData;
	Data FinalOwnWealthData, DeltaOwnWealthData, DeltaUtility;
	double price, AgentsInitialWealth, AgentsFinalWealth;
	
	for (int i = 1; i <= NumberOfAgents; ++i)
	{
		//	First compute agent wealth one commodity at a time...
		//
		AgentsInitialWealth = 0.0;
		AgentsFinalWealth = 0.0;
		for (int j = 1; j <= NumberOfCommodities; ++j)
		{
			price = LnMRSsData.GetData(j)->GetExpAverage();
			AgentsInitialWealth += Agents[i]->GetEndowment(j) * price;
			AgentsFinalWealth += Agents[i]->GetAllocation(j) * price;
		}
		InitialMarketWealthData.AddDatuum(AgentsInitialWealth);
		FinalMarketWealthData.AddDatuum(AgentsFinalWealth);
		DeltaMarketWealthData.AddDatuum(AgentsFinalWealth - AgentsInitialWealth);
		
		//	The following line is a minor optimization...variable name should include 'own' to be mnemonic...
		//
		AgentsFinalWealth = Agents[i]->Wealth(Agents[i]->GetCurrentMRSs());
		FinalOwnWealthData.AddDatuum(AgentsFinalWealth);
		DeltaOwnWealthData.AddDatuum(AgentsFinalWealth - Agents[i]->GetInitialWealth());
		DeltaUtility.AddDatuum(Agents[i]->Utility() - Agents[i]->GetInitialUtility());
	}	//	for i...
	std::cout << "Average initial wealth (@ market prices) = " << InitialMarketWealthData.GetAverage() << "; standard deviation = " << InitialMarketWealthData.GetStdDev() << std::endl;
	std::cout << "Average final wealth (@ market prices)   = " << FinalMarketWealthData.GetAverage() << "; standard deviation = " << FinalMarketWealthData.GetStdDev() << std::endl;
	std::cout << "Average change in market wealth          = " << DeltaMarketWealthData.GetAverage() << "; standard deviation = " << DeltaMarketWealthData.GetStdDev() << std::endl;
	std::cout << "Average final wealth (@ own prices)      = " << FinalOwnWealthData.GetAverage() << "; standard deviation = " << FinalOwnWealthData.GetStdDev() << std::endl;
	std::cout << "Average change in own wealth             = " << DeltaOwnWealthData.GetAverage() << "; standard deviation = " << DeltaOwnWealthData.GetStdDev() << std::endl;
	std::cout << "Minimum increase in utility              = " << DeltaUtility.GetMin() << "; maximum increase = " << DeltaUtility.GetMax() << std::endl;
	std::cout << "Final sum of utilities                   = " << ComputeSumOfUtilities() << std::endl;
	std::cout << std::endl;
	if (PrintFinalCommodityList)
		for (int i = 1; i <= NumberOfCommodities; ++i)
		{
			std::cout << "Commodity " << i << ": volume = " << VolumeStats[i];
			std::cout << "; avg. MRS = " << LnMRSsData.GetData(i)->GetExpAverage() << "; s.d. = " << LnMRSsData.GetData(i)->GetStdDev() << std::endl;
		}
}	//	AgentPopulation::ConvergenceStatistics()

void	AgentPopulation::CompareTwoAgents(AgentPtr Agent1, AgentPtr Agent2)
{
	for (int j = 1; j <= NumberOfCommodities; ++j)
	{
		if (Agent1->GetAlpha(j) != Agent2->GetAlpha(j))
			std::cout << "Bad alpha copying!" << std::endl;
		if (Agent1->GetEndowment(j) != Agent2->GetEndowment(j))
			std::cout << "Bad endowment copying!" << std::endl;
		if (Agent1->GetAllocation(j) != Agent2->GetAllocation(j))
			std::cout << "Bad allocation copying!" << std::endl;
		if (Agent1->GetInitialMRS(j) != Agent2->GetInitialMRS(j))
			std::cout << "Bad InitialMRSs copying!" << std::endl;
		if (Agent1->GetCurrentMRS(j) != Agent2->GetCurrentMRS(j))
			std::cout << "Bad CurrentMRSs copying!" << " " << j << std::endl;
	}
	if (Agent1->GetInitialUtility() != Agent2->GetInitialUtility())
		std::cout << "Bad InitialUtility copying!" << std::endl;
	if (Agent1->GetInitialWealth() != Agent2->GetInitialWealth())
		std::cout << "Bad InitialWealth copying!" << std::endl;
}	//	AgentPopulation::CompareTwoAgents()

void AgentPopulation::ShockAgentPreferences()
{
	AgentPtr ActiveAgent;
	double oldPref, newPref, pref;
	
	int CommodityToShock = RNG.IntegerInRange(1, NumberOfCommodities);
	bool sign = RNG.IntegerInRange(0,1);
	double shock = RNG.RealInRange(MinShock, MaxShock);
	
	for (int AgentIndex = 1; AgentIndex <= NumberOfAgents; ++AgentIndex)
	{
		ActiveAgent = Agents[AgentIndex];
		oldPref = ActiveAgent->GetAlpha(CommodityToShock);
		if (sign)
			newPref = oldPref * shock;
		else
			newPref = oldPref/shock;
		ActiveAgent->SetAlpha(CommodityToShock, newPref);
		
		for (int CommodityIndex = 1; CommodityIndex <= NumberOfCommodities; ++CommodityIndex)
		{
			pref = ActiveAgent->GetAlpha(CommodityIndex);
			ActiveAgent->SetAlpha(CommodityIndex, pref/(1.0 - oldPref + newPref));
		}	//	for (CommodityIndex...
	}	//	for (AgentIndex...
}	//	AgentPopulation::ShockAgentPreferenes()

/*================= End of Methods =================*/

void	InitMiscellaneous()
{
	//	Print some configuration information...
	//
	std::cout << "B I L A T E R A L   E X C H A N G E" << std::endl <<std::endl;
	std::cout << "Robert Axtell" << std::endl;
	std::cout << "Brookings Institution -> George Mason University" << std::endl;
	std::cout << "Version " << Version << std::endl << std::endl;
	
	std::cout << "Number of agents = " << NumberOfAgents << std::endl;
	std::cout << "Number of commodities = " << NumberOfCommodities << std::endl << std::endl;
	
	std::cout << "Agent valuations must differ by more than epsilon = " << trade_eps << std::endl;
	std::cout << "  in order for trade to occur" << std::endl << std::endl;
	
	std::cout << "A time period is defined by " << 2 * PairwiseInteractionsPerPeriod << " agents being active (";
	std::cout << PairwiseInteractionsPerPeriod << " pairings)" << std::endl << std::endl;

	std::cout << "Termination is checked each " << CheckTerminationPeriod << " period";
	if (CheckTerminationPeriod > 1)
		std::cout << "s";
	std::cout << std::endl;
	std::cout << "  once the first " << CheckTerminationThreshhold << " periods have gone by" << std::endl;
	std::cout << "Termination occurs once the ";
	switch (termination_criterion)
	{
		case -1:
			std::cout << "L2 norm of the standard deviation" << std::endl;
			std::cout << "  in agent MRSs falls below ";
			break;
		case 0:
			std::cout << "maximum standard deviation" << std::endl;
			std::cout << "  in agent MRSs falls below ";
			break;
		case 1:
			std::cout << "relative increase in the" << std::endl;
			std::cout << "  sum of utilities is less than ";
			break;
		case 2:
			std::cout << "increase in the sum" << std::endl;
			std::cout << "  of utilities is less than ";
			break;
	}
	std::cout << termination_eps << std::endl << std::endl;
	
	std::cout << "Agents are being processed ";
	if (DefaultSerialExecution)							//	Can't use 'Population.Serial' here because the population has not been initialized
		std::cout << "serially" << std::endl;
	else
	{
		std::cout << "in parallel" << std::endl;
		std::cout << "Number of agents rearranged in the agent list each period: " << AgentsToRandomize << std::endl;
	}
	std::cout << std::endl;
	
	if (!UseRandomSeed)
		std::cout << "Seed: " << NonRandomSeed << std::endl;
	std::cout << std::endl;
	
	std::cout << "Number of equilibrations to be performed: " << RequestedEquilibrations << std::endl;
	if (RequestedEquilibrations > 1)
	{
		std::cout << "Agents are ";
		if (SameAgentInitialCondition)
			std::cout << "not ";
		std::cout << "given new preferences and endowments for each equilibration" << std::endl;
	}	
	std::cout << std::endl;
	
	MemoryState.WriteMemoryRequirements();
}  //	InitMiscellaneous()

int	main()
{
	//	First, initialize variously...
	//
	InitMiscellaneous();
	
	AgentPopulationPtr PopulationPtr = new AgentPopulation;
	
	//	Equilibrate the agent economy once...
	//
	int EquilibrationNumber = 1;
	long sum = PopulationPtr->Equilibrate(EquilibrationNumber);
	long sum2 = sum*sum;
	
	//	Equilibrate again if the user has requested this...
	//
	long interactions;
	
	EquilibrationNumber = 2;
	while (EquilibrationNumber <= RequestedEquilibrations)
	{
		if (SameAgentInitialCondition)
			PopulationPtr->Reset();
		else
			PopulationPtr->Init();
		interactions = PopulationPtr->Equilibrate(EquilibrationNumber);
		sum += interactions;
		sum2 += interactions*interactions;
		++EquilibrationNumber;
	}
	double avg = (double) sum/(EquilibrationNumber-1);
	std::cout << "Average number of interactions: " << avg;
	if (EquilibrationNumber > 2)
		std::cout << "; std. dev.: " << sqrt((sum2 - (EquilibrationNumber - 1) * avg * avg)/(EquilibrationNumber - 2));
}
