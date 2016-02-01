// Copyright 2015 <Stefan McCabe>
/*
Bilateral Exchange Model
Originally written by Rob Axtell
Extended by Stefan McCabe

12/11/2015
*/

//  Rob:
//  Additional work to do on this code:
//
//  A. Multithread the main 'Equilibrate()' method
//  B. Rationalize use of Init()s, Reset()s and constructors
//  C. Use my 'timer()' object instead of what is here...
//  D. Generalize to CES preferences?
//  E. Shrink below 1000 lines?

// I used this thread-safe random number generator in another project; it will probably be useful later.
//
// int IntegerInRange(int min, int max) {
//     static thread_local mt19937* generator = nullptr;
//     if (!generator) generator = new mt19937(clock() + hash<thread::id>()(this_thread::get_id()));
//     uniform_int_distribution<int> distribution(min, max);
//     return distribution(*generator);
// }

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <gflags/gflags.h>
#define ELPP_THREAD_SAFE
#include "./easylogging++.h"
#include <libconfig.h++>
#include <vector>
using namespace std;

// Initialize the logger. This should come immediately after the #includes are finished.
INITIALIZE_EASYLOGGINGPP

// Declare gflags.
DEFINE_bool(test, true, "Test GFlags");
DEFINE_string(file, "parameters.cfg", "Specify the configuation file containing the model parameters.");

// addressing the RNG first

bool UseRandomSeed;
unsigned int NonRandomSeed;


mt19937 rng(0);
// We seed 0 by default; in main we will seed the RNG properly.

double Version = 1.0;

int NumberOfAgents;
int NumberOfCommodities;

int PairwiseInteractionsPerPeriod;

double alphaMin;
double alphaMax;
unsigned int wealthMin;
unsigned int wealthMax;

bool DefaultSerialExecution;   // FALSE means parallel
int AgentsToRandomize;
//  This is done only for parallel interactions

int RandomizationMethod;            // 0 means swap pairs
                                    // 1 means "threaded"
int RequestedEquilibrations;
bool SameAgentInitialCondition;

double trade_eps;
double exp_trade_eps;

int termination_criterion;              // -2 means termination based on time
                                        // -1 means use L2 norm on MRSs
                                        //  0 means use L∞ norm on MRSs
                                        //  1 means use relative increase in V
                                        //  2 means use absolute increase in V
long long TerminationTime;

double termination_eps;    //  10^-1
long long CheckTerminationThreshold; // was 100000
//  1.5 x 10^9 periods ~ 60 x 10^9 interactions

int CheckTerminationPeriod;
//  Typical values:
//  A = 100,        N = 10,     use 1
//  A = 100,        N = 50,     use 25
//  A = 100,        N = 100,    use 100
//  A = 100,        N = 500,    use 2500
//  A = 100,        N = 1000,   use 10,000
//  A = 100,        N = 5000,   use 250,000
//  A = 100,        N = 10,000, use 1,000,000
//  A = 1000,       N = 2,      use 1
//  A = 10,000,     N = 2,      use 10
//  A = 100,000,    N = 2,      use 100
//  A = 1,000,000,  N = 2,      use 1000
//  A = 10,000,000, N = 2,      use 10,000

bool ShockPreferences;
int ShockPeriod;
double MinShock;
double MaxShock;

bool debug;
bool PrintEndowments;
bool PrintIntermediateOutput;
int IntermediateOutputPrintPeriod;  // 20x10^6
bool PrintConvergenceStats;
bool PrintFinalCommodityList;

// DISTRIBUTIONS
uniform_int_distribution<unsigned long> randomAgent;  
uniform_int_distribution<unsigned long> randomCommodity;
uniform_int_distribution<int> randomBinary;
uniform_real_distribution<double> randomShock;
uniform_real_distribution<double> randomAlpha;
uniform_real_distribution<double> randomWealth;

// TODO: A lot of things depend on this typedef but I would like to make 
// to delay declaring the size of the array so that I can take number 
// of commodites as a parameter (user input). For now I'm admitting defeat and still
// hard-coding the number of commodities. 
// typedef double CommodityArray[NumberOfCommodities+1];
//  Size: NumberOfCommodities * 8 bytes
// New typedef *seems* to work, provided that I resize everything in the constructors.
 typedef vector<double> CommodityArray;

typedef CommodityArray *CommodityArrayPtr;

double inline Dot(CommodityArrayPtr vector1, CommodityArrayPtr vector2);

class MemoryObject {
    long long start;

public:
    MemoryObject();
    void WriteMemoryRequirements();
} MemoryState;


class Data {        //      Size:
    int N;          //      4 bytes
    double min;     //      8 bytes
    double max;     //      8 bytes
    double sum;     //      8 bytes
    double sum2;    //      8 bytes
public:
    Data();
    void Init();
    void AddDatuum(double Datuum);
    int GetN() { return N; }
    double GetMin() { return min; }
    double GetMax() { return max; }
    double GetDelta() { return max - min; }
    double GetAverage() {
        if (N > 0) {
            return sum/N;
        } else {
            return 0.0;
        }}
        double GetExpAverage() { return exp(GetAverage()); }
        double GetVariance();
        double GetStdDev() { return sqrt(GetVariance()); }
};  //  Total:  36 bytes

typedef Data *DataPtr;

class CommodityData {   //  Size:
    // Data data[NumberOfCommodities+1];
    vector<Data> data;
    //Data data = new Data[NumberOfCommodities + 1];
                        // NumberOfCommodities * 8 bytes
public:
    CommodityData();
    void Clear();
    DataPtr GetData(size_t index) { return &data[index]; }
    double L2StdDev();
    double LinfStdDev();
};                       //  Total:  NumberOfCommodities * 8 bytes
                        //  Example: NumberOfCommodities = 100, size = 800 bytes


class Agent {                       //  Size:
    CommodityArray alphas;          //      NumberOfCommodities * 8 bytes
    CommodityArray endowment;       //      NumberOfCommodities * 8 bytes
    CommodityArray initialMRSs;     //      NumberOfCommodities * 8 bytes
    CommodityArray allocation;      //      NumberOfCommodities * 8 bytes
    CommodityArray currentMRSs;     //      NumberOfCommodities * 8 bytes
    Agent();
    double initialUtility;          //      8 bytes
    double initialWealth;           //      8 bytes

public:
    Agent(int size);
    void Init();
    void Reset();
    double GetAlpha(size_t CommodityIndex) {
        return alphas[CommodityIndex];
    }
    void SetAlpha(size_t CommodityIndex, double alpha) {
        alphas[CommodityIndex] = alpha;
    }
    double GetEndowment(size_t CommodityIndex) {
        return endowment[CommodityIndex];
    }
    double MRS(size_t CommodityIndex, size_t Numeraire);
    void ComputeMRSs();
    double GetInitialMRS(size_t CommodityIndex) {
        return initialMRSs[CommodityIndex];
    }
    double GetAllocation(size_t CommodityIndex) {
        return allocation[CommodityIndex];
    }
    void IncreaseAllocation(size_t CommodityIndex, double amount) {
        allocation[CommodityIndex] += amount;
    }
    double GetCurrentMRS(size_t CommodityIndex) {
        return currentMRSs[CommodityIndex];
    }
    CommodityArrayPtr GetCurrentMRSs() {
        return &currentMRSs;
    }
    double Utility();
    double GetInitialUtility() {
        return initialUtility;
    }
    double Wealth(CommodityArrayPtr prices) {
        return Dot(&allocation, prices);
    }
    double GetInitialWealth() {
        return initialWealth;
    }
};
                //  Total: 20 + NumberOfCommodities * 40 bytes
                //  Example: NumberOfCommodities = 2, size = 20 + 80 = 100
                //  Example: NumberOfCommodities = 100, size = 20 + 4000 = 4020
typedef Agent *AgentPtr;

class AgentPopulation {                             //  Size:
    //AgentPtr Agents[NumberOfAgents + 1];            //  NumberOfAgents * 4 bytes
    vector<AgentPtr> Agents;
    CommodityArray Volume;
    CommodityData AlphaData, EndowmentData, LnMRSsData;
    bool LnMRSsDataUpToDate;
    void ComputeLnMRSsDistribution();
    double LastSumOfUtilities;
    double ComputeSumOfUtilities();
    double ComputeIncreaseInSumOfUtilities() {
        return ComputeSumOfUtilities() - LastSumOfUtilities;
    }
    double ComputeRelativeIncreaseInSumOfUtilities() {
        return ComputeIncreaseInSumOfUtilities() / LastSumOfUtilities;
    }
    void GetSerialAgentPair(AgentPtr& Agent1, AgentPtr& Agent2);
    size_t ActiveAgentIndex;
    void RandomizeAgents(int NumberToRandomize);
    void RandomizeAgents2(int NumberToRandomize);
    void GetParallelAgentPair(AgentPtr& Agent1, AgentPtr& Agent2);
    void(AgentPopulation::*GetAgentPair) (AgentPtr& Agent1, AgentPtr& Agent2);

public:
    AgentPopulation(int size);                              //  Constructor...
    void Init();
    void Reset();
    long long Equilibrate(int NumberOfEquilibrationsSoFar);
    void ConvergenceStatistics(CommodityArray VolumeStats);
    void CompareTwoAgents(AgentPtr Agent1, AgentPtr Agent2);
    void ShockAgentPreferences();
};
//  Example:    NumberOfAgents = 100, size = 400 bytes

typedef AgentPopulation *AgentPopulationPtr;

double inline Dot(CommodityArrayPtr vector1, CommodityArrayPtr vector2) {
    double sum = 0.0;
    for (size_t i = 1; i <= static_cast<size_t>(NumberOfCommodities); ++i) {
        sum += (*vector1)[i] * (*vector2)[i];
    }
    return sum;
}   //  Dot

/*===========
    ==Methods==
    ===========*/

    MemoryObject::MemoryObject():
    start(0)
{}  //  MemoryObject::MemoryObject()


void MemoryObject::WriteMemoryRequirements() {
    LOG(DEBUG) << "Size of Agent in memory: " << sizeof(Agent) << " bytes";
    LOG(DEBUG) << "Size of Data in memory: " << sizeof(Data) << " bytes";
    LOG(DEBUG) << "Size of CommodityData in memory: " << sizeof(CommodityData) << " bytes";
    LOG(DEBUG) << "Size of AgentPopulation in memory: " << sizeof(AgentPopulation) << " bytes";
    LOG(DEBUG) << "Total bytes required (approximate): " << \
    (sizeof(rng) + sizeof(MemoryObject) + sizeof(CommodityData) \
        + sizeof(AgentPopulation) + sizeof(Agent) * static_cast<size_t>(NumberOfAgents)) << " bytes";
}

Data::Data():
N(0), min(1000000.0), max(0.0), sum(0.0), sum2(0.0)
{}  //  Data::Data()

void Data::Init() {
    N = 0;
    min = 1000000.0;    //  10^6
    max = 0.0;
    sum = 0.0;
    sum2 = 0.0;
}   //  Data::Data()

void Data::AddDatuum(double Datuum) {
    N = N + 1;
    if (Datuum < min) {
        min = Datuum;
    }
    if (Datuum > max) {
        max = Datuum;
    }
    sum += Datuum;
    sum2 += Datuum * Datuum;
}   //  Data::AddDatuum()

double Data::GetVariance() {
    double avg, arg;

    if (N > 1) {
        avg = GetAverage();
        arg = sum2 - N * avg * avg;
        return arg / (N - 1);
    } else {
        return 0.0;
    }
}   //  Data::GetVariance()

CommodityData::CommodityData() {
data.resize(static_cast<size_t>(NumberOfCommodities+1));
//    data = 
}

void CommodityData::Clear() {
    //  Initialize data objects
    //
    for (size_t i = 1; i <= static_cast<size_t>(NumberOfCommodities); ++i) {
        data[i].Init();
    }
}   //  CommodityData::Clear()

double CommodityData::L2StdDev() {
    double sum = 0.0;

    // Instead of computing the standard deviation for each commodity (and calling sqrt() N times),
    // find the largest variance and from it get the std dev
    for (size_t i = 1; i <= static_cast<size_t>(NumberOfCommodities); ++i) {
        sum += data[i].GetVariance();
    }
    return sqrt(sum);
}   //  CommodityData::L2StdDev

double  CommodityData::LinfStdDev() {
    double var, max = 0.0;

    // Instead of computing the standard deviation for each commodity (and calling sqrt() N times),
    // find the largest variance and from it get the std dev
    for (size_t i = 1; i <= static_cast<size_t>(NumberOfCommodities); ++i) {
        var = data[i].GetVariance();
        if (var > max) {
            max = var;
        }
    }
    return sqrt(max);
}   //  CommodityData::LinfStdDev


Agent::Agent(int size): initialUtility(0.0), initialWealth(0.0) {
    alphas.resize(static_cast<size_t>(size+1));
    endowment.resize(static_cast<size_t>(size+1));   
    initialMRSs.resize(static_cast<size_t>(size+1)); 
    allocation.resize(static_cast<size_t>(size+1));  
    currentMRSs.resize(static_cast<size_t>(size+1));
    Init();
}

void Agent::Init() {
    size_t CommodityIndex;

    //  First generate and normalize the exponents...
    //
    double sum = 0.0;
    for (CommodityIndex = 1; CommodityIndex <= static_cast<size_t>(NumberOfCommodities); ++CommodityIndex) {
        alphas[CommodityIndex] = randomAlpha(rng);
        sum += alphas[CommodityIndex];
    }
    //  Next, fill up the rest of the agent fields...
    //
    for (CommodityIndex = 1; CommodityIndex <= static_cast<size_t>(NumberOfCommodities); ++CommodityIndex) {
        alphas[CommodityIndex] = alphas[CommodityIndex] / sum;
        endowment[CommodityIndex] = randomWealth(rng);
        allocation[CommodityIndex] = endowment[CommodityIndex];
        initialMRSs[CommodityIndex] = MRS(CommodityIndex, 1);
        currentMRSs[CommodityIndex] = initialMRSs[CommodityIndex];
    }
    initialUtility = Utility();
    initialWealth = Wealth(&initialMRSs);
}   //  Agent:Init()

void Agent::Reset() {
    for (size_t CommodityIndex = 1; CommodityIndex <= static_cast<size_t>(NumberOfCommodities); ++CommodityIndex) {
        allocation[CommodityIndex] =  endowment[CommodityIndex];
        currentMRSs[CommodityIndex] = initialMRSs[CommodityIndex];
    }
}   //  Agent::Reset

double Agent::MRS (size_t CommodityIndex, size_t Numeraire) {
    return (alphas[CommodityIndex] * allocation[Numeraire]) / (alphas[Numeraire] * allocation[CommodityIndex]);
}   //  Agent::MRS()

void Agent::ComputeMRSs() {
    size_t i;
    currentMRSs[1] = 1.0;
    for (i = 2; i <= static_cast<size_t>(NumberOfCommodities); ++i) {
        currentMRSs[i] = MRS(i, 1);
    }
}   //      Agent::ComputeMRSs()

double Agent::Utility() {
    double product = 1.0;

    for (size_t i = 1; i <= static_cast<size_t>(NumberOfCommodities); ++i) {
        product *= pow(allocation[i], alphas[i]);
    }

    return product;
}   //  Agent::Utility()

void AgentPopulation::ComputeLnMRSsDistribution() {
    LnMRSsData.Clear();

    for (size_t i = 1; i <= static_cast<size_t>(NumberOfAgents); ++i) {
        Agents[i]->ComputeMRSs();
        for (size_t j = 1; j <= static_cast<size_t>(NumberOfCommodities); ++j) {
            LnMRSsData.GetData(j)->AddDatuum(log(Agents[i]->GetCurrentMRS(j)));
        }
    }       //  for i...
    LnMRSsDataUpToDate = true;
}   //  AgentPopulation::ComputeLnMRSsDistribution()

double AgentPopulation::ComputeSumOfUtilities() {
    double sum = 0.0;

    for (size_t i = 1; i <= static_cast<size_t>(NumberOfAgents); ++i) {
        sum += Agents[i]->Utility();
    }
    return sum;
}   //  AgentPopulation::ComputeSumOfUtilities()

void inline AgentPopulation::GetSerialAgentPair(AgentPtr& Agent1, AgentPtr& Agent2) {
    Agent1 = Agents[randomAgent(rng)];
    do {
        Agent2 = Agents[randomAgent(rng)];
    } while (Agent2 == Agent1);
}   //  AgentPopulation::GetSerialAgentPair()

void inline AgentPopulation::RandomizeAgents(int NumberToRandomize) {
    // Note:  If NumberToRandomize < 2 then it is increased to 2, so that 1 pair is swapped each period
    //
    size_t AgentIndex1, AgentIndex2;
    int PairsToRandomize;
    AgentPtr TempPtr;

    PairsToRandomize = NumberToRandomize / 2;
    if (PairsToRandomize < 1) {
        PairsToRandomize = 1;
    }

    for (int i = 1; i <= PairsToRandomize; ++i) {
        AgentIndex1 = randomAgent(rng);
        do {
            AgentIndex2 = randomAgent(rng);
        } while (AgentIndex1 == AgentIndex2);

        TempPtr = Agents[AgentIndex1];
        Agents[AgentIndex1] = Agents[AgentIndex2];
        Agents[AgentIndex2] = TempPtr;
    }
}   //  AgentPopulation::RandomizeAgents()

void inline AgentPopulation::RandomizeAgents2(int NumberToRandomize) {
    // Note:  The NumberToRandomize need not be even, just > 1
    size_t AgentIndex1, AgentIndex2;
    AgentPtr TempPtr;

    AgentIndex1 = randomAgent(rng);
    TempPtr = Agents[AgentIndex1];

    if (NumberToRandomize == 1) {
        NumberToRandomize = 2;
    }
    for (int i = 1; i < NumberToRandomize; ++i) {
        AgentIndex2 = randomAgent(rng);
        Agents[AgentIndex1] = Agents[AgentIndex2];
        AgentIndex1 = AgentIndex2;
    }
    Agents[AgentIndex2] = TempPtr;
}   //  AgentPopulation::RandomizeAgents2()

void inline AgentPopulation::GetParallelAgentPair (AgentPtr& Agent1, AgentPtr& Agent2) {
    if (ActiveAgentIndex > static_cast<size_t>(NumberOfAgents)) {
        ActiveAgentIndex = 1;
        RandomizeAgents(AgentsToRandomize);
    }
    Agent1 = Agents[ActiveAgentIndex];
    Agent2 = Agents[++ActiveAgentIndex];
    ActiveAgentIndex++;
}   //  AgentPopulation::GetParallelAgentPair()

AgentPopulation::AgentPopulation(int size):
AlphaData(), EndowmentData(), LnMRSsData(), LnMRSsDataUpToDate(true), LastSumOfUtilities(0.0), ActiveAgentIndex(0), GetAgentPair(NULL) {
    Volume.resize(static_cast<size_t>(size+1));
    Agents.resize(static_cast<size_t>(NumberOfAgents+1));
    for (size_t i = 1; i <= static_cast<size_t>(NumberOfAgents); i++) {
        Agents[i] = new Agent(NumberOfCommodities);
    }

    if (DefaultSerialExecution) {
        GetAgentPair = &AgentPopulation::GetSerialAgentPair;
    } else
    GetAgentPair = &AgentPopulation::GetParallelAgentPair;
}   // Constructor...

void AgentPopulation::Init() {
    AgentPtr ActiveAgent;
    size_t CommodityIndex;
    Data InitialOwnWealthData;

    ActiveAgentIndex = 1;
    AlphaData.Clear();
    EndowmentData.Clear();
    LnMRSsData.Clear();

    //  Cycle through the agents...
    //
    for (size_t AgentIndex = 1; AgentIndex <= static_cast<size_t>(NumberOfAgents); ++AgentIndex) {
        ActiveAgent = Agents[AgentIndex];
        ActiveAgent->Init();

        //  Next, fill up the rest of the agent fields...
        //
        for (CommodityIndex = 1; CommodityIndex <= static_cast<size_t>(NumberOfCommodities); ++CommodityIndex) {
            AlphaData.GetData(CommodityIndex)->AddDatuum(ActiveAgent->GetAlpha(CommodityIndex));
            EndowmentData.GetData(CommodityIndex)->AddDatuum(ActiveAgent->GetEndowment(CommodityIndex));
            LnMRSsData.GetData(CommodityIndex)->AddDatuum(log(ActiveAgent->GetInitialMRS(CommodityIndex)));
        }   //  for (CommodityIndex...
            InitialOwnWealthData.AddDatuum(ActiveAgent->GetInitialWealth());
    }   //  for (AgentIndex...

        LnMRSsDataUpToDate = true;
        LastSumOfUtilities = ComputeSumOfUtilities();

    //  Finally, display stats on the instantiated population...
    //
        if (PrintEndowments) {
            LOG(INFO) << "Initial endowments:";
            for (CommodityIndex = 1; CommodityIndex <= static_cast<size_t>(NumberOfCommodities); ++CommodityIndex) {
                LOG(INFO) << "Commodity " << CommodityIndex << ": <exp.> = " << AlphaData.GetData(CommodityIndex)->GetAverage() << "; s.d. = " << AlphaData.GetData(CommodityIndex)->GetStdDev() <<
                "; <endow.> = " << EndowmentData.GetData(CommodityIndex)->GetAverage() << "; s.d. = " << EndowmentData.GetData(CommodityIndex)->GetStdDev() << "; <MRS> = " << 
                LnMRSsData.GetData(CommodityIndex)->GetExpAverage() << "; s.d. = " << LnMRSsData.GetData(CommodityIndex)->GetStdDev();
            }
        }
        LOG(INFO) << "Average initial wealth (@ own prices) = " << InitialOwnWealthData.GetAverage() << "; standard deviation = " << InitialOwnWealthData.GetStdDev();
        LOG(INFO) << "Initial sum of utilities = " << LastSumOfUtilities;
}   //  AgentPopulation::Init()

void AgentPopulation::Reset() {
    for (size_t AgentIndex = 1; AgentIndex <= static_cast<size_t>(NumberOfAgents); ++AgentIndex) {
        Agents[AgentIndex]->Reset();
    }
}   //  AgentPopulation::Reset

long long AgentPopulation::Equilibrate(int NumberOfEquilibrationsSoFar) {
    char OutputStr[63];
    bool Converged = false;
    long long theTime = 0;
    long long TotalInteractions = 0;
    AgentPtr Agent1, Agent2, LargerMRSAgent, SmallerMRSAgent;
    size_t Commodity1, Commodity2;
    double MRSratio12, MRSratio;
    double LAgentalpha1, LAgentalpha2, LAgentx1, LAgentx2, SAgentalpha1, SAgentalpha2, SAgentx1, SAgentx2;
    double num, denom, delta1, delta2;
    double Agent1PreTradeUtility, Agent2PreTradeUtility;

    LOG(INFO) << "Equilibration #" << NumberOfEquilibrationsSoFar << " starting";
    //  Next, initialize some variables...
    //
    for (size_t i = 1; i <= static_cast<size_t>(NumberOfCommodities); ++i) {
        Volume[i] = 0.0;
    }

    //  Start up the exchange process here...
    //
    do {
        LnMRSsDataUpToDate = false;  //  Since these data are gonna change...

        ++theTime;

        if ((ShockPreferences) && (theTime % ShockPeriod == 0)) {
            ShockAgentPreferences();
        }
        for (int i = 1; i <= PairwiseInteractionsPerPeriod; ++i) {
            //  First, select the agents to be active...
            //
            (this->*GetAgentPair) (Agent1, Agent2);

            if (debug) {
                Agent1PreTradeUtility = Agent1->Utility();
                Agent2PreTradeUtility = Agent2->Utility();
            }
            //  Next, select the commodities to trade...
            //
            if (NumberOfCommodities == 2) {
                Commodity1 = 1;
                Commodity2 = 2;
            } else {
                Commodity1 = randomCommodity(rng);
                do {
                    Commodity1 = randomCommodity(rng);
                } while (Commodity2 == Commodity1);
            }
            //  Compare MRSs...

            MRSratio12 = Agent1->MRS(Commodity2, Commodity1) / Agent2->MRS(Commodity2, Commodity1);
            if (MRSratio12 > 1.0) {
                MRSratio = MRSratio12;
            } else {
                MRSratio = 1.0/MRSratio12;
            }
            // cout << MRSratio << endl;
            if (MRSratio >= exp_trade_eps)  {  //  do exchange
                if (MRSratio12 > 1.0) {
                    LargerMRSAgent = Agent1;
                    SmallerMRSAgent = Agent2;
                } else {
                    LargerMRSAgent = Agent2;
                    SmallerMRSAgent = Agent1;
                }
                //  Here are the guts of bilateral Walrasian exchange
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
            }   //  if (MRSratio...

                if (debug) {
                    if (Agent1->Utility() < Agent1PreTradeUtility) {
                        LOG(WARNING) << "!!!Utility decreasing trade by agent #1!!!  Actual utility change = " << Agent1->Utility() - Agent1PreTradeUtility;
                    }
                    if (Agent2->Utility() < Agent2PreTradeUtility) {
                        LOG(WARNING) << "!!!Utility decreasing trade by agent #2!!!  Actual utility change = " << Agent2->Utility() - Agent2PreTradeUtility;
                    }
                }
        }   //  for i...

        //  Check for termination...
        //
        if ((theTime > CheckTerminationThreshold) && (theTime % CheckTerminationPeriod == 0))  {
            switch (termination_criterion) {
                case -2:
                if (theTime >= TerminationTime) {
                    Converged = true;
                }
                break;
                case -1:  //  Termination based on L2 norm of MRS distribution
                ComputeLnMRSsDistribution();
                if (LnMRSsData.L2StdDev() < termination_eps) {
                    Converged = true;
                }
                break;
                case 0:  //  Termination based on L∞ norm of MRS distribution
                ComputeLnMRSsDistribution();
                if (LnMRSsData.LinfStdDev() < termination_eps) {
                    Converged = true;
                }
                break;
                case 1:  //  Termination based on relative increase in V
                if (ComputeRelativeIncreaseInSumOfUtilities() < termination_eps) {
                    Converged = true;
                }
                break;
                case 2:  //  Termination based on absolute increase in V
                if (ComputeIncreaseInSumOfUtilities() < termination_eps) {
                    Converged = true;
                }
                break;
            }   //  switch...
        }
        //  Display stats if the time is right...
        //
        if (PrintIntermediateOutput) {
            if (theTime % IntermediateOutputPrintPeriod == 0) {
                LOG(INFO) << "Through time " << theTime << ", " << TotalInteractions << " total exchanges; ";
                switch (termination_criterion) {
                    case -2:
                    if (!LnMRSsDataUpToDate) {
                        ComputeLnMRSsDistribution();
                    }
                    LOG(INFO) << "current L2 s.d. in MRS = " << LnMRSsData.L2StdDev();
                    break;
                    case -1:
                    if (!LnMRSsDataUpToDate) {
                        ComputeLnMRSsDistribution();
                    }
                    LOG(INFO) << "current L2 s.d. in MRS = " << LnMRSsData.L2StdDev();
                    break;
                    case 0:
                    if (!LnMRSsDataUpToDate) {
                        ComputeLnMRSsDistribution();
                    }
                    LOG(INFO) << "current max s.d. in MRS = " << LnMRSsData.LinfStdDev();
                    break;
                    case 1:
                    LOG(INFO) << "relative increase in ∑U = " << ComputeRelativeIncreaseInSumOfUtilities();
                    break;
                    case 2:
                    LOG(INFO) << "increase in ∑U = " << ComputeIncreaseInSumOfUtilities();
                    break;
                }   // switch...
            }   //  theTime...
        }

        //  Store the sum of utilities if it will be needed next period for either termincation check or printing
        if (termination_criterion > 0) {
            if (((theTime > CheckTerminationThreshold) && ((theTime + 1) % CheckTerminationPeriod == 0)) || ((PrintIntermediateOutput) && ((theTime + 1) % IntermediateOutputPrintPeriod == 0))) {
                LastSumOfUtilities = ComputeSumOfUtilities();
            }
        }
    }  while (!Converged);
    
    //
    //  Agents are either equilibrated or user has asked for termination; display stats for the former case
    //
    if (!Converged) {
        LOG(INFO) << "Terminated by user!";
        return 0;
    } else {    //  the economy has converged...
        LOG(INFO) << "Equilibrium achieved at time " << theTime << " via " << TotalInteractions << " interactions";
        //LOG(INFO) << sprintf(OutputStr, "Equilibration #%d ended", NumberOfEquilibrationsSoFar);
        LOG(INFO) << "Equilibration #" << NumberOfEquilibrationsSoFar << " ended";
        if (PrintConvergenceStats) {
            ConvergenceStatistics(Volume);
        }
        return TotalInteractions;
    }
}   //  AgentPopulation::Equilibrate

void AgentPopulation::ConvergenceStatistics(CommodityArray VolumeStats) {
    Data InitialMarketWealthData, FinalMarketWealthData, DeltaMarketWealthData;
    Data FinalOwnWealthData, DeltaOwnWealthData, DeltaUtility;
    double price, AgentsInitialWealth, AgentsFinalWealth;

    for (size_t i = 1; i <= static_cast<size_t>(NumberOfAgents); ++i) {
        //  First compute agent wealth one commodity at a time...
        //
        AgentsInitialWealth = 0.0;
        AgentsFinalWealth = 0.0;
        for (size_t j = 1; j <= static_cast<size_t>(NumberOfCommodities); ++j) {
            price = LnMRSsData.GetData(j)->GetExpAverage();
            AgentsInitialWealth += Agents[i]->GetEndowment(j) * price;
            AgentsFinalWealth += Agents[i]->GetAllocation(j) * price;
        }
        InitialMarketWealthData.AddDatuum(AgentsInitialWealth);
        FinalMarketWealthData.AddDatuum(AgentsFinalWealth);
        DeltaMarketWealthData.AddDatuum(AgentsFinalWealth - AgentsInitialWealth);

        //  The following line is a minor optimization...variable name should include 'own' to be mnemonic...
        //
        AgentsFinalWealth = Agents[i]->Wealth(Agents[i]->GetCurrentMRSs());
        FinalOwnWealthData.AddDatuum(AgentsFinalWealth);
        DeltaOwnWealthData.AddDatuum(AgentsFinalWealth - Agents[i]->GetInitialWealth());
        DeltaUtility.AddDatuum(Agents[i]->Utility() - Agents[i]->GetInitialUtility());
    }   //  for i...

    LOG(INFO) << "Average initial wealth (@ market prices) = " << InitialMarketWealthData.GetAverage() << "; standard deviation = " << InitialMarketWealthData.GetStdDev();
    LOG(INFO) << "Average final wealth (@ market prices)   = " << FinalMarketWealthData.GetAverage() << "; standard deviation = " << FinalMarketWealthData.GetStdDev();
    LOG(INFO) << "Average change in market wealth          = " << DeltaMarketWealthData.GetAverage() << "; standard deviation = " << DeltaMarketWealthData.GetStdDev();
    LOG(INFO) << "Average final wealth (@ own prices)      = " << FinalOwnWealthData.GetAverage() << "; standard deviation = " << FinalOwnWealthData.GetStdDev();
    LOG(INFO) << "Average change in own wealth             = " << DeltaOwnWealthData.GetAverage() << "; standard deviation = " << DeltaOwnWealthData.GetStdDev();
    LOG(INFO) << "Minimum increase in utility              = " << DeltaUtility.GetMin() << "; maximum increase = " << DeltaUtility.GetMax();
    LOG(INFO) << "Final sum of utilities                   = " << ComputeSumOfUtilities();
    if (PrintFinalCommodityList) {
        for (size_t i = 1; i <= static_cast<size_t>(NumberOfCommodities); ++i) {
            LOG(INFO) << "Commodity " << i << ": volume = " << VolumeStats[i] << "; avg. MRS = " << LnMRSsData.GetData(i)->GetExpAverage() << 
            "; s.d. = " << LnMRSsData.GetData(i)->GetStdDev();
        }
    }
}   //  AgentPopulation::ConvergenceStatistics()

void AgentPopulation::CompareTwoAgents(AgentPtr Agent1, AgentPtr Agent2) {
    for (size_t j = 1; j <= static_cast<size_t>(NumberOfCommodities); ++j) {
        if (Agent1->GetAlpha(j) != Agent2->GetAlpha(j)) {
            LOG(WARNING) << "Bad alpha copying!";
        }
        if (Agent1->GetEndowment(j) != Agent2->GetEndowment(j)) {
            LOG(WARNING) << "Bad endowment copying!";
        }
        if (Agent1->GetAllocation(j) != Agent2->GetAllocation(j)) {
            LOG(WARNING) << "Bad allocation copying!";
        }
        if (Agent1->GetInitialMRS(j) != Agent2->GetInitialMRS(j)) {
            LOG(WARNING) << "Bad InitialMRSs copying!";
        }
        if (Agent1->GetCurrentMRS(j) != Agent2->GetCurrentMRS(j)) {
            LOG(WARNING) << "Bad CurrentMRSs copying!";
        }
    }
    if (Agent1->GetInitialUtility() != Agent2->GetInitialUtility()) {
        LOG(WARNING) << "Bad InitialUtility copying!";
    }
    if (Agent1->GetInitialWealth() != Agent2->GetInitialWealth()) {
        LOG(WARNING) << "Bad InitialWealth copying!";
    }
}   //  AgentPopulation::CompareTwoAgents()

void AgentPopulation::ShockAgentPreferences() {
    LOG(DEBUG) << "Shocking agent preferences...";
    AgentPtr ActiveAgent;
    double oldPref, newPref, pref;

    size_t CommodityToShock = randomCommodity(rng);
    bool sign = randomBinary(rng);
    double shock = randomShock(rng);

    if (sign) {
        LOG(DEBUG) << "Shocking commodity " << CommodityToShock << " * " << shock;
    } else {
        LOG(DEBUG) << "Shocking commodity " << CommodityToShock << " * 1/" << shock;
    }
    for (size_t AgentIndex = 1; AgentIndex <= static_cast<size_t>(NumberOfAgents); ++AgentIndex) {
        ActiveAgent = Agents[AgentIndex];
        oldPref = ActiveAgent->GetAlpha(CommodityToShock);
        if (sign) {
            newPref = oldPref * shock;
        } else {
            newPref = oldPref/shock;
        }
        ActiveAgent->SetAlpha(CommodityToShock, newPref);

        for (size_t CommodityIndex = 1; CommodityIndex <= static_cast<size_t>(NumberOfCommodities); ++CommodityIndex) {
            pref = ActiveAgent->GetAlpha(CommodityIndex);
            ActiveAgent->SetAlpha(CommodityIndex, pref/(1.0 - oldPref + newPref));
        }   //  for (CommodityIndex...
    }   //  for (AgentIndex...
}   //  AgentPopulation::ShockAgentPreferenes()

/*================= End of Methods =================*/

void InitMiscellaneous() {
    LOG(INFO) << "Model version: " << Version;
    LOG(INFO) << "Number of agents: " << NumberOfAgents;
    LOG(INFO) << "Number of commodities: " << NumberOfCommodities;
    LOG(INFO) << "Trade epsilon: " << trade_eps;
    LOG(INFO) << "Time period: " << 2 * PairwiseInteractionsPerPeriod;
    LOG(INFO) << "Termination period check rate: " << CheckTerminationPeriod;
    LOG(INFO) << "Termination period check threshold: " << CheckTerminationThreshold;
    switch (termination_criterion) {
        case -1:
        LOG(INFO) << "Termination criterion: L2 norm of standard deviation of agent MRSes";
        break;
        case 0:
        LOG(INFO) << "Termination criterion: Maximum standard deviation of agent MRSes";
        break;
        case 1:
        LOG(INFO) << "Termination criterion: relative increase in sum of agent utilities";
        break;
        case 2:
        LOG(INFO) << "Termination criterion: increase in the sum of agent utilities";
        break;
    }
    LOG(INFO) << "Termination threshold: " << termination_eps;
    if (DefaultSerialExecution) {
        LOG(INFO) << "Parallel activation: FALSE";
    } else {
        LOG(INFO) << "Parallel activation: TRUE";
        LOG(INFO) << "Agent randomization size: " << AgentsToRandomize;
    }
    LOG(INFO) << "Number of equilibrations: " << RequestedEquilibrations;
    LOG(INFO) << "Vary agent initial conditions: " << !SameAgentInitialCondition;

    MemoryState.WriteMemoryRequirements();
} // InitMiscellaneous()

void SeedRNG() {
// Seed the random number generator.
    if (!UseRandomSeed) {
        random_device rd;
        unsigned int seed = rd();
        rng.seed(seed);
        LOG(INFO) << "Using random seed " << seed;
    } else {
        LOG(INFO) << "Using fixed seed " << NonRandomSeed;
        rng.seed(NonRandomSeed);
    }
    // initialize the distributions now that we know the relevant ranges
    randomAgent = uniform_int_distribution<unsigned long>(1, static_cast<size_t>(NumberOfAgents));  //  why are we 1-indexing?
    randomCommodity = uniform_int_distribution<unsigned long>(1, static_cast<size_t>(NumberOfCommodities));
    randomBinary = uniform_int_distribution<int>(0, 1);
    randomShock = uniform_real_distribution<double>(MinShock, MaxShock);
    randomAlpha = uniform_real_distribution<double>(alphaMin, alphaMax);
    randomWealth = uniform_real_distribution<double>(wealthMin, wealthMax);
//  uniform_int_distribution<double> randomDouble(0, 1);
} // SeedRNG()

void ReadConfigFile(string file) {
    // This function sets all relevant model parameters by reading from a config file (libconfig). 
    // If there's some issue with the formatting or reading of the config file, it catches the exception
    // and terminates the program. The config file *must* be proper for the model to run. See parameters.cfg
    // in the repository for an example.

    libconfig::Config config;
    try { 
        config.readFile(file.c_str());
        LOG(INFO) << "Loaded configuration from " << file;
        if (config.lookupValue("debug.enabled", debug) && debug) { LOG(DEBUG) << "debug: " << debug; }
        if (config.lookupValue("number.commodities", NumberOfCommodities) && debug) { LOG(DEBUG) << "NumberOfCommodities: " << NumberOfCommodities; }
        if (config.lookupValue("number.agents", NumberOfAgents) && debug) { LOG(DEBUG) << "NumberOfAgents: " << NumberOfAgents; }
        if (config.lookupValue ("rand.use_seed", UseRandomSeed) && debug) { LOG(DEBUG) << "UseRandomSeed: " << UseRandomSeed; }
        if (config.lookupValue("rand.seed", NonRandomSeed) && debug) { LOG(DEBUG) << "NonRandomSeed: " << NonRandomSeed; }
        if (config.lookupValue("interactions_per_period", PairwiseInteractionsPerPeriod) && debug) { LOG(DEBUG) << "PairwiseInteractionsPerPeriod: " << PairwiseInteractionsPerPeriod; }
        if (config.lookupValue("alpha.min", alphaMin) && debug) { LOG(DEBUG) << "alphaMin: " << alphaMin; }
        if (config.lookupValue("alpha.max", alphaMax) && debug) { LOG(DEBUG) << "alphaMax: " << alphaMax; }
        if (config.lookupValue("wealth.min", wealthMin) && debug) { LOG(DEBUG) << "wealthMin: " << wealthMin; }
        if (config.lookupValue("wealth.max", wealthMax) && debug) { LOG(DEBUG) << "wealthMax: " << wealthMax; }
        if (config.lookupValue("parallel.disabled", DefaultSerialExecution) && debug) { LOG(DEBUG) << "DefaultSerialExecution: " << DefaultSerialExecution; }
        if (config.lookupValue("parallel.agents_to_randomize", AgentsToRandomize) && debug) { LOG(DEBUG) << "AgentsToRandomize: " << AgentsToRandomize; }
        if (config.lookupValue("parallel.randomization_method", RandomizationMethod) && debug) { LOG(DEBUG) << "RandomizationMethod: " << RandomizationMethod; }
        if (config.lookupValue("num_equilibrations", RequestedEquilibrations) && debug) { LOG(DEBUG) << "RequestedEquilibrations: " << RequestedEquilibrations; }
        if (config.lookupValue("same_conditions_each_equilibration", SameAgentInitialCondition) && debug) { LOG(DEBUG) << "SameAgentInitialCondition: " << SameAgentInitialCondition; }
        if (config.lookupValue("trade.eps", trade_eps) && debug) { LOG(DEBUG) << "trade_eps: " << trade_eps; }
        if (config.lookupValue("termination.criterion", termination_criterion) && debug) { LOG(DEBUG) << "termination_criterion: " << termination_criterion; }
        if (config.lookupValue("termination.time", TerminationTime) && debug) { LOG(DEBUG) << "TerminationTime: " << TerminationTime; }
        if (config.lookupValue("termination.eps", termination_eps) && debug) { LOG(DEBUG) << "termination_eps: " << termination_eps; }
        if (config.lookupValue("termination.threshold", CheckTerminationThreshold) && debug) { LOG(DEBUG) << "CheckTerminationThreshold: " << CheckTerminationThreshold; }
        if (config.lookupValue("termination.period", CheckTerminationPeriod) && debug) { LOG(DEBUG) << "CheckTerminationPeriod: " << CheckTerminationPeriod; }
        if (config.lookupValue("shock.enabled", ShockPreferences) && debug) { LOG(DEBUG) << "ShockPreferences: " << ShockPreferences; }
        if (config.lookupValue("shock.period", ShockPeriod) && debug) { LOG(DEBUG) << "ShockPeriod: " << ShockPeriod; }
        if (config.lookupValue("shock.min", MinShock) && debug) { LOG(DEBUG) << "MinShock: " << MinShock; }
        if (config.lookupValue("shock.max", MaxShock) && debug) { LOG(DEBUG) << "MaxShock: " << MaxShock; }
        if (config.lookupValue("debug.print_endowments", PrintEndowments) && debug) { LOG(DEBUG) << "PrintEndowments: " << PrintEndowments; }
        if (config.lookupValue("debug.print_intermediate_output", PrintIntermediateOutput) && debug) { LOG(DEBUG) << "PrintIntermediateOutput: " << PrintIntermediateOutput; }
        if (config.lookupValue("debug.intermediate_output_print_period", IntermediateOutputPrintPeriod) && debug) { LOG(DEBUG) << "IntermediateOutputPrintPeriod: " << IntermediateOutputPrintPeriod; }
        if (config.lookupValue("debug.print_convergence_stats", PrintConvergenceStats) && debug) { LOG(DEBUG) << "PrintConvergenceStats: " << PrintConvergenceStats; }
        if (config.lookupValue("debug.print_final_commodity_list", PrintFinalCommodityList) && debug) { LOG(DEBUG) << "PrintFinalCommodityList: " << PrintFinalCommodityList; }
        // NumberOfAgents = config.lookup("number.agents");
        // NumberOfCommodities = config.lookup("number.commodities");
        exp_trade_eps = exp(trade_eps);

    } catch (...) {
        LOG(ERROR) << "Error reading config file";
        terminate();
    }
} //ReadConfigFile

int main(int argc, char** argv) {
    // Preliminaries: Parse flags, etc.
    string usage = "An agent-based model of bilateral exchange. Usage:\n";
    usage += argv[0];
    gflags::SetUsageMessage(usage);
    gflags::ParseCommandLineFlags(&argc, &argv, true);
    LOG(INFO) << "Opening log file...";

    // Read the config file passed through the -file flag, or read the default parameters.cfg.
    ReadConfigFile(FLAGS_file);

    SeedRNG();

    //  Initialize the model and print preliminaries to the log.
    InitMiscellaneous();
    AgentPopulationPtr PopulationPtr = new AgentPopulation(NumberOfCommodities);

    //  Equilibrate the agent economy once...
    int EquilibrationNumber = 1;
    long long sum = PopulationPtr->Equilibrate(EquilibrationNumber);
    long long sum2 = sum*sum;

    //  Equilibrate again if the user has requested this...
    long long interactions;

    EquilibrationNumber = 2;
    while (EquilibrationNumber <= RequestedEquilibrations) {
        if (SameAgentInitialCondition) {
            PopulationPtr->Reset();
        } else {
            PopulationPtr->Init();
        }
        interactions = PopulationPtr->Equilibrate(EquilibrationNumber);
        sum += interactions;
        sum2 += interactions*interactions;
        ++EquilibrationNumber;
    }

    double avg = static_cast<double>(sum)/(EquilibrationNumber-1);
    LOG(INFO) << "Average number of interactions: " << avg;
    if (EquilibrationNumber > 2) {
        LOG(INFO) << "std. dev.: " << sqrt((sum2 - (EquilibrationNumber - 1) * avg * avg)/(EquilibrationNumber - 2));
    }
} // main()
