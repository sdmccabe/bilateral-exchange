// Copyright 2015 <Stefan McCabe>
/*
Bilateral Exchange Model
Originally written by Rob Axtell
Extended by Stefan McCabe

12/11/2015
*/

#ifndef MAIN_H_
#define MAIN_H_
#endif  // MAIN_H_

#include "./RNG.h"
#include "tbb/tbb.h"
#include "tbb/concurrent_vector.h"
#include <mutex>


// Global variables specifying model parameters. See parameters.cfg for documentation.
bool UseRandomSeed;
unsigned int NonRandomSeed;
double Version = 1.0;
unsigned int NumberOfAgents;
unsigned int NumberOfCommodities;
unsigned int PairwiseInteractionsPerPeriod;
double alphaMin;
double alphaMax;
unsigned int wealthMin;
unsigned int wealthMax;
bool DefaultSerialExecution;
int AgentsToRandomize;
int RandomizationMethod;
int RequestedEquilibrations;
bool SameAgentInitialCondition;
double trade_eps;
double exp_trade_eps;
int termination_criterion;
long long TerminationTime;
double termination_eps;
long long CheckTerminationThreshold;
int CheckTerminationPeriod;
bool ShockPreferences;
int ShockPeriod;
double MinShock;
double MaxShock;
bool debug;
bool PrintEndowments;
bool PrintIntermediateOutput;
int IntermediateOutputPrintPeriod;
bool PrintConvergenceStats;
bool PrintFinalCommodityList;
int activationMethod;
const char* outputFilename = NULL;
bool fileAppend;
bool writeToFile;
std::ofstream outfile;
int NumberOfThreads;
bool ForkAndJoin;
bool ShuffleAfterJoin;
bool DumpAgentInformation;

RNGptr Rand;


typedef tbb::concurrent_vector<double> CommodityArray;
typedef CommodityArray *CommodityArrayPtr;

// Functions
double inline Dot(CommodityArrayPtr vector1, CommodityArrayPtr vector2);
void ReadConfigFile(std::string file);
void InitMiscellaneous();
void OpenFile(const char * filename);
void WriteHeader();
void WriteLine();

// Classes and methods
class MemoryObject {
    long long start;

 public:
    MemoryObject();
    void WriteMemoryRequirements();
} MemoryState;


class Data {
    int N;
    double min;
    double max;
    double sum;
    double sum2;
 public:
    Data();
    void Init();
    void AddDatum(double Datum);
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
};

typedef Data *DataPtr;

class CommodityData {
    tbb::concurrent_vector<Data> data;

 public:
    CommodityData();
    void Clear();
    DataPtr GetData(size_t index) { return &data[index]; }
    double L2StdDev();
    double LinfStdDev();
};

class Agent {
    size_t id;
    
    CommodityArray alphas;
    CommodityArray endowment;
    CommodityArray initialMRSs;
    CommodityArray allocation;
    CommodityArray currentMRSs;
    
    Agent();
    double initialUtility;
    double initialWealth;
    double lambda;  // for Poisson activation
    double nextTime;  // for Poisson activation
    double initiative; // exogenous Poisson variable, distributed U(0,1)
    long long interactions = 0;
    long long activations = 0;
 public:
    explicit Agent(int size, size_t x);
    void Init();
    void Reset();
    
    
    std::mutex m;
    void MarkActivated() {
        //std::lock_guard<std::mutex> lock(m);
        activations++;
    }
    long long GetNumberOfActivations() {
        return activations;
    }
    void MarkSuccessfulTrade() {
        interactions++;
    }
    long long GetNumberOfTrades() {
        return interactions;
    }
    size_t GetId() {
        return id;
    }
    double GetAlpha(size_t CommodityIndex) {
        //std::lock_guard<std::mutex> lock(m);
        return alphas[CommodityIndex];
    }
    void SetAlpha(size_t CommodityIndex, double alpha) {
        //std::lock_guard<std::mutex> lock(m);
        alphas[CommodityIndex] = alpha;
    }
    void SetLambda(double lam) {
        //std::lock_guard<std::mutex> lock(m);
        lambda = lam;
    }
    double GetLambda() {
        //std::lock_guard<std::mutex> lock(m);
        return lambda;
    }
    void SetNextTime(double nextT) {
        //std::lock_guard<std::mutex> lock(m);
        nextTime = nextT;
    }
    double GetNextTime() {
        //std::lock_guard<std::mutex> lock(m);
        return nextTime;
    }
    double GetEndowment(size_t CommodityIndex) {
        //std::lock_guard<std::mutex> lock(m);
        return endowment[CommodityIndex];
    }
    double MRS(size_t CommodityIndex, size_t Numeraire);
    void ComputeMRSs();
    double GetInitialMRS(size_t CommodityIndex) {
        //std::lock_guard<std::mutex> lock(m);
        return initialMRSs[CommodityIndex];
    }
    double GetAllocation(size_t CommodityIndex) {
        //std::lock_guard<std::mutex> lock(m);
        return allocation[CommodityIndex];
    }
    void IncreaseAllocation(size_t CommodityIndex, double amount) {
        //std::lock_guard<std::mutex> lock(m);
        allocation[CommodityIndex] += amount;
    }
    double GetCurrentMRS(size_t CommodityIndex) {
        //std::lock_guard<std::mutex> lock(m);
        return currentMRSs[CommodityIndex];
    }
    CommodityArrayPtr GetCurrentMRSs() {
        //std::lock_guard<std::mutex> lock(m);
        return &currentMRSs;
    }
    double Utility();
    double GetInitialUtility() {
        //std::lock_guard<std::mutex> lock(m);
        return initialUtility;
    }
    double Wealth(CommodityArrayPtr prices) {
        //std::lock_guard<std::mutex> lock(m);
        return Dot(&allocation, prices);
    }
    double GetInitialWealth() {
        //std::lock_guard<std::mutex> lock(m);
        return initialWealth;
    }
    double GetInitiative() {
        return initiative;
    }
};

typedef std::shared_ptr<Agent> AgentPtr;

class AgentPopulation {
    //tbb::concurrent_vector<AgentPtr> Agents;
    tbb::concurrent_vector<AgentPtr> Agents;
    tbb::concurrent_vector<size_t> AgentIndices;
    tbb::concurrent_vector<std::pair<double,AgentPtr>> PoissonActivations;
    size_t AgentIndex = 0;
    bool PoissonUpToDate;
    Data InitialOwnWealthData;
    CommodityArray Volume;
    CommodityData AlphaData, EndowmentData, LnMRSsData;
    bool LnMRSsDataUpToDate;
    tbb::concurrent_vector<std::tuple<long long, double, size_t, double, size_t>> results;
    void ComputeLnMRSsDistribution();
    double LastSumOfUtilities;
    double ComputeSumOfUtilities();
    double ComputeIncreaseInSumOfUtilities() {
        return ComputeSumOfUtilities() - LastSumOfUtilities;
    }
    double ComputeRelativeIncreaseInSumOfUtilities() {
        return ComputeIncreaseInSumOfUtilities() / LastSumOfUtilities;
    }

    void GetRandomAgentPair(AgentPtr& Agent1, AgentPtr& Agent2);
    void GetUniformAgentPair(AgentPtr& Agent1, AgentPtr& Agent2);
    void GetFixedAgentPair(AgentPtr& Agent1, AgentPtr& Agent2);
    void GetPoissonAgentPair(AgentPtr& Agent1, AgentPtr& Agent2);
    void SetPoissonAgentDistribution();
    void Trade(AgentPtr a1, AgentPtr a2);

    void(AgentPopulation::*GetAgentPair) (AgentPtr& Agent1, AgentPtr& Agent2);
    
    std::tuple<long long, double, size_t, double, size_t> ParallelTrade (AgentPtr a1, AgentPtr a2);
    void TradeInFork (tbb::concurrent_vector<AgentPtr> a);
    
    void(AgentPopulation::*GetAgentPairInFork) (AgentPtr& Agent1, AgentPtr& Agent2, tbb::concurrent_vector<AgentPtr> a, tbb::concurrent_vector<AgentPtr>::iterator &ait);
    void GetRandomAgentPairInFork(AgentPtr& Agent1, AgentPtr& Agent2, tbb::concurrent_vector<AgentPtr> a, tbb::concurrent_vector<AgentPtr>::iterator &ait);
    void GetUniformAgentPairInFork(AgentPtr& Agent1, AgentPtr& Agent2, tbb::concurrent_vector<AgentPtr> a, tbb::concurrent_vector<AgentPtr>::iterator &ait);
    void GetFixedAgentPairInFork(AgentPtr& Agent1, AgentPtr& Agent2, tbb::concurrent_vector<AgentPtr> a, tbb::concurrent_vector<AgentPtr>::iterator &ait);
    void GetPoissonAgentPairInFork(AgentPtr& Agent1, AgentPtr& Agent2, tbb::concurrent_vector<AgentPtr> a, tbb::concurrent_vector<AgentPtr>::iterator &ait);
    void SetPoissonAgentDistributionInFork(tbb::concurrent_vector<AgentPtr> a, tbb::concurrent_vector<AgentPtr>::iterator &ait);

    void IntermediateOutput();
    bool TestConvergence();
    void DumpAgentInfo();

 public:
    explicit AgentPopulation(int size);
    bool Converged;
    long long theTime; 
    long long TotalInteractions;

    void IncreaseTotalInteractions(long long x) {
        TotalInteractions += x;
    }
    void Init();
    void Reset();
    long long Equilibrate(int NumberOfEquilibrationsSoFar);
    long long ParallelEquilibrate(int NumberOfEquilibrationsSoFar);
    long long ForkAndJoinEquilibrate(int NumberOfEquilibrationsSoFar);
    void ConvergenceStatistics(CommodityArray VolumeStats);
    void CompareTwoAgents(AgentPtr Agent1, AgentPtr Agent2);
    void ShockAgentPreferences();
    std::string WriteWealthInfo();
    std::string WriteUtilityInfo();
    void WriteLine();
    std::mutex m;
};

typedef AgentPopulation *AgentPopulationPtr;
