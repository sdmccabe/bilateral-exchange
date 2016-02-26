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
#include <memory>
#include <random>
#include <string>
#include <vector>
#include <fstream>

// Global variables specifying model parameters. See parameters.cfg for documentation.
bool UseRandomSeed;
unsigned int NonRandomSeed;
double Version = 1.0;
int NumberOfAgents;
int NumberOfCommodities;
int PairwiseInteractionsPerPeriod;
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

RNGptr Rand;


typedef std::vector<double> CommodityArray;
typedef CommodityArray *CommodityArrayPtr;

// Functions
double inline Dot(CommodityArrayPtr vector1, CommodityArrayPtr vector2);
void ReadConfigFile(std::string file);
void InitMiscellaneous();
//void SeedRNG();
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
    std::vector<Data> data;
    
public:
    CommodityData();
    void Clear();
    DataPtr GetData(size_t index) { return &data[index]; }
    double L2StdDev();
    double LinfStdDev();
};

class Agent {
    CommodityArray alphas;
    CommodityArray endowment;
    CommodityArray initialMRSs;
    CommodityArray allocation;
    CommodityArray currentMRSs;
    Agent();
    double initialUtility;
    double initialWealth;
    double lambda; // for Poisson activation
    double nextTime; // for Poisson activation

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
    void SetLambda(double lam) {
        lambda = lam;
    }
    double GetLambda() {
        return lambda;
    }
    void SetNextTime(double nextT) {
        nextTime = nextT;
    }
    double GetNextTime() {
        return nextTime;
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




//typedef Agent *AgentPtr;
typedef std::shared_ptr<Agent> AgentPtr;

class AgentPopulation {
    std::vector<AgentPtr> Agents;
    std::vector<size_t> AgentIndices;
    std::vector<std::pair<double,AgentPtr>> PoissonActivations;
    size_t AgentIndex = 0;
    bool PoissonUpToDate;
    Data InitialOwnWealthData;
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

    void GetRandomAgentPair(AgentPtr& Agent1, AgentPtr& Agent2);
    void GetUniformAgentPair(AgentPtr& Agent1, AgentPtr& Agent2);
    void GetFixedAgentPair(AgentPtr& Agent1, AgentPtr& Agent2);
    void GetPoissonAgentPair(AgentPtr& Agent1, AgentPtr& Agent2);
    void SetPoissonAgentDistribution();
    void TradeInParallel(AgentPtr a1, AgentPtr a2);
    

    void(AgentPopulation::*GetAgentPair) (AgentPtr& Agent1, AgentPtr& Agent2);

public:
    AgentPopulation(int size); 
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
    void ConvergenceStatistics(CommodityArray VolumeStats);
    void CompareTwoAgents(AgentPtr Agent1, AgentPtr Agent2);
    void ShockAgentPreferences();
    std::string WriteWealthInfo();
    std::string WriteUtilityInfo();
    void WriteLine();
};

typedef AgentPopulation *AgentPopulationPtr;
