// Copyright 2015 <Stefan McCabe>
/*
Bilateral Exchange Model
Originally written by Rob Axtell
Extended by Stefan McCabe

12/11/2015
*/
#ifndef MAIN_H_
#define MAIN_H_
#endif

#include <gflags/gflags.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
#include "./easylogging++.h"
#define ELPP_THREAD_SAFE
#include <libconfig.h++>

//using namespace std;

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


// Random number distributions. This will probably be replaced with some sort
// of thread-safe random-number generator.
std::uniform_int_distribution<unsigned long> randomAgent;  
std::uniform_int_distribution<unsigned long> randomCommodity;
std::uniform_int_distribution<int> randomBinary;
std::uniform_real_distribution<double> randomShock;
std::uniform_real_distribution<double> randomAlpha;
std::uniform_real_distribution<double> randomWealth;


typedef std::vector<double> CommodityArray;
typedef CommodityArray *CommodityArrayPtr;

double inline Dot(CommodityArrayPtr vector1, CommodityArrayPtr vector2);

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

typedef Agent *AgentPtr;

class AgentPopulation {
    std::vector<AgentPtr> Agents;
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
    AgentPopulation(int size); 
    void Init();
    void Reset();
    long long Equilibrate(int NumberOfEquilibrationsSoFar);
    void ConvergenceStatistics(CommodityArray VolumeStats);
    void CompareTwoAgents(AgentPtr Agent1, AgentPtr Agent2);
    void ShockAgentPreferences();
};

typedef AgentPopulation *AgentPopulationPtr;

