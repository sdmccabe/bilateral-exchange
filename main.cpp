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

#include "./main.h"
#include <assert.h>
#include <gflags/gflags.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <libconfig.h++>
#include <random>
#include <string>
#include <vector>
#define ELPP_THREAD_SAFE
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wundef"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wctor-dtor-privacy"
#include "./easylogging++.h"
#pragma GCC diagnostic pop


//using namespace std;

// Initialize the logger. This should come immediately after the #includes are finished.
INITIALIZE_EASYLOGGINGPP

// Declare gflags.
DEFINE_string(file, "parameters.cfg", "Specify the configuation file containing the model parameters.");

// addressing the RNG first. Seeded to 0 for initialization, 
// will be properly seeded in main().
//std::mt19937 rng(0);
std::knuth_b rng(0);

/*===========
    ==Methods==
    ===========*/
    double inline Dot(CommodityArrayPtr vector1, CommodityArrayPtr vector2) {
        double sum = 0.0;
        for (size_t i = 0; i < static_cast<size_t>(NumberOfCommodities); ++i) {
            sum += (*vector1)[i] * (*vector2)[i];
        }
        return sum;
    }

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
        + sizeof(AgentPopulation) + sizeof(Agent) * static_cast<unsigned long>(NumberOfAgents)) << " bytes";
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
    // Constructor resizes the data vector to the appropriate size.
    // It is initially zero-initialized, which is obviously problematic.
    data.resize(static_cast<size_t>(NumberOfCommodities));
}

void CommodityData::Clear() {
    //  Initialize data objects
    //
    for (auto& commodity : data) {
        commodity.Init();
    }
}   //  CommodityData::Clear()

double CommodityData::L2StdDev() {
    double sum = 0.0;

    // Instead of computing the standard deviation for each commodity (and calling sqrt() N times),
    // find the largest variance and from it get the std dev
    for (auto& commodity : data) {
        sum += commodity.GetVariance();
    }
    return sqrt(sum);
}   //  CommodityData::L2StdDev

double  CommodityData::LinfStdDev() {
    double var, max = 0.0;

    // Instead of computing the standard deviation for each commodity (and calling sqrt() N times),
    // find the largest variance and from it get the std dev
    for (auto& commodity : data) {
        var = commodity.GetVariance();
        if (var > max) {
            max = var;
        }
    }
    return sqrt(max);
}   //  CommodityData::LinfStdDev


Agent::Agent(int size): initialUtility(0.0), initialWealth(0.0) {
    // Constructor resizes the data vector to the appropriate size.
    // It is initially zero-initialized, which is obviously problematic.
    // This allows me to replace arrays with vectors.

    alphas.resize(static_cast<size_t>(size));
    endowment.resize(static_cast<size_t>(size));   
    initialMRSs.resize(static_cast<size_t>(size)); 
    allocation.resize(static_cast<size_t>(size));  
    currentMRSs.resize(static_cast<size_t>(size));
    Init();
}

void Agent::Init() {
    size_t CommodityIndex;

    //  First generate and normalize the exponents...
    //
    double sum = 0.0;
    for (auto& alpha : alphas) {
        alpha = randomAlpha(rng);
        sum += alpha;
    }
    //  Next, fill up the rest of the agent fields...
    //
    for (auto& alpha : alphas) {
        auto i = static_cast<size_t>(&alpha - &alphas[0]);
        alpha = alpha / sum;
        endowment[i] = randomWealth(rng);
        allocation[i] = endowment[i];
        initialMRSs[i] = MRS(i, 1);
        currentMRSs[i] = initialMRSs[i];
    }
    initialUtility = Utility();
    initialWealth = Wealth(&initialMRSs);
}   //  Agent:Init()

void Agent::Reset() {
    for (auto& currentAgentMRS : currentMRSs) {
        auto i = static_cast<size_t>(&currentAgentMRS - &currentMRSs[0]);
        allocation[i] = endowment[i];
        currentAgentMRS = initialMRSs[i];
    }
}   //  Agent::Reset

// MRS = marginal rate of substitution
double Agent::MRS(size_t CommodityIndex, size_t Numeraire) {
    return (alphas[CommodityIndex] * allocation[Numeraire]) / (alphas[Numeraire] * allocation[CommodityIndex]);
}   //  Agent::MRS()

void Agent::ComputeMRSs() {
    for (auto& currentAgentMRS : currentMRSs) {
        auto i = static_cast<size_t>(&currentAgentMRS - &currentMRSs[0]);
        if (i == 0) {
            currentAgentMRS = 1.0;
        } else {
            currentMRSs[i] = MRS(i, 0);
        }
    }
}   //      Agent::ComputeMRSs()

double Agent::Utility() {
    double product = 1.0;

    for (size_t i = 0; i < allocation.size(); ++i) {
        product *= pow(allocation[i], alphas[i]);
    }

    return product;
}   //  Agent::Utility()

void AgentPopulation::ComputeLnMRSsDistribution() {
    LnMRSsData.Clear();

    for (auto& agent : Agents) {
        agent->ComputeMRSs();
        for (size_t j = 0; j < static_cast<size_t>(NumberOfCommodities); ++j) {
            LnMRSsData.GetData(j)->AddDatuum(log(agent->GetCurrentMRS(j)));
        }
    }       //  for i...
    LnMRSsDataUpToDate = true;
}   //  AgentPopulation::ComputeLnMRSsDistribution()

double AgentPopulation::ComputeSumOfUtilities() {
    double sum = 0.0;

    for (auto& agent : Agents) {
        sum += agent->Utility();
    }
    return sum;
}   //  AgentPopulation::ComputeSumOfUtilities()

void AgentPopulation::GetRandomAgentPair(AgentPtr& Agent1, AgentPtr& Agent2) {
    Agent1 = Agents[randomAgent(rng)];
    do {
        Agent2 = Agents[randomAgent(rng)];
    } while (Agent2 == Agent1);
}   //  AgentPopulation::GetRandomAgentPair()

void AgentPopulation::GetUniformAgentPair(AgentPtr& Agent1, AgentPtr& Agent2) {
    if (NumberOfAgents % 2 > 0 && AgentIndex == 0) {
        LOG(DEBUG) << "Warning: Uniform activation requires an even number of agents.";
    }
    Agent1 = Agents[AgentIndices[AgentIndex++]];
    Agent2 = Agents[AgentIndices[AgentIndex++]];
    if (AgentIndex >= static_cast<size_t>(NumberOfAgents)) {
        //LOG(DEBUG) << "Rolling over uniform indices...";
        AgentIndex = 0;
        std::shuffle(AgentIndices.begin(), AgentIndices.end(), rng);
    }
}

void AgentPopulation::GetFixedAgentPair(AgentPtr& Agent1, AgentPtr& Agent2) {
    if (NumberOfAgents % 2 > 0 && AgentIndex == 0) {
        LOG(DEBUG) << "Warning: Fixed activation requires an even number of agents.";
    }
    Agent1 = Agents[AgentIndices[AgentIndex++]];
    Agent2 = Agents[AgentIndices[AgentIndex++]];
    if (AgentIndex >= static_cast<size_t>(NumberOfAgents)) {
        //LOG(DEBUG) << "Rolling over uniform indices...";
        AgentIndex = 0;
    }
}

void AgentPopulation::GetPoissonAgentPair(AgentPtr& Agent1, AgentPtr& Agent2) {
    if (!PoissonUpToDate) {
        SetPoissonAgentDistribution();
    }
    size_t AgentIndex1 = AgentIndex++;
    size_t AgentIndex2 = AgentIndex;

    Agent1 = PoissonActivations[AgentIndex1].second;
    Agent2 = PoissonActivations[AgentIndex2].second;

    size_t swap = 1;
    while (Agent1 == Agent2) {
        std::swap(PoissonActivations[AgentIndex2], PoissonActivations[AgentIndex2+swap]);
        ++swap;
        Agent2 = PoissonActivations[AgentIndex2].second;
    }
    
    ++AgentIndex;
    if (AgentIndex >= static_cast<size_t>(NumberOfAgents)) {
        PoissonUpToDate = false;
    }
}

//TODO: Verify all of this.
void AgentPopulation::SetPoissonAgentDistribution() {
    // reset data structures
    PoissonActivations.clear();
    PoissonActivations.shrink_to_fit(); // memory leaks are bad
    std::vector<std::pair<double,AgentPtr>>().swap(PoissonActivations); 
    // more ritual to stave off memory leaks
    AgentIndex = 0;

    double totalWealth = 0.0;
    double totalDistanceFromMean = 0.0;
    double denom;
    double totalLambda = 0.0;

    //determine mean wealth
    for (auto &a : Agents) {
        totalWealth += a->Wealth(a->GetCurrentMRSs());
    }
    double meanWealth = totalWealth / static_cast<double>(NumberOfAgents);

    //determine total distance from the mean
    for (auto &a : Agents) {
        double distFromMean = std::abs(a->Wealth(a->GetCurrentMRSs()) - meanWealth);
        totalDistanceFromMean += distFromMean;
    }

    //update lambdas based on distance from the mean
    //this is where the Poisson activation methods are differentiated
    //TODO: "regular Poisson" - inverse distance from mean - don't really care about it
    for (auto &a : Agents) {
            if (activationMethod == 2) { // furthest from mean activate faster
                denom = std::abs(a->Wealth(a->GetCurrentMRSs()) - meanWealth);
                if (denom == 0) {
                    denom = 0.0001;
                }
                double lam = totalDistanceFromMean / denom;
                a->SetLambda(lam);
                totalLambda += lam;
            } else if (activationMethod == 3) { // poor activate faster
                denom = a->Wealth(a->GetCurrentMRSs());
                if (denom == 0) {
                    denom = 0.0001;
                }
                double lam = 1 / denom;
                a->SetLambda(lam);
                totalLambda += lam;
        } else if (activationMethod == 4) { // rich activate faster
            denom = a->Wealth(a->GetCurrentMRSs());
            if (denom == 0) {
                denom = 0.0001;
            }
            double lam = denom;
                //std::cout << lam << " ";
            a->SetLambda(lam);
            totalLambda += lam;
        } else if (activationMethod == 5) { // closest to mean activate faster
            double lam = std::abs(a->Wealth(a->GetCurrentMRSs()) - meanWealth);
            a->SetLambda(lam);
            totalLambda += lam;
        } 
    }

    //normalize lambdas
    for (auto &a : Agents) {
        auto i = static_cast<size_t>(&a - &Agents[0]);
        double lam = a->GetLambda() * static_cast<double>(NumberOfAgents) * 1.1/totalLambda;
        if (lam == 0) {
            lam = 1.0 / static_cast<double>(NumberOfAgents);
        }
        //done normalizing lambdas

        //schedule agents based on lambda
        a->SetLambda(lam);
        a->SetNextTime(-1 * log(randomDouble(rng)) / a->GetLambda());
        
        while (a->GetNextTime() < 1 ) { //no agent is more than 20% of activations
            PoissonActivations.push_back(std::make_pair(a->GetNextTime(), a));
            a->SetNextTime(a->GetNextTime() + -1 * log(randomDouble(rng)) / a->GetLambda());
        }     
    }

    //sort the activations vector
    std::sort(PoissonActivations.begin(), PoissonActivations.end(), [](const std::pair<double, AgentPtr> &left, const std::pair<double, AgentPtr> &right) {
        return left.first < right.first;
    });

    PoissonUpToDate = true;
}

AgentPopulation::AgentPopulation(int size):
AlphaData(), EndowmentData(), LnMRSsData(), LnMRSsDataUpToDate(true), LastSumOfUtilities(0.0), GetAgentPair(NULL) {
    Volume.resize(static_cast<size_t>(size));
    //Agents.resize(static_cast<size_t>(NumberOfAgents));
    AgentIndices.resize(static_cast<size_t>(NumberOfAgents));
    for (size_t i = 0; i < AgentIndices.size(); ++i) {
        AgentIndices[i] = i;
    }
    
    if (activationMethod == 1) { std::shuffle(AgentIndices.begin(), AgentIndices.end(), rng); }
    
    AgentIndex = 0;
    // for (auto& agent : Agents) {
    //     agent = new Agent(NumberOfCommodities);
    // }
    for (size_t i = 0; i < static_cast<size_t>(NumberOfAgents); ++i) {
        Agents.emplace_back(new Agent{NumberOfCommodities});
    }

    PoissonUpToDate = false;
    switch (activationMethod) {
        case -1:
        GetAgentPair = &AgentPopulation::GetFixedAgentPair;
        LOG(INFO) << "Using fixed activation";
        LOG(INFO) << "WARNING: Do not use fixed activation.";
        break;
        case 0:
        GetAgentPair = &AgentPopulation::GetRandomAgentPair;
        LOG(INFO) << "Using random activation";
        break;
        case 1:
        GetAgentPair = &AgentPopulation::GetUniformAgentPair;
        LOG(INFO) << "Using uniform activation";
        break;
        case 2:
        GetAgentPair = &AgentPopulation::GetPoissonAgentPair;
        LOG(INFO) << "Using Poisson activation (λ = 1/|wealth - mean(wealth)|)";
        break;
        case 3:
        GetAgentPair = &AgentPopulation::GetPoissonAgentPair;
        LOG(INFO) << "Using Poisson activation (λ = 1/wealth)";
        break;
        case 4:
        GetAgentPair = &AgentPopulation::GetPoissonAgentPair;
        LOG(INFO) << "Using Poisson activation (λ = wealth)";
        break;
        case 5:
        GetAgentPair = &AgentPopulation::GetPoissonAgentPair;
        LOG(INFO) << "Using Poisson activation (λ = |wealth - mean(wealth)|)";
        break;
        default:
        LOG(ERROR) << "Invalid activation method (or NYI)";
        std::terminate();
        break;
    }
}   // Constructor...

void AgentPopulation::Init() {
    size_t CommodityIndex;
    Data InitialOwnWealthData;

    
    AlphaData.Clear();
    EndowmentData.Clear();
    LnMRSsData.Clear();

    //  Cycle through the agents...
    //
    for (auto& ActiveAgent : Agents) {
        ActiveAgent->Init();

        //  Next, fill up the rest of the agent fields...
        //
        for (CommodityIndex = 0; CommodityIndex < static_cast<size_t>(NumberOfCommodities); ++CommodityIndex) {
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
        for (CommodityIndex = 0; CommodityIndex < static_cast<size_t>(NumberOfCommodities); ++CommodityIndex) {
            LOG(INFO) << "Commodity " << CommodityIndex << ": <exp.> = " << AlphaData.GetData(CommodityIndex)->GetAverage() << "; s.d. = " << AlphaData.GetData(CommodityIndex)->GetStdDev() <<
            "; <endow.> = " << EndowmentData.GetData(CommodityIndex)->GetAverage() << "; s.d. = " << EndowmentData.GetData(CommodityIndex)->GetStdDev() << "; <MRS> = " <<
            LnMRSsData.GetData(CommodityIndex)->GetExpAverage() << "; s.d. = " << LnMRSsData.GetData(CommodityIndex)->GetStdDev();
        }
    }
    LOG(INFO) << "Average initial wealth (@ own prices) = " << InitialOwnWealthData.GetAverage() << "; standard deviation = " << InitialOwnWealthData.GetStdDev();
    LOG(INFO) << "Initial sum of utilities = " << LastSumOfUtilities;
}   //  AgentPopulation::Init()

void AgentPopulation::Reset() {
    for (auto& agent : Agents) {
        agent->Reset();
    }
}   //  AgentPopulation::Reset

long long AgentPopulation::Equilibrate(int NumberOfEquilibrationsSoFar) {
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
    for (auto& vol : Volume) {
        vol = 0.0;
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
                Commodity1 = 0;
                Commodity2 = 1;
            } else {
                Commodity1 = randomCommodity(rng);
                do {
                    Commodity2 = randomCommodity(rng);
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
                default:
                LOG(ERROR) << "Invalid termination criterion";
                std::terminate();
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
                    default:
                    LOG(ERROR) << "Invalid termination criterion";
                    std::terminate();
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

    for (auto& agent : Agents) {
        //  First compute agent wealth one commodity at a time...
        //
        AgentsInitialWealth = 0.0;
        AgentsFinalWealth = 0.0;
        for (size_t j = 0; j < static_cast<size_t>(NumberOfCommodities); ++j) {
            price = LnMRSsData.GetData(j)->GetExpAverage();
            AgentsInitialWealth += agent->GetEndowment(j) * price;
            AgentsFinalWealth += agent->GetAllocation(j) * price;
        }
        InitialMarketWealthData.AddDatuum(AgentsInitialWealth);
        FinalMarketWealthData.AddDatuum(AgentsFinalWealth);
        DeltaMarketWealthData.AddDatuum(AgentsFinalWealth - AgentsInitialWealth);

        //  The following line is a minor optimization...variable name should include 'own' to be mnemonic...
        //
        AgentsFinalWealth = agent->Wealth(agent->GetCurrentMRSs());
        FinalOwnWealthData.AddDatuum(AgentsFinalWealth);
        DeltaOwnWealthData.AddDatuum(AgentsFinalWealth - agent->GetInitialWealth());
        DeltaUtility.AddDatuum(agent->Utility() - agent->GetInitialUtility());
    }   //  for i...

    LOG(INFO) << "Average initial wealth (@ market prices) = " << InitialMarketWealthData.GetAverage() << "; standard deviation = " << InitialMarketWealthData.GetStdDev();
    LOG(INFO) << "Average final wealth (@ market prices)   = " << FinalMarketWealthData.GetAverage() << "; standard deviation = " << FinalMarketWealthData.GetStdDev();
    LOG(INFO) << "Average change in market wealth          = " << DeltaMarketWealthData.GetAverage() << "; standard deviation = " << DeltaMarketWealthData.GetStdDev();
    LOG(INFO) << "Average final wealth (@ own prices)      = " << FinalOwnWealthData.GetAverage() << "; standard deviation = " << FinalOwnWealthData.GetStdDev();
    LOG(INFO) << "Average change in own wealth             = " << DeltaOwnWealthData.GetAverage() << "; standard deviation = " << DeltaOwnWealthData.GetStdDev();
    LOG(INFO) << "Minimum increase in utility              = " << DeltaUtility.GetMin() << "; maximum increase = " << DeltaUtility.GetMax();
    LOG(INFO) << "Final sum of utilities                   = " << ComputeSumOfUtilities();
    if (PrintFinalCommodityList) {
        for (size_t i = 0; i < static_cast<size_t>(NumberOfCommodities); ++i) {
            LOG(INFO) << "Commodity " << i << ": volume = " << VolumeStats[i] << "; avg. MRS = " << LnMRSsData.GetData(i)->GetExpAverage() <<
            "; s.d. = " << LnMRSsData.GetData(i)->GetStdDev();
        }
    }
}   //  AgentPopulation::ConvergenceStatistics()

void AgentPopulation::CompareTwoAgents(AgentPtr Agent1, AgentPtr Agent2) {
    for (size_t j = 0; j < static_cast<size_t>(NumberOfCommodities); ++j) {
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
    if (debug) { LOG(DEBUG) << "Shocking agent preferences..."; }
    double oldPref, newPref, pref;

    size_t CommodityToShock = randomCommodity(rng);
    bool sign = randomBinary(rng);
    double shock = randomShock(rng);
    if (debug) {
        if (sign) {
            LOG(DEBUG) << "Shocking commodity " << CommodityToShock << " * " << shock;
        } else {
            LOG(DEBUG) << "Shocking commodity " << CommodityToShock << " * 1/" << shock;
        }
    }
    for (auto& ActiveAgent : Agents) {
        oldPref = ActiveAgent->GetAlpha(CommodityToShock);
        if (sign) {
            newPref = oldPref * shock;
        } else {
            newPref = oldPref/shock;
        }
        ActiveAgent->SetAlpha(CommodityToShock, newPref);

        for (size_t CommodityIndex = 0; CommodityIndex < static_cast<size_t>(NumberOfCommodities); ++CommodityIndex) {
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
        case -2:
        LOG(INFO) << "Termination criterion: After " << TerminationTime << " time steps";
        break;
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
        default:
        LOG(ERROR) << "Invalid termination criterion";
        std::terminate();
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
        std::random_device rd;
        unsigned int seed = rd();
        rng.seed(seed);
        LOG(INFO) << "Using random seed " << seed;
    } else {
        LOG(INFO) << "Using fixed seed " << NonRandomSeed;
        rng.seed(NonRandomSeed);
    }
    // initialize the distributions now that we know the relevant ranges
    randomAgent = std::uniform_int_distribution<unsigned long>(0, static_cast<size_t>(NumberOfAgents-1));
    randomCommodity = std::uniform_int_distribution<unsigned long>(0, static_cast<size_t>(NumberOfCommodities-1));
    randomBinary = std::uniform_int_distribution<int>(0, 1);
    randomShock = std::uniform_real_distribution<double>(MinShock, MaxShock);
    randomAlpha = std::uniform_real_distribution<double>(alphaMin, alphaMax);
    randomWealth = std::uniform_real_distribution<double>(wealthMin, wealthMax);
    randomDouble = std::uniform_real_distribution<double>(0, 1);
} // SeedRNG()

void ReadConfigFile(std::string file) {
    // This function sets all relevant model parameters by reading from a config file (libconfig). 
    // If there's some issue with the formatting or reading of the config file, it catches the exception
    // and terminates the program. The config file *must* be proper for the model to run. See parameters.cfg
    // in the repository for an example.

    libconfig::Config config;
    try { 
        LOG(INFO) << "Loading configuration from " << file << "...";
        config.readFile(file.c_str());
        LOG(INFO) << "Loaded configuration from " << file;
        if (config.lookupValue("debug.enabled", debug) && debug) { LOG(DEBUG) << "debug: " << debug; }
        if (config.lookupValue("number.commodities", NumberOfCommodities) && debug) { LOG(DEBUG) << "NumberOfCommodities: " << NumberOfCommodities; }
        if (config.lookupValue("number.agents", NumberOfAgents) && debug) { LOG(DEBUG) << "NumberOfAgents: " << NumberOfAgents; }
        if (config.lookupValue("rand.use_seed", UseRandomSeed) && debug) { LOG(DEBUG) << "UseRandomSeed: " << UseRandomSeed; }
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
        if (config.lookupValue("activation.method", activationMethod) && debug) { LOG(DEBUG) << "Activation Method: " << activationMethod; } // TODO: make this an enum
        exp_trade_eps = exp(trade_eps);

    } catch (...) {
        LOG(ERROR) << "Error reading config file";
        std::terminate();
    }
} //ReadConfigFile

int main(int argc, char** argv) {
    // Preliminaries: Parse flags, etc.
    std::string usage = "An agent-based model of bilateral exchange. Usage:\n";
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
    long long sum;
    if (DefaultSerialExecution) {
        sum = PopulationPtr->Equilibrate(EquilibrationNumber);
    } else {
        //sum = PopulationPtr->ParallelEquilibrate(EquilibrationNumber);
    }
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
}  // main()
