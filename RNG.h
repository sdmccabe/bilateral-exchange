// Copyright 2015 <Stefan McCabe>
/*
Bilateral Exchange Model
Originally written by Rob Axtell
Extended by Stefan McCabe

12/11/2015
*/

#ifndef RNG_H_
#define RNG_H_
#endif  // RNG_H_

#include <memory>
#include <random>
#include <string>
#include <vector>
#include <fstream>

class RNG {
    unsigned int seed;
    std::mt19937 rng;

    unsigned long num_agents, num_com;
    double alpha_min, alpha_max, shock_min, shock_max, wealth_min, wealth_max;

 public:
    RNG (bool randSeed, unsigned int s, unsigned int numagents, unsigned int numcom, double shockmin, double shockmax, double minalpha, double maxalpha, double minwealth, double maxwealth);
    unsigned int GetSeed() { return seed; }
    void SetSeed(unsigned int s) {
        seed = s;
        rng.seed(s);
    }
    std::mt19937 GetGenerator() { return rng; }

    int ValueInRange(int min, int max);
    unsigned long ValueInRange(unsigned long min, unsigned long max);
    double ValueInRange(double min, double max);
    
    unsigned long RandomAgent() { return ValueInRange(0L, num_agents-1); }
    unsigned long RandomCommodity() { return ValueInRange(0L, num_com-1); }
    int RandomBinary() { return ValueInRange(0, 1); }
    double RandomShock() { return ValueInRange(shock_min, shock_max); }
    double RandomAlpha() { return ValueInRange(alpha_min, alpha_max); }
    double RandomWealth() { return ValueInRange(wealth_min, wealth_max); }
    double RandomDouble() { return ValueInRange(0.0, 1.0); }
};

typedef RNG *RNGptr;
