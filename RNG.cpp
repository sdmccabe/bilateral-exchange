// Copyright 2015 <Stefan McCabe>
/*
Bilateral Exchange Model
Originally written by Rob Axtell
Extended by Stefan McCabe

12/11/2015
*/

// The key mechanism here - the use of thread_local mt19337 pointers - 
// was suggested on Stack Exchange.
// https://stackoverflow.com/questions/21237905/how-do-i-generate-thread-safe-uniform-random-numbers

#include "./RNG.h"
#include <thread>
#define ELPP_THREAD_SAFE
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wundef"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wctor-dtor-privacy"
#include "./easylogging++.h"
#pragma GCC diagnostic pop

int RNG::ValueInRange(int min, int max) {
    static thread_local std::mt19937* generator = nullptr;
    if (!generator) generator = new std::mt19937(std::clock() + std::hash<std::thread::id>()(std::this_thread::get_id()));
    std::uniform_int_distribution<int> distribution(min, max);
    return distribution(*generator);
}

unsigned long RNG::ValueInRange(unsigned long min, unsigned long max) {
    static thread_local std::mt19937* generator = nullptr;
    if (!generator) generator = new std::mt19937(std::clock() + std::hash<std::thread::id>()(std::this_thread::get_id()));
    std::uniform_int_distribution<unsigned long> distribution(min, max);
    return distribution(*generator);
}

double RNG::ValueInRange(double min, double max) {
    static thread_local std::mt19937* generator = nullptr;
    if (!generator) generator = new std::mt19937(std::clock() + std::hash<std::thread::id>()(std::this_thread::get_id()));
    std::uniform_real_distribution<double> distribution(min, max);
    return distribution(*generator);
}

RNG::RNG (bool randSeed, unsigned int s, unsigned int numagents, unsigned int numcom, double shockmin, double shockmax, double minalpha, double maxalpha, double minwealth, double maxwealth) {
	// Seed the random number generator.
    if (!randSeed) {
        std::random_device rd;
        seed = rd();
        rng.seed(seed);
        LOG(INFO) << "Using random seed " << seed;
    } else {
    	seed = s;
        LOG(INFO) << "Using fixed seed " << seed;
        rng.seed(s);
    }

	alpha_min = minalpha;
	alpha_max = maxalpha;
	shock_min = shockmin;
	shock_max = shockmax;
	num_agents = static_cast<unsigned long>(numagents);
	num_com = static_cast<unsigned long>(numcom);
	wealth_min = minwealth;
	wealth_max = maxwealth;
}
