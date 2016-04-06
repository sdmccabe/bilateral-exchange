#Bilateral exchange model

This is an agent-based model (ABM) of economic exchange, originally presented in Axtell (2005) and extended for my master's thesis in McCabe (2016). 

Dependencies:

- GCC 5+
- Intel Threading Building Blocks
- libconfig++

Usage:
The model takes a configuration file as an argument. If no argument is passed, it attempts to read the file "parameters.cfg" in the local directory. See the source for parameters.cfg for more information on parameters, etc.

A Dockerfile is included for running the model (and its associated batch scripts) in a Docker session with Ubuntu 15.10 in case GCC, TBB, or libconfig++ are not available on your system. 


##References
Axtell, R. L. (2005). The complexity of exchange. The Economic Journal, 115(504), F193–F210. doi:10.1111/j.1468-0297.2005.01001.x

McCabe, S. D. (2016). Communicating sequential agents: An analysis of concurrent agent scheduling (Master’s thesis). George Mason University, Fairfax, VA.
