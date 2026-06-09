#ifndef _SOLVER_HPP_
#define _SOLVER_HPP_

#include <iostream>
#include <fstream>
#include <sstream>
#include "parameters.h"

enum HostType
{
    P = 0, // promiscuous host
    C = 1 // choosy host
};

enum PhageType
{
    G1 = 0,  // single resistant phage or plasmid
    G2 = 1 // double resistant phage or plasmid
};

class Solver
{
    public:

        // parameter objedt
        Parameters params;

        // reference to data file
        std::ostream &data_file;

        double popsize[2] = {0.0,0.0};

        // host type x phage type
        double popsize_infected[2][2] = {
            {0.0,0.0}
            ,{0.0,0.0}
        };

        double popsize_superinfected[2] = {0,0};

        double N = 0.0;

        unsigned long time_step = 0;

        Solver(Parameters const &par
                ,std::ostream &data_output);

        // update the population size
        void update_N();

        double dSdt(HostType host_idx) const;

        double dIdt(HostType host_idx, 
                PhageType phage_idx) const;

        double dIGBdt(HostType host_idx) const;

        // fecundity of susceptible individuals
        double b(HostType host_idx) const;
        // death rate of susceptible individuals
        double d(HostType host_idx) const;

        // fecundity rate of infected individuals
        double b(HostType host_idx, PhageType phage_idx) const;
        // death rate of infected individuals
        double d(HostType host_idx, PhageType phage_idx) const;

        // fecundity of superinfected individuals
        double b(HostType host_idx, PhageType phage1_idx, 
                PhageType phage2_idx) const;
        
        // death rate of superinfected individuals
        double dBG(HostType host_idx) const;

        void write_data();
        void write_data_headers();
        void write_parameters();
        
        double clamp(double const val);
};

#endif
