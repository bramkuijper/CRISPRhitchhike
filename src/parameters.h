#ifndef _PARAMETERS_HPP_
#define _PARAMETERS_HPP_

#include <string>
#include <limits>


class Parameters
{
    public:
        // FGE loss rate
        // first index is host_type, (P = 0, C = 1)
        // next index is phage_type, (G1 = 0, G2 = 1)
        double gamma[2][2]{{1.0,1.0},{1.0,1.0}};

        // force of infection
        // by G1 and G2 respectively
        double psi[2]{1.0,1.0};

        // fecundity of susceptibles
        double FS{1.0};

        // fecundity of individuals infected with 
        // G1 or G2
        double F[2]{1.0,1.0};
       
        // fecundity upon superinfection  
        double FGB{0.0};
        
        // probability of resisting infection by G1
        double pi{1.0};
        
        // cost of resistance
        double c{1.0};
        // strenght of density dependence
        double kappa{0.001};
        
        // probability of superinfection
        double sigma{0.0};
        
        // death rate of P and C susceptibles
        double dS[2]{1.0,1.0};
        
        // death rate of infected individuals
        // first index host P or C
        // second index FGE G1 or G2
        double dI[2][2]{{1.0,1.0},{1.0,1.0}};
        
        // death rate fo superinfected individuals
        double dBG[2]{1.0,1.0};

        // whether to write output to a file or return 
        // everything as a bunch of vectors
        std::string base_name{"sim_fixed_resistance"};

        // initial population size
        // S_p, S_c, I_pg1, I_pg2, I_cg1, I_cg2
        double init_popsize[2] = {100,100};
        double init_popsize_infected[2][2] = {{1,1},{1,1}};
        double init_popsize_superinfected[2] = {0,0};

        // double minimum boundery
        // this is to output files to string without having
        // later c++ code (i.e., stod) choke on reading back in
        // the output
        double dbl_min_val{std::numeric_limits<double>::lowest()*2};
        double dbl_max_val{std::numeric_limits<double>::max()*2};
        
        // Euler's constant
        double eul{0.001};

        // maximum number of time steps to solve over
        long unsigned max_ecol_time{50'000'000};

        // print every n lines
        long unsigned print_interval{1000};

        // threshold below which we assume an equation has vanished
        double vanish_threshold = 1e-07;

        // whether there should be demographical feedback in the model
        bool demog_feedback = false;
        bool debug = false;
};

#endif
