#include <Rcpp.h>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>
#include <regex>
#include "solver.h"
#include "parameters.h"


// helper function that parses ofstream into a dataframe
void stream_parser(std::stringstream &data_stream, 
                   Rcpp::DataFrame &df)
{
    // first put ostringstream back to beginning
    data_stream.seekp(0);

    // string containing the headers
    std::string headers{};

    data_stream.seekp(0);
    
    // then obtain first line, containing the headers
    std::getline(data_stream, headers);
    
    std::vector <std::string> all_lines{};
    
    std::string tmp;
    
//    // then the next lines which should just be data until we have a line which immediately ends
    while(!data_stream.eof())
    {
        std::getline(data_stream,tmp);
   
        // stop when the string is empty     
        // this is the end of the data
        if (tmp == "")
        {
            break;
        }
        
        all_lines.push_back(tmp);
    }
    
    // parse the headers
    std::istringstream header_to_parse(headers);
    
    std::vector <std::string> header_vec;
    
    // now split the headers into columns
    while(!header_to_parse.eof())
    {
        std::getline(header_to_parse, tmp, ';');
        
        header_vec.push_back(tmp);
    }
    
    // allocate a 2d vector for the data 
    std::vector < std::vector <double> > the_data(
        all_lines.size(), std::vector<double>(header_vec.size(), 0.0)
    );
    
    // now go through all lines and accumulate data
    // in the 2d vector
    for (unsigned row_idx{0}; row_idx < all_lines.size(); ++row_idx)
    {
        std::istringstream line_to_parse(all_lines[row_idx]);
        
        unsigned col_idx{0};
        
        while (!line_to_parse.eof())
        {
            std::getline(line_to_parse, tmp, ';');
            
            if (tmp=="")
            {
                continue;
            }
            
            the_data[row_idx][col_idx] = std::stod(tmp);
            ++col_idx;
        }
    }
    
    // move data from 2d vector to data frame
    // probably can all be done way more efficiently but it works
    for (unsigned col_idx{0}; col_idx < header_vec.size(); ++col_idx)
    {
        Rcpp::NumericVector vx(all_lines.size());
        
        for (unsigned row_idx{0}; row_idx < all_lines.size(); ++row_idx)
        {
            vx[row_idx] = the_data[row_idx][col_idx];
        }
        
        df[header_vec[col_idx]] = vx;
    }
} // end stream_parser()    

//' Run a single simulation of (S1) in Domingues, Kuijper et al.
//' Notation may differ from what is in the main manuscript 
//' P is the non-immune host that accepts all foreign genetic elements
//' C is the CRISPR-immune host that accepts only G2 foreign genetic elements
//' most importantly, G1 is a worse foreign genetic element than G2 
//'
//' @param kappa strength of density-dependence
//' @param gammaPG1 FGE loss rate in P (promiscuous) hosts of G1 FGEs
//' @param gammaPG2 FGE loss rate in P (promiscuous) hosts of G2 FGEs
//' @param gammaCG1 FGE loss rate in C (choosy) hosts of G1 FGEs
//' @param gammaCG2 FGE loss rate in C (choosy) hosts of G2 FGEs
//' @param dSP mortality rate of susceptible P (promiscuous) hosts 
//' @param dSC mortality rate of susceptible C (choosy) hosts 
//' @param dPG1 mortality rate in P (promiscuous) hosts of G1 FGEs
//' @param dPG2 mortality rate in P (promiscuous) hosts of G2 FGEs
//' @param dCG1 mortality rate in C (choosy) hosts of G1 FGEs
//' @param dCG2 mortality rate in C (choosy) hosts of G2 FGEs
//' @param psiG1 force of infection by G1 FGEs
//' @param psiG2 force of infection by G2 FGEs
//' @param FG1 fecundity benefit of the G1 FGE
//' @param FG2 fecundity benefit of the G2 FGE
//' @param FG1G2 fecundity benefit of the G1 + G2 FGEs when superinfected
//' @param pi probability of resisting infection by the G1 FGE
//' @param c cost of resistance
//' @param sigma probability of superinfection
//' @param init_S initial density of S
//' @param init_fraction_SC initial fraction of S that has the C allele, remainder has the P allele
//' @param init_PG1 initial density of P hosts infected with G1 
//' @param init_PG2 initial density of P hosts infected with G2
//' @param init_CG1 initial density of C hosts infected with G1
//' @param init_CG2 initial density of C hosts infected with G2
//' @param demog_feedback boolean allowing for demographical feedbacks or not
//' @param debug boolean whether to give debug output or not
//' @param eul Euler's constant
//' @param max_time The maximum number of time steps the numerical simulation shoujld run for
//' @returns A \code{data.frame} that contains the various frequencies
//' @export
//' @rawNamespace importFrom(Rcpp, sourceCpp);useDynLib("CRISPRhitchhike");
// [[Rcpp::export]]
Rcpp::DataFrame CRISPRhh(
        double kappa=0.001
        ,double gammaPG1=1
        ,double gammaPG2=1
        ,double gammaCG1=1
        ,double gammaCG2=1
        ,double dSP=5
        ,double dSC=5
        ,double dPG1=5
        ,double dPG2=1
        ,double dCG1=5
        ,double dCG2=1
        ,double psiG1=10
        ,double psiG2=10
        ,double FG1=6
        ,double FG2=10
        ,double FG1G2=0
        ,double pi=0.5
        ,double c=0.02
        ,double sigma=0.0
        ,double init_S=50
        ,double init_fraction_SC=0.5
        ,double init_PG1=1
        ,double init_PG2=1
        ,double init_CG1=1
        ,double init_CG2=1
        ,bool demog_feedback=true
        ,bool debug=false
        ,double eul=0.01
        ,unsigned long max_time=5e08
        )
{
    Parameters pars;
    // set parameters
    pars.kappa = kappa;
    
    pars.gamma[P][G1] = gammaPG1;
    pars.gamma[P][G2] = gammaPG2;
    pars.gamma[C][G1] = gammaCG1;
    pars.gamma[C][G2] = gammaCG2;
    
    pars.psi[G1] = psiG1;
    pars.psi[G2] = psiG2;
    
    pars.F[G1] = FG1;
    pars.F[G2] = FG2;
    
    pars.FGB = FG1G2;
    
    pars.pi = pi;
    pars.c = c;
    
    pars.dS[P] = dSP;
    pars.dS[C] = dSC;
    pars.dI[P][G1] = dPG1;
    pars.dI[P][G2] = dPG2;
    pars.dI[C][G1] = dCG1;
    pars.dI[C][G2] = dCG2;
    
    pars.init_popsize[P] = init_S * (1.0 - init_fraction_SC);
    pars.init_popsize[C] = init_S * init_fraction_SC;
    
    pars.init_popsize_infected[P][G1] = init_PG1;
    pars.init_popsize_infected[P][G2] = init_PG2;
    
    // with the initial population size of CG1
    // we need to make sure that fully resistant individuals
    // cannot be infected by G1, as they are fully resistant...
    pars.init_popsize_infected[C][G1] = init_CG1;
    pars.init_popsize_infected[C][G2] = init_CG2;
    
    pars.demog_feedback = demog_feedback;
    pars.debug = debug;
    
    pars.eul = eul;

    pars.max_ecol_time = max_time;
    
    // stringstream to store model output from the solver
    std::stringstream model_output;
    
    Solver sol{pars, model_output};
    
    if (pars.debug)
    {
        Rcpp::Rcout << model_output.str() << std::endl;
    }

    Rcpp::DataFrame time_data;
    stream_parser(model_output, time_data);
    
    return(time_data);
} // CRISPRhh
