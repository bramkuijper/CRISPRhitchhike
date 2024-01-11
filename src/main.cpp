#include <Rcpp.h>
#include <cstdio>
#include <sstream>
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

//' Run a single simulation of (S1) in Domingues, Kuijper et al 
//' @param kappa strength of density-dependence
//' @returns A \code{data.frame} that contains the various frequencies
//' @export
//' @rawNamespace importFrom(Rcpp, sourceCpp);useDynLib("CRISPRhitchhike");
// [[Rcpp::export]]
Rcpp::DataFrame CRISPRhh(
        double kappa=0.001
        )
{
    Parameters pars;
    pars.kappa = kappa;
    
    std::stringstream model_output;
    
    Solver sol{pars, model_output};

    Rcpp::DataFrame time_data;
    stream_parser(model_output, time_data);
    
    return(time_data);
}
