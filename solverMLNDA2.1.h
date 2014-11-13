#ifndef SOLVERMLNDA_H // Include gard
#define SOLVERMLNDA_H

#include <vector>

void GNDAsolutionFS(int);

void GNDAsolutionKE(int, double , std::vector< std::vector<double> > );

void NDAsolution(int, int , int , std::vector< std::vector<double> >, std::vector< std::vector<double> >, std::vector< std::vector<double> >);

#endif // SOLVERMLNDA_H


