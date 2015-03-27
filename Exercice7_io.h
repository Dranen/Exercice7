#ifndef EXERCICE7_IO_H
#define EXERCICE7_IO_H

#include <vector>
#include <iostream>

//
// this is specialication function for the simplified output
// of std:vector's without the need of writing all the
// elements oneself
//
template <class T>
std::ostream& operator << (std::ostream& o, const std::vector<T>& v)
{
  int len = v.size();
  for(int i = 0; i < (len - 1); ++i)
    o << v[i] << " ";

  if(len > 0)
    o << v[len-1];

  return o;
}

void contexte_general(int& Ninter, double& xL, double& xR, int& eqref, int& ucase)
{
    std::cerr << "Number of grid points? " << std::flush;
    std::cin >> Ninter;
    std::cerr << "xL? " << std::flush;
    std::cin >> xL;
    std::cerr << "xR? " << std::flush;
    std::cin >> xR;
    std::cerr << "Quelle equation (1) ou (2)?"<< std::flush;
    std::cin >> eqref;
    std::cerr << "quel question? 7.2 (taper 0) / 7.3 (taper 1) / 7.4 (taper 2):" << std::flush;
    std::cin >> ucase;
}


#endif // EXERCICE7_IO_H
