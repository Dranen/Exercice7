#ifndef EXERCICE7_IO_H
#define EXERCICE7_IO_H

#include <vector>
#include <iostream>

#include "Exercice7_calcul_donnee.h"

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

void contexte_vitesse(int const& ucase, double const& xL, double const& xR, double const& dx, int const& Npos, double& u2_max, u_squared& u2)
{
    double uL, uR, u, hocean = 0, xocean=0, hplage=0;
    switch(ucase){
      case 0:
        cerr << "value of u? " << flush;
        cin >> u;
        u2 = u_squared(u);
        u2_max = u2(0);
        break;
      case 1:
        cerr << "u_L? " << flush;
        cin >> uL;
        cerr << "u_R? " << flush;
        cin >> uR;
        u2 = u_squared(uL,uR,xL,xR);
        for(int ip = 0; ip < Npos; ++ip)
      if(u2(xL+ip*dx) > u2_max)
        u2_max = u2(xL+ip*dx);
        break;
      case 2:
        cerr << "hOcean? " << endl;
        cin >> hocean;
        cerr << "xOcean? " << endl;
        cin >> xocean;
        cerr << "hPlage? " << endl;
        cin >> hplage;
        u2 = u_squared(xL, xR, hocean, xocean, hplage);
        for(int ip = 0; ip < Npos; ++ip)
      if(u2(xL+ip*dx) > u2_max)
        u2_max = u2(xL+ip*dx);
        break;
     }

     cout << "u2_max="<<u2_max <<endl;
}

void contexte_temporelle(double& dt, double& tfinal)
{
    double CFL;
    cerr << "max CFL coefficent? " << flush;
    cin >> CFL;

    dt = CFL * dx / sqrt(u2_max);

    cout << "dt=" << dt << endl;
    cerr << "tfinal? " << flush;
    cin >> tfinal;
}

#endif // EXERCICE7_IO_H
