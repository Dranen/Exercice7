#ifndef EXERCICE7_IO_H
#define EXERCICE7_IO_H

#include <vector>
#include <iostream>

#include "Exercice7_calcul.h"

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

void contexte_general(std::string& nom ,int& Ninter, double& xL, double& xR, int& eqref, int& ucase)
{
    std::cerr << "Nom de la simulation " << std::flush;
    std::cin >> nom;
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
        std::cerr << "value of u? " << std::flush;
        std::cin >> u;
        u2 = u_squared(u);
        u2_max = u2(0);
        break;
      case 1:
        std::cerr << "u_L? " << std::flush;
        std::cin >> uL;
        std::cerr << "u_R? " << std::flush;
        std::cin >> uR;
        u2 = u_squared(uL,uR,xL,xR);
        for(int ip = 0; ip < Npos; ++ip)
      if(u2(xL+ip*dx) > u2_max)
        u2_max = u2(xL+ip*dx);
        break;
      case 2:
        std::cerr << "hOcean? " << std::endl;
        std::cin >> hocean;
        std::cerr << "xOcean? " << std::endl;
        std::cin >> xocean;
        std::cerr << "hPlage? " << std::endl;
        std::cin >> hplage;
        u2 = u_squared(xL, xR, hocean, xocean, hplage);
        for(int ip = 0; ip < Npos; ++ip)
      if(u2(xL+ip*dx) > u2_max)
        u2_max = u2(xL+ip*dx);
        break;
     }

     std::cout << "u2_max="<<u2_max <<std::endl;
}

void contexte_temporelle(double& dt, double& tfinal, double const& dx, double const& u2_max)
{
    double CFL;
    std::cerr << "max CFL coefficent? " << std::flush;
    std::cin >> CFL;

    dt = CFL * dx / sqrt(u2_max);

    std::cout << "dt=" << dt << std::endl;
    std::cerr << "tfinal? " << std::flush;
    std::cin >> tfinal;
}

void contexte_bord(std::string bord, boundary_condition& bc, double& boundaryvalue, double& A, double& omega)
{
    int n;
    std::cerr << bord << " boundary condition (fixed = 0, free = 1, excited = 2, outgoing = 3)? " << std::flush;
    std::cin >> n;

    switch(n)
      {
      case 0:
        bc = fixed;
        break;

      case 1:
        bc = libre;
        break;

      case 2:
        bc = excited;
        break;

      case 3:
        bc = outgoing;
        break;

      default:
        std::cerr << "no valid " << bord << " boundary condition!";
        throw "no valid boundary condition!";
      }

  if(bc == fixed)
  {
    std::cerr << "for fixed left boundary condition: value of f at the left boundary? " << std::flush;
    std::cin >> boundaryvalue;
  }

  if(bc == excited)
    {
      std::cerr << "Amplitude of the excitation? " << std::flush;
      std::cin >> A;
      std::cerr << "Angular frequency omega of the excitation? " << std::flush;
      std::cin >> omega;
    }
}

void contexte_scan_frequence(int& nscan, double& omega, double& omega_start, double& omega_stop, double& delta_omega)
{
    omega_start = omega;
    omega_stop = omega;

    std::cerr << "number of frequencies to scan? " << std::flush;
    std::cin >> nscan;

    if (nscan > 1) // do a frequency scan only if you are interested in more than 1 frequency
      {
        std::cerr << "scan from omega-start=omega to omega-end=? " << std::flush;
        std::cin >> omega_stop;
        omega_start = omega;
        delta_omega = (omega_stop - omega_start)/nscan;
      }

    omega = omega_start;
}

#endif // EXERCICE7_IO_H
