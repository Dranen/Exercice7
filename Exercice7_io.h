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
    std::cerr << std::endl << "Number of grid points? " << std::flush;
    std::cin >> Ninter;
    std::cerr << std::endl << "xL? " << std::flush;
    std::cin >> xL;
    std::cerr << std::endl << "xR? " << std::flush;
    std::cin >> xR;
    std::cerr << std::endl << "Quelle equation (1) ou (2) ou (3) ou (4)?"<< std::flush;
    std::cin >> eqref;
    std::cerr << std::endl << "quel question? 7.2 (taper 0) / 7.3 (taper 1) / 7.4 (taper 2):" << std::flush;
    std::cin >> ucase;
}

void contexte_vitesse(int const& ucase, double const& xL, double const& xR, u_squared& u2)
{
    double uL, uR, u, hocean = 0, xocean=0, hplage=0;
    switch(ucase){
      case 0:
        std::cerr << std::endl << "value of u? " << std::flush;
        std::cin >> u;
        u2 = u_squared(u);
        break;
      case 1:
        std::cerr << std::endl << "u_L? " << std::flush;
        std::cin >> uL;
        std::cerr << std::endl << "u_R? " << std::flush;
        std::cin >> uR;
        u2 = u_squared(uL,uR,xL,xR);
        break;
      case 2:
        std::cerr << std::endl << "hOcean? " << std::endl;
        std::cin >> hocean;
        std::cerr << std::endl << "xOcean? " << std::endl;
        std::cin >> xocean;
        std::cerr << std::endl << "hPlage? " << std::endl;
        std::cin >> hplage;
        u2 = u_squared(xL, xR, hocean, xocean, hplage);
        break;
     }
}

void contexte_temporelle(double& tfinal, double& CFL)
{
    std::cerr << std::endl << "max CFL coefficent? " << std::flush;
    std::cin >> CFL;

    std::cerr << std::endl << "tfinal? " << std::flush;
    std::cin >> tfinal;
}

void contexte_bord(std::string bord, boundary_condition& bc, double& boundaryvalue, double& A, double& omega)
{
    int n;
    std::cerr << std::endl << bord << " boundary condition (fixed = 0, free = 1, excited = 2, outgoing = 3)? " << std::flush;
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
        std::cerr << std::endl << "no valid " << bord << " boundary condition!";
        throw "no valid boundary condition!";
      }

  if(bc == fixed)
  {
    std::cerr << std::endl << "for fixed left boundary condition: value of f at the left boundary? " << std::flush;
    std::cin >> boundaryvalue;
  }

  if(bc == excited)
    {
      std::cerr << std::endl << "Amplitude of the excitation? " << std::flush;
      std::cin >> A;
      std::cerr << std::endl << "Angular frequency omega of the excitation? " << std::flush;
      std::cin >> omega;
    }
}

void contexte_mode(int& nscan, double& omega_stop, int& Ninter_stop, double& CFL_stop, mode& choix, int& ech_t)
{
    int n;
    std::cerr << std::endl << " mode (Unique = 0, scan frequence = 1, convergence CFL = 2, convergence dx = 3)? " << std::flush;
    std::cin >> n;

    switch(n)
    {
      case 0:
        choix = Unique;
        break;

      case 1:
        choix = frequence;
        break;

      case 2:
        choix = convergence_CFL;
        break;

      case 3:
        choix = convergence_Ninter;
        break;

      default:
        std::cerr << std::endl << "no valid mode";
        throw "no valid mode";
    }

    if (choix == Unique)
    {
        std::cerr << std::endl << "Echantillonage x ?" << std::flush;
        std::cin >> ech_t;
    }
    else if(choix == frequence)
    {
        std::cerr << std::endl << "number of frequencies to scan? " << std::flush;
        std::cin >> nscan;

        if (nscan > 1) // do a frequency scan only if you are interested in more than 1 frequency
        {
            std::cerr << std::endl << "scan from omega-start=omega to omega-end=? " << std::flush;
            std::cin >> omega_stop;
        }
        else
        {
            nscan = 1;
            choix = Unique;
        }
    }
    else if(choix == convergence_CFL)
    {
        std::cerr << std::endl << "number of CFL to scan? " << std::flush;
        std::cin >> nscan;

        if (nscan > 1)
        {
            std::cerr << std::endl << "scan from CFL-start=CFL to CFL-end=? " << std::flush;
            std::cin >> CFL_stop;
        }
        else
        {
            nscan = 1;
            choix = Unique;
        }
    }
    else if(choix == convergence_Ninter)
    {
        std::cerr << std::endl << "number of Ninter to scan? " << std::flush;
        std::cin >> nscan;

        if (nscan > 1)
        {
            std::cerr << std::endl << "scan from Ninter-start=Ninter to Ninter-end=? " << std::flush;
            std::cin >> Ninter_stop;
        }
        else
        {
            nscan = 1;
            choix = Unique;
        }
    }
    else
    {
        nscan = 1;
        choix = Unique;
    }
}

#endif // EXERCICE7_IO_H
