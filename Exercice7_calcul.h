#ifndef EXERCICE7_CALCUL_DONNEE_H
#define EXERCICE7_CALCUL_DONNEE_H

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <vector>

typedef enum {fixed,free,excited,outgoing} boundary_condition;


// function object for the position dependent velocity
class u_squared {
public:

  // standard constructor
  u_squared()
  {}

  // constant u: set all private variables to zero but u2_ and uconst_
  u_squared(double u) :
  u2(u*u),ucase(0),g(9.81)
  {}
  u_squared(double u_L,double u_R,double xL_, double xR_) :
  u2_L(u_L*u_L),u2_R(u_R*u_R),xL(xL_),xR(xR_),hocean(0),ucase(1),g(9.81)
  {}

  // variable u: use all possible parameters
  u_squared(double xL_, double xR_, double hOcean_, double xOcean_, double hPlage_) :
    xL(xL_),xR(xR_),hocean(hOcean_),xocean(xOcean_),hplage(hPlage_),ucase(2) ,g(9.81)
  {}

  // return the value of u^2(x) at position x for all possible cases:
  double operator()(double x)
  {
    switch(ucase)
    {
    case 0:
      break;

    case 1:
      u2 = u2_L+(u2_R-u2_L)*(x-xL)/(xR-xL);
      break;

    case 2:

    if (x <= xocean)
          u2 = g * hocean;
    else
      {
            double tmp = sin(M_PI * (x - xocean) / (2.0 * (xR - xocean)));
            u2 = g * (hocean + (hplage - hocean) * tmp*tmp );
      }
     }

    double result = u2;

    assert(result > 0);

    return result;

  }

private:
  double u2_L,u2_R,u2,  xL, xR, hocean,xocean,hplage, g;
  int ucase;
};


// calculate the energy of the wave
double get_energy(const std::vector<double>& f, const double dx)
{
  int npos = f.size();

  double erg = 0.;
  for(int ip = 0; ip < (npos - 1); ++ip)
    erg +=0; //TODO: calculer l'energie.

  return erg * dx;
}


inline void calcul_condition_bord(bool gauche, boundary_condition const& bc, double const& boundaryvalue, double const& A, double const& omega, double const& t, std::vector<double> *fnow, std::vector<double> *fnext, int i, double beta)
{
    switch(bc)
      {
      case fixed:
        (*fnext)[i] = boundaryvalue;
        break;

      case free:
        if(gauche)
        {
            (*fnext)[i] = (*fnext)[i+1];
        }
        else
        {
            (*fnext)[i] = (*fnext)[i-1];
        }
        break;

      case excited:
        fnext->[i] = A * sin(omega * t);
        break;

      case outgoing:
        if(gauche)
        {
            (*fnext)[i] = (*fnow)[i]-beta*((*fnow)[i+1]-(*fnow)[i]);
        }
        else
        {
            (*fnext)[i] = (*fnow)[i]-beta*((*fnow)[i]-(*fnow)[i-1]);
        }
        break;

      default:
        std::cerr << "No valid boundary condition given" <<std::endl;
        break;
      }
}

#endif // EXERCICE7_CALCUL_DONNEE_H
