#ifndef EXERCICE7_CALCUL_DONNEE_H
#define EXERCICE7_CALCUL_DONNEE_H

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <chrono>
#include <fstream>

#include "Exercice7_io.h"

template <class T>
std::ostream& operator << (std::ostream& o, const std::vector<T>& v);

typedef enum {fixed,libre,excited,outgoing} boundary_condition;
typedef enum {Unique, frequence, convergence_CFL, convergence_Ninter} mode;


// function object for the position dependent velocity
class u_squared {
public:

  // standard constructor
  u_squared() :
  u2(1), ucase(0)
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
    erg += f[ip]*f[ip] + f[ip+1]*f[ip+1];

  return erg * dx;
}


inline void calcul_condition_bord(bool gauche, boundary_condition const& bc, double const& boundaryvalue, double const& A, double const& omega, double const& t, std::vector<double> *fnow, std::vector<double> *fnext, int i, double beta)
{
    switch(bc)
      {
      case fixed:
        (*fnext)[i] = boundaryvalue;
        break;

      case libre:
        if(gauche)
        {
            (*fnext)[i] = (*fnow)[i+1];
        }
        else
        {
            (*fnext)[i] = (*fnow)[i-1];
        }
        break;

      case excited:
        (*fnext)[i] = A * sin(omega * t);
        break;

      case outgoing:
        if(gauche)
        {
            (*fnext)[i] = (*fnow)[i]+beta*((*fnow)[i+1]-(*fnow)[i]);
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

void simulation(std::vector<double> *u_1, std::vector<double> *beta, std::vector<double> *coeff, double const& dx, double const& dt, int const& eqref, int const& ucase, int const& Npos, double const& tfinal, double const& omega, double const& A, boundary_condition const& left_bc, double const& leftboundaryvalue, boundary_condition const& right_bc, double const& rightboundaryvalue, mode const& sortie, std::ofstream& w_ofs, std::ofstream& energy_ofs, std::ofstream& maxenergy_ofs, int ech_t)
{
    std::vector<double> *fpast = new std::vector<double>(Npos);
    std::vector<double> *fnow = new std::vector<double>(Npos);
    std::vector<double> *fnext = new std::vector<double>(Npos);
    std::vector<double> *temp = NULL;
    int count_t = 0;


     // initialize the first two time slices of f
     double xo=3e5;
     double x;
     for(int ip = 0; ip < Npos; ++ip)
     {
        if(ucase==2){
          x = static_cast<double>(ip)*dx-xo;
          (*fnow)[ip] = exp(-x*x/2.5e9);
          (*fpast)[ip]=(*fnow)[ip];
        }else{
          (*fpast)[ip] = 0.0;
          (*fnow)[ip] = 0.0;
        }
     }


    // output the positions and the phase velocity squared for plotting purposes
    double t = 0;
    double energy = get_energy(*fnow,dx);
    double maxenergy = energy;


        if(sortie == Unique)
        {
          w_ofs << t << " " << *fnow << std::endl;
          energy_ofs << t << " " << energy << std::endl;
        }

    // time loop --------------------------------
    std::chrono::high_resolution_clock::time_point end;
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    do
    {
        count_t++;
          //cout << " time evolution: " << t << endl;
          if(eqref == 1)
          {
              #pragma omp for simd
              for(int ip = 1; ip < (Npos - 1); ++ip)
              {
                  (*fnext)[ip]=-(*fpast)[ip]+2*(1-(*coeff)[ip])*(*fnow)[ip]+(*coeff)[ip]*((*fnow)[ip-1]+(*fnow)[ip+1]);
              }
          }
          else if(eqref == 2)
          {
              #pragma omp for simd
              for(int ip = 1; ip < (Npos - 1); ++ip)
              {
                  (*fnext)[ip] = 0.5*(*beta)[ip]*((*u_1)[ip+1]-(*u_1)[ip-1])*((*fnow)[ip+1]-(*fnow)[ip-1])*dt/dx + (*beta)[ip]*(*beta)[ip]*((*fnow)[ip+1]-2*(*fnow)[ip]+(*fnow)[ip-1]) + 2*(*fnow)[ip] - (*fpast)[ip];
              }
          }
          else
          {
                #pragma omp for simd
                for(int ip = 1; ip < (Npos - 1); ++ip)
                {
                    (*fnext)[ip] = 0.25*((*coeff)[ip+1]-(*coeff)[ip-1])*((*fnow)[ip+1]-(*fnow)[ip-1])*dt/dx + (*coeff)[ip]*((*fnow)[ip+1]-2*(*fnow)[ip]+(*fnow)[ip-1]) + 2*(*fnow)[ip] - (*fpast)[ip];
                }
          }


          // apply boundary conditions
          //
          // left boundary (at xL)
          calcul_condition_bord(true, left_bc, leftboundaryvalue, A, omega, t, fnow, fnext, 0, (*beta)[0]);
          // right boundary (at xr)
          calcul_condition_bord(false, right_bc, rightboundaryvalue, A, omega, t, fnow, fnext, Npos-1, (*beta)[Npos-1]);

          t += dt;

          temp = fpast;
          fpast = fnow;
          fnow  = fnext;
          fnext = temp;

          double energy = get_energy((*fnow),dx);

          if(energy > maxenergy)
            maxenergy = energy;

          if(sortie == Unique)
          {
              if(count_t >= ech_t)
              {
                   w_ofs << t << " " << *fnow << std::endl;
                   energy_ofs << t << " " << energy << std::endl;
                   count_t = 0;
              }
          }

    } while(t < tfinal); // end of time loop -----------------------------
    end = std::chrono::high_resolution_clock::now();
    #pragma omp critical
    {
        std::cerr << std::endl << std::setprecision(std::numeric_limits<double>::digits10 + 1) << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::flush;
    }

    #pragma omp critical
    {
        if(sortie == convergence_CFL || sortie == convergence_Ninter)
        {
            maxenergy_ofs << dt << " " << dx << " " << maxenergy << std::endl;
        }
        else
        {
            maxenergy_ofs << omega << " " << maxenergy << std::endl;
        }
    }

    #pragma omp critical
    {
          if(sortie == convergence_CFL || sortie == convergence_Ninter)
          {
            w_ofs << dt << " " << dx << " " << *fnow << std::endl;
            energy_ofs << dt << " " << dx << " " << energy << std::endl;
          }
    }

    delete fpast;
    delete fnow;
    delete fnext;
}

void simulation_par(std::vector<double> *u_1, std::vector<double> *beta, std::vector<double> *coeff, double const& dx, double const& dt, int const& eqref, int const& ucase, int const& Npos, double const& tfinal, double const& omega, double const& A, boundary_condition const& left_bc, double const& leftboundaryvalue, boundary_condition const& right_bc, double const& rightboundaryvalue, mode const& sortie, std::ofstream& w_ofs, std::ofstream& energy_ofs, std::ofstream& maxenergy_ofs, int ech_t)
{
    std::vector<double> *fpast = new std::vector<double>(Npos);
    std::vector<double> *fnow = new std::vector<double>(Npos);
    std::vector<double> *fnext = new std::vector<double>(Npos);
    std::vector<double> *temp = NULL;
    int count_t = 0;


     // initialize the first two time slices of f
     double xo=3e5;
     double x;
     for(int ip = 0; ip < Npos; ++ip)
     {
        if(ucase==2){
          x = static_cast<double>(ip)*dx-xo;
          (*fnow)[ip] = exp(-x*x/2.5e9);
          (*fpast)[ip]=(*fnow)[ip];
        }else{
          (*fpast)[ip] = 0.0;
          (*fnow)[ip] = 0.0;
        }
     }


    // output the positions and the phase velocity squared for plotting purposes
    double t = 0;
    double energy = get_energy(*fnow,dx);
    double maxenergy = energy;


        if(sortie == Unique)
        {
          w_ofs << t << " " << *fnow << std::endl;
          energy_ofs << t << " " << energy << std::endl;
        }

    // time loop --------------------------------
    std::chrono::high_resolution_clock::time_point end;
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    do
    {
        count_t++;
          //cout << " time evolution: " << t << endl;
          if(eqref == 1)
          {
              #pragma omp parallel for simd
              for(int ip = 1; ip < (Npos - 1); ++ip)
              {
                  (*fnext)[ip]=-(*fpast)[ip]+2*(1-(*coeff)[ip])*(*fnow)[ip]+(*coeff)[ip]*((*fnow)[ip-1]+(*fnow)[ip+1]);
              }
          }
          else if(eqref == 2)
          {
              #pragma omp parallel for simd
              for(int ip = 1; ip < (Npos - 1); ++ip)
              {
                  (*fnext)[ip] = 0.5*(*beta)[ip]*((*u_1)[ip+1]-(*u_1)[ip-1])*((*fnow)[ip+1]-(*fnow)[ip-1])*dt/dx + (*beta)[ip]*(*beta)[ip]*((*fnow)[ip+1]-2*(*fnow)[ip]+(*fnow)[ip-1]) + 2*(*fnow)[ip] - (*fpast)[ip];
              }
          }
          else
          {
                #pragma omp parallel for simd
                for(int ip = 1; ip < (Npos - 1); ++ip)
                {
                    (*fnext)[ip] = 0.25*((*coeff)[ip+1]-(*coeff)[ip-1])*((*fnow)[ip+1]-(*fnow)[ip-1])*dt/dx + (*coeff)[ip]*((*fnow)[ip+1]-2*(*fnow)[ip]+(*fnow)[ip-1]) + 2*(*fnow)[ip] - (*fpast)[ip];
                }
          }


          // apply boundary conditions
          //
          // left boundary (at xL)
          calcul_condition_bord(true, left_bc, leftboundaryvalue, A, omega, t, fnow, fnext, 0, (*beta)[0]);
          // right boundary (at xr)
          calcul_condition_bord(false, right_bc, rightboundaryvalue, A, omega, t, fnow, fnext, Npos-1, (*beta)[Npos-1]);

          t += dt;

          temp = fpast;
          fpast = fnow;
          fnow  = fnext;
          fnext = temp;

          double energy = get_energy((*fnow),dx);

          if(energy > maxenergy)
            maxenergy = energy;

          if(count_t >= ech_t)
          {
               w_ofs << t << " " << *fnow << std::endl;
               energy_ofs << t << " " << energy << std::endl;
               count_t = 0;
          }

    } while(t < tfinal); // end of time loop -----------------------------
    end = std::chrono::high_resolution_clock::now();
    #pragma omp critical
    {
        std::cerr << std::endl << std::setprecision(std::numeric_limits<double>::digits10 + 1) << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::flush;
    }

    #pragma omp critical
    {
        if(sortie == convergence_CFL || sortie == convergence_Ninter)
        {
            maxenergy_ofs << dt << " " << dx << " " << maxenergy << std::endl;
        }
        else
        {
            maxenergy_ofs << omega << " " << maxenergy << std::endl;
        }
    }

    #pragma omp critical
    {
          if(sortie == convergence_CFL || sortie == convergence_Ninter)
          {
            w_ofs << dt << " " << dx << " " << *fnow << std::endl;
            energy_ofs << dt << " " << dx << " " << energy << std::endl;
          }
    }

    delete fpast;
    delete fnow;
    delete fnext;
}

#endif // EXERCICE7_CALCUL_DONNEE_H
