#include <iostream>
#include <vector> 
#include <fstream> 
#include <cmath> 
#include <cassert> 
#include <cstdlib> 
#include <chrono>
#include <iomanip>

#include "Exercice7_io.h"
#include "Exercice7_calcul.h"

using namespace std;

//
// The main program
//
int main() 
{ 
  int Ninter,eqref, ucase, Npos, nscan;
  double xL, xR, dx, u2_max=0, tfinal = 0, dt, leftboundaryvalue = 0, rightboundaryvalue = 0;
  double A=1., omega_start=0., omega_stop, delta_omega = 0.0, t, energy, maxenergy;
  u_squared u2;
  std::vector<double> *fpast = NULL, *fnow = NULL, *fnext = NULL, *coeff = NULL, *temp = NULL, *u_1 = NULL, *u_2 = NULL, *beta = NULL;
  boundary_condition left_bc, right_bc;
  std::string nom, provisoire;

  //initialisation de parametres généraux
  contexte_general(nom, Ninter, xL, xR, eqref, ucase);


  // open output file streams
  provisoire = nom + "_pos_usquared.dat";
  ofstream p_u_ofs(provisoire);
  p_u_ofs.precision(15);

  provisoire = nom + "_wave.dat";
  ofstream w_ofs(provisoire);
  w_ofs.precision(15);

  provisoire = nom + "_energy.dat";
  ofstream energy_ofs (provisoire);
  energy_ofs.precision(15);

  provisoire = nom + "_maxenergy.dat";
  ofstream maxenergy_ofs (provisoire);
  maxenergy_ofs.precision(15);


  //initialisation des différentes valeurs

  Npos = Ninter+1;
  dx=(xR - xL) / Ninter;

  contexte_vitesse(ucase, xL, xR, dx, Npos, u2_max, u2);
  contexte_temporelle(dt, tfinal, dx, u2_max);

  coeff = new vector<double>(Npos);
  u_1 = new vector<double>(Npos);
  u_2 = new vector<double>(Npos);
  beta = new vector<double>(Npos);

  
  contexte_bord("left", left_bc, leftboundaryvalue, A, omega_start);
  contexte_bord("right", right_bc, rightboundaryvalue, A, omega_start);

  contexte_scan_frequence(nscan, omega_stop);
  vector<double> omega(nscan);
  if(nscan > 1)
  {
      for(int jscan = 0; jscan < nscan; ++jscan)
      {
         omega[jscan] = omega_start + static_cast<double>(jscan)*((omega_stop-omega_start)/static_cast<double>(nscan-1));
      }
  }
  else
  {
      omega[0] = omega_start;
  }


  //calcul des vitesse en fonctions des positions et calcul du parametre beta
  for(int ip = 0; ip < Npos; ++ip)
  {
      (*u_2)[ip] = u2(xL + dx * ip);
      (*u_1)[ip] = sqrt(((*u_2)[ip]));
      (*beta)[ip] = (*u_1)[ip]*dt/dx;
      // coeff is beta^2
      (*coeff)[ip] = (*beta)[ip] * (*beta)[ip];
  }


  //Sortie sur la vitesse en fonction de la position ------

  double pos;
  for(int ip = 0; ip < Npos; ++ip)
    {
      pos = xL + dx * ip;
      p_u_ofs << pos << " " << u2(pos) << std::endl;
    }
  p_u_ofs.close();
  


  // frequency scan ------------------------------------------
  #pragma omp parallel for default(shared) private(energy, maxenergy, temp, fpast, fnow, fnext, t) schedule(dynamic) if(nscan > 1)
  for(int jscan = 0; jscan < nscan; ++jscan) 
  {
      fpast = new vector<double>(Npos);
      fnow = new vector<double>(Npos);
      fnext = new vector<double>(Npos);
      

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
      t = 0;
      energy = get_energy(*fnow,dx);
      maxenergy = energy;
      
      if(jscan == 0) 
      {
        w_ofs << t << " " << *fnow << endl;
        energy_ofs << t << " " << energy << endl;
      }
      

      // time loop --------------------------------
      chrono::high_resolution_clock::time_point end;
      chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
      do
      {
            //cout << " time evolution: " << t << endl;
            if(eqref == 1)
            {
                #pragma omp parallel for if(Npos > 1000 && nscan == 1)
                for(int ip = 1; ip < (Npos - 1); ++ip)
                {
                    (*fnext)[ip]=-(*fpast)[ip]+2*(1-(*coeff)[ip])*(*fnow)[ip]+(*coeff)[ip]*((*fnow)[ip-1]+(*fnow)[ip+1]);
                }
            }
            else
            {
                #pragma omp parallel for if(Npos > 1000 && nscan == 1)
                for(int ip = 1; ip < (Npos - 1); ++ip)
                {
                    (*fnext)[ip] = 0.5*(*beta)[ip]*((*u_1)[ip+1]-(*u_1)[ip-1])*((*fnow)[ip+1]-(*fnow)[ip-1])*dt/dx + (*beta)[ip]*(*beta)[ip]*((*fnow)[ip+1]-2*(*fnow)[ip]+(*fnow)[ip-1]) + 2*(*fnow)[ip] - (*fpast)[ip];
                }
            }


            // apply boundary conditions
            //
            // left boundary (at xL)
            calcul_condition_bord(true, left_bc, leftboundaryvalue, A, omega[jscan], t, fnow, fnext, 0, (*beta)[0]);
            // right boundary (at xr)
            calcul_condition_bord(false, right_bc, rightboundaryvalue, A, omega[jscan], t, fnow, fnext, Npos-1, (*beta)[Npos-1]);

            t += dt;

            temp = fpast;
            fpast = fnow;
            fnow  = fnext;
            fnext = temp;

            double energy = get_energy((*fnow),dx);

            if(energy > maxenergy)
              maxenergy = energy;

            if(nscan == 1)
              {
                w_ofs << t<< " " << *fnow << endl;
                energy_ofs << t << " " << energy << endl;
              }

      } while(t < tfinal); // end of time loop -----------------------------
      end = std::chrono::high_resolution_clock::now();
      #pragma omp critical
      {
          cout << "Scan frequency: " << jscan << endl;
          cerr << endl << setprecision(std::numeric_limits<double>::digits10 + 1) << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms" << endl;

          maxenergy_ofs << omega[jscan] << " " << maxenergy << endl;
      }

      delete fpast;
      delete fnow;
      delete fnext;
      
    } // end of frequency scan loop -------------------

  w_ofs.close();
  energy_ofs.close();
  maxenergy_ofs.close();

  delete coeff;
  delete u_1;
  delete u_2;
  delete beta;
}
