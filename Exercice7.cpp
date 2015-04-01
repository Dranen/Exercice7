#include <iostream>
#include <vector> 
#include <fstream> 
#include <cmath> 
#include <cassert> 
#include <cstdlib> 

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
  double A=1., omega=0., omega_start, omega_stop, delta_omega = 0.0, t, energy, maxenergy;
  u_squared u2;
  std::vector<double> *fpast, *fnow, *fnext, *coeff, *temp;
  boundary_condition left_bc, right_bc;

  // open output file streams

  ofstream p_u_ofs("pos_usquared.dat");
  p_u_ofs.precision(15);

  ofstream w_ofs("wave.dat");
  w_ofs.precision(15);

  ofstream energy_ofs ("energy.dat");
  energy_ofs.precision(15);

  ofstream maxenergy_ofs ("maxenergy.dat");
  maxenergy_ofs.precision(15);

  //initialisation des différentes valeurs

  contexte_general(Ninter, xL, XR, eqref, ucase);

  Npos = Ninter+1;
  dx=(xR - xL) / Ninter;

  contexte_vitesse(ucase, xL, xR, dx, Npos, u2_max, u2);
  contexte_temporelle(dt, tfinal);

  fpast = new vector<double>(Npos);
  fnow = new vector<double>(Npos);
  fnext = new vector<double>(Npos);
  coeff = new vector<double>(Npos);
  
  // coeff is beta^2
  for(int ip = 0; ip < Npos; ++ip)
    coeff->[ip] = u2(xL+ip*dx) * dt * dt / (dx * dx);
  
  contexte_bord("left", left_bc, leftboundaryvalue, A, omega);
  contexte_bord("right", right_bc, rightboundaryvalue, A, omega);

  contexte_scan_frequence(nscan, omega, omega_start, omega_stop, delta_omega);



  //Sortie sur la vitesse en fonction de la position ------

  double pos;
  for(int ip = 0; ip < Npos; ++ip)
    {
      pos = xL + dx * ip;
      p_u_ofs << pos << " " << u2(pos) << std::endl;
    }
  p_u_ofs.close();
  


  // frequency scan ------------------------------------------

  for(int jscan = 0; jscan < nscan; ++jscan) 
  {   cout << "Scan frequency: " << jscan << endl;
      

       // initialize the first two time slices of f
       double xo=3e5;
       double x;
       for(int ip = 0; ip < Npos; ++ip)
       {
          if(ucase==2){
            x = ((double) ip)*dx-xo;
            fnow->[ip] = exp(-x*x/2.5e9);
            fpast->[ip]=fnow->[ip];
          }else{
            fpast->[ip] = 0.0;
            fnow->[ip] = 0.0;
          }
       }


      // output the positions and the phase velocity squared for plotting purposes
      t = 0;
      energy = get_energy(fnow,dx);
      maxenergy = energy;
      
      if(jscan == 0) 
      {
        w_ofs << t << " " << fnow << endl;
        energy_ofs << t << " " << energy << endl;
      }
      

      // time loop --------------------------------
      do
      {
            //cout << " time evolution: " << t << endl;
            if(eqref == 1)
            {
                for(int ip = 1; ip < (Npos - 1); ++ip)
                {
                    //TODO: mettre le code du shema numerique des equations (1).
                }
            }
            else if(eqref == 2)
            {
                for(int ip = 1; ip < (Npos - 1); ++ip)
                {
                    //TODO: mettre le code du shema numerique des equations (2).
                }
            }


            // apply boundary conditions
            //
            // left boundary (at xL)
            calcul_condition_bord(left_bc, leftboundaryvalue, A, omega, t, fnow, fnext, 0);
            // right boundary (at xr)
            calcul_condition_bord(right_bc, rightboundaryvalue, A, omega, t, fnow, fnext, Ninter);

            t += dt;

            temp = fpast;
            fpast = fnow;
            fnow  = fnext;
            fnext = temp;

            double energy = get_energy(fnow,dx);

            if(energy > maxenergy)
              maxenergy = energy;

            if(nscan == 1)
              {
                w_ofs << t<< " " << fnow << endl;
                energy_ofs << t << " " << energy << endl;
              }

      } while(t < tfinal); // end of time loop -----------------------------
      
      maxenergy_ofs << omega << " " << maxenergy << endl;
      
      omega += delta_omega;
      
    } // end of frequency scan loop -------------------

  delete fpast;
  delete fnow;
  delete fnext;
  
  w_ofs.close();
  energy_ofs.close();
  maxenergy_ofs.close();
}
