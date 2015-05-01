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
  int eqref, ucase, Ninter_start, Ninter_stop, nscan, ech_t = 1;
  double xL, xR, u2_max=0, tfinal = 0, leftboundaryvalue = 0, rightboundaryvalue = 0;
  double A=1., omega_start=0., omega_stop, CFL_start, CFL_stop;
  u_squared u2;
  std::vector<double>  *coeff = NULL, *u_1 = NULL, *beta = NULL;
  boundary_condition left_bc, right_bc;
  std::string nom, provisoire;
  mode choix;

  //initialisation de parametres généraux
  contexte_general(nom, Ninter_start, xL, xR, eqref, ucase);


  // open output file streams
  provisoire = nom + "_pos_usquared.dat";
  ofstream p_u_ofs(provisoire);
  p_u_ofs.precision(15);

  provisoire = nom + "_wave.dat";
  ofstream w_ofs(provisoire);
  w_ofs.precision(15);

  provisoire = nom + "_energy.dat";
  ofstream energy_ofs (provisoire);
  energy_ofs.precision(20);

  provisoire = nom + "_maxenergy.dat";
  ofstream maxenergy_ofs (provisoire);
  maxenergy_ofs.precision(15);


  //initialisation des différentes valeurs

  contexte_vitesse(ucase, xL, xR, u2);
  contexte_temporelle(tfinal, CFL_start);

  contexte_bord("left", left_bc, leftboundaryvalue, A, omega_start);
  contexte_bord("right", right_bc, rightboundaryvalue, A, omega_start);

  contexte_mode(nscan, omega_stop, Ninter_stop, CFL_stop, choix, ech_t);

   if(choix != convergence_Ninter)
   {
       Ninter_stop = Ninter_start;
   }
   if(choix != convergence_CFL)
   {
       CFL_stop = CFL_start;
   }
   if(choix != frequence)
   {
       omega_stop = omega_start;
   }

  //Sortie sur la vitesse en fonction de la position ------

    double pos, u_2 = 0;
    for(int ip = 0; ip < Ninter_stop+1; ++ip)
    {
        pos = xL + (xR - xL) / (Ninter_stop+1) * ip;
        u_2 = u2(pos);
        p_u_ofs << pos << " " << u_2 << endl;
        if(u_2 > u2_max)
        {
            u2_max = u_2;
        }
    }
    cout << "u2_max="<<u2_max << endl;
    p_u_ofs.close();
  
    if(choix == Unique)
    {
        int Npos = Ninter_start+1;
        double dx =(xR - xL) / Ninter_start;
        double dt = CFL_start * dx / sqrt(u2_max);

        cout << "dt=" << dt << endl;

        coeff = new vector<double>(Npos);
        u_1 = new vector<double>(Npos);
        beta = new vector<double>(Npos);

        //calcul des vitesse en fonctions des positions et calcul du parametre beta
        for(int ip = 0; ip < Npos; ++ip)
        {
            (*u_1)[ip] = sqrt((u2(xL + dx * ip)));
            (*beta)[ip] = (*u_1)[ip]*dt/dx;
            // coeff is beta^2
            (*coeff)[ip] = (*beta)[ip] * (*beta)[ip];
        }

        if(Npos > 1000)
        {
            simulation_par(u_1,beta,coeff,dx,dt,eqref,ucase,Npos,tfinal,omega_start,A,left_bc,leftboundaryvalue,right_bc,rightboundaryvalue,choix,w_ofs,energy_ofs,maxenergy_ofs, ech_t);
        }
        else
        {
            simulation(u_1,beta,coeff,dx,dt,eqref,ucase,Npos,tfinal,omega_start,A,left_bc,leftboundaryvalue,right_bc,rightboundaryvalue,choix,w_ofs,energy_ofs,maxenergy_ofs, ech_t);
        }

        delete coeff;
        delete u_1;
        delete beta;
    }
    else if(choix == frequence)
    {
        int Npos = Ninter_start+1;
        double dx =(xR - xL) / Ninter_start;
        double dt = CFL_start * dx / sqrt(u2_max);

        cout << "dt=" << dt << endl;

        coeff = new vector<double>(Npos);
        u_1 = new vector<double>(Npos);
        beta = new vector<double>(Npos);

        //calcul des vitesse en fonctions des positions et calcul du parametre beta
        for(int ip = 0; ip < Npos; ++ip)
        {
            (*u_1)[ip] = sqrt((u2(xL + dx * ip)));
            (*beta)[ip] = (*u_1)[ip]*dt/dx;
            // coeff is beta^2
            (*coeff)[ip] = (*beta)[ip] * (*beta)[ip];
        }


        vector<double> omega(nscan);
        for(int jscan = 0; jscan < nscan; ++jscan)
        {
            omega[jscan] = omega_start + static_cast<double>(jscan)*((omega_stop-omega_start)/static_cast<double>(nscan-1));
        }

        #pragma omp parallel for default(shared) schedule(dynamic)
        for(int jscan = 0; jscan < nscan; ++jscan)
        {
            #pragma omp critical
            {
                cout << "Scan frequency: " << jscan << endl;
            }
            simulation(u_1,beta,coeff,dx,dt,eqref,ucase,Npos,tfinal,omega[jscan],A,left_bc,leftboundaryvalue,right_bc,rightboundaryvalue,choix,w_ofs,energy_ofs,maxenergy_ofs, ech_t);
        }

        delete coeff;
        delete u_1;
        delete beta;
    }
    else if(choix == convergence_CFL)
    {
        int Npos = Ninter_start+1;
        double dx = (xR - xL) / Ninter_start;

        vector<double> dt(nscan);
        for(int jscan = 0; jscan < nscan; ++jscan)
        {
            dt[jscan] = (CFL_start + static_cast<double>(jscan)*((CFL_stop-CFL_start)/static_cast<double>(nscan-1))) * dx / sqrt(u2_max);
        }

        coeff = new vector<double>(Npos);
        u_1 = new vector<double>(Npos);
        beta = new vector<double>(Npos);

        //#pragma omp parallel for default(shared) private(coeff, u_1, beta) schedule(dynamic)
        for(int jscan = 0; jscan < nscan; ++jscan)
        {
            //calcul des vitesse en fonctions des positions et calcul du parametre beta
            for(int ip = 0; ip < Npos; ++ip)
            {
                (*u_1)[ip] = sqrt((u2(xL + dx * ip)));
                (*beta)[ip] = (*u_1)[ip]*dt[jscan]/dx;
                // coeff is beta^2
                (*coeff)[ip] = (*beta)[ip] * (*beta)[ip];
            }

            #pragma omp critical
            {
                cout << "dt=" << dt[jscan] << endl;
            }
            if(Npos > 1000)
            {
                simulation_par(u_1,beta,coeff,dx,dt[jscan],eqref,ucase,Npos,tfinal,omega_start,A,left_bc,leftboundaryvalue,right_bc,rightboundaryvalue,choix,w_ofs,energy_ofs,maxenergy_ofs, ech_t);
            }
            else
            {
                simulation(u_1,beta,coeff,dx,dt[jscan],eqref,ucase,Npos,tfinal,omega_start,A,left_bc,leftboundaryvalue,right_bc,rightboundaryvalue,choix,w_ofs,energy_ofs,maxenergy_ofs, ech_t);
            }
        }
        delete coeff;
        delete u_1;
        delete beta;
    }
    else if(choix == convergence_Ninter)
    {
        vector<double> dt(nscan);
        vector<double> dx(nscan);
        vector<int> Npos(nscan);
        for(int jscan = 0; jscan < nscan; ++jscan)
        {
            Npos[jscan] = (Ninter_start + static_cast<double>(jscan)*((Ninter_stop-Ninter_start)/static_cast<double>(nscan-1))) + 1;
            dx[jscan] = (xR - xL)/(Npos[jscan]-1);
            dt[jscan] =  CFL_start * dx[jscan] / sqrt(u2_max);
        }

        //#pragma omp parallel for default(shared) private(coeff, u_1, beta) schedule(dynamic)
        for(int jscan = 0; jscan < nscan; ++jscan)
        {
            coeff = new vector<double>(Npos[jscan]);
            u_1 = new vector<double>(Npos[jscan]);
            beta = new vector<double>(Npos[jscan]);

            //calcul des vitesse en fonctions des positions et calcul du parametre beta
            for(int ip = 0; ip < Npos[jscan]; ++ip)
            {
                (*u_1)[ip] = sqrt((u2(xL + dx[jscan] * ip)));
                (*beta)[ip] = (*u_1)[ip]*dt[jscan]/dx[jscan];
                // coeff is beta^2
                (*coeff)[ip] = (*beta)[ip] * (*beta)[ip];
            }

            #pragma omp critical
            {
                cout << "dt=" << dt[jscan] << endl;
                cout << "dx=" << dx[jscan] << endl;
                cout << "jscan=" << jscan << endl;
            }
            if(Npos[jscan] > 1000)
            {
                simulation_par(u_1,beta,coeff,dx[jscan],dt[jscan],eqref,ucase,Npos[jscan],tfinal,omega_start,A,left_bc,leftboundaryvalue,right_bc,rightboundaryvalue,choix,w_ofs,energy_ofs,maxenergy_ofs, ech_t);
            }
            else
            {
                simulation(u_1,beta,coeff,dx[jscan],dt[jscan],eqref,ucase,Npos[jscan],tfinal,omega_start,A,left_bc,leftboundaryvalue,right_bc,rightboundaryvalue,choix,w_ofs,energy_ofs,maxenergy_ofs, ech_t);
            }
            delete coeff;
            delete u_1;
            delete beta;

        }
    }

  w_ofs.close();
  energy_ofs.close();
  maxenergy_ofs.close();
}
