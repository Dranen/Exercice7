#include <iostream>
#include <vector> 
#include <fstream> 
#include <cmath> 
#include <cassert> 
#include <cstdlib> 

#include "Exercice7_io.h"
#include "Exercice7_calcul_donnee.h"

using namespace std;

//
// The main program
//
int main() 
{ 
  int Ninter,eqref, ucase, Npos;
  double xL, xR, dx, u2_max=0, tfinal = 0, dt;
  u_squared u2;

  contexte_general(Ninter, xL, XR, eqref, ucase);

  Npos = Ninter+1;
  dx=(xR - xL) / Ninter;

  std::vector<double> fpast(Npos), fnow(Npos), fnext(Npos), coeff(Npos);

  contexte_vitesse(ucase, xL, xR, dx, Npos, u2_max, u2);
  contexte_temporelle(dt, tfinal);

  
  // coeff is beta^2
  for(int ip = 0; ip < Npos; ++ip)
    coeff[ip] = u2(xL+ip*dx) * dt * dt / (dx * dx);
  
  typedef enum {fixed,free,excited,outgoing} boundary_condition;
  
  boundary_condition left_bc;
  boundary_condition right_bc;
  
  int n;
  cerr << "left boundary condition (fixed = 0, free = 1, excited = 2, outgoing = 3)? " << flush;
  cin >> n;
  
  switch(n)
    {
    case 0:
      left_bc = fixed;
      break;

    case 1:
      left_bc = free;
      break;

    case 2:
      left_bc = excited;
      break;

    case 3:
      left_bc = outgoing;
      break;
      
    default:
      cerr << "no valid left boundary condition!";
      return 1;
      
    }

  cerr << "right boundary condition (fixed = 0, free = 1, excited = 2, outgoing = 3)? " << flush;
  cin >> n;
  switch(n)
    {
    case 0:
      right_bc = fixed;
      break;

    case 1:
      right_bc = free;
      break;

    case 2:
      right_bc = excited;
      break;

    case 3:
      right_bc = outgoing;
      break;
      
    default:
      cerr << "no valid left boundary condition!";
      return 1;
      
    }

  double leftboundaryvalue = 0, rightboundaryvalue = 0;
  if(left_bc == fixed)
    {
      cerr << "for fixed left boundary condition: value of f at the left boundary? " << flush;
      cin >> leftboundaryvalue;
    }
  
  if(right_bc == fixed)
    {
      cerr << "for fixed right boundary condition: value of f at the right boundary? " << flush;
      cin >> rightboundaryvalue;
    }

  double A = 1.;
  double omega = 0;
 
  if(left_bc == excited || right_bc == excited)
    {
      cerr << "Amplitude of the excitation? " << flush;
      cin >> A;
      cerr << "Angular frequency omega of the excitation? " << flush;
      cin >> omega;
    }
  
  // open output file streams
  
  ofstream p_u_ofs("pos_usquared.dat");
  p_u_ofs.precision(15);
  
  ofstream w_ofs("wave.dat");
  w_ofs.precision(15);
  
  ofstream energy_ofs ("energy.dat");
  energy_ofs.precision(15);
  
  ofstream maxenergy_ofs ("maxenergy.dat");
  maxenergy_ofs.precision(15);
  
  // prepare frequency scan
  
  double omega_start = omega;
  double omega_stop = omega;
  double delta_omega = 0.0;
  int nscan;
  cerr << "number of frequencies to scan? " << flush;
  cin >> nscan;
    
  if (nscan > 1) // do a frequency scan only if you are interested in more than 1 frequency 
    {
      cerr << "scan from omega-start=omega to omega-end=? " << flush;
      cin >> omega_stop;
      omega_start = omega;
//       omega_stop = omega + 2.0;
      delta_omega = (omega_stop - omega_start)/nscan;
    }
  
  omega = omega_start;
  
  // frequency scan ------------------------------------------
  
  for(int jscan = 0; jscan < nscan; ++jscan) 
    {
      
      cout << "Scan frequency: " << jscan << endl;
      
      //
      // initialize the first two time slices of f
      //
      
      double xo=3e5;
      for(int ip = 0; ip < Npos; ++ip) 
	{
	  if(ucase==2){
            double x = ((double) ip)*dx-xo;
	    fnow[ip] = exp(-x*x/2.5e9);
	    fpast[ip]=fnow[ip];
	  }else{
	    fpast[ip] = 0.0;
	    fnow[ip] = 0.0;
	  }
	}
      
      //
      // output the positions and the phase velocity squared
      // for plotting purposes
      //
      
      if(jscan == 0) 
	{
	  for(int ip = 0; ip < Npos; ++ip) 
	    {
	      double pos = xL + dx * ip;
	      p_u_ofs << pos << " " << u2(pos) << std::endl;
	    }
	  p_u_ofs.close();
	}
      
      double t = 0;
      double energy = get_energy(fnow,dx);
      double maxenergy = energy;
      
      if(jscan == 0) 
	{
	  w_ofs << t << " " << fnow << endl;
	  energy_ofs << t << " " << energy << endl;
	}
      
      // time loop --------------------------------
      
      do {
        //cout << " time evolution: " << t << endl;
	// do update
	for(int ip = 1; ip < (Npos - 1); ++ip) 
	  {
	    //TODO: mettre le code du shema numerique des equations (1) et (2).
	  }
	
	//
	// apply boundary conditions
	//
	// left boundary (at xL)
	switch(left_bc) 
	  {
	    
	  case fixed:
	    fnow[0] = leftboundaryvalue;
	    fnext[0] = fnow[0];
	    break;
	    
	  case free:
	    //TODO: Mettre le code du calcul de la condition de bord gauche libre.
	    break;
	    
	  case excited:
	        fnext[0] = A * sin(omega * t);
	    break;
	    
	  case outgoing:
	    //TODO: Mettre le code du calcul de la condition de bord gauche sortante.
	    break;
	    
	  default:
	    std::cerr << "No valid boundary condition given" <<std::endl;
	    break;
	  }
	
	// right boundary (at xr)
	switch(right_bc) 
	  {
	    // note that Ninter=Npos-1
	  case fixed:
	    fnow[Ninter] = rightboundaryvalue;
	    fnext[Ninter] = fnow[Ninter];
	    break;
	    
	  case free:
	    //TODO: Mettre le code du calcul de la condition de bord droit libre.
	    break;
	    
	  case excited:
	    fnext[Ninter] = A * sin(omega * t);
	    break;
	    
	  case outgoing:
	    //TODO: Mettre le code du calcul de la condition de bord droit sortante.
	    break;
	    
	  default:
	    std::cerr << "No valid boundary condition given" <<std::endl;
	    break;
	  }
	
	t += dt;
	
	fpast = fnow;
	fnow  = fnext;
	
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
  
  w_ofs.close();
  energy_ofs.close();
  maxenergy_ofs.close();
}
