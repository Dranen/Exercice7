#include <iostream>
#include <vector> 
#include <fstream> 
#include <cmath> 
#include <cassert> 
#include <cstdlib> 

using namespace std;

// 
// function object for the position dependent velocity
//
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

//
// calculate the energy of the wave 
//

double get_energy(const std::vector<double>& f, const double dx) 
{
  int npos = f.size();
  
  double erg = 0.;	
  for(int ip = 0; ip < (npos - 1); ++ip)
    erg +=0; //TODO: calculer l'energie.
  
  return erg * dx;
}

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

//
// The main program
//
int main() 
{
    //double pi = 3.14159265358979;
  
    int Ninter,eqref;
  cerr << "number of grid points? " << flush;
  cin >> Ninter;
  int Npos = Ninter+1;
  
  double xL, xR, hocean = 0, xocean=0,hplage=0;
  
  int ucase;
  cerr << "xL? " << flush;
  cin >> xL;
  cerr << "xR? " << flush;
  cin >> xR;
  cerr << "Quelle equation (1) ou (2)?"<<flush;
  cin >> eqref;
  cerr << "quel question? 7.2 (taper 0) / 7.3 (taper 1) / 7.4 (taper 2):" << flush;
  cin >> ucase;
  
  double dx=(xR - xL) / Ninter;
  
  double u,uL,uR;
  u_squared u2;
  double u2_max = 0;  
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
  
  
//   double u0 = u2(xL);
  
  double CFL;
  cerr << "max CFL coefficent? " << flush;
  cin >> CFL;
  
  double dt = CFL * dx / sqrt(u2_max);
  
  cout << "dt=" << dt << endl;
  double tfinal = 0;
  cerr << "tfinal? " << flush;
  cin >> tfinal;
  
  std::vector<double> fpast(Npos), fnow(Npos), fnext(Npos), coeff(Npos);
  
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
      cerr << "for fixed right boundary condition: value of f at the right boundary? " << flush;
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
