#include <iostream>
#include <cmath>


const double speedOfLight = 2.99792458e8;
const double electronCharge = 1.602e-19;
const double pi = 3.14159;
const double fineStructureConstant = (1.0/137.0);

static long double g_beta;

//For a C++ compiler 
//     Change TMath::Power(,) to pow(,)
//     Change int tammFrankCalculator() to int main()
//

//============================================================================================================================================//

//Index of Refraction of Radiator Gasses

double cf4Index = 1.00054;//CF4
//double cf4Index = 1.03; //Aerogel
//double cf4Index = 1.33; //Water

//Index of Refraction for CF4 at various lambda
double refractiveIndex(double lambda)
{
  //Conversion from nm to m
  lambda = lambda * TMath::Power(10,9);

  const double A = 0.124523, lambda_n2 = 2.61154e-4; //CF4
  //  const double A = , lambda_n2 = ; //Water
  double cf4index_lambda;

  cf4index_lambda = ((A * TMath::Power(10,-6) / (lambda_n2 - pow(lambda,-2))) + 1);  
  return (cf4index_lambda);
}

//============================================================================================================================================//  

// Particle Mass

enum Particletype
{
    ANTI_PROTON = 0,
    EP = 1,
    EN = 2,
    MUP = 3,
    MUN = 4,
    OPTICALPHOTON = 5,
    PIP = 6,
    PI0 = 7,
    PROTON = 8
};

double getParticleMass()
{
  using namespace std;

  cout << "\n 0: anti_proton \n 1: e+ \n 2: e- \n 3: mu+ \n 4: mu- \n 5: OpticalPhoton \n 6: pi+ \n 7: pi0 \n 8: proton \n\n Enter Particle Type: ";
  int inputparticle;
  cin >> inputparticle;
  Particletype particle = static_cast<Particletype>(inputparticle);

  double particlemass = 0;

  if (particle == ANTI_PROTON)
    particlemass = 1.672621777e-27;
  else if (particle == EP)
    particlemass = 9.10938291e-31;
  else if (particle == EN)
    particlemass = 9.10938291e-31;
  else if (particle == MUP)
    particlemass = 1.883531475e-28;
  else if (particle == MUN)
    particlemass = 1.883531475e-28;
  else if (particle == OPTICALPHOTON)
    particlemass = 0;
  else if (particle == PIP)
    particlemass = 2.497872e-28;
  else if (particle == PI0)
    particlemass = 2.406176e-28;
  else if (particle == PROTON)
    particlemass = 1.672621777e-27;
  else
    cout << "\n Error: wrong particle type entered \n" << endl;

  return particlemass;
}

//============================================================================================================================================//

//Particle Velocity:

void getParticleVelocity(double particlemass, double particleenergy)
{
  using namespace std;

  long double particlevelocity;

//Converts from eV to Joules
  particleenergy = particleenergy * 1.602e-19;

//Calculates beta and particle velocity
  g_beta = sqrt(1.0 - TMath::Power( (1.0/((particleenergy/(particlemass*pow(speedOfLight,2))) + 1 )) , 2 ));
  particlevelocity = g_beta * speedOfLight;
  cout << "\n Energy:   " << particleenergy << " J \n" << " Mass:     " << particlemass << " kg \n" << " Beta:     " << g_beta << "\n" << " Velocity: " << particlevelocity << " m/s \n";
}

//============================================================================================================================================//  

//Cerenkov Integration:
long double Integrand(Double_t *x, Double_t *par)

{
  //Using Sellenmeier approximation:
  return (1*(1.0/(x[0]*x[0]))*TMath::Power(TMath::Sin(TMath::ACos(1.0/(refractiveIndex(x[0])*par[0]))),2.0));
  
  //For Constant Index:
  //return (1*(1.0/(x[0]*x[0]))*TMath::Power(TMath::Sin(TMath::ACos(1.0/(cf4Index*par[0]))),2.0));
}


Double_t IntegrationSpecial()

{
  double lower = 120e-9, upper = 200e-9;

  // Create the function and wrap it
  // Extend the function...                                                                                                                                                                               
  TF1 f("myfunc",Integrand, lower,upper,1);
  f.SetParameter(0,g_beta);
  ROOT::Math::WrappedTF1 wf1(f);

  // Create the Integrator                                                                                                                                                                                         
  ROOT::Math::GaussLegendreIntegrator ig;

  // Set parameters of the integration                                                                                                                                                                             
  ig.SetFunction(wf1);

  //Set accuracy of integration here: Setting number of points even one order of magnitude higher will take probably an hour or so more.                                                                           
  ig.SetRelTolerance(0.00001);
  ig.SetNumberPoints(10000);

  // cout << ig.Integral(62e-9,1e-2) << endl;                                                                                                                                                                      
  //Range of Integration...[should always be same as function]
  return ig.Integral(lower,upper);
}

//============================================================================================================================================// 

//Cerenkov Photon Creation Rate:

double getCerenkovPhotonNumber()
{
  using namespace std;

  double cerenkovanglerad, cerenkovangledeg, photonnumber;

  // Calculates Cerenkov angle using cf4Index = 1.00054  
  cerenkovanglerad = acos(1/(cf4Index*g_beta));
  cerenkovangledeg = 180*cerenkovanglerad/pi;
  cout << "\n Cerenkov Angle: " << cerenkovanglerad << " radians \n Cerenkov Angle: " << cerenkovangledeg << " degrees" << endl;  

  // Calculates Photon Number
  photonnumber = 2*pi*fineStructureConstant*IntegrationSpecial();
  cout << " Photon Number:  " << photonnumber << " photons/m " << endl;

  return photonnumber;
} 

//============================================================================================================================================//  

int tammFrankCalculator()
{
  using namespace std;

  double particlemass = 0; 
  while(particlemass == 0)
    particlemass = getParticleMass();
  
  cout << "\n Enter particle energy in GeV: ";
  double particleenergy;
  cin >> particleenergy;

  //Converts from GeV to eV
  particleenergy = particleenergy * TMath::Power(10,9); 

  //Gets beta
  getParticleVelocity(particlemass, particleenergy);
  
  //Calculating Photon Number
  cout << "Calculating..." << endl;
  double photonnumber;
  photonnumber = getCerenkovPhotonNumber();

  cout << "\n//======================================================================//\n" << endl;

  return 0;
}
