//PROBLEM
//Two compartment model + extravascular absportion + nonlinear clearance

$PARAM
TVCL = 1, TVVC = 20, TVQ =2, TVVP = 10, TVKA = 1, TVVMAX = 10, TVKM   =  2

$CMT  @annotated 
// similar to $MODEL in NONMEM
// if @annotated is not used, simpley declare compartment names with a space

EV1    : First extravascular compartment (mass)
CENT   : Central compartment (mass)
PERIPH : Peripheral compartment (mass) 

$MAIN
//akin to $PK block in NONMEM
  
double CL   = TVCL*exp(ETA(1));
double VC   = TVVC*exp(ETA(2));
double Q    = TVQ;
double VP   = TVVP;
double KA   = TVKA;
double VMAX = TVVMAX;
double KM   = TVKM;

  
$GLOBAL 
// declare variables you want to either use further in computation
// OR derive variables you want calculated from computations

//int TEST = 1; // declaring variables needs ";" at the end
#define CP (CENT/VC) // # define does not require ";"
#define CT (PERIPH/VP)
#define CLNL (VMAX/(KM+CP))

$OMEGA @block
// generates a block matrix with the formulation (1,1), (1,2), (2,2)....
// @correlation will specify values to be corrleations instead as variances

0.1 0.02 0.3

$SIGMA 
// sigma variances go here

0.01

$ODE
//similar to $DES in NONMEM
//differential equations go here
  
dxdt_EV1 = -KA*EV1;
dxdt_CENT = KA*EV1 - (CL+CLNL+Q)*CP  + Q*CT;
dxdt_PERIPH = Q*CP - Q*CT;

$TABLE
// Use $TABLE to interact with parameters, compartment values, 
// and other user-defined variables after the system advances to the next time
// double CP = (CENT/VC);

double DV = CP * (1 + EPS(1));

$CAPTURE
//This is a block to identify variables that should be captured in the simulated output

CP, DV