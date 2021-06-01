//===================================================================================
// Joao C.C. Henriques 2021 | IST Ocean Energy Systems and Technology Course
// email: joaochenriques@tecnico.ulisboa.pt
//===================================================================================

model PTOsim_2021
  
  //=================================================================================
  // Constants (prefix "c_" )
  //=================================================================================

  constant Real c_pi( unit = "Pa" ) = Modelica.Constants.pi;
  constant Real c_gamma( unit = "1" ) = 1.4 "nitrogen specific heat ratio";
  constant Real c_patm( unit = "Pa" ) = 101325.0 "atmospheric pressure";
  constant Real c_beta( unit = "Pa" ) = 1.0E9 "atmospheric pressure [Pa]";
  constant Real c_rho_oil( unit = "kg/m3" ) = 850.0 "hydraulic fluid density";
  constant Real c_floater_mass( unit = "kg" ) = 140000.0 "floater mass";

  constant Real c_Racc( unit = "m" ) = 0.05 "acc. inlet/outlet hose radius";
  constant Real c_Cd( unit = "1" ) = 0.65 "acc. inlet/outlet discharge coef.";
  constant Real c_Inertia( unit = "kgm2" ) = 2.0 "motor/generator inertia";

  constant Real c_p0A( unit = "Pa" ) = 80 * c_patm "initial pressure chamber A";
  constant Real c_p0D( unit = "Pa" ) = 80 * c_patm "initial pressure chamber D";

  // floater excitation force
  //-------------------------
  // First frequency
  constant Real c_Aexc1( unit = "N" ) = 6000E3 "exc. force amplitude";
  constant Real c_Texc1( unit = "s" ) = 8.0 "exc. force period";
  constant Real c_omega1_exc( unit = "rad/s" ) = 2.0*c_pi/c_Texc1 "exc. force freq.";

  // Second frequency
  constant Real c_Aexc2( unit = "N" ) = 2500E3 "exc. force amplitude";
  constant Real c_Texc2( unit = "s" ) = 13.0 "exc. force period";
  constant Real c_omega2_exc( unit = "rad/s" ) = 2.0*c_pi/c_Texc2 "exc. force freq.";
  constant Real c_phase2_exc( unit = "rad" ) = c_pi * 0.37 "exc. force freq.";
  
  // Third frequency
  constant Real c_Aexc3( unit = "N" ) = 3500E3 "exc. force amplitude";
  constant Real c_Texc3( unit = "s" ) = 11.0 "exc. force period";
  constant Real c_omega3_exc( unit = "rad/s" ) = 2.0*c_pi/c_Texc3 "exc. force freq.";
  constant Real c_phase3_exc( unit = "rad" ) = c_pi * 0.67 "exc. force freq.";
  
  //=================================================================================
  // Parameters (prefix "p_" )
  //=================================================================================

  // hydraulic cylinder
  parameter Real p_Rp( unit = "m" ) = 0.4 "cylinder radius";
  parameter Real p_lp( unit = "m" ) = 16.0 "cylinder height";

  // accumulators data
  // p_Vacc should be a positive constant times p_Vacc
  parameter Real p_Vacc( unit = "m3" ) = 2.0 * p_Vh0 "total acc. volume";
  parameter Real p_V0( unit = "m3" ) = 0.5 * p_Vacc "initial hydraulic fluid volume";
  parameter Real p_p0H( unit = "Pa" ) = 1.1 * c_p0A "initial pressure chamber A";
  parameter Real p_p0L( unit = "Pa" ) = 0.9 * c_p0D "initial pressure chamber D";

  // PTO system
  parameter Real p_Dv( unit = "m3/rad" ) = 0.0008 "hydraulic motor volumetric disp.";
  parameter Real p_a( unit = "Nm/(rad2/s2)" ) = 0.1 "generator control parameter";

  //=================================================================================
  // do not change these parameters
  //=================================================================================

  parameter Real p_Vg0( unit = "m3" ) = p_Vacc - p_V0 "initial gas volume";
  parameter Real p_Aacc( unit = "m2" ) = c_pi*c_Racc^2 "acc. inlet/outlet area";
  parameter Real p_Ap( unit = "m2" ) = c_pi * p_Rp^2 "cylinder/piston area";
  parameter Real p_Vh0( unit = "m3" ) = 0.5 * p_lp * p_Ap "cylinder chamber volume";

  //=================================================================================
  // Variables (no prefix)
  // Equations of differential equations must have initial conditions ( start)
  //=================================================================================

  // algebraic equations
  Real Fexc( unit = "N" ) "excitation force";
  Real smoothing( unit = "1" ) "initial smoothing ramp";
  Real y( start = 0, unit = "m", fixed = true ) "piston position";
  Real v( start = 0, unit = "m/s", fixed = true ) "piston velocity";
  
  // Control manifold (Integers have no units)
  Integer sAB "connection A->B open [1/0]";
  Integer tilde_sAB "connection A->B closed [1/0]";
  Integer sCA "connection C->A open [1/0]";
  Integer sDB"connection D->B open [1/0]";
  Integer tilde_sDB "connection D->B closed [1/0]";
  Integer sCD "connection C->D open [1/0]";
  
  // cylinder chamber A
  Real VA( unit = "m3" ) "volume"; 
  Real mA( unit = "kg" ) "hydraulic fluid mass"; 
  Real wVA( unit = "kg/s" ) "piston pump mass flow rate";
  Real wA( unit = "kg/s" ) "inlet/outlet mass flow rate"; 
  Real rhoA( start = c_rho_oil, unit = "kg/m3", fixed = true ) "hyd. fluid density";
  Real pA( start = p_p0H, unit = "Pa", fixed = true ) "hydraulic fluid pressure";
  
  // cylinder chamber D
  Real VD( unit = "m3" ) "volume"; 
  Real mD( unit = "kg" ) "hydraulic fluid mass"; 
  Real wVD( unit = "kg/s" ) "piston pump mass flow rate";
  Real wD( unit = "kg/s" ) "inlet/outlet mass flow rate"; 
  Real rhoD( start = c_rho_oil, unit = "kg/m3", fixed = true ) "hyd. fluid density";
  Real pD( start = c_p0D, unit = "Pa", fixed = true ) "hydraulic fluid pressure";

  // PTO system
  Real wM( unit = "m3/s" ) "hydraulic motor flow rate";
  Real TM( unit = "Nm" ) "hydraulic motor torque"; 
  Real PM( unit = "W" ) "hydraulic motor power"; 
  Real TG( unit = "Nm" ) "generator electromagnetic shaft torque"; 
  Real PG( unit = "W" ) "generator electromagnetic shaft power";
  Real Omega( start = 100, unit = "rad/s", fixed = true ) "PTO rotational speed";
  Real EG( start = 0.0, unit = "J", fixed = true ) "generator electromagnetic shaft power";
  
  // High-pressure accumulator
  Real wH( unit = "kg/s" ) "inlet/outlet mass flow rate";
  Real mH( unit = "kg" ) "hydraulic fluid mass";
  Real VH( start = p_V0, unit = "m3", fixed = true ) "hydraulic fluid volume";
  Real rhoH( start = c_rho_oil, unit = "kg/m3", fixed = true ) "hyd. fluid density";
  Real VgH( start = p_Vg0, unit = "m3/s", fixed = true ) "gas volume";
  Real pH( start = p_p0H, unit = "Pa", fixed = true ) "accumulator pressure";
  
  // Low-pressure accumulator
  Real wL( unit = "kg/s" ) "inlet/outlet mass flow rate";
  Real mL( unit = "kg" ) "hydraulic fluid mass";
  Real VL( start = p_V0, unit = "m3", fixed = true ) "hydraulic fluid volume";
  Real rhoL( start = c_rho_oil, unit = "kg/m3", fixed = true ) "hyd. fluid density";
  Real VgL( start = p_Vg0, unit = "m3", fixed = true ) "gas volume";
  Real pL( start = p_p0L, unit = "Pa", fixed = true ) "accumulator pressure";
  
  // numerical model checking
  Real mT( unit = "kg" ) "total system mass";

equation

  //=================================================================================
  // control manifold chamber A
  if (pA > pH) then
    sAB = 1;
    tilde_sAB = 0;
    sCA = 0;
  elseif (pA < pL) then
    sAB = 0;
    tilde_sAB = 1;
    sCA = 1;
  else
    sAB = 0;
    tilde_sAB = 1;
    sCA = 0;
  end if;
  
  //=================================================================================
  // control manifold chamber D
  if (pD > pH) then
    sDB = 1;
    tilde_sDB = 0;
    sCD = 0;
  elseif (pD < pL) then
    sDB = 0;
    tilde_sDB = 1;
    sCD = 1;
  else
    sDB = 0;
    tilde_sDB = 1;
    sCD = 0;
  end if;
  
  //=================================================================================
  // Cylinder/piston motion
  smoothing = min( time / ( 2.0 * c_Texc1 ), 1.0 ); // 2 periods smoothing ramp 
  
  Fexc = c_Aexc1 * sin( c_omega1_exc * time ) 
       + c_Aexc2 * sin( c_omega2_exc * time + c_phase2_exc ) 
       + c_Aexc3 * sin( c_omega3_exc * time + c_phase3_exc );

  c_floater_mass * der(v) = smoothing * Fexc - p_Ap * (pA - pD);
  der(y) = v;

  //=================================================================================
  // Cylinder's chamber A
  VA = p_Vh0 - p_Ap * y;
  mA = rhoA * VA;  

  wVA = -rhoA * p_Ap * v;
  wA = +sAB * c_Cd * p_Aacc * sqrt( 2 * rhoA * max(pA - pH,0.0) ) 
       -sCA * c_Cd * p_Aacc * sqrt( 2 * rhoL * max(pL - pA,0.0) );
         
  der(rhoA) / rhoA + wVA / mA = -wA / mA;
  der(pA) / c_beta = der(rhoA) / rhoA;
  
  //=================================================================================
  // Cylinder's chamber D
  VD = p_Vh0 + p_Ap * y;
  mD = rhoD * VD;

  wVD = rhoD * p_Ap * v;
  wD = +sDB * c_Cd * p_Aacc * sqrt( 2 * rhoD * max(pD - pH,0.0) ) 
       -sCD * c_Cd * p_Aacc * sqrt( 2 * rhoL * max(pL - pD,0.0) );
         
  der(rhoD) / rhoD + wVD / mD = -wD / mD;
  der(pD) / c_beta = der(rhoD) / rhoD;
  
  //=================================================================================
  // PTO system
  wM = ( sAB * rhoA + sDB * rhoD + tilde_sAB * tilde_sDB * rhoH ) * p_Dv * Omega;
  
  TM = p_Dv * max( pH - pL, 0.0 ); // non-reversible motor
  PM = TM * Omega;
  TG = p_a * Omega^2;
  PG = TG * Omega;
  
  c_Inertia * der(Omega) = TM - TG;
  der(EG) = PG;
    
  //=================================================================================
  // high-pressure accumulator
  wH = sAB * wA + sDB * wD - wM;
  mH = rhoH * VH;
  
  der(mH) = wH;
  der(VgH) + der(VH) = 0.0;  
  der(pH) / pH + c_gamma * der(VgH) / VgH = 0.0;
  der(pH) / c_beta = der(rhoH) / rhoH;

  //=================================================================================
  // low-pressure accumulator
  wL = wM + ( sCA * wA + sCD * wD ); // NOTE: wA and wD have negative signs
  mL = rhoL * VL;
  
  der(mL) = wL;
  der(VgL) + der(VL) = 0.0;  
  der(pL) / pL + c_gamma * der(VgL) / VgL = 0.0;
  der(pL) / c_beta = der(rhoL) / rhoL;
  
  //=================================================================================
  // check total system mass (must be constant a time constant)
  mT = mA + mD + mH + mL;
  
end PTOsim_2021;
