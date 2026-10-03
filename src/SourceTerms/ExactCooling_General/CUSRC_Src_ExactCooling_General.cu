#include "CUFLU.h"
#include "Global.h"

#ifdef EXACT_COOLING_GENERAL

// external functions and GPU-related set-up
#ifdef __CUDACC__

#include "CUDA_CheckError.h"
#include "CUDA_ConstMemory.h"
#include "CUFLU_Shared_FluUtility.cu"

#endif // #ifdef __CUDACC__


// local function prototypes
#ifndef __CUDACC__

void Src_SetAuxArray_ExactCooling_General( double [], int [] );
void Src_SetConstMemory_ExactCooling_General( const double AuxArray_Flt[], const int AuxArray_Int[],
                                              double *&DevPtr_Flt, int *&DevPtr_Int );
void Src_SetCPUFunc_ExactCooling_General( SrcFunc_t & );
#ifdef GPU
void Src_SetGPUFunc_ExactCooling_General( SrcFunc_t & );
#endif
void Src_WorkBeforeMajorFunc_ExactCooling_General( const int lv, const double TimeNew, const double TimeOld, const double dt,
                                                    double AuxArray_Flt[], int AuxArray_Int[] );
void Src_End_ExactCooling_General();
double Mis_GetTimeStep_ExactCooling_General( const int lv, const double dTime_dt );
#endif // #ifdef __CUDACC__

GPU_DEVICE static
double ExactCooling_GetLambda( const double Temp, const int k, const double Tks[], const double Lks[], const double aks[] );
GPU_DEVICE static
double ExactCooling_GetTcool( const double Temp, const double n, const double n_e, const double n_H, const double Lambda, const double kb, const double gamma );
GPU_DEVICE static
double TEF_PiecewisePowerLaw( const double Temp, const int k, const int N, const double Tks[], const double Lks[], const double aks[], double Yks[] );
GPU_DEVICE static
double TEF_inverse_PiecewisePowerLaw( const double Y, const int N, const double Yks[], const double Lks[], const double aks[], const double Tks[] );

/********************************************************
1. General exact-cooling source term
   --> Enabled by the compilation option "EXACT_COOLING_GENERAL" and the runtime option "SRC_EXACTCOOLING_GENERAL"

2. This file is shared by both CPU and GPU
   
   CUSRC_Src_ExactCooling_General.cu -> CPU_Src_ExactCooling_General.cpp

3. Four steps are required to implement a source term

   I.   Set auxiliary arrays
   II.  Implement the source-term function
   III. [Optional] Add the work to be done every time
        before calling the major source-term function
   IV.  Set initialization functions

4. The source-term function must be thread-safe and
   not use any global variable
********************************************************/



// =======================
// I. Set auxiliary arrays
// =======================

//-------------------------------------------------------------------------------------------------------
// Function    :  Src_SetAuxArray_ExactCooling_General
// Description :  Set the auxiliary arrays AuxArray_Flt/Int[]
//
// Note        :  1. Invoked by Src_Init_ExactCooling_General()
//                2. AuxArray_Flt/Int[] have the size of SRC_NAUX_EXACTCOOLING_GENERAL defined in Macro.h (default = 20)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//                4. Change the values to implement your piecewise power law cooling funciton
//
// Parameter   :  AuxArray_Flt/Int : Floating-point/Integer arrays to be filled up
//
// Return      :  AuxArray_Flt/Int[]
//-------------------------------------------------------------------------------------------------------
#ifndef __CUDACC__
void Src_SetAuxArray_ExactCooling_General( double AuxArray_Flt[], int AuxArray_Int[] )
{
   // ==========================
   // Gas settings
   // ==========================
   const double X_H   = 0.7;                                                         // hydrogen mass fraction
   const double Z     = 0.018;                                                       // metal mass fraction

   // mean molecular weights (for full ionization regime)
   const double mu_e  = 2.0 / ( 1.0 + X_H );                                        // mean electron molecular weight
   const double mu_H  = 1.0 / X_H;                                                  // mean proton molecular weight
   const double mu    = 1.0 / ( 2.0 * X_H + 0.75 * ( 1.0 - X_H - Z ) + 0.5 * Z );   // mean total molecularweight
   
   // ===========================================================
   // Piecewise power law cooling function
   // Power law :  Lambda(T) = Lk * (T/Tk)^ak for Tk <= T < Tk+1
   // ===========================================================
   // example : 4 temperature points (N = 4) with 3 power laws
   const int    N     = 4;                                                          // number of temperature points
   
   // Temperature points (K)
   const double T1    = 1.0e4;
   const double T2    = 1.0e5;
   const double T3    = pow(10.0, 7.4);
   const double T4    = 1.0e9;

#  ifdef GAMER_DEBUG
   if ( T1 <= 0.0 )
      Aux_Error( ERROR_INFO, "T1 (%24.16e) <= 0.0 !!\n", T1 );
#  endif   

   // Power law slopes
   const double a1    =  0.74;
   const double a2    = -0.70;
   const double a3    =  0.50;

   // Normalizations (erg cm^3 s^-1)
   const double L1    =  9.93e-23;
   const double L2    =  5.51e-22;
   const double L3    =  1.15e-23;


   // store in the aux arrays
   AuxArray_Flt[0] = X_H;                                                                                                  
   AuxArray_Flt[1] = Z;       

   AuxArray_Flt[2] = mu_e;                                                                    
   AuxArray_Flt[3] = mu_H;                                                                                     
   AuxArray_Flt[4] = mu;   

   AuxArray_Flt[5] = Const_mp / UNIT_M;           // proton mass (g)
   AuxArray_Flt[6] = Const_kB / UNIT_E;           // Boltzmann constant (erg/K)
   AuxArray_Flt[7] = GAMMA;                       // ratio of specific heats

   AuxArray_Int[0] = N;            

   // Temperature points (K)
   AuxArray_Flt[8]  = T1;             
   AuxArray_Flt[9]  = T2;             
   AuxArray_Flt[10] = T3;   
   AuxArray_Flt[11] = T4;             

   // Power law slopes
   AuxArray_Flt[12] = a1;         
   AuxArray_Flt[13] = a2;              
   AuxArray_Flt[14] = a3;             

   // Normalizations (erg cm^3 s^-1)
   AuxArray_Flt[15] = L1 / ( UNIT_E * UNIT_L * UNIT_L * UNIT_L / UNIT_T );         
   AuxArray_Flt[16] = L2 / ( UNIT_E * UNIT_L * UNIT_L * UNIT_L / UNIT_T );       
   AuxArray_Flt[17] = L3 / ( UNIT_E * UNIT_L * UNIT_L * UNIT_L / UNIT_T );         

} // FUNCTION : Src_SetAuxArray_ExactCooling_General
#endif // #ifndef __CUDACC__



// ======================================
// II. Implement the source-term function
// ======================================

//-------------------------------------------------------------------------------------------------------
// Function    :  Src_ExactCooling_General
// Description :  Major source-term function
//
// Note        :  1. Invoked by CPU/GPU_SrcSolver_IterateAllCells()
//                2. See Src_SetAuxArray_ExactCooling_General() for the values stored in AuxArray_Flt/Int[]
//                3. Follow Townsend (2009, ApJS, 181, 391) to implement the exact integration of the piecewise power law cooling function
//                4. Shared by both CPU and GPU
//                5. Please change the length of Yks to your N value (line 242)
//
// Parameter   :  fluid             : Fluid array storing both the input and updated values
//                                    --> Including both active and passive variables
//                B                 : Cell-centered magnetic field
//                SrcTerms          : Structure storing all source-term variables
//                dt                : Time interval to advance solution
//                dh                : Grid size
//                x/y/z             : Target physical coordinates
//                TimeNew           : Target physical time to reach
//                TimeOld           : Physical time before update
//                                    --> This function updates physical time from TimeOld to TimeNew
//                MinDens/Pres/Eint : Density, pressure, and internal energy floors
//                PassiveFloor      : Bitwise flag to specify the passive scalars to be floored
//                EoS               : EoS object
//                AuxArray_*        : Auxiliary arrays (see the Note above)
//
// Return      :  fluid[]
//-----------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static void Src_ExactCooling_General( real fluid[], const real B[],
                                      const SrcTerms_t *SrcTerms, const real dt, const real dh,
                                      const double x, const double y, const double z,
                                      const double TimeNew, const double TimeOld,
                                      const real MinDens, const real MinPres, const real MinEint, const long PassiveFloor,
                                      const EoS_t *EoS, const double AuxArray_Flt[], const int AuxArray_Int[] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );
   if ( AuxArray_Int == NULL )   printf( "ERROR : AuxArray_Int == NULL in %s !!\n", __FUNCTION__ );
#  endif
   
   // ==========================
   // (1) read parameters
   // ==========================
   // mean molecular weights (for full ionization regime)
   const double mu_e  = AuxArray_Flt[2];                             // mean electron molecular weight
   const double mu_H  = AuxArray_Flt[3];                             // mean proton molecular weight
   const double mu    = AuxArray_Flt[4];                             // mean total molecular weight

   // Constants
   const double mp    = AuxArray_Flt[5];                             // proton mass (g)
   const double kB    = AuxArray_Flt[6];                             // Boltzmann constant (erg/K)
   const double gamma = AuxArray_Flt[7];                             // ratio of specific heats

   // ===========================================================
   
   // Piecewise power law cooling function
   const int    N     = AuxArray_Int[0];                             // number of temperature points
   
   const double *Tks  = &AuxArray_Flt[8];                            // Temperature points (K)
   const double *aks  = &AuxArray_Flt[8+N];                          // Power law slopes
   const double *Lks  = &AuxArray_Flt[8+2*N-1];                      // Normalizations 

   const double TN    = Tks[N-1];                                    // maximum temperature (K)
   const double LN    = Lks[N-2] * pow( TN / Tks[N-2], aks[N-2] );   


   // ==========================
   // (2) calculate parameters
   // ==========================
   // mass density
   const double rho  = fluid[DENS];           // g cm^-3                  

   // number densities
   const double n_H  = rho / (mu_H * mp);     // cm^-3       
   const double n_e  = rho / (mu_e * mp);     // cm^-3                   
   const double n    = rho / (mu   * mp);     // cm^-3                  

// please change the length of Yks to your N value
   double T_now, T_new, Eint, Eintf, Emag, Pres, Lambda_T, Y_T, t_cool, Y, Yks[4];
   int    k = -1;

   // ===========================================================
   // (3) calculate the current temperature and internal energy
   // ===========================================================
   // current internal energy
   const bool CheckMinEint_No = false;
   Emag  = 0.0;
#  ifdef __CUDACC__
   Eint  = Hydro_Con2Eint( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY],
                           CheckMinEint_No, NULL_REAL, PassiveFloor, Emag, 
                           EoS->GuessHTilde_FuncPtr, EoS->HTilde2Temp_FuncPtr, 
                           EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int,
                           EoS->Table );
#  else
   Eint  = Hydro_Con2Eint( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY],
                           CheckMinEint_No, NULL_REAL, PassiveFloor, Emag,
                           EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                           EoS_AuxArray_Flt, EoS_AuxArray_Int,
                           h_EoS_Table );
#  endif

   // current temperature
#  ifdef __CUDACC__
   T_now = EoS->DensEint2Temp_FuncPtr( fluid[DENS], Eint, NULL, EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int, EoS->Table  );
#  else
   T_now = EoS_DensEint2Temp_CPUPtr  ( fluid[DENS], Eint, NULL, EoS_AuxArray_Flt,        EoS_AuxArray_Int,        h_EoS_Table );
#  endif

   // =====================================================================
   // (4) calculate the new temperature and internal energy after cooling
   // =====================================================================
   // new temperature
   if ( T_now >= Tks[0] && T_now <= Tks[N-1] )
   {
      // k : the index of the temperature interval where T_now falls into
      for (int i = N-2; i >= 0; i--)
      {
         if ( T_now >= Tks[i] )
         {
            k = i+1;
            break;
         }
      }

      // cooling efficiency
      Lambda_T = ExactCooling_GetLambda( T_now, k, Tks, Lks, aks );          // erg cm^3 s^-1

      // Y(T)
      Y_T      = TEF_PiecewisePowerLaw( T_now, k, N, Tks, Lks, aks, Yks );

      // new temperature
      t_cool   = ExactCooling_GetTcool( T_now, n, n_e, n_H, Lambda_T, kB, gamma );  // s
      Y        = Y_T + T_now / TN * LN / Lambda_T * dt / t_cool;
      
      T_new    = TEF_inverse_PiecewisePowerLaw( Y, N, Yks, Lks, aks, Tks );
      
      if (T_new < Tks[0])
         T_new = Tks[0];

   }
   else
   {
#     ifdef GAMER_DEBUG
      printf( "ERROR : T_now = %24.16e is out of range (min: %24.16e, max: %24.16e) at TimeNew = %24.16e !!\n",
              T_now, Tks[0], TN, TimeNew );
#     endif
      
      fluid[ENGY] = NAN;
      return;
   }

   // new internal energy
#  ifdef __CUDACC__
   Pres  = EoS->DensTemp2Pres_FuncPtr( fluid[DENS], T_new, NULL, EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int, EoS->Table  );
   Eintf = EoS->DensPres2Eint_FuncPtr( fluid[DENS], Pres , NULL, EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int, EoS->Table  );
#  else
   Pres  = EoS_DensTemp2Pres_CPUPtr  ( fluid[DENS], T_new, NULL, EoS_AuxArray_Flt,        EoS_AuxArray_Int,        h_EoS_Table );
   Eintf = EoS_DensPres2Eint_CPUPtr  ( fluid[DENS], Pres , NULL, EoS_AuxArray_Flt,        EoS_AuxArray_Int,        h_EoS_Table );
#  endif

   // convert back
   fluid[ENGY] = Hydro_ConEint2Etot  ( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], Eintf, Emag );


} // FUNCTION : Src_ExactCooling_General


//-------------------------------------------------------------------------------------------------------
// Function    :  ExactCooling_GetLambda
// Description :  Calculate the cooling efficiency for the piecewise power law cooling function
//
// Note        :  1. Lambda(T) = Lk * (T/Tk)^ak for Tk <= T < Tk+1
//
// Parameter   :  Temp : Temperature in Kelvin
//                k    : Temperature interval index
//                Tks  : Array of temperature intervals
//                Lks  : Array of cooling rates
//                aks  : Array of power law slopes
//
// Return      :  Lambda(Temp)
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE static
double ExactCooling_GetLambda( const double Temp, const int k, const double Tks[], const double Lks[], const double aks[] )
{
   return Lks[k-1] * pow( Temp / Tks[k-1], aks[k-1]);   // erg cm^3 s^-1
} // FUNCTION : ExactCooling_GetLambda


//-------------------------------------------------------------------------------------------------------
// Function    :  ExactCooling_GetTcool
// Description :  Calculate the single point cooling time for the piecewise power law cooling function
//
// Note        :  t_cool = n * kB * T / ( (gamma - 1) * n_e * n_H * Lambda(T) ) 
//
// Parameter   :  Temp : Temperature in Kelvin
//                n    : Number density in cm^-3
//                n_e  : Electron number density in cm^-3
//                n_H  : Proton number density in cm^-3
//                Lambda : Cooling efficiency in erg cm^3 s^-1
//
// Return      :  t_cool
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE static
double ExactCooling_GetTcool( const double Temp, const double n, const double n_e, const double n_H, const double Lambda, const double kB, const double gamma )
{
   return n * kB * Temp / ( (gamma - 1.0) * n_e * n_H * Lambda );  // s 
} // FUNCTION : ExactCooling_GetTcool


//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_GetTimeStep_ExactCooling_General
// Description :  Estimate the evolution time-step constrained by the user's exact-cooling source term
//
// Note        :  1. Scan all cells at the target refinement level and find the minimum cooling time.
//                2. The cooling time is calculated using exactly the same gas settings and
//                   piecewise power-law cooling function as Src_ExactCooling_General().
//                3. The returned timestep is
//
//                      dt = SRC_EXACTCOOLING_GENERAL_DT * min(tcool)
//
//                   where tcool and dt is in GAMER code time units (s).
//                4. The minimum is first found over OpenMP threads and then over MPI processes.
//                5. Cells already below the cooling temperature floor are ignored.
//
// Parameter   :  lv       : Target refinement level
//                dTime_dt : dTime/dt (== 1.0 if COMOVING is off)
//
// Return      :  dt in GAMER code time units (s)
//-------------------------------------------------------------------------------------------------------
#ifndef __CUDACC__
double Mis_GetTimeStep_ExactCooling_General( const int lv, const double dTime_dt )
{
   (void) dTime_dt;   // 1.0 for non-comoving

   // ==========================
   // (1) read parameters
   // ==========================

   const double mu_e  = Src_ExactCooling_General_AuxArray_Flt[2];          // mean electron molecular weight
   const double mu_H  = Src_ExactCooling_General_AuxArray_Flt[3];          // mean proton molecular weight
   const double mu    = Src_ExactCooling_General_AuxArray_Flt[4];          // mean total molecular weight

   const double mp    = Src_ExactCooling_General_AuxArray_Flt[5];          // proton mass (g)
   const double kB    = Src_ExactCooling_General_AuxArray_Flt[6];          // Boltzmann constant (erg/K)
   const double gamma = Src_ExactCooling_General_AuxArray_Flt[7];          // ratio of specific heats
   
   const int    N     = Src_ExactCooling_General_AuxArray_Int[0];          // number of temperature points
   
   const double *Tks  = &Src_ExactCooling_General_AuxArray_Flt[8];         // Temperature points (K)
   const double *aks  = &Src_ExactCooling_General_AuxArray_Flt[8+N];       // Power law slopes
   const double *Lks  = &Src_ExactCooling_General_AuxArray_Flt[8+2*N-1];   // Normalizations (erg cm^3 s^-1)

   // ==========================
   // (2) OpenMP
   // ==========================
#  ifdef OPENMP
   const int NT = OMP_NTHREAD;
#  else
   const int NT = 1;
#  endif

   double  dt_Cool     = HUGE_NUMBER;
   double *OMP_dt_Cool = new double [NT];

   // ==========================
   // (3) Scan all cells
   // ==========================
#  pragma omp parallel
   {
#     ifdef OPENMP
      const int TID = omp_get_thread_num();
#     else
      const int TID = 0;
#     endif

      OMP_dt_Cool[TID] = HUGE_NUMBER;


#     pragma omp for schedule( static )
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
         for (int k=0; k<PS1; k++)
         {
         for (int j=0; j<PS1; j++)
         {
         for (int i=0; i<PS1; i++)
         {

            // =================================================
            // (3.1) Get fluid variables
            // =================================================

            real fluid[FLU_NIN_S];

            for (int v=0; v<FLU_NIN_S; v++)
               fluid[v] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i];


            // =================================================
            // (3.2) Calculate internal energy
            // =================================================

            double     Eint, Emag;
            const bool CheckMinEint_No = false;
            Emag  = 0.0;
            Eint  = Hydro_Con2Eint( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY],
                                    CheckMinEint_No, NULL_REAL, NULL, Emag, 
                                    EoS.GuessHTilde_FuncPtr, EoS.HTilde2Temp_FuncPtr, 
                                    EoS.AuxArrayDevPtr_Flt, EoS.AuxArrayDevPtr_Int,
                                    EoS.Table );
            
            //=================================================
            // (3.3) Calculate current temperature
            // =================================================
            
            double T_now = EoS.DensEint2Temp_FuncPtr( fluid[DENS], Eint, NULL, EoS.AuxArrayDevPtr_Flt, EoS.AuxArrayDevPtr_Int, EoS.Table   );

            // =================================================
            // (3.4) Skip cells out of temperature range
            // =================================================

            if ( T_now < Tks[0] || T_now > Tks[N-1] )
               continue;

            // =================================================
            // (3.5) Find the temperature interval
            // =================================================

            int k_interval = -1;

            for (int kk=N-2; kk>=0; kk--)
            {
               if ( T_now >= Tks[kk] )
               {
                  k_interval = kk + 1;
                  break;
               }
            }

            // This should never happen because Tks[0] <= T_now <= Tks[N-1].
            if ( k_interval < 1 || k_interval > N-1 )
               continue;


            // =================================================
            // (3.6) Calculate cooling efficiency Lambda(T)
            // =================================================

            const double Lambda = ExactCooling_GetLambda( T_now, k_interval, Tks, Lks, aks );

            //=================================================
            // (3.7) Calculate number densities
            // =================================================

            const double rho = fluid[DENS];
            const double n_H = rho / (mu_H * mp);
            const double n_e = rho / (mu_e * mp);
            const double n   = rho / (mu   * mp);

            // =================================================
            // (3.8) Calculate cooling time
            // =================================================

            const double tcool = ExactCooling_GetTcool( T_now, n, n_e, n_H, Lambda, kB, gamma );

            // =================================================
            // (3.9) Cooling timestep 
            // =================================================

            const double dt_cell = SrcTerms.ExactCooling_General_dt * tcool;

            // =================================================
            // (3.10) Store the minimum timestep
            // =================================================

            if ( Aux_IsFinite(dt_cell) && dt_cell > 0.0 )
               OMP_dt_Cool[TID] = fmin( OMP_dt_Cool[TID], dt_cell );

         }}} // i,j,k
      } // PID

   } // OpenMP parallel region


   // ==========================================================
   // (4) Minimum over OpenMP threads
   // ==========================================================

   for (int TID=0; TID<NT; TID++)
      dt_Cool = fmin( dt_Cool, OMP_dt_Cool[TID] );


   delete [] OMP_dt_Cool;


   // ==========================================================
   // (5) Minimum over MPI processes
   // ==========================================================

#  ifndef SERIAL
   MPI_Allreduce
   (
      MPI_IN_PLACE,
      &dt_Cool,
      1,
      MPI_DOUBLE,
      MPI_MIN,
      MPI_COMM_WORLD
   );
#  endif


   return dt_Cool;

} // FUNCTION : Mis_GetTimeStep_ExactCooling_General
# endif


//-------------------------------------------------------------------------------------------------------
// Function    : TEF_PiecewisePowerLaw
// Description : Temporal evolution function (TEF) for Piecewise Power Law Cooling Function
//
// Note        :
//
// Parameter   : Temp          : Temperature in Kelvin
//               k             : Temperature interval index
//               N             : Number of temperature intervals
//               Tks           : Array of temperature intervals
//               Lks           : Array of cooling rates
//               aks           : Array of power law slopes
//               Yks           : Array of Yk values
//
// Return      : TEF(Temp)
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE static
double TEF_PiecewisePowerLaw( const double Temp, const int k, const int N, const double Tks[], const double Lks[], const double aks[], double Yks[] )
{
   const double TN = Tks[N-1];
   const double LN = Lks[N-2] * pow( TN / Tks[N-2], aks[N-2] );

   const double Tk = Tks[k-1];
   const double Lk = Lks[k-1];
   const double ak = aks[k-1];

   // Yk
   Yks[N-1] = 0.0;   // YN = 0

   if ( ak == 1.0 )
   {
      for ( int i = N-1; i > 0; i-- )
      {
         const double Ti = Tks[i-1];
         const double Li = Lks[i-1];
         Yks[i-1] = Yks[i] - LN / Li * Ti / TN * log( Ti / Tks[i]);
      }
   }
   else
   {
      for ( int i = N-1; i > 0; i-- )
      {
         const double Ti = Tks[i-1];
         const double Li = Lks[i-1];
         const double ai = aks[i-1];
         Yks[i-1] = Yks[i] - 1.0 / (1.0 - ai) * LN / Li * Ti / TN * (1.0 - pow( Ti / Tks[i], (ai - 1.0)));
      }
   }

   // Y(T)
   double Y_T;
   if ( ak == 1.0 )
   {
      Y_T = Yks[k-1] + LN / Lk * Tk / TN * log( Tk / Temp );
   }
   else
   {
      Y_T = Yks[k-1] + 1.0 / (1.0 - ak) * LN / Lk * Tk / TN * (1.0 - pow( Tk / Temp, (ak - 1.0)));
   }

   return Y_T;

} // FUNCTION : TEF_PiecewisePowerLaw


//-------------------------------------------------------------------------------------------------------
// Function    : TEF_inverse_PiecewisePowerLaw
// Description : Inverse Temporal evolution function (TEF) for Piecewise Power Law Cooling Function
//
// Note        :
//
// Parameter   : Y             : Input parameter for the inverse TEF
//               N             : Number of temperature intervals
//               Yks           : Array of Yk values
//               Lks           : Array of cooling rates
//               aks           : Array of power law slopes
//               Tks           : Array of temperature intervals
//               
//
// Return      : inverseTEF(Y)
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE static
double TEF_inverse_PiecewisePowerLaw( const double Y, const int N, const double Yks[], const double Lks[], const double aks[], const double Tks[] )
{
   // cooling floor
   if ( Y >= Yks[0] )
       return Tks[0];
   
   int k = -1;
   for ( int i = 1; i < N; i++ )
   {
      if ( Yks[i] <= Y && Y < Yks[i-1] )
      {
         k = i;
         break;
      }
   }
   if ( k < 1 )
   {
   #  ifdef GAMER_DEBUG   
      printf( "ERROR : Invalid Y = %e in %s !!\n", Y, __FUNCTION__ );
   #  endif
      return Tks[0];
   }

   const double Tk = Tks[k-1];
   const double Lk = Lks[k-1];
   const double ak = aks[k-1];
   const double Yk = Yks[k-1];

   const double TN = Tks[N-1];
   const double LN = Lks[N-2] * pow( TN / Tks[N-2], aks[N-2] );

   // inverseY(Y)
   double inverseY;
   if ( ak == 1.0 )
   {
      inverseY = Tk * exp(-Lk / LN * TN / Tk * (Y - Yk));
   }
   else
   {
      inverseY = Tk * pow( (1.0 - (1.0 - ak) * Lk / LN * TN / Tk * (Y - Yk)), 1.0 / (1.0 - ak));
   }
   
   return inverseY;
}  // FUNCTION : TEF_inverse_PiecewisePowerLaw



// ==================================================
// III. [Optional] Add the work to be done every time
//      before calling the major source-term function
// ==================================================

//-------------------------------------------------------------------------------------------------------
// Function    :  Src_WorkBeforeMajorFunc_ExactCooling_General
// Description :  Specify work to be done every time before calling the major source-term function
//
// Note        :  1. Invoked by Src_WorkBeforeMajorFunc()
//                   --> By linking to "Src_WorkBeforeMajorFunc_User_Ptr" in Src_Init_ExactCooling_General()
//                2. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  lv               : Target refinement level
//                TimeNew          : Target physical time to reach
//                TimeOld          : Physical time before update
//                                   --> The major source-term function will update the system from TimeOld to TimeNew
//                dt               : Time interval to advance solution
//                                   --> Physical coordinates : TimeNew - TimeOld == dt
//                                       Comoving coordinates : TimeNew - TimeOld == delta(scale factor) != dt
//                AuxArray_Flt/Int : Auxiliary arrays
//                                   --> Can be used and/or modified here
//                                   --> Must call Src_SetConstMemory_ExactCooling_General() after modification
//
// Return      :  AuxArray_Flt/Int[]
//-------------------------------------------------------------------------------------------------------
#ifndef __CUDACC__
void Src_WorkBeforeMajorFunc_ExactCooling_General( const int lv, const double TimeNew, const double TimeOld, const double dt,
                                                   double AuxArray_Flt[], int AuxArray_Int[] )
{

// nothing to do here

} // FUNCTION : Src_WorkBeforeMajorFunc_ExactCooling_General
#endif // #ifndef __CUDACC__



// ================================
// IV. Set initialization functions
// ================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE SrcFunc_t SrcFunc_Ptr = Src_ExactCooling_General;

//-----------------------------------------------------------------------------------------
// Function    :  Src_SetCPU/GPUFunc_ExactCooling_General
// Description :  Return the function pointer of the CPU/GPU source-term function
//
// Note        :  1. Invoked by Src_Init_ExactCooling_General()
//                2. Call-by-reference
//
// Parameter   :  SrcFunc_CPU/GPUPtr : CPU/GPU function pointer to be set
//
// Return      :  SrcFunc_CPU/GPUPtr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void Src_SetGPUFunc_ExactCooling_General( SrcFunc_t &SrcFunc_GPUPtr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &SrcFunc_GPUPtr, SrcFunc_Ptr, sizeof(SrcFunc_t) )  );
} // FUNCTION : Src_SetGPUFunc_ExactCooling_General

#else

void Src_SetCPUFunc_ExactCooling_General( SrcFunc_t &SrcFunc_CPUPtr )
{
   SrcFunc_CPUPtr = SrcFunc_Ptr;
} // FUNCTION : Src_SetCPUFunc_ExactCooling_General

#endif // #ifdef __CUDACC__ ... else ...



#ifdef __CUDACC__
//-------------------------------------------------------------------------------------------------------
// Function    :  Src_SetConstMemory_ExactCooling_General
// Description :  Set the constant memory variables on GPU
//
// Note        :  1. Adopt the suggested approach for CUDA version >= 5.0
//                2. Invoked by Src_Init_ExactCooling_General() and, if necessary, Src_WorkBeforeMajorFunc_ExactCooling_General()
//                3. SRC_NAUX_EXACTCOOLING_GENERAL is defined in Macro.h
//
// Parameter   :  AuxArray_Flt/Int : Auxiliary arrays to be copied to the constant memory
//                DevPtr_Flt/Int   : Pointers to store the addresses of constant memory arrays
//
// Return      :  c_Src_ExactCooling_General_AuxArray_Flt[], c_Src_ExactCooling_General_AuxArray_Int[], DevPtr_Flt, DevPtr_Int
//---------------------------------------------------------------------------------------------------
void Src_SetConstMemory_ExactCooling_General( const double AuxArray_Flt[], const int AuxArray_Int[],
                                              double *&DevPtr_Flt, int *&DevPtr_Int )
{

// copy data to constant memory
   CUDA_CHECK_ERROR(  cudaMemcpyToSymbol( c_Src_ExactCooling_General_AuxArray_Flt, AuxArray_Flt, SRC_NAUX_EXACTCOOLING_GENERAL*sizeof(double) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyToSymbol( c_Src_ExactCooling_General_AuxArray_Int, AuxArray_Int, SRC_NAUX_EXACTCOOLING_GENERAL*sizeof(int   ) )  );

// obtain the constant-memory pointers
   CUDA_CHECK_ERROR(  cudaGetSymbolAddress( (void **)&DevPtr_Flt, c_Src_ExactCooling_General_AuxArray_Flt) );
   CUDA_CHECK_ERROR(  cudaGetSymbolAddress( (void **)&DevPtr_Int, c_Src_ExactCooling_General_AuxArray_Int) );

} // FUNCTION : Src_SetConstMemory_ExactCooling_General
#endif // #ifdef __CUDACC__



#ifndef __CUDACC__

//-----------------------------------------------------------------------------------------
// Function    :  Src_Init_ExactCooling_General
// Description :  Initialize a user-specified source term
//
// Note        :  1. Set auxiliary arrays by invoking Src_SetAuxArray_*()
//                   --> Copy to the GPU constant memory and store the associated addresses
//                2. Set the source-term function by invoking Src_SetCPU/GPUFunc_*()
//                3. Set the function pointers "Src_WorkBeforeMajorFunc_ExactCooling_General_Ptr" and "Src_End_ExactCooling_General_Ptr"
//                4. Invoked by Src_Init()
//                   --> Enable it by linking to the function pointer "Src_Init_ExactCooling_General_Ptr"
//                5. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  None
//
// Return      :  None
//-----------------------------------------------------------------------------------------
void Src_Init_ExactCooling_General()
{

// set the auxiliary arrays
   Src_SetAuxArray_ExactCooling_General( Src_ExactCooling_General_AuxArray_Flt, Src_ExactCooling_General_AuxArray_Int );

// copy the auxiliary arrays to the GPU constant memory and store the associated addresses
#  ifdef GPU
   Src_SetConstMemory_ExactCooling_General( Src_ExactCooling_General_AuxArray_Flt, Src_ExactCooling_General_AuxArray_Int,
                                            SrcTerms.ExactCooling_General_AuxArrayDevPtr_Flt, SrcTerms.ExactCooling_General_AuxArrayDevPtr_Int );
#  else
   SrcTerms.ExactCooling_General_AuxArrayDevPtr_Flt = Src_ExactCooling_General_AuxArray_Flt;
   SrcTerms.ExactCooling_General_AuxArrayDevPtr_Int = Src_ExactCooling_General_AuxArray_Int;
#  endif

// set the major source-term function
   Src_SetCPUFunc_ExactCooling_General( SrcTerms.ExactCooling_General_CPUPtr );

#  ifdef GPU
   Src_SetGPUFunc_ExactCooling_General( SrcTerms.ExactCooling_General_GPUPtr );
   SrcTerms.ExactCooling_General_FuncPtr = SrcTerms.ExactCooling_General_GPUPtr;
#  else
   SrcTerms.ExactCooling_General_FuncPtr = SrcTerms.ExactCooling_General_CPUPtr;
#  endif

   if ( OPT__INIT == INIT_BY_RESTART )
      for (int i=0; i<NLEVEL; i++)   SrcTerms.ExactCooling_General_TCoolInit[i] = true;
   else
      for (int i=0; i<NLEVEL; i++)   SrcTerms.ExactCooling_General_TCoolInit[i] = false;

// initialize the cooling function
   const int      N   = Src_ExactCooling_General_AuxArray_Int[0];   // number of temperature points
   
   const double X_H   = Src_ExactCooling_General_AuxArray_Flt[0];   // hydrogen mass fraction
   const double Z     = Src_ExactCooling_General_AuxArray_Flt[1];   // metal mass fraction
   const double mu_e  = Src_ExactCooling_General_AuxArray_Flt[2];   // mean electron molecular weight
   const double mu_H  = Src_ExactCooling_General_AuxArray_Flt[3];   // mean proton molecular weight
   const double mu    = Src_ExactCooling_General_AuxArray_Flt[4];   // mean total molecular weight

   const double mp    = Src_ExactCooling_General_AuxArray_Flt[5];   // proton mass (g)
   const double kB    = Src_ExactCooling_General_AuxArray_Flt[6];   // Boltzmann constant (erg/K)
   const double gamma = Src_ExactCooling_General_AuxArray_Flt[7];   // ratio of specific heats

   // please add T5, a4, L4,... if your cooling function have more than 3 power law
   const double T1    = Src_ExactCooling_General_AuxArray_Flt[8];   // temperature point 1 (K)
   const double T2    = Src_ExactCooling_General_AuxArray_Flt[9];   // temperature point 2 (K)
   const double T3    = Src_ExactCooling_General_AuxArray_Flt[10];  // temperature point 3 (K)
   const double T4    = Src_ExactCooling_General_AuxArray_Flt[11];  // temperature point 4 (K)
   
   const double a1    = Src_ExactCooling_General_AuxArray_Flt[12];  // slope 1
   const double a2    = Src_ExactCooling_General_AuxArray_Flt[13];  // slope 2
   const double a3    = Src_ExactCooling_General_AuxArray_Flt[14];  // slope 3

   const double L1    = Src_ExactCooling_General_AuxArray_Flt[15];  // normalization 1 (erg cm^3 s^-1)
   const double L2    = Src_ExactCooling_General_AuxArray_Flt[16];  // normalization 2 (erg cm^3 s^-1)
   const double L3    = Src_ExactCooling_General_AuxArray_Flt[17];  // normalization 3 (erg cm^3 s^-1)

} // FUNCTION : Src_Init_ExactCooling_General



//-----------------------------------------------------------------------------------------
// Function    :  Src_End_ExactCooling_General
// Description :  Free the resources used by a user-specified source term
//
// Note        :  1. Invoked by Src_End()
//                   --> Enable it by linking to the function pointer "Src_End_User_Ptr"
//                2. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  None
//
// Return      :  None
//-----------------------------------------------------------------------------------------
void Src_End_ExactCooling_General()
{


} // FUNCTION : Src_End_ExactCooling_General

#endif // #ifndef __CUDACC__

#endif // #ifdef EXACT_COOLING_GENERAL