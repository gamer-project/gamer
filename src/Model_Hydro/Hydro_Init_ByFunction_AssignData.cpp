#include "GAMER.h"

#if ( MODEL == HYDRO )

// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void Init_Function_User( real fluid[], const double x, const double y, const double z, const double Time,
                                const int lv, double AuxArray[] );

// this function pointer may be overwritten by various test problem initializers
void (*Init_Function_User_Ptr)( real fluid[], const double x, const double y, const double z, const double Time,
                                const int lv, double AuxArray[] ) = Init_Function_User;

extern bool (*Flu_ResetByUser_Func_Ptr)( real fluid[], const double x, const double y, const double z, const double Time,
                                         const int lv, double AuxArray[] );

#ifdef MHD
// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void Init_Function_BField_User( real magnetic[], const double x, const double y, const double z, const double Time,
                                       const int lv, double AuxArray[] );

// this function pointer may be overwritten by various test problem initializers
void (*Init_Function_BField_User_Ptr)( real magnetic[], const double x, const double y, const double z, const double Time,
                                       const int lv, double AuxArray[] ) = Init_Function_BField_User;
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Function_User
// Description :  Function to initialize the fluid field
//
// Note        :  1. Invoked by Hydro_Init_ByFunction_AssignData() using the function pointer "Init_Function_User_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   (unless OPT__INIT_GRID_WITH_OMP is disabled)
//                   --> Please ensure that everything here is thread-safe
//                3. Even when DUAL_ENERGY is adopted, one does NOT need to set the dual-energy variable here
//                   --> It will be set automatically in Hydro_Init_ByFunction_AssignData()
//                4. For MHD, do NOT add magnetic energy (i.e., 0.5*B^2) to fluid[ENGY] here
//                   --> It will be added automatically later
//
// Parameter   :  fluid    : Fluid field to be initialized
//                x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void Init_Function_User( real fluid[], const double x, const double y, const double z, const double Time,
                         const int lv, double AuxArray[] )
{

   const real P0   = 1.0;
   const real Rho0 = 1.0;
   const real Vx   = 1.25e-1;
   const real Vy   = 2.30e-1;
   const real Vz   = 3.70e-1;

   fluid[DENS] = Rho0 + 0.2*exp( -( SQR(1.1*x-0.5*amr->BoxSize[0])+SQR(2.2*y-0.5*amr->BoxSize[1])+SQR(3.3*z-0.5*amr->BoxSize[2]) ) / SQR(1.8*amr->BoxSize[2]) );
   fluid[MOMX] = fluid[DENS]*Vx + 0.1;
   fluid[MOMY] = fluid[DENS]*Vy + 0.2;
   fluid[MOMZ] = fluid[DENS]*Vz + 0.3;
   fluid[ENGY] = (P0)/(GAMMA-1.0)*(2.0+sin(2.0*M_PI*(4.5*x+5.5*y*6.5*z)/amr->BoxSize[2]))
                 + 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];

// set passive scalars

} // FUNCTION : Init_Function_User



#ifdef MHD
//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Function_BField_User
// Description :  Function to initialize the magnetic field
//
// Note        :  1. Invoked by Hydro_Init_ByFunction_AssignData() using the function pointer "Init_Function_BField_User_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   (uless OPT__INIT_GRID_WITH_OMP is disabled)
//                   --> Please ensure that everything here is thread-safe
//
// Parameter   :  magnetic : Array to store the output magnetic field
//                x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  magnetic
//-------------------------------------------------------------------------------------------------------
void Init_Function_BField_User( real magnetic[], const double x, const double y, const double z, const double Time,
                                const int lv, double AuxArray[] )
{

   magnetic[MAGX] = 1.0;
   magnetic[MAGY] = 2.0;
   magnetic[MAGZ] = 3.0;

} // FUNCTION : Init_Function_BField_User
#endif // #ifdef MHD



//-------------------------------------------------------------------------------------------------------
// Function    :  HydroInit_ByFunction_AssignData
// Description :  Construct the initial condition in HYDRO
//
// Note        :  1. Work for the option "OPT__INIT == INIT_BY_FUNCTION"
//                2. The function pointers "Init_Function_User_Ptr/Flu_ResetByUser_Func_Ptr/Init_Function_BField_User_Ptr"
//                   point to "Init_Function_User()/Flu_ResetByUser_Func()/Init_Function_BField_User()" by default
//                   but may be overwritten by various test problem initializers
//                3. One can disable OpenMP in this routine by setting OPT__INIT_GRID_WITH_OMP = 0
//                   --> Useful if "Init_Function_User_Ptr/Flu_ResetByUser_Func_Ptr/Init_Function_BField_User_Ptr"
//                       do not support OpenMP
//                       (e.g., they may not be thread-safe or may involve a random number generator for which
//                       all threads would share the same random seed when adopting OpenMP)
//
// Parameter   :  lv : Target refinement level
//-------------------------------------------------------------------------------------------------------
void Hydro_Init_ByFunction_AssignData( const int lv )
{

// check
   if ( Init_Function_User_Ptr == NULL )  Aux_Error( ERROR_INFO, "Init_Function_User_Ptr == NULL !!\n" );

#  ifdef MHD
   if ( Init_Function_BField_User_Ptr == NULL )    Aux_Error( ERROR_INFO, "Init_Function_BField_User_Ptr == NULL !!\n" );
#  endif


// set the number of OpenMP threads
#  ifdef OPENMP
   const int OMP_NThread = ( OPT__INIT_GRID_WITH_OMP ) ? OMP_NTHREAD : 1;
#  endif

   const int    NSub     = ( INIT_SUBSAMPLING_NCELL <= 0 ) ? 1 : INIT_SUBSAMPLING_NCELL;
   const double dh       = amr->dh[lv];
   const double dh_sub   = dh / NSub;
   const double _NSub3   = 1.0/CUBE(NSub);
   const real   Gamma_m1 = GAMMA - (real)1.0;
   const real  _Gamma_m1 = (real)1.0 / Gamma_m1;
#  ifdef MHD
   const double _NSub2   = 1.0/SQR(NSub);
#  endif


#  pragma omp parallel for schedule( runtime ) num_threads( OMP_NThread )
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
//    1. set the magnetic field
#     ifdef MHD
      real magnetic_1v, magnetic_sub[NCOMP_MAG];

//    loop over B_X/Y/Z to set one component at a time
//    --> because different components are defined at different cell faces
      for (int v=0; v<NCOMP_MAG; v++)
      {
         int    ijk_end[3], sub_end[3], idx=0;
         double dxyz0[3];

         for (int d=0; d<3; d++)
         {
            ijk_end[d] = ( d == v ) ? PS1+1 : PS1;
            sub_end[d] = ( d == v ) ? 1     : NSub;
            dxyz0  [d] = ( d == v ) ? 0.0   : 0.5*dh_sub;
         }

         for (int k=0; k<ijk_end[2]; k++)    {  const double z0 = amr->patch[0][lv][PID]->EdgeL[2] + k*dh + dxyz0[2];
         for (int j=0; j<ijk_end[1]; j++)    {  const double y0 = amr->patch[0][lv][PID]->EdgeL[1] + j*dh + dxyz0[1];
         for (int i=0; i<ijk_end[0]; i++)    {  const double x0 = amr->patch[0][lv][PID]->EdgeL[0] + i*dh + dxyz0[0];

            magnetic_1v = (real)0.0;

            for (int kk=0; kk<sub_end[2]; kk++)    {  const double z = z0 + kk*dh_sub;
            for (int jj=0; jj<sub_end[1]; jj++)    {  const double y = y0 + jj*dh_sub;
            for (int ii=0; ii<sub_end[0]; ii++)    {  const double x = x0 + ii*dh_sub;

               Init_Function_BField_User_Ptr( magnetic_sub, x, y, z, Time[lv], lv, NULL );

               magnetic_1v += magnetic_sub[v];

            }}}

            amr->patch[ amr->MagSg[lv] ][lv][PID]->magnetic[v][ idx ++ ] = magnetic_1v*_NSub2;
         }}} // i,j,k
      } // for (int v=0; v<NCOMP_MAG; v++)
#     endif // #ifdef MHD



//    2. set the fluid field
      real fluid[NCOMP_TOTAL], fluid_sub[NCOMP_TOTAL];

      for (int k=0; k<PS1; k++)  {  const double z0 = amr->patch[0][lv][PID]->EdgeL[2] + k*dh + 0.5*dh_sub;
      for (int j=0; j<PS1; j++)  {  const double y0 = amr->patch[0][lv][PID]->EdgeL[1] + j*dh + 0.5*dh_sub;
      for (int i=0; i<PS1; i++)  {  const double x0 = amr->patch[0][lv][PID]->EdgeL[0] + i*dh + 0.5*dh_sub;

         for (int v=0; v<NCOMP_TOTAL; v++)   fluid[v] = (real)0.0;

         for (int kk=0; kk<NSub; kk++)    {  const double z = z0 + kk*dh_sub;
         for (int jj=0; jj<NSub; jj++)    {  const double y = y0 + jj*dh_sub;
         for (int ii=0; ii<NSub; ii++)    {  const double x = x0 + ii*dh_sub;

            Init_Function_User_Ptr( fluid_sub, x, y, z, Time[lv], lv, NULL );

//          modify the initial condition if required
            if ( OPT__RESET_FLUID  &&  Flu_ResetByUser_Func_Ptr != NULL)
               Flu_ResetByUser_Func_Ptr( fluid_sub, x, y, z, Time[lv], lv, NULL );

            for (int v=0; v<NCOMP_TOTAL; v++)   fluid[v] += fluid_sub[v];

         }}}

         for (int v=0; v<NCOMP_TOTAL; v++)   fluid[v] *= _NSub3;


//       add the magnetic energy
#        ifdef MHD
         const real EngyB = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, amr->MagSg[lv] );
         fluid[ENGY] += EngyB;
#        else
         const real EngyB = NULL_REAL;
#        endif


//       check minimum density and pressure
         fluid[DENS] = FMAX( fluid[DENS], (real)MIN_DENS );
         fluid[ENGY] = Hydro_CheckMinPresInEngy( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY],
                                                 Gamma_m1, _Gamma_m1, MIN_PRES, EngyB );

//       calculate the dual-energy variable (entropy or internal energy)
#        if   ( DUAL_ENERGY == DE_ENPY )
         fluid[ENPY] = Hydro_Fluid2Entropy( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY], Gamma_m1, EngyB );
#        elif ( DUAL_ENERGY == DE_EINT )
#        error : DE_EINT is NOT supported yet !!
#        endif

//       floor and normalize passive scalars
#        if ( NCOMP_PASSIVE > 0 )
         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  fluid[v] = FMAX( fluid[v], TINY_NUMBER );

         if ( OPT__NORMALIZE_PASSIVE )
            Hydro_NormalizePassive( fluid[DENS], fluid+NCOMP_FLUID, PassiveNorm_NVar, PassiveNorm_VarIdx );
#        endif

         for (int v=0; v<NCOMP_TOTAL; v++)   amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i] = fluid[v];

      }}} // i,j,k
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

} // FUNCTION : Hydro_Init_ByFunction_AssignData



#endif // #if ( MODEL == HYDRO )
