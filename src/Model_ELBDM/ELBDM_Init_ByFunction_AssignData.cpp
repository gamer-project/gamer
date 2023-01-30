#include "GAMER.h"

#if ( MODEL == ELBDM )

// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void Init_Function_User_Template( real fluid[], const double x, const double y, const double z, const double Time,
                                         const int lv, double AuxArray[] );

// this function pointer must be set by a test problem initializer
void (*Init_Function_User_Ptr)( real fluid[], const double x, const double y, const double z, const double Time,
                                const int lv, double AuxArray[] ) = NULL;

extern int (*Flu_ResetByUser_Func_Ptr)( real fluid[], const double x, const double y, const double z, const double Time,
                                        const int lv, double AuxArray[] );




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Function_User_Template
// Description :  Function template to initialize the ELBDM field
//
// Note        :  1. Invoked by ELBDM_Init_ByFunction_AssignData() using the function pointer
//                   "Init_Function_User_Ptr", which must be set by a test problem initializer
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   (unless OPT__INIT_GRID_WITH_OMP is disabled)
//                   --> Please ensure that everything here is thread-safe
//
// Parameter   :  fluid    : ELBDM field array to be initialized
//                x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void Init_Function_User_Template( real fluid[], const double x, const double y, const double z, const double Time,
                                  const int lv, double AuxArray[] )
{

   const double C1[3]   = { 0.5*amr->BoxSize[0]+100.0,
                            0.5*amr->BoxSize[1]+200.0,
                            0.5*amr->BoxSize[2]+300.0 };
   const double C2[3]   = { 20.0, 40.0, 10.0 };
   const double Height1 =   5.0;
   const double Height2 =   8.0;
   const double Width1  =  64.0;
   const double Width2  = 512.0;

   fluid[REAL] = 1.0 + Height1*exp(  -( SQR(x-C1[0]) + SQR(y-C1[1]) + SQR(z-C1[2]) ) /SQR(Width1)  );
   fluid[IMAG] = 1.0 + Height2*exp(  -( SQR(x-C2[0]) + SQR(y-C2[1]) + SQR(z-C2[2]) ) /SQR(Width2)  );
   fluid[DENS] = fluid[REAL]*fluid[REAL] + fluid[IMAG]*fluid[IMAG];

// ELBDM does not support passive scalars yet ...

} // FUNCTION : Init_Function_User_Template



//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_Init_ByFunction_AssignData
// Description :  Construct the initial condition in ELBDM
//
// Note        :  1. Work for the option "OPT__INIT == INIT_BY_FUNCTION"
//                2. Function pointers "Init_Function_User_Ptr/Flu_ResetByUser_Func_Ptr" must be set by a
//                   test problem initializer
//                3. One can disable OpenMP in this routine by setting OPT__INIT_GRID_WITH_OMP = 0
//                   --> Useful if "Init_Function_User_Ptr/Flu_ResetByUser_Func_Ptr" do not support OpenMP
//                       (e.g., they may not be thread-safe or may involve a random number generator for which
//                       all threads would share the same random seed when adopting OpenMP)
//
// Parameter   :  lv : Target refinement level
//-------------------------------------------------------------------------------------------------------
void ELBDM_Init_ByFunction_AssignData( const int lv )
{

// check
   if ( Init_Function_User_Ptr == NULL )  Aux_Error( ERROR_INFO, "Init_Function_User_Ptr == NULL !!\n" );

   if ( OPT__RESET_FLUID_INIT  &&  Flu_ResetByUser_Func_Ptr == NULL )
      Aux_Error( ERROR_INFO, "Flu_ResetByUser_Func_Ptr == NULL for OPT__RESET_FLUID_INIT !!\n" );


// set the number of OpenMP threads
#  ifdef OPENMP
   const int OMP_NThread = ( OPT__INIT_GRID_WITH_OMP ) ? OMP_NTHREAD : 1;
#  endif


   const int    NSub   = ( INIT_SUBSAMPLING_NCELL <= 0 ) ? 1 : INIT_SUBSAMPLING_NCELL;
   const double dh     = amr->dh[lv];
   const double dh_sub = dh / NSub;
   const double _NSub3 = 1.0/(NSub*NSub*NSub);

   real   fluid[NCOMP_TOTAL], fluid_sub[NCOMP_TOTAL];
   double x, y, z, x0, y0, z0;


#  pragma omp parallel for private( fluid, fluid_sub, x, y, z, x0, y0, z0 ) schedule( runtime ) num_threads( OMP_NThread )
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   for (int k=0; k<PS1; k++)  {  z0 = amr->patch[0][lv][PID]->EdgeL[2] + k*dh + 0.5*dh_sub;
   for (int j=0; j<PS1; j++)  {  y0 = amr->patch[0][lv][PID]->EdgeL[1] + j*dh + 0.5*dh_sub;
   for (int i=0; i<PS1; i++)  {  x0 = amr->patch[0][lv][PID]->EdgeL[0] + i*dh + 0.5*dh_sub;

      for (int v=0; v<NCOMP_TOTAL; v++)   fluid[v] = 0.0;

      for (int kk=0; kk<NSub; kk++)    {  z = z0 + kk*dh_sub;
      for (int jj=0; jj<NSub; jj++)    {  y = y0 + jj*dh_sub;
      for (int ii=0; ii<NSub; ii++)    {  x = x0 + ii*dh_sub;

         Init_Function_User_Ptr( fluid_sub, x, y, z, Time[lv], lv, NULL );

//       modify the initial condition if required
         if ( OPT__RESET_FLUID_INIT )
            Flu_ResetByUser_Func_Ptr( fluid_sub, x, y, z, Time[lv], lv, NULL );

         for (int v=0; v<NCOMP_TOTAL; v++)   fluid[v] += fluid_sub[v];

      }}}

//    ensure density = real_part^2 + imaginary_part^2
      fluid[REAL] *= _NSub3;
      fluid[IMAG] *= _NSub3;
      fluid[DENS]  = fluid[REAL]*fluid[REAL] + fluid[IMAG]*fluid[IMAG];

//    check minimum density (but keep phase fixed)
      if ( fluid[DENS] < (real)MIN_DENS )
      {
         const real Rescale = SQRT( (real)MIN_DENS/fluid[DENS] );

         fluid[REAL] *= Rescale;
         fluid[IMAG] *= Rescale;
         fluid[DENS]  = (real)MIN_DENS;
      }

//    floor and normalize passive scalars (actually passive scalars are NOT supported by ELBDM yet)
      /*
#     if ( NCOMP_PASSIVE > 0 )
      for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  fluid[v] = FMAX( fluid[v], TINY_NUMBER );

      if ( OPT__NORMALIZE_PASSIVE )
         Hydro_NormalizePassive( fluid[DENS], fluid+NCOMP_FLUID, PassiveNorm_NVar, PassiveNorm_VarIdx );
#     endif
      */

      for (int v=0; v<NCOMP_TOTAL; v++)   amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i] = fluid[v];

   }}}

} // FUNCTION : ELBDM_Init_ByFunction_AssignData



#endif // #if ( MODEL == ELBDM )
