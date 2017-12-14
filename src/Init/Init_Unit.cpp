#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Unit
// Description :  Initialize code units
//
// Note        :  1. Work with OPT__UNIT == true. If OPT__UNIT == false, all units are set to unity as if
//                   the simulation is dimensionless
//                2. One can specify three and only three units among UNIT_L/M/T/V/D in the input parameter file
//                   --> Other units are calculated from these three basic units
//                3. Some physical constants are also converted to code units by this function
//                   --> NEWTON_G, ELBDM_MASS, ELBDM_PLANCK_CONST, ...
//                4. For cosmological simulations, basic units are fixed to
//                       UNIT_L = Mpc/h
//                       UNIT_D = current matter density = 3*OMEGA_M0*H0^2/(8*pi*G)
//                       UNIT_T = 1/H0 = s*Mpc/(100*h*km)
//                   --> UNIT_V = UNIT_L/UNIT_T = Mpc/h * 100*h*km/(s*Mpc) = 100*km/s (**independent of h**)
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_Unit()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// set default units for cosmological simulations
#  ifdef COMOVING
   if ( OMEGA_M0 < 0.0  ||  OMEGA_M0 > 1.0 )
      Aux_Error( ERROR_INFO, "incorrect OMEGA_M0 (0.0 <= OMEGA_M0 <= 1.0) !!\n" );

   if ( HUBBLE0 <= 0.0 )
      Aux_Error( ERROR_INFO, "HUBBLE0 = %14.7e <= 0.0 !!\n", HUBBLE0 );

   OPT__UNIT = true;    // always turn on OPT__UNIT for COMOVING
   UNIT_L    = Const_Mpc/HUBBLE0;
   UNIT_T    = Const_s*Const_Mpc/(100.0*HUBBLE0*Const_km);
   UNIT_D    = 3.0*OMEGA_M0*SQR( 100.0*HUBBLE0*Const_km/(Const_s*Const_Mpc) )/( 8.0*M_PI*Const_NewtonG );
   UNIT_V    = -1.0;    // UNIT_V/M/E/P will be calculated automatically below
   UNIT_M    = -1.0;
   UNIT_E    = -1.0;
   UNIT_P    = -1.0;

   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "NOTE : code units in cosmological simulations are fixed to\n" );
      Aux_Message( stdout, "   UNIT_L = Mpc/h\n" );
      Aux_Message( stdout, "   UNIT_D = current matter density = 3*OMEGA_M0*H0^2/(8*pi*G)\n" );
      Aux_Message( stdout, "   UNIT_T = 1/H0 = s*Mpc/(100*h*km)\n" );
      Aux_Message( stdout, "   UNIT_V = UNIT_L/UNIT_T = 100 km/s\n" );
   }
#  endif // #ifdef COMOVING


// check user-provided units
   if ( OPT__UNIT )
   {
//    check if there are exactly three basic units
      int NBasicUnit = 0;
      if ( UNIT_L > 0.0 )    NBasicUnit ++;
      if ( UNIT_M > 0.0 )    NBasicUnit ++;
      if ( UNIT_T > 0.0 )    NBasicUnit ++;
      if ( UNIT_V > 0.0 )    NBasicUnit ++;
      if ( UNIT_D > 0.0 )    NBasicUnit ++;

      if ( NBasicUnit != 3 )  Aux_Error( ERROR_INFO, "Number of basic units set in Input__Parameter = %d != 3 !!\n", NBasicUnit );


//    set all code units
//    (1) length unit
      if ( UNIT_L <= 0.0 )
      {
         if ( UNIT_T > 0.0  &&  UNIT_V > 0.0 )
         {
            UNIT_L = UNIT_V*UNIT_T;
            if ( MPI_Rank == 0 )    Aux_Message( stdout, "NOTE : UNIT_L is set to %13.7e\n", UNIT_L );
         }

         else if ( UNIT_M > 0.0  &&  UNIT_D > 0.0 )
         {
            UNIT_L = pow( UNIT_M/UNIT_D, 1.0/3.0 );
            if ( MPI_Rank == 0 )    Aux_Message( stdout, "NOTE : UNIT_L is set to %13.7e\n", UNIT_L );
         }

         else
            Aux_Error( ERROR_INFO, "cannot determine the length unit !!\n" );
      }

//    (2) mass unit
      if ( UNIT_M <= 0.0 )
      {
         if ( UNIT_L > 0.0  &&  UNIT_D > 0.0 )
         {
            UNIT_M = UNIT_D*CUBE(UNIT_L);
            if ( MPI_Rank == 0 )    Aux_Message( stdout, "NOTE : UNIT_M is set to %13.7e\n", UNIT_M );
         }

         else
            Aux_Error( ERROR_INFO, "cannot determine the mass unit !!\n" );
      }

//    (3) time unit
      if ( UNIT_T <= 0.0 )
      {
         if ( UNIT_L > 0.0  &&  UNIT_V > 0.0 )
         {
            UNIT_T = UNIT_L/UNIT_V;
            if ( MPI_Rank == 0 )    Aux_Message( stdout, "NOTE : UNIT_T is set to %13.7e\n", UNIT_T );
         }

         else
            Aux_Error( ERROR_INFO, "cannot determine the time unit !!\n" );
      }

//    (4) velocity unit
      if ( UNIT_V <= 0.0 )
      {
         if ( UNIT_L > 0.0  &&  UNIT_T > 0.0 )
         {
            UNIT_V = UNIT_L/UNIT_T;
            if ( MPI_Rank == 0 )    Aux_Message( stdout, "NOTE : UNIT_V is set to %13.7e\n", UNIT_V );
         }

         else
            Aux_Error( ERROR_INFO, "cannot determine the velocity unit !!\n" );
      }

//    (5) mass density unit
      if ( UNIT_D <= 0.0 )
      {
         if ( UNIT_L > 0.0  &&  UNIT_M > 0.0 )
         {
            UNIT_D = UNIT_M/CUBE(UNIT_L);
            if ( MPI_Rank == 0 )    Aux_Message( stdout, "NOTE : UNIT_D is set to %13.7e\n", UNIT_D );
         }

         else
            Aux_Error( ERROR_INFO, "cannot determine the density unit !!\n" );
      }

//    (6) energy unit (which cannot be set by users)
      if ( true )
      {
         if ( UNIT_M > 0.0  &&  UNIT_V > 0.0 )
         {
            UNIT_E = UNIT_M*SQR(UNIT_V);
            if ( MPI_Rank == 0 )    Aux_Message( stdout, "NOTE : UNIT_E is set to %13.7e\n", UNIT_E );
         }

         else
            Aux_Error( ERROR_INFO, "cannot determine the energy unit !!\n" );
      }

//    (7) energy density unit (which cannot be set by users)
      if ( true )
      {
         if ( UNIT_D > 0.0  &&  UNIT_V > 0.0 )
         {
            UNIT_P = UNIT_D*SQR(UNIT_V);
            if ( MPI_Rank == 0 )    Aux_Message( stdout, "NOTE : UNIT_P is set to %13.7e\n", UNIT_P );
         }

         else
            Aux_Error( ERROR_INFO, "cannot determine the energy density unit !!\n" );
      }


//    convert physical constants to code units
//    (1) gravitational constant
#     ifdef GRAVITY
      NEWTON_G = Const_NewtonG / ( 1.0/UNIT_D/SQR(UNIT_T) );
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "NOTE : NEWTON_G is set to %13.7e internally\n", NEWTON_G );
#     endif

#     if ( MODEL == ELBDM )
//    (2) particle mass in ELBDM
//    note that the input value of ELBDM_MASS is always in eV/c^2 when OPT__UNIT is on (i.e., independent of UNIT_*)
      ELBDM_MASS = ELBDM_MASS*Const_eV/SQR(Const_c) / UNIT_M;
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "NOTE : ELBDM_MASS is set to %13.7e internally\n", ELBDM_MASS );

//    (3) reduced Planck constant
      ELBDM_PLANCK_CONST = Const_Planck / ( SQR(UNIT_L)*UNIT_M/UNIT_T );
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "NOTE : ELBDM_PLANCK_CONST is set to %13.7e internally\n", ELBDM_PLANCK_CONST );
#     endif

   } // if ( OPT__UNIT )



// set all code units to unity if OPT__UNIT == false (just in case any of UNIT_* is misused)
   else
   {
      UNIT_L = UNIT_M = UNIT_T = UNIT_V = UNIT_D = UNIT_E = UNIT_P = 1.0;

//    check if physical constants are properly set by users
#     ifdef GRAVITY
      if ( NEWTON_G <= 0.0 )              Aux_Error( ERROR_INFO, "NEWTON_G = %14.7e <= 0.0 !!\n", NEWTON_G );
#     endif

#     if ( MODEL == ELBDM )
      if ( ELBDM_PLANCK_CONST <= 0.0 )    Aux_Error( ERROR_INFO, "ELBDM_PLANCK_CONST = %14.7e <= 0.0 !!\n", ELBDM_PLANCK_CONST );
#     endif
   } // if ( OPT__UNIT ) ... else ...


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_Unit


