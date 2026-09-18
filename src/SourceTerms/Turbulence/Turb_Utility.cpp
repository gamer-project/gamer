#include "GAMER.h"
#include <algorithm>

#ifdef TURBULENCE

extern void Src_SetAuxArray_Turbulence( double [], int [] );
extern void Src_SetConstMemory_Turbulence( const double AuxArray_Flt[], const int AuxArray_Int[],
                                           double *&DevPtr_Flt, int *&DevPtr_Int );
extern void Src_PassData2GPU_Turbulence( int IdxTable );

/********************************************************************************************************
Turbulence structure:

1. Stores turbulence Fourier modes, amplitude, random phases, and sin, cos Fourier basis

2. Use Helmholtz decomposition to separate compressive and solenoial components

3. Update random phases by Ornstein-Uhlenbeck process:
   x(t+dt) = f * x(t) + sigma * sqrt (1 - f^2) * z_n
   where f = exp( -dt/tau ), tau is correlation time, z_n is Gaussian random number
   correlation <x(t+dt), x(t)> = sigma^2 * f

4. References: Federrath et al. (2010), A&A 512, A81 (https://doi.org/10.1051/0004-6361/200912437)
               TurbGen (https://github.com/chfeder/turbulence_generator)

*********************************************************************************************************/



//-------------------------------------------------------------------------------------------------------
// Function    :  Turb_Init_Modes
// Description :  Initialize turbulence modes, amplitudes
//
// Note        :  1. Invoked by Src_Init_Turbulence()
//
// Parameter   :  None
//
// Return      :  Turb
//-------------------------------------------------------------------------------------------------------
void Turb_Init_Modes()
{
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

   if ( Turb == NULL )     Aux_Error( ERROR_INFO, "Turb == NULL !!\n" );

// assign values to structure members
   Turb->Tdecay   = BOX_SIZE / SRC_TURB_KDRIV / SRC_TURB_VEL;
   Turb->dt       = Turb->Tdecay / SRC_TURB_UPDATE_STEP;
   Turb->SetRNGState( SRC_TURB_RSEED_INIT );

   const double ZetaNorm = sqrt(3.0) / sqrt( 1.0 - 2.0*SRC_TURB_ZETA + 3.0*SQR( SRC_TURB_ZETA ) );
   const double EnergyInputRate = CUBE( SRC_TURB_AMPL_FACTOR*0.15*SRC_TURB_VEL ) / BOX_SIZE;

// OUvar ~ a_rms
   Turb->OUvar = sqrt( EnergyInputRate/Turb->Tdecay );

// initialize k modes
   double kmin   = (SRC_TURB_KMIN - __DBL_EPSILON__) * 2*M_PI / BOX_SIZE;
   double kmax   = (SRC_TURB_KMAX + __DBL_EPSILON__) * 2*M_PI / BOX_SIZE;

   if ( kmax < kmin )
      Aux_Error( ERROR_INFO, "Turbulence: kmax ( %13.7e ) < kmin ( %13.7e )!!\n", kmax, kmin );

   double kmid   = 0.5*(kmin + kmax);
   int    Nmax   = 2*(int)SRC_TURB_KMAX + 1;
   int    nmodes = 0;

   for (int k = 0; k < Nmax; k++) {
   for (int j = 0; j < Nmax; j++) {
   for (int i = 0; i < Nmax; i++) {
      double kx = 2*M_PI/amr->BoxSize[0] * (i - (Nmax - 1)/2.0);
      double ky = 2*M_PI/amr->BoxSize[1] * (j - (Nmax - 1)/2.0);
      double kz = 2*M_PI/amr->BoxSize[2] * (k - (Nmax - 1)/2.0);

      double kmag = sqrt( SQR(kx) + SQR(ky) + SQR(kz) );
      if ( kmag >= kmin && kmag <= kmax ) nmodes++;

   }}}

   if ( nmodes == 0 )  Aux_Error( ERROR_INFO, "number of turbulence modes = 0 !!\n" );

   if ( nmodes > SRC_TURB_MAX_NMODE )  Aux_Error( ERROR_INFO, "number of turbulence modes ( %d ) exceeds maximum mode (%d) !!\n"
                                                              "try lowering SRC_TURB_KMAX !!\n" , nmodes, SRC_TURB_MAX_NMODE );

   if ( MPI_Rank == 0 ) Aux_Message( stdout, "   initialize %d turbulence modes\n", nmodes );

   Turb->NMode = nmodes;

   if ( Turb->Amplitude == NULL)  Turb->Amplitude = new double [nmodes];

   for (int i = 0; i < 3; i++)
   {
      if ( Turb->Kmode[i] == NULL ) Turb->Kmode[i] = new double [nmodes];
   }

// get amplitude of each mode
   int n = 0;
   for (int k = 0; k < Nmax; k++)
   for (int j = 0; j < Nmax; j++)
   for (int i = 0; i < Nmax; i++)
   {
      double amp = 0;
      double kx  = 2*M_PI/amr->BoxSize[0] * (i - (Nmax - 1)/2.0);
      double ky  = 2*M_PI/amr->BoxSize[1] * (j - (Nmax - 1)/2.0);
      double kz  = 2*M_PI/amr->BoxSize[2] * (k - (Nmax - 1)/2.0);

      double kmag = sqrt( SQR(kx) + SQR(ky) + SQR(kz) );
      if ( kmag >= kmin && kmag <= kmax ) {
//       constant
         if ( SRC_TURB_SPEC_FORM == 0)
            amp = 1.0*kmin/kmag;
//       parabolic
         else if ( SRC_TURB_SPEC_FORM == 1 )
            amp = sqrt( fabs(-4 * SQR( (kmag - kmid)/(kmax - kmin) ) + 1) )*kmid/kmag;
//       power law
         else if ( SRC_TURB_SPEC_FORM == 2 )
            amp = sqrt( pow(kmag/kmin, SRC_TURB_POW) )*kmin/kmag;
         else
            Aux_Error( ERROR_INFO, "Unknown TURB_SPEC_FORM = %d!!\n", SRC_TURB_SPEC_FORM );

         Turb->Kmode[0][n] = kx;
         Turb->Kmode[1][n] = ky;
         Turb->Kmode[2][n] = kz;

         Turb->Amplitude[n] = amp*2*ZetaNorm;
         n++;
      } // if ( kmag >= kmin && kmag <= kmax )
   } // for i, j, k

// print turbulence information
   if ( MPI_Rank == 0 )
   {
       Aux_Message( stdout, "Turbulence parameters:\n" );
       Aux_Message( stdout, "   velocity dispersion    = %13.7e\n", SRC_TURB_VEL         );
       Aux_Message( stdout, "   amplitude factor       = %13.7e\n", SRC_TURB_AMPL_FACTOR );
       Aux_Message( stdout, "   energy injection rate  = %13.7e\n", EnergyInputRate      );
       Aux_Message( stdout, "   kmin                   = %13.7e\n", kmin                 );
       Aux_Message( stdout, "   kmax                   = %13.7e\n", kmax                 );
       Aux_Message( stdout, "   correlation time       = %13.7e\n", Turb->Tdecay         );
       Aux_Message( stdout, "   update pattern dt      = %13.7e\n", Turb->dt             );
       Aux_Message( stdout, "   OU variance            = %13.7e\n", Turb->OUvar          );
       Aux_Message( stdout, "   solenoidal weight norm = %13.7e\n", ZetaNorm             );
       Aux_Message( stdout, "\n");

      if ( OPT__VERBOSE )
      for (int n = 0; n < Turb->NMode; n++)
      {
         Aux_Message( stdout, "    mode = %3d, amplitude = %13.7e\n", n, Turb->Amplitude[n] );
      }
   }

} // FUNCTION : Turb_Init_Modes



//-------------------------------------------------------------------------------------------------------
// Function    :  Turb_Init_Field
// Description :  Initialize turbulence OU phases, and fill in AccTable
//
// Note        :  1. Invoked by Init_GAMER() after Time[lv] is initialized
//                2. When restart, load OU phases, times, and rng state.
//
// Parameter   :  None
//
// Return      :  Turb
//-------------------------------------------------------------------------------------------------------
void Turb_Init_Field()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

   if ( Turb == NULL )     Aux_Error( ERROR_INFO, "Turb == NULL !!\n" );

// initialize OU noise
   for (int t = 0; t < 2; t++)
      if ( Turb->OUphase[t] == NULL ) Turb->OUphase[t] = new double [6*Turb->NMode];

// when restart, load turbulence field
   if ( OPT__INIT == INIT_BY_RESTART && !OPT__RESTART_RESET && !SRC_TURB_RESET )
   {
//    load with rank 0
      if ( MPI_Rank == 0 )
      {
#ifndef  SUPPORT_HDF5
         Aux_Error( ERROR_INFO, "restart turbulence field must enable SUPPORT_HDF5!!\n" );
#else
         const char FileName[] = "RESTART";

         if ( !Aux_CheckFileExist(FileName) )
            Aux_Error( ERROR_INFO, "restart HDF5 file \"%s\" does not exist !!\n", FileName );

         if ( !H5Fis_hdf5(FileName) )
            Aux_Error( ERROR_INFO, "restart HDF5 file \"%s\" is not in the HDF5 format !!\n", FileName );

         hid_t  H5_FileID, H5_GroupID_Turb, H5_SetID_Turb;
         herr_t H5_Status;

         H5_FileID = H5Fopen( FileName, H5F_ACC_RDONLY, H5P_DEFAULT );
         if ( H5_FileID < 0 )         Aux_Error( ERROR_INFO, "failed to open the restart HDF5 file \"%s\" !!\n", FileName );

         H5_GroupID_Turb = H5Gopen( H5_FileID, "Turbulence", H5P_DEFAULT );
         if ( H5_GroupID_Turb < 0 )   Aux_Error( ERROR_INFO, "failed to open the group \"%s\" !!\n"
                                                             "set SRC_TURB_RESET or OPT__RESTART_RESET to reset turbulence when restarting !!\n", "Turbulence" );

         int RS_NMode;
         H5_SetID_Turb = H5Dopen ( H5_GroupID_Turb, "NMode", H5P_DEFAULT);
         H5_Status     = H5Dread ( H5_SetID_Turb, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &RS_NMode );
         H5_Status     = H5Dclose( H5_SetID_Turb );

//       check
         if ( Turb->NMode != RS_NMode )
             Aux_Error( ERROR_INFO, "number of modes (%d) != number of modes from RESTART (%d) !!\n"
                                    "enable SRC_TURB_RESET to change the spectrum !\n" , Turb->NMode, RS_NMode );

         H5_SetID_Turb = H5Dopen ( H5_GroupID_Turb, "TimeLast", H5P_DEFAULT);
         H5_Status     = H5Dread ( H5_SetID_Turb, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Turb->TimeLast );
         H5_Status     = H5Dclose( H5_SetID_Turb );

         H5_SetID_Turb = H5Dopen ( H5_GroupID_Turb, "TimeNext", H5P_DEFAULT);
         H5_Status     = H5Dread ( H5_SetID_Turb, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Turb->TimeNext );
         H5_Status     = H5Dclose( H5_SetID_Turb );

         H5_SetID_Turb = H5Dopen ( H5_GroupID_Turb, "RNGState", H5P_DEFAULT);
         H5_Status     = H5Dread ( H5_SetID_Turb, H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Turb->RNGState );
         H5_Status     = H5Dclose( H5_SetID_Turb );

         H5_SetID_Turb = H5Dopen ( H5_GroupID_Turb, "OUArrLast", H5P_DEFAULT);
         H5_Status     = H5Dread ( H5_SetID_Turb, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Turb->OUphase[Turb->IdxLast] );
         if ( H5_Status < 0 ) Aux_Error( ERROR_INFO, "Failed to load OUphase from RESTART !!\n" );
         H5_Status     = H5Dclose( H5_SetID_Turb );

         H5_SetID_Turb = H5Dopen ( H5_GroupID_Turb, "OUArrNext", H5P_DEFAULT);
         H5_Status     = H5Dread ( H5_SetID_Turb, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Turb->OUphase[Turb->IdxNext] );
         if ( H5_Status < 0 ) Aux_Error( ERROR_INFO, "Failed to load OUphase from RESTART !!\n" );
         H5_Status     = H5Dclose( H5_SetID_Turb );

//       close file
         H5_Status = H5Gclose( H5_GroupID_Turb  );
         H5_Status = H5Fclose( H5_FileID        );
#        endif // ifdef SUPPORT_HDF5
      } // if ( MPI_Rank == 0 )

      MPI_Bcast( &Turb->TimeLast, 1, MPI_DOUBLE,   0, MPI_COMM_WORLD );
      MPI_Bcast( &Turb->TimeNext, 1, MPI_DOUBLE,   0, MPI_COMM_WORLD );
      MPI_Bcast( &Turb->RNGState, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD );
      MPI_Bcast( Turb->OUphase[Turb->IdxLast], Turb->NMode*6, MPI_DOUBLE, 0, MPI_COMM_WORLD );
      MPI_Bcast( Turb->OUphase[Turb->IdxNext], Turb->NMode*6, MPI_DOUBLE, 0, MPI_COMM_WORLD );

   } // if ( OPT__INIT == INIT_BY_RESTART && !OPT__RESTART_RESET && !SRC_TURB_RESET )
   else
   {
//    loop through two sets
      for (int t = 0; t < 2 ; t++)
      {
//       construct OU phase vector
         for (int n = 0; n < Turb->NMode; n++)
         {
            double kk       = 0;
            double k_dot_Nr = 0;
            double k_dot_Ni = 0;
            double Nr[3], Ni[3];
            for (int d = 0; d < 3; d++)
            {
//             get random number Nr and Ni
               Turb->GetRNG( Nr[d], Ni[d] );

               kk       += SQR( Turb->Kmode[d][n] );
               k_dot_Nr += Turb->Kmode[d][n]*Nr[d];
               k_dot_Ni += Turb->Kmode[d][n]*Ni[d];
            }

//          Helmholtz decomposition
            for (int d = 0; d < 3; d++)
            {
               Turb->OUphase[t][2*3*n+2*d  ] = SRC_TURB_ZETA*Nr[d] + (1 - 2*SRC_TURB_ZETA)*Turb->Kmode[d][n]*k_dot_Nr/kk;
               Turb->OUphase[t][2*3*n+2*d+1] = SRC_TURB_ZETA*Ni[d] + (1 - 2*SRC_TURB_ZETA)*Turb->Kmode[d][n]*k_dot_Ni/kk;
            }
         } // for (int n = 0; n < Turb->NMode; n++)
      } // for t

//    perform Ornstein-Uhlenbeck process to update OUphase[Next]
      double coeff1 = exp( -Turb->dt/Turb->Tdecay );
      double coeff2 = sqrt( 1 - SQR(coeff1) );
      for (int n = 0; n < Turb->NMode; n++)
      {
         for (int d = 0; d < 3; d++)
         {
            Turb->OUphase[Turb->IdxNext][2*3*n+2*d  ] = coeff1 * Turb->OUphase[Turb->IdxLast][2*3*n+2*d  ] + coeff2 * Turb->OUphase[Turb->IdxNext][2*3*n+2*d  ];
            Turb->OUphase[Turb->IdxNext][2*3*n+2*d+1] = coeff1 * Turb->OUphase[Turb->IdxLast][2*3*n+2*d+1] + coeff2 * Turb->OUphase[Turb->IdxNext][2*3*n+2*d+1];
         }
      }

//    set next update time
      Turb->TimeLast = floor( Time[0]/Turb->dt )*Turb->dt;
      Turb->TimeNext = Turb->TimeLast + Turb->dt;

//    be careful about round-off errors
      if (   (  Turb->TimeNext <= Time[0]  )                                             ||
             (  Time[0] != 0.0 && fabs( (Time[0]-Turb->TimeNext)/Time[0] ) < 1.0e-8   )  ||
             (  Time[0] == 0.0 && fabs(  Time[0]-Turb->TimeNext          ) < 1.0e-12  )      )
      {
         Turb->TimeLast  = Turb->TimeNext;
         Turb->TimeNext += Turb->dt;
      }

   } // !( OPT__INIT == INIT_BY_RESTART && !OPT__RESTART_RESET && !SRC_TURB_RESET )

// initialize acc table, store values on box corner
   const long NPoint = SRC_TURB_TABLE_SIZE + 1;
   const double dx   = amr->BoxSize[0] / SRC_TURB_TABLE_SIZE;
   const double dy   = amr->BoxSize[1] / SRC_TURB_TABLE_SIZE;
   const double dz   = amr->BoxSize[2] / SRC_TURB_TABLE_SIZE;

   for (int d = 0; d < 3; d++)
   {
      Turb->Sin[d] = new double [ NPoint*Turb->NMode ];
      Turb->Cos[d] = new double [ NPoint*Turb->NMode ];
   }

// pre-compute Fourier basis
#  pragma omp parallel for schedule( runtime )
   for (int n = 0; n < Turb->NMode; n++)
   {
      for (int i = 0; i < SRC_TURB_TABLE_SIZE; i++)
      {
         Turb->Sin[0][ n*NPoint + i ] = sin( Turb->Kmode[0][n]*i*dx );
         Turb->Cos[0][ n*NPoint + i ] = cos( Turb->Kmode[0][n]*i*dx );
      }
      for (int j = 0; j < SRC_TURB_TABLE_SIZE; j++)
      {
         Turb->Sin[1][ n*NPoint + j ] = sin( Turb->Kmode[1][n]*j*dy );
         Turb->Cos[1][ n*NPoint + j ] = cos( Turb->Kmode[1][n]*j*dy );
      }
      for (int k = 0; k < SRC_TURB_TABLE_SIZE; k++)
      {
         Turb->Sin[2][ n*NPoint + k ] = sin( Turb->Kmode[2][n]*k*dz );
         Turb->Cos[2][ n*NPoint + k ] = cos( Turb->Kmode[2][n]*k*dz );
      }
//    apply periodicity
      for (int d = 0; d < 3; d++)
      {
         Turb->Sin[d][ n*NPoint + SRC_TURB_TABLE_SIZE ] = Turb->Sin[d][ n*NPoint ];
         Turb->Cos[d][ n*NPoint + SRC_TURB_TABLE_SIZE ] = Turb->Cos[d][ n*NPoint ];
      }
   }

// fill in both tables
   Turb_FillinTable( Turb->IdxLast );
   Turb_FillinTable( Turb->IdxNext );
#  ifdef GPU
   Src_PassData2GPU_Turbulence( Turb->IdxLast );
   Src_PassData2GPU_Turbulence( Turb->IdxNext );
#  endif

// update AuxArray
   Src_SetAuxArray_Turbulence( Src_Turb_AuxArray_Flt, Src_Turb_AuxArray_Int );
#  ifdef GPU
   Src_SetConstMemory_Turbulence( Src_Turb_AuxArray_Flt, Src_Turb_AuxArray_Int,
                                  SrcTerms.Turb_AuxArrayDevPtr_Flt, SrcTerms.Turb_AuxArrayDevPtr_Int );
#  endif

// immediately check update in case DumpTime = Turb->TimeNext
   Turb_CheckUpdate();

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Turb_Init_Field



//-------------------------------------------------------------------------------------------------------
// Function    :  Turb_CheckUpdate
// Description :  Check if Time[0] = Turb->TimeNext -> update Turb field.
//
// Note        :  1. Invoked by  Turb_Init_Field()
//                2. Invoked during main loop after output data
//
// Parameter   :  None
//
// Return      :  Turb
//-------------------------------------------------------------------------------------------------------
void Turb_CheckUpdate()
{
   if (   ( Time[0] != 0.0 && fabs( (Time[0]-Turb->TimeNext)/Time[0] ) < 1.0e-8  )
       || ( Time[0] == 0.0 && fabs(  Time[0]-Turb->TimeNext          ) < 1.0e-12 )   )
   {
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "Time ( %13.7e ) = Turbulence TimeNext ( %13.7e ): Update turbulence pattern ...", Time[0], Turb->TimeNext );

      const double coeff1 = exp( -Turb->dt/Turb->Tdecay );
      const double coeff2 = sqrt( 1 - SQR(coeff1) );

//    swap last and next indices
      std::swap( Turb->IdxLast, Turb->IdxNext );

//    construct OU phase vector
      for (int n = 0; n < Turb->NMode; n++)
      {
         double kk       = 0;
         double k_dot_Nr = 0;
         double k_dot_Ni = 0;
         double Nr[3], Ni[3];

         for (int d = 0; d < 3; d++)
         {
//          get random number Nr and Ni
            Turb->GetRNG( Nr[d], Ni[d] );

            kk       += SQR( Turb->Kmode[d][n] );
            k_dot_Nr += Turb->Kmode[d][n]*Nr[d];
            k_dot_Ni += Turb->Kmode[d][n]*Ni[d];
         }

         for (int d = 0; d < 3; d++)
         {
//          Helmholtz decomposition
            Nr[d] = SRC_TURB_ZETA*Nr[d] + (1 - 2*SRC_TURB_ZETA)*Turb->Kmode[d][n]*k_dot_Nr/kk;
            Ni[d] = SRC_TURB_ZETA*Ni[d] + (1 - 2*SRC_TURB_ZETA)*Turb->Kmode[d][n]*k_dot_Ni/kk;

//          Update OU phases to time_new
            Turb->OUphase[Turb->IdxNext][2*3*n+2*d  ] = coeff1 * Turb->OUphase[Turb->IdxLast][2*3*n+2*d  ] + coeff2 * Nr[d];
            Turb->OUphase[Turb->IdxNext][2*3*n+2*d+1] = coeff1 * Turb->OUphase[Turb->IdxLast][2*3*n+2*d+1] + coeff2 * Ni[d];
         }

      } // for (int n = 0; n < Turb->NMode; n++)

      if ( MPI_Rank == 0 )    Aux_Message( stdout, " done\n" );

//    update turbulence time
      Turb->TimeLast = Turb->TimeNext;
      Turb->TimeNext = round( Time[0]/Turb->dt + 1.0 )*Turb->dt;

//    update new table
      Turb_FillinTable( Turb->IdxNext );
#     ifdef GPU
      Src_PassData2GPU_Turbulence( Turb->IdxNext );
#     endif

//    update AuxArray
      Src_SetAuxArray_Turbulence( Src_Turb_AuxArray_Flt, Src_Turb_AuxArray_Int );
#     ifdef GPU
      Src_SetConstMemory_Turbulence( Src_Turb_AuxArray_Flt, Src_Turb_AuxArray_Int,
                                     SrcTerms.Turb_AuxArrayDevPtr_Flt, SrcTerms.Turb_AuxArrayDevPtr_Int );
#     endif

   } // if ( hasUpdate > 0 )

} // FUNCTION : Turb_CheckUpdate



//-------------------------------------------------------------------------------------------------------
// Function    :  Turb_FillinTable
// Description :  Fillin h_SrcTurb_AccTable host array
//
// Note        :  1. Invoked by Turb_Init(), Turb_CheckUpdate()
//
// Parameter   :  IdxTable : Turbulence table index-> IdxLast or IdxNext
//
// Return      :  h_SrcTurb_AccTable
//-------------------------------------------------------------------------------------------------------
void Turb_FillinTable( int IdxTable )
{
   const long NPoint = SRC_TURB_TABLE_SIZE + 1;

#  pragma omp parallel for schedule( runtime )
   for (int k = 0; k < NPoint; k++)  {
   for (int j = 0; j < NPoint; j++)  {
   for (int i = 0; i < NPoint; i++)  {
      double Acc[3] = {0};
      long  idx = IDX321( i, j, k, NPoint, NPoint );

      for (int n = 0; n < Turb->NMode; n++)
      {
         double sinx = Turb->Sin[0][ n*NPoint + i ];
         double cosx = Turb->Cos[0][ n*NPoint + i ];
         double siny = Turb->Sin[1][ n*NPoint + j ];
         double cosy = Turb->Cos[1][ n*NPoint + j ];
         double sinz = Turb->Sin[2][ n*NPoint + k ];
         double cosz = Turb->Cos[2][ n*NPoint + k ];
         double amp  = Turb->Amplitude[n];

         const double real_part = ( cosx*cosy - sinx*siny ) * cosz - ( sinx*cosy + cosx*siny ) * sinz;
         const double imag_part = ( cosy*sinz + siny*cosz ) * cosx + ( cosy*cosz - siny*sinz ) * sinx;

         for (int d=0; d<3; d++)
            Acc[d] += amp*( Turb->OUphase[IdxTable][2*3*n+2*d]*real_part - Turb->OUphase[IdxTable][2*3*n+2*d+1]*imag_part );

      } // for (int n = 0; n < Turb->NMode; n++)

      for (int d=0; d<3; d++)
         h_SrcTurb_AccTable[IdxTable][3*idx + d] = (real)Acc[d];

   }}} // for i, j, k

} // FUNCTION : Turb_FillinTable


#endif // ifdef TURBULENCE
