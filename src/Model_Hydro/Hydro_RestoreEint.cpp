#include "GAMER.h"

#if ( MODEL == HYDRO )


static real (*EintBk)[PS1][PS1][PS1] = NULL;
static int  EintBk_NP                = 0;

void Hydro_RestoreEint_MemAlloc( const int NP );




//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_RestoreEint_Backup
// Description :  Backup the internal energy of all real patches, which is later used in Hydro_RestoreEint_Check()
//                to restore the original internal energy if unphysical data are detected
//
// Note        :  1. Currently only backup *real* patches
//                2. Internal energy includes cosmic-ray energy if present
//                3. Hydro_RestoreEint_Backup() and Hydro_RestoreEint_Check() must be called
//                   as a paired operation for the same level and fluid/B-field sandglasses
//
// Parameter   :  lv    : Target refinement level
//                FluSg : Sandglass of the fluid   data
//                MagSg : Sandglass of the B field data
//
// Return      :  EintBk[] and EintBk_NP
//-------------------------------------------------------------------------------------------------------
void Hydro_RestoreEint_Backup( const int lv, const int FluSg, const int MagSg )
{

   const int  NP              = amr->NPatchComma[lv][1]; // only target *real* patches
   const bool CheckMinEint_No = false;


// allocate memory if needed
   Hydro_RestoreEint_MemAlloc( NP );


// backup the internal energy of all real patches
#  pragma omp parallel for schedule( runtime )
   for (int PID=0; PID<NP; PID++)
   {
      const real (*fluid)[PS1][PS1][PS1] = amr->patch[FluSg][lv][PID]->fluid;

      for (int k=0; k<PS1; k++)
      for (int j=0; j<PS1; j++)
      for (int i=0; i<PS1; i++)
      {
#        ifdef MHD
         const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, MagSg );
#        else
         const real Emag = NULL_REAL;
#        endif
         const real Eint = Hydro_Con2Eint( fluid[DENS][k][j][i],
                                           fluid[MOMX][k][j][i],
                                           fluid[MOMY][k][j][i],
                                           fluid[MOMZ][k][j][i],
                                           fluid[ENGY][k][j][i],
                                           CheckMinEint_No, NULL_REAL, PassiveFloorMask, Emag,
                                           EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                           EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
         EintBk[PID][k][j][i] = Eint;

//       perform this extra check only in the debug mode, since we have assumed that the *input* internal energy
//       is always physical
//       --> excluding floating-point rounding errors (CK_UNPHY_RND_NO) to avoid false alarms
#        ifdef GAMER_DEBUG
         real fluid_ck[NCOMP_TOTAL];
         for (int v=0; v<NCOMP_TOTAL; v++)   fluid_ck[v] = amr->patch[FluSg][lv][PID]->fluid[v][k][j][i];

         if (  Hydro_IsUnphysical( UNPHY_MODE_CONS, fluid_ck, Emag,
                                   EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                   EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table,
                                   PassiveFloorMask, ERROR_INFO, UNPHY_VERBOSE, CK_UNPHY_RND_NO )  )
            Aux_Error( ERROR_INFO, "unphysical input data detected !!\n" );
#        endif
      } // i, j, k
   } // for (int PID=0; PID<NP; PID++)

} // FUNCTION : Hydro_RestoreEint_Backup



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_RestoreEint_Check
// Description :  Restore the original internal energy if the current conserved variables imply an
//                unphysical internal energy
//
// Note        :  1. Currently only check *real* patches
//                2. Internal energy includes cosmic-ray energy if present
//                3. Buffer patches will be updated as well if any unphysical internal energy is detected
//                4. Restore the original internal energy by modifying only the total energy
//                   --> This routine assumes that density, momentum, passive scalars, and magnetic fields
//                       are already physical, which are thus kept unchanged
//
// Parameter   :  lv    : Target refinement level
//                FluSg : Sandglass of the fluid   data
//                MagSg : Sandglass of the B field data
//
// Return      :  amr->patch[FluSg][lv][*]->fluid[ENGY]
//-------------------------------------------------------------------------------------------------------
void Hydro_RestoreEint_Check( const int lv, const int FluSg, const int MagSg )
{

   const int  NP              = amr->NPatchComma[lv][1]; // only target *real* patches
   const bool CheckMinEint_No = false;

   bool CheckFailed_AnyCell = false;


// checks
#  ifdef GAMER_DEBUG
   if ( NP > EintBk_NP )   Aux_Error( ERROR_INFO, "NP (%d) > EintBk_NP (%d) !!\n", NP, EintBk_NP );
   if ( EintBk == NULL )   Aux_Error( ERROR_INFO, "EintBk == NULL !!\n" );
#  endif


// check and restore the internal energy of all real patches
#  pragma omp parallel for reduction ( ||:CheckFailed_AnyCell ) schedule( runtime )
   for (int PID=0; PID<NP; PID++)
   {
      for (int k=0; k<PS1; k++)
      for (int j=0; j<PS1; j++)
      for (int i=0; i<PS1; i++)
      {
         real fluid[NCOMP_TOTAL];
         for (int v=0; v<NCOMP_TOTAL; v++)   fluid[v] = amr->patch[FluSg][lv][PID]->fluid[v][k][j][i];

#        ifdef MHD
         const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, MagSg );
#        else
         const real Emag = NULL_REAL;
#        endif

//       check unphysical internal energy
//       --> including floating-point rounding errors (CK_UNPHY_RND_YES) to be more cautious
         const bool CheckFailed_ThisCell =
            Hydro_IsUnphysical( UNPHY_MODE_CONS, fluid, Emag,
                                EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table,
                                PassiveFloorMask, ERROR_INFO, UNPHY_SILENCE, CK_UNPHY_RND_YES );

//       restore the original internal energy if needed
         if ( CheckFailed_ThisCell )
         {
//          recompute Etot directly from the backed-up internal energy
//          --> avoid computing the difference between the backed-up and new internal energy,
//              since the latter may be non-finite
            const real EtotNew = Hydro_ConEint2Etot( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ],
                                                     EintBk[PID][k][j][i], Emag );

            amr->patch[FluSg][lv][PID]->fluid[ENGY][k][j][i] = EtotNew;

            CheckFailed_AnyCell = true;

//          perform this extra check only in the debug mode, since we have assumed that the *input* internal energy
//          is always physical
//          --> excluding floating-point rounding errors (CK_UNPHY_RND_NO) to avoid false alarms
#           ifdef GAMER_DEBUG
            fluid[ENGY] = EtotNew;

            if (  Hydro_IsUnphysical( UNPHY_MODE_CONS, fluid, Emag,
                                      EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                      EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table,
                                      PassiveFloorMask, ERROR_INFO, UNPHY_VERBOSE, CK_UNPHY_RND_NO )  )
               Aux_Error( ERROR_INFO, "unphysical output data detected !!\n" );
#           endif
         }
      } // i, j, k
   } // for (int PID=0; PID<NP; PID++)


// since only real patches are checked and corrected above, update buffer-patch energy
// if any correction has been applied on any MPI rank
#  ifndef SERIAL
   MPI_Allreduce( MPI_IN_PLACE, &CheckFailed_AnyCell, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD );

   if ( CheckFailed_AnyCell )
      Buf_GetBufferData( lv, FluSg, NULL_INT, NULL_INT, DATA_GENERAL, _ENGY, _NONE, Flu_ParaBuf, USELB_YES );
#  endif


// print a debug message
   if ( CheckFailed_AnyCell  &&  OPT__VERBOSE  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "\nWARNING : unphysical internal energy detected in %s() !!\n", __FUNCTION__ );

} // FUNCTION : Hydro_RestoreEint_Check



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_RestoreEint_MemAlloc
// Description :  Allocate EintBk[]
//
// Note        :  1. Invoked by Hydro_RestoreEint_Backup()
//
// Parameter   :  NP : Number of patches to be allocated
//-------------------------------------------------------------------------------------------------------
void Hydro_RestoreEint_MemAlloc( const int NP )
{

// allocate memory only if needed
   if ( EintBk_NP < NP )
   {
      delete [] EintBk;

      EintBk_NP = NP;
      EintBk    = new real [EintBk_NP][PS1][PS1][PS1];
   }

} // FUNCTION : Hydro_RestoreEint_MemAlloc



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_RestoreEint_MemFree
// Description :  Free EintBk[]
//
// Note        :  1. Invoked by End_MemFree()
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Hydro_RestoreEint_MemFree()
{

   delete [] EintBk;
   EintBk    = NULL;
   EintBk_NP = 0;

} // FUNCTION : Hydro_RestoreEint_MemFree



#endif // #if ( MODEL == HYDRO )
