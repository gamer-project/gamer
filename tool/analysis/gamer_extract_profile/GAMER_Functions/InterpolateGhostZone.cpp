#include "ExtractProfile.h"

static int Table_01( const int SibID, const int Side, const char dim, const int w01, const int w02,
                     const int w10, const int w11, const int w12, const int w20, const int w21 );




//-------------------------------------------------------------------------------------------------------
// Function    :  InterpolateGhostZone
// Description :  Fill up the ghost-zone values by spatial and temporal interpolation 
//
// Note        :  1. Work for the function "Prepare_PatchData"
//                2. Use the input parameter "NVar_Flu" and "TFluVarIdxList" to control the target variables
//                3. Use the input parameter "IntScheme" to control the interpolation scheme
//                4. Invoke the function "Interpolate" for spatial interpolation
//                5. Data preparation order: FLU -> POTE
// 
// Parameter   :  lv             : Targeted "coarse-grid" refinement level 
//                PID            : Patch ID at level "lv" used for interpolation 
//                IntData        : Array to store the interpolation result 
//                SibID          : Sibling index (0~25) used to determine the interpolation region
//                GhostSize      : Number of ghost zones
//                IntScheme      : Interpolation scheme
//                                 --> currently supported schemes include
//                                     INT_MINMOD1D : MinMod-1D
//                                     INT_MINMOD3D : MinMod-3D
//                                     INT_VANLEER  : vanLeer
//                                     INT_CQUAD    : conservative quadratic
//                                     INT_QUAD     : quadratic
//                                     INT_CQUAR    : conservative quartic
//                                     INT_QUAR     : quartic
//                NTSib          : Number of targeted sibling patches along different sibling directions 
//                TSib           : Targeted sibling indices along different sibling directions
//                TVar           : Targeted variables to be prepared
//                                 --> Supported variables in different models:
//                                     HYDRO : _DENS, _MOMX, _MOMY, _MOMZ, _ENGY, _FLUID, _VELX, _VELY, _VELZ, _PRES,
//                                             [, _POTE] [, _PAR_DENS] [, _PASSIVE]
//                                     MHD   : 
//                                     ELBDM : _DENS, _REAL, _IMAG [, _POTE] [, _PAR_DENS]
//                NVar_Tot       : Total number of variables to be prepared
//                NVar_Flu       : Number of fluid variables to be prepared
//                TFluVarIdxList : List recording the target fluid and passive variable indices ( = [0 ... NCOMP_FLUID+NCOMP_PASSIVE-1] )
//                IntPhase       : true --> Perform interpolation on rho/phase instead of real/imag parts in ELBDM
//-------------------------------------------------------------------------------------------------------
void InterpolateGhostZone( const int lv, const int PID, real IntData[], const int SibID, const int GhostSize, 
                           const IntScheme_t IntScheme, const int NTSib[], int *TSib[], const long TVar, const int NVar_Tot,
                           const int NVar_Flu, const int TFluVarIdxList[], const bool IntPhase )
{

// check
#  ifdef GAMER_DEBUG
// nothing to do if GhostSize == 0
   if ( GhostSize == 0 )   
   {
      Aux_Message( stderr, "WARNING : GhostSize == 0 !!\n" );
      return;
   }
#  endif // #ifdef GAMER_DEBUG


   const bool PrepPot     = ( TVar & _POTE     ) ? true : false;
   const bool PrepParDens = ( TVar & _PAR_DENS ) ? true : false;

// set up parameters for the adopted interpolation scheme
   int NSide, CGhost, CSize[3], FSize[3], CSize3D, FSize3D;

   Int_Table( IntScheme, NSide, CGhost );

   const int  GhostSize_Padded = GhostSize + (GhostSize&1);
   const int  CGrid            = (GhostSize+1)/2 + 2*CGhost;   // number of coarse grids required for interpolation
   const int  La               = CGrid - CGhost;               // number of coarse grids in the central region

   for (int d=0; d<3; d++)    
   {
      CSize[d] = TABLE_01( SibID, 'x'+d, CGrid, PATCH_SIZE+2*CGhost, CGrid );
      FSize[d] = TABLE_01( SibID, 'x'+d, GhostSize_Padded, 2*PATCH_SIZE, GhostSize_Padded );
   }

   CSize3D = CSize[0]*CSize[1]*CSize[2];
   FSize3D = FSize[0]*FSize[1]*FSize[2];


// we assume that we only need ONE coarse-grid patch in each sibling direction
   if ( La > PATCH_SIZE )  Aux_Error( ERROR_INFO, "La (%d) > PATCH_SIZE (%d) !!\n", La, PATCH_SIZE );


// coarse-grid array stored all data required for interpolation (including the ghost zones in each side)
   real *CData_Ptr = NULL;
   real *CData     = new real [ NVar_Tot*CSize3D ];


// 0. initialize CData as some strange numbers to ensure that the patches with fluid==NULL are indeed useless
// ------------------------------------------------------------------------------------------------------------
   for (int t=0; t<NVar_Tot*CSize3D; t++)    CData[t] = __FLT_MAX__;


// a. fill up the central region of CData
// ------------------------------------------------------------------------------------------------------------
   int i1, i2, j1, j2, k1, k2, Idx, TFluVarIdx, Disp1[3], Disp2[3], Loop1[3];

   for (int d=0; d<3; d++)
   {
      Loop1[d] = TABLE_01( SibID, 'x'+d, La, PATCH_SIZE, La );
      Disp1[d] = TABLE_01( SibID, 'x'+d, PATCH_SIZE-La, 0, 0 );
      Disp2[d] = TABLE_01( SibID, 'x'+d, 0, CGhost, CGhost );
   }

   
   CData_Ptr = CData;

// skip the useless patch
   if ( amr.patch[lv][PID]->fluid != NULL )
   {
//    a1. fluid data            
      for (int v=0; v<NVar_Flu; v++)         
      {  
         TFluVarIdx = TFluVarIdxList[v]; 

         for (int k=0; k<Loop1[2]; k++)   {  k1 = k + Disp1[2];   k2 = k + Disp2[2];
         for (int j=0; j<Loop1[1]; j++)   {  j1 = j + Disp1[1];   j2 = j + Disp2[1];
                                             Idx = IDX321( Disp2[0], j2, k2, CSize[0], CSize[1] );
         for (i1=Disp1[0]; i1<Disp1[0]+Loop1[0]; i1++)   {

            CData_Ptr[Idx] = amr.patch[lv][PID]->fluid[TFluVarIdx][k1][j1][i1];

            Idx ++;
         }}}

         CData_Ptr += CSize3D;
      }


//    a2. potential data
      if ( PrepPot )
      {
         for (int k=0; k<Loop1[2]; k++)   {  k1 = k + Disp1[2];   k2 = k + Disp2[2];
         for (int j=0; j<Loop1[1]; j++)   {  j1 = j + Disp1[1];   j2 = j + Disp2[1];
                                             Idx = IDX321( Disp2[0], j2, k2, CSize[0], CSize[1] );
         for (i1=Disp1[0]; i1<Disp1[0]+Loop1[0]; i1++)   {

            CData_Ptr[Idx] = amr.patch[lv][PID]->pot[k1][j1][i1];

            Idx ++;
         }}}

         CData_Ptr += CSize3D;
      }


//    a3. particle density data
      if ( PrepParDens )
      {
         for (int k=0; k<Loop1[2]; k++)   {  k1 = k + Disp1[2];   k2 = k + Disp2[2];
         for (int j=0; j<Loop1[1]; j++)   {  j1 = j + Disp1[1];   j2 = j + Disp2[1];
                                             Idx = IDX321( Disp2[0], j2, k2, CSize[0], CSize[1] );
         for (i1=Disp1[0]; i1<Disp1[0]+Loop1[0]; i1++)   {

            CData_Ptr[Idx] = amr.patch[lv][PID]->par_dens[k1][j1][i1];

            Idx ++;
         }}}

         CData_Ptr += CSize3D;
      }
   } // if ( amr.patch[lv][PID]->fluid != NULL )


// b. fill up the ghost zone of CData
// ------------------------------------------------------------------------------------------------------------
   int Loop2[3], Disp3[3], Disp4[3], Side, SibPID;

   for (int CSib=0; CSib<NTSib[SibID]; CSib++)
   {
      Side   = TSib[SibID][CSib];
      SibPID = amr.patch[lv][PID]->sibling[Side];

      for (int d=0; d<3; d++)
      {
         Loop2[d] = Table_01( SibID, Side, 'x'+d, La, CGhost, CGhost, PATCH_SIZE, CGhost, CGhost, La );
         Disp3[d] = Table_01( SibID, Side, 'x'+d, 0, La, 0, CGhost, CGhost+PATCH_SIZE, 0, CGhost );
         Disp4[d] = Table_01( SibID, Side, 'x'+d, PATCH_SIZE-La, 0, PATCH_SIZE-CGhost, 0, 0, 
                              PATCH_SIZE-CGhost, 0 );
      }


//    b1. if the targeted sibling patch exists --> just copy data from the nearby patch at the same level
      if ( SibPID >= 0 )
      {
         CData_Ptr = CData;

//       skip the useless patch
         if ( amr.patch[lv][SibPID]->fluid != NULL )
         {
//          b1-1. fluid data            
            for (int v=0; v<NVar_Flu; v++)         
            {  
               TFluVarIdx = TFluVarIdxList[v]; 

               for (int k=0; k<Loop2[2]; k++)   {  k1 = k + Disp3[2];   k2 = k + Disp4[2];
               for (int j=0; j<Loop2[1]; j++)   {  j1 = j + Disp3[1];   j2 = j + Disp4[1];
                                                   Idx = IDX321( Disp3[0], j1, k1, CSize[0], CSize[1] );
               for (i2=Disp4[0]; i2<Disp4[0]+Loop2[0]; i2++)   {

                  CData_Ptr[Idx] = amr.patch[lv][SibPID]->fluid[TFluVarIdx][k2][j2][i2];

                  Idx ++;
               }}}

               CData_Ptr += CSize3D;
            }


//          b1-2. potential data
            if ( PrepPot )
            {
               for (int k=0; k<Loop2[2]; k++)   {  k1 = k + Disp3[2];   k2 = k + Disp4[2];
               for (int j=0; j<Loop2[1]; j++)   {  j1 = j + Disp3[1];   j2 = j + Disp4[1];
                                                   Idx = IDX321( Disp3[0], j1, k1, CSize[0], CSize[1] );
               for (i2=Disp4[0]; i2<Disp4[0]+Loop2[0]; i2++)   {

                  CData_Ptr[Idx] = amr.patch[lv][SibPID]->pot[k2][j2][i2];

                  Idx ++;
               }}}

               CData_Ptr += CSize3D;
            }


//          b1-3. particle density data
            if ( PrepParDens )
            {
               for (int k=0; k<Loop2[2]; k++)   {  k1 = k + Disp3[2];   k2 = k + Disp4[2];
               for (int j=0; j<Loop2[1]; j++)   {  j1 = j + Disp3[1];   j2 = j + Disp4[1];
                                                   Idx = IDX321( Disp3[0], j1, k1, CSize[0], CSize[1] );
               for (i2=Disp4[0]; i2<Disp4[0]+Loop2[0]; i2++)   {

                  CData_Ptr[Idx] = amr.patch[lv][SibPID]->par_dens[k2][j2][i2];

                  Idx ++;
               }}}

               CData_Ptr += CSize3D;
            }
         } // if ( amr.patch[lv][SibPID]->fluid != NULL )
      } // if ( SibPID >= 0 )


//    b2. if the targeted sibling patch does not exist --> something is wrong !!
      else if ( SibPID == -1 )
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d Lv %d, PID %d, Side %d) !!\n", "SibPID", SibPID, lv, PID, Side );

      else
         Aux_Error( ERROR_INFO, "SibPID == %d (PID %d, Side %d) !!\n", SibPID, PID, Side );

   } // for (int CSib=0; CSib<NTSib[SibID]; CSib++)



// c. interpolation : CData --> IntData
// ------------------------------------------------------------------------------------------------------------
   const bool PhaseUnwrapping_Yes    = true;
   const bool PhaseUnwrapping_No     = false;
   const bool EnsureMonotonicity_Yes = true;
   const bool EnsureMonotonicity_No  = false;
   int CStart[3], CRange[3], FStart[3], NVar_SoFar;

   for (int d=0; d<3; d++)
   {
      CStart[d] = CGhost;
      CRange[d] = CSize[d] - 2*CGhost;
      FStart[d] = 0;
   }


// determine the variables which must be positive
   bool Monotonicity[NVar_Flu];

   for (int v=0; v<NVar_Flu; v++)
   {
      TFluVarIdx = TFluVarIdxList[v]; 

#     if ( MODEL == HYDRO )
//    we now apply monotonic interpolation to ALL fluid variables (which helps alleviate the issue of negative density/pressure)
      /*
      if ( TFluVarIdx == DENS  ||  TFluVarIdx == ENGY )  Monotonicity[v] = EnsureMonotonicity_Yes;
      else                                               Monotonicity[v] = EnsureMonotonicity_No;
      */
                                                         Monotonicity[v] = EnsureMonotonicity_Yes;

#     elif ( MODEL == MHD )
#     warning : WAIT MHD !!!

#     elif ( MODEL == ELBDM )
      if ( TFluVarIdx == DENS )                          Monotonicity[v] = EnsureMonotonicity_Yes;
      else                                               Monotonicity[v] = EnsureMonotonicity_No;

#     else
#     warning : WARNING : DO YOU WANT TO ENSURE THE POSITIVITY OF INTERPOLATION ??          
#     endif // MODEL
   }


// interpolation
#  if ( MODEL == ELBDM )
   real *CData_Dens = NULL; 
   real *CData_Real = NULL;
   real *CData_Imag = NULL;
   real *FData_Dens = NULL;
   real *FData_Real = NULL;
   real *FData_Imag = NULL;

// c1. interpolation on phase in ELBDM
   if ( IntPhase )
   {
      int DensIdx=-1, RealIdx=-1, ImagIdx=-1;

      for (int v=0; v<NVar_Flu; v++)
      {
         TFluVarIdx = TFluVarIdxList[v]; 

         if      ( TFluVarIdx == DENS )   DensIdx = v;
         else if ( TFluVarIdx == REAL )   RealIdx = v;
         else if ( TFluVarIdx == IMAG )   ImagIdx = v;
      }

//    check
#     ifdef GAMER_DEBUG
      if ( RealIdx == -1  ||  ImagIdx == -1 )
         Aux_Error( ERROR_INFO, "real and/or imag parts are not found for phase interpolation in ELBDM !!\n" );
#     endif


//    determine the array index to store density
      CData_Dens = CData   + ( (DensIdx==-1) ? ImagIdx : DensIdx )*CSize3D;
      CData_Real = CData   + RealIdx*CSize3D;
      CData_Imag = CData   + ImagIdx*CSize3D;
      FData_Dens = IntData + ( (DensIdx==-1) ? ImagIdx : DensIdx )*FSize3D;
      FData_Real = IntData + RealIdx*FSize3D;
      FData_Imag = IntData + ImagIdx*FSize3D;

//    get the wrapped phase (store in the REAL component) and density (store in the IMAG component)
      real Re, Im;

      for (int t=0; t<CSize3D; t++)    
      {
         Re = CData_Real[t];
         Im = CData_Imag[t];

         CData_Real[t] = ATAN2( Im, Re );
         if ( DensIdx == -1 )
         CData_Dens[t] = Re*Re + Im*Im;
      }

//    interpolate density 
      Interpolate( CData_Dens, CSize, CStart, CRange, FData_Dens, FSize, FStart, 1, IntScheme, 
                   PhaseUnwrapping_No, EnsureMonotonicity_Yes );

//    interpolate phase
      Interpolate( CData_Real, CSize, CStart, CRange, FData_Real, FSize, FStart, 1, IntScheme, 
                   PhaseUnwrapping_Yes, EnsureMonotonicity_No );
   } // if ( IntPhase )

// c2. interpolation on real/imag parts in ELBDM
   else // if ( IntPhase )
   {
      for (int v=0; v<NVar_Flu; v++)
      Interpolate( CData+CSize3D*v, CSize, CStart, CRange, IntData+FSize3D*v, FSize, FStart, 1, 
                   IntScheme, PhaseUnwrapping_No, Monotonicity[v] );
   } // if ( IntPhase ) ... else ...

// retrieve real and imaginary parts when phase interpolation is adopted
   if ( IntPhase )
   {
      real Amp, Phase, Rho;

      for (int t=0; t<FSize3D; t++)
      {
         Phase = FData_Real[t];
         Rho   = FData_Dens[t];

//       be careful about the negative density introduced from the round-off errors
//       --> also note that some data are not loaded in GAMER_ExtractProfile (these data are actually useless)
         if ( Rho < (real)0.0 )
         {
            FData_Dens[t] = (real)0.0;
            Rho           = (real)0.0;
         }

         Amp           = SQRT( Rho );
         FData_Real[t] = Amp*COS( Phase );
         FData_Imag[t] = Amp*SIN( Phase );
      }
   }

#  else // #if ( MODEL == ELBDM )

// c3. interpolation on original variables
   for (int v=0; v<NVar_Flu; v++)
      Interpolate( CData+CSize3D*v, CSize, CStart, CRange, IntData+FSize3D*v, FSize, FStart, 1, 
                   IntScheme, PhaseUnwrapping_No, Monotonicity[v] );

#  endif // #if ( MODEL == ELBDM ) ... else 

   NVar_SoFar = NVar_Flu;

// c4. interpolation on potential
   if ( PrepPot )
   {
      Interpolate( CData+CSize3D*NVar_SoFar, CSize, CStart, CRange, IntData+FSize3D*NVar_SoFar, FSize, FStart, 1, 
                   IntScheme, PhaseUnwrapping_No, EnsureMonotonicity_No );
      NVar_SoFar ++;
   }

// c5. interpolation on particle density
   if ( PrepParDens )
   {
      Interpolate( CData+CSize3D*NVar_SoFar, CSize, CStart, CRange, IntData+FSize3D*NVar_SoFar, FSize, FStart, 1, 
                   IntScheme, PhaseUnwrapping_No, EnsureMonotonicity_No );
      NVar_SoFar ++;
   }

   delete [] CData;

} // FUNCTION : InterpolateGhostZone



// ============
// |  Tables  | 
// ============

//-------------------------------------------------------------------------------------------------------
// Function    :  Table_01
// Description :  return the loop size and displacement required by the function "InterpolateGhostZone"
// 
// Parameter   :  SibID       : Sibling index ONE (0~25)
//                Side        : Sibling index TWO (0~25)
//                dim         : Targeted spatial direction (x/y/z)
//                w01 ... w21 : Returned values  
//-------------------------------------------------------------------------------------------------------
int Table_01( const int SibID, const int Side, const char dim, const int w01, const int w02,
              const int w10, const int w11, const int w12, const int w20, const int w21 )
{

   switch ( dim )
   {
      case 'x':
      {
         switch ( SibID )
         {
            case 0: case 6: case 8: case 14: case 15: case 18: case 20: case 22: case 24:
            {
               switch ( Side )
               {
                  case 0: case 6: case 8: case 14: case 15: case 18: case 20: case 22: case 24:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n", 
                                "SibID", SibID, "Side", Side );

                  case 2: case 3: case 4: case 5: case 10: case 11: case 12: case 13:
                     return w01;

                  case 1: case 7: case 9: case 16: case 17: case 19: case 21: case 23: case 25:
                     return w02;
               }
            }

            case 2: case 3: case 4: case 5: case 10: case 11: case 12: case 13:
            {
               switch ( Side )
               {
                  case 0: case 6: case 8: case 14: case 15: case 18: case 20: case 22: case 24:
                     return w10;

                  case 2: case 3: case 4: case 5: case 10: case 11: case 12: case 13:
                     return w11;

                  case 1: case 7: case 9: case 16: case 17: case 19: case 21: case 23: case 25:
                     return w12;
               }
            }

            case 1: case 7: case 9: case 16: case 17: case 19: case 21: case 23: case 25:
            {
               switch ( Side )
               {
                  case 0: case 6: case 8: case 14: case 15: case 18: case 20: case 22: case 24:
                     return w20;

                  case 2: case 3: case 4: case 5: case 10: case 11: case 12: case 13:
                     return w21;

                  case 1: case 7: case 9: case 16: case 17: case 19: case 21: case 23: case 25:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n", 
                                "SibID", SibID, "Side", Side );
               }
            }

         } // switch ( SibID )
      } // case 'x':


      case 'y':
      {
         switch ( SibID )
         {
            case 2: case 6: case 7: case 10: case 12: case 18: case 19: case 22: case 23:
            {
               switch ( Side )
               {
                  case 2: case 6: case 7: case 10: case 12: case 18: case 19: case 22: case 23:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n", 
                                "SibID", SibID, "Side", Side );

                  case 0: case 1: case 4: case 5: case 14: case 15: case 16: case 17:
                     return w01;

                  case 3: case 8: case 9: case 11: case 13: case 20: case 21: case 24: case 25:
                     return w02;
               }
            }

            case 0: case 1: case 4: case 5: case 14: case 15: case 16: case 17:
            {
               switch ( Side )
               {
                  case 2: case 6: case 7: case 10: case 12: case 18: case 19: case 22: case 23:
                     return w10;

                  case 0: case 1: case 4: case 5: case 14: case 15: case 16: case 17:
                     return w11;

                  case 3: case 8: case 9: case 11: case 13: case 20: case 21: case 24: case 25:
                     return w12;
               }
            }

            case 3: case 8: case 9: case 11: case 13: case 20: case 21: case 24: case 25:
            {
               switch ( Side )
               {
                  case 2: case 6: case 7: case 10: case 12: case 18: case 19: case 22: case 23:
                     return w20;

                  case 0: case 1: case 4: case 5: case 14: case 15: case 16: case 17:
                     return w21;

                  case 3: case 8: case 9: case 11: case 13: case 20: case 21: case 24: case 25:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n", 
                                "SibID", SibID, "Side", Side );
               }
            }

         } // switch ( SibID )
      } // case 'y':


      case 'z':
      {
         switch ( SibID )
         {
            case 4: case 10: case 11: case 14: case 16: case 18: case 19: case 20: case 21:
            {
               switch ( Side )
               {
                  case 4: case 10: case 11: case 14: case 16: case 18: case 19: case 20: case 21:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n", 
                                "SibID", SibID, "Side", Side );

                  case 0: case 1: case 2: case 3: case 6: case 7: case 8: case 9:
                     return w01;

                  case 5: case 12: case 13: case 15: case 17: case 22: case 23: case 24: case 25:
                     return w02;
               }
            }

            case 0: case 1: case 2: case 3: case 6: case 7: case 8: case 9:
            {
               switch ( Side )
               {
                  case 4: case 10: case 11: case 14: case 16: case 18: case 19: case 20: case 21:
                     return w10;

                  case 0: case 1: case 2: case 3: case 6: case 7: case 8: case 9:
                     return w11;

                  case 5: case 12: case 13: case 15: case 17: case 22: case 23: case 24: case 25:
                     return w12;
               }
            }

            case 5: case 12: case 13: case 15: case 17: case 22: case 23: case 24: case 25:
            {
               switch ( Side )
               {
                  case 4: case 10: case 11: case 14: case 16: case 18: case 19: case 20: case 21:
                     return w20;

                  case 0: case 1: case 2: case 3: case 6: case 7: case 8: case 9:
                     return w21;

                  case 5: case 12: case 13: case 15: case 17: case 22: case 23: case 24: case 25:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n", 
                                "SibID", SibID, "Side", Side );
               }
            }

         } // switch ( SibID )
      } // case 'z':


      default:
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %c !!\n", "dim", dim );
         exit(1);

   } // switch ( dim )

   return NULL_INT;

} // FUNCTION : Table_01
