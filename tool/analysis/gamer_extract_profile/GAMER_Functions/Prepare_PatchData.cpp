#include "ExtractProfile.h"

void InterpolateGhostZone( const int lv, const int PID, real IntData[], const int SibID, const int GhostSize, 
                           const IntScheme_t IntScheme, const int NTSib[], int *TSib[], const int TVar, const int NVar_Tot,
                           const int NVar_Flu, const int TFluVarIdxList[], const bool IntPhase );
static void SetTargetSibling( int NTSib[], int *TSib[] );
static int Table_01( const int SibID, const char dim, const int Count, const int GhostSize );
static int Table_02( const int lv, const int PID, const int Side );




//-------------------------------------------------------------------------------------------------------
// Function    :  Prepare_PatchData
// Description :  Prepare the uniform data including ghost zones for the target patches or patch groups
//
// Note        :  1. Use the input parameter "TVar" to control the target variables
//                   --> TVar can be any combination of the symbolic constants defined in "Macro.h" 
//                       (e.g., "TVar = _DENS", "TVar = _MOMX|ENGY", or "TVar = _FLUID|_POTE")
//                2. If "GhostSize != 0" --> the function "InterpolateGhostZone" will be used to fill up the
//                   ghost-zone values by spatial interpolation if the corresponding sibling patches do
//                   NOT exist
//                3. Use "patch group" as the preparation unit
//                   --> The data of all patches within the same patch group will be prepared 
//                4. Data preparation order: FLU --> POTE
//
// Parameter   :  lv             : Targeted refinement level
//                h_Input_Array  : Host array to store the prepared data
//                GhostSize      : Number of ghost zones to be prepared
//                NPG            : Number of patch groups prepared at a time
//                PID0_List      : List recording the patch indices with LocalID==0 to be prepared
//                TVar           : Targeted variables to be prepared
//                                 --> Supported variables in different models:
//                                     HYDRO : _DENS, _MOMX, _MOMY, _MOMZ, _ENGY, _FLUID, _VELX, _VELY, _VELZ, _PRES,
//                                             [, _POTE] [, _PAR_DENS] [, _PASSIVE]
//                                     MHD   : 
//                                     ELBDM : _DENS, _REAL, _IMAG [, _POTE] [, _PAR_DENS]
//                IntScheme      : Interpolation scheme
//                                 --> currently supported schemes include
//                                     INT_MINMOD1D : MinMod-1D
//                                     INT_MINMOD3D : MinMod-3D
//                                     INT_VANLEER  : vanLeer
//                                     INT_CQUAD    : conservative quadratic
//                                     INT_QUAD     : quadratic
//                                     INT_CQUAR    : conservative quartic
//                                     INT_QUAR     : quartic
//                NSide          : Number of sibling directions to prepare data
//                                 --> NSIDE_06 (=  6) : prepare only sibling directions 0~5
//                                     NSIDE_26 (= 26) : prepare all sibling directions 0~25
//                IntPhase       : true --> Perform interpolation on rho/phase instead of real/imag parts in ELBDM
//                                      --> TVar must contain _REAL and _IMAG
//-------------------------------------------------------------------------------------------------------
void Prepare_PatchData( const int lv, real *h_Input_Array, const int GhostSize, const int NPG, const int *PID0_List, 
                        const int TVar, const IntScheme_t IntScheme, const NSide_t NSide, const bool IntPhase )
{

// check
#  ifdef GAMER_DEBUG

   if ( TVar & ~(_TOTAL|_POTE|_PAR_DENS) )   Aux_Error( ERROR_INFO, "unsupported parameter %s = %d !!\n", "TVar", TVar );

   if ( IntPhase )
   {
#     if ( MODEL == ELBDM )
      if (  !(TVar & _REAL)  ||  !(TVar & _IMAG)  )
      Aux_Error( ERROR_INFO, "real and/or imag parts are not found for phase interpolation in ELBDM !!\n" );
#     else
      Aux_Error( ERROR_INFO, "\"interpolation on phase\" is useful only in ELBDM !!\n" );
#     endif
   }
#  endif // #ifdef GAMER_DEBUG


   const int  PGSize1D    = 2*( PATCH_SIZE + GhostSize );    // size of a single patch group including the ghost zone
   const int  PGSize3D    = PGSize1D*PGSize1D*PGSize1D;
   const bool PrepPot     = ( TVar & _POTE     ) ? true : false;
   const bool PrepParDens = ( TVar & _PAR_DENS ) ? true : false;

// TFluVarIdxList : List recording the targeted fluid and passive variable indices ( = [0 ... NCOMP_TOTAL-1] )
   int NTSib[26], *TSib[26], NVar_Flu, NVar_Tot, TFluVarIdxList[NCOMP_TOTAL];

// set up the targeted sibling indices for the function "InterpolateGhostZone"
   SetTargetSibling( NTSib, TSib );

// determine the components to be prepared
   NVar_Flu = 0;
   for (int v=0; v<NCOMP_TOTAL; v++)
      if ( TVar & (1<<v) )    TFluVarIdxList[ NVar_Flu++ ] = v;

   NVar_Tot = NVar_Flu;

   if ( PrepPot )       NVar_Tot ++; 
   if ( PrepParDens )   NVar_Tot ++; 


   if ( NVar_Tot == 0 )    
   {
      Aux_Message( stderr, "WARNING : no targeted variable is found !!\n" );
      return;
   }


// start to prepare data
//#  pragma omp parallel
   {
      int  J, K, I2, J2, K2, Idx1, Idx2, PID0, TFluVarIdx;

//    Array : array to store the prepared data of one patch group (including the ghost-zone data) 
      real *Array_Ptr = NULL;
      real *Array     = new real [ NVar_Tot*PGSize3D ];

      
//    prepare eight nearby patches (one patch group) at a time 
//#     pragma omp for schedule( runtime )
      for (int TID=0; TID<NPG; TID++)
      {
         PID0 = PID0_List[TID];

#        ifdef GAMER_DEBUG
         if ( PID0 < 0  ||  PID0 >= amr.num[lv] )  
            Aux_Error( ERROR_INFO, "PID0 (%d) is not within the correct range [%d ... %d] !!\n", PID0, 0, amr.num[lv]-1 );

         if ( PID0%8 != 0 )
            Aux_Error( ERROR_INFO, "PID0 (%d) does not have LocalID == 0 !!\n", PID0 );
#        endif



//       0. initialize Array as some strange numbers to ensure that the patches with fluid==NULL are indeed useless
// ------------------------------------------------------------------------------------------------------------
         for (int t=0; t<NVar_Tot*PGSize3D; t++)   Array[t] = __FLT_MAX__;
         

         
//       a. fill up the central region of Array (ghost zone is not filled up yet)
// ------------------------------------------------------------------------------------------------------------
         for (int LocalID=0; LocalID<8; LocalID++ )
         {
            const int PID    = PID0 + LocalID;
            const int Disp_i = TABLE_02( LocalID, 'x', GhostSize, GhostSize+PATCH_SIZE );
            const int Disp_j = TABLE_02( LocalID, 'y', GhostSize, GhostSize+PATCH_SIZE );
            const int Disp_k = TABLE_02( LocalID, 'z', GhostSize, GhostSize+PATCH_SIZE );

//          skip the useless patch
            if ( amr.patch[lv][PID]->fluid == NULL )  continue;

            Array_Ptr = Array;
            
//          a1. fluid data
            for (int v=0; v<NVar_Flu; v++)
            {  
               TFluVarIdx = TFluVarIdxList[v];

               for (int k=0; k<PATCH_SIZE; k++)    {  K    = k + Disp_k;
               for (int j=0; j<PATCH_SIZE; j++)    {  J    = j + Disp_j;
                                                      Idx1 = IDX321( Disp_i, J, K, PGSize1D, PGSize1D );
               for (int i=0; i<PATCH_SIZE; i++)    {
               
                  Array_Ptr[ Idx1 ++ ] = amr.patch[lv][PID]->fluid[TFluVarIdx][k][j][i];
               
               }}}

               Array_Ptr += PGSize3D;
            }


//          a2. potential data
            if ( PrepPot )
            {
               for (int k=0; k<PATCH_SIZE; k++)    {  K    = k + Disp_k;
               for (int j=0; j<PATCH_SIZE; j++)    {  J    = j + Disp_j;
                                                      Idx1 = IDX321( Disp_i, J, K, PGSize1D, PGSize1D );
               for (int i=0; i<PATCH_SIZE; i++)    {
               
                  Array_Ptr[ Idx1 ++ ] = amr.patch[lv][PID]->pot[k][j][i];
               
               }}}

               Array_Ptr += PGSize3D;
            }


//          a3. particle density data
            if ( PrepParDens )
            {
               for (int k=0; k<PATCH_SIZE; k++)    {  K    = k + Disp_k;
               for (int j=0; j<PATCH_SIZE; j++)    {  J    = j + Disp_j;
                                                      Idx1 = IDX321( Disp_i, J, K, PGSize1D, PGSize1D );
               for (int i=0; i<PATCH_SIZE; i++)    {
               
                  Array_Ptr[ Idx1 ++ ] = amr.patch[lv][PID]->par_dens[k][j][i];
               
               }}}

               Array_Ptr += PGSize3D;
            }
         } // for (int LocalID=0; LocalID<8; LocalID++ )


//       b. fill up the ghost zone of Array
// ------------------------------------------------------------------------------------------------------------
         for (int Side=0; Side<NSide; Side++)
         {
//          nothing to do if no ghost zone is required
            if ( GhostSize == 0 )   break;


            const int SibPID0 = Table_02( lv, PID0, Side );    // the 0th patch of the sibling patch group

//          (b1) if the targeted sibling patch exists --> just copy data from the nearby patches at the same level
            if ( SibPID0 >= 0 )
            {
               const int Loop_i  = TABLE_01( Side, 'x', GhostSize, PATCH_SIZE, GhostSize );
               const int Loop_j  = TABLE_01( Side, 'y', GhostSize, PATCH_SIZE, GhostSize );
               const int Loop_k  = TABLE_01( Side, 'z', GhostSize, PATCH_SIZE, GhostSize );
               const int Disp_i2 = TABLE_01( Side, 'x', PATCH_SIZE-GhostSize, 0, 0 );
               const int Disp_j2 = TABLE_01( Side, 'y', PATCH_SIZE-GhostSize, 0, 0 );
               const int Disp_k2 = TABLE_01( Side, 'z', PATCH_SIZE-GhostSize, 0, 0 );
                  
//###OPTIMIZATION: simplify TABLE_03 and TABLE_04
               for (int Count=0; Count<TABLE_04( Side ); Count++)
               {
                  const int SibPID = TABLE_03( Side, Count ) + SibPID0;
                  const int Disp_i = Table_01( Side, 'x', Count, GhostSize );
                  const int Disp_j = Table_01( Side, 'y', Count, GhostSize );
                  const int Disp_k = Table_01( Side, 'z', Count, GhostSize );

//                skip the useless patch
                  if ( amr.patch[lv][SibPID]->fluid == NULL )  continue;
                  
                  Array_Ptr = Array;

//                (b1-1) fluid data
                  for (int v=0; v<NVar_Flu; v++)
                  {
                     TFluVarIdx = TFluVarIdxList[v];

                     for (int k=0; k<Loop_k; k++)  {  K = k + Disp_k;   K2 = k + Disp_k2;
                     for (int j=0; j<Loop_j; j++)  {  J = j + Disp_j;   J2 = j + Disp_j2;
                                                      Idx1 = IDX321( Disp_i, J, K, PGSize1D, PGSize1D );
                     for (I2=Disp_i2; I2<Disp_i2+Loop_i; I2++) {

                        Array_Ptr[ Idx1 ++ ] = amr.patch[lv][SibPID]->fluid[TFluVarIdx][K2][J2][I2];

                     }}}

                     Array_Ptr += PGSize3D;
                  }


//                (b1-2) potential data
                  if ( PrepPot )
                  {
                     for (int k=0; k<Loop_k; k++)  {  K = k + Disp_k;   K2 = k + Disp_k2;
                     for (int j=0; j<Loop_j; j++)  {  J = j + Disp_j;   J2 = j + Disp_j2;
                                                      Idx1 = IDX321( Disp_i, J, K, PGSize1D, PGSize1D );
                     for (I2=Disp_i2; I2<Disp_i2+Loop_i; I2++) {

                        Array_Ptr[ Idx1 ++ ] = amr.patch[lv][SibPID]->pot[K2][J2][I2];

                     }}}

                     Array_Ptr += PGSize3D;
                  }


//                (b1-2) particle density data
                  if ( PrepParDens )
                  {
                     for (int k=0; k<Loop_k; k++)  {  K = k + Disp_k;   K2 = k + Disp_k2;
                     for (int j=0; j<Loop_j; j++)  {  J = j + Disp_j;   J2 = j + Disp_j2;
                                                      Idx1 = IDX321( Disp_i, J, K, PGSize1D, PGSize1D );
                     for (I2=Disp_i2; I2<Disp_i2+Loop_i; I2++) {

                        Array_Ptr[ Idx1 ++ ] = amr.patch[lv][SibPID]->par_dens[K2][J2][I2];

                     }}}

                     Array_Ptr += PGSize3D;
                  }
               } // for (int Count=0; Count<TABLE_04( Side ); Count++)
            } // if ( SibPID0 >= 0 )  


//          (b2) if the targeted sibling patch does not exist --> interpolate from patches at level lv-1
            else if ( SibPID0 == -1 )
            {
//             interpolation should never be applied to the base level
#              ifdef GAMER_DEBUG
               if ( lv == 0 )    Aux_Error( ERROR_INFO, "performing interpolation at the base level !!\n" );
#              endif


//             allocate the array to store the interpolation result
               const int GhostSize_Padded = GhostSize + (GhostSize&1);

               int FSize[3];
               for (int d=0; d<3; d++)  FSize[d] = TABLE_01( Side, 'x'+d, GhostSize_Padded, 2*PATCH_SIZE, GhostSize_Padded );

               real *IntData_Ptr = NULL;
               real *IntData     = new real [ NVar_Tot*FSize[0]*FSize[1]*FSize[2] ];


//             determine the parameters for the spatial and temporal interpolations
               const int FaPID    = amr.patch[lv][PID0]->father;
               const int FaSibPID = amr.patch[lv-1][FaPID]->sibling[Side];

#              ifdef GAMER_DEBUG
               if ( FaSibPID < 0 )  Aux_Error( ERROR_INFO, "FaSibPID = %d < 0 (lv %d, PID0 %d, FaPID %d, sib %d) !!\n",
                                               FaSibPID, lv, PID0, FaPID, Side ); 
#              endif


//             perform interpolation and store the results in IntData
               InterpolateGhostZone( lv-1, FaSibPID, IntData, Side, GhostSize, IntScheme, NTSib, TSib, TVar, NVar_Tot, 
                                     NVar_Flu, TFluVarIdxList, IntPhase );


//             properly copy data from IntData array to Array
               const int NUseless = GhostSize & 1;
               const int Loop_i   = TABLE_01( Side, 'x', GhostSize, 2*PATCH_SIZE, GhostSize );
               const int Loop_j   = TABLE_01( Side, 'y', GhostSize, 2*PATCH_SIZE, GhostSize );
               const int Loop_k   = TABLE_01( Side, 'z', GhostSize, 2*PATCH_SIZE, GhostSize );
               const int Disp_i1  = TABLE_01( Side, 'x', 0, GhostSize, GhostSize+2*PATCH_SIZE );
               const int Disp_j1  = TABLE_01( Side, 'y', 0, GhostSize, GhostSize+2*PATCH_SIZE );
               const int Disp_k1  = TABLE_01( Side, 'z', 0, GhostSize, GhostSize+2*PATCH_SIZE );
               const int Disp_i2  = TABLE_01( Side, 'x', NUseless, 0, 0 );
               const int Disp_j2  = TABLE_01( Side, 'y', NUseless, 0, 0 );
               const int Disp_k2  = TABLE_01( Side, 'z', NUseless, 0, 0 );

               Array_Ptr   = Array;
               IntData_Ptr = IntData;

               for (int v=0; v<NVar_Tot; v++)
               {
                  for (int k=0; k<Loop_k; k++)  {  K = k + Disp_k1;  K2 = k + Disp_k2;
                  for (int j=0; j<Loop_j; j++)  {  J = j + Disp_j1;  J2 = j + Disp_j2;
                                                   Idx1 = IDX321( Disp_i1, J,  K,  PGSize1D, PGSize1D );
                                                   Idx2 = IDX321( Disp_i2, J2, K2, FSize[0], FSize[1] );
                  for (int i=0; i<Loop_i; i++)  {

                     Array_Ptr[ Idx1 ++ ] = IntData_Ptr[ Idx2 ++ ];

                  }}}

                  Array_Ptr   += PGSize3D;
                  IntData_Ptr += FSize[0]*FSize[1]*FSize[2];
               }

               delete [] IntData;

            } // else if ( SibPID0 == -1 )


            else
               Aux_Error( ERROR_INFO, "SibPID0 == %d (PID0 %d, Side %d) !!\n", SibPID0, PID0, Side );

         } // for (int Side=0; Side<NSide; Side++)
            

//       c. copy data from Array to h_Input_Array
// ------------------------------------------------------------------------------------------------------------
//       if ( PrepUnit == UNIT_PATCH ) // separate the prepared patch group data into individual patches
         if ( true )                   // always prepare individual patches 
         {
            const int PSize1D = PATCH_SIZE + 2*GhostSize;  // size of a single patch including the ghost zone
            const int PSize3D = PSize1D*PSize1D*PSize1D;
            real *InArray_Ptr = NULL;

            for (int LocalID=0; LocalID<8; LocalID++)
            {
               const int N      = 8*TID + LocalID;
               const int Disp_i = TABLE_02( LocalID, 'x', 0, PATCH_SIZE );
               const int Disp_j = TABLE_02( LocalID, 'y', 0, PATCH_SIZE );
               const int Disp_k = TABLE_02( LocalID, 'z', 0, PATCH_SIZE );

               Array_Ptr   = Array;
               InArray_Ptr = h_Input_Array + N*NVar_Tot*PSize3D;
               Idx2        = 0;
               
               for (int v=0; v<NVar_Tot; v++)
               {
                  for (int k=Disp_k; k<Disp_k+PSize1D; k++)
                  for (int j=Disp_j; j<Disp_j+PSize1D; j++)
                  {
                     Idx1 = IDX321( Disp_i, j, k, PGSize1D, PGSize1D );

                     for (int i=0; i<PSize1D; i++)    InArray_Ptr[ Idx2 ++ ] = Array_Ptr[ Idx1 ++ ];
                  }

                  Array_Ptr += PGSize3D;
               }
            }
         } // if ( PatchByPatch )

         /*
         else if ( PrepUnit == UNIT_PATCHGROUP )
            memcpy( h_Input_Array + TID*NVar_Tot*PGSize3D, Array, NVar_Tot*PGSize3D*sizeof(real) );

         else
            Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "PrepUnit", PrepUnit );
         */

      } // for (int TID=0; TID<NPG; TID++)

      delete [] Array;

   } // OpenMP parallel region


// free memroy
   for (int s=0; s<26; s++)   delete [] TSib[s];

} // FUNCTION : Prepare_PatchData



// ============
// |  Tables  | 
// ============

//-------------------------------------------------------------------------------------------------------
// Function    :  Table_01 
// Description :  Return the displacement for the function "Prepare_PatchData"
//
// Parameter   :  SibID     : Sibling index (0~25) 
//                dim       : Targeted spatial direction (x/y/z)
//                Count     : Patch counter (0~3)
//                GhostSize : Number of ghost zones
//-------------------------------------------------------------------------------------------------------
int Table_01( const int SibID, const char dim, const int Count, const int GhostSize )
{

   switch ( dim )
   {
      case 'x':
      {
         switch ( SibID )
         {
            case 0:case 6: case 8: case 14: case 15: case 18: case 20: case 22: case 24:
               return 0;
               
            case 2: case 3:
            {
               switch ( Count )
               {
                  case 0: case 2:   return GhostSize;
                  case 1: case 3:   return GhostSize + PATCH_SIZE;
                  default:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n", 
                                "SibID", SibID, "Count", Count );
               }
            }
               
            case 4: case 5:
            {
               switch ( Count )
               {
                  case 0: case 1:   return GhostSize;
                  case 2: case 3:   return GhostSize + PATCH_SIZE;
                  default:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n", 
                                "SibID", SibID, "Count", Count );
               }
            }
               
            case 10: case 11: case 12: case 13:
            {
               switch ( Count )
               {
                  case 0:  return GhostSize;
                  case 1:  return GhostSize + PATCH_SIZE;
                  default:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n", 
                                "SibID", SibID, "Count", Count );
               }
            }
               
            case 1: case 7: case 9: case 16: case 17: case 19: case 21: case 23: case 25:
               return GhostSize + 2*PATCH_SIZE;

            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "SibID", SibID );

         } // switch ( SibID )
      } // case 'x':


      case 'y':
      {
         switch ( SibID )
         {
            case 2: case 6: case 7: case 10: case 12: case 18: case 19: case 22: case 23:
               return 0;
               
            case 0: case 1:
            {
               switch ( Count )
               {
                  case 0: case 1:   return GhostSize;
                  case 2: case 3:   return GhostSize + PATCH_SIZE;
                  default:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n", 
                                "SibID", SibID, "Count", Count );
               }
            }

            case 4: case 5:
            {
               switch ( Count )
               {
                  case 0: case 2:   return GhostSize;
                  case 1: case 3:   return GhostSize + PATCH_SIZE;
                  default:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n", 
                                "SibID", SibID, "Count", Count );
               }
            }

            case 14: case 15: case 16: case 17:
            {
               switch ( Count ) 
               {
                  case 0:  return GhostSize;
                  case 1:  return GhostSize + PATCH_SIZE;
                  default:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n", 
                                "SibID", SibID, "Count", Count );
               }
            }
                  
            case 3: case 8: case 9: case 11: case 13: case 20: case 21: case 24: case 25:
               return GhostSize + 2*PATCH_SIZE; 

            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "SibID", SibID );

         } // switch ( SibID )
      } // case 'y':


      case 'z':
      {
         switch ( SibID )
         {
            case 4: case 10: case 11: case 14: case 16: case 18: case 19: case 20: case 21:
               return 0;
               
            case 0: case 1:
            {
               switch ( Count )
               {
                  case 0: case 2:   return GhostSize;
                  case 1: case 3:   return GhostSize + PATCH_SIZE;
                  default:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n", 
                                "SibID", SibID, "Count", Count );
               }
            }
               
            case 2: case 3:
            {
               switch ( Count )
               {
                  case 0: case 1:   return GhostSize;
                  case 2: case 3:   return GhostSize + PATCH_SIZE;                     
                  default:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n", 
                                "SibID", SibID, "Count", Count );
               }
            }
               
            case 6: case 7: case 8: case 9:
            {
               switch ( Count )
               {
                  case 0:  return GhostSize;
                  case 1:  return GhostSize + PATCH_SIZE;
                  default:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n", 
                                "SibID", SibID, "Count", Count );
               }
            }
                  
            case 5: case 12: case 13: case 15: case 17: case 22: case 23: case 24: case 25:
               return GhostSize + 2*PATCH_SIZE;

            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "SibID", SibID );

         } // switch ( SibID )
      } // case 'z':


      default:
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %c !!\n", "dim", dim );
         exit(1);
   } // switch ( dim )

   return NULL_INT;

} // FUNCTION : Table_01



//-------------------------------------------------------------------------------------------------------
// Function    :  Table_02 
// Description :  Return the patch ID of the 0th patch (local ID = 0) of the sibling patch group 
//
// Note        :  Work for the function "Prepare_PatchData" 
//
// Parameter   :  lv    : Targeted refinement level 
//                PID   : Targeted patch ID to find its sibling patches
//                Side  : Sibling index (0~25) 
//-------------------------------------------------------------------------------------------------------
int Table_02( const int lv, const int PID, const int Side )
{

   int Sib;

   switch ( Side )
   {
      case 0: 
         Sib = amr.patch[lv][PID  ]->sibling[0];
         if ( Sib >= 0 )  return Sib-1;
         else             return Sib;

      case 1:
         Sib = amr.patch[lv][PID+1]->sibling[1];
         if ( Sib >= 0 )  return Sib;
         else             return Sib;

      case 2:
         Sib = amr.patch[lv][PID  ]->sibling[2];
         if ( Sib >= 0 )  return Sib-2;
         else             return Sib;

      case 3:
         Sib = amr.patch[lv][PID+2]->sibling[3];
         if ( Sib >= 0 )  return Sib;
         else             return Sib;

      case 4:
         Sib = amr.patch[lv][PID  ]->sibling[4];
         if ( Sib >= 0 )  return Sib-3;
         else             return Sib;

      case 5: 
         Sib = amr.patch[lv][PID+3]->sibling[5];
         if ( Sib >= 0 )  return Sib;
         else             return Sib;

      case 6:
         Sib = amr.patch[lv][PID  ]->sibling[6];
         if ( Sib >= 0 )  return Sib-4;
         else             return Sib;

      case 7:
         Sib = amr.patch[lv][PID+1]->sibling[7];
         if ( Sib >= 0 )  return Sib-2;
         else             return Sib;

      case 8:
         Sib = amr.patch[lv][PID+2]->sibling[8];
         if ( Sib >= 0 )  return Sib-1;
         else             return Sib;

      case 9:
         Sib = amr.patch[lv][PID+4]->sibling[9];
         if ( Sib >= 0 )  return Sib;
         else             return Sib;

      case 10: 
         Sib = amr.patch[lv][PID  ]->sibling[10];
         if ( Sib >= 0 )  return Sib-5;
         else             return Sib;

      case 11:
         Sib = amr.patch[lv][PID+2]->sibling[11];
         if ( Sib >= 0 )  return Sib-3;
         else             return Sib;

      case 12:
         Sib = amr.patch[lv][PID+3]->sibling[12];
         if ( Sib >= 0 )  return Sib-2;
         else             return Sib;

      case 13:
         Sib = amr.patch[lv][PID+5]->sibling[13];
         if ( Sib >= 0 )  return Sib;
         else             return Sib;

      case 14:
         Sib = amr.patch[lv][PID  ]->sibling[14];
         if ( Sib >= 0 )  return Sib-6;
         else             return Sib;

      case 15: 
         Sib = amr.patch[lv][PID+3]->sibling[15];
         if ( Sib >= 0 )  return Sib-1;
         else             return Sib;

      case 16:
         Sib = amr.patch[lv][PID+1]->sibling[16];
         if ( Sib >= 0 )  return Sib-3;
         else             return Sib;

      case 17:
         Sib = amr.patch[lv][PID+6]->sibling[17];
         if ( Sib >= 0 )  return Sib;
         else             return Sib;

      case 18:
         Sib = amr.patch[lv][PID  ]->sibling[18];
         if ( Sib >= 0 )  return Sib-7;
         else             return Sib;

      case 19:
         Sib = amr.patch[lv][PID+1]->sibling[19];
         if ( Sib >= 0 )  return Sib-5;
         else             return Sib;

      case 20: 
         Sib = amr.patch[lv][PID+2]->sibling[20];
         if ( Sib >= 0 )  return Sib-6;
         else             return Sib;

      case 21:
         Sib = amr.patch[lv][PID+4]->sibling[21];
         if ( Sib >= 0 )  return Sib-3;
         else             return Sib;

      case 22:
         Sib = amr.patch[lv][PID+3]->sibling[22];
         if ( Sib >= 0 )  return Sib-4;
         else             return Sib;

      case 23:
         Sib = amr.patch[lv][PID+6]->sibling[23];
         if ( Sib >= 0 )  return Sib-2;
         else             return Sib;

      case 24:
         Sib = amr.patch[lv][PID+5]->sibling[24];
         if ( Sib >= 0 )  return Sib-1;
         else             return Sib;

      case 25: 
         Sib = amr.patch[lv][PID+7]->sibling[25];
         if ( Sib >= 0 )  return Sib;
         else             return Sib;

      default:
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "Side", Side );
         exit(1);

   } // switch ( Side )

   return NULL_INT;

} // FUNCTION : Table_02



//-------------------------------------------------------------------------------------------------------
// Function    :  SetTargetSibling 
// Description :  Set the targeted sibling directions for preparing the ghost-zone data at the coarse-grid level
//
// Note        :  1. Work for the function "Prepare_PatchData"
//                2. TSib needs to be deallocated manually
//                3. Sibling directions recorded in TSib must be in ascending numerical order for filling the 
//                   non-periodic ghost-zone data in the function "InterpolateGhostZone"
//                   --> Therefore, this function CANNOT be applied in ""LB_RecordExchangeDataPatchID", in which 
//                       case "SetTargetSibling" and "SetReceiveSibling" must be declared consistently
// 
// Parameter   :  NTSib : Number of targeted sibling patches along different sibling directions 
//                TSib  : Targeted sibling indices along different sibling directions
//-------------------------------------------------------------------------------------------------------
void SetTargetSibling( int NTSib[], int *TSib[] )
{

   for (int t= 0; t< 6; t++)  NTSib[t] = 17;
   for (int t= 6; t<18; t++)  NTSib[t] = 11;
   for (int t=18; t<26; t++)  NTSib[t] =  7;

   for (int s=0; s<26; s++)   TSib[s] = new int [ NTSib[s] ];
   
   TSib[ 0][ 0] =  1;
   TSib[ 0][ 1] =  2;
   TSib[ 0][ 2] =  3;
   TSib[ 0][ 3] =  4;
   TSib[ 0][ 4] =  5;
   TSib[ 0][ 5] =  7;
   TSib[ 0][ 6] =  9;
   TSib[ 0][ 7] = 10;
   TSib[ 0][ 8] = 11;
   TSib[ 0][ 9] = 12;
   TSib[ 0][10] = 13;
   TSib[ 0][11] = 16;
   TSib[ 0][12] = 17;
   TSib[ 0][13] = 19;
   TSib[ 0][14] = 21;
   TSib[ 0][15] = 23;
   TSib[ 0][16] = 25;
   
   TSib[ 1][ 0] =  0;
   TSib[ 1][ 1] =  2;
   TSib[ 1][ 2] =  3;
   TSib[ 1][ 3] =  4;
   TSib[ 1][ 4] =  5;
   TSib[ 1][ 5] =  6;
   TSib[ 1][ 6] =  8;
   TSib[ 1][ 7] = 10;
   TSib[ 1][ 8] = 11;
   TSib[ 1][ 9] = 12;
   TSib[ 1][10] = 13;
   TSib[ 1][11] = 14;
   TSib[ 1][12] = 15;
   TSib[ 1][13] = 18;
   TSib[ 1][14] = 20;
   TSib[ 1][15] = 22;
   TSib[ 1][16] = 24;
   
   TSib[ 2][ 0] =  0;
   TSib[ 2][ 1] =  1;
   TSib[ 2][ 2] =  3;
   TSib[ 2][ 3] =  4;
   TSib[ 2][ 4] =  5;
   TSib[ 2][ 5] =  8;
   TSib[ 2][ 6] =  9;
   TSib[ 2][ 7] = 11;
   TSib[ 2][ 8] = 13;
   TSib[ 2][ 9] = 14;
   TSib[ 2][10] = 15;
   TSib[ 2][11] = 16;
   TSib[ 2][12] = 17;
   TSib[ 2][13] = 20;
   TSib[ 2][14] = 21;
   TSib[ 2][15] = 24;
   TSib[ 2][16] = 25;
   
   TSib[ 3][ 0] =  0;
   TSib[ 3][ 1] =  1;
   TSib[ 3][ 2] =  2;
   TSib[ 3][ 3] =  4;
   TSib[ 3][ 4] =  5;
   TSib[ 3][ 5] =  6;
   TSib[ 3][ 6] =  7;
   TSib[ 3][ 7] = 10;
   TSib[ 3][ 8] = 12;
   TSib[ 3][ 9] = 14;
   TSib[ 3][10] = 15;
   TSib[ 3][11] = 16;
   TSib[ 3][12] = 17;
   TSib[ 3][13] = 18;
   TSib[ 3][14] = 19;
   TSib[ 3][15] = 22;
   TSib[ 3][16] = 23;
   
   TSib[ 4][ 0] =  0;
   TSib[ 4][ 1] =  1;
   TSib[ 4][ 2] =  2;
   TSib[ 4][ 3] =  3;
   TSib[ 4][ 4] =  5;
   TSib[ 4][ 5] =  6;
   TSib[ 4][ 6] =  7;
   TSib[ 4][ 7] =  8;
   TSib[ 4][ 8] =  9;
   TSib[ 4][ 9] = 12;
   TSib[ 4][10] = 13;
   TSib[ 4][11] = 15;
   TSib[ 4][12] = 17;
   TSib[ 4][13] = 22;
   TSib[ 4][14] = 23;
   TSib[ 4][15] = 24;
   TSib[ 4][16] = 25;
   
   TSib[ 5][ 0] =  0;
   TSib[ 5][ 1] =  1;
   TSib[ 5][ 2] =  2;
   TSib[ 5][ 3] =  3;
   TSib[ 5][ 4] =  4;
   TSib[ 5][ 5] =  6;
   TSib[ 5][ 6] =  7;
   TSib[ 5][ 7] =  8;
   TSib[ 5][ 8] =  9;
   TSib[ 5][ 9] = 10;
   TSib[ 5][10] = 11;
   TSib[ 5][11] = 14;
   TSib[ 5][12] = 16;
   TSib[ 5][13] = 18;
   TSib[ 5][14] = 19;
   TSib[ 5][15] = 20;
   TSib[ 5][16] = 21;
   
   TSib[ 6][ 0] =  1;
   TSib[ 6][ 1] =  3;
   TSib[ 6][ 2] =  4;
   TSib[ 6][ 3] =  5;
   TSib[ 6][ 4] =  9;
   TSib[ 6][ 5] = 11;
   TSib[ 6][ 6] = 13;
   TSib[ 6][ 7] = 16;
   TSib[ 6][ 8] = 17;
   TSib[ 6][ 9] = 21;
   TSib[ 6][10] = 25;
   
   TSib[ 7][ 0] =  0;
   TSib[ 7][ 1] =  3;
   TSib[ 7][ 2] =  4;
   TSib[ 7][ 3] =  5;
   TSib[ 7][ 4] =  8;
   TSib[ 7][ 5] = 11;
   TSib[ 7][ 6] = 13;
   TSib[ 7][ 7] = 14;
   TSib[ 7][ 8] = 15;
   TSib[ 7][ 9] = 20;
   TSib[ 7][10] = 24;
   
   TSib[ 8][ 0] =  1;
   TSib[ 8][ 1] =  2;
   TSib[ 8][ 2] =  4;
   TSib[ 8][ 3] =  5;
   TSib[ 8][ 4] =  7;
   TSib[ 8][ 5] = 10;
   TSib[ 8][ 6] = 12;
   TSib[ 8][ 7] = 16;
   TSib[ 8][ 8] = 17;
   TSib[ 8][ 9] = 19;
   TSib[ 8][10] = 23;
   
   TSib[ 9][ 0] =  0;
   TSib[ 9][ 1] =  2;
   TSib[ 9][ 2] =  4;
   TSib[ 9][ 3] =  5;
   TSib[ 9][ 4] =  6;
   TSib[ 9][ 5] = 10;
   TSib[ 9][ 6] = 12;
   TSib[ 9][ 7] = 14;
   TSib[ 9][ 8] = 15;
   TSib[ 9][ 9] = 18;
   TSib[ 9][10] = 22;
   
   TSib[10][ 0] =  0;
   TSib[10][ 1] =  1;
   TSib[10][ 2] =  3;
   TSib[10][ 3] =  5;
   TSib[10][ 4] =  8;
   TSib[10][ 5] =  9;
   TSib[10][ 6] = 13;
   TSib[10][ 7] = 15;
   TSib[10][ 8] = 17;
   TSib[10][ 9] = 24;
   TSib[10][10] = 25;
   
   TSib[11][ 0] =  0;
   TSib[11][ 1] =  1;
   TSib[11][ 2] =  2;
   TSib[11][ 3] =  5;
   TSib[11][ 4] =  6;
   TSib[11][ 5] =  7;
   TSib[11][ 6] = 12;
   TSib[11][ 7] = 15;
   TSib[11][ 8] = 17;
   TSib[11][ 9] = 22;
   TSib[11][10] = 23;
   
   TSib[12][ 0] =  0;
   TSib[12][ 1] =  1;
   TSib[12][ 2] =  3;
   TSib[12][ 3] =  4;
   TSib[12][ 4] =  8;
   TSib[12][ 5] =  9;
   TSib[12][ 6] = 11;
   TSib[12][ 7] = 14;
   TSib[12][ 8] = 16;
   TSib[12][ 9] = 20;
   TSib[12][10] = 21;
   
   TSib[13][ 0] =  0;
   TSib[13][ 1] =  1;
   TSib[13][ 2] =  2;
   TSib[13][ 3] =  4;
   TSib[13][ 4] =  6;
   TSib[13][ 5] =  7;
   TSib[13][ 6] = 10;
   TSib[13][ 7] = 14;
   TSib[13][ 8] = 16;
   TSib[13][ 9] = 18;
   TSib[13][10] = 19;
   
   TSib[14][ 0] =  1;
   TSib[14][ 1] =  2;
   TSib[14][ 2] =  3;
   TSib[14][ 3] =  5;
   TSib[14][ 4] =  7;
   TSib[14][ 5] =  9;
   TSib[14][ 6] = 12;
   TSib[14][ 7] = 13;
   TSib[14][ 8] = 17;
   TSib[14][ 9] = 23;
   TSib[14][10] = 25;
   
   TSib[15][ 0] =  1;
   TSib[15][ 1] =  2;
   TSib[15][ 2] =  3;
   TSib[15][ 3] =  4;
   TSib[15][ 4] =  7;
   TSib[15][ 5] =  9;
   TSib[15][ 6] = 10;
   TSib[15][ 7] = 11;
   TSib[15][ 8] = 16;
   TSib[15][ 9] = 19;
   TSib[15][10] = 21;
   
   TSib[16][ 0] =  0;
   TSib[16][ 1] =  2;
   TSib[16][ 2] =  3;
   TSib[16][ 3] =  5;
   TSib[16][ 4] =  6;
   TSib[16][ 5] =  8;
   TSib[16][ 6] = 12;
   TSib[16][ 7] = 13;
   TSib[16][ 8] = 15;
   TSib[16][ 9] = 22;
   TSib[16][10] = 24;
   
   TSib[17][ 0] =  0;
   TSib[17][ 1] =  2;
   TSib[17][ 2] =  3;
   TSib[17][ 3] =  4;
   TSib[17][ 4] =  6;
   TSib[17][ 5] =  8;
   TSib[17][ 6] = 10;
   TSib[17][ 7] = 11;
   TSib[17][ 8] = 14;
   TSib[17][ 9] = 18;
   TSib[17][10] = 20;
   
   TSib[18][ 0] =  1;
   TSib[18][ 1] =  3;
   TSib[18][ 2] =  5;
   TSib[18][ 3] =  9;
   TSib[18][ 4] = 13;
   TSib[18][ 5] = 17;
   TSib[18][ 6] = 25;
   
   TSib[19][ 0] =  0;
   TSib[19][ 1] =  3;
   TSib[19][ 2] =  5;
   TSib[19][ 3] =  8;
   TSib[19][ 4] = 13;
   TSib[19][ 5] = 15;
   TSib[19][ 6] = 24;
   
   TSib[20][ 0] =  1;
   TSib[20][ 1] =  2;
   TSib[20][ 2] =  5;
   TSib[20][ 3] =  7;
   TSib[20][ 4] = 12;
   TSib[20][ 5] = 17;
   TSib[20][ 6] = 23;
   
   TSib[21][ 0] =  0;
   TSib[21][ 1] =  2;
   TSib[21][ 2] =  5;
   TSib[21][ 3] =  6;
   TSib[21][ 4] = 12;
   TSib[21][ 5] = 15;
   TSib[21][ 6] = 22;
   
   TSib[22][ 0] =  1;
   TSib[22][ 1] =  3;
   TSib[22][ 2] =  4;
   TSib[22][ 3] =  9;
   TSib[22][ 4] = 11;
   TSib[22][ 5] = 16;
   TSib[22][ 6] = 21;
   
   TSib[23][ 0] =  0;
   TSib[23][ 1] =  3;
   TSib[23][ 2] =  4;
   TSib[23][ 3] =  8;
   TSib[23][ 4] = 11;
   TSib[23][ 5] = 14;
   TSib[23][ 6] = 20;
   
   TSib[24][ 0] =  1;
   TSib[24][ 1] =  2;
   TSib[24][ 2] =  4;
   TSib[24][ 3] =  7;
   TSib[24][ 4] = 10;
   TSib[24][ 5] = 16;
   TSib[24][ 6] = 19;
   
   TSib[25][ 0] =  0;
   TSib[25][ 1] =  2;
   TSib[25][ 2] =  4;
   TSib[25][ 3] =  6;
   TSib[25][ 4] = 10;
   TSib[25][ 5] = 14;
   TSib[25][ 6] = 18;

} // FUNCTION : SetTargetSibling
