#include "GAMER.h"

static void WriteFile( FILE *File, const int lv, const int PID, const int i, const int j, const int k,
                       const int ii, const int jj, const int kk, const real (*DerField)[ CUBE(PS1) ] );
static void GetDerivedField( real (*Der_FluIn)[NCOMP_TOTAL][ CUBE(DER_NXT)            ],
                             real (*Der_Out  )             [ CUBE(PS1)                ],
                             real (*Der_MagFC)[NCOMP_MAG  ][ (DER_NXT+1)*SQR(DER_NXT) ],
                             real (*Der_MagCC)             [ CUBE(DER_NXT)            ],
                             const int lv, const int PID, const bool PrepFluIn );




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_DumpData_Part
// Description :  Output part of data in the ASCII form
//
// Note        :  1. Used for the runtime option "OPT__OUTPUT_PART"
//                2. For MHD, this function outputs the **cell-centered** magnetic field and energy
//
// Parameter   :  Part     : OUTPUT_XY   : xy plane
//                           OUTPUT_YZ   : yz plane
//                           OUTPUT_XZ   : xz plane
//                           OUTPUT_X    : x  line
//                           OUTPUT_Y    : y  line
//                           OUTPUT_Z    : z  line
//                           OUTPUT_DIAG : diagonal along (+1,+1,+1)
//                           OUTPUT_BOX  : entire box
//
//                BaseOnly : Only output the base-level data
//
//                x        : x coordinate
//                y        : y coordinate
//                z        : z coordinate
//
//                FileName : Name of the output file
//-------------------------------------------------------------------------------------------------------
void Output_DumpData_Part( const OptOutputPart_t Part, const bool BaseOnly, const double x, const double y,
                           const double z, const char *FileName )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d)           ...\n", __FUNCTION__, DumpID );


// check the input parameters
   if ( Part != OUTPUT_XY    &&  Part != OUTPUT_YZ  &&  Part != OUTPUT_XZ  &&
        Part != OUTPUT_X     &&  Part != OUTPUT_Y   &&  Part != OUTPUT_Z   &&
        Part != OUTPUT_DIAG  &&  Part != OUTPUT_BOX )
      Aux_Error( ERROR_INFO, "unsupported option \"Part = %d\" [0 ~ 8] !!\n", Part );

   if (  ( Part == OUTPUT_YZ  ||  Part == OUTPUT_Y  ||  Part == OUTPUT_Z )  &&
         ( x < 0.0  ||  x >= amr->BoxSize[0] )  )
      Aux_Error( ERROR_INFO, "incorrect x (out of range [0<=X<%lf]) !!\n", amr->BoxSize[0] );

   if (  ( Part == OUTPUT_XZ  ||  Part == OUTPUT_X  ||  Part == OUTPUT_Z )  &&
         ( y < 0.0  ||  y >= amr->BoxSize[1] )  )
      Aux_Error( ERROR_INFO, "incorrect y (out of range [0<=Y<%lf]) !!\n", amr->BoxSize[1] );

   if (  ( Part == OUTPUT_XY  ||  Part == OUTPUT_X  ||  Part == OUTPUT_Y )  &&
         ( z < 0.0  ||  z >= amr->BoxSize[2] )  )
      Aux_Error( ERROR_INFO, "incorrect z (out of range [0<=Z<%lf]) !!\n", amr->BoxSize[2] );

   if ( Part == OUTPUT_DIAG  &&  ( amr->BoxSize[0] != amr->BoxSize[1] || amr->BoxSize[0] != amr->BoxSize[2] )  )
      Aux_Error( ERROR_INFO, "simulation domain must be cubic for \"OUTPUT_DIAG\" !!\n" );


// check the synchronization
   for (int lv=1; lv<NLEVEL; lv++)
      if ( NPatchTotal[lv] != 0 )   Mis_CompareRealValue( Time[0], Time[lv], __FUNCTION__, true );


// check if the file already exists
   if ( MPI_Rank == 0 )
   {
      if ( Aux_CheckFileExist(FileName) )
      {
         Aux_Message( stderr, "WARNING : file \"%s\" already exists and will be overwritten !!\n", FileName );

         FILE *Temp = fopen( FileName, "w" );
         fclose( Temp );
      }
   }

   const double dh_min = amr->dh[NLEVEL-1];
   const int    NLv    = ( BaseOnly ) ? 1 : NLEVEL;

   int     ii, jj, kk, scale;
   double  dh, xx, yy, zz;    // xx,yy,zz => physical coordinates of cell left edge
   int    *Corner  = NULL;    // patch corner in scale
   double *EdgeL   = NULL;    // patch corner in physical coord.
   double *EdgeR   = NULL;
   bool    Check_x = false;
   bool    Check_y = false;
   bool    Check_z = false;

   switch ( Part )
   {
      case OUTPUT_XY   :                                      Check_z = true;   break;
      case OUTPUT_YZ   :  Check_x = true;                                       break;
      case OUTPUT_XZ   :                    Check_y = true;                     break;
      case OUTPUT_X    :                    Check_y = true;   Check_z = true;   break;
      case OUTPUT_Y    :  Check_x = true;                     Check_z = true;   break;
      case OUTPUT_Z    :  Check_x = true;   Check_y = true;                     break;

      case OUTPUT_DIAG :
      case OUTPUT_BOX  :
      case OUTPUT_PART_NONE : break; // do nothing
   }


// for the derived fields
   const int Der_NP = 8;
   bool Der_PrepFluIn;

   real (*Der_FluIn)[NCOMP_TOTAL][ CUBE(DER_NXT)            ] = new real [Der_NP][NCOMP_TOTAL ][ CUBE(DER_NXT)            ];
   real (*Der_Out  )             [ CUBE(PS1)                ] = new real         [DER_NOUT_MAX][ CUBE(PS1)                ];
#  ifdef MHD
   real (*Der_MagFC)[NCOMP_MAG  ][ (DER_NXT+1)*SQR(DER_NXT) ] = new real [Der_NP][NCOMP_MAG   ][ (DER_NXT+1)*SQR(DER_NXT) ];
   real (*Der_MagCC)             [ CUBE(DER_NXT)            ] = new real         [NCOMP_MAG   ][ CUBE(DER_NXT)            ];
#  else
   real (*Der_MagFC)[NCOMP_MAG  ][ (DER_NXT+1)*SQR(DER_NXT) ] = NULL;
   real (*Der_MagCC)             [ CUBE(DER_NXT)            ] = NULL;
#  endif


   for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)
   {
      if ( MPI_Rank == TargetMPIRank )
      {
         FILE *File = fopen( FileName, "a" );

//       output header
         if ( TargetMPIRank == 0 )
         {
            fprintf( File, "#%10s %10s %10s %20s %20s %20s", "i", "j", "k", "x", "y", "z" );

            for (int v=0; v<NCOMP_TOTAL; v++)
            fprintf( File, " %*s", StrLen_Flt, FieldLabel[v] );

#           ifdef MHD
            for (int v=0; v<NCOMP_MAG; v++)
            fprintf( File, " %*s", StrLen_Flt, MagLabel[v] );

            fprintf( File, " %*s", StrLen_Flt, "MagEngy" );
#           endif

#           ifdef GRAVITY
            if ( OPT__OUTPUT_POT )     fprintf( File, " %*s", StrLen_Flt, PotLabel );
#           endif

//          derived fields
#           if ( MODEL == HYDRO )
            if ( OPT__OUTPUT_PRES )    fprintf( File, " %*s", StrLen_Flt, "Pressure" );
            if ( OPT__OUTPUT_TEMP )    fprintf( File, " %*s", StrLen_Flt, "Temperature" );
            if ( OPT__OUTPUT_ENTR )    fprintf( File, " %*s", StrLen_Flt, "Entropy" );
            if ( OPT__OUTPUT_CS )      fprintf( File, " %*s", StrLen_Flt, "Sound speed" );
            if ( OPT__OUTPUT_DIVVEL )  fprintf( File, " %*s", StrLen_Flt, "Div(Vel)" );
            if ( OPT__OUTPUT_MACH   )  fprintf( File, " %*s", StrLen_Flt, "Mach" );
#           endif
#           ifdef MHD
            if ( OPT__OUTPUT_DIVMAG )  fprintf( File, " %*s", StrLen_Flt, "Div(Mag)" );
#           endif
#           ifdef SRHD
            if ( OPT__OUTPUT_LORENTZ ) fprintf( File, " %*s", StrLen_Flt, "Lorentz" );
            if ( OPT__OUTPUT_3VELOCITY )
            {
                                       fprintf( File, " %*s", StrLen_Flt, "Velocity X" );
                                       fprintf( File, " %*s", StrLen_Flt, "Velocity Y" );
                                       fprintf( File, " %*s", StrLen_Flt, "Velocity Z" );
            }
            if ( OPT__OUTPUT_ENTHALPY )
                                       fprintf( File, " %*s", StrLen_Flt, "Reduced enthalpy" );
#           endif
            if ( OPT__OUTPUT_USER_FIELD ) {
               for (int v=0; v<UserDerField_Num; v++)
                                       fprintf( File, " %*s", StrLen_Flt, UserDerField_Label[v] );
            }

            fprintf( File, "\n" );
         } // if ( TargetMPIRank == 0 )


//       output data
         for (int lv=0; lv<NLv; lv++)
         {
            dh    = amr->dh   [lv];
            scale = amr->scale[lv];

            for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
            {
               if ( PID % 8 == 0 )  Der_PrepFluIn = true;

//             output the patch data only if it has no son (if the option "BaseOnly" is turned off)
               if ( amr->patch[0][lv][PID]->son == -1  ||  BaseOnly )
               {
                  Corner = amr->patch[0][lv][PID]->corner;
                  EdgeL  = amr->patch[0][lv][PID]->EdgeL;
                  EdgeR  = amr->patch[0][lv][PID]->EdgeR;

                  if ( Part == OUTPUT_DIAG ) // (+1,+1,+1) diagonal
                  {
//                   check whether the patch corner is along the diagonal
                     if ( Corner[0] == Corner[1]  &&  Corner[0] == Corner[2] )
                     {
//                      compute the derived fields
                        GetDerivedField( Der_FluIn, Der_Out, Der_MagFC, Der_MagCC, lv, PID, Der_PrepFluIn );
                        Der_PrepFluIn = false;

//                      write data
                        for (int k=0; k<PS1; k++)
                        {
                           kk = Corner[2] + k*scale;

                           WriteFile( File, lv, PID, k, k, k, kk, kk, kk, Der_Out );
                        }
                     }
                  } // if ( Part == OUTPUT_DIAG )


                  else // x/y/z lines || xy/yz/xz slices
                  {
//                   check whether the patch corner is within the target range
                     if (  !Check_x  ||  ( EdgeL[0]<=x && EdgeR[0]>x )  )
                     if (  !Check_y  ||  ( EdgeL[1]<=y && EdgeR[1]>y )  )
                     if (  !Check_z  ||  ( EdgeL[2]<=z && EdgeR[2]>z )  )
                     {
//                      compute the derived fields
                        GetDerivedField( Der_FluIn, Der_Out, Der_MagFC, Der_MagCC, lv, PID, Der_PrepFluIn );
                        Der_PrepFluIn = false;

//                      write data
//                      --> check whether the cell is within the target range
                        for (int k=0; k<PS1; k++)  {  kk = Corner[2] + k*scale;  zz = kk*dh_min;
                                                      if ( Check_z && ( zz>z || zz+dh<=z ) )    continue;

                        for (int j=0; j<PS1; j++)  {  jj = Corner[1] + j*scale;  yy = jj*dh_min;
                                                      if ( Check_y && ( yy>y || yy+dh<=y ) )    continue;

                        for (int i=0; i<PS1; i++)  {  ii = Corner[0] + i*scale;  xx = ii*dh_min;
                                                      if ( Check_x && ( xx>x || xx+dh<=x ) )    continue;

                           WriteFile( File, lv, PID, i, j, k, ii, jj, kk, Der_Out );

                        }}}
                     } // if patch corner is within the target range

                  } // if ( Part == OUTPUT_DIAG ) ... else ...
               } // if ( amr->patch[0][lv][PID]->son == -1 )
            } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         } // for (int lv=0; lv<NLv; lv++)

         fclose( File );

      } // if ( MPI_Rank == TargetMPIRank )

      MPI_Barrier( MPI_COMM_WORLD );

   } // for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)


   delete [] Der_FluIn;
   delete [] Der_Out;
#  ifdef MHD
   delete [] Der_MagFC;
   delete [] Der_MagCC;
#  endif

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d)           ... done\n", __FUNCTION__, DumpID );

} // FUNCTION : Output_DumpData_Part



//-------------------------------------------------------------------------------------------------------
// Function    :  WriteFile
// Description :  Output data to file
//
// Parameter   :  File     : File pointer
//                lv       : Target refinement level
//                PID      : Patch ID
//                i/j/k    : Cell indices within the patch
//                ii/jj/kk : Cell scale indices in the simulation domain
//                DerField : Array storing the derived fields
//-------------------------------------------------------------------------------------------------------
void WriteFile( FILE *File, const int lv, const int PID, const int i, const int j, const int k,
                const int ii, const int jj, const int kk, const real (*DerField)[ CUBE(PS1) ] )
{

   const double dh_min  = amr->dh[TOP_LEVEL];
   const double scale_2 = 0.5*amr->scale[lv];
   real u[NCOMP_TOTAL];

   for (int v=0; v<NCOMP_TOTAL; v++)   u[v] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i];

// cell indices and coordinates
   fprintf( File, " %10d %10d %10d %20.14e %20.14e %20.14e",
            ii, jj, kk, (ii+scale_2)*dh_min, (jj+scale_2)*dh_min, (kk+scale_2)*dh_min );

// output all variables in the fluid array
   for (int v=0; v<NCOMP_TOTAL; v++)   fprintf( File, BlankPlusFormat_Flt, u[v] );

// magnetic field
#  if ( MODEL == HYDRO )
#  ifdef MHD
   const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, amr->MagSg[lv] );
   real B[3];
   MHD_GetCellCenteredBFieldInPatch( B, lv, PID, i, j, k, amr->MagSg[lv] );
   fprintf( File, BlankPlusFormat_Flt, B[MAGX] );
   fprintf( File, BlankPlusFormat_Flt, B[MAGY] );
   fprintf( File, BlankPlusFormat_Flt, B[MAGZ] );
   fprintf( File, BlankPlusFormat_Flt, Emag    );
#  else
   const real Emag = NULL_REAL;
#  endif
#  endif // # if ( MODEL == HYDRO )

// output potential
#  ifdef GRAVITY
   if ( OPT__OUTPUT_POT )
      fprintf( File, BlankPlusFormat_Flt, amr->patch[ amr->PotSg[lv] ][lv][PID]->pot[k][j][i] );
#  endif

// output derived fields
   const int Der_CellIdx = IDX321( i, j, k, PS1, PS1 );
   int Der_FieldIdx = 0;

#  if ( MODEL == HYDRO )
   const bool CheckMinPres_No = false;
   const bool CheckMinTemp_No = false;
   const bool CheckMinEntr_No = false;
   real Pres=-1.0, Temp=-1.0, Entr=-1.0, Cs=-1.0;  // initial values must be negative

// no need to increase Der_FieldIdx for fields not using DerField[]
   if ( OPT__OUTPUT_PRES ) {
      Pres = Hydro_Con2Pres( u[DENS], u[MOMX], u[MOMY], u[MOMZ], u[ENGY], u+NCOMP_FLUID,
                             CheckMinPres_No, NULL_REAL, Emag, EoS_DensEint2Pres_CPUPtr,
                             EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                             EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL );
      fprintf( File, BlankPlusFormat_Flt, Pres );
   }

   if ( OPT__OUTPUT_TEMP ) {
      Temp = Hydro_Con2Temp( u[DENS], u[MOMX], u[MOMY], u[MOMZ], u[ENGY], u+NCOMP_FLUID,
                             CheckMinTemp_No, NULL_REAL, Emag, EoS_DensEint2Temp_CPUPtr,
                             EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                             EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
      fprintf( File, BlankPlusFormat_Flt, Temp );
   }

#  ifndef SRHD
   if ( OPT__OUTPUT_ENTR ) {
      Entr = Hydro_Con2Entr( u[DENS], u[MOMX], u[MOMY], u[MOMZ], u[ENGY], u+NCOMP_FLUID,
                             CheckMinEntr_No, NULL_REAL, Emag, EoS_DensEint2Entr_CPUPtr,
                             EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
      fprintf( File, BlankPlusFormat_Flt, Entr );
   }
#  endif

#  ifdef SRHD
   real Prim[NCOMP_TOTAL], LorentzFactor=-1.0, HTilde=-1.0;
   if ( OPT__OUTPUT_CS || OPT__OUTPUT_LORENTZ || OPT__OUTPUT_3VELOCITY )
      Hydro_Con2Pri( u, Prim, (real)-HUGE_NUMBER, NULL_BOOL, NULL_INT, NULL,
                     NULL_BOOL, NULL_REAL, EoS_DensEint2Pres_CPUPtr,
                     EoS_DensPres2Eint_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                     EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL, &LorentzFactor );
#  endif

   if ( OPT__OUTPUT_CS ) {
#     ifdef SRHD
      Cs = SQRT(  EoS_DensPres2CSqr_CPUPtr( Prim[0], Prim[4], NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table )  );
#     else
//    compute pressure if it is not done yet
      if ( Pres < 0.0 )
      Pres = Hydro_Con2Pres( u[DENS], u[MOMX], u[MOMY], u[MOMZ], u[ENGY], u+NCOMP_FLUID,
                             CheckMinPres_No, NULL_REAL, Emag, EoS_DensEint2Pres_CPUPtr,
                             EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                             EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL );
      Cs   = SQRT(  EoS_DensPres2CSqr_CPUPtr( u[DENS], Pres, u+NCOMP_FLUID, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                              h_EoS_Table )  );
#     endif
      fprintf( File, BlankPlusFormat_Flt, Cs );
   }

   if ( OPT__OUTPUT_DIVVEL )
      fprintf( File, BlankPlusFormat_Flt, DerField[ Der_FieldIdx ++ ][Der_CellIdx] );

   if ( OPT__OUTPUT_MACH )
      fprintf( File, BlankPlusFormat_Flt, DerField[ Der_FieldIdx ++ ][Der_CellIdx] );
#  endif // #if ( MODEL == HYDRO )

#  ifdef MHD
   if ( OPT__OUTPUT_DIVMAG ) {
      const real DivB = MHD_GetCellCenteredDivBInPatch( lv, PID, i, j, k, amr->MagSg[lv] );
      fprintf( File, BlankPlusFormat_Flt, DivB );
   }
#  endif

#  ifdef SRHD
   if ( OPT__OUTPUT_LORENTZ )
      fprintf( File, BlankPlusFormat_Flt, LorentzFactor );

   if ( OPT__OUTPUT_3VELOCITY )
   {
      fprintf( File, BlankPlusFormat_Flt, Prim[1] / LorentzFactor );
      fprintf( File, BlankPlusFormat_Flt, Prim[2] / LorentzFactor );
      fprintf( File, BlankPlusFormat_Flt, Prim[3] / LorentzFactor );
   }

   if ( OPT__OUTPUT_ENTHALPY )
   {
      HTilde = Hydro_Con2HTilde( u, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                 EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
      fprintf( File, BlankPlusFormat_Flt, HTilde );
   }
#  endif

   if ( OPT__OUTPUT_USER_FIELD ) {
      for (int v=0; v<UserDerField_Num; v++)
      fprintf( File, BlankPlusFormat_Flt, DerField[ Der_FieldIdx ++ ][Der_CellIdx] );
   }

   fprintf( File, "\n" );

} // FUNCTION : WriteFile



//-------------------------------------------------------------------------------------------------------
// Function    :  GetDerivedField
// Description :  Compute the derived fields
//
// Note        :  1. FluIn[] and MagFC[] will be filled in only if "PrepFluIn == true"
//                2. Called by Output_DumpData_Part()
//
// Parameter   :  FluIn     : Array to store the input fluid data for the derived field functions
//                Out       : Array to store the output derived fields
//                MagFC     : Array to store the temporary face-centered B field
//                MagCC     : Array to store the input cell-centered B field for the derived field functions
//                lv        : Target refinement level
//                PID       : Target patch ID
//                PrepFluIn : Whether to fill in FluIn[] and MagCC[]
//                            --> To prepare patches within the same patch group just once
//
// Return      :  FluIn[], Out[], MagFC[] (MagCC[] is not useful outside this function)
//-------------------------------------------------------------------------------------------------------
void GetDerivedField( real (*FluIn)[NCOMP_TOTAL][ CUBE(DER_NXT)            ],
                      real (*Out  )             [ CUBE(PS1)                ],
                      real (*MagFC)[NCOMP_MAG  ][ (DER_NXT+1)*SQR(DER_NXT) ],
                      real (*MagCC)             [ CUBE(DER_NXT)            ],
                      const int lv, const int PID, const bool PrepFluIn )
{

   const double dh      = amr->dh[lv];
   const int    LocalID = PID % 8;


// prepare the input arrays
   if ( PrepFluIn )
   {
      const int  PID0                = PID - LocalID;
      const int  NP                  = 8;
      const bool IntPhase_No         = false;
      const real MinDens_No          = -1.0;
      const real MinPres_No          = -1.0;
      const real MinTemp_No          = -1.0;
      const real MinEntr_No          = -1.0;
      const bool DE_Consistency_No   = false;
#     ifndef MHD
      const int  OPT__MAG_INT_SCHEME = INT_NONE;
#     endif

//    always prepare all fields
      Prepare_PatchData( lv, Time[lv], FluIn[0][0], MagFC[0][0], DER_GHOST_SIZE, 1, &PID0,
                         _TOTAL, _MAG, OPT__FLU_INT_SCHEME, OPT__MAG_INT_SCHEME, UNIT_PATCH, NSIDE_26,
                         IntPhase_No, OPT__BC_FLU, BC_POT_NONE, MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );
   } // if ( PrepFluIn )


// convert B field from face-centered to cell-centered
#  ifdef MHD
   for (int k=0; k<DER_NXT; k++)
   for (int j=0; j<DER_NXT; j++)
   for (int i=0; i<DER_NXT; i++)
   {
      const int IdxCC = IDX321( i, j, k, DER_NXT, DER_NXT );
      real B_CC[NCOMP_MAG];

      MHD_GetCellCenteredBField( B_CC, MagFC[LocalID][MAGX], MagFC[LocalID][MAGY], MagFC[LocalID][MAGZ],
                                 DER_NXT, DER_NXT, DER_NXT, i, j, k );

      MagCC[MAGX][IdxCC] = B_CC[MAGX];
      MagCC[MAGY][IdxCC] = B_CC[MAGY];
      MagCC[MAGZ][IdxCC] = B_CC[MAGZ];
   }
#  endif // #ifdef MHD


// calculate the derived fields
   int OutFieldIdx = 0;

#  if ( MODEL == HYDRO )
   if ( OPT__OUTPUT_DIVVEL )
   {
      const int NFieldOut = 1;

      if ( OutFieldIdx + NFieldOut > DER_NOUT_MAX )
         Aux_Error( ERROR_INFO, "OutFieldIdx (%d) + NFieldOut (%d) > DER_NOUT_MAX (%d) !!\n",
                    OutFieldIdx, NFieldOut, DER_NOUT_MAX );

      Flu_DerivedField_DivVel( Out[OutFieldIdx], FluIn[LocalID][0], NULL,
                               NFieldOut, DER_NXT, DER_NXT, DER_NXT, DER_GHOST_SIZE, dh );

      OutFieldIdx += NFieldOut;
   }

   if ( OPT__OUTPUT_MACH )
   {
      const int NFieldOut = 1;

      if ( OutFieldIdx + NFieldOut > DER_NOUT_MAX )
         Aux_Error( ERROR_INFO, "OutFieldIdx (%d) + NFieldOut (%d) > DER_NOUT_MAX (%d) !!\n",
                    OutFieldIdx, NFieldOut, DER_NOUT_MAX );

      Flu_DerivedField_Mach( Out[OutFieldIdx], FluIn[LocalID][0], MagCC[0],
                             NFieldOut, DER_NXT, DER_NXT, DER_NXT, DER_GHOST_SIZE, dh );

      OutFieldIdx += NFieldOut;
   }
#  endif // #if ( MODEL == HYDRO )

   if ( OPT__OUTPUT_USER_FIELD )
   {
      const int NFieldOut = UserDerField_Num;

      if ( OutFieldIdx + NFieldOut > DER_NOUT_MAX )
         Aux_Error( ERROR_INFO, "OutFieldIdx (%d) + NFieldOut (%d) > DER_NOUT_MAX (%d) !!\n",
                    OutFieldIdx, NFieldOut, DER_NOUT_MAX );

      Flu_DerivedField_User_Ptr( Out[OutFieldIdx], FluIn[LocalID][0], MagCC[0],
                                 NFieldOut, DER_NXT, DER_NXT, DER_NXT, DER_GHOST_SIZE, dh );

      OutFieldIdx += NFieldOut;
   }

} // FUNCTION : GetDerivedField
