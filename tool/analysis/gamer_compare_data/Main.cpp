#include "CompareData.h"



// no need to share with other files
static AMR_t        amr1, amr2;
static char        *FileName_In1=NULL, *FileName_In2=NULL, *FileName_Out=NULL;
static bool         UseCorner=false, MaxErrOnly=false;
static int          NField1=-1, NField2=-1, NMag1=-1, NMag2=-1, NParAttFlt1=-1, NParAttFlt2=-1, NParAttInt1=-1, NParAttInt2=-1, Format1=-1, Format2=-1;
static long         NPar1=-1, NPar2=-1;
static double       TolErr=__FLT_MIN__;
static real_par   **ParFltData1=NULL, **ParFltData2=NULL;
static long_par   **ParIntData1=NULL, **ParIntData2=NULL;
static char  (*FieldLabel1)[MAX_STRING]=NULL, (*MagLabel1)[MAX_STRING]=NULL, (*ParAttFltLabel1)[MAX_STRING]=NULL, (*ParAttIntLabel1)[MAX_STRING]=NULL;
static char  (*FieldLabel2)[MAX_STRING]=NULL, (*MagLabel2)[MAX_STRING]=NULL, (*ParAttFltLabel2)[MAX_STRING]=NULL, (*ParAttIntLabel2)[MAX_STRING]=NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  ReadOption
// Description :  Read options from the command line
//-------------------------------------------------------------------------------------------------------
void ReadOption( int argc, char **argv )
{

   int c;

   while( (c = getopt(argc, argv, "hcmi:j:o:e:")) != -1 )
      switch(c)
      {
         case 'i': FileName_In1  = optarg;
                   break;
         case 'j': FileName_In2  = optarg;
                   break;
         case 'o': FileName_Out  = optarg;
                   break;
         case 'e': TolErr        = atof(optarg);
                   break;
         case 'c': UseCorner     = true;
                   break;
         case 'm': MaxErrOnly    = true;
                   break;
         case 'h':
         case '?': cerr << endl << "usage: " << argv[0]
                        << " [-h (for help)] [-i Input FileName1] [-j Input FileName2] [-o Output FileName]"
                        << endl << "                          "
                        << " [-e Tolerant Error] [-c (compare patches with the same corner coordinates) [off]]"
                        << endl << "                          "
                        << " [-m Only output the maximum error of each field [off]]"
                        << endl
                        << endl << endl;
                   exit( EXIT_SUCCESS );
      }

} // FUNCTION : ReadOption



//-------------------------------------------------------------------------------------------------------
// Function    :  CheckParameter
// Description :  Verify the input parameters
//-------------------------------------------------------------------------------------------------------
void CheckParameter()
{

   Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   if ( FileName_In1 == NULL )
      Aux_Error( ERROR_INFO, "please provide the name of the input file one (-i FileName1) !!\n" );

   if ( FileName_In2 == NULL )
      Aux_Error( ERROR_INFO, "please provide the name of the input file two (-j FileName2) !!\n" );

   if ( FileName_Out == NULL )
      Aux_Error( ERROR_INFO, "please provide the name of the output file (-o FileName) !!\n" );

   if ( TolErr == __FLT_MIN__ )
      Aux_Message( stderr, "WARNING : please provide the tolerant error (-e Tolerant Error) !!\n" );

   if ( UseCorner == false )
      Aux_Message( stderr, "WARNING : please make sure that the order of patches stored in the input files are the same.\n"
                           "          Otherwise, you should turn on the option \"-c\" !!\n" );


   Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : CheckParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  CompareGridData
// Description :  Compare the grid data between two files
//
// Return      :  EXIT_SUCCESS : no error
//                EXIT_FAILURE : with errors
//-------------------------------------------------------------------------------------------------------
int CompareGridData()
{

   Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// check
// verify that the total numbers of patches on each level are consistent
   int NPatch1=-1, NPatch2=-1;
   for (int lv=0; lv<NLEVEL; lv++)
   {
      NPatch1 = amr1.num[lv];
      NPatch2 = amr2.num[lv];

//    exclude non-leaf patches in HDF5 when comparing it with a C-binary file
      if ( lv != NLEVEL-1 )
      {
         if      ( Format1 == 1  &&  Format2 == 2 )   NPatch2 -= amr2.num[lv+1]/8;
         else if ( Format1 == 2  &&  Format2 == 1 )   NPatch1 -= amr1.num[lv+1]/8;
      }

      if ( NPatch1 != NPatch2 )
         Aux_Error( ERROR_INFO, "inconsistent numbers of patches on level %d (%d vs %d) !!\n",
                    lv, NPatch1, NPatch2 );
   }

// verify that the simulation domains are consistent
   for (int d=0; d<3; d++)
   {
      if ( amr1.nx0_tot[d] != amr2.nx0_tot[d] )
         Aux_Error( ERROR_INFO, "inconsistent root level grid size ([%d]: %d vs %d) !!\n",
                    d, amr1.nx0_tot[d], amr2.nx0_tot[d] );
   }

// verify that the numbers of fields are consistent
   if ( NField1 != NField2 )
      Aux_Error( ERROR_INFO, "inconsistent numbers of fields (%d vs %d) !!\n", NField1, NField2 );

   if ( NMag1 != NMag2 )
      Aux_Error( ERROR_INFO, "inconsistent numbers of magnetic components (%d vs %d) !!\n", NMag1, NMag2 );

// verify that the field labels are consistent
   if ( Format1 == 2  &&  Format2 == 2 )
   {
      for (int v=0; v<NField1; v++)
      {
         if ( FieldLabel1 == NULL  ||  FieldLabel2 == NULL )
            Aux_Error( ERROR_INFO, "FieldLabel[%d] == NULL !!\n", v );

         if (  strcmp( FieldLabel1[v], FieldLabel2[v] )  )
            Aux_Error( ERROR_INFO, "inconsistent field labels ([%d]: \"%s\" vs \"%s\") !!\n",
                       v, FieldLabel1[v], FieldLabel2[v] );
      }

      for (int v=0; v<NMag1; v++)
      {
         if ( MagLabel1 == NULL  ||  MagLabel2 == NULL )
            Aux_Error( ERROR_INFO, "MagLabel[%d] == NULL !!\n", v );

         if (  strcmp( MagLabel1[v], MagLabel2[v] )  )
            Aux_Error( ERROR_INFO, "inconsistent magnetic field labels ([%d]: \"%s\" vs \"%s\") !!\n",
                       v, MagLabel1[v], MagLabel2[v] );
      }
   } // if ( Format1 == 2  &&  Format2 == 2 )


// enable UseCorner when comparing files with different formats
   if ( UseCorner == false  &&  Format1 != Format2 )
   {
      UseCorner = true;

      Aux_Message( stderr, "WARNING : Format1 (%d) != Format2 (%d) !!\n"
                           "          --> The option \"-c\" is turned on automatically\n",
                   Format1, Format2 );
   }


// only compare non-leaf patches for C-binary files
   const bool LeafOnly = ( Format1 == 1 || Format2 == 1 ) ? true : false;


// compare data
   const int MagIdx0 = 1000;

   double Data1, Data2, AbsErr, RelErr;
   int    PID1, PID2;
   int   *Cr1, *Cr2;
   bool   ErrorDetected = false;

   double MaxErrCC_Data[NField1][4], MaxErrFC_Data[NMag1][4];
   int    MaxErrCC_Info[NField1][7], MaxErrFC_Info[NMag1][7];

   for (int v=0; v<NField1; v++)    MaxErrCC_Data[v][3] = -__FLT_MIN__;
   for (int v=0; v<NMag1;   v++)    MaxErrFC_Data[v][3] = -__FLT_MIN__;

   FILE *File = fopen( FileName_Out, "w" );

   if ( Format1 == 2  &&  Format2 == 2 )
   {
      fprintf( File, "# Field list:\n" );
      for (int v=0; v<NField1; v++)    fprintf( File, "# [%4d]: %s\n", v,         FieldLabel1[v] );
      for (int v=0; v<NMag1; v++)      fprintf( File, "# [%4d]: %s\n", v+MagIdx0, MagLabel1  [v] );
      fprintf( File, "\n" );
   }

   fprintf( File, "#%5s%8s%8s  (%3s,%3s,%3s )%6s%16s%16s%16s%16s\n",
                  "Level", "PID1", "PID2", "i", "j", "k", "Field", "Data1", "Data2", "AbsErr", "RelErr" );

   for (int lv=0; lv<NLEVEL; lv++)
   {
      Aux_Message( stdout, "  Comparing level %d ... ", lv );

      for (PID1=0; PID1<amr1.num[lv]; PID1++)
      {
//       only compare leaf patches when enabling LeafOnly
         if ( !LeafOnly  ||  amr1.patch[lv][PID1]->son == -1 )
         {
//          set the targeted patch ID in the second input
            if ( UseCorner )
            {
               Cr1 = amr1.patch[lv][PID1]->corner;

               for (PID2=0; PID2<amr2.num[lv]; PID2++)
               {
                  Cr2 = amr2.patch[lv][PID2]->corner;

                  if ( Cr1[0] == Cr2[0]  &&  Cr1[1] == Cr2[1]  &&  Cr1[2] == Cr2[2] )  break;
               }
            }

            else
            {
               PID2 = PID1;

//             assert the two patches have the same corner coordinates
               Cr1 = amr1.patch[lv][PID1]->corner;
               Cr2 = amr2.patch[lv][PID2]->corner;

               for (int d=0; d<3; d++)
                  if ( Cr1[d] != Cr2[d] )
                     Aux_Error( ERROR_INFO, "patch %d in the two snapshots do not match: d %d, corner %d != %d !!\n"
                                            "        --> Try enabling the option \"-c\"\n",
                                PID1, PID2, d, Cr1[d], Cr2[d] );
            }


//          a. cell-centered fields
            for (int k=0; k<PS1; k++)
            for (int j=0; j<PS1; j++)
            for (int i=0; i<PS1; i++)
            for (int v=0; v<NField1; v++)
            {
               Data1  = amr1.patch[lv][PID1]->field[v][k][j][i];
               Data2  = amr2.patch[lv][PID2]->field[v][k][j][i];
               AbsErr = Data1 - Data2;
               RelErr = AbsErr / ( 0.5*(Data1+Data2) );

               if ( Data1 == 0.0  &&  Data2 == 0.0 )  continue;

               if ( fabs(RelErr) >= TolErr  ||  !isfinite(RelErr) )
               {
                  ErrorDetected = true;

                  if ( MaxErrOnly ) {
                     if ( fabs(RelErr) > fabs(MaxErrCC_Data[v][3]) ) {
                        MaxErrCC_Info[v][0] = lv;
                        MaxErrCC_Info[v][1] = PID1;
                        MaxErrCC_Info[v][2] = PID2;
                        MaxErrCC_Info[v][3] = i;
                        MaxErrCC_Info[v][4] = j;
                        MaxErrCC_Info[v][5] = k;
                        MaxErrCC_Info[v][6] = v;

                        MaxErrCC_Data[v][0] = Data1;
                        MaxErrCC_Data[v][1] = Data2;
                        MaxErrCC_Data[v][2] = AbsErr;
                        MaxErrCC_Data[v][3] = RelErr;
                     }
                  }

                  else
                     fprintf( File, "%6d%8d%8d  (%3d,%3d,%3d )%6d%16.7e%16.7e%16.7e%16.7e\n",
                              lv, PID1, PID2, i, j, k, v, Data1, Data2, AbsErr, RelErr );
               }
            } // i,j,k,v


//          b. face-centered B field
            for (int t=0; t<PS1P1*SQR(PS1); t++)
            for (int v=0; v<NMag1; v++)
            {
               Data1  = amr1.patch[lv][PID1]->mag[v][t];
               Data2  = amr2.patch[lv][PID2]->mag[v][t];
               AbsErr = Data1 - Data2;
               RelErr = AbsErr / ( 0.5*(Data1+Data2) );

               if ( Data1 == 0.0  &&  Data2 == 0.0 )  continue;

               if ( fabs(RelErr) >= TolErr  ||  !isfinite(RelErr) )
               {
                  ErrorDetected = true;

                  int size_i, size_j, size_k, size_ij, i, j, k;

                  switch ( v )
                  {
                     case 0:  size_i = PS1P1;  size_j = PS1;    size_k = PS1;    break;
                     case 1:  size_i = PS1;    size_j = PS1P1;  size_k = PS1;    break;
                     case 2:  size_i = PS1;    size_j = PS1;    size_k = PS1P1;  break;
                  }

                  size_ij = size_i*size_j;
                  i       = t % size_i;
                  j       = t % size_ij / size_i;
                  k       = t / size_ij;

                  if ( MaxErrOnly ) {
                     if ( fabs(RelErr) > fabs(MaxErrFC_Data[v][3]) ) {
                        MaxErrFC_Info[v][0] = lv;
                        MaxErrFC_Info[v][1] = PID1;
                        MaxErrFC_Info[v][2] = PID2;
                        MaxErrFC_Info[v][3] = i;
                        MaxErrFC_Info[v][4] = j;
                        MaxErrFC_Info[v][5] = k;
                        MaxErrFC_Info[v][6] = v+MagIdx0;

                        MaxErrFC_Data[v][0] = Data1;
                        MaxErrFC_Data[v][1] = Data2;
                        MaxErrFC_Data[v][2] = AbsErr;
                        MaxErrFC_Data[v][3] = RelErr;
                     }
                  }

                  else
                     fprintf( File, "%6d%8d%8d  (%3d,%3d,%3d )%6d%16.7e%16.7e%16.7e%16.7e\n",
                              lv, PID1, PID2, i, j, k, v+MagIdx0, Data1, Data2, AbsErr, RelErr );
               }
            } // for v, t

            amr1.patch[lv][PID1]->check = true;
            amr2.patch[lv][PID2]->check = true;

         } // if ( amr1[lv][PID1]->son == -1 )

      } // for (PID1=0; PID1<amr1.num[lv]; PID1++)

      Aux_Message( stdout, "done\n" );

   } // for (int lv=0; lv<NLEVEL; lv++)


// output the maximum errors
   if ( MaxErrOnly )
   {
      for (int v=0; v<NField1; v++)
         if ( MaxErrCC_Data[v][3] != -__FLT_MIN__ )
            fprintf( File, "%6d%8d%8d  (%3d,%3d,%3d )%6d%16.7e%16.7e%16.7e%16.7e\n",
                     MaxErrCC_Info[v][0], MaxErrCC_Info[v][1], MaxErrCC_Info[v][2], MaxErrCC_Info[v][3],
                     MaxErrCC_Info[v][4], MaxErrCC_Info[v][5], MaxErrCC_Info[v][6],
                     MaxErrCC_Data[v][0], MaxErrCC_Data[v][1], MaxErrCC_Data[v][2], MaxErrCC_Data[v][3] );

      for (int v=0; v<NMag1; v++)
         if ( MaxErrFC_Data[v][3] != -__FLT_MIN__ )
            fprintf( File, "%6d%8d%8d  (%3d,%3d,%3d )%6d%16.7e%16.7e%16.7e%16.7e\n",
                     MaxErrFC_Info[v][0], MaxErrFC_Info[v][1], MaxErrFC_Info[v][2], MaxErrFC_Info[v][3],
                     MaxErrFC_Info[v][4], MaxErrFC_Info[v][5], MaxErrFC_Info[v][6],
                     MaxErrFC_Data[v][0], MaxErrFC_Data[v][1], MaxErrFC_Data[v][2], MaxErrFC_Data[v][3] );
   }


   fclose( File );


// verify that all target patches have been checked
   for (int lv=0; lv<NLEVEL; lv++)
   {
      for (int PID=0; PID<amr1.num[lv]; PID++)
      {
         if ( LeafOnly  &&  amr1.patch[lv][PID]->son != -1 )   continue;

         if ( amr1.patch[lv][PID]->check == false )
         {
            ErrorDetected = true;
            Aux_Message( stderr, "WARNING : patch %5d on level %d in input 1 has NOT been checked !!\n", PID, lv );
         }
      }

      for (int PID=0; PID<amr2.num[lv]; PID++)
      {
         if ( LeafOnly  &&  amr2.patch[lv][PID]->son != -1 )   continue;

         if ( amr2.patch[lv][PID]->check == false )
         {
            ErrorDetected = true;
            Aux_Message( stderr, "WARNING : patch %5d on level %d in input 2 has NOT been checked !!\n", PID, lv );
         }
      }
   }

   if ( ! ErrorDetected )  Aux_Message( stdout, "\n*** No error is detected in the grid data ***\n\n" );


   Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

   return ( ErrorDetected ) ? EXIT_FAILURE : EXIT_SUCCESS;

} // FUNCTION : CompareGridData



//-------------------------------------------------------------------------------------------------------
// Function    :  CompareParticleData
// Description :  Compare the particle data between two files
//
// Return      :  EXIT_SUCCESS : no error
//                EXIT_FAILURE : with errors
//-------------------------------------------------------------------------------------------------------
int CompareParticleData()
{

   Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// check if particle information is consistent in both files
   if ( NPar1 != NPar2 )
      Aux_Error( ERROR_INFO, "inconsistent numbers of particles (%ld vs %ld) !!\n", NPar1, NPar2 );

   if ( NPar1 > 0 )
   {
      if ( NParAttFlt1 != NParAttFlt2 )
         Aux_Error( ERROR_INFO, "inconsistent numbers of particle floating-point attributes (%d vs %d) !!\n", NParAttFlt1, NParAttFlt2 );

      if ( NParAttInt1 != NParAttInt2 )
         Aux_Error( ERROR_INFO, "inconsistent numbers of particle integer attributes (%d vs %d) !!\n", NParAttInt1, NParAttInt2 );

      if (  NParAttFlt1 != 0  &&  ( ParFltData1 == NULL || ParFltData2 == NULL )  )
         Aux_Error( ERROR_INFO, "particle floating-point arrays have not been allocated yet !!\n" );

      if (  NParAttInt1 != 0  &&  ( ParIntData1 == NULL || ParIntData2 == NULL )  )
         Aux_Error( ERROR_INFO, "particle integer arrays have not been allocated yet !!\n" );

      if ( Format1 == 2  &&  Format2 == 2 )
      {
         for (int v=0; v<NParAttFlt1; v++)
         {
            if ( ParAttFltLabel1 == NULL  ||  ParAttFltLabel2 == NULL )
               Aux_Error( ERROR_INFO, "ParAttFltLabel[%d] == NULL !!\n", v );

            if (  strcmp( ParAttFltLabel1[v], ParAttFltLabel2[v] )  )
               Aux_Error( ERROR_INFO, "inconsistent particle floating-point attribute labels ([%d]: \"%s\" vs \"%s\") !!\n",
                          v, ParAttFltLabel1[v], ParAttFltLabel2[v] );
         }

         for (int v=0; v<NParAttInt1; v++)
         {
            if ( ParAttIntLabel1 == NULL  ||  ParAttIntLabel2 == NULL )
               Aux_Error( ERROR_INFO, "ParAttIntLabel[%d] == NULL !!\n", v );

            if (  strcmp( ParAttIntLabel1[v], ParAttIntLabel2[v] )  )
               Aux_Error( ERROR_INFO, "inconsistent particle integer attribute labels ([%d]: \"%s\" vs \"%s\") !!\n",
                          v, ParAttIntLabel1[v], ParAttIntLabel2[v] );
         }
      } // if ( Format1 == 2  &&  Format2 == 2 )
   } // if ( NPar1 > 0 )


// sort particles by their x and y coordinates
   long *IdxTable1 = new long [NPar1];
   long *IdxTable2 = new long [NPar2];

// do not call SortParticle() if NPar==0 since ParDataX[0/1/2][] is ill-defined
   if ( NPar1 > 0 )
   {
      SortParticle( NPar1, ParFltData1[1], ParFltData1[2], ParFltData1[3], ParFltData1[4], IdxTable1 );
      SortParticle( NPar2, ParFltData2[1], ParFltData2[2], ParFltData2[3], ParFltData2[4], IdxTable2 );
   }


// compare data
   double     AbsErr, RelErr;
   long       AbsIntErr;
   long       ParID1, ParID2;
   real_par   FltData1, FltData2;
   long_par   IntData1, IntData2;
   bool       ErrorDetected = false;

   double     MaxFltErrPar_Data[NParAttFlt1][4];
   long       MaxFltErrPar_Info[NParAttFlt1][3];
   long       MaxIntErrPar_Data[NParAttInt1][3];
   long       MaxIntErrPar_Info[NParAttInt1][3];

   for (int v=0; v<NParAttFlt1; v++)   MaxFltErrPar_Data[v][3] = -__FLT_MIN__;
   for (int v=0; v<NParAttInt1; v++)   MaxIntErrPar_Data[v][2] = -0L;

   FILE *File = fopen( FileName_Out, "a" );

   fprintf( File, "\n\n" );
   fprintf( File, "#=============================================================================================================\n" );

   if ( Format1 == 2  &&  Format2 == 2 )
   {
      fprintf( File, "# Floating-point attribute list:\n" );
      for (int v=0; v<NParAttFlt1; v++)   fprintf( File, "# [%4d]: %s\n", v, ParAttFltLabel1[v] );
      fprintf( File, "\n" );
   }

   fprintf( File, "# %12s  %12s  %4s  %14s  %14s  %14s  %14s\n",
                  "ParID1", "ParID2", "Att", "Data1", "Data2", "AbsErr", "RelErr" );

   for (int v=0; v<NParAttFlt1; v++)
   for (long p=0; p<NPar1; p++)
   {
      ParID1   = IdxTable1[p];
      ParID2   = IdxTable2[p];
      FltData1 = ParFltData1[v][ParID1];
      FltData2 = ParFltData2[v][ParID2];

      AbsErr = FltData1 - FltData2;
      RelErr = AbsErr / ( 0.5*(FltData1+FltData2) );

      if ( FltData1 == 0.0  &&  FltData2 == 0.0 )  continue;

      if ( fabs(RelErr) >= TolErr  ||  !isfinite(RelErr) )
      {
         ErrorDetected = true;

         if ( MaxErrOnly ) {
            if ( fabs(RelErr) > fabs(MaxFltErrPar_Data[v][3]) ) {
               MaxFltErrPar_Info[v][0] = ParID1;
               MaxFltErrPar_Info[v][1] = ParID2;
               MaxFltErrPar_Info[v][2] = v;

               MaxFltErrPar_Data[v][0] = FltData1;
               MaxFltErrPar_Data[v][1] = FltData2;
               MaxFltErrPar_Data[v][2] = AbsErr;
               MaxFltErrPar_Data[v][3] = RelErr;
            }
         }

         else
            fprintf( File, "  %12ld  %12ld  %4d  %14.7e  %14.7e  %14.7e  %14.7e\n",
                     ParID1, ParID2, v, FltData1, FltData2, AbsErr, RelErr );
      }
   }

// output the maximum errors
   if ( MaxErrOnly )
   {
      for (int v=0; v<NParAttFlt1; v++)
         if ( MaxFltErrPar_Data[v][3] != -__FLT_MIN__ )
            fprintf( File, "  %12ld  %12ld  %4d  %14.7e  %14.7e  %14.7e  %14.7e\n",
                     MaxFltErrPar_Info[v][0], MaxFltErrPar_Info[v][1], (int)MaxFltErrPar_Info[v][2],
                     MaxFltErrPar_Data[v][0], MaxFltErrPar_Data[v][1], MaxFltErrPar_Data[v][2], MaxFltErrPar_Data[v][3] );
   }

   fprintf( File, "\n\n" );
   fprintf( File, "#=============================================================================================================\n" );

   if ( Format1 == 2  &&  Format2 == 2 )
   {
      fprintf( File, "# Integer attribute list:\n" );
      for (int v=0; v<NParAttInt1; v++)   fprintf( File, "# [%4d]: %s\n", v, ParAttIntLabel1[v] );
      fprintf( File, "\n" );
   }

   fprintf( File, "# %12s  %12s  %4s  %14s  %14s  %14s\n",
                  "ParID1", "ParID2", "Att", "Data1", "Data2", "AbsErr" );

   for (int v=0; v<NParAttInt1; v++)
   for (long p=0; p<NPar1; p++)
   {
      ParID1   = IdxTable1[p];
      ParID2   = IdxTable2[p];
      IntData1 = ParIntData1[v][ParID1];
      IntData2 = ParIntData2[v][ParID2];

      AbsIntErr = long(IntData1 - IntData2);

      if ( AbsIntErr != 0L )
      {
         ErrorDetected = true;

         if ( MaxErrOnly ) {
            if ( labs(AbsIntErr) > labs(MaxIntErrPar_Data[v][2]) ) {
               MaxIntErrPar_Info[v][0] = ParID1;
               MaxIntErrPar_Info[v][1] = ParID2;
               MaxIntErrPar_Info[v][2] = v;

               MaxIntErrPar_Data[v][0] = (long)IntData1;
               MaxIntErrPar_Data[v][1] = (long)IntData2;
               MaxIntErrPar_Data[v][2] = AbsIntErr;
            }
         }

         else
            fprintf( File, "  %12ld  %12ld  %4d  %14ld  %14ld  %14ld\n",
                     ParID1, ParID2, v, (long)IntData1, (long)IntData2, AbsIntErr );
      }
   }

// output the maximum errors
   if ( MaxErrOnly )
   {
      for (int v=0; v<NParAttInt1; v++)
         if ( MaxIntErrPar_Data[v][2] != -0L )
            fprintf( File, "  %12ld  %12ld  %4d  %14ld  %14ld  %14ld\n",
                     MaxIntErrPar_Info[v][0], MaxIntErrPar_Info[v][1], (int)MaxIntErrPar_Info[v][2],
                     MaxIntErrPar_Data[v][0], MaxIntErrPar_Data[v][1], MaxIntErrPar_Data[v][2] );
   }

   fclose( File );

   delete [] IdxTable1;
   delete [] IdxTable2;

   if ( ! ErrorDetected )  Aux_Message( stdout, "\n*** No error is detected in the particle data ***\n\n" );


   Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

   return ( ErrorDetected ) ? EXIT_FAILURE : EXIT_SUCCESS;

} // FUNCTION : CompareParticleData



//-------------------------------------------------------------------------------------------------------
// Function    :  FreeMemory
// Description :  Free memory
//-------------------------------------------------------------------------------------------------------
void FreeMemory()
{
   Aux_Message( stdout, "%s ... ", __FUNCTION__ );


   Aux_DeallocateArray2D<real_par>( ParFltData1 );
   Aux_DeallocateArray2D<real_par>( ParFltData2 );
   Aux_DeallocateArray2D<long_par>( ParIntData1 );
   Aux_DeallocateArray2D<long_par>( ParIntData2 );

   delete [] FieldLabel1;
   delete [] FieldLabel2;
   delete [] MagLabel1;
   delete [] MagLabel2;
   delete [] ParAttFltLabel1;
   delete [] ParAttFltLabel2;
   delete [] ParAttIntLabel1;
   delete [] ParAttIntLabel2;


   Aux_Message( stdout, "done\n" );

} // FUNCTION : FreeMemory



//-------------------------------------------------------------------------------------------------------
// Function    :  main
// Description :
//-------------------------------------------------------------------------------------------------------
int main( int argc, char ** argv )
{

   ReadOption( argc, argv );

   CheckParameter();

   LoadData( FileName_In1, amr1, Format1, NField1, NMag1, NParAttFlt1, NParAttInt1, NPar1,
             ParFltData1, ParIntData1, FieldLabel1, MagLabel1, ParAttFltLabel1, ParAttIntLabel1 );
   LoadData( FileName_In2, amr2, Format2, NField2, NMag2, NParAttFlt2, NParAttInt2, NPar2,
             ParFltData2, ParIntData2, FieldLabel2, MagLabel2, ParAttFltLabel2, ParAttIntLabel2 );

   int ErrorDetected = false;

   ErrorDetected |= CompareGridData();
   ErrorDetected |= CompareParticleData();

   FreeMemory();

   Aux_Message( stdout, "Program terminated successfully\n" );

   return ( ErrorDetected ) ? EXIT_FAILURE : EXIT_SUCCESS;

} // FUNCTION : main
