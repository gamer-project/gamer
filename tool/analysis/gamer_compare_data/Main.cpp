#include "CompareData.h"



// no need to share with other files
static AMR_t    amr1, amr2;
static char    *FileName_In1=NULL, *FileName_In2=NULL, *FileName_Out=NULL;
static bool     UseCorner=false;
static int      NField1=-1, NField2=-1, NMag1=-1, NMag2=-1, NParAtt1=-1, NParAtt2=-1, Format1=-1, Format2=-1;
static long     NPar1=-1, NPar2=-1;
static double   TolErr=__FLT_MIN__;
static real   **ParData1=NULL, **ParData2=NULL;
static char  (*FieldLabel1)[MAX_STRING]=NULL, (*MagLabel1)[MAX_STRING]=NULL, (*ParAttLabel1)[MAX_STRING]=NULL;
static char  (*FieldLabel2)[MAX_STRING]=NULL, (*MagLabel2)[MAX_STRING]=NULL, (*ParAttLabel2)[MAX_STRING]=NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  ReadOption
// Description :  Read options from the command line
//-------------------------------------------------------------------------------------------------------
void ReadOption( int argc, char **argv )
{

   int c;

   while( (c = getopt(argc, argv, "hci:j:o:e:")) != -1 )
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
         case 'h':
         case '?': cerr << endl << "usage: " << argv[0]
                        << " [-h (for help)] [-i Input FileName1] [-j Input FileName2] [-o Output FileName]"
                        << endl << "                          "
                        << " [-e Tolerant Error] [-c (compare patches with the same corner coordinates) [off]]"
                        << endl
                        << endl << endl;
                   exit( 1 );
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

   if ( TolErr == 1.e-20 )
      Aux_Message( stderr, "WARNING : please provide the tolerant error (-e Tolerant Error) !!\n" );

   if ( UseCorner == false )
      Aux_Message( stderr, "WARNING : please make sure that the order of patches stored in the input files are the same.\n"
                           "          Otherwise, you should turn on the option \"-c\" !!\n" );


   Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : CheckParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  CompareGridData
// Description :  Compare the grid data between two files
//-------------------------------------------------------------------------------------------------------
void CompareGridData()
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
      Aux_Message( stdout, "  Comparing level %d ...", lv );

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

               if ( fabs( RelErr ) >= TolErr  ||  !isfinite( RelErr )  )
               {
                  ErrorDetected = true;

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

               if ( fabs( RelErr ) >= TolErr  ||  !isfinite( RelErr )  )
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

} // FUNCTION : CompareGridData



//-------------------------------------------------------------------------------------------------------
// Function    :  CompareParticleData
// Description :  Compare the particle data between two files
//-------------------------------------------------------------------------------------------------------
void CompareParticleData()
{

   Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// check if particle information is consistent in both files
   if ( NPar1 != NPar2 )
      Aux_Error( ERROR_INFO, "inconsistent numbers of particles (%ld vs %ld) !!\n", NPar1, NPar2 );

   if ( NPar1 > 0 )
   {
      if ( NParAtt1 != NParAtt2 )
         Aux_Error( ERROR_INFO, "inconsistent numbers of particle attributes (%d vs %d) !!\n", NParAtt1, NParAtt2 );

      if ( ParData1 == NULL  ||  ParData2 == NULL )
         Aux_Error( ERROR_INFO, "particle arrays have not been allocated yet !!\n" );

      if ( Format1 == 2  &&  Format2 == 2 )
      {
         for (int v=0; v<NParAtt1; v++)
         {
            if ( ParAttLabel1 == NULL  ||  ParAttLabel2 == NULL )
               Aux_Error( ERROR_INFO, "ParAttLabel[%d] == NULL !!\n", v );

            if (  strcmp( ParAttLabel1[v], ParAttLabel2[v] )  )
               Aux_Error( ERROR_INFO, "inconsistent particle attribute labels ([%d]: \"%s\" vs \"%s\") !!\n",
                          v, ParAttLabel1[v], ParAttLabel2[v] );
         }
      } // if ( Format1 == 2  &&  Format2 == 2 )
   } // if ( NPar1 > 0 )


// sort particles by their x and y coordinates
   long *IdxTable1 = new long [NPar1];
   long *IdxTable2 = new long [NPar2];

// do not call SortParticle() if NPar==0 since ParDataX[0/1/2][] is ill-defined
   if ( NPar1 > 0 )
   {
      SortParticle( NPar1, ParData1[1], ParData1[2], ParData1[3], ParData1[4], IdxTable1 );
      SortParticle( NPar2, ParData2[1], ParData2[2], ParData2[3], ParData2[4], IdxTable2 );
   }


// compare data
   double AbsErr, RelErr;
   long   ParID1, ParID2;
   real   Data1, Data2;
   bool   ErrorDetected = false;

   FILE *File = fopen( FileName_Out, "a" );

   fprintf( File, "\n\n" );
   fprintf( File, "#=============================================================================================================\n" );

   if ( Format1 == 2  &&  Format2 == 2 )
   {
      fprintf( File, "# Attribute list:\n" );
      for (int v=0; v<NParAtt1; v++)   fprintf( File, "# [%4d]: %s\n", v, ParAttLabel1[v] );
      fprintf( File, "\n" );
   }

   fprintf( File, "# %12s  %12s  %4s  %14s  %14s  %14s  %14s\n",
                  "ParID1", "ParID2", "Att", "Data1", "Data2", "AbsErr", "RelErr" );

   for (int v=0; v<NParAtt1; v++)
   for (long p=0; p<NPar1; p++)
   {
      ParID1 = IdxTable1[p];
      ParID2 = IdxTable2[p];
      Data1  = ParData1[v][ParID1];
      Data2  = ParData2[v][ParID2];

      AbsErr = Data1 - Data2;
      RelErr = AbsErr / ( 0.5*(Data1+Data2) );

      if ( Data1 == 0.0  &&  Data2 == 0.0 )  continue;

      if ( fabs( RelErr ) >= TolErr  ||  !isfinite( RelErr )  )
      {
         ErrorDetected = true;

         fprintf( File, "  %12ld  %12ld  %4d  %14.7e  %14.7e  %14.7e  %14.7e\n",
                  ParID1, ParID2, v, Data1, Data2, AbsErr, RelErr );
      }
   }

   fclose( File );

   delete [] IdxTable1;
   delete [] IdxTable2;

   if ( ! ErrorDetected )  Aux_Message( stdout, "\n*** No error is detected in the particle data ***\n\n" );


   Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : CompareParticleData



//-------------------------------------------------------------------------------------------------------
// Function    :  FreeMemory
// Description :  Free memory
//-------------------------------------------------------------------------------------------------------
void FreeMemory()
{
   Aux_Message( stdout, "%s ... ", __FUNCTION__ );


   Aux_DeallocateArray2D( ParData1 );
   Aux_DeallocateArray2D( ParData2 );

   delete [] FieldLabel1;
   delete [] FieldLabel2;
   delete [] MagLabel1;
   delete [] MagLabel2;
   delete [] ParAttLabel1;
   delete [] ParAttLabel2;


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

   LoadData( FileName_In1, amr1, Format1, NField1, NMag1, NParAtt1, NPar1, ParData1,
             FieldLabel1, MagLabel1, ParAttLabel1 );
   LoadData( FileName_In2, amr2, Format2, NField2, NMag2, NParAtt2, NPar2, ParData2,
             FieldLabel2, MagLabel2, ParAttLabel2 );

   CompareGridData();
   CompareParticleData();

   FreeMemory();

   Aux_Message( stdout, "Program terminated successfully\n" );

   return 0;

} // FUNCTION : main
