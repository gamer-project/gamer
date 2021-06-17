#include "CompareData.h"



AMR_t    amr1, amr2;
char    *FileName_In1=NULL, *FileName_In2=NULL, *FileName_Out=NULL;
bool     WithPot1=false, WithPot2=false, WithMagCC1=false, WithMagCC2=false, WithMagFC1=false, WithMagFC2=false;
bool     UseCorner=false, WithPar1=false, WithPar2=false;
int      WithParDens1=0, WithParDens2=0, NParVarOut1=-1, NParVarOut2=-1, Format1=-1, Format2=-1;
long     NPar1=-1, NPar2=-1;
double   TolErr=__FLT_MIN__;
real   **ParData1=NULL, **ParData2=NULL;




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

   Aux_Message( stdout, "%s ... \n", __FUNCTION__ );


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

   Aux_Message( stdout, "%s ... \n", __FUNCTION__ );


// verify that the total number of patches at each level of amr1 and amr2 are the same
   for (int lv=0; lv<NLEVEL; lv++)
   {
      if ( amr1.num[lv] != amr2.num[lv] )
         Aux_Error( ERROR_INFO, "amr1.num[%d] (%d) != amr2.num[%d] (%d) !!\n",
                    lv, amr1.num[lv], lv, amr2.num[lv] );
   }


// verify that the simulation domains of amr1 and amr2 are the same
   for (int d=0; d<3; d++)
   {
      if ( amr1.nx0_tot[d] != amr2.nx0_tot[d] )
         Aux_Error( ERROR_INFO, "amr1.nx0_tot[%d] (%d) != amr2.nx0_tot[%d] (%d) !!\n",
                    d, amr1.nx0_tot[d], d, amr2.nx0_tot[d] );
   }


// check whether both inputs store the potential data or not
   if ( WithPot1 != WithPot2 )
      Aux_Message( stderr, "WARNING : one of the input files does NOT store the potential data !!\n"
                           "          --> Potential data will not be compared ...\n" );


// check whether both inputs store the particle (or total) density data or not
   if ( WithParDens1 != WithParDens2 )
      Aux_Message( stderr, "WARNING : one of the input files does NOT store the particle (or total) density data !!\n"
                           "          --> Particle density data will not be compared ...\n" );


// check whether both inputs store the B field
   if ( WithMagCC1 != WithMagCC2 )
      Aux_Message( stderr, "WARNING : one of the input files does NOT store the cell-centered B field data !!\n"
                           "          --> Cell-centered B field data will not be compared ...\n" );

   if ( WithMagFC1 != WithMagFC2 )
      Aux_Error( ERROR_INFO, "one of the input files does NOT store the face-centered B field data !!\n" );


// compare data
   double Data1, Data2, AbsErr, RelErr;
   int PID1, PID2;
   int *Cr1, *Cr2;

   FILE *File = fopen( FileName_Out, "w" );

   fprintf( File, "%5s%8s%8s  (%3s,%3s,%3s )%6s%16s%16s%16s%16s\n",
                  "Level", "PID1", "PID2", "i", "j", "k", "Comp", "Data1", "Data2", "AbsErr", "RelErr" );


   for (int lv=0; lv<NLEVEL; lv++)
   {
      Aux_Message( stdout, "  Comparing level %d ... ", lv );

      for (PID1=0; PID1<amr1.num[lv]; PID1++)
      {
//       only compare patches without son
         if ( amr1.patch[lv][PID1]->son == -1 )
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
               PID2 = PID1;

            for (int k=0; k<PATCH_SIZE; k++)
            for (int j=0; j<PATCH_SIZE; j++)
            for (int i=0; i<PATCH_SIZE; i++)
            {
//             fluid data
               for (int v=0; v<NCOMP_TOTAL; v++)
               {
                  Data1  = amr1.patch[lv][PID1]->fluid[v][k][j][i];
                  Data2  = amr2.patch[lv][PID2]->fluid[v][k][j][i];
                  AbsErr = Data1 - Data2;
                  RelErr = AbsErr / ( 0.5*(Data1+Data2) );

                  if ( Data1 == 0.0  &&  Data2 == 0.0 )  continue;

                  if ( fabs( RelErr ) >= TolErr  ||  !isfinite( RelErr )  )
                     fprintf( File, "%5d%8d%8d  (%3d,%3d,%3d )%6d%16.7e%16.7e%16.7e%16.7e\n",
                              lv, PID1, PID2, i, j, k, v, Data1, Data2, AbsErr, RelErr );
               }

//             gravitational potential
               if ( WithPot1 && WithPot2 )
               {
                  Data1  = amr1.patch[lv][PID1]->pot[k][j][i];
                  Data2  = amr2.patch[lv][PID2]->pot[k][j][i];
                  AbsErr = Data1 - Data2;
                  RelErr = AbsErr / ( 0.5*(Data1+Data2) );

                  if ( Data1 == 0.0  &&  Data2 == 0.0 )  continue;

                  if ( fabs( RelErr ) >= TolErr  ||  !isfinite( RelErr )  )
                     fprintf( File, "%5d%8d%8d  (%3d,%3d,%3d )%6d%16.7e%16.7e%16.7e%16.7e\n",
                              lv, PID1, PID2, i, j, k, 1000, Data1, Data2, AbsErr, RelErr );
               }

//             particle density
               if ( WithParDens1 > 0  &&  WithParDens1 == WithParDens2 )
               {
                  Data1  = amr1.patch[lv][PID1]->par_dens[k][j][i];
                  Data2  = amr2.patch[lv][PID2]->par_dens[k][j][i];
                  AbsErr = Data1 - Data2;
                  RelErr = AbsErr / ( 0.5*(Data1+Data2) );

                  if ( Data1 == 0.0  &&  Data2 == 0.0 )  continue;

                  if ( fabs( RelErr ) >= TolErr  ||  !isfinite( RelErr )  )
                     fprintf( File, "%5d%8d%8d  (%3d,%3d,%3d )%6d%16.7e%16.7e%16.7e%16.7e\n",
                              lv, PID1, PID2, i, j, k, 1001, Data1, Data2, AbsErr, RelErr );
               }

//             cell-centered B field
               if ( WithMagCC1 && WithMagCC2 )
               {
                  for (int v=0; v<NCOMP_MAG; v++)
                  {
                     Data1  = amr1.patch[lv][PID1]->mag_cc[v][k][j][i];
                     Data2  = amr2.patch[lv][PID2]->mag_cc[v][k][j][i];
                     AbsErr = Data1 - Data2;
                     RelErr = AbsErr / ( 0.5*(Data1+Data2) );

                     if ( Data1 == 0.0  &&  Data2 == 0.0 )  continue;

                     if ( fabs( RelErr ) >= TolErr  ||  !isfinite( RelErr )  )
                        fprintf( File, "%5d%8d%8d  (%3d,%3d,%3d )%6d%16.7e%16.7e%16.7e%16.7e\n",
                                 lv, PID1, PID2, i, j, k, 2000+v, Data1, Data2, AbsErr, RelErr );
                  }
               }
            } // i,j,k

//          face-centered B field
//          Bx
            if ( WithMagFC1 && WithMagFC2 )
            for (int k=0; k<PS1;   k++)
            for (int j=0; j<PS1;   j++)
            for (int i=0; i<PS1P1; i++)
            {
               Data1  = amr1.patch[lv][PID1]->mag_fc[0][ IDX321_BX(i,j,k) ];
               Data2  = amr2.patch[lv][PID2]->mag_fc[0][ IDX321_BX(i,j,k) ];
               AbsErr = Data1 - Data2;
               RelErr = AbsErr / ( 0.5*(Data1+Data2) );

               if ( Data1 == 0.0  &&  Data2 == 0.0 )  continue;

               if ( fabs( RelErr ) >= TolErr  ||  !isfinite( RelErr )  )
                  fprintf( File, "%5d%8d%8d  (%3d,%3d,%3d )%6d%16.7e%16.7e%16.7e%16.7e\n",
                           lv, PID1, PID2, i, j, k, 2003, Data1, Data2, AbsErr, RelErr );
            }

//          By
            if ( WithMagFC1 && WithMagFC2 )
            for (int k=0; k<PS1;   k++)
            for (int j=0; j<PS1P1; j++)
            for (int i=0; i<PS1;   i++)
            {
               Data1  = amr1.patch[lv][PID1]->mag_fc[1][ IDX321_BY(i,j,k) ];
               Data2  = amr2.patch[lv][PID2]->mag_fc[1][ IDX321_BY(i,j,k) ];
               AbsErr = Data1 - Data2;
               RelErr = AbsErr / ( 0.5*(Data1+Data2) );

               if ( Data1 == 0.0  &&  Data2 == 0.0 )  continue;

               if ( fabs( RelErr ) >= TolErr  ||  !isfinite( RelErr )  )
                  fprintf( File, "%5d%8d%8d  (%3d,%3d,%3d )%6d%16.7e%16.7e%16.7e%16.7e\n",
                           lv, PID1, PID2, i, j, k, 2004, Data1, Data2, AbsErr, RelErr );
            }

//          Bz
            if ( WithMagFC1 && WithMagFC2 )
            for (int k=0; k<PS1P1; k++)
            for (int j=0; j<PS1;   j++)
            for (int i=0; i<PS1;   i++)
            {
               Data1  = amr1.patch[lv][PID1]->mag_fc[2][ IDX321_BZ(i,j,k) ];
               Data2  = amr2.patch[lv][PID2]->mag_fc[2][ IDX321_BZ(i,j,k) ];
               AbsErr = Data1 - Data2;
               RelErr = AbsErr / ( 0.5*(Data1+Data2) );

               if ( Data1 == 0.0  &&  Data2 == 0.0 )  continue;

               if ( fabs( RelErr ) >= TolErr  ||  !isfinite( RelErr )  )
                  fprintf( File, "%5d%8d%8d  (%3d,%3d,%3d )%6d%16.7e%16.7e%16.7e%16.7e\n",
                           lv, PID1, PID2, i, j, k, 2005, Data1, Data2, AbsErr, RelErr );
            }


            amr1.patch[lv][PID1]->check = true;
            amr2.patch[lv][PID2]->check = true;

         } // if ( amr1[lv][PID1]->son == -1 )

      } // for (PID1=0; PID1<amr1.num[lv]; PID1++)

      Aux_Message( stdout, "done\n" );

   } // for (int lv=0; lv<NLEVEL; lv++)


   fclose( File );


// verify that all patches without son have been checked
   for (int lv=0; lv<NLEVEL; lv++)
   for (int PID=0; PID<amr1.num[lv]; PID++)
   {
      if ( amr1.patch[lv][PID]->son == -1  &&  amr1.patch[lv][PID]->check == false )
         Aux_Message( stderr, "WARNING : patch %5d at level %d in input 1 has NOT been checked !!\n", PID, lv );

      if ( amr2.patch[lv][PID]->son == -1  &&  amr2.patch[lv][PID]->check == false )
         Aux_Message( stderr, "WARNING : patch %5d at level %d in input 2 has NOT been checked !!\n", PID, lv );
   }


   Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : CompareGridData



//-------------------------------------------------------------------------------------------------------
// Function    :  CompareParticleData
// Description :  Compare the particle data between two files
//-------------------------------------------------------------------------------------------------------
void CompareParticleData()
{

   Aux_Message( stdout, "%s ... \n", __FUNCTION__ );


// nothing to do if there are not particle data in both files
   if ( !WithPar1  &&  !WithPar2 )
   {
      Aux_Message( stdout, "   No particle data are found in both files\n" );
      return;
   }


// nothing to do if there are no particles at all
// --> we must return here since ParDataX[...] is ill-defined if NParX == 0
   if ( NPar1 == 0  &&  NPar2 == 0 )
   {
      Aux_Message( stdout, "   No particle data are found in both files\n" );
      return;
   }


// check if particle information is consistent in both files
   if ( WithPar1 != WithPar2 )
      Aux_Error( ERROR_INFO, "One of the input files does NOT have particles !!\n" );

   if ( NPar1 != NPar2 )
      Aux_Error( ERROR_INFO, "Numbers of particles in the two files are inconsistent (%ld vs %ld) !!\n",
                 NPar1, NPar2 );

   if ( NParVarOut1 != NParVarOut2 )
      Aux_Error( ERROR_INFO, "Numbers of particle attributes in the two files are inconsistent (%d vs %d) !!\n",
                 NParVarOut1, NParVarOut2 );

   if ( ParData1 == NULL  ||  ParData2 == NULL )
      Aux_Error( ERROR_INFO, "Particle arrays have not been allocated yet !!\n" );


// sort particles by their x and y coordinates
   long *IdxTable1 = new long [NPar1];
   long *IdxTable2 = new long [NPar2];

   SortParticle( NPar1, ParData1[1], ParData1[2], ParData1[3], ParData1[4], IdxTable1 );
   SortParticle( NPar2, ParData2[1], ParData2[2], ParData2[3], ParData2[4], IdxTable2 );


// compare data
   double AbsErr, RelErr;
   long   ParID1, ParID2;
   real   Data1, Data2;

   FILE *File = fopen( FileName_Out, "a" );

   fprintf( File, "\n\n" );
   fprintf( File, "=============================================================================================================\n" );
   fprintf( File, "  %12s  %12s  %4s  %14s  %14s  %14s  %14s\n",
                  "ParID1", "ParID2", "Var", "Data1", "Data2", "AbsErr", "RelErr" );

   for (int v=0; v<NParVarOut1; v++)
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
         fprintf( File, "  %12ld  %12ld  %4d  %14.7e  %14.7e  %14.7e  %14.7e\n",
                  ParID1, ParID2, v, Data1, Data2, AbsErr, RelErr );
   }

   fclose( File );

   delete [] IdxTable1;
   delete [] IdxTable2;


   Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : CompareParticleData



//-------------------------------------------------------------------------------------------------------
// Function    :  main
// Description :
//-------------------------------------------------------------------------------------------------------
int main( int argc, char ** argv )
{

   ReadOption( argc, argv );

   CheckParameter();

   LoadData( amr1, FileName_In1, WithPot1, WithParDens1, WithPar1, NParVarOut1, NPar1, ParData1, WithMagCC1, WithMagFC1, Format1 );
   LoadData( amr2, FileName_In2, WithPot2, WithParDens2, WithPar2, NParVarOut2, NPar2, ParData2, WithMagCC2, WithMagFC2, Format2 );

   CompareGridData();
   CompareParticleData();

   Aux_DeallocateArray2D( ParData1 );
   Aux_DeallocateArray2D( ParData2 );

   Aux_Message( stdout, "Program terminated successfully\n" );

   return 0;

} // FUNCTION : main
