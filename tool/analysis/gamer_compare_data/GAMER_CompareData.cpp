#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <stdarg.h>
#include "TypeDef.h"

using namespace std;

#define ERROR_INFO         __FILE__, __LINE__, __FUNCTION__
#define CUBE( a )       ( (a)*(a)*(a) )

void ReadOption( int argc, char **argv );
void CheckParameter();
void LoadData( GAMER_t &patch, const char *FileName_In, bool &WithPot, int &WithParDens, bool &WithPar, int &NParVarOut,
               long &NPar, real **&ParData );
void Aux_Error( const char *File, const int Line, const char *Func, const char *Format, ... );
void Aux_AllocateArray2D( real** &Array, const int J, const int I );
void Aux_DeallocateArray2D( real** &Array );
void CompareVar( const char *VarName, const int RestartVar, const int RuntimeVar, const bool Fatal );
void Mis_Heapsort( const long N, real Array[], long IdxTable[] );
void Heapsort_SiftDown( const long L, const long R, real Array[], long IdxTable[] );
void Load_Parameter_After_2000( FILE *File, const int FormatVersion, bool &WithPot, int &WithParDens,
                                const long HeaderOffset_Makefile, const long HeaderOffset_Constant,
                                const long HeaderOffset_Parameter, int *NX0_Tot, int *NGPU_X,
                                bool &WithPar, int &NParVarOut );
void CompareGridData();
void CompareParticleData();
void SortParticle( const long NPar, const real *PosX, const real *PosY, const real *PosZ, const real *VelX, long *IdxTable );

GAMER_t  patch1, patch2;
char    *FileName_In1=NULL, *FileName_In2=NULL, *FileName_Out=NULL;
bool     WithPot1, WithPot2;

double   TolErr     = __FLT_MIN__;
bool     UseCorner  = false;

bool     WithPar1=false, WithPar2=false;
int      WithParDens1=0, WithParDens2=0;
long     NPar1=-1, NPar2=-1;
int      NParVarOut1=-1, NParVarOut2=-1;
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

   cout << "CheckParameter ... " << endl;


   if ( FileName_In1 == NULL )
   {
      cerr << "ERROR : please provide the name of the input file one (-i FileName1) !!" << endl;
      exit( 1 );
   }

   if ( FileName_In2 == NULL )
   {
      cerr << "ERROR : please provide the name of the input file two (-j FileName2) !!" << endl;
      exit( 1 );
   }

   if ( FileName_Out == NULL )
   {
      cerr << "ERROR : please provide the name of the output file (-o FileNmae) !!" << endl;
      exit( 1 );
   }

   if ( TolErr == 1.e-20 )
   {
      cerr << "WARNING : please provide the tolerant error (-e Tolerant Error) !!" << endl;
   }

   if ( UseCorner == false )
   {
      cerr << "WARNING : please make sure that the order of patches storing in the input files are the same."
           << endl
           << "          Otherwise, you should turn on the option \"-c\" !!" << endl;
   }


   cout << "CheckParameter ... done" << endl;

} // FUNCTION : CheckParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  CompareGridData
// Description :  Compare the grid data between two files
//-------------------------------------------------------------------------------------------------------
void CompareGridData()
{

   cout << "CompareGridData ... " << endl;


// verify that the total number of patches at each level of patch1 and patch2 are the same
   for (int lv=0; lv<NLEVEL; lv++)
   {
      if ( patch1.num[lv] != patch2.num[lv] )
      {
         fprintf( stderr, "ERROR : patch1.num[%d] (%d) != patch2.num[%d] (%d) !!\n",
                  lv, patch1.num[lv], lv, patch2.num[lv] );
         exit( 1 );
      }
   }


// verify that the simulation domains of patch1 and patch2 are the same
   for (int d=0; d<3; d++)
   {
      if ( patch1.nx0_tot[d] != patch2.nx0_tot[d] )
      {
         fprintf( stderr, "ERROR : patch1.nx0_tot[%d] (%d) != patch2.nx0_tot[%d] (%d) !!\n",
                  d, patch1.nx0_tot[d], d, patch2.nx0_tot[d] );
         exit( 1 );
      }
   }


// verify that the domain decomposition of patch1 and patch2 are the same
   if ( UseCorner == false )
   {
      for (int d=0; d<3; d++)
      {
         if ( patch1.ngpu_x[d] != patch2.ngpu_x[d] )
         {
            fprintf( stderr, "WARNING : patch1.ngpu_x[%d] (%d) != patch2.ngpu_x[%d] (%d) !!\n",
                     d, patch1.ngpu_x[d], d, patch2.ngpu_x[d] );
            fprintf( stderr, "          --> The option \"-c\" is turned on automatically\n" );

            UseCorner = true;
            break;
         }
      }
   }


// check that whether both inputs store the potential data or not
   if ( WithPot1 != WithPot2 )
   {
      fprintf( stderr, "WARNING : one of the input files does NOT store the potential data !!\n" );
      fprintf( stderr, "          --> Potential data will not be compared ...\n" );
   }


// check that whether both inputs store the particle (or total) density data or not
   if ( WithParDens1 != WithParDens2 )
   {
      fprintf( stderr, "WARNING : one of the input files does NOT store the particle (or total) density data !!\n" );
      fprintf( stderr, "          --> Particle density data will not be compared ...\n" );
   }


// compare data
   double Data1, Data2, AbsErr, RelErr;
   int PID1, PID2;
   int *Cr1, *Cr2;

   FILE *File = fopen( FileName_Out, "w" );

   fprintf( File, "%5s%8s%8s  (%3s,%3s,%3s )%6s%16s%16s%16s%16s\n",
                  "Level", "PID1", "PID2", "i", "j", "k", "Comp", "Data1", "Data2", "AbsErr", "RelErr" );


   for (int lv=0; lv<NLEVEL; lv++)
   {
      cout << "  Comparing level " << lv << " ... " << flush;

      for (PID1=0; PID1<patch1.num[lv]; PID1++)
      {
//       only compare patches without son
         if ( patch1.ptr[lv][PID1]->son == -1 )
         {

//          set the targeted patch ID in the second input
            if ( UseCorner )
            {
               Cr1 = patch1.ptr[lv][PID1]->corner;

               for (PID2=0; PID2<patch2.num[lv]; PID2++)
               {
                  Cr2 = patch2.ptr[lv][PID2]->corner;

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
                  Data1  = patch1.ptr[lv][PID1]->fluid[v][k][j][i];
                  Data2  = patch2.ptr[lv][PID2]->fluid[v][k][j][i];
                  AbsErr = Data1 - Data2;
                  RelErr = AbsErr / Data2;

                  if ( Data1 == 0.0  &&  Data2 == 0.0 )  continue;

                  if ( fabs( RelErr ) >= TolErr  ||  !isfinite( RelErr )  )
                     fprintf( File, "%5d%8d%8d  (%3d,%3d,%3d )%6d%16.7e%16.7e%16.7e%16.7e\n",
                              lv, PID1, PID2, i, j, k, v, Data1, Data2, AbsErr, RelErr );
               }

//             gravitational potential
               if ( WithPot1 && WithPot2 )
               {
                  Data1  = patch1.ptr[lv][PID1]->pot[k][j][i];
                  Data2  = patch2.ptr[lv][PID2]->pot[k][j][i];
                  AbsErr = Data1 - Data2;
                  RelErr = AbsErr / Data2;

                  if ( Data1 == 0.0  &&  Data2 == 0.0 )  continue;

                  if ( fabs( RelErr ) >= TolErr  ||  !isfinite( RelErr )  )
                     fprintf( File, "%5d%8d%8d  (%3d,%3d,%3d )%6d%16.7e%16.7e%16.7e%16.7e\n",
                              lv, PID1, PID2, i, j, k, NCOMP_TOTAL, Data1, Data2, AbsErr, RelErr );
               }

//             particle density
               if ( WithParDens1 > 0  &&  WithParDens1 == WithParDens2 )
               {
                  Data1  = patch1.ptr[lv][PID1]->par_dens[k][j][i];
                  Data2  = patch2.ptr[lv][PID2]->par_dens[k][j][i];
                  AbsErr = Data1 - Data2;
                  RelErr = AbsErr / Data2;

                  if ( Data1 == 0.0  &&  Data2 == 0.0 )  continue;

                  if ( fabs( RelErr ) >= TolErr  ||  !isfinite( RelErr )  )
                     fprintf( File, "%5d%8d%8d  (%3d,%3d,%3d )%6d%16.7e%16.7e%16.7e%16.7e\n",
                              lv, PID1, PID2, i, j, k, NCOMP_TOTAL+1, Data1, Data2, AbsErr, RelErr );
               }
            }

            patch1.ptr[lv][PID1]->check = true;
            patch2.ptr[lv][PID2]->check = true;

         } // if ( patch1[lv][PID1]->son == -1 )

      } // for (PID1=0; PID1<patch1.num[lv]; PID1++)

      cout << "done" << endl;

   } // for (int lv=0; lv<NLEVEL; lv++)


   fclose( File );


// verify that all patches without son have been checked
   for (int lv=0; lv<NLEVEL; lv++)
   for (int PID=0; PID<patch1.num[lv]; PID++)
   {
      if ( patch1.ptr[lv][PID]->son == -1  &&  patch1.ptr[lv][PID]->check == false )
         fprintf( stderr, "WARNING : patch %5d at level %d in input 1 has NOT been checked !!\n", PID, lv );

      if ( patch2.ptr[lv][PID]->son == -1  &&  patch2.ptr[lv][PID]->check == false )
         fprintf( stderr, "WARNING : patch %5d at level %d in input 2 has NOT been checked !!\n", PID, lv );
   }


   cout << "CompareGridData ... done" << endl;

} // FUNCTION : CompareGridData



//-------------------------------------------------------------------------------------------------------
// Function    :  CompareParticleData
// Description :  Compare the particle data between two files
//-------------------------------------------------------------------------------------------------------
void CompareParticleData()
{

   cout << "CompareParticleData ... " << endl;


// nothing to do if there are not particle data in both files
   if ( !WithPar1  &&  !WithPar2 )
   {
      fprintf( stdout, "   No particle data are found in both files\n" );
      return;
   }


// nothing to do if there are no particles at all
// --> we must return here since ParDataX[...] is ill-defined if NParX == 0
   if ( NPar1 == 0  &&  NPar2 == 0 )
   {
      fprintf( stdout, "   No particle data are found in both files\n" );
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
      RelErr = AbsErr / Data2;

      if ( Data1 == 0.0  &&  Data2 == 0.0 )  continue;

      if ( fabs( RelErr ) >= TolErr  ||  !isfinite( RelErr )  )
         fprintf( File, "  %12ld  %12ld  %4d  %14.7e  %14.7e  %14.7e  %14.7e\n",
                  ParID1, ParID2, v, Data1, Data2, AbsErr, RelErr );
   }

   fclose( File );

   delete [] IdxTable1;
   delete [] IdxTable2;


   cout << "CompareParticleData ... done" << endl;

} // FUNCTION : CompareParticleData



//-------------------------------------------------------------------------------------------------------
// Function    :  SortParticle
// Description :  Sort particles by their x, y, z, and vx
//-------------------------------------------------------------------------------------------------------
void SortParticle( const long NPar, const real *PosX, const real *PosY, const real *PosZ, const real *VelX, long *IdxTable )
{

   long NParSameX, NParSameY, NParSameZ;

// 0. back up the PosX array since we don't want to modify it during the sorting
   real *PosX_Sorted = new real [NPar];

   memcpy( PosX_Sorted, PosX, sizeof(real)*NPar );


// 1. sort by x
   Mis_Heapsort( NPar, PosX_Sorted, IdxTable );


// 2. sort by y
   for (long x=0; x<NPar-1; x++)
   {
      NParSameX = 1;    // 1 --> itself

//    find all particles with the same x coordinate
      while ( x+NParSameX<NPar  &&  PosX_Sorted[x] == PosX_Sorted[x+NParSameX] )    NParSameX ++;

      if ( NParSameX > 1 )
      {
         long *SortByX_IdxTable = new long [NParSameX];
         long *SortByY_IdxTable = new long [NParSameX];
         real *PosY_Sorted      = new real [NParSameX];

         for (long y=0; y<NParSameX; y++)
         {
            SortByX_IdxTable[y] = IdxTable[ x + y ];
            PosY_Sorted     [y] = PosY[ SortByX_IdxTable[y] ];
         }

         Mis_Heapsort( NParSameX, PosY_Sorted, SortByY_IdxTable );

         for (long y=0; y<NParSameX; y++)    IdxTable[ x + y ] = SortByX_IdxTable[ SortByY_IdxTable[y] ];


//       3. sort by z
         for (long y=0; y<NParSameX-1; y++)
         {
            NParSameY = 1;    // 1 --> itself

//          find all particles with the same x and y coordinates
            while ( y+NParSameY<NParSameX  &&  PosY_Sorted[y] == PosY_Sorted[y+NParSameY] )  NParSameY ++;

            if ( NParSameY > 1 )
            {
               long *SortByZ_IdxTable = new long [NParSameY];
               real *PosZ_Sorted      = new real [NParSameY];

               for (long z=0; z<NParSameY; z++)
               {
                  SortByY_IdxTable[z] = IdxTable[ x + y + z ];
                  PosZ_Sorted     [z] = PosZ[ SortByY_IdxTable[z] ];
               }

               Mis_Heapsort( NParSameY, PosZ_Sorted, SortByZ_IdxTable );

               for (long z=0; z<NParSameY; z++)    IdxTable[ x + y + z ] = SortByY_IdxTable[ SortByZ_IdxTable[z] ];


//             4. sort by vx
               for (long z=0; z<NParSameY-1; z++)
               {
                  NParSameZ = 1;    // 1 --> itself

//                find all particles with the same x, y, and z coordinates
                  while ( z+NParSameZ<NParSameY  &&  PosZ_Sorted[z] == PosZ_Sorted[z+NParSameZ] )  NParSameZ ++;

                  if ( NParSameZ > 1 )
                  {
                     long *SortByVx_IdxTable = new long [NParSameZ];
                     real *VelX_Sorted       = new real [NParSameZ];

                     for (long w=0; w<NParSameZ; w++)
                     {
                        SortByZ_IdxTable[w] = IdxTable[ x + y + z + w ];
                        VelX_Sorted     [w] = VelX[ SortByZ_IdxTable[w] ];
                     }

                     Mis_Heapsort( NParSameZ, VelX_Sorted, SortByVx_IdxTable );

                     for (long w=0; w<NParSameZ; w++)    IdxTable[ x + y + z + w ] = SortByZ_IdxTable[ SortByVx_IdxTable[w] ];

                     delete [] SortByVx_IdxTable;
                     delete [] VelX_Sorted;

                     z += NParSameZ - 1;
                  } // if ( NParSameZ > 1 )
               } // for (long z=0; z<NParSameY-1; z++)

               delete [] SortByZ_IdxTable;
               delete [] PosZ_Sorted;

               y += NParSameY - 1;
            } // if ( NParSameY > 1 )
         } // for (long y=0; y<NParSameX-1; y++)

         delete [] SortByX_IdxTable;
         delete [] SortByY_IdxTable;
         delete [] PosY_Sorted;

         x += NParSameX - 1;
      } // if ( NParSameX > 1 )
   } // for (long x=0; x<NPar-1; x++)

   delete [] PosX_Sorted;

} // FUNCTION : SortParticle



//-------------------------------------------------------------------------------------------------------
// Function    :  LoadData
// Description :  Load the input data from the file "FileName_In"
//
// Parameter   :  patch       : Targeted GAMER pointer
//                FileName_In : Name of the input file
//                WithPot     : true --> the loaded data contain potential field
//                WithParDens : > 0  --> the loaded data contain particle density on grids
//                WithPar     : true --> the loaded data contain particles
//                NParVarOut  : Number of particle attributes stored in the file
//                NPar        : NUmber of particles
//                ParData     : Particle data array (allocated here --> must be dallocated manually later)
//-------------------------------------------------------------------------------------------------------
void LoadData( GAMER_t &patch, const char *FileName_In, bool &WithPot, int &WithParDens, bool &WithPar, int &NParVarOut,
               long &NPar, real **&ParData )
{

   fprintf( stdout, "Loading data %s ...\n", FileName_In );


// check if the target file exists
   FILE *File = fopen( FileName_In, "rb" );

   if ( File == NULL )  Aux_Error( ERROR_INFO, "input file \"%s\" does not exist !!\n", FileName_In );


// record the size of different data types
   const int size_bool   = sizeof( bool   );
   const int size_int    = sizeof( int    );
   const int size_uint   = sizeof( uint   );
   const int size_long   = sizeof( long   );
   const int size_ulong  = sizeof( ulong  );
   const int size_real   = sizeof( real   );
   const int size_double = sizeof( double );

   if ( size_int != size_uint )
      Aux_Error( ERROR_INFO, "sizeof(int) = %d != sizeof(uint) = %d !!\n", size_int, size_uint );

   if ( size_long != size_ulong )
      Aux_Error( ERROR_INFO, "sizeof(long) = %d != sizeof(ulong) = %d !!\n", size_long, size_ulong );


// a. load the information of data format
// =================================================================================================
   long FormatVersion, CheckCode;
   long HeaderSize_Format, HeaderSize_Makefile, HeaderSize_Constant, HeaderSize_Parameter, HeaderSize_SimuInfo;
   long HeaderOffset_Format, HeaderOffset_Makefile, HeaderOffset_Constant, HeaderOffset_Parameter, HeaderOffset_SimuInfo;
   long HeaderSize_Total;


// load the file format
   fread( &FormatVersion, sizeof(long), 1, File );


// load the size of each header (in bytes)
   fread( &HeaderSize_Format,    sizeof(long), 1, File );
   fread( &HeaderSize_Makefile,  sizeof(long), 1, File );
   fread( &HeaderSize_Constant,  sizeof(long), 1, File );
   fread( &HeaderSize_Parameter, sizeof(long), 1, File );
   fread( &HeaderSize_SimuInfo,  sizeof(long), 1, File );
   fread( &CheckCode,            sizeof(long), 1, File );


// determine the offset of each header (in bytes)
   HeaderOffset_Format    = 0;    // it must be zero
   HeaderOffset_Makefile  = HeaderOffset_Format    + HeaderSize_Format;
   HeaderOffset_Constant  = HeaderOffset_Makefile  + HeaderSize_Makefile;
   HeaderOffset_Parameter = HeaderOffset_Constant  + HeaderSize_Constant;
   HeaderOffset_SimuInfo  = HeaderOffset_Parameter + HeaderSize_Parameter;

   HeaderSize_Total       = HeaderOffset_SimuInfo  + HeaderSize_SimuInfo;


// verify the input data format version
   fprintf( stdout, "   The format version of the input file = %ld\n", FormatVersion );

   if ( FormatVersion < 2000 )
   {
      fprintf( stderr, "ERROR : unsupported data format version (only support version >= 2000) !!\n" );
      exit( EXIT_FAILURE );
   }


// check if the size of different data types are consistent
   int size_bool_restart, size_int_restart, size_long_restart, size_real_restart, size_double_restart;

   fread( &size_bool_restart,   sizeof(int), 1, File );
   fread( &size_int_restart,    sizeof(int), 1, File );
   fread( &size_long_restart,   sizeof(int), 1, File );
   fread( &size_real_restart,   sizeof(int), 1, File );
   fread( &size_double_restart, sizeof(int), 1, File );

   if ( size_bool_restart != size_bool )
      Aux_Error( ERROR_INFO, "sizeof(bool) is inconsistent : RESTART file = %d, runtime = %d !!\n",
                 size_bool_restart, size_bool );

   if ( size_int_restart != size_int )
      Aux_Error( ERROR_INFO, "sizeof(int) is inconsistent : RESTART file = %d, runtime = %d !!\n",
                 size_int_restart, size_int );

   if ( size_long_restart != size_long )
      Aux_Error( ERROR_INFO, "sizeof(long) is inconsistent : RESTART file = %d, runtime = %d !!\n",
                 size_long_restart, size_long );

   if ( size_real_restart != size_real )
      Aux_Error( ERROR_INFO, "sizeof(real) is inconsistent : RESTART file = %d, runtime = %d !!\n",
                 size_real_restart, size_real );

   if ( size_double_restart != size_double )
      Aux_Error( ERROR_INFO, "sizeof(double) is inconsistent : RESTART file = %d, runtime = %d !!\n",
                 size_double_restart, size_double );


// b. load all simulation parameters
// =================================================================================================
   Load_Parameter_After_2000( File, FormatVersion, WithPot, WithParDens, HeaderOffset_Makefile, HeaderOffset_Constant,
                              HeaderOffset_Parameter, patch.nx0_tot, patch.ngpu_x, WithPar, NParVarOut );


// c. load the simulation information
// =================================================================================================
   fprintf( stdout, "   Loading simulation information ... \n" );

// verify the check code
   long checkcode;

   fseek( File, HeaderOffset_SimuInfo, SEEK_SET );

   fread( &checkcode, sizeof(long), 1, File );

   if ( checkcode != CheckCode )
      Aux_Error( ERROR_INFO, "incorrect check code in the RESTART file (input %ld <-> expeect %ld) !!\n",
                 checkcode, CheckCode );


// load information necessary for restart
   long   AdvanceCounter[NLEVEL];
   int    NPatchTotal[NLEVEL], NDataPatch_Total[NLEVEL], DumpID;
   long   Step, FileOffset_Particle;
   double Time[NLEVEL];

   fread( &DumpID,          sizeof(int),         1, File );
   fread( Time,             sizeof(double), NLEVEL, File );
   fread( &Step,            sizeof(long),        1, File );
   fread( NPatchTotal,      sizeof(int),    NLEVEL, File );
   fread( NDataPatch_Total, sizeof(int),    NLEVEL, File );
   fread( AdvanceCounter,   sizeof(long),   NLEVEL, File );
   fseek( File, sizeof(double), SEEK_CUR );

   if ( WithPar ) {
   fread( &NPar,                sizeof(long),    1, File );
   fread( &FileOffset_Particle, sizeof(long),    1, File ); }


// verify the size of the RESTART file
   fprintf( stdout, "      Verifying the size of the RESTART file ...\n" );

   long ExpectSize, InputSize, PatchDataSize, DataSize[NLEVEL];
   int  NGridVar = NCOMP_TOTAL;  // number of grid variables

   if ( WithPot )       NGridVar ++;
   if ( WithParDens )   NGridVar ++;

   PatchDataSize = CUBE(PATCH_SIZE)*NGridVar*sizeof(real);
   ExpectSize    = HeaderSize_Total;

   for (int lv=0; lv<NLEVEL; lv++)
   {
      DataSize[lv]  = 0;
      DataSize[lv] += NPatchTotal[lv]*4*sizeof(int);        // 4 = corner(3) + son(1)
      DataSize[lv] += NDataPatch_Total[lv]*PatchDataSize;

      ExpectSize   += DataSize[lv];
   }

   if ( WithPar )
   {
      for (int lv=0; lv<NLEVEL; lv++)
      {
//       2 = NPar + starting particle index stored in each leaf patch
         const long ParInfoSize = (long)NDataPatch_Total[lv]*2*sizeof(long);

         DataSize[lv] += ParInfoSize;
         ExpectSize   += ParInfoSize;
      }

      ExpectSize += (long)NParVarOut*NPar*sizeof(real);
   }

   fseek( File, 0, SEEK_END );
   InputSize = ftell( File );

   if ( InputSize != ExpectSize )
      Aux_Error( ERROR_INFO, "size of the file <%s> is incorrect --> input = %ld <-> expect = %ld !!\n",
                 FileName_In, InputSize, ExpectSize );

   fprintf( stdout, "      Verifying the size of the RESTART file ... passed\n" );
   fprintf( stdout, "   Loading simulation information ... done\n" );


// d. load the simulation data
// =================================================================================================
   int  LoadCorner[3], LoadSon, PID;

   fseek( File, HeaderSize_Total, SEEK_SET );

   for (int lv=0; lv<NLEVEL; lv++)
   {
      fprintf( stdout, "   Loading grid data at level %2d ... ", lv );

      for (int LoadPID=0; LoadPID<NPatchTotal[lv]; LoadPID++)
      {
//       d1. load the patch information
         fread(  LoadCorner, sizeof(int), 3, File );
         fread( &LoadSon,    sizeof(int), 1, File );

//       d2. create the patch and load the physical data only if it is a leaf patch
         if ( LoadSon == -1 )
         {
//          skip particle info
            if ( WithPar ) fseek( File, 2*sizeof(long), SEEK_CUR );

            PID = patch.num[lv];

            patch.pnew( lv, LoadCorner[0], LoadCorner[1], LoadCorner[2], -1, true );
            patch.ptr[lv][PID]->son = -1;    // set the SonPID as -1 to indicate that it's a leaf patch

//          d2-1. load the fluid variables
            fread( patch.ptr[lv][PID]->fluid,    sizeof(real), PATCH_SIZE*PATCH_SIZE*PATCH_SIZE*NCOMP_TOTAL, File );

//          d2-2. load the gravitational potential
            if ( WithPot )
            fread( patch.ptr[lv][PID]->pot,      sizeof(real), PATCH_SIZE*PATCH_SIZE*PATCH_SIZE,             File );

//          d2-3. load the particle density on grids
            if ( WithParDens )
            fread( patch.ptr[lv][PID]->par_dens, sizeof(real), PATCH_SIZE*PATCH_SIZE*PATCH_SIZE,             File );
         }
      } // for (int LoadPID=0; LoadPID<NPatchTotal[lv]; LoadPID++)

      fprintf( stdout, "done\n" );
   } // for (int lv=0; lv<NLEVEL; lv++)

   printf( "   DumpID      = %d\n",  DumpID      );
   printf( "   Step        = %ld\n", Step        );
   printf( "   Time        = %lf\n", Time[0]     );
   printf( "   WithPot     = %d\n",  WithPot     );
   printf( "   WithParDens = %d\n",  WithParDens );
   printf( "   WithPar     = %d\n",  WithPar     );
   if ( WithPar ) {
   printf( "   NParVarOut  = %d\n",  NParVarOut  );
   printf( "   NPar        = %ld\n", NPar        ); }
   for (int lv=0; lv<NLEVEL; lv++)
   printf( "   NPatch[%2d] = %d\n",  lv, patch.num[lv] );

   fclose( File );


// load particles
   if ( WithPar )
   {
      const long ParDataSize1v = NPar*sizeof(real);

      Aux_AllocateArray2D( ParData, NParVarOut, NPar );

      File = fopen( FileName_In, "rb" );

//    one must not call fread when NPar == 0, for which ParData[0...NParVarOut-1] will be ill-defined
      if ( NPar > 0 )
      for (int v=0; v<NParVarOut; v++)
      {
         fseek( File, FileOffset_Particle + v*ParDataSize1v, SEEK_SET );

         fread( ParData[v], sizeof(real), NPar, File );
      }

      fclose( File );
   } // if ( WithPar )


   fprintf( stdout, "Loading data %s ... done\n", FileName_In );

} // FUNCTION : LoadData



//-------------------------------------------------------------------------------------------------------
// Function    :  Load_Parameter_After_2000
// Description :  Load all simulation parameters from the RESTART file with format version >= 2000
//
// Note        :  All floating-point variables are declared as double after version 2000
//
// Parameter   :  File           : RESTART file pointer
//                FormatVersion  : Format version of the RESTART file
//                WithPot        : Whether or not the RESTART file stores the potential data
//                WithParDens    : Whether or not the RESTART file stores the particle density data (deposited onto grids)
//                HeaderOffset_X : Offsets of different headers
//                NX0_Tot        : Total number of base-level cells along each direction
//                NGPU_X         : Number of MPI processes along each direction
//                WithPar        : Whether or not the RESTART file stores particles
//                NParVarOut     : Number of particle attributes stored in the file
//
// Return      :  LoadPot (END_T and END_STEP may also be set to the original values)
//-------------------------------------------------------------------------------------------------------
void Load_Parameter_After_2000( FILE *File, const int FormatVersion, bool &WithPot, int &WithParDens,
                                const long HeaderOffset_Makefile, const long HeaderOffset_Constant,
                                const long HeaderOffset_Parameter, int *NX0_Tot, int *NGPU_X,
                                bool &WithPar, int &NParVarOut )
{

   fprintf( stdout, "   Loading simulation parameters ...\n" );


// a. load the simulation options and parameters defined in the Makefile
// =================================================================================================
   bool gravity, individual_timestep, comoving, gpu, gamer_optimization, gamer_debug, timing, timing_solver;
   bool intel, float8, serial, overlap_mpi, openmp;
   int  model, pot_scheme, flu_scheme, lr_scheme, rsolver, load_balance, nlevel, max_patch, gpu_arch, ncomp_passive;
   bool store_pot_ghost, unsplit_gravity, particle;

   fseek( File, HeaderOffset_Makefile, SEEK_SET );

   fread( &model,                      sizeof(int),                     1,             File );
   fread( &gravity,                    sizeof(bool),                    1,             File );
   fread( &pot_scheme,                 sizeof(int),                     1,             File );
   fread( &individual_timestep,        sizeof(bool),                    1,             File );
   fread( &comoving,                   sizeof(bool),                    1,             File );
   fread( &flu_scheme,                 sizeof(int),                     1,             File );
   fread( &lr_scheme,                  sizeof(int),                     1,             File );
   fread( &rsolver,                    sizeof(int),                     1,             File );
   fread( &gpu,                        sizeof(bool),                    1,             File );
   fread( &gamer_optimization,         sizeof(bool),                    1,             File );
   fread( &gamer_debug,                sizeof(bool),                    1,             File );
   fread( &timing,                     sizeof(bool),                    1,             File );
   fread( &timing_solver,              sizeof(bool),                    1,             File );
   fread( &intel,                      sizeof(bool),                    1,             File );
   fread( &float8,                     sizeof(bool),                    1,             File );
   fread( &serial,                     sizeof(bool),                    1,             File );
   fread( &load_balance,               sizeof(int),                     1,             File );
   fread( &overlap_mpi,                sizeof(bool),                    1,             File );
   fread( &openmp,                     sizeof(bool),                    1,             File );
   fread( &gpu_arch,                   sizeof(int),                     1,             File );
   fread( &nlevel,                     sizeof(int),                     1,             File );
   fread( &max_patch,                  sizeof(int),                     1,             File );
   fread( &store_pot_ghost,            sizeof(bool),                    1,             File );
   fread( &unsplit_gravity,            sizeof(bool),                    1,             File );
   if ( FormatVersion >= 2100 )
   fread( &particle,                   sizeof(bool),                    1,             File );
   else
   particle = false;
   if ( FormatVersion >= 2120 )
   fread( &ncomp_passive,              sizeof(int),                     1,             File );
   else
   ncomp_passive = 0;


// b. load the symbolic constants defined in "Macro.h, CUPOT.h, and CUFLU.h"
// =================================================================================================
   bool   enforce_positive, char_reconstruction, hll_no_ref_state, hll_include_all_waves, waf_dissipate;
   bool   use_psolver_10to14;
   int    ncomp_fluid, patch_size, flu_ghost_size, pot_ghost_size, gra_ghost_size, check_intermediate;
   int    flu_block_size_x, flu_block_size_y, pot_block_size_x, pot_block_size_z, gra_block_size_z;
   int    par_nvar, par_npassive;
   double min_value, max_error;

   fseek( File, HeaderOffset_Constant, SEEK_SET );

   fread( &ncomp_fluid,                sizeof(int),                     1,             File );
   fread( &patch_size,                 sizeof(int),                     1,             File );
   fread( &min_value,                  sizeof(double),                  1,             File );
   fread( &flu_ghost_size,             sizeof(int),                     1,             File );
   fread( &pot_ghost_size,             sizeof(int),                     1,             File );
   fread( &gra_ghost_size,             sizeof(int),                     1,             File );
   fread( &enforce_positive,           sizeof(bool),                    1,             File );
   fread( &char_reconstruction,        sizeof(bool),                    1,             File );
   fread( &check_intermediate,         sizeof(int),                     1,             File );
   fread( &hll_no_ref_state,           sizeof(bool),                    1,             File );
   fread( &hll_include_all_waves,      sizeof(bool),                    1,             File );
   fread( &waf_dissipate,              sizeof(bool),                    1,             File );
   fread( &max_error,                  sizeof(double),                  1,             File );
   fread( &flu_block_size_x,           sizeof(int),                     1,             File );
   fread( &flu_block_size_y,           sizeof(int),                     1,             File );
   fread( &use_psolver_10to14,         sizeof(bool),                    1,             File );
   fread( &pot_block_size_x,           sizeof(int),                     1,             File );
   fread( &pot_block_size_z,           sizeof(int),                     1,             File );
   fread( &gra_block_size_z,           sizeof(int),                     1,             File );
   if ( particle ) {
   fread( &par_nvar,                   sizeof(int),                     1,             File );
   fread( &par_npassive,               sizeof(int),                     1,             File ); }



// c. load the simulation parameters recorded in the file "Input__Parameter"
// =================================================================================================
   bool   opt__adaptive_dt, opt__dt_user, opt__flag_rho, opt__flag_rho_gradient, opt__flag_pres_gradient;
   bool   opt__flag_engy_density, opt__flag_user, opt__fixup_flux, opt__fixup_restrict, opt__overlap_mpi;
   bool   opt__gra_p5_gradient, opt__int_time, opt__output_test_error, opt__output_base, opt__output_pot;
   bool   opt__output_baseps, opt__timing_balance, opt__int_phase, opt__corr_unphy;
   int    nx0_tot[3], mpi_nrank, mpi_nrank_x[3], omp_nthread, regrid_count, opt__output_par_dens;
   int    flag_buffer_size, max_level, opt__lr_limiter, opt__waf_limiter, flu_gpu_npgroup, gpu_nstream;
   int    sor_max_iter, sor_min_iter, mg_max_iter, mg_npre_smooth, mg_npost_smooth, pot_gpu_npgroup;
   int    opt__flu_int_scheme, opt__pot_int_scheme, opt__rho_int_scheme;
   int    opt__gra_int_scheme, opt__ref_flu_int_scheme, opt__ref_pot_int_scheme;
   int    opt__output_total, opt__output_part, opt__output_mode, output_step, opt__corr_unphy_scheme;
   long   end_step;
   double lb_wli_max, gamma, minmod_coeff, ep_coeff, elbdm_mass, elbdm_planck_const, newton_g, sor_omega;
   double mg_tolerated_error, output_part_x, output_part_y, output_part_z;
   double box_size, end_t, omega_m0, dt__fluid, dt__gravity, dt__phase, dt__max_delta_a, output_dt;

   fseek( File, HeaderOffset_Parameter, SEEK_SET );

   fread( &box_size,                   sizeof(double),                  1,             File );
   fread(  nx0_tot,                    sizeof(int),                     3,             File );
   fread( &mpi_nrank,                  sizeof(int),                     1,             File );
   fread(  mpi_nrank_x,                sizeof(int),                     3,             File );
   fread( &omp_nthread,                sizeof(int),                     1,             File );
   fread( &end_t,                      sizeof(double),                  1,             File );
   fread( &end_step,                   sizeof(long),                    1,             File );
   fread( &omega_m0,                   sizeof(double),                  1,             File );
   fread( &dt__fluid,                  sizeof(double),                  1,             File );
   fread( &dt__gravity,                sizeof(double),                  1,             File );
   fread( &dt__phase,                  sizeof(double),                  1,             File );
   fread( &dt__max_delta_a,            sizeof(double),                  1,             File );
   fread( &opt__adaptive_dt,           sizeof(bool),                    1,             File );
   fread( &opt__dt_user,               sizeof(bool),                    1,             File );
   fread( &regrid_count,               sizeof(int),                     1,             File );
   fread( &flag_buffer_size,           sizeof(int),                     1,             File );
   fread( &max_level,                  sizeof(int),                     1,             File );
   fread( &opt__flag_rho,              sizeof(bool),                    1,             File );
   fread( &opt__flag_rho_gradient,     sizeof(bool),                    1,             File );
   fread( &opt__flag_pres_gradient,    sizeof(bool),                    1,             File );
   fread( &opt__flag_engy_density,     sizeof(bool),                    1,             File );
   fread( &opt__flag_user,             sizeof(bool),                    1,             File );
   fread( &lb_wli_max,                 sizeof(double),                  1,             File );
   fread( &gamma,                      sizeof(double),                  1,             File );
   fread( &minmod_coeff,               sizeof(double),                  1,             File );
   fread( &ep_coeff,                   sizeof(double),                  1,             File );
   fread( &opt__lr_limiter,            sizeof(int),                     1,             File );
   fread( &opt__waf_limiter,           sizeof(int),                     1,             File );
   fread( &elbdm_mass,                 sizeof(double),                  1,             File );
   fread( &elbdm_planck_const,         sizeof(double),                  1,             File );
   fread( &flu_gpu_npgroup,            sizeof(int),                     1,             File );
   fread( &gpu_nstream,                sizeof(int),                     1,             File );
   fread( &opt__fixup_flux,            sizeof(bool),                    1,             File );
   fread( &opt__fixup_restrict,        sizeof(bool),                    1,             File );
   fread( &opt__overlap_mpi,           sizeof(bool),                    1,             File );
   fread( &newton_g,                   sizeof(double),                  1,             File );
   fread( &sor_omega,                  sizeof(double),                  1,             File );
   fread( &sor_max_iter,               sizeof(int),                     1,             File );
   fread( &sor_min_iter,               sizeof(int),                     1,             File );
   fread( &mg_max_iter,                sizeof(int),                     1,             File );
   fread( &mg_npre_smooth,             sizeof(int),                     1,             File );
   fread( &mg_npost_smooth,            sizeof(int),                     1,             File );
   fread( &mg_tolerated_error,         sizeof(double),                  1,             File );
   fread( &pot_gpu_npgroup,            sizeof(int),                     1,             File );
   fread( &opt__gra_p5_gradient,       sizeof(bool),                    1,             File );
   fread( &opt__int_time,              sizeof(bool),                    1,             File );
   fread( &opt__int_phase,             sizeof(bool),                    1,             File );
   fread( &opt__flu_int_scheme,        sizeof(int),                     1,             File );
   fread( &opt__pot_int_scheme,        sizeof(int),                     1,             File );
   fread( &opt__rho_int_scheme,        sizeof(int),                     1,             File );
   fread( &opt__gra_int_scheme,        sizeof(int),                     1,             File );
   fread( &opt__ref_flu_int_scheme,    sizeof(int),                     1,             File );
   fread( &opt__ref_pot_int_scheme,    sizeof(int),                     1,             File );
   fread( &opt__output_total,          sizeof(int),                     1,             File );
   fread( &opt__output_part,           sizeof(int),                     1,             File );
   fread( &opt__output_test_error,     sizeof(bool),                    1,             File );
   fread( &opt__output_base,           sizeof(bool),                    1,             File );
   fread( &opt__output_pot,            sizeof(bool),                    1,             File );
   fread( &opt__output_mode,           sizeof(int),                     1,             File );
   fread( &output_step,                sizeof(int),                     1,             File );
   fread( &output_dt,                  sizeof(double),                  1,             File );
   fread( &output_part_x,              sizeof(double),                  1,             File );
   fread( &output_part_y,              sizeof(double),                  1,             File );
   fread( &output_part_z,              sizeof(double),                  1,             File );
   fread( &opt__timing_balance,        sizeof(bool),                    1,             File );
   fread( &opt__output_baseps,         sizeof(bool),                    1,             File );
   fread( &opt__corr_unphy,            sizeof(bool),                    1,             File );
   fread( &opt__corr_unphy_scheme,     sizeof(int),                     1,             File );
   if ( particle )
   fread( &opt__output_par_dens,       sizeof(int),                     1,             File );
   else
   opt__output_par_dens = 0;



   fprintf( stdout, "   Loading simulation parameters ... done\n" );


// d. check parameters (before loading any size-dependent parameters)
// =================================================================================================
   fprintf( stdout, "   Checking loaded parameters ...\n" );


   const bool Fatal    = true;
   const bool NonFatal = false;

#  ifdef FLOAT8
   if ( !float8 )
   {
      fprintf( stderr, "ERROR : %s : RESTART file (%s) != runtime (%s) !!\n", "FLOAT8", "OFF", "ON" );
      exit( 1 );
   }
#  else
   if (  float8 )
   {
      fprintf( stderr, "ERROR : %s : RESTART file (%s) != runtime (%s) !!\n", "FLOAT8", "ON", "OFF" );
      exit( 1 );
   }
#  endif

   CompareVar( "MODEL",                   model,                  MODEL,                        Fatal );
   CompareVar( "NLEVEL",                  nlevel,                 NLEVEL,                       Fatal );
   CompareVar( "NCOMP_FLUID",             ncomp_fluid,            NCOMP_FLUID,                  Fatal );
   CompareVar( "NCOMP_PASSIVE",           ncomp_passive,          NCOMP_PASSIVE,                Fatal );
   CompareVar( "PATCH_SIZE",              patch_size,             PATCH_SIZE,                   Fatal );


   fprintf( stdout, "   Checking loaded parameters ... done\n" );


// set the returned variables
   WithPot = opt__output_pot;
   for (int d=0; d<3; d++)
   {
      NX0_Tot[d]  = nx0_tot    [d];
      NGPU_X [d]  = mpi_nrank_x[d];
   }

   WithParDens = opt__output_par_dens;
   WithPar     = particle;
   if ( WithPar )
   {
      if ( FormatVersion > 2131 )   NParVarOut = par_nvar;           // after version 2131, par_nvar = PAR_NATT_STORED
      else                          NParVarOut = 7 + par_npassive;   // mass, position x/y/z, velocity x/y/z, and passive variables
   }
   else
                                    NParVarOut = -1;

} // FUNCTION : Load_Parameter_After_2000



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Error
// Description :  Output the error messages and force the program to be terminated
//
// Note        :  Use the variable argument lists provided in "cstdarg"
//
// Parameter   :  File     : Name of the file where error occurs
//                Line     : Line number where error occurs
//                Func     : Name of the function where error occurs
//                Format   : Output format
//                ...      : Arguments in vfprintf
//-------------------------------------------------------------------------------------------------------
void Aux_Error( const char *File, const int Line, const char *Func, const char *Format, ... )
{

// flush all previous messages
   fflush( stdout ); fflush( stdout ); fflush( stdout );
   fflush( stderr ); fflush( stderr ); fflush( stderr );


// output error messages
   va_list Arg;
   va_start( Arg, Format );

   fprintf ( stderr, "********************************************************************************\n" );
   fprintf ( stderr, "ERROR : " );
   vfprintf    ( stderr, Format, Arg );
   fprintf ( stderr, "        file <%s>, line <%d>, function <%s>\n", File, Line, Func );
   fprintf ( stderr, "********************************************************************************\n" );

   va_end( Arg );


// terminate the program
   exit( EXIT_FAILURE );

} // FUNCTION : Aux_Error



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_AllocateArray2D
// Description :  Allocate a continuous 2D array with size [J][I]
//
// Note        :  1. Call-by-reference
//                2. Free memory by Aux_DeallocateArray2D
//
// Parameter   :  Array : Pointer to be allocated
//                J/I   : Array dimensions
//-------------------------------------------------------------------------------------------------------
void Aux_AllocateArray2D( real** &Array, const int J, const int I )
{

   if ( J < 0  ||  I < 0 )    Aux_Error( ERROR_INFO, "incorrect array size (J = %d, I = %d) !!\n", J, I );

   if ( J == 0  ||  I == 0 )
   {
      Array = NULL;
      return;
   }

   Array    = new real* [J  ];
   Array[0] = new real  [J*I];

   for (int j=1; j<J; j++)    Array[j] = Array[j-1] + I;

} // FUNCTION : Aux_AllocateArray2D



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_DeallocateArray2D
// Description :  Free an array previously allocated by Aux_AllocateArray2D
//
// Note        :  1. Call-by-reference
//                2. Pointer is reset to NULL
//
// Parameter   :  Array : Pointer to be deallocated
//-------------------------------------------------------------------------------------------------------
void Aux_DeallocateArray2D( real** &Array )
{

   if ( Array == NULL )    return;

   delete [] Array[0];
   delete [] Array;

   Array = NULL;

} // FUNCTION : Aux_DeallocateArray2D



//-------------------------------------------------------------------------------------------------------
// Function    :  CompareVar
// Description :  Compare the input variables
//
// Note        :  This function is overloaded to work with different data types
//
// Parameter   :  VarName     : Name of the targeted variable
//                RestartVar  : Variable loaded from the RESTART file
//                RuntimeVar  : Variable loaded from the Input__Parameter
//                Fatal       : Whether or not the difference between RestartVar and RuntimeVar is fatal
//                              --> true  : terminate the program if the input variables are different
//                                  false : display warning message if the input variables are different
//-------------------------------------------------------------------------------------------------------
void CompareVar( const char *VarName, const int RestartVar, const int RuntimeVar, const bool Fatal )
{

   if ( RestartVar != RuntimeVar )
   {
      if ( Fatal )
      {
         fprintf( stderr, "ERROR : %s : RESTART file (%d) != runtime (%d) !!\n",
                  VarName, RestartVar, RuntimeVar );
         exit( 1 );
      }
      else
         fprintf( stderr, "WARNING : %s : RESTART file (%d) != runtime (%d) !!\n",
                  VarName, RestartVar, RuntimeVar );
   }

} // FUNCTION : CompareVar (int)



//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_Heapsort
// Description :  Use the Heapsort algorithm to sort the input array into ascending numerical order
//                --> An index table will also be constructed if "IdxTable != NULL"
//
// Note        :  1. Ref : Numerical Recipes Chapter 8.3 - 8.4
//
// Parameter   :  N        :  Size of Array
//                Array    :  Array to be sorted into ascending numerical order
//                IdxTable :  Index table
//-------------------------------------------------------------------------------------------------------
void Mis_Heapsort( const long N, real Array[], long IdxTable[] )
{

// initialize the IdxTable
   if ( IdxTable != NULL )
      for (long t=0; t<N; t++)   IdxTable[t] = t;

// heap creation
   for (long L=N/2-1; L>=0; L--)    Heapsort_SiftDown( L, N-1, Array, IdxTable );

// retirement-and-promotion
   real Buf;
   for (long R=N-1; R>0; R--)
   {
      Buf      = Array[R];
      Array[R] = Array[0];
      Array[0] = Buf;

      if ( IdxTable != NULL )
      {
         Buf         = IdxTable[R];
         IdxTable[R] = IdxTable[0];
         IdxTable[0] = Buf;
      }

      Heapsort_SiftDown( 0, R-1, Array, IdxTable );
   }

} // FUNCTION : Mis_Heapsort



//-------------------------------------------------------------------------------------------------------
// Function    :  Heapsort_SiftDown
// Description :  Sift-down process for the Heapsort algorithm
//
// Note        :  1. Ref : Numerical Recipes Chapter 8.3 - 8.4
//
// Parameter   :  L        :  Left  range of the sift-down
//                R        :  Right range of the sift-down
//                Array    :  Array to be sorted into ascending numerical order
//                IdxTable :  Index table
//-------------------------------------------------------------------------------------------------------
void Heapsort_SiftDown( const long L, const long R, real Array[], long IdxTable[] )
{

   long Idx_up    = L;
   long Idx_down  = 2*Idx_up + 1;
   real Target    = Array[Idx_up];
   long TargetIdx = ( IdxTable == NULL ) ? -1 : IdxTable[Idx_up];

   while ( Idx_down <= R )
   {
//    find the better employee
      if ( Idx_down < R  &&  Array[Idx_down+1] > Array[Idx_down] )   Idx_down ++;

//    terminate the sift-down process if the target (supervisor) is better than both its employees
      if ( Target >= Array[Idx_down] )    break;

//    otherwise, promote the better employee
      Array[Idx_up] = Array[Idx_down];
      if ( IdxTable != NULL )    IdxTable[Idx_up] = IdxTable[Idx_down];

//    prepare the next sift-down operation
      Idx_up   = Idx_down;
      Idx_down = 2*Idx_up + 1;
   }

// put target at its best position
   Array[Idx_up] = Target;
   if ( IdxTable != NULL )    IdxTable[Idx_up] = TargetIdx;

} // FUNCTION : Heapsort_SiftDown



//-------------------------------------------------------------------------------------------------------
// Function    :  main
// Description :
//-------------------------------------------------------------------------------------------------------
int main( int argc, char ** argv )
{

   ReadOption( argc, argv );

   CheckParameter();

   LoadData( patch1, FileName_In1, WithPot1, WithParDens1, WithPar1, NParVarOut1, NPar1, ParData1 );
   LoadData( patch2, FileName_In2, WithPot2, WithParDens2, WithPar2, NParVarOut2, NPar2, ParData2 );

   CompareGridData();
   CompareParticleData();

   Aux_DeallocateArray2D( ParData1 );
   Aux_DeallocateArray2D( ParData2 );

   cout << "Program terminated successfully" << endl;

   return 0;

} // FUNCTION : main
