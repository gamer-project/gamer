#include "ExtractProfile.h"



// general-purpose variables
AMR_t       amr;
int         DumpID, NX0_TOT[3];
long int    Step;
double      Time[NLEVEL];
double      GAMMA;
double      Center_Map[3];             // the coordinates of the "mapped" sphere center for the periodic B.C.
int        *BaseP             = NULL;

char       *FileName_In       = NULL;
char       *FileName_Tree     = NULL;
tree_t     *tree              = NULL;
char       *Suffix            = NULL;
double      Center[3]         = { WRONG, WRONG, WRONG };    // the coordinates of the sphere center
double      MaxRadius         = WRONG;    // the maximum radius of the sphere
double      UseMaxRhoPos_R    = WRONG;    // maximum radius used by UseMaxRhoPos
bool        Periodic          = false;    // periodic B.C. (not the case by default)
bool        OutputPot         = false;    // output potential as well
int         OutputParDens     = 0;        // output particle density on grids (0/1/2: off/particle density/total density)
bool        InputScale        = false;    // true --> input cell scales instead of coordinates
bool        UseTree           = false;    // true --> use the tree file to improve the I/O performance
bool        NeedGhost         = false;    // treu --> ghost zones are required (--> need to construct the octree structure)

// variables for the mode "shell average"
double      **Average         = NULL;     // the shell average in each shell ( HYDRO : rho, v_r, v_t, engy, pres, pot )
                                          //                                 ( ELBDM : rho, real, imag, pot, Ek )
double      **RMS             = NULL;     // the standard deviation (root-mean-squre) in each shell
real        **Max             = NULL;     // the maximum value in each shell
real        **Min             = NULL;     // the minimum value in each shell
long int    *NCount           = NULL;     // the total number of grids in each shell
double      *Volume           = NULL;     // the total volume in each shell
double      ShellWidth;                   // the width of each shell (or minimum shell widith if LogBin is enabled)

bool        Mode_ShellAve     = false;    // true --> evaluate the shell average
int         NShell            = WRONG;    // number of shells
int         UseMaxRhoPos      = 8;        // number of highest-density cells for setting the sphere center (<=0 -> disable)
bool        UseMaxRhoPos_Par  = false;    // if OutputParDens is on, use particle (or total) density for UseMaxRhoPos
int         NIn               = NULL_INT; // total number of input grid variables
int         NOut              = NULL_INT; // total number of output variables
bool        OutputCenter      = false;    // output the center coords
double      INT_MONO_COEFF    = 2.0;
#if ( MODEL == ELBDM )
bool        ELBDM_IntPhase    = true;
bool        ELBDM_GetVir      = false;    // analyze the ELBDM virial condition (output vel, Ek, and virial surface terms)
double      ELBDM_ETA         = NULL_REAL;
IntScheme_t IntScheme         = INT_CQUAR;
double    (*ELBDM_Mom)[3]     = NULL;     // momentum used for subtracting the CM Ek in ELBDM_GetVir
double     *ELBDM_RhoUr2      = NULL;     // Rho*(vr^2 + wr^2) in the virial surface terms
double     *ELBDM_dRho_dr     = NULL;     // dRho/dr in the virial surface terms
double     *ELBDM_LapRho      = NULL;     // Lap(Rho) in the virial surface terms
#else
IntScheme_t IntScheme         = INT_CQUAD;
#endif
bool        GetAvePot         = false;    // calculate the spherically symmetric gravitational potential
double      NewtonG           = WRONG;    // gravitational constant (must be set if it cannot be determined from the input data)


// variables for the mode "maximum density"
bool        Mode_MaxRho       = false; // true --> get the cells exceeding the density threshold
double      RhoThres          = WRONG; // the density threshold (for the MaxRho mode)
double      LogBin            = 0.0;   // log scale radial bin (<=1.0 --> disable --> use linear scale)
double      GetNShell         = 2.0;   // set the minimum shell width to "GetNShell" times the cell size around center
                                       // (<=0 -> disable)





//-------------------------------------------------------------------------------------------------------
// Function    :  GetMaxRho
// Description :  Output the mesh information if the density exceeds a given threshold "RhoThres"
//-------------------------------------------------------------------------------------------------------
void GetMaxRho()
{

   cout << "GetMaxRho ..." << endl;


   const double dh_min = amr.dh[NLEVEL-1];
   double Radius, x, x1, x2, y, y1, y2, z, z1, z2, scale;   // (x,y,z) : relative coordinates to the vector "Center"
   double xx, yy, zz;
   real   Pass[NCOMP_PASSIVE];

#  if   ( MODEL == HYDRO )
   real Dens, VelX, VelY, VelZ, Engy, Pres, Pot, ParDens;

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   real Dens, Real, Imag, Pot, ParDens;

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL


// open the file
   char FileName[50] = "MaxRho";

   if ( Suffix != NULL )   strcat( FileName, Suffix );

   if ( NULL != fopen(FileName,"r") )
         fprintf( stderr, "Warning : the file \"%s\" already exists and will be overwritten !!\n", FileName );

   FILE *File = fopen( FileName, "w" );

#  if   ( MODEL == HYDRO )
   fprintf( File, "%10s %10s %10s %2s %12s %12s %12s %12s %12s %12s",
            "x", "y", "z", "Lv", "Dens", "VelX", "VelY", "VelZ", "Engy", "Pres" );
   for (int v=0; v<NCOMP_PASSIVE; v++)    fprintf( File, " %10s%02d", "Passive", v );
   if      ( OutputPot          )         fprintf( File, " %12s",     "Pot" );
   if      ( OutputParDens == 1 )         fprintf( File, " %12s",     "ParDens" );
   else if ( OutputParDens == 2 )         fprintf( File, " %12s",     "TotalDens" );
   fprintf( File, "\n" );

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   fprintf( File, "%10s %10s %10s %2s %12s %12s %12s",
            "x", "y", "z", "Lv", "Dens", "Real", "Imag" );
   for (int v=0; v<NCOMP_PASSIVE; v++)    fprintf( File, " %10s%02d", "Passive", v );
   if      ( OutputPot          )         fprintf( File, " %12s",     "Pot" );
   if      ( OutputParDens == 1 )         fprintf( File, " %12s",     "ParDens" );
   else if ( OutputParDens == 2 )         fprintf( File, " %12s",     "TotalDens" );
   fprintf( File, "\n" );

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL


   for (int lv=0; lv<NLEVEL; lv++)
   {
      scale = (double)amr.scale[lv];
      cout << "   Level " << lv << " ... ";

      for (int PID=0; PID<amr.num[lv]; PID++)
      {
         if ( amr.patch[lv][PID]->fluid == NULL  ||  amr.patch[lv][PID]->son != -1 )   continue;

         for (int k=0; k<PATCH_SIZE; k++) {  zz = amr.patch[lv][PID]->corner[2] + (k+0.5)*scale;
                                             z1 = zz - Center    [2];
                                             z2 = zz - Center_Map[2];
                                             z  = ( fabs(z1) <= fabs(z2) ) ?  z1 : z2;
         for (int j=0; j<PATCH_SIZE; j++) {  yy = amr.patch[lv][PID]->corner[1] + (j+0.5)*scale;
                                             y1 = yy - Center    [1];
                                             y2 = yy - Center_Map[1];
                                             y  = ( fabs(y1) <= fabs(y2) ) ?  y1 : y2;
         for (int i=0; i<PATCH_SIZE; i++) {  xx = amr.patch[lv][PID]->corner[0] + (i+0.5)*scale;
                                             x1 = xx - Center    [0];
                                             x2 = xx - Center_Map[0];
                                             x  = ( fabs(x1) <= fabs(x2) ) ?  x1 : x2;

            Radius = sqrt( x*x + y*y + z*z );

//          output the mesh information if it's within the targeted sphere and the density exceeds "RhoThres"
            if ( Radius < MaxRadius  &&  amr.patch[lv][PID]->fluid[DENS][k][j][i] >= RhoThres )
            {
#              if   ( MODEL == HYDRO )
               Dens    = amr.patch[lv][PID]->fluid[DENS][k][j][i];
               VelX    = amr.patch[lv][PID]->fluid[MOMX][k][j][i] / Dens;
               VelY    = amr.patch[lv][PID]->fluid[MOMY][k][j][i] / Dens;
               VelZ    = amr.patch[lv][PID]->fluid[MOMZ][k][j][i] / Dens;
               Engy    = amr.patch[lv][PID]->fluid[ENGY][k][j][i];

               for (int v=0; v<NCOMP_PASSIVE; v++)
               Pass[v] = amr.patch[lv][PID]->fluid[ NCOMP_FLUID + v ][k][j][i];

               if ( OutputPot )
               Pot     = amr.patch[lv][PID]->pot        [k][j][i];

               if ( OutputParDens )
               ParDens = amr.patch[lv][PID]->par_dens   [k][j][i];

               Pres = (GAMMA-1.0) * ( Engy - 0.5*Dens*(VelX*VelX + VelY*VelY + VelZ*VelZ) );


               fprintf( File, "%10.4e %10.4e %10.4e %2d %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e",
                        xx*dh_min, yy*dh_min, zz*dh_min, lv, Dens, VelX, VelY, VelZ, Engy, Pres );

               for (int v=0; v<NCOMP_PASSIVE; v++)    fprintf( File, " %12.5e", Pass[v] );
               if ( OutputPot     )                   fprintf( File, " %12.5e", Pot );
               if ( OutputParDens )                   fprintf( File, " %12.5e", ParDens );


               fprintf( File, "\n" );

#              elif ( MODEL == MHD )
#              warning : WAIT MHD !!!

#              elif ( MODEL == ELBDM )
               Dens    = amr.patch[lv][PID]->fluid[DENS][k][j][i];
               Real    = amr.patch[lv][PID]->fluid[REAL][k][j][i];
               Imag    = amr.patch[lv][PID]->fluid[IMAG][k][j][i];

               for (int v=0; v<NCOMP_PASSIVE; v++)
               Pass[v] = amr.patch[lv][PID]->fluid[ NCOMP_FLUID + v ][k][j][i];

               if ( OutputPot )
               Pot     = amr.patch[lv][PID]->pot        [k][j][i];

               if ( OutputParDens )
               ParDens = amr.patch[lv][PID]->par_dens   [k][j][i];


               fprintf( File, "%10.4e %10.4e %10.4e %2d %12.5e %12.5e %12.5e",
                        xx*dh_min, yy*dh_min, zz*dh_min, lv, Dens, Real, Imag );

               for (int v=0; v<NCOMP_PASSIVE; v++)    fprintf( File, " %12.5e", Pass[v] );
               if ( OutputPot     )                   fprintf( File, " %12.5e", Pot );
               if ( OutputParDens )                   fprintf( File, " %12.5e", ParDens );


               fprintf( File, "\n" );

#              else
#              error : ERROR : unsupported MODEL !!
#              endif // MODEL

            } // if ( Radius < MaxRadius )
         }}} // k, j, i
      } // for (int PID=0; PID<amr.num[lv]; PID++)

      cout << "done" << endl;

   } // for (int lv=0; lv<NLEVEL; lv++)


   fclose( File );

} // FUNCTION : GetMaxRho



//-------------------------------------------------------------------------------------------------------
// Function    :  GetRMS
// Description :  Evaluate the standard deviations from the average values
//-------------------------------------------------------------------------------------------------------
void GetRMS()
{

   cout << "Evaluating the RMS ..." << endl;


#  if   ( MODEL == HYDRO )
   const int NGhost = 0;
   real rho, vx, vy, vz, vr, vt, egy, pres, Pot, ParDens;

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   const int NGhost  = (ELBDM_GetVir) ? 1 : 0;
   const real _Eta   = 1.0/ELBDM_ETA;
   const real _2Eta  = 0.5*_Eta;
   const real _2Eta2 = _Eta*_2Eta;
   real Dens, Real, Imag, Pot, _Dens, ParDens;
   real Ek_Lap, Ek_Gra, GradR[3], GradI[3], LapR, LapI, _2dh, _dh2;
   real v[3], w[3], vr, vr_abs, vt_abs, wr, wr_abs, wt_abs, GradD[3];

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL

#  if ( MODEL != ELBDM )
   const bool ELBDM_IntPhase = false;
#  endif

   const int ArraySize = PATCH_SIZE + 2*NGhost;
   const int NPG       = 1;
   const NSide_t NSide = NSIDE_26;

   int    ShellID, Var, i, j, k, im, jm, km, ip, jp, kp, POTE, PAR_DENS, NextIdx;
   long   TVar;
   double Radius, scale, dv;
   double x, x1, x2, y, y1, y2, z, z1, z2;   // (x,y,z) : relative coordinates to the vector "Center"
   real   pass[NCOMP_PASSIVE];

   real *Field1D = new real [NPG*8*NIn*ArraySize*ArraySize*ArraySize];
   real (*Field)[NIn][ArraySize][ArraySize][ArraySize] = ( real(*)[NIn][ArraySize][ArraySize][ArraySize] )Field1D;


// determine the target variables
   TVar    = _TOTAL;
   NextIdx = NCOMP_TOTAL;

   if ( OutputPot     )    {  TVar |= _POTE;       POTE     = NextIdx ++;  }
   if ( OutputParDens )    {  TVar |= _PAR_DENS;   PAR_DENS = NextIdx ++;  }


   for (int lv=0; lv<NLEVEL; lv++)
   {
      scale = (double)amr.scale[lv];
      dv    = CUBE(scale);
#     if ( MODEL == ELBDM )
      _2dh  = 0.5/amr.dh[lv];
      _dh2  = 1.0/SQR(amr.dh[lv]);
#     endif

      cout << "   Level " << lv << " ... ";

      for (int PID0=0; PID0<amr.num[lv]; PID0+=8)
      {
//       skip useless patches
         if ( amr.patch[lv][PID0+0]->fluid == NULL  &&  amr.patch[lv][PID0+1]->fluid == NULL  &&
              amr.patch[lv][PID0+2]->fluid == NULL  &&  amr.patch[lv][PID0+3]->fluid == NULL  &&
              amr.patch[lv][PID0+4]->fluid == NULL  &&  amr.patch[lv][PID0+5]->fluid == NULL  &&
              amr.patch[lv][PID0+6]->fluid == NULL  &&  amr.patch[lv][PID0+7]->fluid == NULL    )  continue;

         if ( amr.patch[lv][PID0+0]->son != -1  &&  amr.patch[lv][PID0+1]->son != -1  &&
              amr.patch[lv][PID0+2]->son != -1  &&  amr.patch[lv][PID0+3]->son != -1  &&
              amr.patch[lv][PID0+4]->son != -1  &&  amr.patch[lv][PID0+5]->son != -1  &&
              amr.patch[lv][PID0+6]->son != -1  &&  amr.patch[lv][PID0+7]->son != -1    )    continue;


//       prepare data with ghost zones
         Prepare_PatchData( lv, Field[0][0][0][0], NGhost, NPG, &PID0, TVar, IntScheme, NSide, ELBDM_IntPhase );


//       evaluate the shell average
         for (int PID=PID0, p=0; PID<PID0+8; PID++, p++)
         {
            if ( amr.patch[lv][PID]->fluid == NULL  ||  amr.patch[lv][PID]->son != -1 )   continue;

            for (int kk=0; kk<PATCH_SIZE; kk++) {  z1 = amr.patch[lv][PID]->corner[2] + (kk+0.5)*scale - Center    [2];
                                                   z2 = amr.patch[lv][PID]->corner[2] + (kk+0.5)*scale - Center_Map[2];
                                                   z  = ( fabs(z1) <= fabs(z2) ) ? z1 : z2;
                                                   k  = kk + NGhost;   km = k - 1;   kp = k + 1;
            for (int jj=0; jj<PATCH_SIZE; jj++) {  y1 = amr.patch[lv][PID]->corner[1] + (jj+0.5)*scale - Center    [1];
                                                   y2 = amr.patch[lv][PID]->corner[1] + (jj+0.5)*scale - Center_Map[1];
                                                   y  = ( fabs(y1) <= fabs(y2) ) ? y1 : y2;
                                                   j  = jj + NGhost;   jm = j - 1;   jp = j + 1;
            for (int ii=0; ii<PATCH_SIZE; ii++) {  x1 = amr.patch[lv][PID]->corner[0] + (ii+0.5)*scale - Center    [0];
                                                   x2 = amr.patch[lv][PID]->corner[0] + (ii+0.5)*scale - Center_Map[0];
                                                   x  = ( fabs(x1) <= fabs(x2) ) ? x1 : x2;
                                                   i  = ii + NGhost;   im = i - 1;   ip = i + 1;

               Radius = sqrt( x*x + y*y + z*z );

               if ( Radius < MaxRadius )
               {
                  if ( LogBin > 1.0 )  ShellID = ( Radius < ShellWidth ) ? 0 : int( log(Radius/ShellWidth)/log(LogBin) ) + 1;
                  else                 ShellID = int( Radius / ShellWidth );

                  if ( ShellID >= NShell )
                  {
                     cerr << "ERROR : ShellID >= NShell !!" << endl;
                     exit( 1 );
                  }

#                 if   ( MODEL == HYDRO )
//                evaluate the values on the shell
                  rho     = Field[p][DENS    ][k][j][i];
                  vx      = Field[p][MOMX    ][k][j][i] / rho;
                  vy      = Field[p][MOMY    ][k][j][i] / rho;
                  vz      = Field[p][MOMZ    ][k][j][i] / rho;
                  egy     = Field[p][ENGY    ][k][j][i];

                  for (int u=0, uu=NCOMP_FLUID; u<NCOMP_PASSIVE; u++, uu++)
                  pass[u] = Field[p][uu      ][k][j][i];

                  if ( OutputPot )
                  Pot     = Field[p][POTE    ][k][j][i];

                  if ( OutputParDens )
                  ParDens = Field[p][PAR_DENS][k][j][i];

                  pres = (GAMMA-1.0) * ( egy - 0.5*rho*(vx*vx + vy*vy + vz*vz) );
                  vr   = ( x*vx + y*vy + z*vz ) / Radius;
                  vt   = sqrt( fabs(vx*vx + vy*vy + vz*vz - vr*vr) );


//                evalute the square of deviation on the shell
                  Var = 0;
                  RMS[ShellID][Var] += dv*pow( double(rho    )-Average[ShellID][Var], 2.0 );    Var++;
                  RMS[ShellID][Var] += dv*pow( double(vr     )-Average[ShellID][Var], 2.0 );    Var++;
                  RMS[ShellID][Var] += dv*pow( double(vt     )-Average[ShellID][Var], 2.0 );    Var++;
                  RMS[ShellID][Var] += dv*pow( double(egy    )-Average[ShellID][Var], 2.0 );    Var++;
                  RMS[ShellID][Var] += dv*pow( double(pres   )-Average[ShellID][Var], 2.0 );    Var++;

                  for (int v=0; v<NCOMP_PASSIVE; v++) {
                  RMS[ShellID][Var] += dv*pow( double(pass[v])-Average[ShellID][Var], 2.0 );    Var++; }

                  if ( OutputPot ) {
                  RMS[ShellID][Var] += dv*pow( double(Pot    )-Average[ShellID][Var], 2.0 );    Var++; }

                  if ( OutputParDens ) {
                  RMS[ShellID][Var] += dv*pow( double(ParDens)-Average[ShellID][Var], 2.0 );    Var++; }

#                 elif ( MODEL == MHD )
#                 warning : WAIT MHD !!!

#                 elif ( MODEL == ELBDM )
                  Dens    = Field[p][DENS    ][k][j][i];
                  Real    = Field[p][REAL    ][k][j][i];
                  Imag    = Field[p][IMAG    ][k][j][i];

                  for (int u=0, uu=NCOMP_FLUID; u<NCOMP_PASSIVE; u++, uu++)
                  pass[u] = Field[p][uu      ][k][j][i];

                  if ( OutputPot )
                  Pot     = Field[p][POTE    ][k][j][i];

                  if ( OutputParDens )
                  ParDens = Field[p][PAR_DENS][k][j][i];

                  if ( ELBDM_GetVir )
                  {
                     GradD[0] = _2dh*( Field[p][DENS][k ][j ][ip] - Field[p][DENS][k ][j ][im] );
                     GradD[1] = _2dh*( Field[p][DENS][k ][jp][i ] - Field[p][DENS][k ][jm][i ] );
                     GradD[2] = _2dh*( Field[p][DENS][kp][j ][i ] - Field[p][DENS][km][j ][i ] );

                     GradR[0] = _2dh*( Field[p][REAL][k ][j ][ip] - Field[p][REAL][k ][j ][im] );
                     GradR[1] = _2dh*( Field[p][REAL][k ][jp][i ] - Field[p][REAL][k ][jm][i ] );
                     GradR[2] = _2dh*( Field[p][REAL][kp][j ][i ] - Field[p][REAL][km][j ][i ] );

                     GradI[0] = _2dh*( Field[p][IMAG][k ][j ][ip] - Field[p][IMAG][k ][j ][im] );
                     GradI[1] = _2dh*( Field[p][IMAG][k ][jp][i ] - Field[p][IMAG][k ][jm][i ] );
                     GradI[2] = _2dh*( Field[p][IMAG][kp][j ][i ] - Field[p][IMAG][km][j ][i ] );

                     LapR     = ( Field[p][REAL][k ][j ][ip] + Field[p][REAL][k ][jp][i ] + Field[p][REAL][kp][j ][i ] +
                                  Field[p][REAL][k ][j ][im] + Field[p][REAL][k ][jm][i ] + Field[p][REAL][km][j ][i ] -
                                  6.0*Real )*_dh2;
                     LapI     = ( Field[p][IMAG][k ][j ][ip] + Field[p][IMAG][k ][jp][i ] + Field[p][IMAG][kp][j ][i ] +
                                  Field[p][IMAG][k ][j ][im] + Field[p][IMAG][k ][jm][i ] + Field[p][IMAG][km][j ][i ] -
                                  6.0*Imag )*_dh2;

                     Ek_Lap = -_2Eta2*( Real*LapR + Imag*LapI );
                     Ek_Gra = +_2Eta2*( SQR(GradR[0]) + SQR(GradR[1]) + SQR(GradR[2]) +
                                        SQR(GradI[0]) + SQR(GradI[1]) + SQR(GradI[2])   );

                     _Dens  = 1.0 / Dens;

                     for (int d=0; d<3; d++)
                     {
                        v[d] = _Eta*_Dens*( Real*GradI[d] - Imag*GradR[d] );
                        w[d] = _2Eta*_Dens*GradD[d];
                     }

                     vr     = ( x*v[0] + y*v[1] + z*v[2] ) / Radius;
                     vr_abs = fabs( vr );
                     vt_abs = sqrt(  fabs( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] - vr*vr )  );

                     wr     = ( x*w[0] + y*w[1] + z*w[2] ) / Radius;
                     wr_abs = fabs( wr );
                     wt_abs = sqrt(  fabs( w[0]*w[0] + w[1]*w[1] + w[2]*w[2] - wr*wr )  );
                  } // if ( ELBDM_GetVir )


//                evalute the square of deviation on the shell
                  Var = 0;
                  RMS[ShellID][Var] += dv*pow( double(Dens   )-Average[ShellID][Var], 2.0 );  Var++;
                  RMS[ShellID][Var] += dv*pow( double(Real   )-Average[ShellID][Var], 2.0 );  Var++;
                  RMS[ShellID][Var] += dv*pow( double(Imag   )-Average[ShellID][Var], 2.0 );  Var++;

                  for (int v=0; v<NCOMP_PASSIVE; v++) {
                  RMS[ShellID][Var] += dv*pow( double(pass[v])-Average[ShellID][Var], 2.0 );  Var++; }

                  if ( OutputPot   ) {
                  RMS[ShellID][Var] += dv*pow( double(Pot    )-Average[ShellID][Var], 2.0 );  Var++; }

                  if ( OutputParDens ) {
                  RMS[ShellID][Var] += dv*pow( double(ParDens)-Average[ShellID][Var], 2.0 );  Var++; }

                  if ( ELBDM_GetVir ) {
                  RMS[ShellID][Var] += dv*pow( double(Ek_Lap )-Average[ShellID][Var], 2.0 );  Var++;
                  RMS[ShellID][Var] += dv*pow( double(Ek_Gra )-Average[ShellID][Var], 2.0 );  Var++;
                  RMS[ShellID][Var] += dv*pow( double(vr     )-Average[ShellID][Var], 2.0 );  Var++;
                  RMS[ShellID][Var] += dv*pow( double(vr_abs )-Average[ShellID][Var], 2.0 );  Var++;
                  RMS[ShellID][Var] += dv*pow( double(vt_abs )-Average[ShellID][Var], 2.0 );  Var++;
                  RMS[ShellID][Var] += dv*pow( double(wr     )-Average[ShellID][Var], 2.0 );  Var++;
                  RMS[ShellID][Var] += dv*pow( double(wr_abs )-Average[ShellID][Var], 2.0 );  Var++;
                  RMS[ShellID][Var] += dv*pow( double(wt_abs )-Average[ShellID][Var], 2.0 );  Var++; }

#                 else
#                 error : ERROR : unsupported MODEL !!
#                 endif // MODEL

               } // if ( Radius < MaxRadius )
            }}} // kk, jj, ii
         } // for (int PID=PID0, p=0; PID<PID0+8; PID++, p++)
      } // for (int PID0=0; PID0<amr.num[lv]; PID0+=8)

      cout << "done" << endl;

   } // for (int lv=0; lv<NLEVEL; lv++)


// get the root-mean-square at each level
   for (int n=0; n<NShell; n++)
   for (int v=0; v<NOut; v++)    RMS[n][v] = sqrt( RMS[n][v]/Volume[n] );

   delete [] Field1D;

} // FUNCTION : GetRMS



//-------------------------------------------------------------------------------------------------------
// Function    :  ShellAverage
// Description :  Get the shell average of all variables
//-------------------------------------------------------------------------------------------------------
void ShellAverage()
{

   cout << "Evaluating the shell average ..." << endl;


#  if   ( MODEL == HYDRO )
   const int NGhost = 0;
   real px, py, pz, pr, pt, vr, vt, pres, rho, egy, Pot, ParDens;

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   const int  NGhost = (ELBDM_GetVir) ? 1 : 0;
   const real _Eta   = 1.0/ELBDM_ETA;
   const real _2Eta  = 0.5*_Eta;
   const real _2Eta2 = _Eta*_2Eta;
   const real _Eta2  = _Eta*_Eta;
   real Dens, Real, Imag, Pot, _Dens, ParDens;
   real Ek_Lap, Ek_Gra, GradR[3], GradI[3], LapR, LapI, _2dh, _dh2, dv_Eta;
   real v[3], w[3], vr, vr1, vr2, vt1, vt2, wr, wr1, wr2, wt1, wt2, GradD[3], dR_dr, dI_dr;

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL

#  if ( MODEL != ELBDM )
   const bool ELBDM_IntPhase = false;
#  endif

   const int ArraySize = PATCH_SIZE + 2*NGhost;
   const int NPG       = 1;
   const NSide_t NSide = NSIDE_26;

   int    ShellID, Var, i, j, k, im, jm, km, ip, jp, kp, POTE, PAR_DENS, NextIdx;
   long   TVar;
   double Radius, scale, dv;
   double x, x1, x2, y, y1, y2, z, z1, z2;   // (x,y,z) : relative coordinates to the vector "Center"
   real   pass[NCOMP_PASSIVE];

   real *Field1D = new real [NPG*8*NIn*ArraySize*ArraySize*ArraySize];
   real (*Field)[NIn][ArraySize][ArraySize][ArraySize] = ( real(*)[NIn][ArraySize][ArraySize][ArraySize] )Field1D;


// determine the target variables
   TVar    = _TOTAL;
   NextIdx = NCOMP_TOTAL;

   if ( OutputPot     )    {  TVar |= _POTE;       POTE     = NextIdx ++;  }
   if ( OutputParDens )    {  TVar |= _PAR_DENS;   PAR_DENS = NextIdx ++;  }


   for (int lv=0; lv<NLEVEL; lv++)
   {
      scale  = (double)amr.scale[lv];
      dv     = CUBE(scale);
#     if ( MODEL == ELBDM )
      _2dh   = 0.5/amr.dh[lv];
      _dh2   = 1.0/SQR(amr.dh[lv]);
      dv_Eta = dv/ELBDM_ETA*CUBE( amr.dh[NLEVEL-1] );
#     endif

      cout << "   Level " << lv << " ... ";

      for (int PID0=0; PID0<amr.num[lv]; PID0+=8)
      {
//       skip useless patches
         if ( amr.patch[lv][PID0+0]->fluid == NULL  &&  amr.patch[lv][PID0+1]->fluid == NULL  &&
              amr.patch[lv][PID0+2]->fluid == NULL  &&  amr.patch[lv][PID0+3]->fluid == NULL  &&
              amr.patch[lv][PID0+4]->fluid == NULL  &&  amr.patch[lv][PID0+5]->fluid == NULL  &&
              amr.patch[lv][PID0+6]->fluid == NULL  &&  amr.patch[lv][PID0+7]->fluid == NULL    )  continue;

         if ( amr.patch[lv][PID0+0]->son != -1  &&  amr.patch[lv][PID0+1]->son != -1  &&
              amr.patch[lv][PID0+2]->son != -1  &&  amr.patch[lv][PID0+3]->son != -1  &&
              amr.patch[lv][PID0+4]->son != -1  &&  amr.patch[lv][PID0+5]->son != -1  &&
              amr.patch[lv][PID0+6]->son != -1  &&  amr.patch[lv][PID0+7]->son != -1    )    continue;


//       prepare data with ghost zones
         Prepare_PatchData( lv, Field[0][0][0][0], NGhost, NPG, &PID0, TVar, IntScheme, NSide, ELBDM_IntPhase );


//       evaluate the shell average
         for (int PID=PID0, p=0; PID<PID0+8; PID++, p++)
         {
            if ( amr.patch[lv][PID]->fluid == NULL  ||  amr.patch[lv][PID]->son != -1 )   continue;

            for (int kk=0; kk<PATCH_SIZE; kk++) {  z1 = amr.patch[lv][PID]->corner[2] + (kk+0.5)*scale - Center    [2];
                                                   z2 = amr.patch[lv][PID]->corner[2] + (kk+0.5)*scale - Center_Map[2];
                                                   z  = ( fabs(z1) <= fabs(z2) ) ? z1 : z2;
                                                   k  = kk + NGhost;   km = k - 1;   kp = k + 1;
            for (int jj=0; jj<PATCH_SIZE; jj++) {  y1 = amr.patch[lv][PID]->corner[1] + (jj+0.5)*scale - Center    [1];
                                                   y2 = amr.patch[lv][PID]->corner[1] + (jj+0.5)*scale - Center_Map[1];
                                                   y  = ( fabs(y1) <= fabs(y2) ) ? y1 : y2;
                                                   j  = jj + NGhost;   jm = j - 1;   jp = j + 1;
            for (int ii=0; ii<PATCH_SIZE; ii++) {  x1 = amr.patch[lv][PID]->corner[0] + (ii+0.5)*scale - Center    [0];
                                                   x2 = amr.patch[lv][PID]->corner[0] + (ii+0.5)*scale - Center_Map[0];
                                                   x  = ( fabs(x1) <= fabs(x2) ) ? x1 : x2;
                                                   i  = ii + NGhost;   im = i - 1;   ip = i + 1;

               Radius = sqrt( x*x + y*y + z*z );

               if ( Radius < MaxRadius )
               {
                  if ( LogBin > 1.0 )  ShellID = ( Radius < ShellWidth ) ? 0 : int( log(Radius/ShellWidth)/log(LogBin) ) + 1;
                  else                 ShellID = int( Radius / ShellWidth );

                  if ( ShellID >= NShell )
                  {
                     cerr << "ERROR : ShellID >= NShell !!" << endl;
                     exit( 1 );
                  }

#                 if   ( MODEL == HYDRO )
//                evaluate the values on the shell
                  rho     = Field[p][DENS    ][k][j][i];
                  px      = Field[p][MOMX    ][k][j][i];
                  py      = Field[p][MOMY    ][k][j][i];
                  pz      = Field[p][MOMZ    ][k][j][i];
                  egy     = Field[p][ENGY    ][k][j][i];

                  for (int u=0, uu=NCOMP_FLUID; u<NCOMP_PASSIVE; u++, uu++)
                  pass[u] = Field[p][uu      ][k][j][i];

                  if ( OutputPot )
                  Pot     = Field[p][POTE    ][k][j][i];

                  if ( OutputParDens )
                  ParDens = Field[p][PAR_DENS][k][j][i];

                  pr   = ( x*px + y*py + z*pz ) / Radius;
                  pt   = sqrt( fabs(px*px + py*py + pz*pz - pr*pr) );
                  pres = (GAMMA-1.0) * ( egy - 0.5*(px*px + py*py + pz*pz)/rho );
                  vr   = pr/rho;
                  vt   = pt/rho;


//                sum up values at the same shell
                  Var = 0;
                  Average[ShellID][Var++] += (double)(dv*rho     );
                  Average[ShellID][Var++] += (double)(dv*pr      );
                  Average[ShellID][Var++] += (double)(dv*pt      );
                  Average[ShellID][Var++] += (double)(dv*egy     );
                  Average[ShellID][Var++] += (double)(dv*pres    );

                  for (int v=0; v<NCOMP_PASSIVE; v++)
                  Average[ShellID][Var++] += (double)(dv*pass[v] );

                  if ( OutputPot )
                  Average[ShellID][Var++] += (double)(dv*Pot     );

                  if ( OutputParDens )
                  Average[ShellID][Var++] += (double)(dv*ParDens );


//                store the maximum and minimum values
                  Var = 0;
                  if ( rho     > Max[ShellID][Var] )  Max[ShellID][Var] = rho;       Var++;
                  if ( vr      > Max[ShellID][Var] )  Max[ShellID][Var] = vr;        Var++;
                  if ( vt      > Max[ShellID][Var] )  Max[ShellID][Var] = vt;        Var++;
                  if ( egy     > Max[ShellID][Var] )  Max[ShellID][Var] = egy;       Var++;
                  if ( pres    > Max[ShellID][Var] )  Max[ShellID][Var] = pres;      Var++;

                  for (int v=0; v<NCOMP_PASSIVE; v++) {
                  if ( pass[v] > Max[ShellID][Var] )  Max[ShellID][Var] = pass[v];   Var++; }

                  if ( OutputPot ) {
                  if ( Pot     > Max[ShellID][Var] )  Max[ShellID][Var] = Pot;       Var++; }

                  if ( OutputParDens ) {
                  if ( ParDens > Max[ShellID][Var] )  Max[ShellID][Var] = ParDens;   Var++; }


                  Var = 0;
                  if ( rho     < Min[ShellID][Var] )  Min[ShellID][Var] = rho;      Var++;
                  if ( vr      < Min[ShellID][Var] )  Min[ShellID][Var] = vr;       Var++;
                  if ( vt      < Min[ShellID][Var] )  Min[ShellID][Var] = vt;       Var++;
                  if ( egy     < Min[ShellID][Var] )  Min[ShellID][Var] = egy;      Var++;
                  if ( pres    < Min[ShellID][Var] )  Min[ShellID][Var] = pres;     Var++;

                  for (int v=0; v<NCOMP_PASSIVE; v++) {
                  if ( pass[v] < Min[ShellID][Var] )  Min[ShellID][Var] = pass[v];  Var++; }

                  if ( OutputPot ) {
                  if ( Pot     < Min[ShellID][Var] )  Min[ShellID][Var] = Pot;      Var++; }

                  if ( OutputParDens ) {
                  if ( ParDens < Min[ShellID][Var] )  Min[ShellID][Var] = ParDens;  Var++; }

#                 elif ( MODEL == MHD )
#                 warning : WAIT MHD !!!

#                 elif ( MODEL == ELBDM )
//                evaluate the values on the shell
                  Dens    = Field[p][DENS    ][k][j][i];
                  Real    = Field[p][REAL    ][k][j][i];
                  Imag    = Field[p][IMAG    ][k][j][i];

                  for (int u=0, uu=NCOMP_FLUID; u<NCOMP_PASSIVE; u++, uu++)
                  pass[u] = Field[p][uu      ][k][j][i];

                  if ( OutputPot )
                  Pot     = Field[p][POTE    ][k][j][i];

                  if ( OutputParDens )
                  ParDens = Field[p][PAR_DENS][k][j][i];

                  if ( ELBDM_GetVir )
                  {
                     GradD[0] = _2dh*( Field[p][DENS][k ][j ][ip] - Field[p][DENS][k ][j ][im] );
                     GradD[1] = _2dh*( Field[p][DENS][k ][jp][i ] - Field[p][DENS][k ][jm][i ] );
                     GradD[2] = _2dh*( Field[p][DENS][kp][j ][i ] - Field[p][DENS][km][j ][i ] );

                     GradR[0] = _2dh*( Field[p][REAL][k ][j ][ip] - Field[p][REAL][k ][j ][im] );
                     GradR[1] = _2dh*( Field[p][REAL][k ][jp][i ] - Field[p][REAL][k ][jm][i ] );
                     GradR[2] = _2dh*( Field[p][REAL][kp][j ][i ] - Field[p][REAL][km][j ][i ] );

                     GradI[0] = _2dh*( Field[p][IMAG][k ][j ][ip] - Field[p][IMAG][k ][j ][im] );
                     GradI[1] = _2dh*( Field[p][IMAG][k ][jp][i ] - Field[p][IMAG][k ][jm][i ] );
                     GradI[2] = _2dh*( Field[p][IMAG][kp][j ][i ] - Field[p][IMAG][km][j ][i ] );

                     LapR     = ( Field[p][REAL][k ][j ][ip] + Field[p][REAL][k ][jp][i ] + Field[p][REAL][kp][j ][i ] +
                                  Field[p][REAL][k ][j ][im] + Field[p][REAL][k ][jm][i ] + Field[p][REAL][km][j ][i ] -
                                  6.0*Real )*_dh2;
                     LapI     = ( Field[p][IMAG][k ][j ][ip] + Field[p][IMAG][k ][jp][i ] + Field[p][IMAG][kp][j ][i ] +
                                  Field[p][IMAG][k ][j ][im] + Field[p][IMAG][k ][jm][i ] + Field[p][IMAG][km][j ][i ] -
                                  6.0*Imag )*_dh2;

                     Ek_Lap = -_2Eta2*( Real*LapR + Imag*LapI );
                     Ek_Gra = +_2Eta2*( SQR(GradR[0]) + SQR(GradR[1]) + SQR(GradR[2]) +
                                        SQR(GradI[0]) + SQR(GradI[1]) + SQR(GradI[2])   );

                     _Dens  = 1.0 / Dens;

                     for (int d=0; d<3; d++)
                     {
                        v[d] = _Eta*_Dens*( Real*GradI[d] - Imag*GradR[d] );
                        w[d] = _2Eta*_Dens*GradD[d];
                     }

                     vr  = ( x*v[0] + y*v[1] + z*v[2] ) / Radius;
                     vr2 = vr*vr;
                     vt2 = fabs( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] - vr2 );
                     vr1 = fabs( vr );
                     vt1 = sqrt( vt2 );

                     wr  = ( x*w[0] + y*w[1] + z*w[2] ) / Radius;
                     wr2 = wr*wr;
                     wt2 = fabs( w[0]*w[0] + w[1]*w[1] + w[2]*w[2] - wr2 );
                     wr1 = fabs( wr );
                     wt1 = sqrt( wt2 );
                  } // if ( ELBDM_GetVir )


//                sum up values at the same shell
                  Var = 0;
                  Average[ShellID][Var++] += (double)(dv*Dens     );
                  Average[ShellID][Var++] += (double)(dv*Real     );
                  Average[ShellID][Var++] += (double)(dv*Imag     );

                  for (int v=0; v<NCOMP_PASSIVE; v++)
                  Average[ShellID][Var++] += (double)(dv*pass[v] );

                  if ( OutputPot )
                  Average[ShellID][Var++] += (double)(dv*Pot      );

                  if ( OutputParDens )
                  Average[ShellID][Var++] += (double)(dv*ParDens  );

                  if ( ELBDM_GetVir ) {
                  Average[ShellID][Var++] += (double)(dv*Ek_Lap   );
                  Average[ShellID][Var++] += (double)(dv*Ek_Gra   );
                  Average[ShellID][Var++] += (double)(dv*Dens*vr  );
                  Average[ShellID][Var++] += (double)(dv*Dens*vr2 );
                  Average[ShellID][Var++] += (double)(dv*Dens*vt2 );
                  Average[ShellID][Var++] += (double)(dv*Dens*wr  );
                  Average[ShellID][Var++] += (double)(dv*Dens*wr2 );
                  Average[ShellID][Var++] += (double)(dv*Dens*wt2 );

                  for (int d=0; d<3; d++)    ELBDM_Mom[ShellID][d] += (double)( Real*GradI[d] - Imag*GradR[d] )*dv_Eta;

                  dR_dr                   = ( x*GradR[0] + y*GradR[1] + z*GradR[2] ) / Radius;
                  dI_dr                   = ( x*GradI[0] + y*GradI[1] + z*GradI[2] ) / Radius;
                  ELBDM_RhoUr2 [ShellID] += (double)dv*( SQR(dR_dr) + SQR(dI_dr) )*_Eta2;
                  ELBDM_dRho_dr[ShellID] += (double)dv*( x*GradD[0] + y*GradD[1] + z*GradD[2] ) / Radius;
                  ELBDM_LapRho [ShellID] += (double)dv*( Field[p][DENS][k ][j ][ip] + Field[p][DENS][k ][jp][i ] +
                                                         Field[p][DENS][kp][j ][i ] + Field[p][DENS][k ][j ][im] +
                                                         Field[p][DENS][k ][jm][i ] + Field[p][DENS][km][j ][i ] -
                                                         6.0*Dens )*_dh2; }


//                store the maximum and minimum values
                  Var = 0;
                  if ( Dens    > Max[ShellID][Var] )  Max[ShellID][Var] = Dens;     Var++;
                  if ( Real    > Max[ShellID][Var] )  Max[ShellID][Var] = Real;     Var++;
                  if ( Imag    > Max[ShellID][Var] )  Max[ShellID][Var] = Imag;     Var++;

                  for (int v=0; v<NCOMP_PASSIVE; v++) {
                  if ( pass[v] > Max[ShellID][Var] )  Max[ShellID][Var] = pass[v];  Var++; }

                  if ( OutputPot   ) {
                  if ( Pot     > Max[ShellID][Var] )  Max[ShellID][Var] = Pot;      Var++; }

                  if ( OutputParDens ) {
                  if ( ParDens > Max[ShellID][Var] )  Max[ShellID][Var] = ParDens;  Var++; }

                  if ( ELBDM_GetVir ) {
                  if ( Ek_Lap  > Max[ShellID][Var] )  Max[ShellID][Var] = Ek_Lap;   Var++;
                  if ( Ek_Gra  > Max[ShellID][Var] )  Max[ShellID][Var] = Ek_Gra;   Var++;
                  if ( vr      > Max[ShellID][Var] )  Max[ShellID][Var] = vr;       Var++;
                  if ( vr1     > Max[ShellID][Var] )  Max[ShellID][Var] = vr1;      Var++;
                  if ( vt1     > Max[ShellID][Var] )  Max[ShellID][Var] = vt1;      Var++;
                  if ( wr      > Max[ShellID][Var] )  Max[ShellID][Var] = wr;       Var++;
                  if ( wr1     > Max[ShellID][Var] )  Max[ShellID][Var] = wr1;      Var++;
                  if ( wt1     > Max[ShellID][Var] )  Max[ShellID][Var] = wt1;      Var++; }

                  Var = 0;
                  if ( Dens    < Min[ShellID][Var] )  Min[ShellID][Var] = Dens;     Var++;
                  if ( Real    < Min[ShellID][Var] )  Min[ShellID][Var] = Real;     Var++;
                  if ( Imag    < Min[ShellID][Var] )  Min[ShellID][Var] = Imag;     Var++;

                  for (int v=0; v<NCOMP_PASSIVE; v++) {
                  if ( pass[v] < Min[ShellID][Var] )  Min[ShellID][Var] = pass[v];  Var++; }

                  if ( OutputPot     ) {
                  if ( Pot     < Min[ShellID][Var] )  Min[ShellID][Var] = Pot;      Var++; }

                  if ( OutputParDens ) {
                  if ( ParDens < Min[ShellID][Var] )  Min[ShellID][Var] = ParDens;  Var++; }

                  if ( ELBDM_GetVir ) {
                  if ( Ek_Lap  < Min[ShellID][Var] )  Min[ShellID][Var] = Ek_Lap;   Var++;
                  if ( Ek_Gra  < Min[ShellID][Var] )  Min[ShellID][Var] = Ek_Gra;   Var++;
                  if ( vr      < Min[ShellID][Var] )  Min[ShellID][Var] = vr;       Var++;
                  if ( vr1     < Min[ShellID][Var] )  Min[ShellID][Var] = vr1;      Var++;
                  if ( vt1     < Min[ShellID][Var] )  Min[ShellID][Var] = vt1;      Var++;
                  if ( wr      < Min[ShellID][Var] )  Min[ShellID][Var] = wr;       Var++;
                  if ( wr1     < Min[ShellID][Var] )  Min[ShellID][Var] = wr1;      Var++;
                  if ( wt1     < Min[ShellID][Var] )  Min[ShellID][Var] = wt1;      Var++; }

#                 else
#                 error : ERROR : unsupported MODEL !!
#                 endif // MODEL

                  Volume[ShellID] += dv;
                  NCount[ShellID] ++;

               } // if ( Radius < MaxRadius )
            }}} // kk, jj, ii
         } // for (int PID=PID0, p=0; PID<PID0+8; PID++, p++)
      } // for (int PID0=0; PID0<amr.num[lv]; PID0+=8)

      cout << "done" << endl;

   } // for (int lv=0; lv<NLEVEL; lv++)


// get the average values
   for (int n=0; n<NShell; n++)
   for (int v=0; v<NOut; v++)    Average[n][v] /= Volume[n];

#  if ( MODEL == ELBDM )
   if ( ELBDM_GetVir )
   {
      for (int n=0; n<NShell; n++)
      {
         ELBDM_RhoUr2 [n] /= Volume[n];
         ELBDM_dRho_dr[n] /= Volume[n];
         ELBDM_LapRho [n] /= Volume[n];
      }
   }
#  endif


// get the average velocity
#  if   ( MODEL == HYDRO )
   for (int n=0; n<NShell; n++)
   {
//    <v> = <Rho*v>/<Rho>
      Average[n][1] /= Average[n][0];
      Average[n][2] /= Average[n][0];
   }

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   if ( ELBDM_GetVir )
   {
      const int Idx_v = NCOMP_TOTAL + 2 + ( (OutputPot)?1:0 ) + ( (OutputParDens)?1:0 );
      for (int n=0; n<NShell; n++)
      {
//       <v> = <Rho*v>/<Rho>, <|v|> = sqrt( <Rho*v^2>/<Rho> )
         for (int t=0; t<6; t++)    Average[n][Idx_v+t] /= Average[n][0];

         Average[n][Idx_v+1] = sqrt( Average[n][Idx_v+1] );
         Average[n][Idx_v+2] = sqrt( Average[n][Idx_v+2] );
         Average[n][Idx_v+4] = sqrt( Average[n][Idx_v+4] );
         Average[n][Idx_v+5] = sqrt( Average[n][Idx_v+5] );
      }
   }

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL


   delete [] Field1D;

} // FUNCTION : ShellAverage




//-------------------------------------------------------------------------------------------------------
// Function    :  ReadOption
// Description :  Read the command line options
//-------------------------------------------------------------------------------------------------------
void ReadOption( int argc, char **argv )
{

   int c;

   while ( (c = getopt(argc, argv, "hpsSMPVTDcgi:o:n:x:y:z:r:t:m:a:L:R:u:I:e:G:")) != -1 )
   {
      switch ( c )
      {
         case 'i': FileName_In      = optarg;
                   break;
         case 'e': FileName_Tree    = optarg;
                   break;
         case 'o': Suffix           = optarg;
                   break;
         case 'n': NShell           = atoi(optarg);
                   break;
         case 'x': Center[0]        = atof(optarg);
                   break;
         case 'y': Center[1]        = atof(optarg);
                   break;
         case 'z': Center[2]        = atof(optarg);
                   break;
         case 'r': MaxRadius        = atof(optarg);
                   break;
         case 'R': UseMaxRhoPos_R   = atof(optarg);
                   break;
         case 'p': Periodic         = true;
                   break;
         case 'm': UseMaxRhoPos     = atoi(optarg);
                   break;
         case 'D': UseMaxRhoPos_Par = true;
                   break;
         case 'S': Mode_ShellAve    = true;
                   break;
         case 'M': Mode_MaxRho      = true;
                   break;
         case 't': RhoThres         = atof(optarg);
                   break;
         case 's': InputScale       = true;
                   break;
         case 'L': LogBin           = atof(optarg);
                   break;
         case 'a': GetNShell        = atof(optarg);
                   break;
         case 'I': IntScheme        = (IntScheme_t)atoi(optarg);
                   break;
         case 'T': UseTree          = true;
                   break;
#        if ( MODEL == ELBDM )
         case 'P': ELBDM_IntPhase   = false;
                   break;
         case 'V': ELBDM_GetVir     = true;
                   break;
#        endif
         case 'u': INT_MONO_COEFF   = atof(optarg);
                   break;
         case 'c': OutputCenter     = true;
                   break;
         case 'g': GetAvePot        = true;
                   break;
         case 'G': NewtonG          = atof(optarg);
                   break;
         case 'h':
         case '?': cerr << endl << "usage: " << argv[0]
                        << " [-h (for help)] [-i input fileName] [-o suffix to the output file [none]]"
                        << endl << "                             "
                        << " [-n # of shells [0.5*BoxSize]] [-(x,y,z) sphere center [box center]]"
                        << endl << "                             "
                        << " [-r sphere radius [0.5*BoxSize]] [-p periodic B.C. [off]]"
                        << endl << "                             "
                        << " [-m # of highest-density cells for determining the sphere center (<=0: off) [8]]"
                        << endl << "                             "
                        << " [-D use particle (or total) density instead of just grid density to determine the center [off]]"
                        << endl << "                             "
                        << " [-R maximum radius for determining the sphere center (<=0: default) [sphere radius set by -r]]"
                        << endl << "                             "
                        << " [-a NCell (set minimum shell width equal to NCell*cell_size around center, <=0: off) [2.0]]"
                        << endl << "                             "
                        << " [-s [input cell scales instead of physical coordinates to specify the range] [off]]"
                        << endl << "                             "
                        << " [-L LogBinSize (size of the log scale radial bin, <=1.0: off -> use linear scale) [0.0]]"
                        << endl << "                             "
                        << " [-t density threshold [off]] [-c (output the center coords) [off]]"
                        << endl << "                             "
#                       if ( MODEL == ELBDM )
                        << " [-I interpolation scheme (1,2,3,4,5,6,7->MinMod-3D,MinMod-1D,vanLeer,CQuad,Quad,CQuar,Quar) [6]]"
                        << endl << "                             "
                        << " [-P (do not interpolate on phase) [interpolate on phase]]"
                        << endl << "                             "
                        << " [-V (analyze the ELBDM virial condition) [off]]"
#                       else
                        << " [-I interpolation scheme (1,2,3,4,5,6,7->MinMod-3D,MinMod-1D,vanLeer,CQuad,Quad,CQuar,Quar) [4]]"
#                       endif
                        << endl << "                             "
                        << " [-u coefficient for monotonic interpolation (1->4) [2.0])"
                        << endl << "                             "
                        << " [-g (calculate the spherically symmetric gravitational potential [false]] [-G Newton constant]"
                        << endl << "                             "
                        << " [-T (load the tree file) [off]] [-e tree file]"
                        << endl << "                             "
                        << " [-S (turn on the mode \"shell average\") [off]]"
                        << endl << "                             "
                        << " [-M (turn on the mode \"maximum density\") [off]]"
                        << endl << endl;
                   exit( 1 );

      } // switch ( c )
   } // while ..


// set variables
#  if ( MODEL == ELBDM )
   if ( ELBDM_GetVir )  NeedGhost = true;
#  endif


// initial check
   if ( fopen( FileName_In, "rb" ) == NULL )
   {
      fprintf( stderr, "ERROR : the input file \"%s\" does not exist (-i input filename) !!\n", FileName_In );
      exit( 1 );
   }

   if ( UseTree && fopen( FileName_Tree, "rb" ) == NULL )
   {
      fprintf( stderr, "ERROR : the input tree file \"%s\" does not exist (-e tree file) !!\n", FileName_Tree );
      exit( 1 );
   }

   if ( NeedGhost  &&  UseTree )
   {
      fprintf( stderr, "ERROR : currently the option \"UseTree\" does not work when ghost zones are required !!\n" );
      exit( 1 );
   }

   if ( !Mode_ShellAve  &&  !Mode_MaxRho )
   {
      fprintf( stderr, "ERROR : no operation mode is turned on (-S, -M) !!\n" );
      exit( 1 );
   }

   if ( Mode_MaxRho  &&  RhoThres == WRONG )
   {
      fprintf( stderr, "ERROR : please provide the density threshold (-t density threshold) !!\n" );
      exit( 1 );
   }

   if ( INT_MONO_COEFF < 1.0  ||  INT_MONO_COEFF > 4.0 )
      Aux_Error( ERROR_INFO, "INT_MONO_COEFF (%14.7e) is not within the correct range (1<=INT_MONO_COEFF<=4) !!\n",
                 INT_MONO_COEFF );

   if ( IntScheme != INT_MINMOD3D  &&  IntScheme != INT_MINMOD1D  &&
        IntScheme != INT_VANLEER   &&  IntScheme != INT_CQUAD     &&
        IntScheme != INT_QUAD      &&  IntScheme != INT_CQUAR     &&
        IntScheme != INT_QUAR )
      Aux_Error( ERROR_INFO, "unsupported interpolation scheme (%d) !!\n", IntScheme );

#  if ( MODEL == ELBDM )
   if ( ELBDM_IntPhase  &&  IntScheme == INT_MINMOD1D )
      Aux_Error( ERROR_INFO, "MinMod1D does not support interpolation on phase !!\n" );

   if ( ELBDM_GetVir  &&  !GetAvePot )
   {
      GetAvePot = true;

      fprintf( stdout, "NOTE : option \"GetAvePot (-g)\" is turned on automatically for \"ELBDM_GetVir (-V)\"\n" );
   }
#  endif

   if ( GetAvePot  &&  NewtonG == WRONG )
      Aux_Error( ERROR_INFO, "please set \"NewtonG (-G)\" for the option \"GetAvePot\" !!\n" );

#  ifndef DENS
   if ( GetAvePot )
      Aux_Error( ERROR_INFO, "DENS field is NOT defined for the option \"GetAvePot\" !!\n" );
#  endif

   if ( UseMaxRhoPos_Par  &&  UseMaxRhoPos <= 0 )
         Aux_Error( ERROR_INFO, "UseMaxRhoPos_Par (-D) is useless when UseMaxRhoPos (-m) <= 0 !!\n" );

#  if ( MODEL == HYDRO )
      Aux_Message( stderr, "WARNING : pressure is computed by constant-gamma EoS !!\n" );
#  endif

} // FUNCTION : ReadOption



//-------------------------------------------------------------------------------------------------------
// Function    :  GetMinShellWidth
// Description :  Return the size (scale) of the cell closest to the input sphere center
//
// Parameter   :  TCen     : Center of the targeted sphere
//                TCen_Map : Center of the mapped targeted sphere (useful only if Periodic is enabled)
//-------------------------------------------------------------------------------------------------------
double GetMinShellWidth( const double TCen[], const double TCen_Map[] )
{

   cout << "   GetMinShellWidth ..." << endl;


   double x, x1, x2, xx, y, y1, y2, yy, z, z1, z2, zz;
   double r, r_min, scale, ShellWidth_min;
   int    Lv_min;

   r_min  = __FLT_MAX__;
   Lv_min = 0;


   for (int lv=0; lv<NLEVEL; lv++)
   {
      scale = amr.scale[lv];

      for (int PID=0; PID<amr.num[lv]; PID++)
      {
         if ( amr.patch[lv][PID]->fluid == NULL  ||  amr.patch[lv][PID]->son != -1 )   continue;

         for (int k=0; k<PATCH_SIZE; k++) {  zz = amr.patch[lv][PID]->corner[2] + (k+0.5)*scale;
                                             z1 = zz - TCen    [2];
                                             z2 = zz - TCen_Map[2];
                                             z  = ( fabs(z1) <= fabs(z2) ) ?  z1 : z2;
         for (int j=0; j<PATCH_SIZE; j++) {  yy = amr.patch[lv][PID]->corner[1] + (j+0.5)*scale;
                                             y1 = yy - TCen    [1];
                                             y2 = yy - TCen_Map[1];
                                             y  = ( fabs(y1) <= fabs(y2) ) ?  y1 : y2;
         for (int i=0; i<PATCH_SIZE; i++) {  xx = amr.patch[lv][PID]->corner[0] + (i+0.5)*scale;
                                             x1 = xx - TCen    [0];
                                             x2 = xx - TCen_Map[0];
                                             x  = ( fabs(x1) <= fabs(x2) ) ?  x1 : x2;

            r = sqrt( x*x + y*y + z*z );

            if ( r < r_min )
            {
               r_min  = r;
               Lv_min = lv;
            }

         }}} // i,j,k
      } // for (int PID=0; PID<amr.num[lv]; PID++)
   } // for (int lv=0; lv<NLEVEL; lv++)


   ShellWidth_min = amr.scale[Lv_min];

   printf( "      Minimum distance = %13.7e\n", r_min );
   printf( "      Level            = %d\n",     Lv_min );
   printf( "      Shell scale      = %13.7e\n", ShellWidth_min );


   cout << "   GetMinShellWidth ... done" << endl;

   return ShellWidth_min;

} // FUNCTION : GetMinShellWidth



//-------------------------------------------------------------------------------------------------------
// Function    :  Output_ShellAve
// Description :  Output data for the mode "shell average"
//-------------------------------------------------------------------------------------------------------
void Output_ShellAve()
{

   cout << "Output_ShellAve ... " << flush;


// set output file names
   const int NOutMax = NCOMP_TOTAL + 20;  // 20 is added arbitrarily
         int Var     = NCOMP_TOTAL;

#  if   ( MODEL == HYDRO )
   char FileName[NOutMax][50] = { "AveRho", "AveV_R", "AveV_T", "AveEgy", "AvePre" };

   for (int v=0; v<NCOMP_PASSIVE; v++)   sprintf( FileName[ NCOMP_FLUID + v ], "AvePassive%02d", v );

   if      ( OutputPot )            sprintf( FileName[Var++], "%s", "AvePot"     );
   if      ( OutputParDens == 1 )   sprintf( FileName[Var++], "%s", "AveParDens" );
   else if ( OutputParDens == 2 )   sprintf( FileName[Var++], "%s", "AveTotDens" );

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   char FileName[NOutMax][50] = { "AveDens", "AveReal", "AveImag" };

   for (int v=0; v<NCOMP_PASSIVE; v++)   sprintf( FileName[ NCOMP_FLUID + v ], "AvePassive%02d", v );

   if      ( OutputPot )            sprintf( FileName[Var++], "%s", "AvePote" );
   if      ( OutputParDens == 1 )   sprintf( FileName[Var++], "%s", "AveParDens" );
   else if ( OutputParDens == 2 )   sprintf( FileName[Var++], "%s", "AveTotDens" );
   if ( ELBDM_GetVir )
   {
                        sprintf( FileName[Var++], "%s", "AveEk-L"  );
                        sprintf( FileName[Var++], "%s", "AveEk-G"  );
                        sprintf( FileName[Var++], "%s", "AveVr-N"  );   // N = net = inflow - outflow
                        sprintf( FileName[Var++], "%s", "AveVr-A"  );   // A = absolute = |Vr|
                        sprintf( FileName[Var++], "%s", "AveVt-A"  );
                        sprintf( FileName[Var++], "%s", "AveWr-N"  );
                        sprintf( FileName[Var++], "%s", "AveWr-A"  );
                        sprintf( FileName[Var++], "%s", "AveWt-A"  );
   }

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL

   if ( Suffix != NULL )
      for (int v=0; v<NOut; v++)    strcat( FileName[v], Suffix );


// output data
   const int    PAR_DENS    = ( OutputParDens ) ? ( (OutputPot)?NCOMP_TOTAL+1:NCOMP_TOTAL ) : -1;
   const double dh_min      = amr.dh[NLEVEL-1];
   const double LogBinCoeff = pow( (double)LogBin, -0.5 );  // coefficient to calculate the average radius of each log bin
   double AccMass           = 0.0;                          // accumulated mass
   double AccVolume         = 0.0;                          // accumulated volume
   bool   OutputAccMass     = false;                        // output accumulative mass
   bool   OutputAvePot      = false;                        // output spherically symmetric gravitational potential
#  if ( MODEL == ELBDM )
   double ELBDM_AccEk       = 0.0;                          // accumulated Ek
   double ELBDM_AccEk_NoCM  = 0.0;                          // accumulated Ek with CM Ek subtracted
   double ELBDM_AccMom[3]   = { 0.0, 0.0, 0.0 };            // accumulated momentum
   bool   ELBDM_OutputAccEk = false;
   double *ELBDM_Ek_L       = NULL;
   double *ELBDM_Ek_G       = NULL;
#  endif
   double r, dr;



// get the correct volume (scale -> dh)
   for (int n=0; n<NShell; n++)  Volume[n] *= CUBE(dh_min);


// allocate arrays
#  if ( MODEL == ELBDM )
   if ( ELBDM_GetVir )
   {
      ELBDM_Ek_L = new double [NShell];
      ELBDM_Ek_G = new double [NShell];
   }
#  endif


   for (int v=0; v<NOut; v++)
   {
//    determine whether or not to output the accumulated quantities
#     ifdef DENS
      if (  v == DENS  ||  v == PAR_DENS  ||  ( v >= NCOMP_FLUID && v < NCOMP_TOTAL )  )
      {
         OutputAccMass = true;
         OutputAvePot  = GetAvePot && ( v==DENS || v==PAR_DENS );    // currently we do not calculate potential for the passive scalars
         AccMass       = 0.0;
         AccVolume     = 0.0;
      }
      else
      {
         OutputAccMass = false;
         OutputAvePot  = false;
      }
#     endif

#     if ( MODEL == ELBDM )
      if (   (  v == ( NCOMP_TOTAL+((OutputPot)?1:0)+((OutputParDens)?1:0)   )  )  ||
             (  v == ( NCOMP_TOTAL+((OutputPot)?1:0)+((OutputParDens)?1:0)+1 )  )    )
      {
         ELBDM_OutputAccEk = true;
         ELBDM_AccEk       = 0.0;

         for (int d=0; d<3; d++)
         ELBDM_AccMom[d]   = 0.0;

         AccMass           = 0.0;
      }
      else
         ELBDM_OutputAccEk = false;
#     endif


//    header
      if ( NULL != fopen(FileName[v],"r") )
         fprintf( stderr, "Warning : the file \"%s\" already exists and will be overwritten !!\n", FileName[v] );

      FILE *File = fopen( FileName[v], "w" );

      fprintf( File, "#%19s  %10s  %13s  %13s  %13s  %13s", "Radius", "NCount", "Ave", "RMS", "Max", "Min" );

      if ( OutputAccMass )
      fprintf( File, "  %13s  %13s", "AccMass", "AveDens" );

#     if ( MODEL == ELBDM )
      if ( ELBDM_OutputAccEk )
      fprintf( File, "  %13s  %13s", "AccEk", "AccEk-NoCM" );
#     endif

      if ( OutputAvePot )
      fprintf( File, "  %13s", "AvePot" );

      fprintf( File, "\n" );


//    data
      for (int n=0; n<NShell; n++)
      {
//       skip bins without any data
         if ( NCount[n] == 0 )   continue;

//       determine the average radius of the target bin
         GetR( n, r, dr );

//       print data
         fprintf( File, "%20.14e  %10ld  %13.6e  %13.6e  %13.6e  %13.6e",
                  r, NCount[n], Average[n][v], RMS[n][v], Max[n][v], Min[n][v] );

//       get the accumulated mass
         if ( OutputAccMass )
         {
            AccMass   += Average[n][v]*Volume[n];
            AccVolume += Volume[n];

            fprintf( File, "  %13.6e  %13.6e", AccMass, AccMass/AccVolume );
         }


//       get the accumulated kinematic energy for ELBDM
#        if ( MODEL == ELBDM )
         if ( ELBDM_OutputAccEk )
         {
            ELBDM_AccEk += Average[n][v]*Volume[n];

//          subtract the CM Ek
            for (int d=0; d<3; d++)    ELBDM_AccMom[d] += ELBDM_Mom[n][d];
            AccMass += Average[n][0]*Volume[n];
            ELBDM_AccEk_NoCM = ELBDM_AccEk - 0.5/AccMass*( SQR(ELBDM_AccMom[0]) +
                                                           SQR(ELBDM_AccMom[1]) +
                                                           SQR(ELBDM_AccMom[2]) );

            fprintf( File, "  %13.6e  %13.6e", ELBDM_AccEk, ELBDM_AccEk_NoCM );

//          record Ek for analyzing the virial condition
            if      (  v == ( NCOMP_TOTAL+((OutputPot)?1:0)+((OutputParDens)?1:0)   )  )  ELBDM_Ek_L[n] = ELBDM_AccEk_NoCM;
            else if (  v == ( NCOMP_TOTAL+((OutputPot)?1:0)+((OutputParDens)?1:0)+1 )  )  ELBDM_Ek_G[n] = ELBDM_AccEk_NoCM;
            else
               Aux_Error( ERROR_INFO, "incorrect target component (%d) !!\n", v );
         } // if ( ELBDM_OutputAccEk )
#        endif // if ( MODEL == ELBDM )


//       calculate the spherically symmetric gravitational potential
         if ( OutputAvePot )
         {
            double rr, drr, AvePot=0.0;

            for (int m=0; m<NShell; m++)
            {
//             skip bins without any data
               if ( NCount[m] == 0 )   continue;

               GetR( m, rr, drr );

               if ( rr < r )  AvePot += Average[m][v]*drr*rr*rr/r;
               else           AvePot += Average[m][v]*drr*rr;
            } // for (int m=0; m<NShell; m++)

            AvePot *= -4.0*M_PI*NewtonG;

            fprintf( File, "  %13.6e", AvePot );
         } // if ( OutputAvePot )

         fprintf( File, "\n" );
      } // for (int n=0; n<NShell; n++)

      fclose( File );

   } // for (int v=0; v<NOut; v++)


// output the surface terms in the ELBDM virial condition
#  if ( MODEL == ELBDM )
   if ( ELBDM_GetVir )
   {
      const int Idx_vr_abs = NCOMP_TOTAL + ( (OutputPot)?4:3 ) + ( (OutputParDens)?1:0 );
      const int Idx_wr_abs = NCOMP_TOTAL + ( (OutputPot)?7:6 ) + ( (OutputParDens)?1:0 );

      double Ep, AvePot, mr, mdr, tr, tdr;


//    open file and header
      char FileName_VirSurf[50] = "VirSurf";

      if ( Suffix != NULL )   strcat( FileName_VirSurf, Suffix );

      if ( NULL != fopen(FileName_VirSurf,"r") )
         fprintf( stderr, "Warning : the file \"%s\" already exists and will be overwritten !!\n", FileName_VirSurf );

      FILE *File_VirSurf = fopen( FileName_VirSurf, "w" );

      fprintf( File_VirSurf, "#%19s  %10s  %13s  %13s  %13s  %13s  %13s  %13s  %13s  %13s  %13s  %13s\n",
               "Radius", "NCount", "Vr2", "Wr2", "Ur2", "dRho_dr", "Lap(Rho)", "2*Ek-L", "2*Ek-G", "Sum-L", "Sum-G", "Ep" );


//    data
      for (int n=0; n<NShell; n++)
      {
//       skip bins without any data
         if ( NCount[n] == 0 )   continue;

//       determine the average radius of the target bin
         GetR( n, r, dr );

//       calculate the surface terms
         const double DensRVr2 = -4.0*M_PI*CUBE( r )*SQR( Average[n][Idx_vr_abs] )*Average[n][0];
         const double DensRWr2 = -4.0*M_PI*CUBE( r )*SQR( Average[n][Idx_wr_abs] )*Average[n][0];
         const double DensRUr2 = -4.0*M_PI*CUBE( r )            *ELBDM_RhoUr2 [n];
         const double dRho_dr  = -M_PI*SQR( r )/SQR( ELBDM_ETA )*ELBDM_dRho_dr[n];
         const double LapRho   = M_PI*CUBE( r )/SQR( ELBDM_ETA )*ELBDM_LapRho [n];

//       calculate the cumulative potential energy ( 0.5*Integrate(Rho*Phi) ) --> excluding contribution outside the radius
         Ep = 0.0;

         for (int m=0; m<=n; m++)
         {
            AvePot = 0.0;

//          skip bins without any data
            if ( NCount[m] == 0 )   continue;

            GetR( m, mr, mdr );

            for (int t=0; t<=n; t++)
            {
//             skip bins without any data
               if ( NCount[t] == 0 )   continue;

               GetR( t, tr, tdr );

               if ( tr < mr )    AvePot += Average[t][0]*tdr*tr*tr/mr;
               else              AvePot += Average[t][0]*tdr*tr;
            }

            AvePot *= -4.0*M_PI*NewtonG;
            Ep     += 0.5*4.0*M_PI*SQR(mr)*mdr*Average[m][0]*AvePot;
         }

//       print data
         fprintf( File_VirSurf, "%20.14e  %10ld  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e\n",
                  r, NCount[n], DensRVr2, DensRWr2, DensRUr2, dRho_dr, LapRho, 2.0*ELBDM_Ek_L[n], 2.0*ELBDM_Ek_G[n],
                  DensRUr2+dRho_dr+LapRho+2.0*ELBDM_Ek_L[n], DensRUr2+3.0*dRho_dr+LapRho+2.0*ELBDM_Ek_G[n], Ep );
      } // for (int n=0; n<NShell; n++)

      fclose( File_VirSurf );

   } // if ( ELBDM_GetVir )
#  endif // #if ( MODEL == ELBDM )


// output the center coords.
   if ( OutputCenter )
   {
      const char FileName_Cen[] = "CenterCoords";
      bool OutputHeader         = true;

//    check whether the file already exists
      FILE *File_Check = fopen( FileName_Cen, "r" );
      if ( File_Check != NULL )
      {
         OutputHeader = false;   // do not output header if file already exists

         fclose( File_Check );
      }

//    output
      FILE *File_Cen = fopen( FileName_Cen, "a" );

      if ( OutputHeader )     fprintf( File_Cen, "%7s   %7s   %20s   %20s   %20s   %20s\n",
                                       "DumpID", "Step", "Time", "x", "y", "z" );

      fprintf( File_Cen, "%7d   %7ld   %20.14e   %20.14e   %20.14e   %20.14e\n",
               DumpID, Step, Time[0], Center[0]*dh_min, Center[1]*dh_min, Center[2]*dh_min );

      fclose( File_Cen );

   } // if ( OutputCenter )


// free memory
#  if ( MODEL == ELBDM )
   if ( ELBDM_GetVir )
   {
      delete [] ELBDM_Ek_L;
      delete [] ELBDM_Ek_G;
   }
#  endif


   cout << "done" << endl;

} // FUNCTION : Output_ShellAve



//-------------------------------------------------------------------------------------------------------
// Function    :  CheckParameter
// Description :  Verify parameters
//-------------------------------------------------------------------------------------------------------
void CheckParameter()
{

   cout << "CheckParameter ... " << flush;


#  if ( MODEL != HYDRO  &&  MODEL != ELBDM )
#     error : ERROR : unsupported MODEL !!
#  endif

#  ifndef DENS
#     error : ERROR : symbolic constant DENS is not defined !!
#  endif

   if (  Center[0] > amr.BoxScale[0]  ||  Center[0] < 0.0  )
   {
      fprintf( stderr, "ERROR : the x coordinate of the sphere center lies outside the simulation box !!\n" );
      exit( 1 );
   }

   if (  Center[1] > amr.BoxScale[1]  ||  Center[1] < 0.0  )
   {
      fprintf( stderr, "ERROR : the y coordinate of the sphere center lies outside the simulation box !!\n" );
      exit( 1 );
   }

   if (  Center[2] > amr.BoxScale[2]  ||  Center[2] < 0.0  )
   {
      fprintf( stderr, "ERROR : the z coordinate of the sphere center lies outside the simulation box !!\n" );
      exit( 1 );
   }

   if ( MaxRadius > 0.5*amr.BoxScale[0] )
      fprintf( stderr, "Warning : the sphere radius exceeds the half box size in the x direction !!\n" );

   if ( MaxRadius > 0.5*amr.BoxScale[1] )
      fprintf( stderr, "Warning : the sphere radius exceeds the half box size in the y direction !!\n" );

   if ( MaxRadius > 0.5*amr.BoxScale[2] )
      fprintf( stderr, "Warning : the sphere radius exceeds the half box size in the z direction !!\n" );


// checks for the mode "shell average"
   if ( Mode_ShellAve )
   {
      if ( NShell <= 0 )
      {
         fprintf( stderr, "ERROR : please provide the correct number of shells (-n # of shells) !!\n" );
         exit( 1 );
      }
   }


   cout << "done" << endl;

} // FUNCTION : CheckParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  End
// Description :  Deallocate memory
//-------------------------------------------------------------------------------------------------------
void End()
{

   int NPatch;

   for (int lv=NLEVEL-1; lv>=0; lv--)
   {
      NPatch = amr.num[lv];

      for (int PID=0; PID<NPatch; PID++)  amr.pdelete( lv, PID );
   }

   if ( NCount        != NULL )  delete [] NCount;
   if ( Volume        != NULL )  delete [] Volume;

   if ( Average[0]    != NULL )  delete [] Average[0];
   if ( RMS    [0]    != NULL )  delete [] RMS    [0];
   if ( Max    [0]    != NULL )  delete [] Max    [0];
   if ( Min    [0]    != NULL )  delete [] Min    [0];

   if ( Average       != NULL )  delete [] Average;
   if ( RMS           != NULL )  delete [] RMS;
   if ( Max           != NULL )  delete [] Max;
   if ( Min           != NULL )  delete [] Min;

#  if ( MODEL == ELBDM )
   if ( ELBDM_Mom     != NULL )  delete [] ELBDM_Mom;
   if ( ELBDM_RhoUr2  != NULL )  delete [] ELBDM_RhoUr2;
   if ( ELBDM_dRho_dr != NULL )  delete [] ELBDM_dRho_dr;
   if ( ELBDM_LapRho  != NULL )  delete [] ELBDM_LapRho;
#  endif

   if ( BaseP         != NULL )  delete [] BaseP;
   if ( tree          != NULL )  delete tree;

} // FUNCTION : End



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ShellAve
// Description :  Allocate memory and initialize all variables for the mode "shell average"
//-------------------------------------------------------------------------------------------------------
void Init_ShellAve()
{

   cout << "Init_ShellAve ... " << flush;


// allocate memory
   NCount     = new long int [NShell];
   Volume     = new double   [NShell];

   Average    = new double  *[NShell];
   RMS        = new double  *[NShell];
   Max        = new real    *[NShell];
   Min        = new real    *[NShell];

   Average[0] = new double   [NShell*NOut];
   RMS    [0] = new double   [NShell*NOut];
   Max    [0] = new real     [NShell*NOut];
   Min    [0] = new real     [NShell*NOut];

   for (int t=1; t<NShell; t++)
   {
      Average[t] = Average[0] + t*NOut;
      RMS    [t] = RMS    [0] + t*NOut;
      Max    [t] = Max    [0] + t*NOut;
      Min    [t] = Min    [0] + t*NOut;
   }

#  if ( MODEL == ELBDM )
   if ( ELBDM_GetVir )
   {
      ELBDM_Mom     = new double [NShell][3];
      ELBDM_RhoUr2  = new double [NShell];
      ELBDM_dRho_dr = new double [NShell];
      ELBDM_LapRho  = new double [NShell];
   }
#  endif


// initialize variables for storing average data
   for (int n=0; n<NShell; n++)
   {
      for (int v=0; v<NOut; v++)
      {
         Average[n][v] = 0.0;
         RMS    [n][v] = 0.0;
         Max    [n][v] = -__FLT_MAX__;
         Min    [n][v] = +__FLT_MAX__;
      }

      NCount[n] = 0;
      Volume[n] = 0.0;

#     if ( MODEL == ELBDM )
      if ( ELBDM_GetVir )
      {
         for (int d=0; d<3; d++)    ELBDM_Mom[n][d] = 0.0;

         ELBDM_RhoUr2 [n] = 0.0;
         ELBDM_dRho_dr[n] = 0.0;
         ELBDM_LapRho [n] = 0.0;
      }
#     endif
   } // for (int n=0; n<NShell; n++)


   cout << "done" << endl;

} // FUNCTION : Init_ShellAve



//-------------------------------------------------------------------------------------------------------
// Function    :  SetMaxRhoPos
// Description :  Set the sphere center to the center of mass of the AveN highest-density cells
//
// Parameter   :  AveN  : Number of highest-density cells for determining the sphere center
//-------------------------------------------------------------------------------------------------------
void SetMaxRhoPos( const int AveN )
{

   cout << "SetMaxRhoPos ..." << endl;


// check
   if ( AveN <= 0 )
   {
      fprintf( stderr, "ERROR : AveN (%d) <= 0 in %s !!\n", AveN, __FUNCTION__ );
      exit( 1 );
   }

   if ( UseMaxRhoPos_Par  &&  OutputParDens == 0 )
         Aux_Error( ERROR_INFO, "UseMaxRhoPos_Par (-D) is useless when no particle (or total) density is stored on grids !!\n" );


   const double dh_min = amr.dh[NLEVEL-1];
   real  MinRho        = -1.0;
   int   MinRho_ID     = 0;

   real   MaxRho[AveN], Rho;
   double scale, x, y, z, dx1, dx2, dx, dy1, dy2, dy, dz1, dz2, dz, Radius, MaxRho_Pos[AveN][3];
   real   PeakOutside_MaxRho = -1.0, PeakOutside_Pos[3], PeakOutside_R;
   bool   PeakOutside_First  = true;
   int    PeakOutside_Lv;
   int    MaxRho_Lv[AveN];


// initialize the array recording the maximum density
   for (int t=0; t<AveN; t++)    MaxRho[t] = MinRho;


// begin to search the maximum density
   for (int lv=0; lv<NLEVEL; lv++)
   {
      scale = (double)amr.scale[lv];

      for (int PID=0; PID<amr.num[lv]; PID++)
      {
         if ( amr.patch[lv][PID]->fluid == NULL  ||  amr.patch[lv][PID]->son != -1 )   continue;

         for (int k=0; k<PATCH_SIZE; k++) {  z   = amr.patch[lv][PID]->corner[2] + (k+0.5)*scale;
                                             dz1 = z - Center    [2];
                                             dz2 = z - Center_Map[2];
                                             dz  = ( fabs(dz1) <= fabs(dz2) ) ? dz1 : dz2;
         for (int j=0; j<PATCH_SIZE; j++) {  y   = amr.patch[lv][PID]->corner[1] + (j+0.5)*scale;
                                             dy1 = y - Center    [1];
                                             dy2 = y - Center_Map[1];
                                             dy  = ( fabs(dy1) <= fabs(dy2) ) ? dy1 : dy2;
         for (int i=0; i<PATCH_SIZE; i++) {  x   = amr.patch[lv][PID]->corner[0] + (i+0.5)*scale;
                                             dx1 = x - Center    [0];
                                             dx2 = x - Center_Map[0];
                                             dx  = ( fabs(dx1) <= fabs(dx2) ) ? dx1 : dx2;

            if ( UseMaxRhoPos_Par )    Rho = amr.patch[lv][PID]->par_dens   [k][j][i];
            else                       Rho = amr.patch[lv][PID]->fluid[DENS][k][j][i];

            Radius = sqrt( dx*dx + dy*dy + dz*dz );

            if ( Rho > MinRho  &&  Radius <= UseMaxRhoPos_R )
            {
//             record the density and coordinates
               MaxRho    [MinRho_ID]    = Rho;
               MaxRho_Lv [MinRho_ID]    = lv;
               MaxRho_Pos[MinRho_ID][0] = x;
               MaxRho_Pos[MinRho_ID][1] = y;
               MaxRho_Pos[MinRho_ID][2] = z;

//             find the new density threshold
               MinRho = __FLT_MAX__;

               for (int t=0; t<AveN; t++)
               {
                  if ( MaxRho[t] < MinRho )
                  {
                     MinRho    = MaxRho[t];
                     MinRho_ID = t;
                  }
               }
            } // if ( Rho > MinRho  &&  Radius <= UseMaxRhoPos_R )
         }}} // i, j, k
      } // for (int PID=0; PID<amr.num[lv]; PID++)
   } // for (int lv=0; lv<NLEVEL; lv++)


// output cells with high density but lie outside the target region (just to be more cautious)
   if ( UseMaxRhoPos_R < MaxRadius )
   for (int lv=0; lv<NLEVEL; lv++)
   {
      scale = (double)amr.scale[lv];

      for (int PID=0; PID<amr.num[lv]; PID++)
      {
         if ( amr.patch[lv][PID]->fluid == NULL  ||  amr.patch[lv][PID]->son != -1 )   continue;

         for (int k=0; k<PATCH_SIZE; k++) {  z   = amr.patch[lv][PID]->corner[2] + (k+0.5)*scale;
                                             dz1 = z - Center    [2];
                                             dz2 = z - Center_Map[2];
                                             dz  = ( fabs(dz1) <= fabs(dz2) ) ? dz1 : dz2;
         for (int j=0; j<PATCH_SIZE; j++) {  y   = amr.patch[lv][PID]->corner[1] + (j+0.5)*scale;
                                             dy1 = y - Center    [1];
                                             dy2 = y - Center_Map[1];
                                             dy  = ( fabs(dy1) <= fabs(dy2) ) ? dy1 : dy2;
         for (int i=0; i<PATCH_SIZE; i++) {  x   = amr.patch[lv][PID]->corner[0] + (i+0.5)*scale;
                                             dx1 = x - Center    [0];
                                             dx2 = x - Center_Map[0];
                                             dx  = ( fabs(dx1) <= fabs(dx2) ) ? dx1 : dx2;

            if ( UseMaxRhoPos_Par )    Rho = amr.patch[lv][PID]->par_dens   [k][j][i];
            else                       Rho = amr.patch[lv][PID]->fluid[DENS][k][j][i];

            Radius = sqrt( dx*dx + dy*dy + dz*dz );

            if ( Rho > MinRho  &&  Radius > UseMaxRhoPos_R )
            {
               if ( PeakOutside_First )
               {
                  fprintf( stderr, "   WARNING : high density (> %13.7e) found outside the target region for DumpID %2d\n",
                           MinRho, DumpID );
                  fprintf( stderr, "   (%2s, %13s, %13s, %13s, %13s, %13s)\n", "Lv", "x", "y", "z", "Radius", "Density" );

                  PeakOutside_First = false;
               }

               fprintf( stderr, "   (%2d, %13.7e, %13.7e, %13.7e, %13.7e, %13.7e)\n",
                        lv, x*dh_min, y*dh_min, z*dh_min, Radius*dh_min, Rho );

//             record the information of the maximum density lying outside the target region
               if ( Rho > PeakOutside_MaxRho )
               {
                  PeakOutside_Lv     = lv;
                  PeakOutside_Pos[0] = x;
                  PeakOutside_Pos[1] = y;
                  PeakOutside_Pos[2] = z;
                  PeakOutside_R      = Radius;
                  PeakOutside_MaxRho = Rho;
               }
            } // if ( Rho > MinRho  &&  Radius > UseMaxRhoPos_R )
         }}} // i, j, k
      } // for (int PID=0; PID<amr.num[lv]; PID++)
   } // for (int lv=0; lv<NLEVEL; lv++)

   if ( PeakOutside_First )
      fprintf( stdout, "   GOOD : No high density points (> %13.7e) were found outside the target region\n", MinRho );

   else
   {
      fprintf( stderr, "   The one with peak density:\n" );
      fprintf( stderr, "   (%2d, %13.7e, %13.7e, %13.7e, %13.7e, %13.7e)\n",
               PeakOutside_Lv, PeakOutside_Pos[0]*dh_min, PeakOutside_Pos[1]*dh_min, PeakOutside_Pos[2]*dh_min,
               PeakOutside_R*dh_min, PeakOutside_MaxRho );
   }


// sorting
   real   TempRho;
   double TempPos[3];
   int  TempLv;
   for (int n=0; n<AveN-1; n++)
   for (int m=0; m<AveN-1-n; m++)
   {
      if ( MaxRho[m] < MaxRho[m+1] )
      {
         TempRho = MaxRho   [m];
         TempLv  = MaxRho_Lv[m];
         for (int dim=0; dim<3; dim++)    TempPos        [dim] = MaxRho_Pos[m  ][dim];

         MaxRho   [m] = MaxRho   [m+1];
         MaxRho_Lv[m] = MaxRho_Lv[m+1];
         for (int dim=0; dim<3; dim++)    MaxRho_Pos[m  ][dim] = MaxRho_Pos[m+1][dim];

         MaxRho   [m+1] = TempRho;
         MaxRho_Lv[m+1] = TempLv;
         for (int dim=0; dim<3; dim++)    MaxRho_Pos[m+1][dim] = TempPos        [dim];
      }
   }


// check if we do find some cells with non-zero densities
   if ( MaxRho[0] <= 0.0 )
      Aux_Error( ERROR_INFO, "Maximum density found = %14.7e <= 0.0 !!\n", MaxRho[0] );


// set the new sphere center as the center of mass of AveN highest-density cells
   double PosSum[3] = { 0.0, 0.0, 0.0 };
   double TotalMass = 0.0;
   double dv, Mass;

   for (int t=0; t<AveN; t++)
   {
      dv   = pow( amr.dh[ MaxRho_Lv[t] ], 3.0 );
      Mass = MaxRho[t]*dv;

      for (int dim=0; dim<3; dim++)    PosSum[dim] += MaxRho_Pos[t][dim]*Mass;

      TotalMass += Mass;
   }

   for (int dim=0; dim<3; dim++)    Center[dim] = PosSum[dim] / TotalMass;


// print information
   printf( "   ===================================================================================\n" );
   printf( "   Maximum density:\n" );
   printf( "   %7s  %2s  (%13s, %13s, %13s)  %13s\n", "Ranking", "Lv", "x", "y", "z", "Density" );

   for (int t=0; t<AveN; t++)
      printf( "   %7d  %2d  (%13.7e, %13.7e, %13.7e)  %13.7e\n",
              t, MaxRho_Lv[t], MaxRho_Pos[t][0]*dh_min, MaxRho_Pos[t][1]*dh_min, MaxRho_Pos[t][2]*dh_min,
              MaxRho[t] );

   printf( "\n   New center: (%13.7e, %13.7e, %13.7e)\n", Center[0]*dh_min, Center[1]*dh_min, Center[2]*dh_min );
   printf( "   ===================================================================================\n" );


// remove the allocated patches (because the candidate box may be wrong due to the incorrect sphere center)
   int NPatch;

   for (int lv=NLEVEL-1; lv>=0; lv--)
   {
      NPatch = amr.num[lv];

      for (int PID=0; PID<NPatch; PID++)  amr.pdelete( lv, PID );
   }


// reload the patches with the new sphere center
   LoadData();


   cout << "SetMaxRhoPos ... done" << endl;

} // FUNCTION : SetMaxRhoPos



//-------------------------------------------------------------------------------------------------------
// Function    :  TakeNote
// Description :  Output the simulation parameters
//-------------------------------------------------------------------------------------------------------
void TakeNote( int argc, char **argv )
{

   cout << "TakeNote ... " << endl;


   const double dh_min = amr.dh[NLEVEL-1];

   printf( "======================================================================================\n"   );

// record the command-line options
   fprintf( stdout, "\n" );
   fprintf( stdout, "Command-line arguments :\n" );
   for (int v=0; v<argc; v++)    fprintf( stdout, " %s", argv[v] );
   fprintf( stdout, "\n\n" );

#  if   ( MODEL == HYDRO )
   printf( "MODEL            = %14s\n",     "HYDRO"                   );
#  elif ( MODEL == MHD )
   printf( "MODEL            = %14s\n",     "MHD"                     );
#  elif ( MODEL == ELBDM )
   printf( "MODEL            = %14s\n",     "ELBDM"                   );
#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL

#  ifdef FLOAT8
   printf( "FLOAT8           = %14s\n",     "ON"                      );
#  else
   printf( "FLOAT8           = %14s\n",     "OFF"                     );
#  endif

#  ifdef GAMER_DEBUG
   printf( "GAMER_DEBUG      = %14s\n",     "ON"                      );
#  else
   printf( "GAMER_DEBUG      = %14s\n",     "OFF"                     );
#  endif

#  ifdef SUPPORT_HDF5
   printf( "SUPPORT_HDF5     = %14s\n",     "ON"                      );
#  else
   printf( "SUPPORT_HDF5     = %14s\n",     "OFF"                     );
#  endif

   printf( "NCOMP_FLUID      = %14d\n",     NCOMP_FLUID               );
   printf( "NCOMP_PASSIVE    = %14d\n",     NCOMP_PASSIVE             );
   printf( "DumpID           = %14d\n",     DumpID                    );
   printf( "Time             = %14.7e\n",   Time[0]                   );
   printf( "Step             = %14ld\n",    Step                      );
   printf( "NX0_TOT[0]       = %14d\n",     NX0_TOT[0]                );
   printf( "NX0_TOT[1]       = %14d\n",     NX0_TOT[1]                );
   printf( "NX0_TOT[2]       = %14d\n",     NX0_TOT[2]                );
   printf( "NIn              = %14d\n",     NIn                       );
   printf( "NOut             = %14d\n",     NOut                      );
   printf( "Center_x         = %14.7e\n",   Center[0]*dh_min          );
   printf( "Center_y         = %14.7e\n",   Center[1]*dh_min          );
   printf( "Center_z         = %14.7e\n",   Center[2]*dh_min          );
   printf( "MappedCenter_x   = %14.7e\n",   Center_Map[0]*dh_min      );
   printf( "MappedCenter_y   = %14.7e\n",   Center_Map[1]*dh_min      );
   printf( "MappedCenter_z   = %14.7e\n",   Center_Map[2]*dh_min      );
   printf( "Radius           = %14.7e\n",   MaxRadius*dh_min          );
   printf( "UseMaxRhoPos     = %14d\n",     UseMaxRhoPos              );
   printf( "UseMaxRhoPos_R   = %14.7e\n",   UseMaxRhoPos_R*dh_min     );
   printf( "UseMaxRhoPos_Par = %14d\n",     UseMaxRhoPos_Par          );
   printf( "Periodic         = %14d\n",     Periodic                  );
   printf( "InputScale       = %14d\n",     InputScale                );
   printf( "LogBin           = %14.7e\n",   LogBin                    );
   printf( "GetNShell        = %14.7e\n",   GetNShell                 );
   printf( "NShell           = %14d\n",     NShell                    );
   printf( "ShellWidth       = %14.7e\n",   ShellWidth*dh_min         );
   printf( "RhoThres         = %14.7e\n",   RhoThres                  );
   printf( "OutputPot        =  %s\n",      (OutputPot)?"YES":"NO"    );
   printf( "OutputParDens    =  %d\n",      OutputParDens             );
   printf( "OutputCenter     =  %s\n",      (OutputCenter)?"YES":"NO" );
   printf( "NeedGhost        =  %s\n",      (NeedGhost)?"YES":"NO"    );
   printf( "UseTree          =  %s\n",      (UseTree)?"YES":"NO"      );
   printf( "FileName_Tree    =  %s\n",      FileName_Tree             );
   printf( "IntScheme        =  %s\n",      ( IntScheme == INT_MINMOD3D ) ? "MINMOD3D" :
                                            ( IntScheme == INT_MINMOD1D ) ? "MINMOD1D" :
                                            ( IntScheme == INT_VANLEER  ) ? "VANLEER"  :
                                            ( IntScheme == INT_CQUAD    ) ? "CQUAD"    :
                                            ( IntScheme == INT_QUAD     ) ? "QUAD"     :
                                            ( IntScheme == INT_CQUAR    ) ? "CQUAR"    :
                                            ( IntScheme == INT_QUAR     ) ? "QUAR"     :
                                            "UNKNOWN" );
#  if ( MODEL == ELBDM )
   printf( "ELBDM_IntPhase   =  %s\n",      (ELBDM_IntPhase)?"YES":"NO" );
   printf( "ELBDM_GetVir     =  %s\n",      (ELBDM_GetVir  )?"YES":"NO" );
   printf( "ELBDM_ETA        = %14.7e\n",   ELBDM_ETA );
#  endif
   printf( "INT_MONO_COEFF   = %14.7e\n",   INT_MONO_COEFF );
   printf( "GetAvePot        =  %s\n",      (GetAvePot)?"YES":"NO"      );
   if ( GetAvePot )
   printf( "NewtonG          = %14.7e\n",   NewtonG                     );
   printf( "======================================================================================\n"   );


   cout << "TakeNote ... done" << endl;

} // FUNCTION : TakeNote



//-------------------------------------------------------------------------------------------------------
// Function    :  LoadTree
// Description :  Load the octree structure and verify it is consistent with the data file
//-------------------------------------------------------------------------------------------------------
void LoadTree()
{

   fprintf( stdout, "%s ...\n", __FUNCTION__ ); fflush( stdout );


// 1. load the tree
   int     nlv, data_id, nblock_tot;
   long    step;
   double  time, box_size[3];
   int    *nblock = NULL;
   double *dh     = NULL;

   FILE *FileTree = fopen( FileName_Tree, "rb" );

   if ( FileTree == NULL )
   {
      fprintf( stderr, "ERROR : the input tree file \"%s\" does not exist !!\n", FileName_Tree );
      exit( 1 ) ;
   }

   fread( &nlv,        sizeof(int),    1,   FileTree );
   fread( &data_id,    sizeof(int),    1,   FileTree );
   fread( &nblock_tot, sizeof(int),    1,   FileTree );
   fread( &step,       sizeof(long),   1,   FileTree );
   fread( &time,       sizeof(double), 1,   FileTree );
   fread(  box_size,   sizeof(double), 3,   FileTree );

   nblock = new int    [nlv];
   dh     = new double [nlv];
   fread(  nblock,     sizeof(int),    nlv, FileTree );
   fread(  dh,         sizeof(double), nlv, FileTree );

   tree = new tree_t( nlv, data_id, time, step, nblock, dh, box_size );

   fread( tree->block[0], sizeof(block_t), tree->nblock_tot, FileTree );

   fclose( FileTree );

   delete [] nblock;
   delete [] dh;


// 2. verify the loaded tree is consistent with the data file
   FILE *FileGAMER = fopen( FileName_In, "rb" );

   if ( FileGAMER == NULL )
   {
      fprintf( stderr, "ERROR : the input data file \"%s\" does not exist !!\n", FileName_In );
      exit( 1 ) ;
   }

   const int Offset_Float8 = 256 + 5*sizeof(int) + 9*sizeof(bool);
   const int Offset_NLv    = Offset_Float8 + 1*sizeof(int) + 6*sizeof(bool);
   long    FormatVersion, HeaderSize, CheckCode_Ref, Step, CheckCode;
   int     NLv, DataID;
   bool    Float8;
   double *Time   = NULL;
   int    *NBlock = NULL;

// load
   fread( &FormatVersion, sizeof(long),     1, FileGAMER );
   fread( &HeaderSize,    sizeof(long),     1, FileGAMER );
   fread( &CheckCode_Ref, sizeof(long),     1, FileGAMER );

   fseek( FileGAMER, Offset_Float8, SEEK_SET );
   fread( &Float8,        sizeof(bool),     1, FileGAMER );

   fseek( FileGAMER, Offset_NLv,    SEEK_SET );
   fread( &NLv,           sizeof(int),      1, FileGAMER );

   fseek( FileGAMER, HeaderSize,    SEEK_SET );
   fread( &CheckCode,     sizeof(long),     1, FileGAMER );

   Time   = new double [NLv];
   NBlock = new int    [NLv];

   fread( &DataID,        sizeof(int),      1, FileGAMER );
   fread( Time,           sizeof(double), NLv, FileGAMER );
   fread( &Step,          sizeof(long),     1, FileGAMER );
   fread( NBlock,         sizeof(int),    NLv, FileGAMER );

// check
   fprintf( stdout, "   File format = %ld\n", FormatVersion );

   if ( FormatVersion < 1200 )
   {
      fprintf( stderr, "ERROR : unsupported data format version (only support version >= 1200) !!\n" );
      exit( 1 ) ;
   }

   if ( CheckCode != CheckCode_Ref )
   {
      fprintf( stderr, "ERROR : incorrect check code in the data file (input %ld <-> expeect %ld) !!\n",
               CheckCode, CheckCode_Ref );
      exit( 1 );
   }

#  ifdef FLOAT8
   if ( !Float8 )
   {
      fprintf( stderr, "ERROR : %s : data file (%s) != runtime (%s) !!\n", "FLOAT8", "OFF", "ON" );
      exit( 1 ) ;
   }
#  else
   if (  Float8 )
   {
      fprintf( stderr, "ERROR : %s : data file (%s) != runtime (%s) !!\n", "FLOAT8", "ON", "OFF" );
      exit( 1 ) ;
   }
#  endif

   if ( NLv != tree->nlv )
   {
      fprintf( stderr, "ERROR : inconsistent \"%s\" (data: %d <-> tree: %d) !!\n",
               "NLv", NLv, tree->nlv );
      exit( 1 );
   }

   if ( DataID != tree->data_id )
   {
      fprintf( stderr, "ERROR : inconsistent \"%s\" (data: %d <-> tree: %d) !!\n",
               "DataID", DataID, tree->data_id );
      exit( 1 );
   }

   if ( Time[0] != tree->time )
   {
      fprintf( stderr, "ERROR : inconsistent \"%s\" (data: %13.7e <-> tree: %13.7e) !!\n",
               "Time", Time[0], tree->time );
      exit( 1 );
   }

   if ( Step != tree->step )
   {
      fprintf( stderr, "ERROR : inconsistent \"%s\" (data: %ld <-> tree: %ld) !!\n",
               "Step", Step, tree->step );
      exit( 1 );
   }

   for (int lv=0; lv<NLv; lv++)
   {
      if ( NBlock[lv] != tree->nblock[lv] )
      {
         fprintf( stderr, "ERROR : inconsistent \"%s\" at level %d (data: %d <-> tree: %d) !!\n",
                  "NBlock", lv, NBlock[lv], tree->nblock[lv] );
         exit( 1 );
      }
   }

   fclose( FileGAMER );

   delete [] Time;
   delete [] NBlock;

   fprintf( stdout, "   Check ... passed\n" );


   fprintf( stdout, "%s ... done\n", __FUNCTION__ ); fflush( stdout );

} // FUNCTION : LoadTree



//-------------------------------------------------------------------------------------------------------
// Function    :  GetR
// Description :  Get the radius and width of the target shell
//-------------------------------------------------------------------------------------------------------
void GetR( const int n, double &r, double &dr )
{

   const double dr0 = ShellWidth*amr.dh[NLEVEL-1];

   if ( LogBin > 1.0 )
   {
      r = dr0*pow( (double)LogBin, n-0.5 );

      if ( n == 0)   dr = dr0;
      else           dr = dr0*(  pow( (double)LogBin, n ) - pow( (double)LogBin, n-1 )  );
   }

   else
   {
      r  = (n+0.5)*dr0;
      dr = dr0;
   }

} // FUNCTION : GetR



//-------------------------------------------------------------------------------------------------------
// Function    :  main
// Description :
//-------------------------------------------------------------------------------------------------------
int main( int argc, char ** argv )
{

   ReadOption( argc, argv );

   if ( UseTree )    LoadTree();

   LoadData();

   if ( UseMaxRhoPos > 0 )    SetMaxRhoPos( UseMaxRhoPos );

   TakeNote( argc, argv );

   CheckParameter();

   if ( Mode_MaxRho )   GetMaxRho();

   if ( Mode_ShellAve )
   {
      Init_ShellAve();

      ShellAverage();

      GetRMS();

      Output_ShellAve();
   }

   End();


   cout << "Program terminated successfully" << endl;

   return 0;

} // FUNCTION : main
