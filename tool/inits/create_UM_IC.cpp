#include <cstdio>
#include <cmath>


//#define FLOAT8

#ifdef FLOAT8
typedef double real;
#else
typedef float  real;
#endif


#define SQR( a )     ( (a)*(a) )
#define PATCH_SIZE   8
#define PS2          ( 2*PATCH_SIZE )



int main()
{

// UM_IC parameters
   const char  *Filename          = "UM_IC";
   const int    OPT__UM_IC_LEVEL  = 1;
   const int    OPT__UM_IC_NLEVEL = 3;
   const int    OPT__UM_IC_NVAR   = 5;
   const int    RefineRegion[OPT__UM_IC_NLEVEL-1][6] = {  { 2, 2, 2, 2, 2, 2 },
                                                          { 1, 1, 1, 1, 1, 1 }  };

// other parameters
   const int    N0       = 64;   // grid size on level "OPT__UM_IC_LEVEL"
   const double BOX_SIZE = 1.0;  // simulation box size


// convert RefineRegion[] to the number of cells to be skipped on each side
   int NCellSkip[OPT__UM_IC_NLEVEL][6];

   for (int t=0; t<OPT__UM_IC_NLEVEL; t++)
   for (int s=0; s<6; s++)
   {
      if ( t == 0 )  NCellSkip[t][s] = 0;
      else           NCellSkip[t][s] = RefineRegion[t-1][s]*PS2 + NCellSkip[t-1][s]*2;
   }


// set the left-edge coordinates on each level
   double dh[OPT__UM_IC_NLEVEL], EdgeL[OPT__UM_IC_NLEVEL][3];

   for (int t=0; t<OPT__UM_IC_NLEVEL; t++)
   {
      dh[t] = BOX_SIZE / ( N0*(1<<t) );

      for (int d=0; d<3; d++)    EdgeL[t][d] = NCellSkip[t][2*d]*dh[t];
   }


// determine the number of cells on each level along each direction
   long UM_Size[OPT__UM_IC_NLEVEL][3];

   for (int t=0; t<OPT__UM_IC_NLEVEL; t++)
   for (int d=0; d<3; d++)
      UM_Size[t][d] = N0*(1<<t) - NCellSkip[t][2*d] - NCellSkip[t][2*d+1];


// allocate data
   real *data[OPT__UM_IC_NLEVEL];

   for (int t=0; t<OPT__UM_IC_NLEVEL; t++)
      data[t] = new real [ OPT__UM_IC_NVAR*UM_Size[t][0]*UM_Size[t][1]*UM_Size[t][2] ];


// assign data
   const double Gamma  = 5.0 / 3.0;
   const double Amp    = 1.0e4;
   const double Sigma  = 1.0e-1;
   const double BoxCen = 0.5*BOX_SIZE;

   double x, y, z, dx, dy, dz, Dens, Pres;

   for (int t=0; t<OPT__UM_IC_NLEVEL; t++)
   {
      real (*data_lv)[ UM_Size[t][2] ][ UM_Size[t][1] ][ UM_Size[t][0] ]
         = ( real (*)[ UM_Size[t][2] ][ UM_Size[t][1] ][ UM_Size[t][0] ] )data[t];

      for (int k=0; k<UM_Size[t][2]; k++)  {  z = EdgeL[t][2] + (k+0.5)*dh[t];  dz = z - BoxCen;
      for (int j=0; j<UM_Size[t][1]; j++)  {  y = EdgeL[t][1] + (j+0.5)*dh[t];  dy = y - BoxCen;
      for (int i=0; i<UM_Size[t][0]; i++)  {  x = EdgeL[t][0] + (i+0.5)*dh[t];  dx = x - BoxCen;

         Dens = Amp*exp( -0.5*( SQR(dx) + SQR(dy) + SQR(dz) )/SQR(Sigma) );
         Pres = 2.0*Dens;  // arbitrary

         data_lv[0][k][j][i] = Dens;
         data_lv[1][k][j][i] = 0.0;    // zero velocity
         data_lv[2][k][j][i] = 0.0;
         data_lv[3][k][j][i] = 0.0;
         data_lv[4][k][j][i] = Pres / (Gamma-1.0);

      }}}
   }


// output file
   FILE *File = fopen( Filename, "wb" );

   for (int t=0; t<OPT__UM_IC_NLEVEL; t++)
      fwrite( data[t], sizeof(real), OPT__UM_IC_NVAR*UM_Size[t][0]*UM_Size[t][1]*UM_Size[t][2], File );

   fclose( File );


   for (int t=0; t<OPT__UM_IC_NLEVEL; t++)   delete [] data[t];

} // FUNCTION : main
