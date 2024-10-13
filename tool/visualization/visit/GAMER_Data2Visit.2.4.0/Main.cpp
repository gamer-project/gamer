#include "Data2Visit.h"

AMR_t    amr;
char    *FileName_In;
int      NMesh[NLEVEL];      // number of patches (actually allocated) at each level
double   Time[NLEVEL];
int      DumpID;

#if   ( MODEL == HYDRO )
double   GAMMA;

#elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#elif ( MODEL == ELBDM )
double   ELBDM_ETA;

#else
#  error : ERROR : unsupported MODEL !!
#endif // MODEL

double PhySize  [3]   = { NULL_VALUE, NULL_VALUE, NULL_VALUE }; // physical size of the targeted region
double PhyCorner[3]   = { NULL_VALUE, NULL_VALUE, NULL_VALUE }; // physical corner coord. of the targeted region
int    SizeScale[3]   = { NULL_VALUE, NULL_VALUE, NULL_VALUE }; // grid indices of the targeted region
int    CornerScale[3] = { NULL_VALUE, NULL_VALUE, NULL_VALUE }; // corner grid indices of the targeted region
int    TargetMode     = NULL_VALUE;                             // (0,1,2,3) <-> (xy, yz, xz, 3D)
double TPhySize  [3]  = { NULL_VALUE, NULL_VALUE, NULL_VALUE }; // truncated physical size
double TPhyCorner[3]  = { NULL_VALUE, NULL_VALUE, NULL_VALUE }; // truncated physical corner coord.
bool   SepAllLv       = false;                                  // separate data at different levels
#if ( MODEL == HYDRO )
bool   OutputPres     = false;                                  // output pressure
bool   OutputVort     = false;                                  // output vorticity
#endif
bool   OutputPot      = false;                                  // output potential
bool   InputScale     = false;                                  // input cell scales instead of coordinates




//-------------------------------------------------------------------------------------------------------
// Function    :  CreateSilo
// Description :  Create the silo files
//-------------------------------------------------------------------------------------------------------
void CreateSilo()
{

   WriteBox();

   WriteLevel();

   WriteRoot();

} // FUNCTION : CreateSilo



//-------------------------------------------------------------------------------------------------------
// Function    :  WriteBox
// Description :  Create the silo file for the outline of the entire target surface
//-------------------------------------------------------------------------------------------------------
void WriteBox()
{

   cout << "WriteBox ... " << flush;


   DBfile *dbfile = NULL;
   dbfile = DBCreate( "Box.silo", DB_CLOBBER, DB_LOCAL, NULL, DB_HDF5 );

   int    ndims       = 3;
   int    node_dims[] = {2,2,2};
   double x[2]        = { TPhyCorner[0], TPhyCorner[0]+TPhySize[0] };
   double y[2]        = { TPhyCorner[1], TPhyCorner[1]+TPhySize[1] };
   double z[2]        = { TPhyCorner[2], TPhyCorner[2]+TPhySize[2] };
   double *coords[]   = {x, y, z};

   DBPutQuadmesh( dbfile, "Box", NULL, coords, node_dims, ndims, DB_DOUBLE, DB_COLLINEAR, NULL );

   DBClose( dbfile );


   cout << "done" << endl;

} // FUNCTION : WriteBox



//-------------------------------------------------------------------------------------------------------
// Function    :  WriteLevel
// Description :  Create the silo file for each level
//-------------------------------------------------------------------------------------------------------
void WriteLevel()
{

   cout << "WriteLevel ... " << endl;


   const int PS = PATCH_SIZE;    // patch size
   int Size[3];                  // targeted number of grids in one patch

   switch ( TargetMode )
   {
      case 0:
         Size[0] = PS;  Size[1] = PS;  Size[2] = 1;
         break;
      case 1:
         Size[0] = 1;   Size[1] = PS;  Size[2] = PS;
         break;
      case 2:
         Size[0] = PS;  Size[1] = 1;   Size[2] = PS;
         break;
      case 3:
         Size[0] = PS;  Size[1] = PS;  Size[2] = PS;
         break;
   }


   int ndims              = 3;
   int node_dims      [3] = { Size[0]+1, Size[1]+1, Size[2]+1 };
   int node_dims_patch[3] = { 2, 2, 2 };
   int zone_dims      [3] = { Size[0], Size[1], Size[2] };

   int    i1, i2, j1, j2, k1, k2;
   double x[ Size[0]+1 ], y[ Size[1]+1 ] , z[ Size[2]+1 ];
   double x_patch[2], y_patch[2], z_patch[2];

   double *coords      [3] = { x, y, z };
   double *coords_patch[3] = { x_patch, y_patch, z_patch };

#  if   ( MODEL == HYDRO )
   float *Rho  = new float [ Size[0]*Size[1]*Size[2] ];
   float *Vx   = new float [ Size[0]*Size[1]*Size[2] ];
   float *Vy   = new float [ Size[0]*Size[1]*Size[2] ];
   float *Vz   = new float [ Size[0]*Size[1]*Size[2] ];
   float *Egy  = new float [ Size[0]*Size[1]*Size[2] ];
   float *Pres = new float [ Size[0]*Size[1]*Size[2] ];
   float *Wx   = new float [ Size[0]*Size[1]*Size[2] ];
   float *Wy   = new float [ Size[0]*Size[1]*Size[2] ];
   float *Wz   = new float [ Size[0]*Size[1]*Size[2] ];

   char  *Vnames[3] = { "Vx", "Vy", "Vz" };
   char  *Wnames[3] = { "Wx", "Wy", "Wz" };
   float *Vel   [3] = { Vx, Vy, Vz };
   float *Vor   [3] = { Wx, Wy, Wz };

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   float *Dens = new float [ Size[0]*Size[1]*Size[2] ];
   float *Real = new float [ Size[0]*Size[1]*Size[2] ];
   float *Imag = new float [ Size[0]*Size[1]*Size[2] ];
   float *VelX = new float [ Size[0]*Size[1]*Size[2] ];
   float *VelY = new float [ Size[0]*Size[1]*Size[2] ];
   float *VelZ = new float [ Size[0]*Size[1]*Size[2] ];
   float *VelT = new float [ Size[0]*Size[1]*Size[2] ];  // thermal velocity (internal energy density = 0.5*rho*VelT^2)
   float *VelE = new float [ Size[0]*Size[1]*Size[2] ];  // effective total velocity = sqrt(Vt^2+V^2)

   char  *Vnames[3] = { "Vx", "Vy", "Vz" };
   float *Vel   [3] = { VelX, VelY, VelZ };

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL

   float *Pot  = ( OutputPot ) ? ( new float [ Size[0]*Size[1]*Size[2] ] ) : NULL;


   for (int lv=0; lv<NLEVEL; lv++)
   {
      cout << "   Level " << lv << " ... " << flush;


      const double dh  = amr.dh   [lv];
      const int  scale = amr.scale[lv];

      NMesh[lv] = 0;

      DBfile *dbfile = NULL;
      char filename[100];
      sprintf( filename, "Level_%d.silo", lv );

      dbfile = DBCreate( filename, DB_CLOBBER, DB_LOCAL, NULL, DB_HDF5 );


      for (int PID=0; PID<amr.num[lv]; PID++)
      {

//       set the coordinates for each grid cell and the range of target cells for each patch
         switch ( TargetMode )
         {
            case 0:
               for (int i=0; i<Size[0]+1; i++)     x[i] = ( amr.patch[lv][PID]->corner[0]/scale + i )*dh;
               for (int i=0; i<Size[1]+1; i++)     y[i] = ( amr.patch[lv][PID]->corner[1]/scale + i )*dh;
               for (int i=0; i<Size[2]+1; i++)     z[i] = TPhyCorner[2] + i*TPhySize[2];

               x_patch[0] = amr.patch[lv][PID]->corner[0]/scale*dh;
               x_patch[1] = amr.patch[lv][PID]->corner[0]/scale*dh + PS*dh;
               y_patch[0] = amr.patch[lv][PID]->corner[1]/scale*dh;
               y_patch[1] = amr.patch[lv][PID]->corner[1]/scale*dh + PS*dh;
               z_patch[0] = z[0];
               z_patch[1] = z[1];

               i1 = 0;
               i2 = PS;
               j1 = 0;
               j2 = PS;
               k1 = ( CornerScale[2] - amr.patch[lv][PID]->corner[2] ) / scale;
               k2 = k1 + 1;

               break;


            case 1:
               for (int i=0; i<Size[0]+1; i++)     x[i] = TPhyCorner[0] + i*TPhySize[0];
               for (int i=0; i<Size[1]+1; i++)     y[i] = ( amr.patch[lv][PID]->corner[1]/scale + i )*dh;
               for (int i=0; i<Size[2]+1; i++)     z[i] = ( amr.patch[lv][PID]->corner[2]/scale + i )*dh;

               x_patch[0] = x[0];
               x_patch[1] = x[1];
               y_patch[0] = amr.patch[lv][PID]->corner[1]/scale*dh;
               y_patch[1] = amr.patch[lv][PID]->corner[1]/scale*dh + PS*dh;
               z_patch[0] = amr.patch[lv][PID]->corner[2]/scale*dh;
               z_patch[1] = amr.patch[lv][PID]->corner[2]/scale*dh + PS*dh;

               i1 = ( CornerScale[0] - amr.patch[lv][PID]->corner[0] ) / scale;
               i2 = i1 + 1;
               j1 = 0;
               j2 = PS;
               k1 = 0;
               k2 = PS;

               break;


            case 2:
               for (int i=0; i<Size[0]+1; i++)     x[i] = ( amr.patch[lv][PID]->corner[0]/scale + i )*dh;
               for (int i=0; i<Size[1]+1; i++)     y[i] = TPhyCorner[1] + i*TPhySize[1];
               for (int i=0; i<Size[2]+1; i++)     z[i] = ( amr.patch[lv][PID]->corner[2]/scale + i )*dh;

               x_patch[0] = amr.patch[lv][PID]->corner[0]/scale*dh;
               x_patch[1] = amr.patch[lv][PID]->corner[0]/scale*dh + PS*dh;
               y_patch[0] = y[0];
               y_patch[1] = y[1];
               z_patch[0] = amr.patch[lv][PID]->corner[2]/scale*dh;
               z_patch[1] = amr.patch[lv][PID]->corner[2]/scale*dh + PS*dh;

               i1 = 0;
               i2 = PS;
               j1 = ( CornerScale[1] - amr.patch[lv][PID]->corner[1] ) / scale;
               j2 = j1 + 1;
               k1 = 0;
               k2 = PS;

               break;


            case 3:
               for (int i=0; i<Size[0]+1; i++)     x[i] = ( amr.patch[lv][PID]->corner[0]/scale + i )*dh;
               for (int i=0; i<Size[1]+1; i++)     y[i] = ( amr.patch[lv][PID]->corner[1]/scale + i )*dh;
               for (int i=0; i<Size[2]+1; i++)     z[i] = ( amr.patch[lv][PID]->corner[2]/scale + i )*dh;

               x_patch[0] = amr.patch[lv][PID]->corner[0]/scale*dh;
               x_patch[1] = amr.patch[lv][PID]->corner[0]/scale*dh + PS*dh;
               y_patch[0] = amr.patch[lv][PID]->corner[1]/scale*dh;
               y_patch[1] = amr.patch[lv][PID]->corner[1]/scale*dh + PS*dh;
               z_patch[0] = amr.patch[lv][PID]->corner[2]/scale*dh;
               z_patch[1] = amr.patch[lv][PID]->corner[2]/scale*dh + PS*dh;

               i1 = 0;
               i2 = PS;
               j1 = 0;
               j2 = PS;
               k1 = 0;
               k2 = PS;

               break;
         } // switch ( TargetMode )


//       create patch outline
         char PatchName[20];
         sprintf( PatchName, "Patch_%d_%d", lv, NMesh[lv] );
         DBPutQuadmesh( dbfile, PatchName, NULL, coords_patch, node_dims_patch, ndims, DB_DOUBLE, DB_COLLINEAR,
                        NULL );


//       create mesh
         char MeshName[20];
         sprintf( MeshName, "Mesh_%d_%d", lv, NMesh[lv] );
         DBPutQuadmesh( dbfile, MeshName, NULL, coords, node_dims, ndims, DB_DOUBLE, DB_COLLINEAR, NULL );


//       copy variables into 1-D arrays
#        if   ( MODEL == HYDRO )
         const real Gamma_m1 = GAMMA - 1.0;
         const real (* Rho_ptr)[PATCH_SIZE][PATCH_SIZE] = amr.patch[lv][PID]->fluid[DENS];
         const real (*MomX_ptr)[PATCH_SIZE][PATCH_SIZE] = amr.patch[lv][PID]->fluid[MOMX];
         const real (*MomY_ptr)[PATCH_SIZE][PATCH_SIZE] = amr.patch[lv][PID]->fluid[MOMY];
         const real (*MomZ_ptr)[PATCH_SIZE][PATCH_SIZE] = amr.patch[lv][PID]->fluid[MOMZ];
         const real (* Egy_ptr)[PATCH_SIZE][PATCH_SIZE] = amr.patch[lv][PID]->fluid[ENGY];

         real PxVy, PxVz, PyVx, PyVz, PzVx, PzVy;

#        elif ( MODEL == MHD )
#        warning : WAIT MHD !!!

#        elif ( MODEL == ELBDM )
         const real _Eta = 1.0 / ELBDM_ETA;
         const real (*Dens_ptr)[PATCH_SIZE][PATCH_SIZE] = amr.patch[lv][PID]->fluid[DENS];
         const real (*Real_ptr)[PATCH_SIZE][PATCH_SIZE] = amr.patch[lv][PID]->fluid[REAL];
         const real (*Imag_ptr)[PATCH_SIZE][PATCH_SIZE] = amr.patch[lv][PID]->fluid[IMAG];

         real GradR[3], GradI[3], GradF[3], _Dens;    // F^2 = density = R^2 + I^2

#        else
#        error : ERROR : unsupported MODEL !!
#        endif // MODEL

         const real (*Pot_ptr)[PATCH_SIZE][PATCH_SIZE] = ( OutputPot ) ? amr.patch[lv][PID]->pot : NULL;


         int  ii, jj, kk, ID, dID_m, dID_p;
         real _dh;

         for (int k=k1; k<k2; k++)     {  kk = k - k1;
         for (int j=j1; j<j2; j++)     {  jj = j - j1;
         for (int i=i1; i<i2; i++)     {  ii = i - i1;

            ID = kk*Size[1]*Size[0] + jj*Size[0] + ii;

#           if   ( MODEL == HYDRO )
            Rho[ID] =  Rho_ptr[k][j][i];
            Vx [ID] = MomX_ptr[k][j][i] / Rho[ID];
            Vy [ID] = MomY_ptr[k][j][i] / Rho[ID];
            Vz [ID] = MomZ_ptr[k][j][i] / Rho[ID];
            Egy[ID] =  Egy_ptr[k][j][i];

            if ( OutputPres )
               Pres[ID] = Gamma_m1*( Egy[ID] - 0.5*Rho[ID]*( Vx[ID]*Vx[ID] + Vy[ID]*Vy[ID] + Vz[ID]*Vz[ID] ) );

            if ( OutputVort )
            {
//             partial x
               switch ( i )
               {
                  case 0            : dID_m =  0;   dID_p = 1;   _dh = 1.0/dh;   break;
                  case PATCH_SIZE-1 : dID_m = -1;   dID_p = 0;   _dh = 1.0/dh;   break;
                  default           : dID_m = -1;   dID_p = 1;   _dh = 0.5/dh;   break;
               }

               PxVy = _dh*( MomY_ptr[k][j][i+dID_p] / Rho_ptr[k][j][i+dID_p] -
                            MomY_ptr[k][j][i+dID_m] / Rho_ptr[k][j][i+dID_m] );
               PxVz = _dh*( MomZ_ptr[k][j][i+dID_p] / Rho_ptr[k][j][i+dID_p] -
                            MomZ_ptr[k][j][i+dID_m] / Rho_ptr[k][j][i+dID_m] );

               // partial y
               switch ( j )
               {
                  case 0            : dID_m =  0;   dID_p = 1;   _dh = 1.0/dh;   break;
                  case PATCH_SIZE-1 : dID_m = -1;   dID_p = 0;   _dh = 1.0/dh;   break;
                  default           : dID_m = -1;   dID_p = 1;   _dh = 0.5/dh;   break;
               }

               PyVx = _dh*( MomX_ptr[k][j+dID_p][i] / Rho_ptr[k][j+dID_p][i] -
                            MomX_ptr[k][j+dID_m][i] / Rho_ptr[k][j+dID_m][i] );
               PyVz = _dh*( MomZ_ptr[k][j+dID_p][i] / Rho_ptr[k][j+dID_p][i] -
                            MomZ_ptr[k][j+dID_m][i] / Rho_ptr[k][j+dID_m][i] );

               // partial z
               switch ( k )
               {
                  case 0            : dID_m =  0;   dID_p = 1;   _dh = 1.0/dh;   break;
                  case PATCH_SIZE-1 : dID_m = -1;   dID_p = 0;   _dh = 1.0/dh;   break;
                  default           : dID_m = -1;   dID_p = 1;   _dh = 0.5/dh;   break;
               }

               PzVx = _dh*( MomX_ptr[k+dID_p][j][i] / Rho_ptr[k+dID_p][j][i] -
                            MomX_ptr[k+dID_m][j][i] / Rho_ptr[k+dID_m][j][i] );
               PzVy = _dh*( MomY_ptr[k+dID_p][j][i] / Rho_ptr[k+dID_p][j][i] -
                            MomY_ptr[k+dID_m][j][i] / Rho_ptr[k+dID_m][j][i] );

               Wx[ID] = PyVz - PzVy;
               Wy[ID] = PzVx - PxVz;
               Wz[ID] = PxVy - PyVx;
            } // if ( OutputVort )

#           elif ( MODEL == MHD )
#           warning : WAIT MHD !!!

#           elif ( MODEL == ELBDM )
            Dens[ID] = Dens_ptr[k][j][i];
            Real[ID] = Real_ptr[k][j][i];
            Imag[ID] = Imag_ptr[k][j][i];

//          velocity field
            switch ( i )
            {
               case 0            : dID_m =  0;   dID_p = 1;   _dh = 1.0/dh;   break;
               case PATCH_SIZE-1 : dID_m = -1;   dID_p = 0;   _dh = 1.0/dh;   break;
               default           : dID_m = -1;   dID_p = 1;   _dh = 0.5/dh;   break;
            }

            GradR[0] = _dh*(      Real_ptr[k][j][i+dID_p]  -      Real_ptr[k][j][i+dID_m] );
            GradI[0] = _dh*(      Imag_ptr[k][j][i+dID_p]  -      Imag_ptr[k][j][i+dID_m] );
            GradF[0] = _dh*( sqrt(Dens_ptr[k][j][i+dID_p]) - sqrt(Dens_ptr[k][j][i+dID_m]) );

            // partial y
            switch ( j )
            {
               case 0            : dID_m =  0;   dID_p = 1;   _dh = 1.0/dh;   break;
               case PATCH_SIZE-1 : dID_m = -1;   dID_p = 0;   _dh = 1.0/dh;   break;
               default           : dID_m = -1;   dID_p = 1;   _dh = 0.5/dh;   break;
            }

            GradR[1] = _dh*(      Real_ptr[k][j+dID_p][i]  -      Real_ptr[k][j+dID_m][i] );
            GradI[1] = _dh*(      Imag_ptr[k][j+dID_p][i]  -      Imag_ptr[k][j+dID_m][i] );
            GradF[1] = _dh*( sqrt(Dens_ptr[k][j+dID_p][i]) - sqrt(Dens_ptr[k][j+dID_m][i]) );

            // partial z
            switch ( k )
            {
               case 0            : dID_m =  0;   dID_p = 1;   _dh = 1.0/dh;   break;
               case PATCH_SIZE-1 : dID_m = -1;   dID_p = 0;   _dh = 1.0/dh;   break;
               default           : dID_m = -1;   dID_p = 1;   _dh = 0.5/dh;   break;
            }

            GradR[2] = _dh*(      Real_ptr[k+dID_p][j][i]  -      Real_ptr[k+dID_m][j][i] );
            GradI[2] = _dh*(      Imag_ptr[k+dID_p][j][i]  -      Imag_ptr[k+dID_m][j][i] );
            GradF[2] = _dh*( sqrt(Dens_ptr[k+dID_p][j][i]) - sqrt(Dens_ptr[k+dID_m][j][i]) );

            _Dens    = 1.0 / Dens[ID];
            VelX[ID] = _Eta*( Real[ID]*GradI[0] - Imag[ID]*GradR[0] )*_Dens;
            VelY[ID] = _Eta*( Real[ID]*GradI[1] - Imag[ID]*GradR[1] )*_Dens;
            VelZ[ID] = _Eta*( Real[ID]*GradI[2] - Imag[ID]*GradR[2] )*_Dens;
            VelT[ID] = _Eta*(  sqrt( (GradF[0]*GradF[0] + GradF[1]*GradF[1] + GradF[2]*GradF[2])*_Dens )  );
            VelE[ID] = sqrt( SQR(VelX[ID]) + SQR(VelY[ID]) + SQR(VelZ[ID]) + SQR(VelT[ID]) );

#           else
#           error : ERROR : unsupported MODEL !!
#           endif // MODEL

            if ( OutputPot )
            Pot[ID] = Pot_ptr[k][j][i];

         }}} // i,j/k


//       save variables into mesh
         char VarName[20];

#        if   ( MODEL == HYDRO )
         sprintf( VarName, "Rho_%d_%d", lv, NMesh[lv] );
         DBPutQuadvar1( dbfile, VarName, MeshName, Rho,  zone_dims, ndims, NULL, 0, DB_FLOAT, DB_ZONECENT, NULL );

         sprintf( VarName, "Egy_%d_%d", lv, NMesh[lv] );
         DBPutQuadvar1( dbfile, VarName, MeshName, Egy,  zone_dims, ndims, NULL, 0, DB_FLOAT, DB_ZONECENT, NULL );

         if ( OutputPres )
         {
         sprintf( VarName, "Pre_%d_%d", lv, NMesh[lv] );
         DBPutQuadvar1( dbfile, VarName, MeshName, Pres, zone_dims, ndims, NULL, 0, DB_FLOAT, DB_ZONECENT, NULL );
         }

         sprintf( VarName, "Vel_%d_%d", lv, NMesh[lv] );
         DBPutQuadvar( dbfile, VarName, MeshName, 3, Vnames, Vel, zone_dims, ndims, NULL, 0, DB_FLOAT,
                       DB_ZONECENT, NULL );

         if ( OutputVort )
         {
         sprintf( VarName, "Vor_%d_%d", lv, NMesh[lv] );
         DBPutQuadvar( dbfile, VarName, MeshName, 3, Wnames, Vor, zone_dims, ndims, NULL, 0, DB_FLOAT,
                       DB_ZONECENT, NULL );
         }

#        elif ( MODEL == MHD )
#        warning : WAIT MHD !!!

#        elif ( MODEL == ELBDM )
         sprintf( VarName, "Dens_%d_%d", lv, NMesh[lv] );
         DBPutQuadvar1( dbfile, VarName, MeshName, Dens, zone_dims, ndims, NULL, 0, DB_FLOAT, DB_ZONECENT, NULL );

         sprintf( VarName, "Real_%d_%d", lv, NMesh[lv] );
         DBPutQuadvar1( dbfile, VarName, MeshName, Real, zone_dims, ndims, NULL, 0, DB_FLOAT, DB_ZONECENT, NULL );

         sprintf( VarName, "Imag_%d_%d", lv, NMesh[lv] );
         DBPutQuadvar1( dbfile, VarName, MeshName, Imag, zone_dims, ndims, NULL, 0, DB_FLOAT, DB_ZONECENT, NULL );

         sprintf( VarName, "Vel_%d_%d",  lv, NMesh[lv] );
         DBPutQuadvar( dbfile, VarName, MeshName, 3, Vnames, Vel, zone_dims, ndims, NULL, 0, DB_FLOAT,
                       DB_ZONECENT, NULL );

         sprintf( VarName, "VelT_%d_%d", lv, NMesh[lv] );
         DBPutQuadvar1( dbfile, VarName, MeshName, VelT, zone_dims, ndims, NULL, 0, DB_FLOAT, DB_ZONECENT, NULL );

         sprintf( VarName, "VelE_%d_%d", lv, NMesh[lv] );
         DBPutQuadvar1( dbfile, VarName, MeshName, VelE, zone_dims, ndims, NULL, 0, DB_FLOAT, DB_ZONECENT, NULL );

#        else
#        error : ERROR : unsupported MODEL !!
#        endif // MODEL

         if ( OutputPot )
         {
         sprintf( VarName, "Pot_%d_%d", lv, NMesh[lv] );
         DBPutQuadvar1( dbfile, VarName, MeshName, Pot,  zone_dims, ndims, NULL, 0, DB_FLOAT, DB_ZONECENT, NULL );
         }

         NMesh[lv]++;

      } // for (int PID=0; PID<NMesh; PID++)

      DBClose( dbfile );

      cout << "done" << endl;

   } // for (int lv=0; lv<NLEVEL; lv++)


#  if   ( MODEL == HYDRO )
   delete [] Rho;
   delete [] Vx;
   delete [] Vy;
   delete [] Vz;
   delete [] Egy;
   delete [] Pres;
   delete [] Wx;
   delete [] Wy;
   delete [] Wz;

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   delete [] Dens;
   delete [] Real;
   delete [] Imag;
   delete [] VelX;
   delete [] VelY;
   delete [] VelZ;
   delete [] VelT;
   delete [] VelE;

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL

   if ( OutputPot )
   delete [] Pot;


   cout << "WriteLevel ... done" << endl;

} // FUNCTION : WriteLevel



//-------------------------------------------------------------------------------------------------------
// Function    :  WriteRoot
// Description :  Create the root silo file
//-------------------------------------------------------------------------------------------------------
void WriteRoot()
{

   cout << "WriteRoot ... " << flush;


   DBfile    *dbfile  = DBCreate( "Root.silo", DB_CLOBBER, DB_LOCAL, "Root file of each level", DB_HDF5 );
   DBoptlist *OptList = DBMakeOptlist( 2 );

   DBAddOption( OptList, DBOPT_DTIME, &Time[0] );
   DBAddOption( OptList, DBOPT_CYCLE, &DumpID );

   if ( !SepAllLv ) // create multi-mesh object and integrate data at all levels
   {
      char **patchnames  = NULL;
      char **meshnames   = NULL;
      char **varname_pot = NULL;
      char **varnames0   = NULL, **varnames1 = NULL, **varnames2 = NULL, **varnames3 = NULL, **varnames4 = NULL;
#     if ( MODEL == ELBDM )
      char **varnames5   = NULL;
#     endif
      int nmesh = 0, nvar = 0, line = 0;
      int *patchtypes = NULL, *meshtypes = NULL, *vartypes = NULL;

      for (int lv=0; lv<NLEVEL; lv++)  nmesh += NMesh[lv];
      nvar = nmesh;

      patchnames  = (char **) malloc(nmesh * sizeof(char *));
      meshnames   = (char **) malloc(nmesh * sizeof(char *));
      patchtypes  = (int   *) malloc(nmesh * sizeof(int   ));
      meshtypes   = (int   *) malloc(nmesh * sizeof(int   ));
      varname_pot = (char **) malloc(nvar  * sizeof(char *));
      varnames0   = (char **) malloc(nvar  * sizeof(char *));
      varnames1   = (char **) malloc(nvar  * sizeof(char *));
      varnames2   = (char **) malloc(nvar  * sizeof(char *));
      varnames3   = (char **) malloc(nvar  * sizeof(char *));
      varnames4   = (char **) malloc(nvar  * sizeof(char *));
#     if ( MODEL == ELBDM )
      varnames5   = (char **) malloc(nvar  * sizeof(char *));
#     endif
      vartypes    = (int   *) malloc(nvar  * sizeof(int   ));

      for (int lv=0; lv<NLEVEL; lv++)
      {
         for (int Mesh=0; Mesh<NMesh[lv]; Mesh++)
         {
            char tmp[100];

            sprintf( tmp, "Level_%d.silo:Patch_%d_%d", lv, lv, Mesh );
            patchnames[line] = strdup(tmp);
            patchtypes[line] = DB_QUAD_RECT;

            sprintf( tmp, "Level_%d.silo:Mesh_%d_%d", lv, lv, Mesh );
            meshnames[line] = strdup(tmp);
            meshtypes[line] = DB_QUAD_RECT;

#           if   ( MODEL == HYDRO )
            sprintf( tmp, "Level_%d.silo:Rho_%d_%d",  lv, lv, Mesh );
            varnames0[line] = strdup(tmp);

            sprintf( tmp, "Level_%d.silo:Vel_%d_%d",  lv, lv, Mesh );
            varnames1[line] = strdup(tmp);

            sprintf( tmp, "Level_%d.silo:Egy_%d_%d",  lv, lv, Mesh );
            varnames2[line] = strdup(tmp);

            if ( OutputPres )
            {
            sprintf( tmp, "Level_%d.silo:Pre_%d_%d",  lv, lv, Mesh );
            varnames3[line] = strdup(tmp);
            }

            if ( OutputVort )
            {
            sprintf( tmp, "Level_%d.silo:Vor_%d_%d",  lv, lv, Mesh );
            varnames4[line] = strdup(tmp);
            }

#           elif ( MODEL == MHD )
#           warning : WAIT MHD !!!

#           elif ( MODEL == ELBDM )
            sprintf( tmp, "Level_%d.silo:Dens_%d_%d",  lv, lv, Mesh );
            varnames0[line] = strdup(tmp);

            sprintf( tmp, "Level_%d.silo:Real_%d_%d",  lv, lv, Mesh );
            varnames1[line] = strdup(tmp);

            sprintf( tmp, "Level_%d.silo:Imag_%d_%d",  lv, lv, Mesh );
            varnames2[line] = strdup(tmp);

            sprintf( tmp, "Level_%d.silo:Vel_%d_%d",   lv, lv, Mesh );
            varnames3[line] = strdup(tmp);

            sprintf( tmp, "Level_%d.silo:VelT_%d_%d",  lv, lv, Mesh );
            varnames4[line] = strdup(tmp);

            sprintf( tmp, "Level_%d.silo:VelE_%d_%d",  lv, lv, Mesh );
            varnames5[line] = strdup(tmp);

#           else
#           error : ERROR : unsupported MODEL !!
#           endif // MODEL

            if ( OutputPot )
            {
            sprintf( tmp, "Level_%d.silo:Pot_%d_%d",  lv, lv, Mesh );
            varname_pot[line] = strdup(tmp);
            }

            vartypes[line]  = DB_QUADVAR;

            line++;
         } // for (int Mesh=0; Mesh<NMesh[lv]; Mesh++)
      } // for (int lv=0; lv<NLEVEL; lv++)

      if ( nmesh != 0 )
      {
         char MultiPatchName[30];
         sprintf( MultiPatchName, "%s", "Patch" );
         DBPutMultimesh( dbfile, MultiPatchName, nmesh, patchnames, patchtypes, OptList );

         char MultiMeshName[30];
         sprintf( MultiMeshName, "%s", "Mesh" );
         DBPutMultimesh( dbfile, MultiMeshName, nmesh, meshnames, meshtypes, OptList );

         char MultiVarName[30];
#        if   ( MODEL == HYDRO )
         sprintf( MultiVarName, "%s", "Rho" );
         DBPutMultivar( dbfile, MultiVarName, nvar, varnames0, vartypes, OptList );

         sprintf( MultiVarName, "%s", "Vel" );
         DBPutMultivar( dbfile, MultiVarName, nvar, varnames1, vartypes, OptList );

         sprintf( MultiVarName, "%s", "Egy" );
         DBPutMultivar( dbfile, MultiVarName, nvar, varnames2, vartypes, OptList );

         if ( OutputPres )
         {
         sprintf( MultiVarName, "%s", "Pre" );
         DBPutMultivar( dbfile, MultiVarName, nvar, varnames3, vartypes, OptList );
         }

         if ( OutputVort )
         {
         sprintf( MultiVarName, "%s", "Vor" );
         DBPutMultivar( dbfile, MultiVarName, nvar, varnames4, vartypes, OptList );
         }

#        elif ( MODEL == MHD )
#        warning : WAIT MHD !!!

#        elif ( MODEL == ELBDM )
         sprintf( MultiVarName, "%s", "Dens" );
         DBPutMultivar( dbfile, MultiVarName, nvar, varnames0, vartypes, OptList );

         sprintf( MultiVarName, "%s", "Real" );
         DBPutMultivar( dbfile, MultiVarName, nvar, varnames1, vartypes, OptList );

         sprintf( MultiVarName, "%s", "Imag" );
         DBPutMultivar( dbfile, MultiVarName, nvar, varnames2, vartypes, OptList );

         sprintf( MultiVarName, "%s", "Vel" );
         DBPutMultivar( dbfile, MultiVarName, nvar, varnames3, vartypes, OptList );

         sprintf( MultiVarName, "%s", "VelT" );
         DBPutMultivar( dbfile, MultiVarName, nvar, varnames4, vartypes, OptList );

         sprintf( MultiVarName, "%s", "VelE" );
         DBPutMultivar( dbfile, MultiVarName, nvar, varnames5, vartypes, OptList );

#        else
#        error : ERROR : unsupported MODEL !!
#        endif // MODEL

         if ( OutputPot )
         {
         sprintf( MultiVarName, "%s", "Pot" );
         DBPutMultivar( dbfile, MultiVarName, nvar, varname_pot, vartypes, OptList );
         }

      } // if ( nmesh != 0 )

      free( patchnames  );
      free( meshnames   );
      free( patchtypes  );
      free( meshtypes   );
      free( varname_pot );
      free( varnames0   );
      free( varnames1   );
      free( varnames2   );
      free( varnames3   );
      free( varnames4   );
#     if ( MODEL == ELBDM )
      free( varnames5   );
#     endif
      free( vartypes    );

   } // if ( !SepAllLv )


   else // create multi-mesh object for each level
   {
      for (int lv=0; lv<NLEVEL; lv++)
      {
         char **patchnames  = NULL;
         char **meshnames   = NULL;
         char **varname_pot = NULL;
         char **varnames0   = NULL, **varnames1 = NULL, **varnames2 = NULL, **varnames3 = NULL, **varnames4 = NULL;
#        if ( MODEL == ELBDM )
         char **varnames5   = NULL;
#        endif
         int nmesh = NMesh[lv], nvar = NMesh[lv];
         int *patchtypes = NULL, *meshtypes = NULL, *vartypes = NULL;
         int line = 0;

         patchnames  = (char **) malloc(nmesh * sizeof(char *));
         meshnames   = (char **) malloc(nmesh * sizeof(char *));
         patchtypes  = (int   *) malloc(nmesh * sizeof(int   ));
         meshtypes   = (int   *) malloc(nmesh * sizeof(int   ));
         varname_pot = (char **) malloc(nvar  * sizeof(char *));
         varnames0   = (char **) malloc(nvar  * sizeof(char *));
         varnames1   = (char **) malloc(nvar  * sizeof(char *));
         varnames2   = (char **) malloc(nvar  * sizeof(char *));
         varnames3   = (char **) malloc(nvar  * sizeof(char *));
         varnames4   = (char **) malloc(nvar  * sizeof(char *));
#        if ( MODEL == ELBDM )
         varnames5   = (char **) malloc(nvar  * sizeof(char *));
#        endif
         vartypes    = (int   *) malloc(nvar  * sizeof(int   ));


         for (int Mesh=0; Mesh<nmesh; Mesh++)
         {
            char tmp[100];

            sprintf( tmp, "Level_%d.silo:Patch_%d_%d", lv, lv, Mesh );
            patchnames[line] = strdup(tmp);
            patchtypes[line] = DB_QUAD_RECT;

            sprintf( tmp, "Level_%d.silo:Mesh_%d_%d", lv, lv, Mesh );
            meshnames[line] = strdup(tmp);
            meshtypes[line] = DB_QUAD_RECT;

#           if   ( MODEL == HYDRO )
            sprintf( tmp, "Level_%d.silo:Rho_%d_%d",  lv, lv, Mesh );
            varnames0[line] = strdup(tmp);

            sprintf( tmp, "Level_%d.silo:Vel_%d_%d",  lv, lv, Mesh );
            varnames1[line] = strdup(tmp);

            sprintf( tmp, "Level_%d.silo:Egy_%d_%d",  lv, lv, Mesh );
            varnames2[line] = strdup(tmp);

            if ( OutputPres )
            {
            sprintf( tmp, "Level_%d.silo:Pre_%d_%d",  lv, lv, Mesh );
            varnames3[line] = strdup(tmp);
            }

            if ( OutputVort )
            {
            sprintf( tmp, "Level_%d.silo:Vor_%d_%d",  lv, lv, Mesh );
            varnames4[line] = strdup(tmp);
            }

#           elif ( MODEL == MHD )
#           warning : WAIT MHD !!!

#           elif ( MODEL == ELBDM )
            sprintf( tmp, "Level_%d.silo:Dens_%d_%d",  lv, lv, Mesh );
            varnames0[line] = strdup(tmp);

            sprintf( tmp, "Level_%d.silo:Real_%d_%d",  lv, lv, Mesh );
            varnames1[line] = strdup(tmp);

            sprintf( tmp, "Level_%d.silo:Imag_%d_%d",  lv, lv, Mesh );
            varnames2[line] = strdup(tmp);

            sprintf( tmp, "Level_%d.silo:Vel_%d_%d",   lv, lv, Mesh );
            varnames3[line] = strdup(tmp);

            sprintf( tmp, "Level_%d.silo:VelT_%d_%d",  lv, lv, Mesh );
            varnames4[line] = strdup(tmp);

            sprintf( tmp, "Level_%d.silo:VelE_%d_%d",  lv, lv, Mesh );
            varnames5[line] = strdup(tmp);

#           else
#           error : ERROR : unsupported MODEL !!
#           endif // MODEL

            if ( OutputPot )
            {
            sprintf( tmp, "Level_%d.silo:Pot_%d_%d",  lv, lv, Mesh );
            varname_pot[line] = strdup(tmp);
            }

            vartypes[line]  = DB_QUADVAR;

            line++;
         } // for (int Mesh=0; Mesh<nmesh; Mesh++)


         if ( nmesh != 0 )
         {
            char MultiPatchName[30];
            sprintf( MultiPatchName, "Patch_%d", lv );
            DBPutMultimesh( dbfile, MultiPatchName, nmesh, patchnames, patchtypes, OptList );

            char MultiMeshName[30];
            sprintf( MultiMeshName, "Mesh_%d", lv );
            DBPutMultimesh( dbfile, MultiMeshName, nmesh, meshnames, meshtypes, OptList );

            char MultiVarName[30];
#           if   ( MODEL == HYDRO )
            sprintf( MultiVarName, "Rho_%d", lv );
            DBPutMultivar( dbfile, MultiVarName, nvar, varnames0, vartypes, OptList );

            sprintf( MultiVarName, "Vel_%d", lv );
            DBPutMultivar( dbfile, MultiVarName, nvar, varnames1, vartypes, OptList );

            sprintf( MultiVarName, "Egy_%d", lv );
            DBPutMultivar( dbfile, MultiVarName, nvar, varnames2, vartypes, OptList );

            if ( OutputPres )
            {
            sprintf( MultiVarName, "Pre_%d", lv );
            DBPutMultivar( dbfile, MultiVarName, nvar, varnames3, vartypes, OptList );
            }

            if ( OutputVort )
            {
            sprintf( MultiVarName, "Vor_%d", lv );
            DBPutMultivar( dbfile, MultiVarName, nvar, varnames4, vartypes, OptList );
            }

#           elif ( MODEL == MHD )
#           warning : WAIT MHD !!!

#           elif ( MODEL == ELBDM )
            sprintf( MultiVarName, "Dens_%d", lv );
            DBPutMultivar( dbfile, MultiVarName, nvar, varnames0, vartypes, OptList );

            sprintf( MultiVarName, "Real_%d", lv );
            DBPutMultivar( dbfile, MultiVarName, nvar, varnames1, vartypes, OptList );

            sprintf( MultiVarName, "Imag_%d", lv );
            DBPutMultivar( dbfile, MultiVarName, nvar, varnames2, vartypes, OptList );

            sprintf( MultiVarName, "Vel_%d", lv );
            DBPutMultivar( dbfile, MultiVarName, nvar, varnames3, vartypes, OptList );

            sprintf( MultiVarName, "VelT_%d", lv );
            DBPutMultivar( dbfile, MultiVarName, nvar, varnames4, vartypes, OptList );

            sprintf( MultiVarName, "VelE_%d", lv );
            DBPutMultivar( dbfile, MultiVarName, nvar, varnames5, vartypes, OptList );

#           else
#           error : ERROR : unsupported MODEL !!
#           endif // MODEL

            if ( OutputPot )
            {
            sprintf( MultiVarName, "Pot_%d", lv );
            DBPutMultivar( dbfile, MultiVarName, nvar, varname_pot, vartypes, OptList );
            }

         } // if ( nmesh != 0 )

         free( patchnames  );
         free( meshnames   );
         free( patchtypes  );
         free( meshtypes   );
         free( varname_pot );
         free( varnames0   );
         free( varnames1   );
         free( varnames2   );
         free( varnames3   );
         free( varnames4   );
#        if ( MODEL == ELBDM )
         free( varnames5   );
#        endif
         free( vartypes    );
      }
   } // if ( !SepAllLv ) ... else ...


// write the Box into the Root file
   char *name   = "Box.silo:Box";
   int  type[1] = { DB_QUAD_RECT };
   DBPutMultimesh( dbfile, "Box", 1, &name, type, OptList );


   DBClose( dbfile );
   DBFreeOptlist( OptList );


   cout << "done" << endl;

} // FUNCTION : WriteRoot



//-------------------------------------------------------------------------------------------------------
// Function    :  ReadOption
// Description :  Read the command-line options
//-------------------------------------------------------------------------------------------------------
void ReadOption( int argc, char **argv )
{

   cout << "ReadOption ... " << endl;


   double Temp_Start[3] = { NULL_VALUE, NULL_VALUE, NULL_VALUE };
   double Temp_Size [3] = { NULL_VALUE, NULL_VALUE, NULL_VALUE };
   int c;

   while( (c = getopt(argc, argv, "hdpwsi:x:y:z:X:Y:Z:m:")) != -1 )
   switch(c)
   {
      case 'i': FileName_In   = optarg;
                break;
      case 'x': Temp_Start[0] = atof(optarg);
                break;
      case 'y': Temp_Start[1] = atof(optarg);
                break;
      case 'z': Temp_Start[2] = atof(optarg);
                break;
      case 'X': Temp_Size[0]  = atof(optarg);
                break;
      case 'Y': Temp_Size[1]  = atof(optarg);
                break;
      case 'Z': Temp_Size[2]  = atof(optarg);
                break;
      case 'm': TargetMode    = atoi(optarg);
                break;
      case 'd': SepAllLv      = true;
                break;
#     if ( MODEL == HYDRO )
      case 'p': OutputPres    = true;
                break;
      case 'w': OutputVort    = true;
                break;
#     endif
      case 's': InputScale    = true;
                break;
      case 'h':
      case '?': cerr << "usage: " << argv[0]
                     << " [-h (for help)] [-i Input filename] [-d (separate data at different levels) [off]]"
                     << endl << "                         "
                     << " [-(x,y,z) starting coordinate [0]] [-(X/Y/Z) size [BoxSize]]"
                     << endl << "                         "
                     << " [-m (0,1,2,3) --> (xy slice, yz slice, xz slice, 3D)]"
                     << endl << "                         "
#                    if ( MODEL == HYDRO )
                     << " [-p (output pressure) [off]] [-w (output vorticity vector) [off]]"
                     << endl << "                         "
#                    endif
                     << " [-s [input cell scales instead of physical coordinates to specify the range] [off]]"
                     << endl;
                exit( 1 );
   }


// set target range
   if ( InputScale )
   {
      for (int d=0; d<3; d++)
      {
         CornerScale[d] = (int)Temp_Start[d];
         SizeScale  [d] = (int)Temp_Size [d];
      }
   }

   else
   {
      for (int d=0; d<3; d++)
      {
         PhyCorner[d] = Temp_Start[d];
         PhySize  [d] = Temp_Size [d];
      }
   }


// check
   if ( TargetMode == NULL_VALUE  ||  TargetMode < 0  ||  TargetMode > 3 )
   {
      fprintf( stderr, "Error : please provide the correct target mode (-m 0/1/2/3) !!\n" );
      exit( 1 );
   }

   FILE *File = fopen( FileName_In, "rb" );
   if ( File == NULL )
   {
      fprintf( stderr, "Error : the input filename \"%s\" does not exist (-i Filename) !!\n", FileName_In );
      exit( 1 );
   }
   fclose( File );


   cout << "ReadOption ... done" << endl;


// record the command-line options
   fprintf( stdout, "---------------------------------------------------------------------------------------------------\n" );
   fprintf( stdout, "Command-line arguments :\n" );
   for (int v=0; v<argc; v++)    fprintf( stdout, " %s", argv[v] );
   fprintf( stdout, "\n" );
   fprintf( stdout, "---------------------------------------------------------------------------------------------------\n" );

} // FUNCTION : ReadOption



//-------------------------------------------------------------------------------------------------------
// Function    :  SetDefaultParameter
// Description :  Set parameters to the default values
//-------------------------------------------------------------------------------------------------------
void SetDefaultParameter()
{

   if ( InputScale )
   {
      for (int d=0; d<3; d++)
      {
         if ( CornerScale[d] == (int)NULL_VALUE )     CornerScale[d] = 0;

         if ( TargetMode == (d+1)%3 )                 SizeScale  [d] = 1;
         else if ( SizeScale[d] == (int)NULL_VALUE )  SizeScale  [d] = amr.BoxScale[d] - CornerScale[d];
      }
   }

   else
   {
      for (int d=0; d<3; d++)
      {
         if ( PhyCorner[d] == NULL_VALUE )      PhyCorner[d] = 0.0;

         if ( TargetMode == (d+1)%3 )           PhySize  [d] = 0.0;
         else if ( PhySize[d] == NULL_VALUE )   PhySize  [d] = amr.BoxSize[d] - PhyCorner[d];
      }
   }

} // FUNCTION : SetDefaultParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  CheckParameter
// Description :  verify parameters
//-------------------------------------------------------------------------------------------------------
void CheckParameter()
{

   cout << "   CheckParameter ... " << endl;


#  if ( MODEL != HYDRO  &&  MODEL != ELBDM )
#  error : ERROR : unsupported MODEL !!
#  endif

   for (int d=0; d<3; d++)
   {
      if ( PhyCorner[d] < 0.0  ||  PhyCorner[d] >= amr.BoxSize[d] )
      {
         fprintf( stderr, "Error : starting coordinate [%d] = %lf lies outside the simulation box !!\n",
                  d, PhyCorner[d] );
         exit( 1 );
      }

      if (  TargetMode != (d+1)%3  &&  ( PhySize[d] <= 0.0 || PhyCorner[d]+PhySize[d] > amr.BoxSize[d] )  )
      {
         fprintf( stderr, "Error : incorrect target size [%d] = %lf !!\n", d, PhySize[d] );
         exit( 1 );
      }
   }

   for (int d=0; d<3; d++)
   {
      if (  CornerScale[d] >= amr.BoxScale[d]  ||  CornerScale[d] < 0  )
      {
         fprintf( stderr, "Error : the starting grid index [%d] = %d lies outside the simulation box !!\n",
                  d, CornerScale[d] );
         exit( 1 );
      }

      if ( CornerScale[d]+SizeScale[d] > amr.BoxScale[d] )
      {
         fprintf( stderr, "Error : the ending grid index [%d] = %d exceeds the simulation box !!\n",
                  d, CornerScale[d]+SizeScale[d] );
         exit( -1 );
      }
   }


   cout << "   CheckParameter ... done" << endl;

} // FUNCTION : CheckParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  TruncateBox
// Description :  Truncate the targeted slice/box to fit the base-level patches
//-------------------------------------------------------------------------------------------------------
void TruncateBox()
{

   cout << "   TruncateBox ... " << endl;


   if ( InputScale ) // cell scale --> physical coordinates
   {
      for (int d=0; d<3; d++)    PhyCorner[d] = (double)CornerScale[d]/(double)amr.BoxScale[d]*amr.BoxSize[d];

      switch ( TargetMode )
      {
         case 0:
            PhySize[0] = (double)SizeScale[0]/(double)amr.BoxScale[0]*amr.BoxSize[0];
            PhySize[1] = (double)SizeScale[1]/(double)amr.BoxScale[1]*amr.BoxSize[1];
            PhySize[2] = 0.0;
            break;

         case 1:
            PhySize[0] = 0.0;
            PhySize[1] = (double)SizeScale[1]/(double)amr.BoxScale[1]*amr.BoxSize[1];
            PhySize[2] = (double)SizeScale[2]/(double)amr.BoxScale[2]*amr.BoxSize[2];
            break;

         case 2:
            PhySize[0] = (double)SizeScale[0]/(double)amr.BoxScale[0]*amr.BoxSize[0];
            PhySize[1] = 0.0;
            PhySize[2] = (double)SizeScale[2]/(double)amr.BoxScale[2]*amr.BoxSize[2];
            break;

         case 3:
            PhySize[0] = (double)SizeScale[0]/(double)amr.BoxScale[0]*amr.BoxSize[0];
            PhySize[1] = (double)SizeScale[1]/(double)amr.BoxScale[1]*amr.BoxSize[1];
            PhySize[2] = (double)SizeScale[2]/(double)amr.BoxScale[2]*amr.BoxSize[2];
            break;

         default:
            fprintf( stderr, "ERROR : incorrect variable (%s = %d) !!\n", "TargetMode", TargetMode );
            exit( 1 );

      } // switch ( TargetMode )
   } // if ( InputScale )

   else // physical coordinates --> cell scale
   {
      for (int d=0; d<3; d++)    CornerScale[d] = int( PhyCorner[d]/amr.BoxSize[d]*amr.BoxScale[d] );

      switch ( TargetMode )
      {
         case 0:
            SizeScale[0] = int( PhySize[0]/amr.BoxSize[0]*amr.BoxScale[0] );
            SizeScale[1] = int( PhySize[1]/amr.BoxSize[1]*amr.BoxScale[1] );
            SizeScale[2] = 1;
            break;

         case 1:
            SizeScale[0] = 1;
            SizeScale[1] = int( PhySize[1]/amr.BoxSize[1]*amr.BoxScale[1] );
            SizeScale[2] = int( PhySize[2]/amr.BoxSize[2]*amr.BoxScale[2] );
            break;

         case 2:
            SizeScale[0] = int( PhySize[0]/amr.BoxSize[0]*amr.BoxScale[0] );
            SizeScale[1] = 1;
            SizeScale[2] = int( PhySize[2]/amr.BoxSize[2]*amr.BoxScale[2] );
            break;

         case 3:
            SizeScale[0] = int( PhySize[0]/amr.BoxSize[0]*amr.BoxScale[0] );
            SizeScale[1] = int( PhySize[1]/amr.BoxSize[1]*amr.BoxScale[1] );
            SizeScale[2] = int( PhySize[2]/amr.BoxSize[2]*amr.BoxScale[2] );
            break;

         default:
            fprintf( stderr, "ERROR : incorrect variable (%s = %d) !!\n", "TargetMode", TargetMode );
            exit( 1 );

      } // switch ( TargetMode )
   } // if ( InputScale ) ... else ...


// truncate the targeted region to fit the root-level grids
   const int PS0 = amr.scale[0]*PATCH_SIZE;

   switch ( TargetMode )
   {
      case 0:
         CornerScale[0] = PS0 * ( CornerScale[0]/PS0 );
         CornerScale[1] = PS0 * ( CornerScale[1]/PS0 );
         SizeScale  [0] = PS0 * ( (SizeScale[0]-1)/PS0 + 1 );
         SizeScale  [1] = PS0 * ( (SizeScale[1]-1)/PS0 + 1 );
         break;

      case 1:
         CornerScale[1] = PS0 * ( CornerScale[1]/PS0 );
         CornerScale[2] = PS0 * ( CornerScale[2]/PS0 );
         SizeScale  [1] = PS0 * ( (SizeScale[1]-1)/PS0 + 1 );
         SizeScale  [2] = PS0 * ( (SizeScale[2]-1)/PS0 + 1 );
         break;

      case 2:
         CornerScale[0] = PS0 * ( CornerScale[0]/PS0 );
         CornerScale[2] = PS0 * ( CornerScale[2]/PS0 );
         SizeScale  [0] = PS0 * ( (SizeScale[0]-1)/PS0 + 1 );
         SizeScale  [2] = PS0 * ( (SizeScale[2]-1)/PS0 + 1 );
         break;

      case 3:
         CornerScale[0] = PS0 * ( CornerScale[0]/PS0 );
         CornerScale[1] = PS0 * ( CornerScale[1]/PS0 );
         CornerScale[2] = PS0 * ( CornerScale[2]/PS0 );
         SizeScale  [0] = PS0 * ( (SizeScale[0]-1)/PS0 + 1 );
         SizeScale  [1] = PS0 * ( (SizeScale[1]-1)/PS0 + 1 );
         SizeScale  [2] = PS0 * ( (SizeScale[2]-1)/PS0 + 1 );
         break;
   }


   for (int d=0; d<3; d++)
   {
      TPhySize  [d] = (double)SizeScale  [d]/amr.BoxScale[d]*amr.BoxSize[d];
      TPhyCorner[d] = (double)CornerScale[d]/amr.BoxScale[d]*amr.BoxSize[d];
   }

   printf( "      Truncation results:\n" );
   printf( "      mode           = %8d\n",     TargetMode     );
   printf( "      corner x       = %14.5lf\n", TPhyCorner [0] );
   printf( "      corner y       = %14.5lf\n", TPhyCorner [1] );
   printf( "      corner z       = %14.5lf\n", TPhyCorner [2] );
   printf( "      size x         = %14.5lf\n", TPhySize   [0] );
   printf( "      size y         = %14.5lf\n", TPhySize   [1] );
   printf( "      size z         = %14.5lf\n", TPhySize   [2] );
   printf( "      corner scale x = %8d\n",     CornerScale[0] );
   printf( "      corner scale y = %8d\n",     CornerScale[1] );
   printf( "      corner scale z = %8d\n",     CornerScale[2] );
   printf( "      size scale x   = %8d\n",     SizeScale  [0] );
   printf( "      size scale y   = %8d\n",     SizeScale  [1] );
   printf( "      size scale z   = %8d\n",     SizeScale  [2] );


   cout << "   TruncateBox ... done" << endl;

} // FUNCTION : TruncateBox



//-------------------------------------------------------------------------------------------------------
// Function    :  main
// Description :
//-------------------------------------------------------------------------------------------------------
int main( int argc, char ** argv )
{

    ReadOption( argc, argv );

    LoadData();

    CreateSilo();

    cout << "Program terminated successfully" << endl;
    return 0;

} // FUNCTION : main
