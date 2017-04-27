#include "GAMER.h"

static void Flu_ResetByUser_Func( real fluid[], const double x, const double y, const double z, const double Time );

extern int      Jet_NJet;
extern double  *Jet_Radius;
extern double  *Jet_HalfHeight;
extern double  *Jet_SrcVel;
extern double  *Jet_SrcDens;
extern double  *Jet_SrcTemp;
extern double (*Jet_Vec)[3];
extern double  *Jet_SrcEint;
extern double (*Jet_Cen)[3];
extern double  *Jet_WaveK;




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_ResetByUser_Func
// Description :  Function to reset the fluid field
//
// Note        :  1. Invoked by "Flu_ResetByUser"
//                2. Input "fluid" array stores the original values
//                3. This function will reset the values inside the jet source region to the initial values
//                4. The source regions of different jets must NOT overlap
//
// Parameter   :  fluid : Fluid array storing both the input (origial) and reset values
//                        --> Including both active and passive variables
//                x/y/z : Target physical coordinates
//                Time  : Target physical time
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void Flu_ResetByUser_Func( real fluid[], const double x, const double y, const double z, const double Time )
{

   const double r[3] = { x, y, z };

   double Jet_dr, Jet_dh, S, Area;
   double Dis_c2m, Dis_c2v, Dis_v2m, Vec_c2m[3], Vec_v2m[3];
   real   MomSin;


// vectors for calculating the distance between cells and the jet sources
   double TempVec[Jet_NJet][3];

   for (int n=0; n<Jet_NJet; n++)
   for (int d=0; d<3; d++)
      TempVec[n][d] = Jet_Cen[n][d] + Jet_Vec[n][d];


// velocity components along different directions
   double _Jet_VecAbs, Jet_SrcVel_xyz[Jet_NJet][3];

   for (int n=0; n<Jet_NJet; n++)
   {
      _Jet_VecAbs = 1.0/sqrt( SQR(Jet_Vec[n][0]) + SQR(Jet_Vec[n][1]) + SQR(Jet_Vec[n][2]) );

      for (int d=0; d<3; d++)    Jet_SrcVel_xyz[n][d] = Jet_SrcVel[n] * Jet_Vec[n][d] * _Jet_VecAbs;
   }


// loop over all cells to add the jet source
   for (int n=0; n<Jet_NJet; n++)
   {
//    exclude cells far away from the jet source
      if ( fabs(r[2]-Jet_Cen[n][2]) > Jet_HalfHeight[n] )   continue;

//    distance: jet center to mesh
      for (int d=0; d<3; d++)    Vec_c2m[d] = r[d] - Jet_Cen[n][d];
      Dis_c2m = sqrt( SQR(Vec_c2m[0]) + SQR(Vec_c2m[1]) + SQR(Vec_c2m[2]) );

//    distance: temporary vector to mesh
      for (int d=0; d<3; d++)    Vec_v2m[d] = r[d] - TempVec[n][d];
      Dis_v2m = sqrt( SQR(Vec_v2m[0]) + SQR(Vec_v2m[1]) + SQR(Vec_v2m[2]) );

//    distance: jet center to temporary vector
      Dis_c2v = sqrt( SQR(Jet_Vec[n][0]) + SQR(Jet_Vec[n][1]) + SQR(Jet_Vec[n][2]) );

//    check whether or not the target cell is within the jet source
      S      = 0.5*( Dis_c2m + Dis_v2m + Dis_c2v );
      Area   = sqrt( S*(S-Dis_c2m)*(S-Dis_v2m)*(S-Dis_c2v) );
      Jet_dr = 2.0*Area/Dis_c2v;
      Jet_dh = sqrt( Dis_c2m*Dis_c2m - Jet_dr*Jet_dr );

      if ( Jet_dh <= Jet_HalfHeight[n]  &&  Jet_dr <= Jet_Radius[n] )
      {
         fluid[DENS] = Jet_SrcDens[n];

//       use a sine function to make the velocity smooth within the jet from +Jet_SrcVel to -Jet_SrcVel
         MomSin      = fluid[DENS]*sin( Jet_WaveK[n]*Jet_dh );
         MomSin     *= SIGN( Vec_c2m[0]*Jet_Vec[n][0] + Vec_c2m[1]*Jet_Vec[n][1] + Vec_c2m[2]*Jet_Vec[n][2] );
         fluid[MOMX] = MomSin*Jet_SrcVel_xyz[n][0];
         fluid[MOMY] = MomSin*Jet_SrcVel_xyz[n][1];
         fluid[MOMZ] = MomSin*Jet_SrcVel_xyz[n][2];

         fluid[ENGY] = Jet_SrcEint[n] + 0.5*( SQR(fluid[MOMX])+SQR(fluid[MOMY])+SQR(fluid[MOMZ]) ) / fluid[DENS];
      } // if (  Jet_dh <= Jet_HalfHeight[n]  &&  Jet_dr <= Jet_Radius[n] )
   } // for (int n=0; n<Jet_NJet; n++)

} // FUNCTION : Flu_ResetByUser_Func



//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_ResetByUser
// Description :  Reset the fluid array
//
// Note        :  1. Work for the option "OPT__RESET_FLUID == true"
//                2. Invokded by both "Init_StartOver" and "EvolveLevel"
//                3. Will not be applied to reset the input uniform array
//                   --> Init_UM does NOT call this function
//                4. Currently does not work with "OPT__OVERLAP_MPI"
//
// Parameter   :  lv    :  Target refinement level
//                FluSg :  Target fluid sandglass
//                TTime :  Target physical time
//-------------------------------------------------------------------------------------------------------
void Flu_ResetByUser( const int lv, const int FluSg, const double TTime )
{

   const double dh = amr->dh[lv];

   real   fluid[NCOMP_TOTAL];
   double x, y, z, x0, y0, z0;


#  pragma omp parallel for private( fluid, x, y, z, x0, y0, z0 ) schedule( runtime )
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      x0 = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh;
      y0 = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh;
      z0 = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh;

      for (int k=0; k<PS1; k++)  {  z = z0 + k*dh;
      for (int j=0; j<PS1; j++)  {  y = y0 + j*dh;
      for (int i=0; i<PS1; i++)  {  x = x0 + i*dh;

         for (int v=0; v<NCOMP_TOTAL; v++)   fluid[v] = amr->patch[FluSg][lv][PID]->fluid[v][k][j][i];

         Flu_ResetByUser_Func( fluid, x, y, z, TTime );

         for (int v=0; v<NCOMP_TOTAL; v++)   amr->patch[FluSg][lv][PID]->fluid[v][k][j][i] = fluid[v];

      }}}
   }

} // FUNCTION : Flu_ResetByUser
