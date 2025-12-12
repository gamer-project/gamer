#include "GAMER.h"

#ifdef FEEDBACK

//-------------------------------------------------------------------------------------------------------
// Function    :  FB_Accretion
// Description :  Accretion for particle, Ref: Federrath et al. 2010.
// Note        :  1. Input and output fluid and particle data are stored in Fluid[] and ParAttFlt[], respectively
//                   --> This function is responsible for updating gas and particles within
//                       ** FB_GHOST_SIZE <= cell indices i,j,k < FB_GHOST_SIZE+PS2 **
//                   --> Updating gas and particles outside this range is fine but will have no effect at all
//                2. Must use ParSortID[] to access ParAttFlt[]
//                   --> ParAttFlt[PAR_MASS/PAR_POSX/etc][ ParSortID[...] ]
//                3. Particles may be outside the target region
//                4. To ensure the consistency of random numbers, one must call the random number generator for
//                   ALL particles, including those too far away to affect the target region
//                5. No need to worry about the periodic boundary condition here
//                   --> Particle positions have been remapped in FB_AdvanceDt()
//                6. CoarseFine[] records the coarse-fine boundaries along the 26 sibling directions, defined as
//                            24  13  25
//                            15  05  17     z+1 plane
//                            22  12  23
//
//                            08  03  09
//                            00  XX  01     z   plane
//                            06  02  07
//                   y
//                   ^        20  11  21
//                   |        14  04  16     z-1 plane
//                   --->x    18  10  19
//                7. Invoked by FB_AdvanceDt()
//                8. Must NOT change particle positions
//                9. Since Fluid[] stores both the input and output data, the order of particles may affect the
//                   final output results
//                   --> For example, particle 2 may use the data updated by particle 1 as the input data
//                   --> Actually, even if we separate Fluid[] to input and output arrays, the final output results
//                       may still depend on the order of particles for non-local feedback since different particles
//                       may update the same cell
//                10. In general, it is recommended to have the maximum feedback radius no larger than half of the patch size
//                    (i.e., PATCH_SIZE/2=4 cells for PATCH_SIZE=8)
//                    --> Increase PATCH_SIZE if necessary
//                11. Linked to FB_User_Ptr in FB_Init_Plummer()
//                12. Must have FB_GHOST_SIZE>=1 for the mass accretion feedback
//
// Parameter   :  lv           : Target refinement level
//                GasDensThres : Minimum gas density for creating sink particles                (--> "SF_CREATE_SINK_MIN_GAS_DENS"  )
//                AccCellNum   : Accretion radius in cells at highest refinement level          (--> "SF_CREATE_SINK_ACC_RADIUS"    )
//                NPar         : Number of particles
//                ParSortID    : Sorted particle IDs
//                ParAttFlt    : Particle floating-point attribute arrays
//                ParAttInt    : Particle integer        attribute arrays
//                Fluid        : Array to store the input/output fluid data
//                               --> Array size is fixed to (FB_NXT)^3=(PS2+2*FB_GHOST_SIZE)^3
//                EdgeL        : Left edge of Fluid[]
//                               --> Right edge is given by EdgeL[]+FB_NXT*dh
//                dh           : Cell size of Fluid[]
//                CoarseFine   : Coarse-fine boundaries along the 26 sibling directions
//
// Return      :  Fluid, ParAttFlt
//-------------------------------------------------------------------------------------------------------
int FB_Accretion( const int lv, const real GasDensThres, const real AccCellNum, const int NPar, const long *ParSortID, 
                  real_par *ParAttFlt[PAR_NATT_FLT_TOTAL], long_par *ParAttInt[PAR_NATT_INT_TOTAL], 
                  real (*Fluid)[FB_NXT][FB_NXT][FB_NXT], const double EdgeL[], const double dh, bool CoarseFine[] )
{

// check
#  ifdef GAMER_DEBUG
   if ( Fluid == NULL )    Aux_Error( ERROR_INFO, "Fluid == NULL !!\n" );
   if ( NPar > 0 )
   {
      if ( ParSortID == NULL )   Aux_Error( ERROR_INFO, "ParSortID == NULL for NPar = %d !!\n", NPar );
      if ( ParAttFlt == NULL )   Aux_Error( ERROR_INFO, "ParAttFlt == NULL for NPar = %d !!\n", NPar );
   }
#  endif // #ifdef GAMER_DEBUG

   if ( FB_GHOST_SIZE < 2*AccCellNum )
      Aux_Error( ERROR_INFO, "FB_GHOST_SIZE should be larger than twice of SF_CREATE_SINK_ACC_RADIUS !!" );

   real GasDens, Eg, Ekin, GasCell2ParDisti, GasRelVel[3]; 
   real Eg2, GasCell2ParDist2;
   real ControlPos[3];
   real Corner_Array[3]; // the corner of the ghost zone

   const int    NGhost         = PS1 / 2; // the number of ghost cell at each side
   const int    MaxRemovalGas  = CUBE( PS1 );
   const double AccRadius      = AccCellNum*dh;
   const double _dh            = 1.0 / dh;
   const double dv             = CUBE( dh );
   const double _dv            = 1.0 / dv;
   const double epsilon        = 1e-5*dh;

   long   (*RemovalIdx)[3]         = new long [MaxRemovalGas][3];

// prepare the corner array
   for (int d=0; d<3; d++)    Corner_Array[d] = EdgeL[d] + 0.5*dh ;

   bool CheckCF = false;
#  if ( FB_GHOST_SIZE > 0 )
   for (int s=0; s<26; s++) {
      if ( CoarseFine[s] ) {
         CheckCF = true;
         break;
      }
   }
#  endif

   for (int t=0; t<NPar; t++)
   {
      const int    p      = ParSortID[t];
      const double xyz[3] = { ParAttFlt[PAR_POSX][p], ParAttFlt[PAR_POSY][p], ParAttFlt[PAR_POSZ][p] }; // particle position

      int idx[3]; // cell idx in FB_NXT^3
      for (int d=0; d<3; d++)    idx[d] = (int)floor( ( xyz[d] - EdgeL[d] )*_dh );

//    skip particles too close to coarse-fine boundaries
      bool Skip = false;

      if ( CheckCF )
      {
//       only need to check the 27 extreme cases
         int ijk[3];
         for (int dk=-1; dk<=1; dk++)  {  if ( Skip )  break;  ijk[2] = idx[2] + dk*NGhost;
         for (int dj=-1; dj<=1; dj++)  {  if ( Skip )  break;  ijk[1] = idx[1] + dj*NGhost;
         for (int di=-1; di<=1; di++)  {  if ( Skip )  break;  ijk[0] = idx[0] + di*NGhost;

            const int CellPatchRelPos = FB_Aux_CellPatchRelPos( ijk );
            if ( CellPatchRelPos != -1  &&  CoarseFine[CellPatchRelPos] )    Skip = true;   // cell is in a coarse patch

         }}}
      }

      if ( Skip )    continue;

      if ( idx[0] < FB_GHOST_SIZE-NGhost  ||  idx[0] >= FB_GHOST_SIZE+PS2+NGhost  ||
           idx[1] < FB_GHOST_SIZE-NGhost  ||  idx[1] >= FB_GHOST_SIZE+PS2+NGhost  ||
           idx[2] < FB_GHOST_SIZE-NGhost  ||  idx[2] >= FB_GHOST_SIZE+PS2+NGhost   ) // we want completed control volume
         continue;

      int NRemove = 0;

//    Check the control volume
      for (int vki=idx[2]-NGhost; vki<=idx[2]+NGhost; vki++)
      for (int vji=idx[1]-NGhost; vji<=idx[1]+NGhost; vji++)
      for (int vii=idx[0]-NGhost; vii<=idx[0]+NGhost; vii++) // loop the nearby cells, to find the cells inside the control volumne (v)
      {
//       Inside the accretion radius:
//       The test cell must be inside the accretion radius.
//       ===========================================================================================================
         ControlPos[0] = Corner_Array[0] + vii*dh;
         ControlPos[1] = Corner_Array[1] + vji*dh;
         ControlPos[2] = Corner_Array[2] + vki*dh;

         GasCell2ParDisti = SQRT(SQR(ControlPos[0] - xyz[0])+SQR(ControlPos[1] - xyz[1])+SQR(ControlPos[2] - xyz[2])); // distance to the sink
         if ( GasCell2ParDisti > AccRadius )                 continue; // check whether it is inside the control volume

//       Density threshold:
//       Check whether the gas density is larger than the threshold
//       ===========================================================================================================
         GasDens = Fluid[DENS][vki][vji][vii];
         if ( GasDens <= GasDensThres )                continue;

//       Negative radial velocity:
//       The gas must be moving toward the particle.
//       ===========================================================================================================
         real RadialVel;

         GasRelVel[0] = Fluid[MOMX][vki][vji][vii]/GasDens - ParAttFlt[PAR_VELX][p];
         GasRelVel[1] = Fluid[MOMY][vki][vji][vii]/GasDens - ParAttFlt[PAR_VELY][p];
         GasRelVel[2] = Fluid[MOMZ][vki][vji][vii]/GasDens - ParAttFlt[PAR_VELZ][p];

         RadialVel = (GasRelVel[0]*(ControlPos[0] - xyz[0]) + GasRelVel[1]*(ControlPos[1] - xyz[1]) + GasRelVel[2]*(ControlPos[2] - xyz[2]))/GasCell2ParDisti;

         if ( SIGN(RadialVel) != -1.0 )    continue;

//       Bound state check:
//       The gas should be bound to the particle.
//       ===========================================================================================================
         real SelfPhi = (real)0.0; // self-potential

         SelfPhi += -NEWTON_G*ParAttFlt[PAR_MASS][p]/(GasCell2ParDisti+epsilon); // potential from the sink

         Eg   = GasDens*dv*SelfPhi; // gravitational potential energy
         Ekin = 0.5*GasDens*dv*( SQR(GasRelVel[0]) + SQR(GasRelVel[1]) + SQR(GasRelVel[2])); // kinetic energy

         if ( ( Eg + Ekin ) >= 0 )      continue;

//       Overlapped accretion radius check:
//       Check whether the cell is overlapped with another accretion radius.
//       ===========================================================================================================
         bool NotMinEg = false;
         for (int tt=0; tt<NPar; tt++) // find the nearby sink
         {
            const int    pp      = ParSortID[tt];
            if ( pp == p )              continue;

            const double xxyyzz[3] = { ParAttFlt[PAR_POSX][pp], ParAttFlt[PAR_POSY][pp], ParAttFlt[PAR_POSZ][pp] }; // particle position
            
            int idxx[3]; // cell idx in FB_NXT^3
            for (int d=0; d<3; d++)    idxx[d] = (int)floor( ( xxyyzz[d] - EdgeL[d] )*_dh );

            GasCell2ParDist2 = SQRT(SQR(ControlPos[0] - xxyyzz[0])+SQR(ControlPos[1] - xxyyzz[1])+SQR(ControlPos[2] - xxyyzz[2])); // distance to the sink
            if ( GasCell2ParDist2 > AccRadius )       continue;

            if ( idxx[0] < FB_GHOST_SIZE-NGhost  ||  idxx[0] >= FB_GHOST_SIZE+PS2+NGhost  ||
                 idxx[1] < FB_GHOST_SIZE-NGhost  ||  idxx[1] >= FB_GHOST_SIZE+PS2+NGhost  ||
                 idxx[2] < FB_GHOST_SIZE-NGhost  ||  idxx[2] >= FB_GHOST_SIZE+PS2+NGhost   ) // we want completed control volume
               continue;

            real SelfPhi2 = (real)0.0; // self-potential

            SelfPhi2 += -NEWTON_G*ParAttFlt[PAR_MASS][pp]/(GasCell2ParDist2+epsilon); // potential from the sink

            Eg2   = GasDens*dv*SelfPhi2;
            if ( Eg2 <= Eg ) // whether the gas cell is more bound to another particle
            {
               NotMinEg = true;
               break;
            }
         } // for (int tt=0; tt<NPar; tt++)

         if ( NotMinEg )                              continue; 

//       Record the index
//       ===========================================================================================================
         RemovalIdx[NRemove][0] = vii;
         RemovalIdx[NRemove][1] = vji;
         RemovalIdx[NRemove][2] = vki;

         NRemove ++;
      } // vii, vji, vki

      int i, j, k;
      real DeltaM;
      real DeltaMTot   = (real)0.0;
      real DeltaMom[3] = { (real)0.0, (real)0.0, (real)0.0 };

      for (int N=0; N<NRemove; N++)
      {
         i = RemovalIdx[N][0];
         j = RemovalIdx[N][1];
         k = RemovalIdx[N][2];

         GasDens = Fluid[DENS][k][j][i];

         DeltaM     = (GasDens - GasDensThres)*dv; // the mass to be accreted
         DeltaMTot  += DeltaM;

         DeltaMom[0] += DeltaM*Fluid[MOMX][k][j][i]/GasDens; // the momentum of DeltaM
         DeltaMom[1] += DeltaM*Fluid[MOMY][k][j][i]/GasDens;
         DeltaMom[2] += DeltaM*Fluid[MOMZ][k][j][i]/GasDens;

//       Update the cells
//       ===========================================================================================================
         for (int v=0; v<NCOMP_TOTAL; v++)
         Fluid[v][k][j][i] *= GasDensThres/GasDens;
      } // for (int N=0; N<NRemove; N++)

//    Update particle mass and velocity
//    ===========================================================================================================
      ParAttFlt[PAR_VELX][p] =  (DeltaMom[0] + ParAttFlt[PAR_MASS][p]*ParAttFlt[PAR_VELX][p])/(DeltaMTot + ParAttFlt[PAR_MASS][p]);  // center-of-mass velocity of the sink after accretion
      ParAttFlt[PAR_VELY][p] =  (DeltaMom[1] + ParAttFlt[PAR_MASS][p]*ParAttFlt[PAR_VELY][p])/(DeltaMTot + ParAttFlt[PAR_MASS][p]);
      ParAttFlt[PAR_VELZ][p] =  (DeltaMom[2] + ParAttFlt[PAR_MASS][p]*ParAttFlt[PAR_VELZ][p])/(DeltaMTot + ParAttFlt[PAR_MASS][p]);
      ParAttFlt[PAR_MASS][p] +=  DeltaMTot;
   } // for (int t=0; t<NPar; t++)

   delete [] RemovalIdx;

   return GAMER_SUCCESS;

} // FUNCTION : FB_SinkAccretion



//-------------------------------------------------------------------------------------------------------
// Function    :  FB_End_Accretion
// Description :  Free the resources used by the user-specified feedback
//
// Note        :  1. Invoked by FB_End()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void FB_End_Accretion()
{

} // FUNCTION : FB_End_Accretion



//-------------------------------------------------------------------------------------------------------
// Function    :  FB_Init_SinkAccretion
// Description :  Initialize the user-specified feedback
//
// Note        :  1. Invoked by FB_Init()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void FB_Init_Accretion()
{

} // FUNCTION : FB_Init_Accretion



#endif // #ifdef FEEDBACK