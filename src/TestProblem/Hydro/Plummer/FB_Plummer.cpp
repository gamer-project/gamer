#include "GAMER.h"

#ifdef FEEDBACK

extern bool   Plummer_FB_Exp;
extern double Plummer_FB_ExpEMin;
extern double Plummer_FB_ExpEMax;
extern double Plummer_FB_ExpMMin;
extern double Plummer_FB_ExpMMax;
extern bool   Plummer_FB_Acc;
extern double Plummer_FB_Like;


// function pointers to be set by FB_Init_Plummer()
extern void (*FB_User_Ptr)( const int lv, const double TimeNew, const double TimeOld, const double dt,
                            const int NPar, const int *ParSortID, real *ParAtt[PAR_NATT_TOTAL],
                            real (*Fluid)[FB_NXT][FB_NXT][FB_NXT], const double EdgeL[], const double dh, bool CoarseFine[],
                            const int TID, RandomNumber_t *RNG );
extern void (*FB_End_User_Ptr)();




//-------------------------------------------------------------------------------------------------------
// Function    :  FB_Plummer
// Description :  Example of explosion and mass accretion feedbacks
//                --> Just an example and is by no means to be physically correct
//
// Note        :  1. Input and output fluid and particle data are stored in Fluid[] and ParAtt[], respectively
//                   --> This function is responsible for updating gas and particles within
//                       ** FB_GHOST_SIZE <= cell indices i,j,k < FB_GHOST_SIZE+PS2 **
//                   --> Updating gas and particles outside this range is fine but will have no effect at all
//                2. Must use ParSortID[] to access ParAtt[]
//                   --> ParAtt[PAR_MASS/PAR_POSX/etc][ ParSortID[...] ]
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
//                10. Linked to FB_User_Ptr in FB_Init_Plummer()
//
// Parameter   :  lv         : Target refinement level
//                TimeNew    : Target physical time to reach
//                TimeOld    : Physical time before update
//                             --> This function updates physical time from TimeOld to TimeNew
//                dt         : Time interval to advance solution
//                NPar       : Number of particles
//                ParSortID  : Sorted particle IDs
//                ParAtt     : Particle attribute arrays
//                Fluid      : Array to store the input/output fluid data
//                             --> Array size is fixed to (FB_NXT)^3=(PS2+2*FB_GHOST_SIZE)^3
//                EdgeL      : Left edge of Fluid[]
//                             --> Right edge is given by EdgeL[]+FB_NXT*dh
//                dh         : Cell size of Fluid[]
//                CoarseFine : Coarse-fine boundaries along the 26 sibling directions
//                TID        : Thread ID
//                RNG        : Random number generator
//                             --> Random number can be obtained by "RNG->GetValue( TID, Min, Max )",
//                                 where Min/Max specify the range of random numbers
//
// Return      :  Fluid, ParAtt
//-------------------------------------------------------------------------------------------------------
void FB_Plummer( const int lv, const double TimeNew, const double TimeOld, const double dt,
                 const int NPar, const int *ParSortID, real *ParAtt[PAR_NATT_TOTAL],
                 real (*Fluid)[FB_NXT][FB_NXT][FB_NXT], const double EdgeL[], const double dh, bool CoarseFine[],
                 const int TID, RandomNumber_t *RNG )
{

// check
#  ifdef GAMER_DEBUG
   if ( Fluid == NULL )    Aux_Error( ERROR_INFO, "Fluid == NULL !!\n" );
   if ( NPar > 0 )
   {
      if ( ParSortID == NULL )   Aux_Error( ERROR_INFO, "ParSortID == NULL for NPar = %d !!\n", NPar );
      if ( ParAtt == NULL )      Aux_Error( ERROR_INFO, "ParAtt == NULL for NPar = %d !!\n", NPar );
   }
#  endif // #ifdef GAMER_DEBUG


   const double _dh  = 1.0 / dh;

   for (int t=0; t<NPar; t++)
   {
      if ( RNG->GetValue(TID,0.0,1.0) < Plummer_FB_Like )
      {
         const int    p       = ParSortID[t];
         const double xyz[3]  = { ParAtt[PAR_POSX][p], ParAtt[PAR_POSY][p], ParAtt[PAR_POSZ][p] };
         const real   MassFac = 1.0 - RNG->GetValue( TID, Plummer_FB_ExpMMin, Plummer_FB_ExpMMax );
         const real   EngyFac = 1.0 + RNG->GetValue( TID, Plummer_FB_ExpEMin, Plummer_FB_ExpEMax );

         int idx[3];
         for (int d=0; d<3; d++)    idx[d] = (int)floor( ( xyz[d] - EdgeL[d] )*_dh );

         if ( OPT__VERBOSE ) {
            Aux_Message( stdout, "\n" );
            Aux_Message( stdout, "---------------------------------------------------\n" );
            Aux_Message( stdout, "Rank %d, Time %21.14e, ParID % 5d: MassFac %13.7e, EngyFac %13.7e\n",
                         MPI_Rank, TimeNew, p, MassFac, EngyFac );
            Aux_Message( stdout, "---------------------------------------------------\n" );
         }

//       explosion
         if ( Plummer_FB_Exp )
         {
//          inject explosion energy to the nearby 27 cells including itself
//          --> skip cells in the ghost zones
//          --> should also skip cells too close to the coarse-fine boundaries using CoarseFine[] (not implemented here)
            for (int dk=-1; dk<=1; dk++)  {  const int k = idx[2] + dk;  if ( k < FB_GHOST_SIZE || k >= FB_GHOST_SIZE+PS2 )   continue;
            for (int dj=-1; dj<=1; dj++)  {  const int j = idx[1] + dj;  if ( j < FB_GHOST_SIZE || j >= FB_GHOST_SIZE+PS2 )   continue;
            for (int di=-1; di<=1; di++)  {  const int i = idx[0] + di;  if ( i < FB_GHOST_SIZE || i >= FB_GHOST_SIZE+PS2 )   continue;

               Fluid[ENGY][k][j][i] *= EngyFac;

            }}}

//          reduce particle mass
//          --> skip particles in the ghost zones
            if ( idx[0] >= FB_GHOST_SIZE  &&  idx[0] < FB_GHOST_SIZE+PS2  &&
                 idx[1] >= FB_GHOST_SIZE  &&  idx[1] < FB_GHOST_SIZE+PS2  &&
                 idx[2] >= FB_GHOST_SIZE  &&  idx[2] < FB_GHOST_SIZE+PS2   )
               ParAtt[PAR_MASS][p] *= MassFac;
         } // if ( Plummer_FB_Exp )

//       mass accretion
         if ( Plummer_FB_Acc )
         {
//          --> should also skip cells too close to the coarse-fine boundaries using CoarseFine[] (not implemented here)
         } // if ( Plummer_FB_Acc )
      } // if ( RNG->GetValue(TID,0.0,1.0) < Plummer_FB_Like )
   } // for (int t=0; t<NPar; t++)

} // FUNCTION : FB_Plummer



//-------------------------------------------------------------------------------------------------------
// Function    :  FB_End_Plummer
// Description :  Free the resources used by the user-specified feedback
//
// Note        :  1. Invoked by FB_End()
//                2. Linked to FB_End_User_Ptr in FB_Init_Plummer()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void FB_End_Plummer()
{


} // FUNCTION : FB_End_Plummer



//-------------------------------------------------------------------------------------------------------
// Function    :  FB_Init_Plummer
// Description :  Initialize the user-specified feedback
//
// Note        :  1. Invoked by FB_Init()
//                   --> Enable it by linking to the function pointer "FB_Init_User_Ptr"
//                2. Set FB_User_Ptr and FB_End_User_Ptr
//
// Parameter   :  None
//
// Return      :  FB_User_Ptr and FB_End_User_Ptr
//-------------------------------------------------------------------------------------------------------
void FB_Init_Plummer()
{

   FB_User_Ptr     = FB_Plummer;
   FB_End_User_Ptr = FB_End_Plummer;

} // FUNCTION : FB_Init_Plummer



#endif // #ifdef FEEDBACK
