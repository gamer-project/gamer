#include "Copyright.h"
#include "GAMER.h"

#ifdef GRAVITY




//-------------------------------------------------------------------------------------------------------
// Function    :  Poi_Prepare_Pot
// Description :  Fill up the h_Pot_Array_P_In array with the coarse-grid (lv-1) potential for the Poisson solver 
//
// Parameter   :  lv                : Targeted refinement level
//                PrepTime          : Targeted physical time to prepare the coarse-grid data
//                h_Pot_Array_P_In  : Host array to store the prepared coarse-grid potential
//                NPG               : Number of patch groups prepared at a time
//                PID0_List         : List recording the patch indicies with LocalID==0 to be udpated
//-------------------------------------------------------------------------------------------------------
void Poi_Prepare_Pot( const int lv, const double PrepTime, real h_Pot_Array_P_In[][POT_NXT][POT_NXT][POT_NXT], 
                      const int NPG, const int *PID0_List )
{

// check
   if ( lv == 0 )    Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "lv", lv );


//###REVISE: support interpolation schemes requiring 2 ghost cells on each side
   const int IntGhost = 1;    // assuming interpolation ghost-zone == 1
   const int CGhost   = (POT_GHOST_SIZE+1)/2 + IntGhost;
   const int CWidth   = PATCH_SIZE + 2*CGhost;


// determine the parameters for the spatial and temporal interpolations
   bool IntTime, IntSg;

   if      (  Mis_CompareRealValue( PrepTime,      Time[lv-1], NULL, false )  )
   {
      IntTime = false;
      IntSg   = amr->PotSg[lv-1];
   }

   else if (  Mis_CompareRealValue( PrepTime, Time_Prev[lv-1], NULL, false )  )
   {
      IntTime = false;
      IntSg   = 1 - amr->PotSg[lv-1];
   }

   else
   {
#     ifndef INDIVIDUAL_TIMESTEP
      Aux_Error( ERROR_INFO, "lv %d: PrepTime %20.14e, Time[lv-1] %20.14e, Time_Prev[lv-1] %20.14e !!\n",
                 lv, PrepTime, Time[lv-1], Time_Prev[lv-1] );
#     endif
      IntTime = ( OPT__INT_TIME ) ? true : false;
      IntSg   = 1 - amr->PotSg[lv-1];
   }


#  pragma omp parallel
   {
//    CPot : store the coarse-grid potential for interpoaltion
      real (*CPot)[CWidth][CWidth] = new real [CWidth][CWidth][CWidth];

      int FaPID, FaSibPID, PID0;
      int N, I, J, K, I2, J2, K2, Disp_i, Disp_j, Disp_k, Disp_i2, Disp_j2, Disp_k2, Loop_i, Loop_j, Loop_k;


//    prepare the coarse-grid potential for eight patches (one patch group) at a time 
#     pragma omp for schedule( runtime )
      for (int TID=0; TID<NPG; TID++)
      {
         PID0  = PID0_List[TID];
         FaPID = amr->patch[0][lv][PID0]->father;

#        ifdef GAMER_DEBUG
         if ( FaPID < 0 )  Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "FaPID", FaPID );
#        endif

//       a. fill up the central region (ghost zone is not filled up yet) of the CPot array
// ------------------------------------------------------------------------------------------------------------
         for (int k=0; k<PATCH_SIZE; k++)    {  K = k + CGhost;
         for (int j=0; j<PATCH_SIZE; j++)    {  J = j + CGhost;
         for (int i=0; i<PATCH_SIZE; i++)    {  I = i + CGhost;

            CPot[K][J][I] = amr->patch[IntSg][lv-1][FaPID]->pot[k][j][i];

         }}}

//       interpolation in time
         if ( IntTime )
         {
            for (int k=0; k<PATCH_SIZE; k++)    {  K = k + CGhost;
            for (int j=0; j<PATCH_SIZE; j++)    {  J = j + CGhost;
            for (int i=0; i<PATCH_SIZE; i++)    {  I = i + CGhost;

               CPot[K][J][I] = 0.5 * ( CPot[K][J][I] + amr->patch[!IntSg][lv-1][FaPID]->pot[k][j][i] );

            }}}
         }


//       b. fill up the ghost zone (no spatial interpolation is required) of the CPot array
// ------------------------------------------------------------------------------------------------------------
         for (int sib=0; sib<26; sib++)   // we always need coarse-grid ghost zones in 26 sibling directions
         {
            FaSibPID = amr->patch[0][lv-1][FaPID]->sibling[sib];

#           ifdef GAMER_DEBUG
            if ( FaSibPID < 0 )  Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "FaSibPID", FaSibPID );
#           endif

            Loop_i   = TABLE_01( sib, 'x', CGhost, PATCH_SIZE, CGhost );
            Loop_j   = TABLE_01( sib, 'y', CGhost, PATCH_SIZE, CGhost );
            Loop_k   = TABLE_01( sib, 'z', CGhost, PATCH_SIZE, CGhost );
            Disp_i   = TABLE_01( sib, 'x', 0, CGhost, CGhost+PATCH_SIZE );
            Disp_j   = TABLE_01( sib, 'y', 0, CGhost, CGhost+PATCH_SIZE );
            Disp_k   = TABLE_01( sib, 'z', 0, CGhost, CGhost+PATCH_SIZE );
            Disp_i2  = TABLE_01( sib, 'x', PATCH_SIZE-CGhost, 0, 0 );
            Disp_j2  = TABLE_01( sib, 'y', PATCH_SIZE-CGhost, 0, 0 );
            Disp_k2  = TABLE_01( sib, 'z', PATCH_SIZE-CGhost, 0, 0 );
           
            for (int k=0; k<Loop_k; k++)  {  K = k + Disp_k;   K2 = k + Disp_k2;
            for (int j=0; j<Loop_j; j++)  {  J = j + Disp_j;   J2 = j + Disp_j2;
            for (int i=0; i<Loop_i; i++)  {  I = i + Disp_i;   I2 = i + Disp_i2;

               CPot[K][J][I] = amr->patch[IntSg][lv-1][FaSibPID]->pot[K2][J2][I2];

            }}}

//          interpolation in time
            if ( IntTime )
            {
               for (int k=0; k<Loop_k; k++)  {  K = k + Disp_k;   K2 = k + Disp_k2;
               for (int j=0; j<Loop_j; j++)  {  J = j + Disp_j;   J2 = j + Disp_j2;
               for (int i=0; i<Loop_i; i++)  {  I = i + Disp_i;   I2 = i + Disp_i2;

                  CPot[K][J][I] = (real)0.5*( CPot[K][J][I] + 
                                              amr->patch[!IntSg][lv-1][FaSibPID]->pot[K2][J2][I2] );

               }}}

            }

         } // for (int sib=0; sib<26; sib++)


//       c. copy data from the CPot array to the h_Pot_Array_P_In array
// ------------------------------------------------------------------------------------------------------------
         for (int LocalID=0; LocalID<8; LocalID++)
         {
            N      = 8*TID + LocalID;
            Disp_i = TABLE_02( LocalID, 'x', 0, PATCH_SIZE/2 ); 
            Disp_j = TABLE_02( LocalID, 'y', 0, PATCH_SIZE/2 ); 
            Disp_k = TABLE_02( LocalID, 'z', 0, PATCH_SIZE/2 ); 

            for (int k=0; k<POT_NXT; k++)    {  K = k + Disp_k;
            for (int j=0; j<POT_NXT; j++)    {  J = j + Disp_j;
            for (int i=0; i<POT_NXT; i++)    {  I = i + Disp_i;

               h_Pot_Array_P_In[N][k][j][i] = CPot[K][J][I];

            }}}
         }

      } // for (int TID=0; TID<NPG; TID++)


      delete [] CPot;

   } // OpenMP parallel region

} // FUNCTION : Poi_Prepare_Pot



#endif // #ifdef GRAVITY
