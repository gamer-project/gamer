#include "GAMER.h"
#include "CUPOT.h"

#if ( defined GRAVITY  &&  !defined GPU  &&  POT_SCHEME == MG )



#define POT_NXT_INT  ( (POT_NXT-2)*2    )    // size of the array "Pot_Array_Int"
#define POT_USELESS  ( POT_GHOST_SIZE%2 )    // # of useless cells in each side of the array "Pot_Array_Int"
#define MAX_NLV      10                      // maximum number of multigrid levels

static void Smoothing( real *Sol_1D, const real *RHS_1D, const real dh, const int NGrid );
static void ComputeDefect( const real *Sol_1D, const real *RHS_1D, real *Def_1D, const real dh,
                           const int NGrid, const bool EstimateError, real &Error );
static void Restrict( const real *FData_1D, real *CData_1D, const int NGrid_F, const int NGrid_C );
static void Prolongate_and_Correct( const real *CData_1D, real *FData_1D, const int NGrid_C,
                                    const int NGrid_F );




//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_PoissonSolver_MG
// Description :  Use CPU to solve the Poisson equation by the multigrid scheme
//
// Note        :  Reference : Numerical Recipes, Chapter 20.6
//
// Parameter   :  Rho_Array         : Array to store the input density
//                Pot_Array_In      : Array to store the input "coarse-grid" potential for interpolation
//                Pot_Array_Out     : Array to store the output potential
//                NPatchGroup       : Number of patch groups evaluated at a time
//                dh_Min            : Grid size of the input data
//                Max_Iter          : Maximum number of iterations for multigrid
//                NPre_Smooth       : Number of pre-smoothing steps for multigrid
//                NPost_Smooth      : Number of post-smoothing steps for multigrid
//                Tolerated_Error   : Maximum tolerated error for multigrid
//                Poi_Coeff         : Coefficient in front of the RHS in the Poisson eq.
//                IntScheme         : Interpolation scheme for potential
//                                    --> currently supported schemes include
//                                        INT_CQUAD : conservative quadratic interpolation
//                                        INT_QUAD  : quadratic interpolation
//-------------------------------------------------------------------------------------------------------
void CPU_PoissonSolver_MG( const real Rho_Array    [][RHO_NXT][RHO_NXT][RHO_NXT],
                           const real Pot_Array_In [][POT_NXT][POT_NXT][POT_NXT],
                                 real Pot_Array_Out[][GRA_NXT][GRA_NXT][GRA_NXT],
                           const int NPatchGroup, const real dh_Min, const int Max_Iter, const int NPre_Smooth,
                           const int NPost_Smooth, const real Tolerated_Error, const real Poi_Coeff,
                           const IntScheme_t IntScheme )
{

   const int  NPatch    = NPatchGroup*8;
   const real Const_8   = (real)1.0/(real)  8.0;
   const real Const_64  = (real)1.0/(real) 64.0;
   const real Const_512 = (real)1.0/(real)512.0;
   const real Mp[3]     = { (real)-3.0/32.0, (real)+30.0/32.0, (real)+5.0/32.0 };
   const real Mm[3]     = { (real)+5.0/32.0, (real)+30.0/32.0, (real)-3.0/32.0 };


// set the depth of multigrid V-cycle
   int  BottomLv, NGrid[MAX_NLV], NBottom_Smooth;
   real dh[MAX_NLV];

   BottomLv = 0;
   NGrid[0] = PATCH_SIZE + 2*POT_GHOST_SIZE;
   dh   [0] = dh_Min;

   while ( (NGrid[BottomLv]+1)/2 >= 3  &&  BottomLv+1<MAX_NLV )
   {
      BottomLv ++;
      NGrid[BottomLv] = ( NGrid[BottomLv-1] + 1 ) / 2;
      dh   [BottomLv] = dh_Min * ( NGrid[0] - 1 ) / ( NGrid[BottomLv] - 1 );
   }

#  ifdef GAMER_DEBUG
   if ( BottomLv == MAX_NLV-1 )  Aux_Message( stderr, "WARNING : BottomLv == MAX_NLV-1 (%d) !!\n", MAX_NLV-1 );
#  endif

   switch ( NGrid[BottomLv] )
   {
      case 3:  NBottom_Smooth = 1;  break;
      case 4:  NBottom_Smooth = 7;  break;
      default:
         Aux_Error( ERROR_INFO, "NGrid at the bottom level != 3 or 4 --> Please specify NBottom_Smooth !!\n" );
   }


#  pragma omp parallel
   {
      int ip, jp, kp, im, jm, km, I, J, K, Ip, Jp, Kp, ii, jj, kk, Iter, x, y, z, Count, Idx;
      real Slope_x, Slope_y, Slope_z, C2_Slope[13], Error;

//    multigrid arrays
      real (**Sol) = new real* [BottomLv+1];    // solution
      real (**RHS) = new real* [BottomLv+1];    // right-hand-side
      real (**Def) = new real* [BottomLv+1];    // defect

      for (int Lv=0; Lv<=BottomLv; Lv++)
      {
         Sol[Lv] = new real [ NGrid[Lv]*NGrid[Lv]*NGrid[Lv] ];
         RHS[Lv] = new real [ NGrid[Lv]*NGrid[Lv]*NGrid[Lv] ];
         Def[Lv] = new real [ NGrid[Lv]*NGrid[Lv]*NGrid[Lv] ];
      }

//    array to store the interpolated "fine-grid" potential (as the initial guess and the B.C.)
      real (*Pot_Array_Int)[POT_NXT_INT][POT_NXT_INT];
      if ( POT_USELESS == 0 )    Pot_Array_Int = ( real(*)[POT_NXT_INT][POT_NXT_INT] )Sol[0];
      else                       Pot_Array_Int = new real [POT_NXT_INT][POT_NXT_INT][POT_NXT_INT];

//    initialize Def as zero (actually we only need to set boundary values as zero)
      for (int Lv=0; Lv<=BottomLv; Lv++)
      for (int t=0; t<NGrid[Lv]*NGrid[Lv]*NGrid[Lv]; t++)
         Def[Lv][t] = (real)0.0;


//    loop over all patches
#     pragma omp for schedule( runtime )
      for (int P=0; P<NPatch; P++)
      {

//       a. interpolation : Pot_Array_In --> Pot_Array_Int
// ------------------------------------------------------------------------------------------------------------
         switch ( IntScheme )
         {
            /*
            case INT_CENTRAL :
            {
               for (int k=1; k<POT_NXT-1; k++)  {  K = (k-1)*2;   Kp = K + 1;    kp = k + 1;    km = k - 1;
               for (int j=1; j<POT_NXT-1; j++)  {  J = (j-1)*2;   Jp = J + 1;    jp = j + 1;    jm = j - 1;
               for (int i=1; i<POT_NXT-1; i++)  {  I = (i-1)*2;   Ip = I + 1;    ip = i + 1;    im = i - 1;

                  Slope_x = (real)0.125 * ( Pot_Array_In[P][k ][j ][ip] - Pot_Array_In[P][k ][j ][im] );
                  Slope_y = (real)0.125 * ( Pot_Array_In[P][k ][jp][i ] - Pot_Array_In[P][k ][jm][i ] );
                  Slope_z = (real)0.125 * ( Pot_Array_In[P][kp][j ][i ] - Pot_Array_In[P][km][j ][i ] );

                  Pot_Array_Int[K ][J ][I ] = Pot_Array_In[P][k][j][i] - Slope_z - Slope_y - Slope_x;
                  Pot_Array_Int[K ][J ][Ip] = Pot_Array_In[P][k][j][i] - Slope_z - Slope_y + Slope_x;
                  Pot_Array_Int[K ][Jp][I ] = Pot_Array_In[P][k][j][i] - Slope_z + Slope_y - Slope_x;
                  Pot_Array_Int[K ][Jp][Ip] = Pot_Array_In[P][k][j][i] - Slope_z + Slope_y + Slope_x;
                  Pot_Array_Int[Kp][J ][I ] = Pot_Array_In[P][k][j][i] + Slope_z - Slope_y - Slope_x;
                  Pot_Array_Int[Kp][J ][Ip] = Pot_Array_In[P][k][j][i] + Slope_z - Slope_y + Slope_x;
                  Pot_Array_Int[Kp][Jp][I ] = Pot_Array_In[P][k][j][i] + Slope_z + Slope_y - Slope_x;
                  Pot_Array_Int[Kp][Jp][Ip] = Pot_Array_In[P][k][j][i] + Slope_z + Slope_y + Slope_x;

               }}}
            }
            break; // INT_CENTRAL
            */


            case INT_CQUAD :
            {
               for (int k=1; k<POT_NXT-1; k++)  {  K = (k-1)*2;   Kp = K + 1;    kp = k + 1;    km = k - 1;
               for (int j=1; j<POT_NXT-1; j++)  {  J = (j-1)*2;   Jp = J + 1;    jp = j + 1;    jm = j - 1;
               for (int i=1; i<POT_NXT-1; i++)  {  I = (i-1)*2;   Ip = I + 1;    ip = i + 1;    im = i - 1;

                  C2_Slope[ 0] = Const_8   * ( Pot_Array_In[P][k ][j ][ip] - Pot_Array_In[P][k ][j ][im] );
                  C2_Slope[ 1] = Const_8   * ( Pot_Array_In[P][k ][jp][i ] - Pot_Array_In[P][k ][jm][i ] );
                  C2_Slope[ 2] = Const_8   * ( Pot_Array_In[P][kp][j ][i ] - Pot_Array_In[P][km][j ][i ] );

                  C2_Slope[ 3] = Const_64  * ( Pot_Array_In[P][km][j ][ip] - Pot_Array_In[P][km][j ][im] );
                  C2_Slope[ 4] = Const_64  * ( Pot_Array_In[P][km][jp][i ] - Pot_Array_In[P][km][jm][i ] );
                  C2_Slope[ 5] = Const_64  * ( Pot_Array_In[P][k ][jm][ip] - Pot_Array_In[P][k ][jm][im] );
                  C2_Slope[ 6] = Const_64  * ( Pot_Array_In[P][k ][jp][ip] - Pot_Array_In[P][k ][jp][im] );
                  C2_Slope[ 7] = Const_64  * ( Pot_Array_In[P][kp][j ][ip] - Pot_Array_In[P][kp][j ][im] );
                  C2_Slope[ 8] = Const_64  * ( Pot_Array_In[P][kp][jp][i ] - Pot_Array_In[P][kp][jm][i ] );

                  C2_Slope[ 9] = Const_512 * ( Pot_Array_In[P][km][jm][ip] - Pot_Array_In[P][km][jm][im] );
                  C2_Slope[10] = Const_512 * ( Pot_Array_In[P][km][jp][ip] - Pot_Array_In[P][km][jp][im] );
                  C2_Slope[11] = Const_512 * ( Pot_Array_In[P][kp][jm][ip] - Pot_Array_In[P][kp][jm][im] );
                  C2_Slope[12] = Const_512 * ( Pot_Array_In[P][kp][jp][ip] - Pot_Array_In[P][kp][jp][im] );


                  Pot_Array_Int[K ][J ][I ] = - C2_Slope[ 0] - C2_Slope[ 1] - C2_Slope[ 2] - C2_Slope[ 3]
                                              - C2_Slope[ 4] - C2_Slope[ 5] + C2_Slope[ 6] + C2_Slope[ 7]
                                              + C2_Slope[ 8] - C2_Slope[ 9] + C2_Slope[10] + C2_Slope[11]
                                              - C2_Slope[12] + Pot_Array_In[P][k][j][i];

                  Pot_Array_Int[K ][J ][Ip] = + C2_Slope[ 0] - C2_Slope[ 1] - C2_Slope[ 2] + C2_Slope[ 3]
                                              - C2_Slope[ 4] + C2_Slope[ 5] - C2_Slope[ 6] - C2_Slope[ 7]
                                              + C2_Slope[ 8] + C2_Slope[ 9] - C2_Slope[10] - C2_Slope[11]
                                              + C2_Slope[12] + Pot_Array_In[P][k][j][i];

                  Pot_Array_Int[K ][Jp][I ] = - C2_Slope[ 0] + C2_Slope[ 1] - C2_Slope[ 2] - C2_Slope[ 3]
                                              + C2_Slope[ 4] + C2_Slope[ 5] - C2_Slope[ 6] + C2_Slope[ 7]
                                              - C2_Slope[ 8] + C2_Slope[ 9] - C2_Slope[10] - C2_Slope[11]
                                              + C2_Slope[12] + Pot_Array_In[P][k][j][i];

                  Pot_Array_Int[K ][Jp][Ip] = + C2_Slope[ 0] + C2_Slope[ 1] - C2_Slope[ 2] + C2_Slope[ 3]
                                              + C2_Slope[ 4] - C2_Slope[ 5] + C2_Slope[ 6] - C2_Slope[ 7]
                                              - C2_Slope[ 8] - C2_Slope[ 9] + C2_Slope[10] + C2_Slope[11]
                                              - C2_Slope[12] + Pot_Array_In[P][k][j][i];

                  Pot_Array_Int[Kp][J ][I ] = - C2_Slope[ 0] - C2_Slope[ 1] + C2_Slope[ 2] + C2_Slope[ 3]
                                              + C2_Slope[ 4] - C2_Slope[ 5] + C2_Slope[ 6] - C2_Slope[ 7]
                                              - C2_Slope[ 8] + C2_Slope[ 9] - C2_Slope[10] - C2_Slope[11]
                                              + C2_Slope[12] + Pot_Array_In[P][k][j][i];

                  Pot_Array_Int[Kp][J ][Ip] = + C2_Slope[ 0] - C2_Slope[ 1] + C2_Slope[ 2] - C2_Slope[ 3]
                                              + C2_Slope[ 4] + C2_Slope[ 5] - C2_Slope[ 6] + C2_Slope[ 7]
                                              - C2_Slope[ 8] - C2_Slope[ 9] + C2_Slope[10] + C2_Slope[11]
                                              - C2_Slope[12] + Pot_Array_In[P][k][j][i];

                  Pot_Array_Int[Kp][Jp][I ] = - C2_Slope[ 0] + C2_Slope[ 1] + C2_Slope[ 2] + C2_Slope[ 3]
                                              - C2_Slope[ 4] + C2_Slope[ 5] - C2_Slope[ 6] - C2_Slope[ 7]
                                              + C2_Slope[ 8] - C2_Slope[ 9] + C2_Slope[10] + C2_Slope[11]
                                              - C2_Slope[12] + Pot_Array_In[P][k][j][i];

                  Pot_Array_Int[Kp][Jp][Ip] = + C2_Slope[ 0] + C2_Slope[ 1] + C2_Slope[ 2] - C2_Slope[ 3]
                                              - C2_Slope[ 4] - C2_Slope[ 5] + C2_Slope[ 6] + C2_Slope[ 7]
                                              + C2_Slope[ 8] + C2_Slope[ 9] - C2_Slope[10] - C2_Slope[11]
                                              + C2_Slope[12] + Pot_Array_In[P][k][j][i];
               }}} // i, j, k
            }
            break; // INT_CQUAD


            case INT_QUAD :
            {
               for (int k=0; k<POT_NXT_INT; k++)
               for (int j=0; j<POT_NXT_INT; j++)
               for (int i=0; i<POT_NXT_INT; i++)   Pot_Array_Int[k][j][i] = (real)0.0;

               for (int k=1; k<POT_NXT-1; k++)  {  K = (k-1)*2;   Kp = K + 1;
               for (int j=1; j<POT_NXT-1; j++)  {  J = (j-1)*2;   Jp = J + 1;
               for (int i=1; i<POT_NXT-1; i++)  {  I = (i-1)*2;   Ip = I + 1;

                  for (int dk=-1; dk<=1; dk++)  {  z = dk+1;  kk = k + dk;
                  for (int dj=-1; dj<=1; dj++)  {  y = dj+1;  jj = j + dj;
                  for (int di=-1; di<=1; di++)  {  x = di+1;  ii = i + di;

                     Pot_Array_Int[K ][J ][I ] += Pot_Array_In[P][kk][jj][ii] * Mm[z] * Mm[y] * Mm[x];
                     Pot_Array_Int[K ][J ][Ip] += Pot_Array_In[P][kk][jj][ii] * Mm[z] * Mm[y] * Mp[x];
                     Pot_Array_Int[K ][Jp][I ] += Pot_Array_In[P][kk][jj][ii] * Mm[z] * Mp[y] * Mm[x];
                     Pot_Array_Int[K ][Jp][Ip] += Pot_Array_In[P][kk][jj][ii] * Mm[z] * Mp[y] * Mp[x];
                     Pot_Array_Int[Kp][J ][I ] += Pot_Array_In[P][kk][jj][ii] * Mp[z] * Mm[y] * Mm[x];
                     Pot_Array_Int[Kp][J ][Ip] += Pot_Array_In[P][kk][jj][ii] * Mp[z] * Mm[y] * Mp[x];
                     Pot_Array_Int[Kp][Jp][I ] += Pot_Array_In[P][kk][jj][ii] * Mp[z] * Mp[y] * Mm[x];
                     Pot_Array_Int[Kp][Jp][Ip] += Pot_Array_In[P][kk][jj][ii] * Mp[z] * Mp[y] * Mp[x];

                  }}}
               }}} // i, j, k
            }
            break; // INT_QUAD


            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "IntScheme", IntScheme );

         } // switch ( IntScheme )



//       b. initialize the level 0 multigrid arrays
// ------------------------------------------------------------------------------------------------------------
         if ( POT_USELESS != 0 )
         {
            Count = 0;
            for (int k=POT_USELESS; k<POT_NXT_INT-POT_USELESS; k++)
            for (int j=POT_USELESS; j<POT_NXT_INT-POT_USELESS; j++)
            for (int i=POT_USELESS; i<POT_NXT_INT-POT_USELESS; i++)
               Sol[0][ Count ++ ] = Pot_Array_Int[k][j][i];
         }

         for (int k=0; k<RHO_NXT; k++)    {  kp = k + 1;
         for (int j=0; j<RHO_NXT; j++)    {  jp = j + 1;
         for (int i=0; i<RHO_NXT; i++)    {  ip = i + 1;

            Idx         = (kp*NGrid[0] + jp)*NGrid[0] + ip;
            RHS[0][Idx] = Poi_Coeff*Rho_Array[P][k][j][i];

         }}}



//       c. use the MG scheme to evaluate potential
// ------------------------------------------------------------------------------------------------------------
         Iter  = 0;
         Error = __FLT_MAX__;

         while ( Iter < Max_Iter  &&  Error > Tolerated_Error )
         {
//          V-cycle : finer --> coarser grids
            for (int Lv=0; Lv<BottomLv; Lv++)
            {
//             pre-smoothing (apply relaxation to compute solution/correction)
               for (int PreStep=0; PreStep<NPre_Smooth; PreStep++)
               Smoothing( Sol[Lv], RHS[Lv], dh[Lv], NGrid[Lv] );

//             compute defect
               ComputeDefect( Sol[Lv], RHS[Lv], Def[Lv], dh[Lv], NGrid[Lv], false, Error );

//             restrict defect (use as the RHS at the next level)
               Restrict( Def[Lv], RHS[Lv+1], NGrid[Lv], NGrid[Lv+1] );

//             initialize the correction at the next level to zero
               for (int t=0; t<NGrid[Lv+1]*NGrid[Lv+1]*NGrid[Lv+1]; t++)   Sol[Lv+1][t] = (real)0.0;
            }


//          calculate the correction at the bottom level
            for (int BottomStep=0; BottomStep<NBottom_Smooth; BottomStep++)
            Smoothing( Sol[BottomLv], RHS[BottomLv], dh[BottomLv], NGrid[BottomLv] );


//          V-cycle : coarser --> finer grids
            for (int Lv=BottomLv-1; Lv>=0; Lv--)
            {
//             prolongate correction (from Lv+1 to Lv) and correct solution/correction at Lv
               Prolongate_and_Correct( Sol[Lv+1], Sol[Lv], NGrid[Lv+1], NGrid[Lv] );

//             post-smoothing (apply relaxation to compute solution/correction again)
               for (int PostStep=0; PostStep<NPost_Smooth; PostStep++)
               Smoothing( Sol[Lv], RHS[Lv], dh[Lv], NGrid[Lv] );
            }

//          estimate error
            ComputeDefect( Sol[0], RHS[0], Def[0], dh[0], NGrid[0], true, Error );
            Iter ++;

//          Aux_Message( stdout, "Patch %3d, Iter %3d: Error = %13.7e\n", P, Iter, Error );

         } // while ( Iter < Max_Iter  &&  Error > Tolerated_Error )


//       Aux_Message( stdout, "Patch %3d : number of iterations = %3d\n", P, Iter );

         if ( Error > Tolerated_Error )
         {
            Aux_Message( stderr, "WARNING : Rank = %2d, Patch %6d exceeds the maximum tolerated error ",
                         MPI_Rank, P );
            Aux_Message( stderr, "(error = %13.7e)\n", Error );
         }


//       d. copy data : Sol[0] --> Pot_Array_Out
// ------------------------------------------------------------------------------------------------------------
         real (*Pot_3D)[RHO_NXT+2][RHO_NXT+2] = ( real(*)[RHO_NXT+2][RHO_NXT+2] )Sol[0];

         for (int k=0; k<GRA_NXT; k++)    {  K = k + POT_GHOST_SIZE - GRA_GHOST_SIZE;
         for (int j=0; j<GRA_NXT; j++)    {  J = j + POT_GHOST_SIZE - GRA_GHOST_SIZE;
         for (int i=0; i<GRA_NXT; i++)    {  I = i + POT_GHOST_SIZE - GRA_GHOST_SIZE;

            Pot_Array_Out[P][k][j][i] = Pot_3D[K][J][I];

         }}}

      } // for (int P=0; P<NPatch; P++)


//    free memory
      for (int Lv=0; Lv<=BottomLv; Lv++)
      {
         delete [] Sol[Lv];
         delete [] RHS[Lv];
         delete [] Def[Lv];
      }
      if ( POT_USELESS != 0 )    delete [] Pot_Array_Int;
      delete [] Sol;
      delete [] RHS;
      delete [] Def;

   } // OpenMP parallel region

} // FUNCTION : CPU_PoissonSolver_MG



//-------------------------------------------------------------------------------------------------------
// Function    :  Smoothing
// Description :  Use Gauss-Seidel method for smoothing
//
// Note        :  B.C. should be stored in the input array "Sol"
//
// Parameter   :  Sol_1D   : 1D array to store the output solution to the Poisson equation
//                RHS_1D   : 1D array storing the RHS of the Poisson equation
//                dh       : Grid size
//                NGrid    : Number of cells in each spatial direction for Sol and RHS
//-------------------------------------------------------------------------------------------------------
void Smoothing( real *Sol_1D, const real *RHS_1D, const real dh, const int NGrid )
{

   const real dh2     = dh*dh;
   const real One_Six = (real)1.0/(real)6.0;

   int i_start, i_start_pass, i_start_k;     // i_start_(pass,k) : record the i_start in the (pass,k) loop
   int ip, jp, kp, im, jm, km;

// typecasting to 3D arrays
   typedef real (*vla)[NGrid][NGrid];
         vla Sol = ( vla )Sol_1D;
   const vla RHS = ( vla )RHS_1D;


// odd-even ordering
   i_start_pass = 1;

   for (int pass=0; pass<2; pass++)
   {
      i_start_k = i_start_pass;

      for (int k=1; k<NGrid-1; k++)
      {
         i_start = i_start_k;
         kp      = k+1;
         km      = k-1;

         for (int j=1; j<NGrid-1; j++)
         {
            jp = j+1;
            jm = j-1;

            for (int i=i_start; i<NGrid-1; i+=2)
            {
               ip = i+1;
               im = i-1;

//             update solution
               Sol[k][j][i] = One_Six*(   Sol[kp][j ][i ] + Sol[km][j ][i ]
                                        + Sol[k ][jp][i ] + Sol[k ][jm][i ]
                                        + Sol[k ][j ][ip] + Sol[k ][j ][im] - dh2*RHS[k][j][i]  );
            } // i
            i_start = 3 - i_start;
         } // j
         i_start_k = 3 - i_start_k;
      } // k
      i_start_pass = 3 - i_start_pass;
   } // for (int pass=0; pass<2; pass++)

} // FUNCTION : Smoothing



//-------------------------------------------------------------------------------------------------------
// Function    :  ComputeDefect
// Description :  Compute negative defect defined as "-(Laplacian(Sol)-RHS)"
//
// Note        :  1. B.C. should be stored in the input array "Sol"
//                2. It is assumed that the boundary values of Def_1D have already been initialized as zero
//
// Parameter   :  Sol_1D         : 1D array storing the input solution to the Poisson equation
//                RHS_1D         : 1D array storing the RHS of the Poisson equation
//                Def_1D         : 1D array to store the output defect
//                dh             : Grid size
//                NGrid          : Number of cells in each spatial direction for Sol and RHS
//                EstimateError  : estimate the L1 error
//                ERROR          : L1 error to be returned
//-------------------------------------------------------------------------------------------------------
void ComputeDefect( const real *Sol_1D, const real *RHS_1D, real *Def_1D, const real dh,
                    const int NGrid, const bool EstimateError, real &Error )
{

   const real _dh2 = (real)-1.0/(dh*dh);

   int ip, jp, kp, im, jm, km;

// typecasting to 3D arrays
   typedef real (*vla)[NGrid][NGrid];
   const vla Sol = ( vla )Sol_1D;
   const vla RHS = ( vla )RHS_1D;
         vla Def = ( vla )Def_1D;


   for (int k=1; k<NGrid-1; k++)
   {
      kp = k + 1;
      km = k - 1;

      for (int j=1; j<NGrid-1; j++)
      {
         jp = j + 1;
         jm = j - 1;

         for (int i=1; i<NGrid-1; i++)
         {
            ip = i + 1;
            im = i - 1;

            Def[k][j][i] = _dh2*(   Sol[kp][j ][i ] + Sol[km][j ][i ]
                                  + Sol[k ][jp][i ] + Sol[k ][jm][i ]
                                  + Sol[k ][j ][ip] + Sol[k ][j ][im] - (real)6.0*Sol[k][j][i] ) + RHS[k][j][i];
         } // i
      } // j
   } // k


// estimate the L1 error
   if ( EstimateError )
   {
      real SumDef = (real)0.0;
      real SumSol = (real)0.0;

      for (int k=1; k<NGrid-1; k++)
      for (int j=1; j<NGrid-1; j++)
      for (int i=1; i<NGrid-1; i++)
      {
         SumDef += FABS( Def[k][j][i] );
         SumSol += FABS( Sol[k][j][i] );
      }

      Error = dh*dh*SumDef/SumSol;
   }

} // FUNCTION : ComputeDefect



//-------------------------------------------------------------------------------------------------------
// Function    :  Restrict
// Description :  Restrict the input fine-grid data to get the coarse-grid data
//
// Note        :  1. It is assumed that the input arrays follow the "finite-difference" fashion, in which the data
//                   are defined in the cell intersections instead of cell averages
//                   --> N^3 cells define a 3D grid with the size equal to (N-1)^3
//                2. Fine-grid and coarse-grid data at boundaries are assumed to be zero (because defect at
//                   boundaries are always zero)
//
// Parameter   :  FData_1D : 1D array storing the input fine-grid data
//                CData_1D : 1D array to store the output coarse-grid data
//                NGrid_F  : Number of fine-grid cells in each spatial direction
//                NGrid_C  : Number of coarse-grid cells in each spatial direction
//-------------------------------------------------------------------------------------------------------
void Restrict( const real *FData_1D, real *CData_1D, const int NGrid_F, const int NGrid_C )
{

   const real Ratio = real(NGrid_F-1) / real(NGrid_C-1);

// typecasting to 3D arrays
   typedef real (*vlaf)[NGrid_F][NGrid_F];
   typedef real (*vlac)[NGrid_C][NGrid_C];
   const vlaf FData = ( vlaf )FData_1D;
         vlac CData = ( vlac )CData_1D;

   int  iim, ii, iip, jjm, jj, jjp, kkm, kk, kkp;
   real x, y, z, Coeff[3][3]; // Coeff[x/y/z][left/center/right]


   for (int k=1; k<NGrid_C-1; k++)
   {
      z           = k*Ratio;
      kk          = int( z + (real)0.5 );
      kkm         = kk - 1;
      kkp         = kk + 1;
      Coeff[2][0] = (real)0.5 * ( kk + (real)0.5 - z );
      Coeff[2][2] = (real)0.5 - Coeff[2][0];
      Coeff[2][1] = (real)0.5;
//    Coeff[2][1] = (real)1.0 - Coeff[2][0] - Coeff[2][2];

      for (int j=1; j<NGrid_C-1; j++)
      {
         y           = j*Ratio;
         jj          = int( y + (real)0.5 );
         jjm         = jj - 1;
         jjp         = jj + 1;
         Coeff[1][0] = (real)0.5 * ( jj + (real)0.5 - y );
         Coeff[1][2] = (real)0.5 - Coeff[1][0];
         Coeff[1][1] = (real)0.5;
//       Coeff[1][1] = (real)1.0 - Coeff[1][0] - Coeff[1][2];

         for (int i=1; i<NGrid_C-1; i++)
         {
            x           = i*Ratio;
            ii          = int( x + (real)0.5 );
            iim         = ii - 1;
            iip         = ii + 1;
            Coeff[0][0] = (real)0.5 * ( ii + (real)0.5 - x );
            Coeff[0][2] = (real)0.5 - Coeff[0][0];
            Coeff[0][1] = (real)0.5;
//          Coeff[0][1] = (real)1.0 - Coeff[0][0] - Coeff[0][2];

//###OPTIMIZATION: follow the same strategy adopted in "Int_Quadratic"
            CData[k][j][i] =    Coeff[2][0] * Coeff[1][0] * Coeff[0][0] * FData[kkm][jjm][iim]
                             +  Coeff[2][0] * Coeff[1][0] * Coeff[0][1] * FData[kkm][jjm][ii ]
                             +  Coeff[2][0] * Coeff[1][0] * Coeff[0][2] * FData[kkm][jjm][iip]
                             +  Coeff[2][0] * Coeff[1][1] * Coeff[0][0] * FData[kkm][jj ][iim]
                             +  Coeff[2][0] * Coeff[1][1] * Coeff[0][1] * FData[kkm][jj ][ii ]
                             +  Coeff[2][0] * Coeff[1][1] * Coeff[0][2] * FData[kkm][jj ][iip]
                             +  Coeff[2][0] * Coeff[1][2] * Coeff[0][0] * FData[kkm][jjp][iim]
                             +  Coeff[2][0] * Coeff[1][2] * Coeff[0][1] * FData[kkm][jjp][ii ]
                             +  Coeff[2][0] * Coeff[1][2] * Coeff[0][2] * FData[kkm][jjp][iip]

                             +  Coeff[2][1] * Coeff[1][0] * Coeff[0][0] * FData[kk ][jjm][iim]
                             +  Coeff[2][1] * Coeff[1][0] * Coeff[0][1] * FData[kk ][jjm][ii ]
                             +  Coeff[2][1] * Coeff[1][0] * Coeff[0][2] * FData[kk ][jjm][iip]
                             +  Coeff[2][1] * Coeff[1][1] * Coeff[0][0] * FData[kk ][jj ][iim]
                             +  Coeff[2][1] * Coeff[1][1] * Coeff[0][1] * FData[kk ][jj ][ii ]
                             +  Coeff[2][1] * Coeff[1][1] * Coeff[0][2] * FData[kk ][jj ][iip]
                             +  Coeff[2][1] * Coeff[1][2] * Coeff[0][0] * FData[kk ][jjp][iim]
                             +  Coeff[2][1] * Coeff[1][2] * Coeff[0][1] * FData[kk ][jjp][ii ]
                             +  Coeff[2][1] * Coeff[1][2] * Coeff[0][2] * FData[kk ][jjp][iip]

                             +  Coeff[2][2] * Coeff[1][0] * Coeff[0][0] * FData[kkp][jjm][iim]
                             +  Coeff[2][2] * Coeff[1][0] * Coeff[0][1] * FData[kkp][jjm][ii ]
                             +  Coeff[2][2] * Coeff[1][0] * Coeff[0][2] * FData[kkp][jjm][iip]
                             +  Coeff[2][2] * Coeff[1][1] * Coeff[0][0] * FData[kkp][jj ][iim]
                             +  Coeff[2][2] * Coeff[1][1] * Coeff[0][1] * FData[kkp][jj ][ii ]
                             +  Coeff[2][2] * Coeff[1][1] * Coeff[0][2] * FData[kkp][jj ][iip]
                             +  Coeff[2][2] * Coeff[1][2] * Coeff[0][0] * FData[kkp][jjp][iim]
                             +  Coeff[2][2] * Coeff[1][2] * Coeff[0][1] * FData[kkp][jjp][ii ]
                             +  Coeff[2][2] * Coeff[1][2] * Coeff[0][2] * FData[kkp][jjp][iip];

//          coefficient adopted in Enzo which seems to give faster convergence rate
//          CData[k][j][i] *= (real)0.52*Ratio;
         } // i
      } // j
   } // k

} // FUNCTION : Restrict



//-------------------------------------------------------------------------------------------------------
// Function    :  Prolongate_and_Correct
// Description :  Prolongate the input coarse-grid correction to correct the fine-grid solution/correction
//
// Note        :  1. It is assumed that the input arrays follow the "finite-difference" fashion, in which the data
//                   are defined in the cell intersections instead of cell averages
//                   --> N^3 cells define a 3D grid with the size equal to (N-1)^3
//                2. Boundary data of FData_1D are not corrected (since solution/correction at boundaries
//                   should be fixed
//
// Parameter   :  CData_1D : 1D array storing the input coarse-grid data
//                FData_1D : 1D array to store the output fine-grid data
//                NGrid_C  : Number of coarse-grid cells in each spatial direction
//                NGrid_F  : Number of fine-grid cells in each spatial direction
//-------------------------------------------------------------------------------------------------------
void Prolongate_and_Correct( const real *CData_1D, real *FData_1D, const int NGrid_C, const int NGrid_F )
{

   const real Ratio = real(NGrid_C-1) / real(NGrid_F-1);

// typecasting to 3D arrays

   typedef real (*vlaf)[NGrid_F][NGrid_F];
   typedef real (*vlac)[NGrid_C][NGrid_C];
         vlaf FData = ( vlaf )FData_1D;
   const vlac CData = ( vlac )CData_1D;

   int  ii, iip, jj, jjp, kk, kkp;
   real x, y, z, Coeff[3][2]; // Coeff[x/y/z][left/right]


   for (int k=1; k<NGrid_F-1; k++)
   {
      z           = k*Ratio;
      kk          = int(z);
      kkp         = kk + 1;
      Coeff[2][0] = kkp - z;
      Coeff[2][1] = (real)1.0 - Coeff[2][0];

      for (int j=1; j<NGrid_F-1; j++)
      {
         y           = j*Ratio;
         jj          = int(y);
         jjp         = jj + 1;
         Coeff[1][0] = jjp - y;
         Coeff[1][1] = (real)1.0 - Coeff[1][0];

         for (int i=1; i<NGrid_F-1; i++)
         {
            x           = i*Ratio;
            ii          = int(x);
            iip         = ii + 1;
            Coeff[0][0] = iip - x;
            Coeff[0][1] = (real)1.0 - Coeff[0][0];

            FData[k][j][i] +=    Coeff[2][0] * Coeff[1][0] * Coeff[0][0] * CData[kk ][jj ][ii ]
                              +  Coeff[2][0] * Coeff[1][0] * Coeff[0][1] * CData[kk ][jj ][iip]
                              +  Coeff[2][0] * Coeff[1][1] * Coeff[0][0] * CData[kk ][jjp][ii ]
                              +  Coeff[2][1] * Coeff[1][0] * Coeff[0][0] * CData[kkp][jj ][ii ]
                              +  Coeff[2][0] * Coeff[1][1] * Coeff[0][1] * CData[kk ][jjp][iip]
                              +  Coeff[2][1] * Coeff[1][1] * Coeff[0][0] * CData[kkp][jjp][ii ]
                              +  Coeff[2][1] * Coeff[1][0] * Coeff[0][1] * CData[kkp][jj ][iip]
                              +  Coeff[2][1] * Coeff[1][1] * Coeff[0][1] * CData[kkp][jjp][iip];
         } // i
      } // j
   } // k

} // FUNCTION : Prolongate_and_Correct



#endif // #if ( defined GRAVITY  &&  !defined GPU  &&  POT_SCHEME == MG )
