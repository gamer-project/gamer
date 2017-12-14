#include "GAMER.h"
#include "CUPOT.h"

#if ( defined GRAVITY  &&  !defined GPU  &&  POT_SCHEME == SOR )



#define POT_NXT_INT  ( (POT_NXT-2)*2    )    // size of the array "Pot_Array_Int"
#define POT_USELESS  ( POT_GHOST_SIZE%2 )    // # of useless cells in each side of the array "Pot_Array_Int"




//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_PoissonSolver_SOR
// Description :  Use CPU to solve the Poisson equation by the SOR scheme
//
// Note        :  1. Reference : Numerical Recipes, Chapter 20.5
//                2. Typically, the number of iterations required to reach round-off errors is 20 ~ 25 (single precision)
//
// Parameter   :  Rho_Array      : Array to store the input density
//                Pot_Array_In   : Array to store the input "coarse-grid" potential for interpolation
//                Pot_Array_Out  : Array to store the output potential
//                NPatchGroup    : Number of patch groups evaluated at a time
//                dh             : Grid size
//                Min_Iter       : Minimum # of iterations for SOR
//                Max_Iter       : Maximum # of iterations for SOR
//                Omega          : Over-relaxation parameter
//                Poi_Coeff      : Coefficient in front of the RHS in the Poisson eq.
//                IntScheme      : Interpolation scheme for potential
//                                 --> currently supported schemes include
//                                     INT_CQUAD : conservative quadratic interpolation
//                                     INT_QUAD  : quadratic interpolation
//-------------------------------------------------------------------------------------------------------
void CPU_PoissonSolver_SOR( const real Rho_Array    [][RHO_NXT][RHO_NXT][RHO_NXT],
                            const real Pot_Array_In [][POT_NXT][POT_NXT][POT_NXT],
                                  real Pot_Array_Out[][GRA_NXT][GRA_NXT][GRA_NXT],
                            const int NPatchGroup, const real dh, const int Min_Iter, const int Max_Iter,
                            const real Omega, const real Poi_Coeff, const IntScheme_t IntScheme )
{

   const int  NPatch    = NPatchGroup*8;
   const real Const     = Poi_Coeff*dh*dh;
   const real Omega_6   = Omega/(real)6.0;
   const real Const_8   = (real)1.0/(real)  8.0;
   const real Const_64  = (real)1.0/(real) 64.0;
   const real Const_512 = (real)1.0/(real)512.0;
   const real Mp[3]     = { (real)-3.0/32.0, (real)+30.0/32.0, (real)+5.0/32.0 };
   const real Mm[3]     = { (real)+5.0/32.0, (real)+30.0/32.0, (real)-3.0/32.0 };

#  pragma omp parallel
   {
      int i_start, i_start_pass, i_start_k;     // i_start_(pass,k) : record the i_start in the (pass,k) loop
      int ip, jp, kp, im, jm, km, I, J, K, Ip, Jp, Kp, ii, jj, kk, Iter, x, y, z;
      real Slope_x, Slope_y, Slope_z, C2_Slope[13], Residual_Total_Old, Residual_Total, Residual;

//    array to store the interpolated "fine-grid" potential (as the initial guess and the B.C.)
      real (*Pot_Array_Int)[POT_NXT_INT][POT_NXT_INT] = new real [POT_NXT_INT][POT_NXT_INT][POT_NXT_INT];


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
               Aux_Error( ERROR_INFO, "ERROR : incorrect parameter %s = %d !!\n", "IntScheme", IntScheme );

         } // switch ( IntScheme )



//       b. use the SOR scheme to evaluate potential (store in the Pot_Array_Int array)
// ------------------------------------------------------------------------------------------------------------
         Residual_Total_Old = __FLT_MAX__;

         for (Iter=0; Iter<Max_Iter; Iter++)
         {
            Residual_Total = (real)0.0;
            i_start_pass   = 1 + POT_USELESS;

//          odd-even ordering
            for (int pass=0; pass<2; pass++)
            {
               i_start_k = i_start_pass;

               for (int k=1+POT_USELESS; k<POT_NXT_INT-1-POT_USELESS; k++)
               {
                  i_start = i_start_k;
                  kp      = k+1;
                  km      = k-1;
                  kk      = k-1-POT_USELESS;

                  for (int j=1+POT_USELESS; j<POT_NXT_INT-1-POT_USELESS; j++)
                  {
                     jp = j+1;
                     jm = j-1;
                     jj = j-1-POT_USELESS;

                     for (int i=i_start; i<POT_NXT_INT-1-POT_USELESS; i+=2)
                     {
                        ip = i+1;
                        im = i-1;
                        ii = i-1-POT_USELESS;

//                      evaluate the residual of potential
                        Residual = (             Pot_Array_Int[kp][j ][i ] + Pot_Array_Int[km][j ][i ]
                                     +           Pot_Array_Int[k ][jp][i ] + Pot_Array_Int[k ][jm][i ]
                                     +           Pot_Array_Int[k ][j ][ip] + Pot_Array_Int[k ][j ][im]
                                     - (real)6.0*Pot_Array_Int[k ][j ][i ] - Const*Rho_Array[P][kk][jj][ii]  );

//                      update potential
                        Pot_Array_Int[k][j][i] += Omega_6*Residual;

//                      sum up the 1-norm of all residuals
                        Residual_Total += FABS( Residual );

                     } // i

                     i_start = 3 - i_start + 2*POT_USELESS;

                  } // j

                  i_start_k = 3 - i_start_k + 2*POT_USELESS;

               } // k

               i_start_pass = 3 - i_start_pass + 2*POT_USELESS;

            } // for (int pass=0; pass<2; pass++)


//          terminate the SOR iteration if the total residual begins to grow
//          we set the minimum number of iterations because usually the total residual will grow at the first step
            if (  Iter+1 >= Min_Iter  &&  Residual_Total > Residual_Total_Old )
            {
               Iter++;
               break;
            }

            Residual_Total_Old = Residual_Total;

         } // for (int Iter=0; Iter<Max_Iter; Iter++)


         if ( Iter == Max_Iter )
            Aux_Message( stderr, "WARNING : Rank = %2d, Patch %6d exceeds Max_Iter in the SOR iteration !!\n",
                         MPI_Rank, P );


//       c. copy data : Pot_Array_Int --> Pot_Array_Out
// ------------------------------------------------------------------------------------------------------------
         for (int k=0; k<GRA_NXT; k++)    {  K = k + POT_GHOST_SIZE + POT_USELESS - GRA_GHOST_SIZE;
         for (int j=0; j<GRA_NXT; j++)    {  J = j + POT_GHOST_SIZE + POT_USELESS - GRA_GHOST_SIZE;
         for (int i=0; i<GRA_NXT; i++)    {  I = i + POT_GHOST_SIZE + POT_USELESS - GRA_GHOST_SIZE;

            Pot_Array_Out[P][k][j][i] = Pot_Array_Int[K][J][I];

         }}}

      } // for (int P=0; P<NPatch; P++)


      delete [] Pot_Array_Int;

   } // OpenMP parallel region

} // FUNCTION : CPU_PoissonSolver_SOR



#endif // #if ( defined GRAVITY  &&  !defined GPU  &&  POT_SCHEME == SOR )
