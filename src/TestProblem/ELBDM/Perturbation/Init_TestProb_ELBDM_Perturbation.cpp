#include "GAMER.h"



// problem-specific global variables
// =======================================================================================
static double   Perturbation_Amplitude;      // maximum amplitude of perturbations
static double   Perturbation_BgAmplitude;    // amplitude of background density
static int      Perturbation_N;              // number of excitations
static int      Perturbation_NDim;           // dimensionality of the problem
static double   Perturbation_Wavelength[3];  // wavelength of the largest perturbation in x, y, z direction

// number of modes in Fourier series
#define FOURIER_MODE_NUMBER         4
// number of combinations of odd and even waves in 3D
// sine = s; cosine = c
// |sss, ssc, scc, ccc, scs, css, csc, ccs| = 8
#define EVEN_ODD_WAVE_COMBINATIONS  8

// array of uniform, random numbers between 0 and 1 used as amplitudes of plane waves
// numbers are kept fixed to ensure reproducibility of test
const real Perturbation_Coefficients[EVEN_ODD_WAVE_COMBINATIONS][FOURIER_MODE_NUMBER][FOURIER_MODE_NUMBER][FOURIER_MODE_NUMBER] = {{{{ 3.70157933e-01,  2.87661055e-01, -1.14617545e-01,-1.30226995e-01}, { 4.83364854e-01, -4.75112608e-01, -7.29319183e-02, 1.38984806e-01}, { 2.16194361e-01,  3.41238342e-01, -2.87929415e-01, 2.25786458e-02}, {-4.83181107e-01, -2.57694103e-01,  2.91195071e-01, 3.34835372e-01}},{{ 3.99972133e-01,  4.51725653e-01,  3.35312718e-01, 2.66909282e-02}, {-4.23703362e-02, -4.11875769e-01, -1.83606956e-01,-1.78758890e-01}, { 4.96922272e-02, -1.50883654e-01, -1.91524518e-01, 4.98346117e-01}, { 3.03750126e-01,  4.06333736e-01,  4.08116864e-01,-2.81812666e-01}},{{ 3.22437194e-01,  2.48851934e-01,  5.39195840e-02, 2.11813153e-01}, { 3.18778971e-01,  3.40759185e-02, -4.80865598e-01,-2.95609488e-01}, { 4.47373859e-01,  4.03135667e-01,  4.16538014e-01,-2.69895849e-01}, { 3.70663193e-01, -1.37924688e-01, -1.42938379e-01,-4.09558646e-01}},{{-3.11079156e-01, -2.76410015e-01, -3.28665676e-01, 4.47260927e-02}, { 1.76360058e-01, -4.66259868e-01, -2.38358332e-01,-1.47708993e-01}, {-4.82777662e-01,  1.12256224e-01, -3.01588320e-01,-4.74030422e-01}, { 1.07085148e-01, -2.34686965e-01, -2.02382251e-01,-3.37103671e-01}}},       {{{-4.08713347e-01, -3.13828583e-01,  2.18916025e-01,-3.32790995e-01}, {-1.00412594e-01, -4.39902472e-01,  3.62668702e-01,-2.70383886e-01}, {-1.03632364e-01, -8.13027821e-03, -2.76968194e-01, 1.77794425e-01}, { 4.84028109e-01,  2.00890594e-01, -4.76499451e-01, 3.66511913e-01}},{{-3.49768123e-01,  2.42690500e-01, -7.02395920e-02, 4.39562475e-03}, { 1.29983783e-01, -2.37280298e-01,  8.77797873e-02,-3.60590530e-01}, {-2.02311492e-01,  8.37045684e-02, -1.04425110e-01, 3.83863934e-01}, { 1.27421036e-02, -4.11458917e-01,  3.61996430e-01, 8.09036036e-02}},{{-1.42854127e-02,  2.27874506e-01,  4.88905822e-01, 3.88393874e-01}, { 1.21405037e-01,  4.25367583e-01,  2.84862906e-01, 3.52336664e-01}, { 2.57797957e-01,  1.64064348e-01,  1.49956765e-01, 4.88134653e-01}, {-9.86264265e-02,  4.99082438e-01,  2.21083792e-01,-3.97632096e-01}},{{-2.63078419e-03, -4.58650858e-01,  1.71625327e-01,-4.04565318e-01}, {-3.05466543e-01, -1.58551752e-01,  4.72060392e-01,-4.70209661e-01}, { 4.56730484e-01, -4.45751177e-01, -3.13578080e-01, 4.92641643e-01}, { 3.03883947e-01, -1.85514920e-01,  3.37478030e-01,-4.07561715e-02}}},       {{{ 2.62097921e-01,  2.97092238e-01, -3.73326843e-01,-3.39537017e-01}, {-4.95683398e-01,  3.21975929e-01,  3.54554653e-01,-1.72553711e-01}, {-3.62993596e-01,  9.64676591e-02,  1.80556501e-01,-4.45890979e-01}, { 3.13441971e-01,  4.34439590e-01,  2.42752851e-01,-1.82196611e-01}},{{ 9.20899806e-03,  8.97535287e-02,  3.19834761e-01, 1.56834140e-01}, { 1.83504262e-01, -3.95596977e-01, -3.79329605e-01,-4.21551606e-01}, { 3.17537613e-01, -4.08168915e-01, -1.14194353e-01, 4.29083816e-01}, { 2.88347641e-01,  1.78310318e-01,  3.27172528e-01,-3.21315615e-01}},{{ 3.60055424e-01,  4.45609527e-01,  4.31798216e-02, 2.70207899e-01}, { 1.84280226e-01,  3.46331783e-01, -7.68265538e-02, 2.59857745e-01}, { 1.60618091e-01, -3.44586908e-01, -1.34722316e-01,-2.08190523e-03}, {-7.01402042e-02,  4.85526042e-02,  4.97900282e-01, 7.61367799e-03}},{{-4.41327900e-01,  3.42778044e-01,  3.42692990e-01, 2.45246921e-01}, {-1.89767451e-01,  2.49228412e-01, -9.12325997e-02, 2.67070264e-02}, {-3.68875850e-01,  2.77455339e-01,  3.63014414e-01,-3.24135262e-01}, { 2.51946115e-01,  4.39481375e-01,  3.56650168e-01, 3.84638290e-01}}},       {{{-4.53928002e-02, -1.69417888e-02, -1.46374118e-01,-2.23106453e-01}, { 6.87297269e-02,  4.51183654e-01,  7.83833733e-02, 1.02898840e-01}, {-4.83877530e-01,  3.69301569e-01, -4.06203833e-01,-3.85813112e-01}, { 3.07595826e-01,  3.38553389e-01, -7.22095211e-02,-3.64892489e-01}},{{ 4.54481059e-02, -3.36558916e-01, -2.00978848e-01, 3.33308470e-01}, { 1.71080805e-01,  9.06790849e-03,  2.29485173e-01, 1.85577715e-01}, {-1.80768843e-01,  2.93846645e-01,  1.31851251e-01,-6.74915994e-02}, {-4.44747160e-02,  1.53734555e-01,  3.79579201e-01, 9.48838933e-02}},{{ 3.12423448e-01, -1.86750768e-01,  2.02667198e-01,-5.76802213e-02}, { 2.00821073e-02, -2.15202092e-01,  3.48253885e-01,-2.01645085e-01}, { 2.82557207e-01, -3.18390960e-01,  4.35763156e-01,-4.10745827e-01}, {-2.24121513e-01, -2.43306191e-01,  4.34910113e-02,-2.01556085e-01}},{{ 4.95275319e-01, -2.08364975e-01, -2.14653283e-02, 2.25898744e-01}, {-2.97036076e-02,  3.61856364e-01, -2.48804357e-01,-4.75634739e-01}, {-2.41975291e-02, -3.75060214e-01,  4.62327114e-01, 4.98309613e-01}, { 4.64745001e-01,  1.70913893e-01, -2.26453867e-01,-1.10852500e-01}}},       {{{ 2.82354296e-01, -2.60326987e-01,  1.09137825e-01, 2.50924407e-01}, { 1.73186938e-01, -2.04184583e-01,  2.19582731e-01,-2.31767254e-01}, { 6.39284787e-02, -9.91275331e-02,  4.56313778e-01,-2.97354470e-01}, { 4.86058379e-01, -4.23742689e-01, -2.29603730e-01,-2.92610904e-01}},{{ 2.71353093e-01,  3.02961538e-02, -3.26585760e-01,-3.28235288e-01}, {-3.79138984e-01,  4.23960923e-01,  4.43697676e-01,-8.72572340e-02}, {-1.62106341e-01, -2.92301229e-01,  5.72712020e-02, 4.85068780e-01}, {-4.38076467e-01, -6.33783327e-02,  3.11416607e-01, 3.26701174e-01}},{{ 3.48273688e-01,  1.25223847e-01, -3.21015071e-01,-2.31706245e-01}, {-1.59149642e-01, -4.70362740e-01,  6.03381845e-02, 3.68068173e-01}, { 8.75896960e-02, -1.69765926e-01,  6.50150737e-02, 8.81910599e-02}, {-5.88508645e-02,  2.43264706e-01, -3.81893714e-01,-1.22887531e-01}},{{ 1.39862438e-01,  1.97928033e-01, -2.92122745e-01,-4.88128859e-01}, { 1.98964559e-01,  3.90963169e-01, -4.20036820e-01, 5.42304392e-02}, {-2.50761287e-02,  3.01265782e-01,  4.54165847e-01,-2.17137378e-01}, { 1.03521342e-01,  4.59106162e-01, -4.79573327e-01, 4.96686363e-01}}},       {{{ 2.11983181e-01,  8.93027270e-02,  3.57318532e-01, 3.25551087e-01}, { 4.00435477e-01,  9.83212895e-02,  4.76031060e-01,-4.50948616e-01}, {-3.63721516e-01, -4.48825789e-01,  3.30491120e-02, 1.78796882e-01}, { 4.03174572e-01,  4.31331858e-01, -3.49696233e-01, 2.36679207e-01}},{{ 3.67147385e-01, -1.81236388e-01,  2.74202109e-01, 1.15034060e-01}, {-4.43555318e-01,  1.80084099e-01, -2.66868662e-01, 8.24194506e-02}, {-1.84774575e-01, -2.13175273e-01, -1.63463928e-02, 4.89946802e-01}, { 8.89487215e-03,  2.37742434e-01, -4.23641968e-01, 3.32801898e-01}},{{-2.52801356e-01, -2.55633512e-01,  4.30348206e-01, 2.35484929e-01}, { 1.74822782e-01, -1.30220469e-01, -3.81384359e-01,-6.83468781e-02}, {-4.32505606e-03, -7.89800007e-02,  3.85830341e-01,-2.95375805e-01}, {-3.66000771e-01, -3.92086527e-01, -4.31152595e-01, 1.35240517e-01}},{{ 4.87098219e-01, -3.18127523e-01,  7.67794535e-02,-5.72398610e-02}, {-4.46582225e-01,  2.38632127e-01,  2.29158497e-01, 2.54726381e-01}, { 4.61098653e-01, -2.40051714e-01,  1.13974696e-01, 2.03922752e-01}, {-1.60903152e-01, -3.53346814e-01,  1.52593709e-01,-4.16927445e-01}}},       {{{ 4.99163616e-01, -4.81549572e-01, -3.53001458e-01, 3.09968565e-01}, {-3.62482142e-01, -1.51464726e-02,  9.94907721e-03, 4.98152728e-01}, {-1.48042408e-01,  1.34836973e-01,  1.24917919e-01, 1.19056334e-01}, { 3.05008661e-01,  2.61492516e-01, -2.55336855e-01, 1.10654974e-01}},{{ 3.04598880e-01, -1.17095954e-01,  6.66073447e-02,-5.77098993e-02}, {-2.12113708e-01,  1.04647314e-01, -3.30773266e-01, 2.98431503e-01}, { 4.42795331e-01,  1.75432930e-01, -4.56080097e-01, 1.79593750e-01}, {-4.50738968e-01, -3.10565110e-01, -1.34886859e-04, 1.98093939e-03}},{{ 2.40201438e-01, -1.14853257e-01,  3.92644134e-01, 6.80367991e-02}, {-3.52991709e-01,  4.12317138e-01, -4.24015245e-01, 1.57561217e-01}, { 3.33282598e-01,  4.53179149e-01,  4.69813280e-01,-1.40855321e-01}, {-2.55796165e-01, -4.58509046e-01, -3.08439026e-01, 3.55076958e-01}},{{ 1.14146174e-01,  1.21676640e-01,  4.23661797e-02,-1.23919474e-01}, { 8.27115567e-02,  2.10506290e-01, -3.85005862e-01, 3.39318514e-01}, { 1.46903203e-01, -2.15374053e-01, -1.71556316e-01,-4.55641350e-01}, {-3.55815584e-01,  4.64417590e-01, -5.71050339e-02, 3.08085723e-01}}},       {{{-1.57677351e-01,  9.22630040e-02,  3.04553538e-01,-4.39615238e-01}, { 2.32706436e-01,  1.17098293e-01, -2.88543970e-01, 1.70200705e-01}, {-4.96956557e-01, -3.36441256e-01,  9.20915790e-02,-1.87133783e-01}, { 1.71003256e-01,  4.83561674e-02, -1.82343048e-01,-4.87830473e-01}},{{-6.60428634e-02, -3.53534901e-01,  4.50634031e-01, 4.33414701e-01}, {-2.49377988e-01, -7.57523106e-02,  4.85388323e-01, 3.25545359e-01}, { 1.06891548e-01,  3.74882682e-01,  4.04337197e-01, 4.69607825e-02}, {-1.63350312e-02, -4.06013837e-01,  1.23478918e-01, 4.10723417e-01}},{{ 1.72497094e-02, -9.92720126e-02, -1.51539312e-01,-3.35211738e-01}, { 7.68460268e-02,  1.67474396e-01, -1.57993138e-02,-3.83150405e-01}, { 3.09872075e-01,  3.81183746e-01, -4.60159700e-01,-3.46941557e-01}, {-2.00808788e-01, -2.57349550e-01,  1.60312747e-01, 4.04826745e-01}},{{-2.86764579e-01,  3.34851608e-01,  2.10534524e-02,-2.57263471e-01}, { 1.93732324e-01, -4.15945686e-01,  4.80815766e-01,-5.61149598e-02}, { 4.03009815e-02, -1.87718629e-01,  2.38502533e-01, 2.23772782e-01}, { 3.56198340e-01, -5.08654700e-02, -5.24900482e-02,-9.56326772e-02}}}};
// =======================================================================================




//-------------------------------------------------------------------------------------------------------
// Function    :  Validate
// Description :  Validate the compilation flags and runtime parameters for this test problem
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Validate()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ...\n", TESTPROB_ID );


// errors
#  if ( MODEL != ELBDM )
   Aux_Error( ERROR_INFO, "MODEL != ELBDM !!\n" );
#  endif

#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   if ( !ELBDM_MATCH_PHASE )
      Aux_Message(stderr, "WARNING: If ELBDM_MATCH_PHASE is not enabled, this test case will fail !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == ELBDM )
//-------------------------------------------------------------------------------------------------------
// Function    :  SetParameter
// Description :  Load and set the problem-specific runtime parameters
//
// Note        :  1. Filename is set to "Input__TestProb" by default
//                2. Major tasks in this function:
//                   (1) load the problem-specific runtime parameters
//                   (2) set the problem-specific derived parameters
//                   (3) reset other general-purpose parameters if necessary
//                   (4) make a note of the problem-specific parameters
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void SetParameter()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );


// (1) load the problem-specific runtime parameters
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

// (1-1) add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., Useless_bool, Eps_double, NoMin_int, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",           &VARIABLE,                   DEFAULT,       MIN,              MAX                  );
// ********************************************************************************************************************************
   ReadPara->Add( "Perturbation_N",            &Perturbation_N,            -1,             1,                FOURIER_MODE_NUMBER  );
   ReadPara->Add( "Perturbation_Amp",          &Perturbation_Amplitude,     NoDef_double,  NoMin_double,     NoMax_double         );
   ReadPara->Add( "Perturbation_BgAmp",        &Perturbation_BgAmplitude,   NoDef_double,  NoMin_double,     NoMax_double         );
   ReadPara->Add( "Perturbation_NDim",         &Perturbation_NDim,          3,             1,                3                    );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values
   for (int i=0; i<3; ++i)
      Perturbation_Wavelength[i] = 2.0 * M_PI / amr->BoxSize[i];

   switch ( Perturbation_NDim ) {
      case 1:
      Perturbation_Wavelength[1] = 0;
      Perturbation_Wavelength[2] = 0;
      break;

      case 2:
      Perturbation_Wavelength[2] = 0;
      break;

      default:
      break;
   }

// (1-3) check the runtime parameters
   if ( Perturbation_Amplitude == NoDef_double )
      Aux_Error( ERROR_INFO, "Runtime parameter \"Perturbation_Amplitude\" is not set !!\n" );

   if ( Perturbation_BgAmplitude == NoDef_double )
      Aux_Error( ERROR_INFO, "Runtime parameter \"Perturbation_Amplitude\" is not set !!\n" );


// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 2.0;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_RESET_PARA( END_STEP, FORMAT_LONG, "" );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_RESET_PARA( END_T, FORMAT_REAL, "" );
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "======================================================================================\n" );
      Aux_Message( stdout, "  test problem ID                   = %d\n",     TESTPROB_ID              );
      Aux_Message( stdout, "  amplitude of background amplitude = %13.7e\n", Perturbation_BgAmplitude );
      Aux_Message( stdout, "  amplitude of perturbation         = %13.7e\n", Perturbation_Amplitude   );
      Aux_Message( stdout, "  number of Fourier modes           = %d\n",     Perturbation_N           );
      Aux_Message( stdout, "  number of dimensions              = %d\n",     Perturbation_NDim        );
      Aux_Message( stdout, "======================================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  1. This function may also be used to estimate the numerical errors when OPT__OUTPUT_USER is enabled
//                   --> In this case, it should provide the analytical solution at the given "Time"
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
//
// Parameter   :  fluid    : Fluid field to be initialized
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] )
{

   double kx, ky, kz, d1, d2, d3, d4, d5, d6, dtotal, omega;

   double Re = Perturbation_BgAmplitude;
   double Im = 0.0;

   for (int n=0; n<Perturbation_N; ++n) {
   for (int m=0; m<Perturbation_N; ++m) {
   for (int l=0; l<Perturbation_N; ++l) {
      kx = Perturbation_Wavelength[0]*( n + 1 );
      ky = Perturbation_Wavelength[1]*( m + 1 );
      kz = Perturbation_Wavelength[2]*( l + 1 );
      d1 = cos( kx*x );
      d2 = cos( ky*y );
      d3 = cos( kz*z );
      d4 = sin( kx*x );
      d5 = sin( ky*y );
      d6 = sin( kz*z );

      dtotal = + Perturbation_Coefficients[0][n][m][l] * d1 * d2 * d3
               + Perturbation_Coefficients[1][n][m][l] * d4 * d2 * d3
               + Perturbation_Coefficients[2][n][m][l] * d1 * d5 * d3
               + Perturbation_Coefficients[2][n][m][l] * d1 * d5 * d3
               + Perturbation_Coefficients[3][n][m][l] * d1 * d2 * d6
               + Perturbation_Coefficients[4][n][m][l] * d4 * d5 * d3
               + Perturbation_Coefficients[5][n][m][l] * d1 * d5 * d6
               + Perturbation_Coefficients[6][n][m][l] * d4 * d2 * d6
               + Perturbation_Coefficients[7][n][m][l] * d4 * d5 * d6;
      omega = 0.5/ELBDM_ETA*( SQR(kx) + SQR(ky) + SQR(kz) );
      Re += Perturbation_Amplitude*dtotal*cos( -Time*omega );
      Im += Perturbation_Amplitude*dtotal*sin( -Time*omega );
   }}}

   fluid[DENS] = SQR( Re ) + SQR( Im );

#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   if ( amr->use_wave_flag[lv] ) {
#  endif
   fluid[REAL] = Re;
   fluid[IMAG] = Im;
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   } else {
   fluid[PHAS] = SATAN2( Im, Re );
   fluid[STUB] = 0.0;
   }
#  endif

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  OutputError
// Description :  Output the L1 error
//
// Note        :  1. Invoke Output_L1Error()
//                2. Use SetGridIC() to provide the analytical solution at any given time
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void OutputError()
{

   const char Prefix[100]     = "Perturbation";
   const OptOutputPart_t Part = OUTPUT_X;

   Output_L1Error( SetGridIC, NULL, Prefix, Part, 0.0, 0.0, 0.0 );

} // FUNCTION : OutputError
#endif // #if ( MODEL == ELBDM )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_ELBDM_Soliton
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_ELBDM_Perturbation()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


#  if ( MODEL == ELBDM )
// validate the compilation flags and runtime parameters
   Validate();

// set the problem-specific runtime parameters
   SetParameter();

   Init_Function_User_Ptr = SetGridIC;
   Output_User_Ptr        = OutputError;
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_Perturbation
