#include <cstdio>

//#define FLOAT8_PAR
#define INT8_PAR



#ifdef FLOAT8_PAR
typedef double real_par;
#else
typedef float  real_par;
#endif

#ifdef INT8_PAR
typedef long long_par;
#else
typedef int  long_par;
#endif


int main()
{

   const int  NUM_PARTICLE      = 20000;
   const int  NUM_ATTRIBUTE_FLT = 7;
   const int  NUM_ATTRIBUTE_INT = 1;
   const bool PAR_IC_ID_ATT     = true; // data format of PAR_IC: (true: [id][attribute], false: [attribute][id]; row-major)

   real_par (*ParIC_Flt)[NUM_PARTICLE] = new real_par [NUM_ATTRIBUTE_FLT][NUM_PARTICLE];
   long_par (*ParIC_Int)[NUM_PARTICLE] = new long_par [NUM_ATTRIBUTE_INT][NUM_PARTICLE];

   for (int p=0; p<NUM_PARTICLE; p++)
   {
//    replace the following lines by your particle initial condition
      ParIC_Flt[0][p] = 1.1;   // mass
      ParIC_Flt[1][p] = 2.2;   // position x
      ParIC_Flt[2][p] = 3.3;   // position y
      ParIC_Flt[3][p] = 4.4;   // position z
      ParIC_Flt[4][p] = 5.5;   // velocity x
      ParIC_Flt[5][p] = 6.6;   // velocity y
      ParIC_Flt[6][p] = 7.7;   // velocity z

      ParIC_Int[0][p] = 1;     // type (generic massive)
   }

   FILE *File = fopen( "PAR_IC", "wb" );

   if ( PAR_IC_ID_ATT )
   {
      for (int p=0; p<NUM_PARTICLE; p++)
      {
         for (int v=0; v<NUM_ATTRIBUTE_FLT; v++) fwrite( &ParIC_Flt[v][p], sizeof(real_par), 1, File );
         for (int v=0; v<NUM_ATTRIBUTE_INT; v++) fwrite( &ParIC_Int[v][p], sizeof(long_par), 1, File );
      }
   }
   else
   {
      for (int v=0; v<NUM_ATTRIBUTE_FLT; v++) fwrite( ParIC_Flt[v], sizeof(real_par), NUM_PARTICLE, File );
      for (int v=0; v<NUM_ATTRIBUTE_INT; v++) fwrite( ParIC_Int[v], sizeof(long_par), NUM_PARTICLE, File );
   }

   fclose( File );

   delete [] ParIC_Flt;
   delete [] ParIC_Int;

} // FUNCTION : main
