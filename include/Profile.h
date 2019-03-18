#ifndef __PROFILE_H__
#define __PROFILE_H__



void Aux_Error( const char *File, const int Line, const char *Func, const char *Format, ... );




//-------------------------------------------------------------------------------------------------------
// Structure   :  Profile_t
// Description :  Data structure for computing a radial profile
//
// Data Member :  NBin        : Total number of radial bins
//                LogBin      : true/false --> log/linear bins
//                LogBinRatio : Ratio of any two adjacent log bin size
//                MaxRadius   : Maximum radius
//                Center      : Target center coordinates
//                Radius      : Radial coordinate at each bin
//                Data        : Profile data at each bin
//                Weight      : Total weighting at each bin
//                NCell       : Number of cells at each bin
//
// Method      :  Profile_t      : Constructor
//               ~Profile_t      : Destructor
//                AllocateMemory : Allocate memory
//                FreeMemory     : Free memory
//-------------------------------------------------------------------------------------------------------
struct Profile_t
{

// data members
// ===================================================================================
   int    NBin;
   bool   LogBin;
   double LogBinRatio;
   double MaxRadius;
   double Center[3];

   double *Radius;
   double *Data;
   double *Weight;
   long   *NCell;


   //===================================================================================
   // Constructor :  Profile_t
   // Description :  Constructor of the structure "Profile_t"
   //
   // Note        :  Initialize the data members
   //
   // Parameter   :  None
   //===================================================================================
   Profile_t()
   {

      NBin   = -1;
      Radius = NULL;
      Data   = NULL;
      Weight = NULL;
      NCell  = NULL;

   } // METHOD : Profile_t



   //===================================================================================
   // Destructor  :  ~Profile_t
   // Description :  Destructor of the structure "Profile_t"
   //
   // Note        :  Free memory
   //===================================================================================
   ~Profile_t()
   {

      FreeMemory();

   } // METHOD : ~Profile_t



   //===================================================================================
   // Method      :  AllocateMemory
   // Description :  Allocate and initialize arrays
   //
   // Note        :  Invoked by Aux_ComputeProfile()
   //
   // Parameter   :  None
   //
   // Return      :  Radius[], Data[], Weight[], NCell[]
   //===================================================================================
   void AllocateMemory()
   {

      if ( NBin < 0 )   Aux_Error( ERROR_INFO, "NBin (%d) < 0 !!\n", NBin );

      Radius = new double [NBin];
      Data   = new double [NBin];
      Weight = new double [NBin];
      NCell  = new long   [NBin];

      for (int b=0; b<NBin; b++)
      {
         Data  [b] = 0.0;
         Weight[b] = 0.0;
         NCell [b] = 0;
      }

   } // METHOD : AllocateMemory



   //===================================================================================
   // Method      :  FreeMemory
   // Description :  Free the memory allocated by AllocateMemory()
   //
   // Note        :  Invoked by ~Profile()
   //
   // Parameter   :  None
   //
   // Return      :  Radius[], Data[], Weight[], NCell[]
   //===================================================================================
   void FreeMemory()
   {

      delete [] Radius;    Radius = NULL;
      delete [] Data;      Data   = NULL;
      delete [] Weight;    Weight = NULL;
      delete [] NCell;     NCell  = NULL;

   } // METHOD : FreeMemory


}; // struct Profile_t



#endif // #ifndef __PROFILE_H__
