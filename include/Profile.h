#ifndef __PROFILE_H__
#define __PROFILE_H__



void Aux_Error( const char *File, const int Line, const char *Func, const char *Format, ... );




//-------------------------------------------------------------------------------------------------------
// Structure   :  Profile_t
// Description :  Data structure for computing a radial profile
//
// Data Member :  NBin        : Total number of radial bins
//                LogBin      : true/false --> log/linear bins
//                LogBinRatio : Ratio of adjacent log bins
//                MaxRadius   : Maximum radius
//                Center      : Target center coordinates
//                Allocated   : Whether or not member arrays such as Radius[] have been allocated
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
   bool   Allocated;

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

      NBin      = -1;
      Allocated = false;
      Radius    = NULL;
      Data      = NULL;
      Weight    = NULL;
      NCell     = NULL;

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
   // Description :  Allocate member arrays
   //
   // Note        :  1. Invoked by Aux_ComputeProfile()
   //                2. No data initialization is done here
   //
   // Parameter   :  None
   //
   // Return      :  Radius[], Data[], Weight[], NCell[], Allocated
   //===================================================================================
   void AllocateMemory()
   {

      if ( NBin < 0 )   Aux_Error( ERROR_INFO, "NBin (%d) < 0 !!\n", NBin );

//    free the previously allocated memory in case NBin has changed
//    --> allows the same Profile_t object to be reused without the need to manually free memory first
      if ( Allocated )  FreeMemory();

      Radius = new double [NBin];
      Data   = new double [NBin];
      Weight = new double [NBin];
      NCell  = new long   [NBin];

      Allocated = true;

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

      Allocated = false;

   } // METHOD : FreeMemory


}; // struct Profile_t



#endif // #ifndef __PROFILE_H__
