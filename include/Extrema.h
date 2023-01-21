#ifndef __EXTREMA_H__
#define __EXTREMA_H__



#include "Macro.h"




//-------------------------------------------------------------------------------------------------------
// Structure   :  Extrema_t
// Description :  Data structure for finding the extrema
//
// Data Member :  Field  : Target field (e.g., _DENS or BIDX(DENS))
//                Radius : Target radius
//                Center : Target center
//
//                Value : Output extreme value
//                Coord : Coordinates of the extrema
//                Rank  : MPI rank of the extrema
//                Level : AMR level of the extrema
//                PID   : Local Patch index of the extrema
//                Cell  : Cell indices within a patch of the extrema
//
// Method      :  Extrema_t     : Constructor
//               ~Extrema_t     : Destructor
//                CreateMPIType : Create a MPI derived datatype for struct Extrema_t
//-------------------------------------------------------------------------------------------------------
struct Extrema_t
{

// data members
// --> must be consistent with CreateMPIType()
// ===================================================================================
// input parameters
   long   Field;
   double Radius;
   double Center[3];

// output parameters
   real   Value;
   double Coord[3];
   int    Rank;
   int    Level;
   int    PID;
   int    Cell[3];


   //===================================================================================
   // Constructor :  Extrema_t
   // Description :  Constructor of the structure "Extrema_t"
   //
   // Note        :  Initialize the data members
   //
   // Parameter   :  None
   //===================================================================================
   Extrema_t()
   {

      Field     = _NONE;
      Radius    = NULL_REAL;
      Center[0] = NULL_REAL;
      Center[1] = NULL_REAL;
      Center[2] = NULL_REAL;

      Value     = NULL_REAL;
      Coord[0]  = NULL_REAL;
      Coord[1]  = NULL_REAL;
      Coord[2]  = NULL_REAL;
      Rank      = -1;
      Level     = -1;
      PID       = -1;
      Cell[0]   = -1;
      Cell[1]   = -1;
      Cell[2]   = -1;

   } // METHOD : Extrema_t



   //===================================================================================
   // Destructor  :  ~Extrema_t
   // Description :  Destructor of the structure "Extrema_t"
   //
   // Note        :  Free memory
   //===================================================================================
   ~Extrema_t()
   {
   } // METHOD : ~Extrema_t



#  ifndef SERIAL
   //===================================================================================
   // Method      :  CreateMPIType
   // Description :  Create a MPI derived datatype for struct Extrema_t
   //
   // Note        :  1. Invoked by Aux_FindExtrema()
   //                2. The returned MPI datatype must be freed manually by calling MPI_Type_free()
   //
   // Parameter   :  MPI_Extrema_t : MPI derived datatype to be returned
   //
   // Return      :  MPI_Extrema_t
   //===================================================================================
   void CreateMPIType( MPI_Datatype *MPI_Extrema_t )
   {

      const int          NBlk         = 9;
      const int          Length[NBlk] = { 1, 1, 3, 1, 3, 1, 1, 1, 3 };
      const MPI_Datatype Type  [NBlk] = { MPI_LONG, MPI_DOUBLE, MPI_DOUBLE, MPI_GAMER_REAL, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_INT };

      MPI_Aint Disp[NBlk];

      Disp[0] = offsetof( Extrema_t, Field  );
      Disp[1] = offsetof( Extrema_t, Radius );
      Disp[2] = offsetof( Extrema_t, Center );
      Disp[3] = offsetof( Extrema_t, Value  );
      Disp[4] = offsetof( Extrema_t, Coord  );
      Disp[5] = offsetof( Extrema_t, Rank   );
      Disp[6] = offsetof( Extrema_t, Level  );
      Disp[7] = offsetof( Extrema_t, PID    );
      Disp[8] = offsetof( Extrema_t, Cell   );

      MPI_Type_create_struct( NBlk, Length, Disp, Type, MPI_Extrema_t );
      MPI_Type_commit( MPI_Extrema_t );

   } // METHOD : CreateMPIType
#  endif // #ifndef SERIAL


}; // struct Extrema_t



#endif // #ifndef __EXTREMA_H__
