#ifndef __EXTREMA_H__
#define __EXTREMA_H__




//-------------------------------------------------------------------------------------------------------
// Structure   :  Extrema_t
// Description :  Data structure for finding the extrema
//
// Data Member :  In_Field    : Target field (e.g., _DENS or BIDX(DENS))
//                In_Radius   : Target radius
//                In_Center   : Target center
//
//                Out_Extrema : Output extreme valu3
//                Out_Coord   : Coordinates of the extrema
//                Out_Rank    : MPI rank of the extrema
//                Out_Level   : AMR level of the extrema
//                Out_PID     : Local Patch index of the extrema
//                Out_Cell    : Cell indices within a patch of the extrema
//
// Method      :  Extrema_t : Constructor
//               ~Extrema_t : Destructor
//-------------------------------------------------------------------------------------------------------
struct Extrema_t
{

// data members
// ===================================================================================
   int    In_Field;
   double In_Radius;
   double In_Center[3];

   real   Out_Extrema;
   double Out_Coord[3];
   int    Out_Rank;
   int    Out_Level;
   int    Out_PID;
   int    Out_Cell[3];


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

      In_Field     = -1;
      In_Radius    = NULL_REAL;
      In_Center[0] = NULL_REAL;
      In_Center[1] = NULL_REAL;
      In_Center[2] = NULL_REAL;

      Out_Extrema  = NULL_REAL;
      Out_Coord[0] = NULL_REAL;
      Out_Coord[1] = NULL_REAL;
      Out_Coord[2] = NULL_REAL;
      Out_Rank     = -1;
      Out_Level    = -1;
      Out_PID      = -1;
      Out_Cell[0]  = -1;
      Out_Cell[1]  = -1;
      Out_Cell[2]  = -1;

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


}; // struct Extrema_t



#endif // #ifndef __EXTREMA_H__
