#include "GAMER.h"

#ifdef SUPPORT_GRACKLE

#ifdef OPENMP
#include <atomic>
#endif


// global variables borrowing from h_Che_Array[]
// --> declared in Init_MemAllocate_Grackle.cpp
extern int Che_NField;
extern int CheIdx_sEint;




//-------------------------------------------------------------------------------------------------------
// Function    :  Grackle_Calculate
// Description :  Calculate various fields using the Grackle API
//
// Note        :  1. Ref: https://grackle.readthedocs.io/en/latest/Interaction.html#calling-the-available-functions
//                2. Currently, the fields have to be calculated in units of patch group
//                3. The output field data type is "real" instead of "real_che"
//                   --> To be consistent with the output usage
//                4. This function can be invoked by
//                   e.g., Output_DumpData_*(), Flag_Real() and Grackle_GetTimeStep_CoolingTime()
//                5. Do not invoke this function by multiple OpenMP threads
//                   (1) Nested parallelization
//                       --> OpenMP parallelization is already implemented in
//                           Grackle_Prepare() and inside Grackle's API functions
//                   (2) Data race
//                       --> Che_FieldData[] is globally shared and thus not thread-safe for concurrent modification
//                           in Grackle_Prepare() or this function
//
// Parameter   :  Out       : Output array; array size = NFieldOut*(NPG*PS2*PS2*PS2)
//                TFields   : Target fields to be calculated:
//                            _GRACKLE_TEMP, _GRACKLE_MU, _GRACKLE_TCOOL
//                            --> Defined in include/Typedef.h
//                lv        : Target refinement level
//                NPG       : Number of patch groups calculated at a time
//                PID0_List : List recording the patch indices with LocalID==0 to be calculated
//
// Return      :  Out[]
//-------------------------------------------------------------------------------------------------------
void Grackle_Calculate( real Out[], const GrackleFieldBIdx_t TFields,
                        const int lv, const int NPG, const int *PID0_List )
{

// 0. Grackle must be enabled
   if ( !GRACKLE_ACTIVATE )
      Aux_Error( ERROR_INFO, "%s() requires GRACKLE_ACTIVATE enabled !!\n", __FUNCTION__ );


// 1. nothing to do if there is no target patch group
   if ( NPG == 0 )   return;

   const int Size1v = NPG * CUBE(PS2); // size of one field


// 2. check the number of output fields
   int NFieldOut = 0;
   const int IdxOut_Undefined = -1;
   const int IdxOut_temp      = ( TFields & _GRACKLE_TEMP  ) ? NFieldOut++ : IdxOut_Undefined;
   const int IdxOut_mu        = ( TFields & _GRACKLE_MU    ) ? NFieldOut++ : IdxOut_Undefined;
   const int IdxOut_tcool     = ( TFields & _GRACKLE_TCOOL ) ? NFieldOut++ : IdxOut_Undefined;

   if ( NFieldOut == 0 )   return;


#  ifdef OPENMP
// 3. check there is no multithreading
   static std::atomic_flag hasOneRunningThread = ATOMIC_FLAG_INIT;
   if ( hasOneRunningThread.test_and_set( std::memory_order_acquire ) )
      Aux_Error( ERROR_INFO, "%s() is invoked concurrently by multiple threads !!\n", __FUNCTION__ );
#  endif


// 4. allocate and prepare the input array for the Grackle solver
   real_che *gr_fields_input = new real_che [ (long)Che_NField*(long)Size1v ];

   Grackle_Prepare( lv, gr_fields_input, NPG, PID0_List );


// 5. allocate the output array for the Grackle solver
   real_che *gr_fields_temperature  = NULL;
   real_che *gr_fields_gamma        = NULL;
   real_che *gr_fields_cooling_time = NULL;

   if ( IdxOut_temp != IdxOut_Undefined  ||  IdxOut_mu != IdxOut_Undefined )
      gr_fields_temperature  = new real_che [Size1v];

   if ( IdxOut_mu != IdxOut_Undefined )
      gr_fields_gamma        = new real_che [Size1v];

   if ( IdxOut_tcool != IdxOut_Undefined )
      gr_fields_cooling_time = new real_che [Size1v];

   typedef real (*vla_out)[Size1v];
   vla_out Out1D = ( vla_out )Out;


// 6. set grid_dimension, grid_start, and grid_end
   const int OptFac = 16;  // optimization factor // TODO: this value requires further optimization for NPG==1
   if ( SQR(PS2)%OptFac != 0 )   Aux_Error( ERROR_INFO, "SQR(PS2) %% OptFac != 0 !!\n" );

   Che_FieldData->grid_dimension[0] = PS2*OptFac;
   Che_FieldData->grid_dimension[1] = 1;
   Che_FieldData->grid_dimension[2] = SQR(PS2)*NPG/OptFac;

   for (int d=0; d<3; d++)
   {
      Che_FieldData->grid_start[d] = 0;
      Che_FieldData->grid_end  [d] = Che_FieldData->grid_dimension[d] - 1;
   }


// 7. invoke Grackle
// --> note that we use the OpenMP implementation in Grackle directly, which applies the parallelization to the first two
//     dimensiones of the input grid
// --> this approach is found to be much more efficient than parallelizing different patches or patch groups here
   if ( IdxOut_temp != IdxOut_Undefined  ||  IdxOut_mu != IdxOut_Undefined )
   {
      if ( calculate_temperature( &Che_Units, Che_FieldData, gr_fields_temperature ) == 0 )
         Aux_Error( ERROR_INFO, "Grackle calculate_temperature() failed !!\n" );
   }

   if ( IdxOut_mu != IdxOut_Undefined )
   {
      if ( calculate_gamma( &Che_Units, Che_FieldData, gr_fields_gamma ) == 0 )
         Aux_Error( ERROR_INFO, "Grackle calculate_gamma() failed !!\n" );
   }

   if ( IdxOut_tcool != IdxOut_Undefined )
   {
      if ( calculate_cooling_time( &Che_Units, Che_FieldData, gr_fields_cooling_time ) == 0 )
         Aux_Error( ERROR_INFO, "Grackle calculate_cooling_time() failed !!\n" );
   }


// 8. fill in the output array
   const real_che *gr_fields_sEint  = gr_fields_input + CheIdx_sEint*Size1v;
   const real_che temperature_units = (real_che) get_temperature_units( &Che_Units );

#  ifdef GAMER_DEBUG
   if ( IdxOut_mu != IdxOut_Undefined  &&  temperature_units <= 0.0 )
      Aux_Error( ERROR_INFO, "Unphysical temperature_units returned by Grackle = %14.7e !!\n", temperature_units );
#  endif // #ifdef GAMER_DEBUG

   for (int i=0; i<Size1v; i++)
   {
//    check
#     ifdef GAMER_DEBUG
      if ( ( IdxOut_temp != IdxOut_Undefined  ||  IdxOut_mu != IdxOut_Undefined )  &&  gr_fields_temperature[i] <= 0.0 )
         Aux_Error( ERROR_INFO, "Unphysical temperature returned by Grackle = %14.7e !!\n", gr_fields_temperature[i] );

      if ( IdxOut_mu != IdxOut_Undefined  &&  gr_fields_gamma[i] <= 1.0 )
         Aux_Error( ERROR_INFO, "Unphysical gamma returned by Grackle = %14.7e !!\n", gr_fields_gamma[i] );

      if ( IdxOut_mu != IdxOut_Undefined  &&  gr_fields_sEint[i] <= 0.0 )
         Aux_Error( ERROR_INFO, "Unphysical sEint = %14.7e !!\n", gr_fields_sEint[i] );

      if ( IdxOut_tcool != IdxOut_Undefined  &&  gr_fields_cooling_time[i] != gr_fields_cooling_time[i] )
         Aux_Error( ERROR_INFO, "Unphysical cooling time returned by Grackle = %14.7e !!\n", gr_fields_cooling_time[i] );
#     endif // #ifdef GAMER_DEBUG

//    temperature
      if ( IdxOut_temp != IdxOut_Undefined )
         Out1D[IdxOut_temp][i]  = (real) gr_fields_temperature[i];

//    mean molecular weight
      if ( IdxOut_mu != IdxOut_Undefined )
         Out1D[IdxOut_mu][i]    = (real) ( gr_fields_temperature[i] / ( gr_fields_sEint[i] * temperature_units * (gr_fields_gamma[i] - (real_che)1.0) ) );

//    cooling time
      if ( IdxOut_tcool != IdxOut_Undefined )
         Out1D[IdxOut_tcool][i] = (real) gr_fields_cooling_time[i];
   }


// 9. free memory
   delete [] gr_fields_input;
   delete [] gr_fields_temperature;
   delete [] gr_fields_gamma;
   delete [] gr_fields_cooling_time;


#  ifdef OPENMP
// 10. clear the flag
   hasOneRunningThread.clear( std::memory_order_release );
#  endif


} // FUNCTION : Grackle_Calculate



#endif // #ifdef SUPPORT_GRACKLE
