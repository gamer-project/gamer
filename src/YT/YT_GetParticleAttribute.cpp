#include "GAMER.h"

#ifdef SUPPORT_LIBYT

void YT_GetPID( const long gid, int *level, int *PID );

#ifdef PARTICLE


//-------------------------------------------------------------------------------------------------------
// Function    :  Get_ParticleAttribute
// Description :  Get the particle attribute.
//
// Note        :  1. This function's pointer will be passed into libyt.
//                2. The argument should be declared like this, in order to match the libyt API.
//                3. This function will be called if yt inline-analysis needs the data.
//
// Parameter   :  long  gid     : Grid GID
//                char *attr    : Request particle attribute
//                void *data    : Data pointer, write the results to here.
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Get_ParticleAttribute( const int list_len, const long *list_gid, const char *ptype, const char *attr, yt_array *data_array )
{

   bool ParAttIsInt;

// get attribute index in GAMER
   FieldIdx_t ParAttrFlt_Idx, ParAttrInt_Idx, ParAttr_Idx;

   ParAttrFlt_Idx = GetParticleAttributeFltIndex( attr, CHECK_OFF );
   ParAttrInt_Idx = GetParticleAttributeIntIndex( attr, CHECK_OFF );

   if ( ParAttrFlt_Idx != Idx_Undefined )
   {
      ParAttIsInt = false;
      ParAttr_Idx = ParAttrFlt_Idx;
   }
   else if ( ParAttrInt_Idx != Idx_Undefined )
   {
      ParAttIsInt = true;
      ParAttr_Idx = ParAttrInt_Idx;
   }
   else
   {
      Aux_Error( ERROR_INFO, "Cannot find the target particle attribute \"%s\" !!\n", attr );
   }

// loop through list_gid
   for (int lid=0; lid<list_len; lid++)
   {
//    parse level and pid from gid
      int level, PID0;
      YT_GetPID( list_gid[lid], &level, &PID0 );

//    write data
#     ifdef LIBYT_USE_PATCH_GROUP
      long ParID;
      long data_idx = 0;
      for (int i=0; i<8; i++) // run through 8 patches
      {
         for (int p=0; p<amr->patch[0][level][PID0 + i]->NPar; p++) // run through particle data in one PID
         {
            ParID = amr->patch[0][level][PID0 + i]->ParList[p];
            if ( ParAttIsInt ) ((long_par *) data_array[lid].data_ptr)[data_idx] = amr->Par->AttributeInt[ParAttr_Idx][ParID];
            else               ((real_par *) data_array[lid].data_ptr)[data_idx] = amr->Par->AttributeFlt[ParAttr_Idx][ParID];
            data_idx += 1;
         }
      }
#     else
      long ParID;
      for (int p=0; p<amr->patch[0][level][PID0]->NPar; p++)
      {
         ParID = amr->patch[0][level][PID0]->ParList[p];
         if ( ParAttIsInt ) ((long_par *) data_array[lid].data_ptr)[p] = amr->Par->AttributeInt[ParAttrInt_Idx][ParID];
         else               ((real_par *) data_array[lid].data_ptr)[p] = amr->Par->AttributeFlt[ParAttrFlt_Idx][ParID];
      }
#     endif // #ifdef LIBYT_USE_PATCH_GROUP
   } // for(int lid=0; lid<list_len; lid++)

} // FUNCITON : Get_ParticleAttribute

#endif // #ifdef PARTICLE

#endif // #ifdef SUPPORT_LIBYT
