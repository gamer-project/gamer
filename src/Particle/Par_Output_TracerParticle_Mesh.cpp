#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Output_TracerParticle_Mesh
// Description :  Prepare and map mesh quantities onto tracer particles at their positions
//                for storage in HDF5 snapshots
//
// Note        :  1. Invoked by Output_DumpData_Total_HDF5()
//                2. The quantities for non-tracer particles are set to __FLT_MAX__
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Par_Output_TracerParticle_Mesh()
{

   const bool       IntPhase_No       = false;
   const bool       DE_Consistency_No = false;
   const real       MinDens_No        = -1.0;
   const real       MinPres_No        = -1.0;
   const real       MinTemp_No        = -1.0;
   const real       MinEntr_No        = -1.0;
   const bool       UseTracers_Yes    = true;
   const bool       TracerVelCorr_No  = false;
#  ifdef GRAVITY
   const OptPotBC_t BC_Pot            = OPT__BC_POT;
#  else
   const OptPotBC_t BC_Pot            = BC_POT_NONE;
#  endif
   const int        ParGhost          = amr->Par->GhostSizeTracer;
   const int        VarSize           = PS1 + 2*ParGhost;
   const real_par  *ParPos[3]         = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
   const long_par  *ParType           = amr->Par->Type;
   const long       ParListSize       = amr->Par->ParListSize;
   const int        ParAttrNum        = amr->Par->Mesh_Attr_Num;
   const long      *ParAttrIdx        = amr->Par->Mesh_Attr_Idx;

   long AllVar_PreparePatch = ( _TOTAL | _DERIVED );

#  ifdef MHD
   AllVar_PreparePatch |= _MAG;
#  endif
#  ifdef GRAVITY
   AllVar_PreparePatch |= _POTE;
#  endif


// nothing to do if no fields have been specified
   if ( !ParAttrNum )   return;


// allocate enough memory for Mesh_Attr
   for (int v=0; v<ParAttrNum; v++)
      amr->Par->Mesh_Attr[v] = (real_par*) realloc( amr->Par->Mesh_Attr[v], ParListSize*sizeof(real_par) );


// initialize the value to __FLT_MAX__ for all types of particles
   for (int  v=0; v<ParAttrNum;  v++)
   for (long p=0; p<ParListSize; p++)
      amr->Par->Mesh_Attr[v][p] = (real_par) __FLT_MAX__;


// OpenMP parallel region
#  pragma omp parallel
   {

   for (int lv=0; lv<NLEVEL; lv++)
   {
//    get the maximum number of particles in a single patch
//    --> must use "NPar" instead of "NParType[(int)PTYPE_TRACER]" since currently
//        both Var_Temp[] and InterpParPos[] still allocate memory for non-tracer particles
      int NParMax = 0;
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)   NParMax = MAX( NParMax, amr->patch[0][lv][PID]->NPar );

//    nothing to do if there is no particle
      if ( NParMax <= 0 )  continue;

      const double TimeNew = amr->FluSgTime[lv][ amr->FluSg[lv] ];
      const double dh      = amr->dh[lv];
      const double _dh     = 1.0/dh;

      bool  GotYou;
      long  ParID;

      real     *Var      = new real     [ 8*CUBE(VarSize) ];   // 8: number of patches per patch group
      real_par *Var_Temp = new real_par [NParMax];

      real_par **InterpParPos = NULL;
      Aux_AllocateArray2D( InterpParPos, 3, NParMax );


//    loop over all **real** patch groups
#     pragma omp for schedule( PAR_OMP_SCHED, PAR_OMP_SCHED_CHUNK )
      for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
      {
//       1. find the patch groups with target tracer particles
//       --> use patch group as the calculation unit since Prepare_PatchData() only works with patch group
//       --> disadvantage: some patches may not have tracer particles ... (they will be skipped later)
         GotYou = false;

         for (int PID=PID0; PID<PID0+8; PID++)
         {
            if ( amr->patch[0][lv][PID]->NParType[(int)PTYPE_TRACER] > 0 )
            {
               GotYou = true;
               break;
            }
         }

//       nothing to do if there are no target tracer particles in the target patch group
         if ( !GotYou )    continue;


         for (int v=0; v<ParAttrNum; v++)
         {
            const long TVar = ParAttrIdx[v];

//          2. prepare the mesh data for the patch group with particles (need NSIDE_26 for ParGhost>0)
            if ( TVar & AllVar_PreparePatch )
               Prepare_PatchData( lv, TimeNew, Var, NULL, ParGhost, 1, &PID0, TVar, _NONE,
                                  OPT__FLU_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_26, IntPhase_No,
                                  OPT__BC_FLU, BC_Pot, MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );

            else
            {
//###REVISE: support derived fields defined in Flu_DerivedField_BuiltIn() and Flu_DerivedField_User()
               Aux_Error( ERROR_INFO, "unsupported target field (%ld) !!\n", TVar );
            }


//          3. map quantities from mesh onto the tracer particles
            for (int PID=PID0, t=0; PID<PID0+8; PID++, t++)
            {
//             skip patches with no tracer particles
               if ( amr->patch[0][lv][PID]->NParType[(int)PTYPE_TRACER] == 0 )   continue;

               double EdgeL[3], EdgeR[3];

               for (int d=0; d<3; d++)
               {
                  EdgeL[d] = amr->patch[0][lv][PID]->EdgeL[d] - dh*ParGhost;
                  EdgeR[d] = amr->patch[0][lv][PID]->EdgeR[d] + dh*ParGhost;
               }

//             retrieve the particle position
               for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
               {
                  ParID = amr->patch[0][lv][PID]->ParList[p];

                  if ( ParType[ParID] != PTYPE_TRACER )   continue;

                  for (int d=0; d<3; d++)   InterpParPos[d][p] = ParPos[d][ParID];
               }

               Par_MapMesh2Particles( EdgeL, EdgeR, _dh, VarSize, Var+t*CUBE(VarSize),
                                      amr->patch[0][lv][PID]->NPar, InterpParPos, ParType,
                                      amr->patch[0][lv][PID]->ParList, UseTracers_Yes, Var_Temp,
                                      TracerVelCorr_No );

//             4. update particles
               for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
               {
                  ParID = amr->patch[0][lv][PID]->ParList[p];

                  if ( ParType[ParID] != PTYPE_TRACER )   continue;

                  amr->Par->Mesh_Attr[v][ParID] = Var_Temp[p];
               }
            } // for (int PID=PID0, t=0; PID<PID0+8; PID++, t++)
         } // for (int v=0; v<ParAttrNum; v++)
      } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)


//    5. free memory
      delete [] Var;
      delete [] Var_Temp;

      Aux_DeallocateArray2D( InterpParPos );
   } // for (int lv=0; lv<NLEVEL; lv++)

   } // end of OpenMP parallel region

} // FUNCTION : Par_Output_TracerParticle_Mesh



#endif // #ifdef PARTICLE
