#include "GAMER.h"

#ifdef PARTICLE

static long Par_Get_MeshIndex( const char *InputLabel );




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Get_MeshIndex
// Description :  Retrieve the bitwise index corresponding to the input field label
//
// Note        :  1. Invoked by Par_Init_Attribute_Mesh()
//                2. Similar to GetFieldIndex(), but supports derived fields
//                3. Return Idx_Undefined if the target field cannot be found
//
// Parameter   :  InputLabel : Target field label
//
// Return      :  Success: bitwise index of the target field
//                Failed : Idx_Undefined
//-------------------------------------------------------------------------------------------------------
long Par_Get_MeshIndex( const char *InputLabel )
{

// predefined active and passive fields
   for (int v=0; v<NCOMP_TOTAL; v++)
      if (  strcmp( FieldLabel[v], InputLabel ) == 0  )
         return BIDX(v);

// derived fields
#  if ( MODEL == HYDRO )
   if (  strcmp( "VelX", InputLabel ) == 0  )   return _VELX;
   if (  strcmp( "VelY", InputLabel ) == 0  )   return _VELY;
   if (  strcmp( "VelZ", InputLabel ) == 0  )   return _VELZ;
   if (  strcmp( "Pres", InputLabel ) == 0  )   return _PRES;
   if (  strcmp( "Temp", InputLabel ) == 0  )   return _TEMP;
   if (  strcmp( "Entr", InputLabel ) == 0  )   return _ENTR;
   if (  strcmp( "Eint", InputLabel ) == 0  )   return _EINT;
#  ifdef MHD
   if (  strcmp( "MagX", InputLabel ) == 0  )   return _MAGX_CC;
   if (  strcmp( "MagY", InputLabel ) == 0  )   return _MAGY_CC;
   if (  strcmp( "MagZ", InputLabel ) == 0  )   return _MAGZ_CC;
   if (  strcmp( "MagE", InputLabel ) == 0  )   return _MAGE_CC;
#  endif
#  ifdef GRAVITY
   if (  strcmp( "Pote", InputLabel ) == 0  )   return _POTE;
#  endif
//###REVISE: support derived fields defined in Flu_DerivedField_BuiltIn() and Flu_DerivedField_User()
#  endif // #if ( MODEL == HYDRO )


   return Idx_Undefined;

} // FUNCTION : Par_Get_MeshIndex



//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_Attribute_Mesh
// Description :  Parse the fields specified in the file Input__Par_Mesh and initialize
//                associated variables
//
// Note        :  1. Invoked by Init_GAMER()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Par_Init_Attribute_Mesh()
{

   char FileName[] = "Input__Par_Mesh";

// retrieve the number of mesh quantities specified in Input__Par_Mesh
   const int Mesh_NAttr = Aux_CountRow( FileName );

   amr->Par->Mesh_Attr_Num = Mesh_NAttr;


// allocate the memory
   amr->Par->Mesh_Attr       = (real_par**) malloc( Mesh_NAttr*sizeof(real_par*) );
   amr->Par->Mesh_Attr_Label = (char    **) malloc( Mesh_NAttr*sizeof(char    *) );
   amr->Par->Mesh_Attr_Idx   = (long    * ) malloc( Mesh_NAttr*sizeof(long     ) );

   for (int v=0; v<Mesh_NAttr; v++)
   {
      amr->Par->Mesh_Attr      [v] = (real_par*) malloc( amr->Par->ParListSize*sizeof(real_par) );
      amr->Par->Mesh_Attr_Label[v] = (char    *) malloc(            MAX_STRING*sizeof(char    ) );
   }


// load data in Input__Par_Mesh and convert the field name to field index
   char FirstItem[MAX_STRING];
   int  NRow=0, NItem;

   char *Line = new char [MAX_STRING];
   FILE *File = fopen( FileName, "r" );

   while ( fgets(Line, MAX_STRING, File) != NULL )
   {
      NItem = sscanf( Line, "%s", FirstItem );

//    skip empty lines and lines starting with #
      if ( NItem > 0  &&  FirstItem[0] != '#' )
      {
//       initialize the bitwise indices and labels of mesh quantities to be mapped from
         const long FieldIdx = Par_Get_MeshIndex( FirstItem );

         if ( FieldIdx == Idx_Undefined )
            Aux_Error( ERROR_INFO, "unknown input field label (%s) in %s !!\n", FirstItem, FileName );

         amr->Par->Mesh_Attr_Idx[NRow] = FieldIdx;

         sprintf( amr->Par->Mesh_Attr_Label[NRow], "Mesh%s", FirstItem );

         NRow++;
      }
   }


   fclose( File );
   delete [] Line;

} // FUNCTION : Par_Init_Attribute_Mesh



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

   const bool      IntPhase_No       = false;
   const bool      DE_Consistency_No = false;
   const real      MinDens_No        = -1.0;
   const real      MinPres_No        = -1.0;
   const real      MinTemp_No        = -1.0;
   const real      MinEntr_No        = -1.0;
   const bool      UseTracers_Yes    = true;
   const bool      TracerVelCorr_No  = false;
   const int       ParGhost          = amr->Par->GhostSizeTracer;
   const int       VarSize           = PS1 + 2*ParGhost;
   const real_par *ParPos[3]         = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
   const real_par *ParType           = amr->Par->Type;
   const long      ParListSize       = amr->Par->ParListSize;
   const int       ParAttrNum        = amr->Par->Mesh_Attr_Num;
   const long     *ParAttrIdx        = amr->Par->Mesh_Attr_Idx;

   long AllVar_PreparePatch = ( _TOTAL | _DERIVED );

#  ifdef MHD
   AllVar_PreparePatch |= _MAG;
#  endif
#  ifdef GRAVITY
   AllVar_PreparePatch |= _POTE;
#  endif


// allocate enough memory for Mesh_Attr
   for (int v=0; v<ParAttrNum; v++)
      amr->Par->Mesh_Attr[v] = (real_par*) realloc( amr->Par->Mesh_Attr[v], ParListSize*sizeof(real_par) );


// initialize the value to __FLT_MAX__ for all types of particles
   for (int v=0; v<ParAttrNum;  v++)
   for (int p=0; p<ParListSize; p++)
      amr->Par->Mesh_Attr[v][p] = (real_par) __FLT_MAX__;


// map quantities from mesh onto the tracer particles
   for (int lv=0; lv<NLEVEL; lv++)
   {
//    get the maximum number of particles in a single patch
//    --> must use "NPar" instead of "NParType[(int)PTYPE_TRACER]" since currently
//        both Vel_Temp[] and InterpParPos[] still allocate memory for non-tracer particles
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
      for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
      {
//       1. find the patch groups with target tracer particles
//       --> use patch group as the calculation unit since Prepare_PatchData() only works with patch group
//       --> disadvantage: some patches may not have tracer particles ... (they will be skipped later)
         GotYou = false;

         for (int PID=PID0; PID<PID0+8; PID++)
         {
            if ( amr->patch[0][lv][PID]->NParType[(int)PTYPE_TRACER] > 0 )    GotYou = true;

            if ( GotYou )  break;
         }

//       nothing to do if there are no target tracer particles in the target patch group
         if ( !GotYou )    continue;


         for (int v=0; v<ParAttrNum; v++)
         {
            const long TVar = ParAttrIdx[v];

//          2. prepare the mesh data for the patch group with particles (need NSIDE_26 for ParGhost>0 )
            if ( TVar & AllVar_PreparePatch )
               Prepare_PatchData( lv, TimeNew, Var, NULL, ParGhost, 1, &PID0, TVar, _NONE,
                                  OPT__FLU_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_26, IntPhase_No,
                                  OPT__BC_FLU, BC_POT_NONE, MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );

            else
            {
//###REVISE: support derived fields defined in Flu_DerivedField_BuiltIn() and Flu_DerivedField_User()
               Aux_Error( ERROR_INFO, "unsupported target field (%ld) !!\n", TVar );
            }


//          3. map quantities from mesh onto the tracer particles
            for (int PID=PID0, P=0; PID<PID0+8; PID++, P++)
            {
//             skip patches with no tracer particles
               if ( amr->patch[0][lv][PID]->NParType[(int)PTYPE_TRACER] == 0 )   continue;

               double EdgeL[3], EdgeR[3];

               for (int d=0; d<3; d++) {
                  EdgeL[d] = amr->patch[0][lv][PID]->EdgeL[d] - dh*ParGhost;
                  EdgeR[d] = amr->patch[0][lv][PID]->EdgeR[d] + dh*ParGhost;
               }

//             retrieve the particle position
               for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
               {
                  ParID = amr->patch[0][lv][PID]->ParList[p];

                  if ( ParType[ParID] != PTYPE_TRACER )
                     continue;

                  for (int d=0; d<3; d++)   InterpParPos[d][p] = ParPos[d][ParID];
               }

               Par_MapMesh2Particles( EdgeL, EdgeR, _dh, VarSize, Var+P*CUBE(VarSize),
                                      amr->patch[0][lv][PID]->NPar, InterpParPos, ParType,
                                      amr->patch[0][lv][PID]->ParList, UseTracers_Yes, Var_Temp,
                                      TracerVelCorr_No );

//             4. update particles
               for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
               {
                  ParID = amr->patch[0][lv][PID]->ParList[p];

                  if ( ParType[ParID] != PTYPE_TRACER )
                     continue;

                  amr->Par->Mesh_Attr[v][ParID] = Var_Temp[p];
               }
            } // for (int PID=PID0, P=0; PID<PID0+8; PID++, P++)
         } // for (int v=0; v<ParAttrNum; v++)
      } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)


//    5. free memory
      delete [] Var;
      delete [] Var_Temp;

      Aux_DeallocateArray2D( InterpParPos );
   } // for (int lv=0; lv<NLEVEL; lv++)

} // FUNCTION : Par_Output_TracerParticle_Mesh



#endif // #ifdef PARTICLE
