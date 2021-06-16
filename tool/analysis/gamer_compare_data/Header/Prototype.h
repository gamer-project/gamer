#ifndef __PROTOTYPE_H__
#define __PROTOTYPE_H__



// CompareData functions
void ReadOption( int argc, char **argv );
void CheckParameter();
void CompareGridData();
void CompareParticleData();

// GAMER functions
void Aux_Error( const char *File, const int Line, const char *Func, const char *Format, ... );
void Aux_Message( FILE *Type, const char *Format, ... );
bool Aux_CheckFileExist( const char *FileName );
void Aux_AllocateArray2D( real** &Array, const int J, const int I );
void Aux_DeallocateArray2D( real** &Array );
void LoadData( AMR_t &patch, const char *FileName, bool &WithPot, int &WithParDens, bool &WithPar,
               int &NParVarOut, long &NPar, real **&ParData, bool &WithMagCC, bool &WithMagFC );
void LoadData_HDF5( AMR_t &amr, const char *FileName, bool &WithPot, int &WithParDens, bool &WithPar,
                    int &NParVarOut, long &NPar, real **&ParData, bool &WithMagCC, bool &WithMagFC );
void SortParticle( const long NPar, const real *PosX, const real *PosY, const real *PosZ, const real *VelX, long *IdxTable );



#endif // __PROTOTYPE_H__
