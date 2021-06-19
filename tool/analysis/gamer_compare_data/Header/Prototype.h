#ifndef __PROTOTYPE_H__
#define __PROTOTYPE_H__



// CompareData functions
void ReadOption( int argc, char **argv );
void CheckParameter();
void CompareGridData();
void CompareParticleData();
void FreeMemory();

// GAMER functions
void Aux_Error( const char *File, const int Line, const char *Func, const char *Format, ... );
void Aux_Message( FILE *Type, const char *Format, ... );
bool Aux_CheckFileExist( const char *FileName );
void Aux_AllocateArray2D( real** &Array, const int J, const int I );
void Aux_DeallocateArray2D( real** &Array );
void LoadData( const char *FileName, AMR_t &amr, int &Format, int &NField, int &NMag, int &NParAtt, long &NPar, real **&ParData,
               char (*&FieldLabel)[MAX_STRING], char (*&MagLabel)[MAX_STRING], char (*&ParAttLabel)[MAX_STRING] );
void LoadData_HDF5( const char *FileName, AMR_t &amr, int &Format, int &NField, int &NMag, int &NParAtt, long &NPar, real **&ParData,
                    char (*&FieldLabel)[MAX_STRING], char (*&MagLabel)[MAX_STRING], char (*&ParAttLabel)[MAX_STRING] );
void SortParticle( const long NPar, const real *PosX, const real *PosY, const real *PosZ, const real *VelX, long *IdxTable );



#endif // __PROTOTYPE_H__
