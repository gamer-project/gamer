#ifndef __PROTOTYPE_H__
#define __PROTOTYPE_H__



// CompareData functions
void ReadOption( int argc, char **argv );
void CheckParameter();
int CompareGridData();
int CompareParticleData();
void FreeMemory();

// GAMER functions
void Aux_Error( const char *File, const int Line, const char *Func, const char *Format, ... );
void Aux_Message( FILE *Type, const char *Format, ... );
bool Aux_CheckFileExist( const char *FileName );
template <typename T> void Aux_AllocateArray2D( T** &Array, const int J, const int I );
template <typename T> void Aux_DeallocateArray2D( T** &Array );
void LoadData( const char *FileName, AMR_t &amr, int &Format, int &NField, int &NMag, int &NParAttFlt, int &NParAttInt,
               long &NPar, real_par **&ParFltData, long_par **&ParIntData, char (*&FieldLabel)[MAX_STRING], char (*&MagLabel)[MAX_STRING],
               char (*&ParAttFltLabel)[MAX_STRING], char (*&ParAttIntLabel)[MAX_STRING] );
void LoadData_HDF5( const char *FileName, AMR_t &amr, int &Format, int &NField, int &NMag, int &NParAttFlt, int &NParAttInt, long &NPar,
                    real_par **&ParFltData, long_par **&ParIntData, char (*&FieldLabel)[MAX_STRING], char (*&MagLabel)[MAX_STRING],
                    char (*&ParAttFltLabel)[MAX_STRING], char (*&ParAttIntLabel)[MAX_STRING] );
void SortParticle( const long NPar, const real_par *PosX, const real_par *PosY, const real_par *PosZ, const real_par *VelX, long *IdxTable );



#endif // __PROTOTYPE_H__
