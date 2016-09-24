#ifndef __PROTOTYPE_H__
#define __PROTOTYPE_H__



// Data2Visit functions
void TruncateBox();
void ReadOption( int, char ** );
void CheckParameter();
void SetDefaultParameter();
void WriteLevel();
void WriteRoot();
void WriteBox();
void CreateSilo();

// GAMER functions
void LoadData();
void LoadData_HDF5();
bool Aux_CheckFileExist( const char *FileName );
void Aux_Message( FILE *Type, const char *Format, ... );
void Aux_Error( const char *File, const int Line, const char *Func, const char *Format, ... );



#endif // __PROTOTYPE_H__
