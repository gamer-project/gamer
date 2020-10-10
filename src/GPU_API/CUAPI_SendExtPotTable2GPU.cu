#include "CUAPI.h"
#include "CUPOT.h"

#if ( defined GPU  &&  defined GRAVITY )


// device pointer
extern real *d_ExtPotTable;




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_SendExtPotTable2GPU
// Description :  Send the external potential table to GPU
//
// Note        :  1. Invoked by Init_LoadExtPotTable()
//
// Parameter   :  h_Table : Host array storing the input table
//
// Return      :  d_ExtPotTable
//-------------------------------------------------------------------------------------------------------
void CUAPI_SendExtPotTable2GPU( const real *h_Table )
{

   const long MemSize = (long)sizeof(real)*EXT_POT_TABLE_NPOINT[0]*EXT_POT_TABLE_NPOINT[1]*EXT_POT_TABLE_NPOINT[2];

   CUDA_CHECK_ERROR(  cudaMemcpyAsync( d_ExtPotTable, h_Table, MemSize, cudaMemcpyHostToDevice )  );

} // FUNCTION : CUAPI_SendExtPotTable2GPU



#endif // #if ( defined GPU  &&  defined GRAVITY )
