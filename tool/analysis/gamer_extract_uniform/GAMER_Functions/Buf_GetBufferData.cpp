#include "ExtractUniform.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Buf_GetBufferData
// Description :  Fill up the data of the buffer patches
//
// Parameter   :  lv          :  The targeted refinement level to get data
//                Package     :  1 -> rho + momentum + energy
//                               2 -> potential
//                               3 -> density [USELESS HERE]
//                               4 -> particle (or total) density
//                ParaBuffer  :  the width of data te be sent in the buffer patches
//-------------------------------------------------------------------------------------------------------
void Buf_GetBufferData( const int lv, const int Package, const int ParaBuffer )
{

   const int TargetSib[26] = { 0,1,2,3,4,5,6,9,7,8,10,13,11,12,14,17,16,15,18,25,19,24,20,23,21,22 };

   int LoopWidth[3], Disp[2][3], NSpecies, ii, jj, kk, sib, TargetRank[2];
   int SendSize[2], RecvSize[2], ID, PID;
   real *SendBuffer[2], *RecvBuffer[2];


// set the number of variables to be sent in a single cell
   switch ( Package )
   {
      case 1:
         NSpecies = NCOMP_TOTAL;
         break;

      case 2: //case 3:
         NSpecies = 1;
         break;

      case 4:
         NSpecies = 1;
         break;

      default:
         fprintf( stderr, "ERROR : \"incorrect parameter %s = %d\" !!\n", "Package", Package );
         fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n", __FILE__, __LINE__,  __FUNCTION__  );
         MPI_Exit();
   }


   for (int s=0; s<26; s+=2)
   {
      for (int d=0; d<3; d++)
         LoopWidth[d] = TABLE_01( TargetSib[s], 'x'+d, ParaBuffer, PATCH_SIZE, ParaBuffer );

//    prepare the SendBuffer and RecvBuffer
      for (int t=0; t<2; t++)
      {
         sib            = TargetSib[s+t];

         for (int d=0; d<3; d++)
            Disp[t][d]  = TABLE_01( sib, 'x'+d, 0, 0, PATCH_SIZE-ParaBuffer );

         TargetRank[t]  = SibRank[sib];

         SendSize[t]    = ParaVar.SendP_NList[lv][sib]*LoopWidth[0]*LoopWidth[1]*LoopWidth[2]*NSpecies;
         RecvSize[t]    = ParaVar.RecvP_NList[lv][sib]*LoopWidth[0]*LoopWidth[1]*LoopWidth[2]*NSpecies;

         SendBuffer[t]  = new real [ SendSize[t] ];
         RecvBuffer[t]  = new real [ RecvSize[t] ];


//       copy the data into the SendBuffer
         switch ( Package )
         {
            case 1:
               for (int P=0; P<ParaVar.SendP_NList[lv][sib]; P++)    {  PID = ParaVar.SendP_IDList[lv][sib][P];
               for (int v=0; v<NSpecies; v++)                        {
               for (int k=0; k<LoopWidth[2]; k++)                    {  kk  = k + Disp[t][2];
               for (int j=0; j<LoopWidth[1]; j++)                    {  jj  = j + Disp[t][1];
               for (int i=0; i<LoopWidth[0]; i++)                    {  ii  = i + Disp[t][0];

                  ID = NSpecies * ( P*LoopWidth[0]*LoopWidth[1]*LoopWidth[2] + k*LoopWidth[0]*LoopWidth[1] +
                                    j*LoopWidth[0] + i );

                     SendBuffer[t][ID+v] = amr.patch[lv][PID]->fluid[v][kk][jj][ii];

               }}}}}
               break;

            case 2:
               for (int P=0; P<ParaVar.SendP_NList[lv][sib]; P++)    {  PID = ParaVar.SendP_IDList[lv][sib][P];
               for (int k=0; k<LoopWidth[2]; k++)                    {  kk  = k + Disp[t][2];
               for (int j=0; j<LoopWidth[1]; j++)                    {  jj  = j + Disp[t][1];
               for (int i=0; i<LoopWidth[0]; i++)                    {  ii  = i + Disp[t][0];

                  ID = NSpecies * ( P*LoopWidth[0]*LoopWidth[1]*LoopWidth[2] + k*LoopWidth[0]*LoopWidth[1] +
                                    j*LoopWidth[0] + i );

                  SendBuffer[t][ID] = amr.patch[lv][PID]->pot[kk][jj][ii];

               }}}}
               break;

/*
            case 3:
               for (int P=0; P<ParaVar.SendP_NList[lv][sib]; P++)    {  PID = ParaVar.SendP_IDList[lv][sib][P];
               for (int k=0; k<LoopWidth[2]; k++)                    {  kk  = k + Disp[t][2];
               for (int j=0; j<LoopWidth[1]; j++)                    {  jj  = j + Disp[t][1];
               for (int i=0; i<LoopWidth[0]; i++)                    {  ii  = i + Disp[t][0];

                  ID = NSpecies * ( P*LoopWidth[0]*LoopWidth[1]*LoopWidth[2] + k*LoopWidth[0]*LoopWidth[1] +
                                    j*LoopWidth[0] + i );

                  SendBuffer[t][ID] = amr.patch[lv][PID]->fluid[kk][jj][ii][0];

               }}}}
               break;
*/

            case 4:
               for (int P=0; P<ParaVar.SendP_NList[lv][sib]; P++)    {  PID = ParaVar.SendP_IDList[lv][sib][P];
               for (int k=0; k<LoopWidth[2]; k++)                    {  kk  = k + Disp[t][2];
               for (int j=0; j<LoopWidth[1]; j++)                    {  jj  = j + Disp[t][1];
               for (int i=0; i<LoopWidth[0]; i++)                    {  ii  = i + Disp[t][0];

                  ID = NSpecies * ( P*LoopWidth[0]*LoopWidth[1]*LoopWidth[2] + k*LoopWidth[0]*LoopWidth[1] +
                                    j*LoopWidth[0] + i );

                  SendBuffer[t][ID] = amr.patch[lv][PID]->par_dens[kk][jj][ii];

               }}}}
               break;

         } // switch ( Package )
      } // for (int t=0; t<2; t++)


//    transfer the buffer data between different ranks
      MPI_ExchangeInfo( TargetRank, SendSize, RecvSize, SendBuffer, RecvBuffer );


//    copy the data in the RecvBuffer back to the amr.patch pointer
      for (int t=0; t<2; t++)
      {
         sib = TargetSib[s+t];

         switch ( Package )
         {
            case 1:
               for (int P=0; P<ParaVar.RecvP_NList[lv][sib]; P++)    {  PID = ParaVar.RecvP_IDList[lv][sib][P];
               for (int v=0; v<NSpecies; v++)                        {
               for (int k=0; k<LoopWidth[2]; k++)                    {  kk  = k + Disp[!t][2];
               for (int j=0; j<LoopWidth[1]; j++)                    {  jj  = j + Disp[!t][1];
               for (int i=0; i<LoopWidth[0]; i++)                    {  ii  = i + Disp[!t][0];

                  ID = NSpecies * ( P*LoopWidth[0]*LoopWidth[1]*LoopWidth[2] + k*LoopWidth[0]*LoopWidth[1] +
                                    j*LoopWidth[0] + i );

                     amr.patch[lv][PID]->fluid[v][kk][jj][ii] = RecvBuffer[t][ID+v];

               }}}}}
               break;

            case 2:
               for (int P=0; P<ParaVar.RecvP_NList[lv][sib]; P++)    {  PID = ParaVar.RecvP_IDList[lv][sib][P];
               for (int k=0; k<LoopWidth[2]; k++)                    {  kk  = k + Disp[!t][2];
               for (int j=0; j<LoopWidth[1]; j++)                    {  jj  = j + Disp[!t][1];
               for (int i=0; i<LoopWidth[0]; i++)                    {  ii  = i + Disp[!t][0];

                  ID = NSpecies * ( P*LoopWidth[0]*LoopWidth[1]*LoopWidth[2] + k*LoopWidth[0]*LoopWidth[1] +
                                    j*LoopWidth[0] + i );

                  amr.patch[lv][PID]->pot[kk][jj][ii] = RecvBuffer[t][ID];

               }}}}
               break;

/*
            case 3:
               for (int P=0; P<ParaVar.RecvP_NList[lv][sib]; P++)    {  PID = ParaVar.RecvP_IDList[lv][sib][P];
               for (int k=0; k<LoopWidth[2]; k++)                    {  kk  = k + Disp[!t][2];
               for (int j=0; j<LoopWidth[1]; j++)                    {  jj  = j + Disp[!t][1];
               for (int i=0; i<LoopWidth[0]; i++)                    {  ii  = i + Disp[!t][0];

                  ID = NSpecies * ( P*LoopWidth[0]*LoopWidth[1]*LoopWidth[2] + k*LoopWidth[0]*LoopWidth[1] +
                                    j*LoopWidth[0] + i );

                  amr.patch[lv][PID]->fluid[kk][jj][ii][0] = RecvBuffer[t][ID];

               }}}}
               break;
*/

            case 4:
               for (int P=0; P<ParaVar.RecvP_NList[lv][sib]; P++)    {  PID = ParaVar.RecvP_IDList[lv][sib][P];
               for (int k=0; k<LoopWidth[2]; k++)                    {  kk  = k + Disp[!t][2];
               for (int j=0; j<LoopWidth[1]; j++)                    {  jj  = j + Disp[!t][1];
               for (int i=0; i<LoopWidth[0]; i++)                    {  ii  = i + Disp[!t][0];

                  ID = NSpecies * ( P*LoopWidth[0]*LoopWidth[1]*LoopWidth[2] + k*LoopWidth[0]*LoopWidth[1] +
                                    j*LoopWidth[0] + i );

                  amr.patch[lv][PID]->par_dens[kk][jj][ii] = RecvBuffer[t][ID];

               }}}}
               break;

         } // switch ( Package )

         delete [] SendBuffer[t];
         delete [] RecvBuffer[t];

      } // for (int t=0; t<2; t++)

   } // for (int s=0; s<26; s+=2)

}
