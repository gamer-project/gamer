#ifndef __CUAPI_H__
#define __CUAPI_H__



#ifndef GAMER_DEBUG
#  define NDEBUG
#endif

#ifndef SERIAL
#  include <mpi.h>
#endif

#include <stdio.h>
#include <unistd.h>
#include "Macro.h"
#include "Typedef.h"
#include "Timer.h"
#include "SrcTerms.h"
#include "EoS.h"
#include "Global.h"
#include "PhysicalConstant.h"

#ifdef SERIAL
#  include "Serial.h"
#endif

#include "CUDA_CheckError.h"



#endif // #ifndef __CUAPI_H__
