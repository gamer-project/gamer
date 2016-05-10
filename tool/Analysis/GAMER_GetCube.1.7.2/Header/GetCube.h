#ifndef __GETCUBE_H__
#define __GETCUBE_H__



#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <unistd.h>
#include <cstdio>
#ifndef SERIAL
#include <mpi.h>
#endif

using namespace std;

#include "TypeDef.h"
#include "Tree.h"
#include "Global.h"
#include "Prototype.h"

#ifdef SUPPORT_HDF5
#  include "HDF5_Typedef.h"
#endif

#ifdef SERIAL
#  include "Serial.h"
#endif

#ifdef OPENMP
#  include <omp.h>
#endif



#endif // #ifndef __GETCUBE_H__
