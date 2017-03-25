#ifndef EXTENDEDKF_GET_TIME_H
#define EXTENDEDKF_GET_TIME_H

#ifdef _WIN32
#include <Windows.h>
#else
#include <sys/time.h>
#include <ctime>
#endif

/* Remove if already defined */
typedef long long int64;
typedef unsigned long long uint64;

/* Returns the amount of milliseconds elapsed since the UNIX epoch. Works on both
 * windows and linux. */

uint64 GetTimeMs64();

#endif //EXTENDEDKF_GET_TIME_H
