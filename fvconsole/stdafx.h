// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#if !defined(LINUX_GCC_BUILD) && !defined(LINUX_ICC_BUILD)
#include "targetver.h"
#include <tchar.h>
#endif

#include <stdio.h>
#include <iostream>
#include <fstream>

#ifndef LINUX_ICC_BUILD
#include <filesystem>
#endif

#ifdef LINUX_GCC_BUILD
#include <math.h>
#endif

// TODO: reference additional headers your program requires here
