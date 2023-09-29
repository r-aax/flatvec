// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#if !defined(LINUX_GCC_BUILD) && !defined(LINUX_ICC_BUILD)
#include "targetver.h"
#define WIN32_LEAN_AND_MEAN // Exclude rarely-used stuff from Windows headers
#endif

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <stack>
#include <queue>
#include <list>
#include <algorithm>
#include <functional>
#include <chrono>

#ifdef LINUX_GCC_BUILD
#include <math.h>
#endif

// TODO: reference additional headers your program requires here
