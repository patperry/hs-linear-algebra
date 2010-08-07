
#include "config.h"

#define xstr(s) str(s)
#define str(s) #s
#define hsc_f77_func(name) printf("\""xstr(F77_FUNC(name))"\"");
