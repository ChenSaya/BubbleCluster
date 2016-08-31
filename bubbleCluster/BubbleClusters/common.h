#pragma once
#include <Eigen/Dense>
extern "C"
{
#include "triangle.h"
}

#ifdef ANSI_DECLARATORS
extern "C" void triangulate(char *, struct triangulateio *, struct triangulateio *,
struct triangulateio *);
extern "C" void trifree(VOIDD *memptr);
#else /* not ANSI_DECLARATORS */
void triangulate();
void trifree();
#endif /* not ANSI_DECLARATORS */