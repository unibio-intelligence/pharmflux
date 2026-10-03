#include <R_ext/Rdynload.h>
extern void R_init_pharmflux_extendr(DllInfo *dll);
void R_init_pharmflux(DllInfo *dll) { R_init_pharmflux_extendr(dll); }
