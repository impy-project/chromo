#include <numpy/random/bitgen.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// called from Fortran
void npynxt_(double* value, int64_t* buffer)
{
    bitgen_t* bg = (bitgen_t*)(*buffer);
    *value = bg->next_double(bg->state);
}