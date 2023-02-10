#include <numpy/random/bitgen.h>
#include <stdint.h>

double random_standard_normal(bitgen_t *bitgen_state);

// called from Fortran
void npynxt_(double* value, int64_t* ptr)
{
    bitgen_t* bg = (bitgen_t*)(*ptr);
    *value = bg->next_double(bg->state);
}

// called from Fortran
void npygas_(double* value, int64_t* ptr)
{
    bitgen_t* bg = (bitgen_t*)(*ptr);
    *value = random_standard_normal(bg);
}
