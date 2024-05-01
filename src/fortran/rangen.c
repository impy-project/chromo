#include <numpy/random/bitgen.h>
#include <stdint.h>

typedef struct {
    // The structure is defined the same as in 
    // "numpy/numpy/random/src/pcg64/pcg64.h"
    // But the order seems swapped compared 
    // with __getstate__() output. 
    // It is not important here, but should be noted in case of problems
    uint64_t high;
    uint64_t low;
} pcg128_t;


typedef struct {
    pcg128_t state;
    pcg128_t inc;
} pcg_state_setseq_128;

typedef struct {
     pcg_state_setseq_128* pcg_state;
     int has_uint32;
     uint32_t uinteger;
} pcg64_state;

double random_standard_normal(bitgen_t *bitgen_state);

// called from Fortran
void npynxt_(double* value, int64_t* ptr)
{
    bitgen_t* bg = (bitgen_t*)(*ptr);
    *value = bg->next_double(bg->state);
}

// called from Fortran
void npynxt_get_state_(uint64_t* state_arr, int64_t* ptr)
{
    bitgen_t* bg = (bitgen_t*)(*ptr);
    pcg64_state *state = (pcg64_state *)bg->state;
    pcg_state_setseq_128* pcg_state128 = state->pcg_state;
    state_arr[0] = pcg_state128->state.high;
    state_arr[1] = pcg_state128->state.low;
    state_arr[2] = pcg_state128->inc.high;
    state_arr[3] = pcg_state128->inc.low;
}

void npynxt_set_state_(uint64_t* state_arr, int64_t* ptr)
{
    bitgen_t* bg = (bitgen_t*)(*ptr);
    pcg64_state *state = (pcg64_state *)bg->state;
    pcg_state_setseq_128* pcg_state128 = state->pcg_state;

    pcg_state128->state.high = state_arr[0];
    pcg_state128->state.low = state_arr[1];
    pcg_state128->inc.high = state_arr[2];
    pcg_state128->inc.low = state_arr[3];
}

// called from Fortran
void npygas_(double* value, int64_t* ptr)
{
    bitgen_t* bg = (bitgen_t*)(*ptr);
    *value = random_standard_normal(bg);
}
