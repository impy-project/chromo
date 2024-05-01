#include <numpy/random/bitgen.h>
#include <stdint.h>
#include <stdio.h>

// Define the structure for a 128-bit integer
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


void npynxt_get_state_(uint64_t* state_arr, int64_t* ptr)
{
    bitgen_t* bg = (bitgen_t*)(*ptr);
    pcg64_state *state = (pcg64_state *)bg->state;
    state_arr[0] = state->pcg_state->state.high;
    state_arr[1] = state->pcg_state->state.low;
    state_arr[2] = state->pcg_state->inc.high;
    state_arr[3] = state->pcg_state->inc.low;

    // printf("State High: %lu\n", state->pcg_state->state.high);
    // printf("State Low: %lu\n", state->pcg_state->state.low);
    // printf("Inc High: %lu\n", state->pcg_state->inc.high);
    // printf("Inc Low: %lu\n", state->pcg_state->inc.low);
}

void npynxt_set_state_(uint64_t* state_arr, int64_t* ptr)
{
    bitgen_t* bg = (bitgen_t*)(*ptr);
    pcg64_state *state = (pcg64_state *)bg->state;

    // printf("Before setting:");
    // printf("State High: %lu\n", state->pcg_state->state.high);
    // printf("State Low: %lu\n", state->pcg_state->state.low);
    // printf("Inc High: %lu\n", state->pcg_state->inc.high);
    // printf("Inc Low: %lu\n", state->pcg_state->inc.low);

    state->pcg_state->state.high = state_arr[0];
    state->pcg_state->state.low = state_arr[1];
    state->pcg_state->inc.high = state_arr[2];
    state->pcg_state->inc.low = state_arr[3];

    // printf("After setting:");
    // printf("State High: %lu\n", state->pcg_state->state.high);
    // printf("State Low: %lu\n", state->pcg_state->state.low);
    // printf("Inc High: %lu\n", state->pcg_state->inc.high);
    // printf("Inc Low: %lu\n", state->pcg_state->inc.low);



    // pcg_state_setseq_128 pstate = state->pcg_state
    // pstate.state
    // pcg_state* state = (pcg_state*)bg->state;
    // printf("%d", *(int64_t*)bg->state->state);

    // Extract the high and low parts of state and inc
    // state_values[0] = state->pcg_state->state.high;
    // state_values[1] = state->pcg_state->state.low;
    // state_values[2] = state->pcg_state->inc.high;
    // state_values[3] = state->pcg_state->inc.low;

    // printf("State High: %lu\n", state->pcg_state->state.high);
    // printf("State Low: %lu\n", state->pcg_state->state.low);
    // printf("Inc High: %lu\n", state->pcg_state->inc.high);
    // printf("Inc Low: %lu\n", state->pcg_state->inc.low);

}

// called from Fortran
void npygas_(double* value, int64_t* ptr)
{
    bitgen_t* bg = (bitgen_t*)(*ptr);
    *value = random_standard_normal(bg);
}
