#include <numpy/random/bitgen.h>
#include <stdint.h>
#include <string.h>

// Maximum size for storing any bit generator state (in uint64_t units)
#define MAX_STATE_SIZE 628

// Generator type IDs (must match Python side)
#define GEN_ID_PCG64     1
#define GEN_ID_PCG64DXSM 2
#define GEN_ID_MT19937   3
#define GEN_ID_PHILOX    4
#define GEN_ID_SFC64     5

// PCG64/PCG64DXSM structures
typedef struct {
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

// MT19937 structure
typedef struct {
    uint64_t key[624];
    int pos;
    int has_uint32;  
    uint32_t uinteger; 
} mt19937_state;

// Philox structure
typedef struct {
    uint64_t counter[4];
    uint64_t key[2];
    int buffer_pos;
    uint64_t buffer[4];
    int has_uint32;        
    uint32_t uinteger;     
} philox_state;

// SFC64 structure
typedef struct {
    uint64_t s[4];
    int has_uint32;
    uint32_t uinteger;
} sfc64_state;

double random_standard_normal(bitgen_t *bitgen_state);

// called from Fortran
void npynxt_(double* value, int64_t* ptr)
{
    bitgen_t* bg = (bitgen_t*)(*ptr);
    *value = bg->next_double(bg->state);
}

// Save generator state using explicit generator ID
void npynxt_get_state_(uint64_t* state_arr, int64_t* ptr, int* gen_id)
{
    bitgen_t* bg = (bitgen_t*)(*ptr);
    
    // Clear the state array
    memset(state_arr, 0, MAX_STATE_SIZE * sizeof(uint64_t));
    
    // Save state based on generator type
    switch (*gen_id) {
        case GEN_ID_PCG64:
        case GEN_ID_PCG64DXSM: {
            pcg64_state *pcg_state = (pcg64_state *)bg->state;
            pcg_state_setseq_128* pcg_state128 = pcg_state->pcg_state;
            state_arr[0] = pcg_state128->state.high;
            state_arr[1] = pcg_state128->state.low;
            state_arr[2] = pcg_state128->inc.high;
            state_arr[3] = pcg_state128->inc.low;
            state_arr[4] = (uint64_t)pcg_state->has_uint32;
            state_arr[5] = (uint64_t)pcg_state->uinteger;
            memcpy(&state_arr[6], "PCG64\0\0\0", 8);
            break;
        }
        
        case GEN_ID_MT19937: {
            mt19937_state *mt_state = (mt19937_state *)bg->state;
            memcpy(state_arr, mt_state->key, sizeof(mt_state->key));
            state_arr[624] = (uint64_t)mt_state->pos;
            state_arr[625] = (uint64_t)mt_state->has_uint32;
            state_arr[626] = (uint64_t)mt_state->uinteger;
            memcpy(&state_arr[627], "MT19937\0", 8);
            break;
        }
        
        case GEN_ID_PHILOX: {
            philox_state *philox = (philox_state *)bg->state;
            memcpy(state_arr, philox->counter, sizeof(philox->counter));
            memcpy(&state_arr[4], philox->key, sizeof(philox->key));
            state_arr[6] = (uint64_t)philox->buffer_pos;
            memcpy(&state_arr[7], philox->buffer, sizeof(philox->buffer));
            state_arr[11] = (uint64_t)philox->has_uint32;
            state_arr[12] = (uint64_t)philox->uinteger;
            memcpy(&state_arr[13], "PHILOX\0\0", 8);
            break;
        }
        
        case GEN_ID_SFC64: {
            sfc64_state *sfc_state = (sfc64_state *)bg->state;
            memcpy(state_arr, sfc_state->s, sizeof(sfc_state->s));
            state_arr[4] = (uint64_t)sfc_state->has_uint32;
            state_arr[5] = (uint64_t)sfc_state->uinteger;
            memcpy(&state_arr[6], "SFC64\0\0\0", 8);
            break;
        }
        
        default:
            // Unknown generator ID - default to PCG64 format
            pcg64_state *pcg_state = (pcg64_state *)bg->state;
            pcg_state_setseq_128* pcg_state128 = pcg_state->pcg_state;
            state_arr[0] = pcg_state128->state.high;
            state_arr[1] = pcg_state128->state.low;
            state_arr[2] = pcg_state128->inc.high;
            state_arr[3] = pcg_state128->inc.low;
            state_arr[4] = (uint64_t)pcg_state->has_uint32;
            state_arr[5] = (uint64_t)pcg_state->uinteger;
            memcpy(&state_arr[6], "PCG64\0\0\0", 8);
            break;
    }
}

// Restore generator state
void npynxt_set_state_(uint64_t* state_arr, int64_t* ptr)
{
    bitgen_t* bg = (bitgen_t*)(*ptr);
    char gen_type[8] = {0};
    
    // Detect generator type from saved marker and restore state
    
    // Check for PCG64 marker
    memcpy(gen_type, &state_arr[6], 8);
    if (memcmp(gen_type, "PCG64", 5) == 0) {
        pcg64_state *pcg_state = (pcg64_state *)bg->state;
        pcg_state_setseq_128* pcg_state128 = pcg_state->pcg_state;
        pcg_state128->state.high = state_arr[0];
        pcg_state128->state.low = state_arr[1];
        pcg_state128->inc.high = state_arr[2];
        pcg_state128->inc.low = state_arr[3];
        pcg_state->has_uint32 = (int)state_arr[4];
        pcg_state->uinteger = (uint32_t)state_arr[5];
        return;
    }
    
    // Check for MT19937 marker
    memcpy(gen_type, &state_arr[627], 8);
    if (memcmp(gen_type, "MT19937", 7) == 0) {
        mt19937_state *mt_state = (mt19937_state *)bg->state;
        memcpy(mt_state->key, state_arr, sizeof(mt_state->key));
        mt_state->pos = (int)state_arr[624];
        mt_state->has_uint32 = (int)state_arr[625];
        mt_state->uinteger = (uint32_t)state_arr[626];
        return;
    }
    
    // Check for Philox marker
    memcpy(gen_type, &state_arr[13], 8);
    if (memcmp(gen_type, "PHILOX", 6) == 0) {
        philox_state *philox = (philox_state *)bg->state;
        memcpy(philox->counter, state_arr, sizeof(philox->counter));
        memcpy(philox->key, &state_arr[4], sizeof(philox->key));
        philox->buffer_pos = (int)state_arr[6];
        memcpy(philox->buffer, &state_arr[7], sizeof(philox->buffer));
        philox->has_uint32 = (int)state_arr[11];
        philox->uinteger = (uint32_t)state_arr[12];
        return;
    }
    
    // Check for SFC64 marker
    memcpy(gen_type, &state_arr[6], 8);
    if (memcmp(gen_type, "SFC64", 5) == 0) {
        sfc64_state *sfc_state = (sfc64_state *)bg->state;
        memcpy(sfc_state->s, state_arr, sizeof(sfc_state->s));
        sfc_state->has_uint32 = (int)state_arr[4];
        sfc_state->uinteger = (uint32_t)state_arr[5];
        return;
    }
    
    // No recognized marker - assume legacy PCG64 format
    pcg64_state *pcg_state = (pcg64_state *)bg->state;
    if (pcg_state->pcg_state != NULL) {
        pcg_state_setseq_128* pcg_state128 = pcg_state->pcg_state;
        pcg_state128->state.high = state_arr[0];
        pcg_state128->state.low = state_arr[1];
        pcg_state128->inc.high = state_arr[2];
        pcg_state128->inc.low = state_arr[3];
    }
}

// called from Fortran
void npygas_(double* value, int64_t* ptr)
{
    bitgen_t* bg = (bitgen_t*)(*ptr);
    *value = random_standard_normal(bg);
}
