//--------------------------------------------------
// sse2_wrapper.h — platform-portable 128-bit integer SIMD wrappers.
// x86/x86_64: uses SSE2 (__m128i via emmintrin.h)
// arm64:      uses NEON (uint32x4_t via arm_neon.h)
//--------------------------------------------------
#ifndef SSE2_WRAPPER_H
#define SSE2_WRAPPER_H

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
// ===== ARM NEON implementation =====
#include <arm_neon.h>
#include <stdint.h>
#include <string.h>  // memcpy for aligned load/store fallback

typedef uint32x4_t b128;

typedef union {
  int _mem[4];
  b128 _b128;
} union_b128_mem;

static inline b128 set_ints_b128(int i3, int i2, int i1, int i0){
  int32_t v[4] = {i0, i1, i2, i3};
  return vreinterpretq_u32_s32(vld1q_s32(v));
}
static inline b128 set_first_int_b128(int i0){
  uint32x4_t z = vdupq_n_u32(0);
  return vsetq_lane_u32((uint32_t)i0, z, 0);
}
static inline b128 set_zero_b128(){ return vdupq_n_u32(0); }
static inline b128 set_all_bytes(char b){
  return vreinterpretq_u32_u8(vdupq_n_u8((uint8_t)b));
}
static inline b128 set_all_shorts(short b){
  return vreinterpretq_u32_s16(vdupq_n_s16(b));
}
static inline b128 set_all_ints(int b){
  return vreinterpretq_u32_s32(vdupq_n_s32(b));
}

// _mm_extract_epi16 / _mm_insert_epi16 equivalents (pos must be 0-7)
#define get_immediate_int16_b128(__A, __IMM) \
  (int)vgetq_lane_s16(vreinterpretq_s16_u32(__A), __IMM)
#define set_immediate_int16_b128(__A, __B, __IMM) \
  vreinterpretq_u32_s16(vsetq_lane_s16((int16_t)(__B), vreinterpretq_s16_u32(__A), __IMM))

static inline int get_int_0_b128(b128 a){ return (int)vgetq_lane_u32(a,0); }
static inline int get_int_1_b128(b128 a){ return (int)vgetq_lane_u32(a,1); }
static inline int get_int_2_b128(b128 a){ return (int)vgetq_lane_u32(a,2); }
static inline int get_int_3_b128(b128 a){ return (int)vgetq_lane_u32(a,3); }

static inline b128 set_int_0_b128(b128 a, int b){ return vsetq_lane_u32((uint32_t)b,a,0); }
static inline b128 set_int_1_b128(b128 a, int b){ return vsetq_lane_u32((uint32_t)b,a,1); }
static inline b128 set_int_2_b128(b128 a, int b){ return vsetq_lane_u32((uint32_t)b,a,2); }
static inline b128 set_int_3_b128(b128 a, int b){ return vsetq_lane_u32((uint32_t)b,a,3); }

static inline int get_int_b128(b128 a, int intpos){
  union_b128_mem u; u._b128 = a; return u._mem[intpos];
}
static inline b128 set_int_b128(b128 a, int b, int intpos){
  union_b128_mem u; u._b128 = a; u._mem[intpos] = b; return u._b128;
}

static inline int get_bit_b128(b128 a, int bitpos){
  return (( get_int_b128(a,bitpos>>5) >> (bitpos & 0x1F)) & 0x1);
}
static inline b128 set_bit_b128(b128 a, int bitpos, int bitvalue){
  union_b128_mem u; u._b128 = a;
  u._mem[bitpos>>5] |= (0x1 << (bitpos & 0x1F));
  (void)bitvalue;
  return u._b128;
}

static inline int equal_b128(b128 a, b128 b){
  uint32x4_t c = vceqq_u32(a,b);
  return vgetq_lane_u32(c,0) | vgetq_lane_u32(c,1) |
         vgetq_lane_u32(c,2) | vgetq_lane_u32(c,3);
}

static inline b128 add_b128(b128 a, b128 b){ return vaddq_u32(a,b); }
static inline b128 sub_b128(b128 a, b128 b){ return vsubq_u32(a,b); }
static inline b128 xor_b128(b128 a, b128 b){
  return vreinterpretq_u32_u8(veorq_u8(vreinterpretq_u8_u32(a),vreinterpretq_u8_u32(b)));
}
static inline b128 or_b128(b128 a, b128 b){
  return vreinterpretq_u32_u8(vorrq_u8(vreinterpretq_u8_u32(a),vreinterpretq_u8_u32(b)));
}
static inline b128 and_b128(b128 a, b128 b){
  return vreinterpretq_u32_u8(vandq_u8(vreinterpretq_u8_u32(a),vreinterpretq_u8_u32(b)));
}
static inline b128 andnot_b128(b128 a, b128 b){
  // SSE2 andnot: (~a) & b  =>  NEON: vbicq (b & ~a)
  return vreinterpretq_u32_u8(vbicq_u8(vreinterpretq_u8_u32(b),vreinterpretq_u8_u32(a)));
}
static inline b128 negate_b128(b128 a){
  return vreinterpretq_u32_u8(vmvnq_u8(vreinterpretq_u8_u32(a)));
}

// Byte shifts: _mm_slli_si128 shifts bytes LEFT (higher addresses)
// vextq_u8(a,b,n) = concat(a,b) shifted right by n bytes
#define shift_bytes_left_b128(__A, __IMM_BYTES) \
  vreinterpretq_u32_u8(vextq_u8(vdupq_n_u8(0), vreinterpretq_u8_u32(__A), (16-(__IMM_BYTES))))
#define shift_bytes_right_b128(__A, __IMM_BYTES) \
  vreinterpretq_u32_u8(vextq_u8(vreinterpretq_u8_u32(__A), vdupq_n_u8(0), (__IMM_BYTES)))

// Per-element 32-bit shifts by shift count in lowest lane of b
static inline b128 shift_each32_bits_right_b128(b128 a, b128 b){
  int32_t sh = -(int32_t)vgetq_lane_u32(b,0);
  return vshlq_u32(a, vdupq_n_s32(sh));
}
static inline b128 shift_each32_bits_left_b128(b128 a, b128 b){
  int32_t sh = (int32_t)vgetq_lane_u32(b,0);
  return vshlq_u32(a, vdupq_n_s32(sh));
}

b128 *alloc_b128(int num_b128);
b128 *calloc_b128(int num_b128);
void  free_b128(b128 *a);

static inline b128 get_from_mem_b128(const b128 *ptr, int pos){
  b128 v; memcpy(&v, ptr+pos, 16); return v;
}
static inline b128 get_b128(const b128 *ptr){ b128 v; memcpy(&v, ptr, 16); return v; }
static inline void set_b128(b128 *ptr, b128 val){ memcpy(ptr, &val, 16); }
static inline void set_b128_NOCACHEPOLLUTE(b128 *ptr, b128 val){ memcpy(ptr, &val, 16); }

static inline int get_bit_in_mem_b128(const b128 *ptr, int bitpos){
  return get_bit_b128(get_from_mem_b128(ptr,(bitpos>>7)), (bitpos & 0x7F));
}
static inline void set_bit_in_mem_b128(b128 *ptr, int bitpos, int bitvalue){
  b128 *b128pos = ptr+(bitpos>>7);
  b128 v = get_b128(b128pos);
  set_b128(b128pos, set_bit_b128(v, (bitpos & 0x7F), bitvalue));
}

static inline void mem_or_b128(b128 *dst, const b128 *src, const int num_b128s){
  for(int i=0;i<num_b128s;i++) dst[i]=or_b128(dst[i],src[i]);
}
static inline void mem_and_b128(b128 *dst, const b128 *src, const int num_b128s){
  for(int i=0;i<num_b128s;i++) dst[i]=and_b128(dst[i],src[i]);
}
static inline void mem_xor_b128(b128 *dst, const b128 *src, const int num_b128s){
  for(int i=0;i<num_b128s;i++) dst[i]=xor_b128(dst[i],src[i]);
}
static inline void mem_andnot_b128(b128 *dst, const b128 *src, const int num_b128s){
  for(int i=0;i<num_b128s;i++) dst[i]=andnot_b128(src[i],dst[i]);
}
static inline void mem_reversed_andnot_b128(b128 *dst, const b128 *src, const int num_b128s){
  for(int i=0;i<num_b128s;i++) dst[i]=andnot_b128(dst[i],src[i]);
}

void print_ints_b128(b128 a);
void print_bits_32(int a);
void print_bits_b128(b128 a);
void print_blocks_b128(b128 a, int block_size);

// No rdtsc on ARM — stub out
typedef unsigned long long ticks;
static inline ticks getticks(void){ return 0ULL; }

// Prefetch and fence stubs for arm64 (hardware prefetcher handles this)
#define _MM_HINT_NTA 0
#define _MM_HINT_T0  1
#define _MM_HINT_T1  2
#define _MM_HINT_T2  3
static inline void _mm_prefetch(const void *p, int hint){ (void)p; (void)hint; }
static inline void _mm_lfence(void){}

#else
// ===== x86/x86_64 SSE2 implementation (original) =====
#include <emmintrin.h>

typedef __m128i b128;

typedef union {
  int _mem[4];
  b128 _b128;
} union_b128_mem;

static __inline b128 set_ints_b128(int i3, int i2, int i1, int i0){
  return _mm_set_epi32(i3,i2,i1,i0);
}
static __inline b128 set_first_int_b128(int i0){ return _mm_cvtsi32_si128(i0); }
static __inline b128 set_zero_b128(){ return _mm_setzero_si128(); }
static __inline b128 set_all_bytes(char b){ return _mm_set1_epi8(b); }
static __inline b128 set_all_shorts(short b){ return _mm_set1_epi16(b); }
static __inline b128 set_all_ints(int b){ return _mm_set1_epi32(b); }

#define get_immediate_int16_b128(__A, __IMM) _mm_extract_epi16(__A,__IMM)
#define set_immediate_int16_b128(__A, __B, __IMM) _mm_insert_epi16(__A,__B,__IMM)

static __inline int get_int_0_b128(b128 a){
  return ((_mm_extract_epi16(a,1)<<16) | _mm_extract_epi16(a,0));
}
static __inline int get_int_1_b128(b128 a){
  return ((_mm_extract_epi16(a,3)<<16) | _mm_extract_epi16(a,2));
}
static __inline int get_int_2_b128(b128 a){
  return ((_mm_extract_epi16(a,5)<<16) | _mm_extract_epi16(a,4));
}
static __inline int get_int_3_b128(b128 a){
  return ((_mm_extract_epi16(a,7)<<16) | _mm_extract_epi16(a,6));
}

static __inline b128 set_int_0_b128(b128 a, int b){
  return _mm_insert_epi16(_mm_insert_epi16(a,(b>>16),1),b,0);
}
static __inline b128 set_int_1_b128(b128 a, int b){
  return _mm_insert_epi16(_mm_insert_epi16(a,(b>>16),3),b,2);
}
static __inline b128 set_int_2_b128(b128 a, int b){
  return _mm_insert_epi16(_mm_insert_epi16(a,(b>>16),5),b,4);
}
static __inline b128 set_int_3_b128(b128 a, int b){
  return _mm_insert_epi16(_mm_insert_epi16(a,(b>>16),7),b,6);
}

static __inline int get_int_b128(b128 a, int intpos){
  union_b128_mem b; b._b128 = a; return b._mem[intpos];
}
static __inline b128 set_int_b128(b128 a, int b, int intpos){
  union_b128_mem c; c._b128 = a; c._mem[intpos] = b; return c._b128;
}

static __inline int get_bit_b128(b128 a, int bitpos){
  return ((get_int_b128(a,bitpos>>5) >> (bitpos & 0x1F)) & 0x1);
}
static __inline b128 set_bit_b128(b128 a, int bitpos, int bitvalue){
  union_b128_mem b; b._b128 = a;
  b._mem[bitpos>>5] |= (0x1 << (bitpos & 0x1F));
  return b._b128;
}

static __inline int equal_b128(b128 a, b128 b){
  b128 c = _mm_cmpeq_epi32(a,b);
  return get_int_0_b128(c)|get_int_1_b128(c)|get_int_2_b128(c)|get_int_3_b128(c);
}

static __inline b128 add_b128(b128 a, b128 b){ return _mm_add_epi32(a,b); }
static __inline b128 sub_b128(b128 a, b128 b){ return _mm_sub_epi32(a,b); }
static __inline b128 xor_b128(b128 a, b128 b){ return _mm_xor_si128(a,b); }
static __inline b128 or_b128 (b128 a, b128 b){ return _mm_or_si128(a,b); }
static __inline b128 and_b128(b128 a, b128 b){ return _mm_and_si128(a,b); }
static __inline b128 andnot_b128(b128 a, b128 b){ return _mm_andnot_si128(a,b); }
static __inline b128 negate_b128(b128 a){ return xor_b128(_mm_set1_epi8(0xff),a); }

#define shift_bytes_left_b128(__A, __IMM_BYTES)  _mm_slli_si128(__A,__IMM_BYTES)
#define shift_bytes_right_b128(__A, __IMM_BYTES) _mm_srli_si128(__A,__IMM_BYTES)

static __inline b128 shift_each32_bits_right_b128(b128 a, b128 b){ return _mm_srl_epi32(a,b); }
static __inline b128 shift_each32_bits_left_b128 (b128 a, b128 b){ return _mm_sll_epi32(a,b); }

b128 *alloc_b128(int num_b128);
b128 *calloc_b128(int num_b128);
void  free_b128(b128 *a);

static __inline b128 get_from_mem_b128(const b128 *ptr,int pos){ return _mm_load_si128(ptr+pos); }
static __inline b128 get_b128(const b128 *ptr){ return _mm_load_si128(ptr); }
static __inline void set_b128(b128 *ptr, b128 val){ _mm_store_si128(ptr,val); }
static __inline void set_b128_NOCACHEPOLLUTE(b128 *ptr, b128 val){ _mm_stream_si128(ptr,val); }

static __inline int get_bit_in_mem_b128(const b128 *ptr, int bitpos){
  return get_bit_b128(get_from_mem_b128(ptr,(bitpos>>7)), (bitpos & 0x7F));
}
static __inline void set_bit_in_mem_b128(b128 *ptr, int bitpos, int bitvalue){
  b128 *b128pos = ptr+(bitpos>>7);
  _mm_store_si128(b128pos, set_bit_b128(_mm_load_si128(b128pos),(bitpos&0x7F),bitvalue));
}

static __inline void mem_or_b128(b128 *dst, const b128 *src, const int n){
  for(int i=0;i<n;i++) dst[i]=or_b128(dst[i],src[i]);
}
static __inline void mem_and_b128(b128 *dst, const b128 *src, const int n){
  for(int i=0;i<n;i++) dst[i]=and_b128(dst[i],src[i]);
}
static __inline void mem_xor_b128(b128 *dst, const b128 *src, const int n){
  for(int i=0;i<n;i++) dst[i]=xor_b128(dst[i],src[i]);
}
static __inline void mem_andnot_b128(b128 *dst, const b128 *src, const int n){
  for(int i=0;i<n;i++) dst[i]=andnot_b128(src[i],dst[i]);
}
static __inline void mem_reversed_andnot_b128(b128 *dst, const b128 *src, const int n){
  for(int i=0;i<n;i++) dst[i]=andnot_b128(dst[i],src[i]);
}

void print_ints_b128(b128 a);
void print_bits_32(int a);
void print_bits_b128(b128 a);
void print_blocks_b128(b128 a, int block_size);

typedef unsigned long long ticks;
static __inline__ ticks getticks(void){
  unsigned a, d;
  asm volatile("rdtsc" : "=a" (a), "=d" (d));
  return ((ticks)a) | (((ticks)d) << 32);
}

#endif  /* ARM vs x86 */

#ifdef __cplusplus
}
#endif

#endif /* SSE2_WRAPPER_H */
