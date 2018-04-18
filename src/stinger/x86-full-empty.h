/* x86 Emulation of Full Empty Bits Using atomic check-and-swap.
 *
 * NOTES:
 * - Using these functions means that the MARKER value defined 
 *   below must be reserved in your application and CANNOT be 
 *   considered a normal value.  Feel free to change the value to
 *   suit your application.
 * - Do not compile this file with optimization.  There are loops
 *   that do not make sense in a serial context.  The compiler will
 *   optimize them out.
 * - Improper use of these functions can and will result in deadlock.
 *
 * author: rmccoll3@gatech.edu
 */

#ifndef	X86_FULL_EMPTY_C
#define	X86_FULL_EMPTY_C

#include    "stinger-atomics.h"

#if !defined(__MTA__) /* x86 only */

#include  <stdint.h>
#define MARKER INT64_MAX

int64_t 
readfe(int64_t * v);

int64_t
writeef(int64_t * v, int64_t new_val);

int64_t
readff(int64_t * v);

int64_t
writeff(int64_t * v, int64_t new_val);

int64_t
writexf(int64_t * v, int64_t new_val);

#endif  /* x86 only */

#endif  /*X86-FULL-EMPTY_C*/

