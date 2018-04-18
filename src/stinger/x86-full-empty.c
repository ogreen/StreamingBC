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

#include  "x86-full-empty.h"

#if !defined(__MTA__) /* x86 only */

int64_t 
readfe(int64_t * v) {
  int64_t val;
  while(1) {
    val = *v;
    while(val == MARKER) {
      val = *v;
    }
    if(val == stinger_int64_cas(v, val, MARKER))
      break;
  }
  return val;
}

int64_t
writeef(int64_t * v, int64_t new_val) {
  int64_t val;
  while(1) {
    val = *v;
    while(val != MARKER) {
      val = *v;
    }
    if(MARKER == stinger_int64_cas(v, MARKER, new_val))
      break;
  }
  return val;
}

int64_t
readff(int64_t * v) {
  int64_t val = *v;
  while(val == MARKER) {
    val = *v;
  }
  return val;
}

int64_t
writeff(int64_t * v, int64_t new_val) {
  int64_t val;
  while(1) {
    val = *v;
    while(val == MARKER) {
      val = *v;
    }
    if(val == stinger_int64_cas(v, val, new_val))
      break;
  }
  return val;
}

int64_t
writexf(int64_t * v, int64_t new_val) {
  *v = new_val;
  return new_val;
}

#endif  /* x86 only */

