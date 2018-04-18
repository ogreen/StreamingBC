#include <stdint.h>

#define MAX_VTXID 0xFFFFF
/** \def MAX_VTXID
*   \brief Maximum number of vertices the physical mapper can support.
*   Note: If the physical mapper produces errors, increase this number.
*/
#define MAX_NODES 0xFFFFF
/** \def MAX_NODES
*   \brief Maximum number internal nodes that the physical mapper can support.
*   Note: If the physical mapper produces errors, increase this number.
*/

typedef struct stinger_physmap stinger_physmap_t;

stinger_physmap_t * 
stinger_physmap_create();

void
stinger_physmap_delete(stinger_physmap_t * map);

int64_t
stinger_physmap_create_mapping (stinger_physmap_t * map, char * string, int64_t length);

int64_t
stinger_physmap_get_mapping (stinger_physmap_t * map, char * string, int64_t length);

int
stinger_physmap_get_key (stinger_physmap_t * map, char ** outbuffer, int64_t * outbufferlength, int64_t vertexID);

