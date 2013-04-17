#ifndef RTERR
#define RTERR
#include <string.h>

const int WHAT_MAX = 200;
class rtErr {
  
 public:
  char what[WHAT_MAX];

  rtErr(char *_what);
  ~rtErr();

};
#endif
