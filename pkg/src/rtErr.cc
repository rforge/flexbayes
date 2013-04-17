#include "rtErr.h"

rtErr::rtErr(char* _what) {
  strncpy(what, _what, WHAT_MAX);
}

rtErr::~rtErr() {};
