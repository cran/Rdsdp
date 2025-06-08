#ifndef DSDP_IO_H
#define DSDP_IO_H

// #define USING_R 1

#ifdef USING_R
  #include <R.h>
  #define dsdp_printf(...) Rprintf(__VA_ARGS__)
  #define dsdp_sprintf(buf, ...) snprintf((buf), sizeof(buf), __VA_ARGS__)
  #define dsdp_vsprintf(buf, fmt, args)  vsnprintf((buf), sizeof(buf), (fmt), (args)) 
  #define dsdp_exit(code)         error("DSDP terminated with exit code %d", (code))
  #define dsdp_error(msg)  Rf_error(msg)
#else
  #include <stdio.h>
  #include <stdlib.h>
  #define dsdp_printf(...) printf(__VA_ARGS__)
  #define dsdp_sprintf sprintf
  #define dsdp_vsprintf vsprintf
  #define dsdp_exit(code) exit(code)
  #define dsdp_error(msg)  do { fprintf(stderr, "%s\n", msg); exit(1); } while(0)
#endif

#endif
