#include "numchol.h"
#include <string.h>

void ShutDown(void){

  /*   sdpdat* sdt = sdat; */
  dsdp_printf("\n Shutdown --  ");

} /* ShutDown */


int ExitProc(int  ccode,
             char *str)
{
  xcode code=(xcode)ccode;

  dsdp_printf("\n Exit -- %d : ",ccode);
  
  switch (code) {
    case OptFound:
      dsdp_printf("optimal solution found");
      return code;
    case OutOfSpc:
      dsdp_printf("out of memory space");
      break;
    default:
      break;
  }
  if (str) dsdp_printf(", %s",str);

  ShutDown();

  dsdp_printf("\n Exiting --  ");

  return 1;
} /* ExitProc */

