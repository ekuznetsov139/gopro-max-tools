// Latest upstream GPMF_parser.c on raw DISP KLVs.  stdin: repeated [u32 nbytes][KLV bytes]; stdout: one row per line.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "GPMF_parser.h"
int main(void){ uint32_t n; static uint32_t buf[4096]; static uint8_t out[8192]; int bad=0;
 while(fread(&n,4,1,stdin)==1){ if(fread(buf,1,n,stdin)!=n) return 2; GPMF_stream ms; 
  if(GPMF_Init(&ms,buf,n)!=GPMF_OK){bad++;continue;} 
  if(GPMF_FindNext(&ms,STR2FOURCC("DISP"),GPMF_RECURSE_LEVELS)!=GPMF_OK){bad++;GPMF_Free(&ms);continue;}
  uint32_t reps=GPMF_Repeat(&ms), el=GPMF_ElementsInStruct(&ms);
  if(GPMF_ScaledData(&ms,out,sizeof(out),0,reps,GPMF_TYPE_UNSIGNED_BYTE)!=GPMF_OK){bad++;GPMF_Free(&ms);continue;}
  for(uint32_t i=0;i<reps*el;i++) printf("%d ",out[i]); printf("\n"); GPMF_Free(&ms);
 } fprintf(stderr,"undecoded %d\n",bad); }
