// stdin: rows of 640 uint8 (or argv[1] = count of a test row).  stdout: per row, u32 length + compressed KLV, from GoPro's compiled GPMFCompress.
#include <stdio.h>
#include <stdint.h>
#include <string.h>
uint32_t GPMFCompress(uint32_t* dst_gpmf, uint32_t *src_gpmf, uint32_t payloadAddition, uint32_t quantize);
int main(void){ static uint32_t src[256], dst[1024]; uint8_t* s=(uint8_t*)src; uint32_t h; static uint8_t row[4096];
 while(fread(&h,4,1,stdin)==1){ if(fread(row,1,h,stdin)!=h) return 2; memset(src,0,sizeof(src)); memset(dst,0,sizeof(dst));
  memcpy(s,"DISP",4); s[4]='B'; s[5]=1; s[6]=h>>8; s[7]=h&255; memcpy(s+8,row,h);
  uint32_t n=GPMFCompress(dst,src,8+((h+3)&~3),1); fwrite(&n,4,1,stdout); fwrite(dst,1,n,stdout); } }
