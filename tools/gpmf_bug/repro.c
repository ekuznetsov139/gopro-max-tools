/* Decodes one compressed GPMF payload: 20 x 50 followed by 20 x 53, uint8.
 * The bytes are the output of GPMFCompress() from gpmf-write for that input.
 * Build from a gpmf-parser checkout:
 *   gcc -I. -Idemo repro.c GPMF_parser.c GPMF_utils.c -lm -o repro
 */
#include <stdio.h>
#include <stdint.h>
#include "GPMF_parser.h"

static const uint8_t payload[] = {
    'D','E','V','C', 0, 1, 0, 32,               /* nest */
    'S','T','R','M', 0, 1, 0, 24,               /* nest */
    'D','I','S','P', '#', 1, 0, 16,             /* compressed, 16 bytes */
    0x42,0x01,0x00,0x28,                        /* uncompressed: 'B', size 1, repeat 40 */
    0x32,0x01,                                  /* first sample 50, quantize 1 */
    0xfc,0x36,0xe3,0xc5,0x00,0x00,0x00,0x00,0x00,0x00,  /* zero run of 19, +3, END */
};

int main(void)
{
    static uint32_t buf[16];
    static uint8_t out[64];
    GPMF_stream ms;
    uint32_t i, n;

    for (i = 0; i < sizeof(payload); i++)
        ((uint8_t *)buf)[i] = payload[i];
    if (GPMF_Init(&ms, buf, sizeof(payload)) != GPMF_OK ||
        GPMF_FindNext(&ms, STR2FOURCC("DISP"), GPMF_RECURSE_LEVELS) != GPMF_OK)
        return 1;
    n = GPMF_Repeat(&ms);
    if (GPMF_ScaledData(&ms, out, sizeof(out), 0, n, GPMF_TYPE_UNSIGNED_BYTE) != GPMF_OK)
        return 2;
    for (i = 0; i < n; i++)
        printf("%d%s", out[i], i + 1 < n ? " " : "\n");
    GPMF_Free(&ms);
    return 0;
}
