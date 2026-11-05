Supporting files for the `GPMF_Decompress` report (see docs/FORMAT.md, "A decoder bug").

- `repro.c`: decodes one small payload, printing what the decoder returns.
- `fix.patch`: the suggested change to gpmf-parser's `GPMF_parser.c`.
- `klv_decode.c`: gpmf-parser's decoder as a filter, for comparing decodings. stdin is repeated `[u32 nbytes][KLV bytes]`, stdout one row per line.
- `enc.c`: gpmf-write's `GPMFCompress` as a filter, for re-encoding rows and comparing with a file's bytes.

Build each against a checkout of gpmf-parser or gpmf-write, for example:

    gcc -I$P -I$P/demo repro.c $P/GPMF_parser.c $P/GPMF_utils.c -lm -o repro
    gcc -I$W enc.c $W/GPMF_writer.c -lpthread -lm -o enc
