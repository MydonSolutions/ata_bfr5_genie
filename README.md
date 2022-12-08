# BFR5 Generation for the ATA

The SETI operations at the ATA have their beamforming information split across 3 different file-sources. These sources will be scraped for the information required to populate a beamforming-recipe (BFR5) file.

## Typical uses

The simplest use of the toplevel `bfr5_gen.jl` script omits the antenna calibration file (populating /calinfo/cal_all with 1+0j instead) and scrapes beam-coordinates from the RA/DEC_OFF%01d key-values of the RAW file's first block:

`julia ./bfr5_gen.jl /path/to/rawfile.0000.raw ./test/ata_telinfo.toml`

Specification of beam coordinates can be made at the CLI too:
`julia ./bfr5_gen.jl /path/to/rawfile.0000.raw ./test/ata_telinfo.toml --beam 3.141,1.618 --beam 3.040,1.608`

Furthermore, consult the --help usage text.