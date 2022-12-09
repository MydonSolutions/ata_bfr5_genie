# BFR5 Generation for the ATA

The SETI operations at the ATA have their beamforming information split across 3 different file-sources. These sources will be scraped for the information required to populate a beamforming-recipe (BFR5) file.

## Installation and Toplevel use

`julia -e 'using Pkg; Pkg.add(url="https://github.com/MydonSolutions/ata_bfr5_genie")'`

An internal function parses command-line arguments and can be invoked via: `julia -e 'using ATA_BFR5_Genie; main_bfr5Gen()' -- -h`. Alternatively, this can be placed in a script and the script provided to the Julia runtime:


```lang=julia
using ATA_BFR5_Genie
main_bfr5Gen()
```

### Typical uses

The simplest use of the toplevel `main_bfr5Gen` function omits the antenna calibration file (populating /calinfo/cal_all with 1+0j instead) and scrapes beam-coordinates from the RA/DEC_OFF%01d key-values of the RAW file's first block:

`julia -e 'using ATA_BFR5_Genie; main_bfr5Gen()' -- /path/to/rawfile.0000.raw ./test/ata_telinfo.toml`

Specification of beam coordinates can be made at the CLI too:
`julia -e 'using ATA_BFR5_Genie; main_bfr5Gen()' -- /path/to/rawfile.0000.raw ./test/ata_telinfo.toml --beam 3.141,1.618 --beam 3.040,1.608`

The optional antenna calibration coefficients file is specified [here](https://github.com/MydonSolutions/ata_antenna_weights_binary).

Furthermore, consult the `--help` text.