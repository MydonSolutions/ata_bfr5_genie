antenna_names = [
    "1a",
    "1b",
    "1c",
    "1d",
    "1e",
    "1f",
    "1g",
    "1h",
    "1j",
    "1k",
    "2a",
    "2b",
    "2c",
    "2d",
    "2e",
    "2f",
    "2g",
    "2h",
    "2j",
    "2k",
    "2l",
    "2m",
    "3c",
    "3d",
    "3e",
    "3f",
    "3g",
    "3h",
    "3j",
    "3l",
    "4e",
    "4f",
    "4g",
    "4h",
    "4j",
    "4k",
    "4l",
    "5b",
    "5c",
    "5e",
    "5g",
    "5h",
    "0a",
]

NANTS = length(antenna_names)
NCHAN = 2048
NPOL = 2

weights = ones(ComplexF64, (NPOL, NCHAN, NANTS))

open("ant_weights.bin", "w") do fio
    write(fio,
        htol(Int32(NANTS)),
        htol(Int32(NCHAN)),
        htol(Int32(NPOL))
    )
    write(fio, join(antenna_names, '\0')*'\0')
    write(fio, htol.(weights))
end
