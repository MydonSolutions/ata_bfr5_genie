using Printf
using ATA_BFR5_Genie
using BeamformerRecipes

rawpath = ARGS[1]
antweightpath = ARGS[2]
telinfopath = ARGS[3]

recipe = collectBfr5(rawpath, antweightpath, telinfopath)

bfr5path = @sprintf("%s.bfr5", rawpath)
to_hdf5(bfr5path, recipe)
println(bfr5path)