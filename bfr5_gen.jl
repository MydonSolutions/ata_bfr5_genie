using Printf
using ATA_BFR5_Genie
using BeamformerRecipes
using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
  "rawpath"
    help = "path to the raw file"
    required = true
  "telinfopath"
    help = "path to the telescope information TOML file"
    required = true
  "outputpath"
    help = "path of the output BFR5 file"
    default = "./output.bfr5"
    required = false
  "--antweightpath"
    help = "path to the antenna-weights file"
    default = nothing
    required = false
  "--beam"
    help = "beam coordinate as 'RA,DEC' (hours,degrees)"
    arg_type = String
    required = false
    action = :append_arg
    default = nothing
end

args = parse_args(s)

recipe = collectBfr5(
  args["rawpath"],
  args["telinfopath"];
  antweights_filepath = args["antweightpath"],
  beam_coords = args["beam"]
)

to_hdf5(args["outputpath"], recipe)
println(args["outputpath"])