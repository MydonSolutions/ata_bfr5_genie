using ATA_BFR5_Genie
using Test

@testset "ATA_BFR5_Genie.jl" begin
	include("./generate_antenna_weights.jl")
	include("./generate_guppiraw.jl")

	recipe = collectBfr5("header.0000.raw", "ant_weights.bin", "ata_telinfo.toml")
	@test true
end

