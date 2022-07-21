using ATA_BFR5_Genie
using Test

@testset "ATA_BFR5_Genie.jl" begin
	include("./generate_antenna_weights.jl")
	include("./generate_guppiraw.jl")

	recipe = collectBfr5("header.0000.raw", "ant_weights.bin", "ata_telinfo.toml", headerentry_limit=28, headers_only=true)
	@test true
end

