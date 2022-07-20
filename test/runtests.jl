using ATA_BFR5_Genie
using Test

@testset "ATA_BFR5_Genie.jl" begin
		
	using Blio: GuppiRaw
	using BeamformerRecipes
	using TOML

	include("./generate_antenna_weights.jl")
	include("./generate_guppiraw.jl")

	fio = open("header.0000.raw", "r")
		header = GuppiRaw.Header(28)
		@assert read!(fio, header)

		diminfo, beaminfo, obsinfo = collectDimBeamObsInfo(header)
		obs_antnames = collectObsAntnames(header)
		obs_chanrange = header["SCHAN"]:(header["SCHAN"] + diminfo.nchan)

		delayinfo = DelayInfo()
		delayinfo.time_array = []
		delayinfo.dut1 = header["DUT1"]
		
		push!(delayinfo.time_array, calculateMJDfromEpoch(calculateEpochGuppiHeader(header, 0.5)))
		
		while read!(fio, header)
			push!(delayinfo.time_array, calculateMJDfromEpoch(calculateEpochGuppiHeader(header, 0.5)))
		end
		
		ntimes = length(delayinfo.time_array)
		delayinfo.delays = zeros(Float64, (diminfo.nants, diminfo.nbeams, ntimes))
		delayinfo.rates = zeros(Float64, (diminfo.nants, diminfo.nbeams, ntimes))
		delayinfo.jds = zeros(Float64, (ntimes))
	close(fio)

	fio = open("ant_weights.bin", "r")
		antcal_weights = collectAntennaWeights(fio, obs_antnames, obs_chanrange)
		# display(antcal_weights):println()
	close(fio)

	calinfo = CalInfo()
	calinfo.refant = obs_antnames[1]
	calinfo.cal_K = zeros(Float32, (diminfo.nants, diminfo.npol))
	calinfo.cal_G = ones(ComplexF32, (diminfo.nants, diminfo.npol))
	calinfo.cal_B = antcal_weights
	calinfo.cal_all = antcal_weights

	fio = open("ata_telinfo.toml", "r")
		telinfo = collectTelinfo(TOML.parse(fio), obs_antnames)
		antenna_positions_xyz = calculateAntennaXyzPositions(telinfo)
		# display(telinfo);println()
	close(fio)


	for (i, midblock_time_unix) in enumerate(delayinfo.time_array)
		delayinfo.delays[:, :, i] = calculateBeamDelays(
			antenna_positions_xyz, 1,
			obsinfo.phase_center_ra, obsinfo.phase_center_dec,
			transpose(hcat(beaminfo.ras, beaminfo.decs)),
			telinfo.longitude, telinfo.latitude, telinfo.altitude,
			midblock_time_unix, delayinfo.dut1
		)
	end

	recipe = BeamformerRecipe(
		diminfo,
		telinfo,
		obsinfo,
		calinfo,
		beaminfo,
		delayinfo
	)

	to_hdf5("out.bfr5", recipe)
	@test true
end

