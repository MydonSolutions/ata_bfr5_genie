module ATA_BFR5_Genie

using Printf: @sprintf

using BeamformerRecipes
using RadioInterferometry: dms2deg, xyz2uvw
using Blio: GuppiRaw
using ERFA: apco13, atciq, atioq
using Geodesy: ECEFfromLLA, wgs84_ellipsoid, LLA
using Dates: unix2datetime
using TOML

export collectBfr5

struct AntennaWeights
	nants::Int32
	nchan::Int32
	npoln::Int32
	names::Array{String, 1}
	weights::Array{ComplexF64, 3}
end

function AntennaWeights(io::IO)::AntennaWeights
	nants = ltoh(read(io, Int32))
	nchan = ltoh(read(io, Int32))
	npoln = ltoh(read(io, Int32))
	return AntennaWeights(
		nants,
		nchan,
		npoln,
		[readuntil(io, '\0') for i in 1:nants],
		reshape(
			collect(
				ltoh(read(io, ComplexF64)) for i in 1:nants*nchan*npoln
			),
			(npoln, nchan, nants)
		)
	)
end

include("collections.jl")
include("calculations.jl")

function collectBfr5(
	guppiraw_filepath::String,
	antweights_filepath::String,
	telinfo_filepath::String;
	headerentry_limit::Integer=256,
	headers_only::Bool=false
)::BeamformerRecipe

	fio = open(guppiraw_filepath, "r")
		header = GuppiRaw.Header(headerentry_limit)
		@assert read!(fio, header, skip_padding=!headers_only)

		diminfo, beaminfo, obsinfo = collectDimBeamObsInfo(header)
		obs_antnames = collectObsAntnames(header)
		obs_chanrange = header["SCHAN"]:(header["SCHAN"] + diminfo.nchan)

		delayinfo = DelayInfo()
		delayinfo.time_array = []
		delayinfo.dut1 = header["DUT1"]
		
		push!(delayinfo.time_array, calculateEpochGuppiHeader(header, 0.5))
		
		while true 
			if !headers_only
				skip(fio, header["BLOCSIZE"])
			end
			try
				if ! read!(fio, header, skip_padding=!headers_only)
					break
				end
			catch err
				println("`read!(..., skip_padding=", skip_padding, ")` error caught at block ", length(delayinfo.time_array), ": ", err)
				break
			end
			push!(delayinfo.time_array, calculateEpochGuppiHeader(header, 0.5))
		end
		delayinfo.jds = map(calculateJDfromEpoch, delayinfo.time_array)
		
		ntimes = length(delayinfo.time_array)
		delayinfo.delays = zeros(Float64, (diminfo.nants, diminfo.nbeams, ntimes))
		delayinfo.rates = zeros(Float64, (diminfo.nants, diminfo.nbeams, ntimes))
	close(fio)

	fio = open(antweights_filepath, "r")
		antcal_weights = collectAntennaWeights(fio, obs_antnames, obs_chanrange)
		# display(antcal_weights):println()
	close(fio)

	calinfo = CalInfo()
	calinfo.refant = obs_antnames[1]
	calinfo.cal_K = zeros(Float32, (diminfo.nants, diminfo.npol))
	calinfo.cal_G = ones(ComplexF32, (diminfo.nants, diminfo.npol))
	calinfo.cal_B = antcal_weights
	calinfo.cal_all = antcal_weights

	fio = open(telinfo_filepath, "r")
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

	return BeamformerRecipe(
		diminfo,
		telinfo,
		obsinfo,
		calinfo,
		beaminfo,
		delayinfo
	)
end

end # ATA_BFR5_Genie module
