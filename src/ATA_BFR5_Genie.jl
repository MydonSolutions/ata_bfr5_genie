module ATA_BFR5_Genie

using Printf: @sprintf

using BeamformerRecipes
using RadioInterferometry: dms2deg, xyz2uvw
using Blio: GuppiRaw
using ERFA: apco13, atciq, atioq
using Geodesy: ECEFfromLLA, wgs84_ellipsoid, LLA
using Dates: unix2datetime
using TOML

export bfr5Collect, main_bfr5Gen

using ArgParse
function main_bfr5Gen()::Cint
  s = ArgParseSettings()
  @add_arg_table s begin
    "rawpath"
      help = "path to the RAW file"
      required = true
    "telinfopath"
      help = "path to the telescope information TOML file"
      required = true
    "outputpath"
      help = "path of the output BFR5 file"
      default = "./output.bfr5"
      required = false
    "--antweightpath"
      help = "path to the antenna-weights file, otherwise all coefficients are 1+0j"
      default = nothing
      required = false
    "--beam"
      help = "beam coordinate as 'RA,DEC' (hours,degrees), repeatable"
      arg_type = String
      required = false
      action = :append_arg
      default = nothing
    "--phase-center"
      help = "phase-center coordinate as 'RA,DEC' (hours,degrees)"
      arg_type = String
      required = false
      default = nothing
  end

  args = parse_args(s)

  recipe = bfr5Collect(
    args["rawpath"],
    args["telinfopath"];
    antweights_filepath = args["antweightpath"],
    beam_coords = args["beam"],
    phase_center = args["phase-center"],
  )

  to_hdf5(args["outputpath"], recipe)
  println(args["outputpath"])

  return 0
end

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

function bfr5Collect(
	guppiraw_filepath::String,
	telinfo_filepath::String;
	headerentry_limit::Integer=256,
	headers_only::Bool=true,
	antweights_filepath::Union{String, Nothing}=nothing,
	beam_coords::Union{Vector{String}, Nothing}=nothing, # RAh,DECdeg strings
	phase_center::Union{String, Nothing}=nothing # RAh,DECdeg string
)::BeamformerRecipe

	guppiraw_stem_match = match(r"(.*)\.(\d{4})\.raw", guppiraw_filepath)
	guppiraw_stempath, guppiraw_stem_index = guppiraw_stem_match[1], parse(Int, guppiraw_stem_match[2])

	fio = open(guppiraw_filepath, "r")
		header = GuppiRaw.Header(headerentry_limit)
		@assert read!(fio, header, skip_padding=headers_only)

		diminfo, beaminfo, obsinfo = collectDimBeamObsInfo(header)
		schan = header["SCHAN"]
		if !isnothing(beam_coords) && length(beam_coords) > 0
			beams = hcat(collect(
				parse.(Float64, split(beam, ","))
				for (i, beam) in enumerate(beam_coords)
			)...)
			diminfo.nbeams = size(beams)[2]
			beaminfo.ras = beams[1, :] .* ((360.0 / 24.0) * (pi/180.0))
			beaminfo.decs = beams[2, :] .* (pi/180.0)
			beaminfo.src_names = collect(@sprintf("BEAM_%01d", i) for i in 0:diminfo.nbeams-1)
		end
		if !isnothing(phase_center)
			coords = parse.(Float64, split(phase_center, ","))
			obsinfo.phase_center_ra = coords[1] * ((360.0 / 24.0) * (pi/180.0))
			obsinfo.phase_center_dec = coords[2] * (pi/180.0)
		end
		obs_antnames = collectObsAntnames(header)

		delayinfo = DelayInfo()
		delayinfo.time_array = []
		delayinfo.dut1 = header["DUT1"]
		
		push!(delayinfo.time_array, calculateEpochGuppiHeader(header, 0.5))
		
		while true 
			if headers_only
				skip(fio, header["BLOCSIZE"])
			end
			try
				if ! read!(fio, header, skip_padding=headers_only)
					guppiraw_stem_index += 1
					next_guppiraw_filepath = @sprintf("%s.%04d.raw", guppiraw_stempath, guppiraw_stem_index)
					if !isfile(next_guppiraw_filepath)
						break
					end

					close(fio)
					println("Covering", next_guppiraw_filepath)
					guppiraw_filepath = next_guppiraw_filepath
					fio = open(guppiraw_filepath, "r")
					@assert read!(fio, header, skip_padding=headers_only)

				end
			catch err
				println("`read!(..., skip_padding=", headers_only, ")` error caught at block ", length(delayinfo.time_array), ": ", err)
				break
			end
			push!(delayinfo.time_array, calculateEpochGuppiHeader(header, 0.5))
		end
		delayinfo.jds = map(calculateJDfromEpoch, delayinfo.time_array)
		
		ntimes = length(delayinfo.time_array)
		delayinfo.delays = zeros(Float64, (diminfo.nants, diminfo.nbeams, ntimes))
		delayinfo.rates = zeros(Float64, (diminfo.nants, diminfo.nbeams, ntimes))
	close(fio)

	if isnothing(antweights_filepath)
		antcal_weights = ones(ComplexF64, (diminfo.nants, diminfo.npol, schan+diminfo.nchan))
	else
		fio = open(antweights_filepath, "r")
			antcal_weights = collectAntennaWeights(fio, obs_antnames, :)
			# display(antcal_weights):println()
		close(fio)
	end

	calinfo = CalInfo()
	calinfo.refant = obs_antnames[1]
	calinfo.cal_K = zeros(Float32, (diminfo.nants, diminfo.npol))
	calinfo.cal_G = ones(ComplexF32, (diminfo.nants, diminfo.npol))
	calinfo.cal_B = antcal_weights
	calinfo.cal_all = antcal_weights
	diminfo.nchan = size(antcal_weights)[3]

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
		) .* 1e9
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
