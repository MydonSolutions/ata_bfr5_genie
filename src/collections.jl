
# Not only filter on antenna-names, but also return `cal` appropriate [nants, npoln, nchan]
function collectAntennaWeights(
	io::IO,
	ant_name_filter::Array{String, 1},
	channel_range::Union{UnitRange, Colon}
)::Array{ComplexF64, 3}
	antweights = AntennaWeights(io)
	ant_indices = [findfirst(x -> x == name, antweights.names) for name in ant_name_filter]

	return permutedims(
		antweights.weights[:, channel_range, ant_indices],
		[3, 1, 2]
	)
end

function collectTelinfo(telinfoDict::Dict)::TelInfo
	collectTelinfo(telinfoDict, Vector{String}())
end

function collectTelinfo(telinfoDict::Dict, ant_name_filter::Array{String, 1})::TelInfo
	telinfo = TelInfo()
	default_diameter = get(telinfoDict, "antenna_diameter", 0.0)

	if length(ant_name_filter) == 0
		for antinfo in values(telinfoDict["antennas"])
			push!(telinfo.antenna_names, antinfo["name"])
			push!(telinfo.antenna_numbers, antinfo["number"])
			push!(telinfo.antenna_diameters, "diameter" in keys(antinfo) ? antinfo["diameter"] : default_diameter)
			telinfo.antenna_positions = cat(telinfo.antenna_positions, antinfo["position"], dims=2)
		end
	else
		name_info_map = Dict{String, Any}(
			map(antinfo -> antinfo["name"] => antinfo, values(telinfoDict["antennas"]))
		)
		for antname in ant_name_filter
			antinfo = name_info_map[antname]
			push!(telinfo.antenna_names, antinfo["name"])
			push!(telinfo.antenna_numbers, antinfo["number"])
			push!(telinfo.antenna_diameters, "diameter" in keys(antinfo) ? antinfo["diameter"] : default_diameter)
			telinfo.antenna_positions = cat(telinfo.antenna_positions, antinfo["position"], dims=2)
		end
	end

	telinfo.antenna_position_frame = telinfoDict["antenna_position_frame"]
	telinfo.latitude = isa(telinfoDict["latitude"], String) ? dms2deg(telinfoDict["latitude"]) : telinfoDict["latitude"]
	telinfo.longitude = isa(telinfoDict["longitude"], String) ? dms2deg(telinfoDict["longitude"]) : telinfoDict["longitude"]
	telinfo.altitude = telinfoDict["altitude"]
	telinfo.telescope_name = telinfoDict["telescope_name"]
	return telinfo
end

function collectDimBeamObsInfo(header::GuppiRaw.Header)
	# infer number of beams
	beaminfo = BeamInfo()
	nbeams = 0
	while @sprintf("RA_OFF%01d",nbeams) in keys(header)
		push!(beaminfo.src_names, @sprintf("BEAM_%01d", nbeams))
		push!(beaminfo.ras, (pi/180.0) * header[@sprintf("RA_OFF%01d", nbeams)] * 360.0 / 24.0) # convert hours to degrees
		push!(beaminfo.decs, (pi/180.0) * header[@sprintf("DEC_OFF%01d", nbeams)])
		nbeams += 1
	end
	if nbeams == 0 
		push!(beaminfo.src_names, "BEAM_BORESIGHT")
		push!(beaminfo.ras, (pi/180.0) * header["RA_STR"] * 360.0 / 24.0) # convert hours to degrees
		push!(beaminfo.decs, (pi/180.0) * header["DEC_STR"])
		nbeams = 1
	end

	diminfo = DimInfo()
	diminfo.npol, diminfo.ntimes, diminfo.nchan, diminfo.nants = size(header)
	diminfo.nbeams = nbeams

	obsinfo = ObsInfo()
	obsinfo.obsid = header["SRC_NAME"]

	# OBSFREQ is the center frequency of the observed data
	obsinfo.freq_array = [header["OBSFREQ"] - (chan - diminfo.nchan/2)*header["CHAN_BW"] for chan in 0:diminfo.nchan-1]
	obsinfo.freq_array /= 1e3 # MHz -> GHz

	obsinfo.phase_center_ra = header["RA_STR"] * 360.0 / 24.0 # convert hours to degrees
	obsinfo.phase_center_dec = header["DEC_STR"]
	obsinfo.phase_center_ra *= pi/180.0
	obsinfo.phase_center_dec *= pi/180.0
	obsinfo.instrument_name = get(header, "TELESCOP", "Unknown")

	return diminfo, beaminfo, obsinfo
end

function collectObsAntnames(header::GuppiRaw.Header)
	obs_antenna_names = Array{String, 1}()
	antnmes_index = 0
	while @sprintf("ANTNMS%02d", antnmes_index) in keys(header)
		append!(obs_antenna_names, split(header[@sprintf("ANTNMS%02d", antnmes_index)], ','))
		antnmes_index += 1
	end

	return map(antlo -> antlo[1:end-1], obs_antenna_names)
end