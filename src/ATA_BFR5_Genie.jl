module ATA_BFR5_Genie

using Printf: @sprintf

using BeamformerRecipes
using RadioInterferometry: dms2deg, xyz2uvw
using Blio: GuppiRaw
using ERFA: apco13, atciq, atioq
using Geodesy: ECEFfromLLA, wgs84_ellipsoid, LLA
using Dates: unix2datetime

export 
	collectAntennaWeights,
	AntennaWeights,
	collectTelinfo,
	collectDimBeamObsInfo,
	collectObsAntnames,
	calculateBeamDelays,
	calculateEpochGuppiHeader,
	calculateMJDfromEpoch,
	calculateAntennaXyzPositions

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

# Not only filter on antenna-names, but also return `cal` appropriate [nants, npoln, nchan]
function collectAntennaWeights(
	io::IO,
	ant_name_filter::Array{String, 1},
	channel_range::UnitRange
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

	for antinfo in values(telinfoDict["antennas"])
		if length(ant_name_filter) == 0 || antinfo["name"] in ant_name_filter
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

function calculateAntennaXyzPositions(telinfo::TelInfo)	
	if telinfo.antenna_position_frame == "ecef"
		lla = LLA(telinfo.longitude, telinfo.latitude, telinfo.altitude)
		center_ecef = ECEFfromLLA(wgs84_ellipsoid)(lla)
		center_xyz = [center_ecef.x, center_ecef.y, center_ecef.z]

		return telinfo.antenna_positions .- center_xyz
	elseif telinfo.antenna_position_frame == "enu"
		return enu2xyz(telinfo.antenna_positions)
	end
	return telinfo.antenna_positions
end

function collectDimBeamObsInfo(header::GuppiRaw.Header)
	# infer number of beams
	beaminfo = BeamInfo()
	nbeams = 0
	while @sprintf("RA_OFF%d",nbeams) in keys(header)
		push!(beaminfo.src_names, @sprintf("BEAM_%d", nbeams))
		push!(beaminfo.ras, header[@sprintf("RA_OFF%d", nbeams)])
		push!(beaminfo.decs, header[@sprintf("DEC_OFF%d", nbeams)])
		nbeams += 1
	end

	diminfo = DimInfo()
	diminfo.npol, diminfo.ntimes, diminfo.nchan, diminfo.nants = size(header)
	diminfo.nbeams = nbeams

	
	obsinfo = ObsInfo()
	obsinfo.obsid = header["SRC_NAME"]
	obsinfo.freq_array = [(header["SCHAN"] + 0.5 + chan)*header["CHAN_BW"] for chan in 0:diminfo.nchan-1]
	obsinfo.phase_center_ra = header["RA_STR"]
	obsinfo.phase_center_dec = header["DEC_STR"]
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

function calculateEpochGuppiHeader(header::GuppiRaw.Header, piperblk_offset_factor::Real = 0.5)::Real
	ntime = size(header, 2)
	header["SYNCTIME"] + (header["PKTIDX"] + header["PIPERBLK"]*piperblk_offset_factor) * ntime/(header["PIPERBLK"] * 1e6 * header["CHAN_BW"])
end

function calculateJDfromEpoch(unix_epoch_seconds::Real)::Real
	(unix_epoch_seconds / 86400) + 2440587.5
end

function calculateEpochfromJD(jd::Real)::Real
	(jd - 2440587.5) * 86400
end

function calculateMJDfromJD(jd::Real)::Real
	jd + 2400000.5
end

function calculateJDfromMJD(mjd::Real)::Real
	mjd - 2400000.5
end

function calculateMJDfromEpoch(unix_epoch_seconds::Real)::Real
	calculateMJDfromJD(calculateJDfromEpoch(unix_epoch_seconds))
end

function calculatePositionDelays(
	positions_xyz::AbstractArray{<:Real, 2},
	reference_position_index::Integer,
	ha_rad::Real, dec_rad::Real, lon_rad::Real
)

	positions_uvw = xyz2uvw(positions_xyz, ha_rad, dec_rad, lon_rad)
	reference_w = positions_uvw[3, reference_position_index]
	
	return (positions_uvw[3, :] .- reference_w) ./ 299792458.0
end

function calculateHaDec(
	ra_rad::Real, dec_rad::Real, astrom
)
	ri, di = atciq(
		ra_rad, dec_rad,
		0, 0, 0, 0,
		astrom
	)
	aob, zob, hour_angle_rad, declination_rad, rob = atioq(
		ri, di,
		astrom
	)
	return hour_angle_rad, declination_rad
end

function calculateBeamDelays(
	antenna_positions_xyz::AbstractArray{<:Real, 2},
	reference_position_index::Integer,
	boresight_ra::Real, boresight_dec::Real,
	beams_radec::AbstractArray{<:Real, 2},
	longitude_rad::Real, latitude_rad::Real, altitude::Real,
	timemjd::Real, dut1::Real
)
	println(unix2datetime(calculateEpochfromJD(calculateJDfromMJD(timemjd))))
	astrom, eo = apco13(
		timemjd, 0,
		dut1,
		longitude_rad, latitude_rad, altitude,
		0, 0,
		0, 0, 0, 0
	)

	boresight_ha, boresight_dec = calculateHaDec(
		boresight_ra, boresight_dec, astrom
	)

	boresight_delays = calculatePositionDelays(
		antenna_positions_xyz, reference_position_index,
		boresight_ha, boresight_dec, longitude_rad
	)
	
	beam_antenna_delays = zeros(
		Float64,
		(
			size(antenna_positions_xyz, 2),
			size(beams_radec, 2)
		)
	)

	for beam_index in 1:size(beams_radec, 2)
		beam_ha, beam_dec = calculateHaDec(
			beams_radec[1, beam_index], beams_radec[2, beam_index], astrom
		)

		beam_delays = calculatePositionDelays(
			antenna_positions_xyz, reference_position_index,
			beam_ha, beam_dec, longitude_rad
		)

		beam_antenna_delays[:, beam_index] = beam_delays .- boresight_delays
	end

	return beam_antenna_delays
end

end # ATA_BFR5_Genie module
