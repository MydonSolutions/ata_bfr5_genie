
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