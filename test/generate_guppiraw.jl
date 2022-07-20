using Blio: GuppiRaw

header = GuppiRaw.Header()
header["CHAN_BW"]		= 0.5
header["OBSNCHAN"]	= 128
header["NANTS"]			= 4
header["NPOL"]			= 2
header["NBITS"]			= 8
header["BLOCSIZE"]	= 128*4*2*2
header["ANTNMS00"]	= "1cB,1eB"
header["ANTNMS01"]	= "1gB,1hB"
header["DEC_STR"]		= 21.9
header["RA_STR"]		= 19.9
header["DEC_OFF0"]	= 21.0
header["RA_OFF0"]		= 19.0
header["DEC_OFF1"]	= 21.1
header["RA_OFF1"]		= 19.1
header["DEC_OFF2"]	= 21.2
header["RA_OFF2"]		= 19.2
header["DIRECTIO"]	= 1
header["DUT1"]			= -0.055208199999999999
header["PIPERBLK"]	= 8192
header["PKTIDX"]		= 1364518798816
header["SCHAN"]			= 928
header["SRC_NAME"]	= "J1934+2153"
header["STTVALID"]	= 1
header["STT_IMJD"]	= 59775
header["STT_OFFS"]	= 0
header["STT_SMJD"]	= 33754
header["SYNCTIME"]	= 1655150143


open("header.0000.raw", "w") do fio
	for i in 1:20
		write(fio, header)
		header["PKTIDX"] += header["PIPERBLK"]
	end
end