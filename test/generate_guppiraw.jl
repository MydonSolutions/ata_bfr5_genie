using Blio: GuppiRaw

n_bits = 8
n_pols = 2
n_time = 1024
n_chan_perant = 128
n_ant = 4

header = GuppiRaw.Header()
header["CHAN_BW"]		= 0.5
header["OBSNCHAN"]	= n_ant*n_chan_perant
header["NANTS"]			= n_ant
header["NPOL"]			= n_pols
header["NBITS"]			= n_bits
header["BLOCSIZE"]	= (n_ant*n_chan_perant*n_time*n_pols*2*n_bits)/8
header["ANTNMS00"]	= "1cB,1eB"
header["ANTNMS01"]	= "1gB,1hB"
header["DEC_STR"]		= 21.9
header["RA_STR"]		= 19.9
header["DEC_OFF0"]	= 21.9
header["RA_OFF0"]		= 19.9
header["DIRECTIO"]	= 1
header["DUT1"]			= -0.055208199999999999
header["PIPERBLK"]	= n_time
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