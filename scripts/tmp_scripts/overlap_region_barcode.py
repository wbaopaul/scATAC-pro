import sys, os, re, math


## get # of overlap reads per barcode, given region files


input_sam = sys.argv[1]
region_file = sys.argv[2]
output_file = sys.argv[3]

fIN = open(input_sam, 'r')


line = fIN.readline()
overlaps = {}
nreads = {}
while (line):
	tmp = line.split('\t')[0:4]
	tmp0 = tmp[0]
	barcode = tmp0.split(':')[0]
	chr0 = tmp[2]
	start0 = int(tmp[3]) 
	end0 = start0 + 100
    
	if not barcode in nreads:
		nreads[barcode] = 0
	
	nreads[barcode] += 1

	if not barcode in overlaps:
        	overlaps[barcode] =0

	fIN_region = open(region_file, 'r')
	rline = fIN_region.readline()
	while (rline):
		rtmp = rline.split('\t')[0:3]
		rchr = rtmp[0]
		rstart = int(rtmp[1])
		rend = int(rtmp[2])
		if(chr0 == rchr):
			if(start0 <= rend & start0 >= rstart):
				
				overlaps[barcode] += 1
				break
			elif(end0 <= rend & end0 >= rstart):
                
				overlaps[barcode] += 1
				break
			
		
		rline = fIN_region.readline()
	
	line = fIN.readline()	


fIN.close()
fIN_region.close()


# start write
fOUT = open(output_file, 'w')

for barcode in nreads.keys():
	fOUT.write((barcode +"\t" + str(nreads[barcode]) +  "\t" + str(overlaps[barcode]) + "\n" ))


fOUT.close()



