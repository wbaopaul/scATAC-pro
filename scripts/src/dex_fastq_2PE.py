# -*- coding: utf-8 -*-



## copy from snaptools with some modificatoin

import os.path
import sys 
import gzip 
import collections
import bz2



def file_type(filename):
    if filename.endswith('.gz'):
        return "gz"
    if filename.endswith('.bz2'):
        return "bz2"
    return "txt"

def dex_fastq(input_fastq1, input_fastq2,
              output_fastq1, output_fastq2,
              index_fastq
              ):
    
    """
    De-multiplex fastq files by adding barcode to the beginning of each read name.
    Required:
    --------
    input_fastq: 
        a fastq format file that contains the sequencing reads;
    output_fastq: 
        a fastq file contains output fastq file;
    index_fastq: 
        a list of fastq file
    """
    # check wheather fastq file exists
    if not os.path.exists(input_fastq1):
        print(('error: ' + input_fastq1 + ' does not exist!'));
        sys.exit(1);
        
    if not os.path.exists(input_fastq2):
        print(('error: ' + input_fastq2 + ' does not exist!'));
        sys.exit(1);

    if os.path.exists(output_fastq1):
        print(('error: ' + output_fastq1 + ' already exists, remove it first!'));
        sys.exit(1);
        
    if os.path.exists(output_fastq2):
        print(('error: ' + output_fastq2 + ' already exists, remove it first!'));
        sys.exit(1);
    
    if not os.path.exists(index_fastq):
        print(('error: ' + index_fastq + ' does not exist!'));
        sys.exit(1);
    
    if file_type(input_fastq1) == "gz":
        fr1 = gzip.open(input_fastq1, 'rb')
    elif file_type(input_fastq1) == "bz2":
        fr1 = bz2.BZ2File(input_fastq1, 'r')
    elif file_type(input_fastq1) == "txt":
        fr1 = open(input_fastq1, 'r')
        
    if file_type(input_fastq2) == "gz":
        fr2 = gzip.open(input_fastq2, 'rb')
    elif file_type(input_fastq2) == "bz2":
        fr2 = bz2.BZ2File(input_fastq2, 'r')
    elif file_type(input_fastq2) == "txt":
        fr2 = open(input_fastq2, 'r')
        
    if file_type(index_fastq) == "gz":
        fix = gzip.open(index_fastq, 'rb')
    elif file_type(index_fastq) == "bz2":
        fix = bz2.BZ2File(index_fastq, 'r')
    elif file_type(index_fastq) == "txt":
        fix = open(index_fastq, 'r')
    

    if output_fastq1.endswith("gz"):
        fout1 = gzip.open(output_fastq1, 'wb')
    elif output_fastq1.endswith("bz2"):
        fout1 = bz2.BZ2File(output_fastq1, 'wb')
    else:
        fout1 = open(output_fastq1, 'wb')            

        
    if output_fastq2.endswith("gz"):
        fout2 = gzip.open(output_fastq2, 'wb')
    elif output_fastq2.endswith("bz2"):
        fout2 = bz2.BZ2File(output_fastq2, 'wb')
    else:
        fout2 = open(output_fastq2, 'wb')  
        
    while True:
        cur_r1_name = fr1.readline().strip()[1:]
        cur_r2_name = fr2.readline().strip()[1:]
        
        if type(cur_r1_name) is bytes:
            cur_r1_name = cur_r1_name.decode();
        if cur_r1_name == "": break        
        cur_r1_read = fr1.readline().strip()
        cur_r1_plus = fr1.readline().strip()
        cur_r1_qual = fr1.readline().strip()
        if type(cur_r1_read) is bytes:
            cur_r1_read = cur_r1_read.decode();
        if type(cur_r1_qual) is bytes:
            cur_r1_qual = cur_r1_qual.decode();
        if type(cur_r1_plus) is bytes:
            cur_r1_plus = cur_r1_plus.decode();
            
        if type(cur_r2_name) is bytes:
            cur_r2_name = cur_r2_name.decode();
        if cur_r2_name == "": break        
        cur_r2_read = fr2.readline().strip()
        cur_r2_plus = fr2.readline().strip()
        cur_r2_qual = fr2.readline().strip()
        if type(cur_r2_read) is bytes:
            cur_r2_read = cur_r2_read.decode();
        if type(cur_r2_qual) is bytes:
            cur_r2_qual = cur_r2_qual.decode();
        if type(cur_r2_plus) is bytes:
            cur_r2_plus = cur_r2_plus.decode();
        
        cur_name = fix.readline().strip()[1:]
        cur_read = fix.readline().strip()
        cur_plus = fix.readline().strip()
        cur_qual = fix.readline().strip()
        if type(cur_name) is bytes:
            cur_name = cur_name.decode();
        if type(cur_read) is bytes:
            cur_read = cur_read.decode();
        if type(cur_plus) is bytes:
            cur_plus = cur_plus.decode();
        if type(cur_qual) is bytes:
            cur_qual = cur_qual.decode();
            
        cur_barcode = cur_read
            
        #if not (cur_name.split()[0] == cur_r1_name.split()[0]): sys.exit("read name does not match")        
        fout1.write(('@' + cur_barcode + ':' + cur_r1_name+"\n").encode('utf-8'))
        fout1.write((cur_r1_read+"\n").encode('utf-8'))
        fout1.write(("+\n").encode('utf-8'))
        fout1.write((cur_r1_qual+"\n").encode('utf-8')) 
        
        fout2.write(('@' + cur_barcode + ':' + cur_r2_name+"\n").encode('utf-8'))
        fout2.write((cur_r2_read+"\n").encode('utf-8'))
        fout2.write(("+\n").encode('utf-8'))
        fout2.write((cur_r2_qual+"\n").encode('utf-8'))    
        
    fout1.close()
    fout2.close()
    fr1.close()
    fr2.close()
    fix.close()

dex_fastq(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]);


