#!/usr/bin/python
import argparse


#HISEQ_HU01:54:C2FA5ACXX:8:2203:19371:16893	99	1	1591985	40	=	1592468	584	Library10
#HISEQ_HU01:70:H7J0DADXX:1:2213:9114:33460	99	1	1592007	40	=	1592377	471	Library1
#---------------------------------------------
#1	1592108	1592377	269	2	40	>Deletion_xxx_00000685<
#
#HISEQ_HU01:46:H0V93ADXX:2:1205:17357:94090	97	1	1591704	0	=	1657198	65574	Library7
#HISEQ_HU01:55:C2EMGACXX:2:2313:1834:73258	81	1	1657510	24	=	1591995	-65617	Library3
#---------------------------------------------
#1	1592096	1657198	65102	2	0	>Deletion_xxx_00000686<
#
#HISEQ_HU01:54:C2FA5ACXX:4:1109:20724:76865	81	1	1659302	0	=	1592899	-66504	Library8
#HISEQ_HU01:69:H7J06ADXX:2:1202:6069:88019	97	1	1593048	0	=	1659164	66188	Library6
#---------------------------------------------
#1	1593149	1659164	66015	2	0	>Deletion_xxx_00000687<

# ---------------------------------------------------------
  
def init():
  parser = argparse.ArgumentParser(prog='Parse_Delly_Calls.py', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('-i', help='Input calls file')
  parser.add_argument('-s', help='Sample definitions file')
  parser.add_argument('-o',  default='INFILE.calls.txt', help='Output tab delimited text file')
  parser.add_argument('-f',  default=False, type=bool, help='Flag to toggle filtering [True|False]')
  parser.add_argument('-qc', default=20, type=int, help='Minimum mapping quality [0-100]')
  parser.add_argument('-rp', default=10, type=int, help='Minimum number of supporting read pairs [1-x]')

  args = parser.parse_args()

  # check for input file
  if not args.i:
    parser.print_help()
    return(0)
    
  if not args.s:
    parser.print_help()
    return(0)

  if args.o == 'INFILE.calls.txt':
    args.o = args.i+'.calls.txt'
  
  print("INPUT:\t"+ args.i)
  print("OUTPUT:\t"+args.o)
  print("")
  
  
  return(args)

# ---------------------------------------------------------
  
def read_sample_names(samples_file):
  samples = {}

  ifile = open(samples_file, 'r')
  header = ifile.readline().strip().split('\t')

  for line in ifile:
    items = line.strip().split('\t')
    if len(items) <= 1:
      continue
    samples[items[0]] = items[1]
    
  ifile.close()
  
  return(samples)

# ---------------------------------------------------------

def make_header_string(samples):
 
  header = "Chromosome\tStart\tStop\tSize\tSupport\tQual\tID"
  
  sorter = sorted(samples)
  for sample_name in sorter:
    header = header + '\t' + sample_name
  
  return(header+'\n')
  
# ---------------------------------------------------------

def make_region_string(items, samples):
  #1	1593149	1659164	66015	2	0	>Deletion_xxx_00000687<  
 
  sorter = sorted(samples)
  line = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(*items)
  for sample_name in sorter:
    line = line + '\t' + str(samples[sample_name])
  
  return(line+'\n')



# ---------------------------------------------------------

def run(args):
  
  sample_ids = read_sample_names(args.s)
  
  
  samples = {}
  #reads = []
  
  for sample_id in sample_ids:
    samples[sample_ids[sample_id]] = 0
  
  
  ifile = open(args.i, 'r')
  ofile = open(args.o, 'w')
  
  ofile.write(make_header_string(samples))
  
  for line in ifile:
    items = line.strip().split('\t')
    
    if line.startswith("HISEQ"):
      #reads.append(items)
      samples[sample_ids[items[-1]]] += 1
      continue
      
    if line.startswith("----"):
      continue
    
    if len(items) <= 1:
      continue
  
    #summary line
    if args.f:
      #int(items[4]) >= args.rp and
      if int(items[5]) >= args.qc:	# Summary must match QC
	if max(samples.values()) >= args.qc:				# At least one sample should match RP
	  ofile.write(make_region_string(items, samples))
    else:
      ofile.write(make_region_string(items, samples))
    #reads=[]
    
    for sample_name in samples:
      samples[sample_name] = 0
  
  
  ofile.close()
  ifile.close()
      
# ---------------------------------------------------------


args = init()
if args is not 0:
  run(args)
