#!/usr/local/bin/python

##################################################################
# A table will be made that contains all SVs from each VCF-file, with the SVs sorted on type and length.
# Each vcf file is put in a different column of the table.
# In the command line a path should be provided of the vcf files to be processed. 
#
# Author: Floor Dussel
##################################################################

import os
import re
import numpy as np
import pandas as pd
import glob
import vcf
import argparse

def categorize(path):
	list_vcf_names = []
	df_all = pd.DataFrame()

	def extractInfoINVDELTRADUP (record):
		chrom = record.CHROM
		pos = record.POS
		end = record.INFO["END"]
		length_calc = end - pos
		length_abs = abs(length_calc)
		length = length_abs + 1
		output = [chrom, pos, end, length]
		return output;

	def extractInfoINS (record):
		chrom = record.CHROM
		pos = record.POS
		end = record.INFO["END"]
		length = 0
		if "INSLEN" in record.INFO:
			length_get = record.INFO["INSLEN"]
			length = abs(length_get)
			#print 'INSLENGTH used as length of insertion'
		elif "SVLEN" in record.INFO:
			length = record.INFO["SVLEN"]
			if length != int:
				length = length[0]
			if length < 0:
				length = length * -1
			#print 'SVLEN used as length of insertion'
		#else:
			#print 'No insertion length found. Note: insertion is not taken into account in the table'

		output = [chrom, pos, end, length]
		return output;

	#BNDs are said to have length = 0 and END is said to be the same as POS
	def extractInfoBND (record):
		chrom = record.CHROM
		pos = record.POS
		end = pos
		length = 0
		output = [chrom, pos, end, length]
		return output;

	def checkLengthSVs(record):
		if record.INFO["SVTYPE"] == "DEL":
			del_output = extractInfoINVDELTRADUP(record)
			length = del_output[3]

			if (length > 1000 and length < 10001): #1000-10001bp, so length of 1-10kb
				del_outputlist10.append(del_output)
			elif (length > 10000 and length <100001): #10-100kbp
				del_outputlist100.append(del_output)
			elif (length > 100000 and length <1000001): #100kb-1mbp
				del_outputlist1000.append(del_output)
			elif (length > 1000000 and length <10000001): #1mbp - 10 mbp
				del_outputlist10000.append(del_output)
			elif (length > 10000000): # >10Mb
				del_outputlist100000.append(del_output)

		if record.INFO["SVTYPE"] == "INS":
			ins_output = extractInfoINS(record)
			length = ins_output[3]

			if (length > 1000 and length < 10001): #1000-10001bp, so length of 1-10kb
				ins_outputlist10.append(ins_output)
			elif (length > 10000 and length <100001): #10-100kbp
				ins_outputlist100.append(ins_output)
			elif (length > 100000 and length <1000001): #100kb-1mbp
				ins_outputlist1000.append(ins_output)
			elif (length > 1000000 and length <10000001): #1mbp - 10 mbp
				ins_outputlist10000.append(ins_output)
			elif (length > 10000000): # >10Mb
				ins_outputlist100000.append(ins_output)

		if record.INFO["SVTYPE"] == "INV":
			inv_output = extractInfoINVDELTRADUP(record)
			length = inv_output[3]
			if (length > 1000 and length < 10001): #1000-10001bp, so length of 1-10kb
				inv_outputlist10.append(inv_output)
			elif (length > 10000 and length <100001): #10-100kbp
				inv_outputlist100.append(inv_output)
			elif (length > 100000 and length <1000001): #100kb-1mbp
				inv_outputlist1000.append(inv_output)
			elif (length > 1000000 and length <10000001): #1mbp - 10 mbp
				inv_outputlist10000.append(inv_output)
			elif (length > 10000000): # >10Mb
				inv_outputlist100000.append(inv_output)

		if  record.INFO["SVTYPE"] == "TRA":
			tra_output = extractInfoINVDELTRADUP(record)
			tra_bnd_outputlist.append(tra_output)

		if record.INFO["SVTYPE"] == "DUP":
			dup_output = extractInfoINVDELTRADUP(record)
			length = dup_output[3]

			if (length > 1000 and length < 10001): #1000-10001bp, so length of 1-10kb
				dup_outputlist10.append(dup_output)
			elif (length > 10000 and length <100001): #10-100kbp
				dup_outputlist100.append(dup_output)
			elif (length > 100000 and length <1000001): #100kb-1mbp
				dup_outputlist1000.append(dup_output)
			elif (length > 1000000 and length <10000001): #1mbp - 10 mbp
				dup_outputlist10000.append(dup_output)
			elif (length > 10000000): # >10Mb
				dup_outputlist100000.append(dup_output)

		if record.INFO["SVTYPE"] == "BND":
			bnd_output = extractInfoBND(record)
			tra_bnd_outputlist.append(bnd_output)

	path_vcf = path
	#path = "/home/cog/fdussel/Documents/*.vcf"
	for vcf_filename in glob.glob(path_vcf):
		#vcf= open(vcf_filename)
		#print vcf_filename
		vcf_file = None
		df_vcf = pd.DataFrame()
		vcf_file_handle = open(vcf_filename, 'r')

		if os.path.getsize(vcf_filename) > 0:
			vcf_file = vcf.Reader(vcf_file_handle)
		else:
			continue

		del_outputlist10 = []
		del_outputlist100 = []
		del_outputlist1000 = []
		del_outputlist10000 = []
		del_outputlist100000 = []
		ins_outputlist10 = []
		ins_outputlist100 = []
		ins_outputlist1000 = []
		ins_outputlist10000 = []
		ins_outputlist100000 = []
		inv_outputlist10 = []
		inv_outputlist100 = []
		inv_outputlist1000 = []
		inv_outputlist10000 = []
		inv_outputlist100000 = []
		tra_bnd_outputlist = []
		dup_outputlist10 = []
		dup_outputlist100 = []
		dup_outputlist1000 = []
		dup_outputlist10000 = []
		dup_outputlist100000 = []

		all_names = vcf_filename.rsplit('/', 1)
		name = all_names[1]
		list_vcf_names.append(name)

		for record in vcf_file:
			checkLengthSVs(record)


			dict_numberSV = {"del_1_10kb" : len(del_outputlist10), "del_10_100kb": len(del_outputlist100), "del_100kb_1mb": len(del_outputlist1000),
						"del_1mb_10mb" : len(del_outputlist10000), "del_>10mb": len(del_outputlist100000),
						"ins_1_10kb" : 	len(ins_outputlist10), "ins_10_100kb" : len(ins_outputlist100), "ins_100kb_1mb": len(ins_outputlist1000),
						"ins_1mb_10mb" : len(ins_outputlist10000), "ins_>10mb" : len(ins_outputlist100000),
			   			"inv_1_10kb" : len(inv_outputlist10), "inv_10_100kb": len(inv_outputlist100), "inv_100kb_1mb" : len(inv_outputlist1000),
						"inv_1mb_10mb" : len(inv_outputlist10000), "inv_>10mb" : len(inv_outputlist100000),
						"tra/bnd ": len(tra_bnd_outputlist),
						"dup_1_10kb" : len(dup_outputlist10), "dup_10_100kb" : len(dup_outputlist100), "dup_100kb_1mb" : len(dup_outputlist1000),
						"dup_1mb_10mb" : len(dup_outputlist10000), "dup_>10mb" : len(dup_outputlist100000)}

			df_vcf = pd.DataFrame({name: dict_numberSV}, index = dict_numberSV.keys())

		vcf_file_handle.close()

		# Check whether data frame is empty
		if not df_vcf.empty:
			df_all = pd.concat([df_all, df_vcf], axis=1)
			df_all = df_all.sort_index()

	#print list_vcf_names

	print df_all

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'Categorizing structural variants on type and length')
	required_named = parser.add_argument_group('Required arguments')
	required_named.add_argument('-p', '--path', help = "input a path for the vcf-file(s) to be processed", required=True)

	args = parser.parse_args()
	categorize(args.path)
