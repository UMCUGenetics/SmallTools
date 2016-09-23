#!/usr/bin/python
import vcf
import os
from optparse import OptionParser

parser = OptionParser()
parser.add_option("--vcf",	dest="vcf_file",     help="Path to VCF to convert",       default=False)
#parser.add_option("--conf",	dest="config_file",  help="Path to DataBase config file", default=False)
(options, args) = parser.parse_args()

if not os.path.exists(options.vcf_file):
	print "Invalid VCF file"
	exit(0)

TABLE_NAME = os.path.basename(options.vcf_file).replace("_design.vcf","")
VCF_READER = vcf.Reader(open(options.vcf_file, 'r'))

SNP_TEMPLATE = """%s tinyint(3) NOT NULL"""
SQL_TEMPLATE = """CREATE TABLE %s (
Sample varchar(100) NOT NULL,
%s,
PRIMARY KEY (Sample)
);
"""

SNPids = []
for vcf_record in VCF_READER:
	SNPids.append(vcf_record.ID)

print SQL_TEMPLATE%(TABLE_NAME,",\n".join([SNP_TEMPLATE%(x) for x in SNPids]))
