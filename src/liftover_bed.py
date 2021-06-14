#!/usr/bin/env python

'''
Title: liftover_bed.py
Date: 20200918
Author: Adam Nunn
Description:
	Takes a position-sorted BED file with unspecific reference coordinates and derives the
	true coordinates from a corresponding, read name-sorted BAM file where reference sequences
	have been mapped against the appropriate genome sequence.
	...
	### BED
	1 123 124
	1 147 148
	3  50  51
	...
	### BAM
	1  99 Chromosome1 1245 ... etc. ATCHACHAHCATG *
	1 147 Chromosome1 1459 ... etc. ATCYHAYCHTACG *

List of functions:
	main()
	check_mate()

Procedure:
	1. Stage initial BAM file for reading
	2. Iterate through each pos in BED file, check for corresponding scaffold in BAM 
	3. On match, find proper paired mate (next line or so)
	4. Get aligned pairs (matches only) and infer read/mate length
	5. Get coordinates from aligned pairs, make adjustments (+/- 1) for R1/R2 and for/rev
	6. Cross-reference pos with mapped coordinates, if present print new coordinates

Usage:
	./liftover_bed.py [-h, --help] [bam] [bed]
eg. ./liftover_bed.py in.bam in.bed
'''

###################
## INIT ENVIRONMENT

import argparse
import pysam

# main process function
def main(BAM,BED):
	new_mate="False"
	with pysam.AlignmentFile(BAM) as bam, open(BED, "r") as bed:

		comp = str.maketrans("ACTG","TGAC")

		# setup vars for initial bam alignments
		try: read = next(bam)
		except StopIteration: 
			print("No alignments detected in bam file")
			raise SystemExit(0)

		rnam = int(read.query_name)

		# iterate through VCF file
		for line in bed:	

			line = line.rstrip()
			line = line.split("\t")
			contig = int(line[0])
			pos = int(line[2])

			# if BAM is behind VCF then iterate 
			while rnam < contig:
				read = next(bam)
				rnam = int(read.query_name)
				new_mate = True
			
			# if BAM is ahead of VCF then current SNP has no equivalent alignment (skip it)
			if (rnam > contig) or read.is_unmapped \
				or read.is_secondary \
					or read.is_qcfail \
						or read.is_duplicate \
							or read.is_supplementary \
								or not read.is_proper_pair:
				new_mate = True
				continue 
			
			# get read mate
			if new_mate:
				mate = check_mate(bam)
				if rnam != int(mate.query_name): continue
				else: new_mate = False

			# get bp alignment (matches only)
			read_pairs = read.get_aligned_pairs(matches_only = True)
			read_pos = [i[0] for i in read_pairs]
			read_ref = [i[1] for i in read_pairs]
			#print(read_pairs)
			#print(read_pos)
			#print(read_ref)

			mate_pairs = mate.get_aligned_pairs(matches_only = True)
			mate_pos = [i[0] for i in mate_pairs]
			mate_ref = [i[1] for i in mate_pairs]
			#print(mate_pairs)
			#print(mate_pos)
			#print(mate_ref)

			# infer read lengths
			read_length = read.infer_read_length()
			mate_length = mate.infer_read_length()
			#print("{}\t{}".format("read length",read_length))
			#print("{}\t{}".format("mate length",mate_length))

			#print("{}\t{}\t{}\t{}\t{}\t{}".format(read.query_name,read.reference_name,read.is_read1,read.is_reverse,SNP.contig,SNP.pos))
			#print("{}\t{}\t{}\t{}\t{}\t{}".format(mate.query_name,mate.reference_name,mate.is_read1,mate.is_reverse,SNP.contig,SNP.pos))

			# re-invent the wheel (read)
			if read.is_read2 and read.is_reverse:
				STRAND = "+"
				read_pos = [i + mate_length for i in read_pos]

			elif read.is_read1 and read.is_reverse:
				STRAND = "-"
				true_pos = [i+1 for i in range(read_length)][::-1]
				read_pos = [true_pos[p] + 1 for p in read_pos]

			elif read.is_read2:
				STRAND = "-"
				true_pos = [i+1 for i in range(read_length)][::-1]
				read_pos = [true_pos[p] + mate_length + 1 for p in read_pos]

			else: STRAND = "+"
			
			# re-invent the wheel (mate)
			if mate.is_read2 and mate.is_reverse:
				mate_pos = [i + read_length for i in mate_pos]

			elif mate.is_read1 and mate.is_reverse:
				true_pos = [i+1 for i in range(mate_length)][::-1]
				mate_pos = [true_pos[p] + 1 for p in mate_pos]

			elif mate.is_read2:
				true_pos = [i+1 for i in range(mate_length)][::-1]
				mate_pos = [true_pos[p] + read_length + 1 for p in mate_pos]

			# combine pairs index for read and mate
			spos = read_pos + mate_pos
			rpos = read_ref + mate_ref

			# reverse complement ref and alt for denovo clusters aligning to crick
			#if STRAND == "-":
				#SNP.ref = SNP.ref.translate(comp)
				#SNP.alts = (x.translate(comp) for x in SNP.alts)

			# cross-reference SNP position, print results
			if pos in read_pos: print('{}\t{}\t{}\tread\t{}\t{}'.format(read.reference_name,read_ref[read_pos.index(pos)]-1,read_ref[read_pos.index(pos)],STRAND,"\t".join(line)))
			elif pos in mate_pos: print('{}\t{}\t{}\tmate\t{}\t{}'.format(mate.reference_name,mate_ref[mate_pos.index(pos)]-1,mate_ref[mate_pos.index(pos)],STRAND,"\t".join(line)))


# recursive mate check function
def check_mate(bam):

	mate = next(bam)
	if mate.is_unmapped \
		or mate.is_secondary \
			or mate.is_qcfail \
				or mate.is_duplicate \
					or mate.is_supplementary \
						or not mate.is_proper_pair: mate = check_mate(bam) 
	
	return mate


# define argparse
usage = ''' Takes a position-sorted BED file with unspecific reference coordinates and derives the
	true coordinates from a corresponding, read name-sorted BAM file where reference sequences
	have been mapped against the appropriate genome sequence. '''

parser = argparse.ArgumentParser(description=usage)

parser.add_argument('bam', metavar='<BAM>', help='[REQUIRED] Path to input BAM file')
parser.add_argument('bed', metavar='<BED>', help='[REQUIRED] Path to input BED file')

args = parser.parse_args()
parser.parse_args()

# run main()
if __name__ == '__main__': main(args.bam,args.bed)
