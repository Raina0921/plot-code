#!/usr/bin/python2
# Sequence Translator (Foward Frames Only)
# James Wright 2019

#used to get cmd line arguments 
import argparse
import re


#Read command line arguments and create help documentation using argparse
parser = argparse.ArgumentParser(
	description='''Read in na FASTA file and translate to protein. james.wright@icr.ac.uk 2019''')
parser.add_argument('fasta', metavar='*.fasta|*.fa', help='FASTA file of proteins sequences to be translated.')
parser.add_argument('--output_fasta', '-o', dest='dout', default='TranslatedProteins.fa', help='Set file to write proteins. Default=TranslatedProteins.fa')
parser.add_argument('--3frame', '-f', dest='frames', default=False, action='store_true', help='Translated in 3 frames. Default=false')
parser.add_argument('--no_isobaric', '-i', dest='iso', default=False, action='store_true', help='Do not make proteins isobaric. Default=false')
parser.add_argument('--longest_orf', '-l', dest='longest', default=False, action='store_true', help='Only store the longest ORF for each Frame. Default=false')
parser.add_argument('--min_orf_size', '-m', type=int, dest='osize', default=10, help='Set minimum amino acid length of ORFs saved. Default=10')
#parser.add_argument('--id_filter', '-x',  dest='flist', default='', help='A file containing a list of IDs to include in the results. Default=NA (include all ids)')

args = parser.parse_args()


#Basic genetic code to use for translation
gcode = {
      'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
      'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
      'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
      'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
      'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
      'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
      'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
      'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
      'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
      'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
      'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
      'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
      'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
      'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
      'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
      'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}


#Function to translate na into aa sequence
def translate ( se ):
      translate = ''.join([gcode.get(se[3*i:3*i+3],'X') for i in range(len(se)//3)])
      return translate


#Function to process each FASTA sequence by selected frame(s) and orf filtering parameters
def process_sequence ( seq, id, nframes, iso, osize, longest, oFa  ):
	olen = 0	#longest ORF length
	lf = 99		#longest ORF frame
	lc = -1		#longest ORF id
	lorf =''	#longest ORF sequence
	lspos = -1 	#longest ORF transcript start position
	lcnt = 1;

	#split id on any white space so acession can be modified
	idsp = id.split(None, 1)

	grid = idsp[0].replace('>','')
	trid = re.sub('_iso\d+', '', grid)

	gid = grid
	if len(idsp) > 1: gid = idsp[1]


	#Loop each forawrd frame (1 or all 3)
	for frame_start in range(nframes):

		#Translate sequence using function
		tseq = translate(seq[frame_start:])

		current_start_pos = frame_start

		#make sequence isobaric (check args for switch off)
		if iso == False:
			tseq = tseq.replace('I', 'L')

		#split sequence on stop codons "_" and loop through each one
		for ocnt, orf in enumerate( tseq.split("_") ):

			#print ("orf::" + orf + "::" + str(len(orf)) + "::" + str(osize))

			#Does the orf meet minimum length requirement from args
			if len(orf) >= osize:

				#if option to only retain longest ORF is set check current orf against longest so far and update stored information
				if longest == True:

					if (len(orf) > olen):
						olen = len(orf)
						lorf = orf
						lf = frame_start
						lc = ocnt
						lspos = current_start_pos
						
				#>alt_3prime.4188_iso1|1 frame=0 offset=0 orf=0 gene=ENSG00000188976.9 geneID=alt_3prime.4188 transcriptID=alt_3prime.4188_iso1
				#>alt_3prime.4188_iso1_Frame0_spos0 orf=0 gene=ENSG00000188976.9

				#Otherwise write each translated ORF into new FASTA
				else:
					oFa.write(idsp[0] + '|' + str(lcnt) + ' frame=' + str(frame_start) + ' offset=' + str(current_start_pos) + ' orf=' + str(ocnt) + ' gene=' + gid  + ' geneID=' + grid + ' transcriptID=' + trid + '\n')
					oFa.write(orf + '\n')
					lcnt+=1

			current_start_pos = current_start_pos + ( len(orf) * 3 ) + 3

	#If selected write longets ORF to FASTA after looping all ORFs in all Frames
	if longest == True:
		#Does the longest ORF meet minimum length requirements
		if len(lorf) >= osize:
			oFa.write(idsp[0] + '|' + str(lcnt) + ' frame=' + str(lf) + ' offset=' + str(lspos) + ' orf=' + str(lc) + ' gene=' + gid  + ' geneID=' + grid + ' transcriptID=' + trid + '\n')
			oFa.write(lorf + '\n')

	return

#################################################

#empty protein sequence, protein id
sequ = ''	
proID = ''

#Set number of frames to be translated
nf = 1
if args.frames == True:
	nf = 3

#print (parser.parse_args())	

#Open FASTA file for writing translated sequences
outFasta = open(args.dout, 'w')

#Open FASTA file using first cmd line argument
fasta = open(args.fasta, 'r')

#loop each line in the file
for line in fasta:
	#if this line starts with ">" then process sequence if not empty
	if line[0] == '>':
		if sequ != '':

			#Process and translate current sequence
			process_sequence(sequ.upper(), proID, nf, args.iso, args.osize, args.longest, outFasta)

		#Extract next sequence identifier (FASTA Header)
		proID = line.rstrip()
		#New sequence
		sequ = ''

	#if not accession line then append aa sequence (with no newline or white space) to seq string
	else: 
		sequ+=line.rstrip()
		
#Close files
fasta.close()

#At end of FASTA file process final sequence
if sequ != '':
	process_sequence(sequ.upper(), proID, nf, args.iso, args.osize, args.longest, outFasta)


outFasta.close()
