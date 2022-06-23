#!usr/bin/env python3

import glob
from Bio.Seq import Seq 
from Bio.SeqUtils import GC 
from Bio.Seq import MutableSeq 
from Bio import SeqIO 
from Bio.SeqRecord import SeqRecord 
from Bio.SeqFeature import SeqFeature, FeatureLocation

#creating the .txt file output containing the number of pseudogenes at end of contigs and pseudogenes in the middle of contigs
pseudo_count = open('pseudogenes_counting.txt', 'w')
header = "Genome" + '\t' + "Pseudogenes at 3' and 5' end of contigs" + '\t' + "Selected pseudogenes" + '\n'
pseudo_count.writelines(header)

x = 0
for file in glob.glob('*.gbff'):


	#creating the .fasta file output without pseudogenes at 3' and 5' ends of contigs
	out = open(file.replace('.gbff','.fasta'), 'w')

	pseudo_end = 0
	pseudo_selected = 0

	print(x)
	x+=1

	locus_tag = []

	#searching pseudogenes over .gbff file
	for seq_record in SeqIO.parse(file, 'genbank'):

		for feature in seq_record.features:

			if feature.type == 'CDS':

				if 'pseudo' in feature.qualifiers:

					#pseudogenes at 3' and 5' ends of contigs
					if 'too short partial abutting assembly gap' in feature.qualifiers['note'][0]:
						pseudo_end+=1

					#pseudogenes selected
					else:
						pseudo_selected+=1

						locus_tag.append(feature.qualifiers['locus_tag'][0])


	counting = file + '\t' + str(pseudo_end) + '\t' + str(pseudo_selected) + '\n'
	pseudo_count.writelines(counting)



	#filtering cds_genomic.fna files

	fna = file.replace('_genomic.gbff','_cds_from_genomic.fna')

	for seq_rec in SeqIO.parse(fna, 'fasta'):

		flag = 0
		
		locus_tag_fasta = seq_rec.description.split('[locus_tag=')[1].split('] [')[0]

		for tag in locus_tag:

			if flag == 0:

				if tag == locus_tag_fasta:

					flag = 1

					SeqIO.write(seq_rec, out, 'fasta')

	out.close()


pseudo_count.close() 
			

			


					

	

	

	





