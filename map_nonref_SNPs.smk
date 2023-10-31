'''
maps haplotypes onto a genome, finds best mapping location for every haplotype,
and compares each position of the haplotype against the genome to see if it
differs from the reference.
'''

fasta_file='/nfs/jbailey5/baileyweb/asimkin/other_people/jmsadler/sample_mixup/summary.fa'
genome='/nfs/jbailey5/baileyweb/bailey_share/UCSC_genome_browsers/private_UCSC_genome_browser/Pf_3D7/Pf_3D7_genome/PlasmoDB-62_Pfalciparum3D7_Genome.fasta'

rule all:
	input:
		mapped_out='nonref_SNPs.bed'

rule map_haps:
	input:
		genome=genome,
		fasta_file=fasta_file
	output:
		alignment='mapped.pslx'
	shell:
		'blat {input.genome} {input.fasta_file} -out=pslx {output.mapped_out}'

rule parse_alignment:
	input:
		alignment='mapped.pslx',
		query_fasta=fasta_file,
		genome_fasta=genome
	output:
		bed_file='nonref_SNPs.bed'
	script:
		'scripts/parse_alignment.py'
