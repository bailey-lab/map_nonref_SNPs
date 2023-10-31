mapped_psl=snakemake.input.alignment
query_seqs=snakemake.input.query_fasta
genome_seqs=snakemake.input.genome_fasta
bed_file=open(snakemake.output.bed_file, 'w')
best_alignments={}

bed_file.write('chrom\tpos\thap\tref\talt\n')
def read_fasta(fasta_file):
	seq, name_list, seq_list, seq_dict='', [], [], {}
	for line in open(fasta_file):
		line=line.strip()
		if '>' in line:
			name_list.append(line[1:])
			if len(seq)>0:
				seq_list.append(seq)
				seq=''
		else:
			seq=seq+line
	seq_list.append(seq)
#	for seq_number, name in enumerate(name_list):
#		seq_dict[name]=seq_list[seq_number]
	return [[name, seq_list[name_number]] for name_number, name in enumerate(name_list)]

def revcom(seq,nuc='DNA'):
	if nuc=='DNA':
		complement={'N':'N','n':'n','A':'T','C':'G','G':'C','T':'A','a':'t','t':'a','c':'g','g':'c', 'U':'A', 'u':'a', '-':'-'}
	else:
		complement={'N':'N','n':'n','A':'U','C':'G','G':'C','U':'A','a':'u','u':'a','c':'g','g':'c','-':'-'}
	return ''.join(reversed([complement[base] for base in seq]))


query_fasta=dict(read_fasta(query_seqs))
genome_fasta=dict(read_fasta(genome_seqs))
new_genome={}
for key in genome_fasta:
	new_key=key.split(' ')[0]
	new_genome[new_key]=genome_fasta[key]
genome_fasta=new_genome

print(genome_fasta.keys())

for line in open(mapped_psl):
	line=tuple(line.strip().split())
	if len(line)>5 and line[0].isdigit():
		score, qname=line[0], line[9]
		if qname not in best_alignments or int(score)>int(best_alignments[qname][0]):
			best_alignments[qname]=set([])
		if len(best_alignments[qname])==0 or int(score)==int(best_alignments[qname][0]):
			best_alignments[qname].add(line)
#			print('best is', best_alignments)
for query in best_alignments:
#	print('\n\nquery is', query)
	for line in best_alignments[query]:
		score, strand, q_name, q_size, q_start, q_end, t_name, t_start, t_end, block_sizes, q_starts, t_starts, q_seqs, t_seqs=line[0], line[8], line[9], line[10], line[11], line[12], line[13], line[15], line[16], line[18], line[19], line[20], line[21], line[22]
		print('line is', line)
		block_sizes, q_starts, t_starts, q_seqs, t_seqs=block_sizes.split(','), q_starts.split(','), t_starts.split(','), q_seqs.split(','), t_seqs.split(',')
		q_block_end, t_block_end=int(q_starts[0]), int(t_starts[0])
		for block_number, block in enumerate(t_starts):
			print(block)
			if len(block)<1:
				break
			t_start=int(t_starts[block_number])
			q_start=int(q_starts[block_number])
			#this code doesn't yet handle reverse complemented indels
			if q_block_end!=int(q_starts[block_number]):
				print('entering q_block')
				if strand=='-':
					seq=revcom(query_fasta[q_name])
				else:
					seq=query_fasta[q_name]
				print(q_block_end, q_starts[block_number])
				bed_file.write(f'{t_name}\t{position}\t{q_name}\tMISSING\t{seq[q_block_end:int(q_starts[block_number])]}\n')
			if t_block_end!=int(t_starts[block_number]):
				print('entering t_block')
				print(t_block_end, t_starts[block_number])
				bed_file.write(f'{t_name}\t{t_block_end}\t{q_name}\t{genome_fasta[t_name][t_block_end:int(t_starts[block_number])]}\tMISSING\n')
			for position_number, q_letter in enumerate(q_seqs[block_number]):
				t_letter=t_seqs[block_number][position_number]
				position=t_start+position_number
				if q_letter!=t_letter:
					bed_file.write(f'{t_name}\t{position}\t{q_name}\t{t_letter.upper()}\t{q_letter.upper()}\n')
			t_block_end=t_start+position_number+1
			q_block_end=q_start+position_number+1
