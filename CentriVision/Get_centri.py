from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import CentriVision.bez as bez

class Get_centri():
	def __init__(self, options):
		# self.a = ""
		for k, v in options:
			setattr(self, str(k), v)
			print(k, ' = ', v)

	def run(self):
		zsl1 = {}
		for line in open(self.gff_file,'r'):
			if line[0] == '#':
				continue
			lt = line.strip('\n').split()
			chro = lt[0]
			if chro == 'Chr_ID':
				continue
			if chro not in zsl1.keys():
				zsl1[chro] = [int(lt[1]),int(lt[2]),abs(int(lt[1])-int(lt[2]))]
			else:
				zsl1[chro].append(int(lt[1]))
				zsl1[chro].append(int(lt[2]))
				zsl1[chro].append(abs(int(lt[1])-int(lt[2])))
		out = open(self.out_fasta,'w')
		for seq_record in SeqIO.parse(self.genome_file, "fasta"):# biopython循环读取
			if seq_record.id not in zsl1.keys():
				continue
			for i in range(int(len(zsl1[seq_record.id])/3)):
				name = seq_record.id+'_'+str(i+1)
				seq0 = str(seq_record.seq)[zsl1[seq_record.id][i*3]-1:zsl1[seq_record.id][i*3+1]]
				out.write('>'+name+'\n'+seq0+'\n')
		out.close()