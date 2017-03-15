import matplotlib.pyplot as plt

from HLAProcessor import *
from Bio.Cluster import kmedoids
from Bio.pairwise2 import align
from numpy import unique
from lib import parse_fasta

def init():
	p = HLAProcessor(
		'data/child10_S4_L001_R1_001.fastq',
		'data/child10_S4_L001_R2_001.fastq',
	)

	return p


def pipe1():
	p = HLAProcessor(
		'data/child10_S4_L001_R1_001.fastq',
		'data/child10_S4_L001_R2_001.fastq',
	)

	nsamples = 1000
	p.sampling(nsamples, random_sampling=True)
	p.seq_span = 80
	p.calc_distance_matrix()

	errors = []
	for nclusters in range(2, 10):
		_, error, nfound = kmedoids(p.dist_matrix, nclusters=nclusters, npass=10)
		print(nfound)

		errors.append(error)

	nclusters = list(range(2, 10))
	plt.plot(nclusters, errors)
	plt.show()


def pipe2():
	p = init()
	p.sampling(1000, random_sampling=True)
	p.calc_distance_matrix()
	p.seq_span = 80
	p.calc_distance_matrix()

	nclusters = 5
	clusterid, error, _ = kmedoids(p.dist_matrix, nclusters=nclusters, npass=20)


def pipe3():
	p = init()
	p.sampling(2000, random_sampling=True)
	p.calc_distance_matrix()
	# print(p.seq_span)
	p.clean(5, 50)

	#ids, _, _ = kmedoids(p.dist_matrix, nclusters=20, npass=20)
	#clusters = unique(ids)

	#for cluster in clusters:
	#	print('Cluster %i: %i' % (cluster, sum(ids == cluster)))
	p.clustering(nclusters=20)
	print(p.clusters)


def pipe4():
	consensus = \
		'CTGAGCTCCCCACTGGCTTTGGCTGGGGACACCCAACCACGTTTCCTGTGG \
		 CAGGGTAAGTATAAGTGTCATTTCTTCAACGGGACGGAGCGGGTGCAGTTC \
		 CTGGAAAGACTCTTCTATAACCAGGAGGAGTTCGTGCGCTTCGACAGCGAC \
		 GTGGGGGAGTACCGGGCGGTGACGGAGCTAGGGCGGCCTGTCGCCGAGTCC \
		 TGGAACAGCCAGAAGGACATCCTGGAGGACAGGCGGGGCCAGGTGGACACC \
		 GTGTGCAGACACAACTACGGGGTTGGTGAGAGCTTTACAGTGCAGC'

	with open('data/fullbase.fasta') as input_file:
		data = input_file.read()
		base = parse_fasta(data, destination = 'dict')

	best_score = 0
	for hla in base:



