import random
import Levenshtein
import matplotlib.pyplot as plt

from Bio import pairwise2
from Bio.Cluster import kmedoids
from lib import parse_fastq, hamming
from numpy import zeros, logical_not, copy, unique
from itertools import islice


class HLAProcessor:
	"""
	HLAProcessor process HLA sequencing data for the IBCh lab. Data consist of
	two packs of reads: left reads and right reads. Left reads pack represent
	left parts of different HLA genes plus noise reads, and the right one 
	represent right parts plus noise.

	HLAProcessor is a command line tool for HLA genes assembling and genotyping.
	Genotyping is carried out with respect to the HLA genes database.

	All methods can be divided into two groups: one pack processign (clustering,
	distance diagrams, cluster consensus) and two pack processing (genes
	reconstruction, HLA genotyping)

	Futher you can find description of defferent HLAProcessor functions and
	abilities.
	"""

	def __init__(self, left_name, right_name, current='left'):
		"""
		HLAProcessor get names of left and right packs in constructor. It's
		possible to investigate content of only one pack (left or right)
		specified by current.
		"""

		self.pack_names = {'left': left_name, 'right': right_name}
		self.current = current if current in ['left', 'right'] else 'left'
		self.alignment_params = {
			'match': 1,
			'miss': -0.5,
			'gap_open': -2,
			'gap_extend': -0.5
		}

		# Pairwise dastance function: f(x, y)
		self.align = lambda x, y: pairwise2.align.globalms(x, y,
			self.alignment_params['match'],
			self.alignment_params['miss'],
			self.alignment_params['gap_open'],
			self.alignment_params['gap_extend'],
			one_alignment_only=True, score_only=True
		)

		self.distance = Levenshtein.distance

		# Length of substrings for computing Levenshtein distance
		self.seq_span = 150


	def __repr__(self):
		return 'Hi, I\'m HLAProcessor!'


	def set_current(self, current):
		self.current = current if current in ['left', 'right'] else 'left'


	def set_alignment_params(self, match, miss, gap_open, gap_extend):
		"""
		Set Affine Gap alignment parameters
		"""

		self.alignment_params = {
			'match': match,
			'miss': miss,
			'gap_open': gap_open,
			'gap_extend': gap_extend
		}


	def sampling(self, n, random_sampling=False):
		"""
		Sampling current pack n times. You can use random sampling algorithm by
		specifying random_sampling. Samples is a list of dictionaries:
		sample = {'name': '...', 'seq': '...', 'quality': [...]}.
		"""

		# It's posiible to use Biopython module SeqIO
		# SeqIO.parse('...', 'fastq')
		# It return list generator of SeqRecord objects

		if not random_sampling:
			self.samples = list( islice(
				parse_fastq(self.pack_names[self.current]), n) )
		else:
			# Non-optimal but fast in implementation
			full = list(parse_fastq(self.pack_names[self.current]))
			random.shuffle(full)

			self.samples = full[:n]


	def calc_distances(self):
		"""
		Compute left lower part of distance matrix
		"""

		n = len(self.samples)
		self.distances = []

		for i, s1 in enumerate(self.samples[1:]):
			for s2 in self.samples[:i]:
				self.distance.append(self.distance(s1, s2))


	def calc_distance_matrix(self):
		"""
		Compute distance matrix as numpy array
		"""

		n = len(self.samples)
		self.dist_matrix = zeros((n, n))

		for i, s1 in enumerate(self.samples):
			row = [
				self.distance(
					s1['seq'][:self.seq_span],
					s2['seq'][:self.seq_span]
				) for s2 in self.samples[(i + 1):]
			]

			self.dist_matrix[i, (i + 1):] = row
			self.dist_matrix[(i + 1):, i] = row


	def neighbours_hist(self, max_dist):
		"""
		Show histogram of neighbours number distribution. We suppose that two
		samples are neighbours if their Levenshtein distance less or equal than
		max_dist. This procedure alow us to exclude noise samples.
		"""

		selector = (self.dist_matrix <= max_dist)
		neighbours = selector.sum(axis=1)

		plt.hist(neighbours, bins=50)
		plt.show()


	def neighbours(self, max_dist):
		return (self.dist_matrix <= max_dist).sum(axis=1)


	def clean(self, max_dist, min_neighbours):
		selector = (self.dist_matrix <= max_dist)
		neighbours = selector.sum(axis=1)

		selector = neighbours < min_neighbours

		# Remove all selected samples
		samples = []
		for i, sample in enumerate(self.samples):
			if not selector[i]:
				samples.append(sample)

		self.samples = samples

		# I can't select from this! kmedoids returns segmentation fault!
		# TODO: try to understand
		self.calc_distance_matrix()

		print('%i samples was cleaned' % sum(selector))


	def clustering(self, nclusters=20, min_elements=None):
		self.ids, _, _ = kmedoids(
			self.dist_matrix,
			nclusters=nclusters,
			npass=20
		)


	# TODO: optimize nested constructions
	def grouping(self, min_elements=None):
		if min_elements is None:
			min_elements = 0.01 * len(self.samples)

		cluster_ids = unique(self.ids)
		self.clusters = {}
		self.cluster_ids = []

		for cluster_id in cluster_ids:
			if sum(self.ids == cluster_id) >= min_elements:
				
				# Make cluster
				self.cluster_ids.append(cluster_id)
				self.clusters[cluster_id] = []
				for i, sample in enumerate(self.samples):
					if self.ids[i] == cluster_id:
						self.clusters[cluster_id].append(sample['seq'])


	def make_consensus(self, cluster_id, n):
		profile = []

		s0 = self.clusters[cluster_id][0]
		length = len(s0)

		for s in self.clusters[cluster_id][1:n]:
			alignment = pairwise2.align.globalxx(
				s0, s, one_alignment_only = True)[0]

			profile.append(alignment[1])
			length = max(length, alignment[4])

		# Alignment of the first string
		profile.append(alignment[0])

		# Extending
		for i, sample in enumerate(profile):
			diff = length - len(s)
		if diff > 0:
			profile[i] = profile[i] + ('-' * diff)


		consensus = []
		symbols = ['A', 'C', 'G', 'T', '-']
		for pack in zip(*profile):
			freqs = {s: pack.count(s) for s in symbols}
			symbol = max(freqs, key=lambda s: freqs[s])

			if symbol == '-':
				symbol = ''
			
			consensus.append(symbol)

		return ''.join(consensus)








