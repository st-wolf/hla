def grouper(iterable, n):
	args = [iter(iterable)] * n
	return zip(*args)


def hamming(xs, ys):
	count = 0
	for x, y in zip(xs, ys):
		count += (x != y)

	return count


def parse_fastq(file_name):
	"""
	Take fastq file name and return parsed file generator.
	"""

	with open(file_name) as input_file:
		for name, seq, _, quality_string in grouper(input_file, 4):

			# Drop last symbol ('\n') from every string
			quality = [ord(x) - ord('!') for x in quality_string[:-1]]

			yield {
				'name': name[:-1],
				'seq': seq[:-1],
				'quality': quality
			}


def parse_fasta(data, destination = 'dict'):
    def parse_to_dict(data):
        fasta = {}

        for item in data.split('>')[1::]:
            item = item.split('\n')
            fasta[item[0]] = ''.join(item[1::])

        return fasta

    def parse_to_list(data):
        fasta = []

        for item in data.split('>')[1::]:
            fasta.append(''.join(item.split('\n')[1::]))

        return fasta

    if destination == 'dict':
        return parse_to_dict(data)
    elif destination == 'list':
        return parse_to_list(data)