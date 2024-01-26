from Bio import SeqIO

def readFastq(filename):
    """Parse Fastq file, tabulate list of sequences and associated BP qualities."""
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

def readGenome(filename):
    """Parse Fasta file."""
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

def max_orf(sequence, rf):
    """Given reading frame, calculate starting position of maximum length ORF in the sequence."""
    seq = sequence[rf-1:]
    len_max_orf = 0
    ind_max_orf = 0

    # Iterate over sequence
    for i in range(0, len(seq), 3):
        if seq[i:i + 3] == 'ATG': # check for 'start' codon
            # iterate over sub-sequence codons following 'start' codon. 
            for j in range(i + 3, len(seq) - 3, 3):
                if seq[j:j + 3] in ['TAA', 'TAG', 'TGA']: # stop iteration after finding 'stop' codon
                    len_orf = j+3-i # update current ORF length
                    # update max orf length and index
                    if len_orf > ind_max_orf:
                        len_max_orf = len_orf
                        ind_max_orf = i
                    break

    return ind_max_orf, len_max_orf

# identify all repeats of length repeat in sequence
def repeat_ns(sequence, repeat):
    """Find all repeated subsequences of length repeat in the given sequence."""
    rep_list = [] # stores repeats

    # iterate over sequence, stop 'repeat' bases before end
    for i in range(len(sequence) - repeat):
        rep_list.append(sequence[i:(i + repeat)])

    return rep_list

def fasta_processing(input_file, rf, n):
    """Wrapper function to perform analysis of input fasta file."""
    fasta_sequences = SeqIO.parse(open(input_file), 'fasta')

    n_records = 0
    len_records = {}
    orf_records = {}
    sub_sum = []

    # iterate over parse results
    for fasta in fasta_sequences:
        name, sequence, description = fasta.id, str(fasta.seq), str(fasta.description).split()
        desc = description[0]
        n_records += 1 # increment number of records
        seq_len = len(sequence) # note length of sequence
        len_records[desc] = seq_len 

        # calculate start position and the length of longest ORF
        orf_start, orf_len = max_orf(sequence, rf)
        orf_records[desc] = [orf_start, orf_len]

        # identify all repeats of length n
        substring_list = repeat_ns(sequence, n)
        for i in range(len(substring_list)):
            sub_sum.append(substring_list[i])

    sorted_record_len = {k: v for k, v in sorted(len_records.items(), key=lambda item: item[1], reverse=True)}
    orf_record_sorted = {k: v for k, v in sorted(orf_records.items(), key=lambda item: item[1][1], reverse=True)}

    return n_records, sorted_record_len, orf_record_sorted, sub_sum