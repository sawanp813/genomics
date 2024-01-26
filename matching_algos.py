import bisect

class MatchIndex(object):
    """Holds a substring index for a text T."""

    def __init__(self, t, k):
        """Create index from all substrings of t of length k."""
        self.k = k  # k-mer length (k)
        self.index = []
        for i in range(len(t) - k + 1):  # for each k-mer
            self.index.append((t[i:i+k], i))  # add (k-mer, offset) pair
        self.index.sort()  # sort index to facilitate binary search

    def query(self, p):
        """Return index hits for first k-mer of p."""
        kmer = p[:self.k]  # query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1))  # binary search
        hits = []

        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits
    
class SubseqIndex(object):
    """ Holds a subsequence index for a text T """
    
    def __init__(self, t, k, ival):
        """ Create index from all subsequences consisting of k characters
            spaced ival positions apart.  E.g., SubseqIndex("ATAT", 2, 2)
            extracts ("AA", 0) and ("TT", 1). """
        self.k = k  # num characters per subsequence extracted
        self.ival = ival  # space between them; 1=adjacent, 2=every other, etc
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):  # for each subseq
            self.index.append((t[i:i+self.span:ival], i))  # add (subseq, offset)
        self.index.sort()  # alphabetize by subseq
    
    def query(self, p):
        """ Return index hits for first subseq of p """
        subseq = p[:self.span:self.ival]  # query with first subseq
        i = bisect.bisect_left(self.index, (subseq, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits
    
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

def index_approximate_match(p, t, t_index, n):
    """Find approximate substring matches (length n) using an index on the context string."""
    segment_length = int(len(p) / (n + 1))
    hits = []
    all_match = set() # avoid repeated adds

    for cod in range(n + 1):
        beg = cod * segment_length
        end = (cod + 1) * segment_length
        segment = p[beg:end]
        
        i = bisect.bisect_left(t_index.index, (segment, -1))
        
        while i < len(t_index.index):
            if segment != t_index.index[i][0]:
                break
            else:
                t_offset = t_index.index[i][1]
                hits.append(t_offset)                    
            i += 1            
        
        for h in hits:
            if beg > h or h-beg + len(p) > len(t):
                continue
                
            errs = 0
            for x in range(0, beg):
                if p[x] != t[h-beg+x]:
                    errs += 1
                if errs > n:
                    break
                    
            for y in range(end, len(p)):
                if p[y] != t[h-beg+y]:
                    errs += 1
                if errs > n:
                    break
                    
            if errs <= n:
                all_match.add(h-beg)
                                
    return hits, list(all_match)

def index_approximate_subseq_match(p, t, t_subseq_index, n):
    """Compute approximate matches between query string and context string."""
    hits = []    
    for position in range(n):
        subseq = p[position:][:22:3]
        
        i = bisect.bisect_left(t_subseq_index.index, (subseq, -1))  # binary search

        while i < len(t_subseq_index.index):  # collect matching index entries
            if t_subseq_index.index[i][0] != subseq:
                break
            hits.append(t_subseq_index.index[i][1])
            i += 1
            
    return hits

def editDistance(x, y):
    """Compute minimum edit distance between two strings via DP."""
    # Create distance matrix
    D = []
    for i in range(len(x)  +1):
        D.append([0] * (len(y) + 1 ))

    # Initialize first row and column of matrix
    for i in range(len(x) + 1):
        D[i][0] = i
    for i in range(len(y) + 1):
        D[0][i] = i

    # Fill in the rest of the matrix
    for i in range(1, len(x) + 1):
        for j in range(1, len(y) + 1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if x[i-1] == y[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)

    return D[-1][-1]

# genome file
# f = 'ERR037900_1.first1000.fastq'
f = 'chr1.GRCh38.excerpt.fasta'
genome = readGenome(f) # just take the first read

# query string
p = 'GGCGCGGTGGCTCACGCCTGTAAT'
print(f"The query string is: {p}.")

# find matches to genome by splitting query string into n parts
n = 2
gen_ind = MatchIndex(genome, 8)
hits, matches = index_approximate_match(p, genome, gen_ind, n)
print("Indices of approximate matches: ", hits)

# Crete index of all three-mers contained in genome string
t_subseq_index = SubseqIndex(genome, 8, 3)
hits = index_approximate_subseq_match(p, genome, t_subseq_index, 2)
print("Indices of approximate matches for 3-mers: ", hits)

d = editDistance(p, genome)
print("Minimum edit distance between read and genome: ", d)
