# Implements Boyer-Moore algorithm (with weak GSR)
from gen_file_parser import readFastq, readFasta

def comp_z(seq):
    """Preprocess string to retrive Z array.
    Inputs:
        seq: string to preprocess
    Outputs:
        Z: Array, Z[i] is the length of the substring beginning at i which is also a prefix
        of seq. 
    """
    seq_len = len(seq)
    assert seq_len > 1, "Input a string with a length greater than 1."
    z = [seq_len] + [0] * (seq_len-1)

    # Initial comparison of s[1:] with prefix
    for i in range(1, seq_len):
        if seq[i] == seq[i-1]:
            z[1] += 1 # z[1] stores length of match of substrings of seq beginning at 0 and 1
        else:
            break
    # indices for z-box
    r = 0
    l = 0
    # NOTE: [l, r] is interval with max r s.t. [l, r] is a prefix  substring
    # if there is a matching base substring
    if z[1] > 0:
        r = z[1]
        l = 1

    # rest of z-array
    for j in range(2, seq_len):
        assert z[j] == 0
        # if outside 'z-box', do explicit check
        if j > r:
            for i in range(j, seq_len):
                if seq[i] == seq[i-j]:
                    z[j] += 1
                else:
                    break
            r = j + z[j] - 1
            l = j
        else: # if within z-box
            # Calculate length of beta
            nbeta = r - j + 1
            zkp = z[j - l] # length of matching substring starting at 0 and 1
            if nbeta > zkp:
                # Beginning of box
                z[j] = zkp
            else:
                # Compare characters past r
                matches = 0

                for i in range(r+1, seq_len):
                    if seq[i] == seq[i - j]:
                        matches += 1
                    else:
                        break
                l = j
                r = r + matches
                z[j] = r - j + 1
    return z


def comp_n(s):
    """Compile the N array from the Z array."""
    return comp_z(s[::-1])[::-1]


def big_l_prime_array(p, n):
    """Compile L' array using p and N array.
    L'[i] = largest index j less than n such that N[j] = |P[i:]|."""
    ip = [0] * len(p)
    for j in range(len(p)-1):
        i = len(p) - n[j]
        if i < len(p):
            ip[i] = j + 1
    return ip


def big_l_array(p, lp):
    """Compile L array using p and L' array.
    L[i] = largest index j less than n such that N[j] >= |P[i:]|."""
    l = [0] * len(p)
    l[1] = lp[1]
    for i in range(2, len(p)):
        l[i] = max(l[i-1], lp[i])
    return l


def small_l_prime_array(n):
    """Compile lp' array using N array."""
    small_lp = [0] * len(n)

    for i in range(len(n)):
        if n[i] == i+1:  # prefix matching a suffix
            small_lp[len(n)-i-1] = i+1

    for i in range(len(n)-2, -1, -1):  # "smear" them out to the left
        if small_lp[i] == 0:
            small_lp[i] = small_lp[i+1]

    return small_lp


def good_suffix_table(p):
    """Return tables needed to apply good suffix rule."""
    n = comp_n(p)
    lp = big_l_prime_array(p, n)
    return lp, big_l_array(p, lp), small_l_prime_array(n)


def good_suffix_mismatch(i, big_l_prime, small_l_prime):
    """Given a mismatch at offset i, and given L/L' and l' arrays,
    return amount to shift as determined by good suffix rule."""
    length = len(big_l_prime)
    assert i < length

    if i == length - 1:
        return 0
    
    i += 1  # i points to leftmost matching position of P

    if big_l_prime[i] > 0:
        return length - big_l_prime[i]
    return length - small_l_prime[i]


def good_suffix_match(small_l_prime):
    """Given a full match of P to T, return amount to shift as
    determined by good suffix rule."""
    return len(small_l_prime) - small_l_prime[1]


def dense_bad_char_tab(p, amap):
    """Given pattern string and list with ordered alphabet characters, create
    and return a list of lists (shape len(p) x |alphabet|) where, for each index i, 
    entries contain latest appearance of a character in the query string."""
    tab = []
    nxt = [0] * len(amap)

    for i in range(0, len(p)):
        c = p[i]
        assert c in amap
        tab.append(nxt[:])
        nxt[amap[c]] = i+1
    return tab

def boyer_moore(p, p_bm, t):
    """ Do Boyer-Moore matching """
    i = 0
    occurrences = []
    comps = 0
    alignments = 0
    # Iterate over all starting points - determined by 'shift' calculated below
    while i < len(t) - len(p) + 1:
        shift = 1 # Assume base case of 1
        mismatched = False
        alignments += 1
        # Iterate over subsequence backwards
        for j in range(len(p)-1, -1, -1):
            comps += 1
            # If mismatch, compute skipped chars from bad char rule and good suffix rule
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                # update shift
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            # append start index for this subsequence
            occurrences.append(i)
            # compute skip based on index
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences, comps, alignments

class BoyerMoore(object):
    """Performs preprocessing for Boyer-Moore algorithm."""
    
    def __init__(self, p, alphabet='ACGT'):
        self.p = p
        self.alphabet = alphabet
        # Create map from alphabet characters to integers
        self.amap = {}
        for i in range(len(self.alphabet)):
            self.amap[self.alphabet[i]] = i
        # Make bad character rule table
        self.bad_char = dense_bad_char_tab(p, self.amap)
        # Create good suffix rule table
        _, self.big_l, self.small_l_prime = good_suffix_table(p)
    
    def bad_character_rule(self, i, c):
        """Returns number of skips given by bad character rule at index i."""
        assert c in self.amap
        ci = self.amap[c]
        assert i > (self.bad_char[i][ci]-1)
        return i - (self.bad_char[i][ci]-1)
    
    def good_suffix_rule(self, i):
        """Given a mismatch at index i, return amount to shift
        as determined by (weak) good suffix rule."""
        length = len(self.big_l)
        assert i < length
        if i == length - 1:
            return 0
        i += 1  # i points to leftmost matching position of P
        if self.big_l[i] > 0:
            return length - self.big_l[i]
        return length - self.small_l_prime[i]
    
    def match_skip(self):
        """Return amount to shift in case where P matches T."""
        return len(self.small_l_prime) - self.small_l_prime[1]
    
f = 'raw_files/chr1.GRCh38.excerpt.fasta'
t = readFasta(f) # just take the first read

# query string
p = 'GGCGCGGTGGCTCACGCCTGTAAT'
print(f"The query string is: {p}.")

ind = BoyerMoore(p, alphabet='ACGT')
occs, comps, alignments = boyer_moore(p, ind, t)
print(f"The number of occurrences of the pattern in the text is {len(occs)}.")
