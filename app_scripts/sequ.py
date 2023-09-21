"""
parse functions return None if the data is invalid for that format. load functions throw an error.
"""
# Author: Dave Curran
# Date: Aug 2022

import weakref # Used to avoid circular references between Sequence & SeqList
from collections import UserList


# TODO
# finish the phylip writing methods.
# could rename `make_unique()` to something like `filter_redundant()`; decide on some function nomenclature about what `remove_`, `fileter_`, or `make_` mean for function names and keep it standard.
# remove_shorter() takes quite a long time on 262k sequences, likely because of the `remove()` call. Provide a 'filter_shorter()' method as well, that just returns a new SeqList instead of removing.
#  - it's probably the _deregister calls that are slow. could maybe have an option to prevent them
# once this is decent, put a copy into a new molecbio.sequ file. eventually i need to replace it.
# translate, inverse, complements
# quality, aln_identity
# some tools for quantifying sequence lengths, possibly even printing a basic histogram. should use seq.nongaps for the lengths
# SeqList methods return the obj, Sequence methods return str(seq). Probably want to standardize them to both return the obj


# # #  Package file I/O functions
def save_fasta(seqlist, filename, line=None, spaces=False, numbers=False):
    """line is an int indicating the length of each line. spaces and numbers are booleans."""
    with open(filename, 'w') as f:
        f.write(seqlist.to_fasta(line=line, spaces=spaces, numbers=numbers))
def save_fasta_sequences(seqlist, filename, line=None, spaces=False, numbers=False):
    """Designed for alignments, removes gaps and saves as FASTA sequences without modifying the passed `seqlist`."""
    new_seqs = seqlist.copy()
    new_seqs.strip_gaps()
    with open(filename, 'w') as f:
        f.write(seqlist.to_fasta(line=line, spaces=spaces, numbers=numbers))
def save_clustal(seqlist, filename, numbers=True, name_len=None):
    with open(filename, 'w') as f:
        f.write(seqlist.to_clustal(numbers=numbers, name_len=name_len))
def save_phylip(seqlist, filename, kind='interleaved', strict=False, per_line=70, chunk_size=10):
    with open(filename, 'w') as f:
        f.write(seqlist.to_phylip(kind=kind, strict=strict, per_line=per_line, chunk_size=chunk_size))

def load(filename, only_these=None):
    """Attempts to parse the sequence file using all known formats."""
    with open(filename) as f:
        first_line = f.readline()
        f.seek(0)
        line_split = first_line.split()
        if first_line.upper().startswith(('CLUSTAL W', 'CLUSTALW')):
            seqs = parse_clustal(f, only_these)
            method = 'Clustal'
        elif len(line_split)==2 and line_split[0].isdigit() and line_split[1].isdigit():
            seqs = parse_phylip(f, only_these)
            method = 'Phylip'
        else:
            seqs = parse_fasta(f, only_these)
            method = 'FASTA'
    if seqs is None:
        raise MolecbioFileFormatError("could not load '{}' as a {}-format file of sequences.".format(filename, method))
    return seqs
def load_fasta(filename, only_these=None):
    """If given, only_these should be a tuple or collection of strings. Only sequences with names starting with at least one of those strings will be returned."""
    with open(filename) as f:
        seqs = parse_fasta(f, only_these)
    if seqs is None:
        raise MolecbioFileFormatError("could not load '{}' as a FASTA file.".format(filename))
    return seqs
def load_clustal(filename, only_these=None):
    with open(filename) as f:
        seqs = parse_clustal(f, only_these)
    if seqs is None:
        raise MolecbioFileFormatError("could not load '{}' as a Clustal alignment file.".format(filename))
    return seqs
def load_phylip(filename, only_these=None, kind='auto', strict=False):
    with open(filename) as f:
        seqs = parse_phylip(f, only_these, kind, strict)
    if seqs is None:
        raise MolecbioFileFormatError("could not load '{}' as a Phylip alignment file.".format(filename))
    return seqs

def parse_fasta(lines, only_these=None):
    if only_these:
        only_these = tuple(only_these)
    seqs = SeqList()
    name, seq_buff = None, []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line[0]=='>':
            if name:
                if not only_these or name.startswith(only_these):
                    seq = Sequence(name, ''.join(seq_buff))
                    seqs.append(seq)
            name = line[1:]
            seq_buff = []
        else:
            seq_buff.append(line)
    if name:
        seq = Sequence(name, ''.join(seq_buff))
        seqs.append(seq)
    if seqs:
        return seqs
    else:
        return None
def parse_clustal(lines, only_these=None):
    if only_these:
        only_these = tuple(only_these)
    for line in lines: # Ensures first line is right
        if not line.upper().startswith(('CLUSTAL W', 'CLUSTALW')):
            return None
        break
    if isinstance(lines, list):
        lines = lines[1:] # allows `lines` to be a file obj or actual list of strings.
    seq_dict = {}
    for line in lines:
        if line.startswith((' ','\t','\n')):
            continue # skips empty + conservation lines
        data = line.split()
        name, seq = data[:2]
        seq_dict.setdefault(name, []).append(seq)
    seqs = SeqList()
    for name, list_seq in seq_dict.items():
        if not only_these or name.startswith(only_these):
            seq = Sequence(name, ''.join(list_seq))
            seqs.append(seq)
    if seqs:
        return seqs
    else:
        return None
def parse_phylip(lines, only_these=None, kind='auto', strict=False):
    """`kind` can be 'auto', 'interleaved', 'i', 'sequential', or 's'. If `strict` is True, the format is from PhyML with name lengths of exactly 10. If False, just expects a space between the name and sequence."""
    if only_these:
        only_these = tuple(only_these)
    for line in lines: # Ensures first line is right
        data = line.split()
        if len(data) != 2 or not data[0].isdigit() or not data[1].isdigit():
            return None
        num_seqs = int(data[0])
        aln_len = int(data[1])
        break
    if isinstance(lines, list):
        lines = lines[1:]
    else:
        lines = lines.readlines()
    lines = [line.strip() for line in lines if line.strip()] # removes blank lines
    if len(lines) % num_seqs != 0:
        return None # Needs to be an equal # of lines associated with each sequence
    kind = kind.lower()
    if kind == 'auto':
        return parse_phylip_lines(lines, only_these, True, strict, num_seqs, aln_len) or parse_phylip_lines(lines, only_these, True, not strict, num_seqs, aln_len) or parse_phylip_lines(lines, only_these, False, strict, num_seqs, aln_len) or parse_phylip_lines(lines, only_these, False, not strict, num_seqs, aln_len)
    elif kind in ('interleaved', 'i'):
        return parse_phylip_lines(lines, only_these, True, strict, num_seqs, aln_len) or parse_phylip_lines(lines, only_these, True, not strict, num_seqs, aln_len)
    elif kind in ('sequential', 's'):
        return parse_phylip_lines(lines, only_these, False, strict, num_seqs, aln_len) or parse_phylip_lines(lines, only_these, False, not strict, num_seqs, aln_len)
    else:
        raise ValueError("unrecognized `kind` argument '{}'. It must be one of: 'auto', 'interleaved', 'i', 'sequential', or 's'.".format(kind))
def parse_phylip_lines(lines, only_these, is_interleaved, strict, num_seqs, aln_len):
    strict_len = 10 # To be compatible with PhyML and similar software
    seqs = SeqList()
    lines_per = int(len(lines) / num_seqs)
    for seq_ind in range(num_seqs):
        if is_interleaved:
            first_ind = seq_ind
        else:
            first_ind = seq_ind * lines_per
        if strict:
            name = lines[first_ind][:strict_len].strip()
            seq_buff = lines[first_ind][strict_len:].split()
        else:
            line_data = lines[first_ind].split()
            name = line_data[0]
            seq_buff = line_data[1:]
        if only_these and not name.startswith(only_these):
            continue
        for ind in range(1, lines_per):
            if is_interleaved:
                next_ind = first_ind + ind * num_seqs
            else:
                next_ind = first_ind + ind
            seq_buff.extend(lines[next_ind].split())
        seq = Sequence(name, sequence=''.join(seq_buff))
        if len(seq) != aln_len:
            return None
        seqs.append(seq)
    if len(seqs) != num_seqs:
        return None
    return seqs


# # #  Module constants
whitespace_name_filter = {ord(' '):'_', ord('\t'):'_', ord('\n'):'_'}
# Restricted characters for phylogenetic software. Use: name.translate(phylo_name_filter)
phylo_name_filter = {ord(' '):'_', ord('\t'):'_', ord('\n'):'_', ord(','):'_', ord(':'):'_', ord('('):None, ord(')'):None, ord('['):None, ord(']'):None, ord('<'):None, ord('>'):None, ord(';'):None, ord('='):None}


# # #  Module classes
class SeqList(UserList):
    """
    The manipulation functions return self, allowing them to be chained."""
    def __init__(self, sequences=[], **kwargs):
        super().__init__()
        # self.names = [] # Read-only property; list of sequence names in order.
        # self.size = int # Read-only property; sum total of all sequence lengths.
        # self.nongaps = int # Read-only property; sum of all alpha characters in sequences.
        # self.gaps = int # Read-only property; sum of all `-` characters in sequences.
        # self.is_alignment = Boolean # Read-only property; checks that all lengths are equal.
        # self.lengths = List # Read-only property; returns a list of all sequence nongap lengths in order.
        # #  Private attributes
        self.names_cache = {} # Reset when SeqList changes, or when a child Sequence name is changed
        self.self_ref = weakref.ref(self) # Used by Sequence children
        # #  Finish initialization
        for seq in sequences:
            self.append(seq)

    # # #  Accession and filtering
    def get(self, name, not_found=None):
        """Returns the most recently added Sequence object with that name. If not found, returns the `not_found` value. Builds a cache just-in-time if needed."""
        if not self.names_cache:
            self.rebuild_cache()
        return self.names_cache.get(name, not_found)
    def get_named(self, names):
        """Returns a new SeqList of all sequences identified with a name in the collection `names`, respecting the order of `names` if applicable. Unfound names are ignored."""
        if not self.names_cache:
            self.rebuild_cache()
        hits = [self.names_cache[name] for name in names if name in self.names_cache]
        return SeqList(sequences=hits)
    def get_where(self, selector):
        """Returns a new SeqList where the function `selector(seqobj)` evaluates to True."""
        matches = SeqList()
        for seq in self.data:
            if selector(seq):
                matches.append(seq)
        return matches
    def copy(self):
        """Returns a deepcopy of self."""
        seqs = SeqList()
        for seq in self.data:
            seqs.append(seq.copy())
        return seqs
    def rebuild_cache(self):
        cache = {}
        for seq in self.data:
            _ = cache.setdefault(seq.name, seq) # Does not overwrite repeated names
        self.names_cache = cache

    # #  Sequence manipulations
    def clean(self, allowed={'-'}):
        """Calls the clean() method for all child Sequences, which removes all numbers, whitespace, and symbols other than `-`."""
        for seq in self.data:
            seq.clean(allowed)
        return self
    def strip_gaps(self):
        """Removes all gap characters in the sequences."""
        for seq in self.data:
            seq.strip_gaps()
        return self
    def strip_spaces(self):
        """Removes all whitespace characters in the sequences."""
        for seq in self.data:
            seq.strip_spaces()
        return self
    def make_unique(self, return_removed=False):
        """Removes any children with identical sequences, keeping only the first object encountered. If `return_removed` is True, returns a SeqList of the removed Sequences."""
        seqset, removed = set(), []
        for seq in self.data[::]:
            if seq.seq in seqset:
                removed.append(seq)
                self.remove(seq)
            else:
                seqset.add(seq.seq)
        if return_removed:
            return SeqList(removed)
        else:
            return self
    def remove_empty(self, return_removed=False):
        """Removes all sequences with no non-gap characters. If `return_removed` is True, returns a SeqList of the removed Sequences."""
        return self.remove_shorter(1, return_removed=return_removed)
    def remove_shorter(self, size, return_removed=False):
        """Removes all sequences with non-gap length < `size`. If `return_removed` is True, returns a SeqList of the removed Sequences."""
        to_remove = []
        for seq in self.data:
            if seq.nongaps < size:
                to_remove.append(seq)
        for seq in to_remove:
            self.remove(seq)
        if return_removed:
            return SeqList(to_remove)
        else:
            return self
    def remove_longer(self, size, return_removed=False):
        """Removes all sequences with non-gap length > `size`. If `return_removed` is True, returns a SeqList of the removed Sequences."""
        to_remove = []
        for seq in self.data:
            if seq.nongaps > size:
                to_remove.append(seq)
        for seq in to_remove:
            self.remove(seq)
        if return_removed:
            return SeqList(to_remove)
        else:
            return self
    def trim_to(self, start=None, end=None):
        """For all Sequences, keeps only the sequence between indices `start` and `end`."""
        for seq in self.data:
            seq.trim_to(start, end)
        return self

    # #  Sequence name manipulations
    def make_names_unique(self, pattern='{}_{}'):
        """Changes the name of child Sequences so that all are unique, by adding sequential integers to any that are identical. If given, `pattern` should be amenable to `pattern.format(old_name, int)`."""
        # return self
        pass
    def clean_names(self, replacements=phylo_name_filter):
        """Removes or replaces characters in sequence `name`s. The default is to be compatible with various phylogenetic software."""
        for seq in self.data:
            seq.clean_name(replacements)
        return self

    # #  Output formats
    def to_fasta(self, line=None, spaces=False, numbers=False):
        """line is an int indicating the length of each line. spaces and numbers are booleans."""
        if line:
            return '\n\n'.join(seq.to_fasta(line, spaces, numbers) for seq in self.data)
        else:
            return '\n\n'.join(str(seq) for seq in self.data)
    def to_clustal(self, numbers=True, name_len=None):
        """Will not work if sequences have different lengths. If given, `name_len` should be an int indicating the max length of name to include. Does not check for uniqueness."""
        strong = (set('STA'), set('NEQK'), set('NHQK'), set('NDEQ'), set('QHRK'),
                  set('MILV'), set('MILF'), set('HY'), set('FYW'))
        weak = (set('CSA'), set('ATV'), set('SAG'), set('STNK'), set('STPA'),
                set('SGND'), set('SNDEQK'), set('NDEQHK'), set('NEQHRK'),
                set('FVLIM'), set('HFY'))
        def _conservation(col_tup):
            col = set(col_tup)
            if '-' in col: return ' '
            if len(col) == 1: return '*'
            for aas in strong:
                if col <= aas: return ':'
            for aas in weak:
                if col <= aas: return '.'
            else:
                return ' '
        names, sequences, seq_nums, seq_lens = [], [], {}, set()
        for seq in self.data:
            name = seq.name.translate(whitespace_name_filter)[:name_len]
            names.append(name)
            sequences.append(seq.seq.upper())
            seq_nums[name] = 0
            seq_lens.add(len(seq))
        if len(seq_lens) != 1:
            raise MolecbioInvalidAlignmentError("cannot save as a Clustal alignment as sequences have different lengths.")
        seq_len = seq_lens.pop()
        max_name = max(len(name) for name in names)
        name_fmt = '{{:<{}}}'.format(max_name)
        cons_pref = ' ' * max_name
        conserv = ''.join(_conservation(col) for col in zip(*sequences))
        buff = ['CLUSTAL W multiple sequence alignment', '\n']
        start_ind = 0
        for end_ind in range(60, seq_len+60, 60):
            for name, sequence in zip(names, sequences):
                line_buff = [name_fmt.format(name)]
                seq = sequence[start_ind:end_ind]
                line_buff.append(seq)
                if numbers:
                    num_chars = len(list(filter(str.isalpha, seq)))
                    seq_nums[name] += num_chars
                    line_buff.append(str(seq_nums[name]))
                buff.append(' '.join(line_buff))
            buff.append('{} {}'.format(cons_pref, conserv[start_ind:end_ind]))
            buff.append('')
            start_ind = end_ind
        return '\n'.join(buff)
    def to_phylip(self, kind='interleaved', strict=False, per_line=70, chunk_size=10):
        """If given, `name_len` should be an int indicating the max length of name to include. Does not check for uniqueness. `line` indicates approximately how long to make each line."""
        # if strict=True, seq names are limited to 10 characters, and edits seq names to be compatible with phylogenetic software, so removes/replaces: ` ()[]<>,:;=`

        per_line = 32
        #chunk_size = 7

        names, aln_seqs, aln_lens = [], [], set()
        for seq in self.data:
            if strict:
                name = '{:10s}'.format(seq.name.translate(phylo_name_filter)[:10])
            else:
                name = seq.name.translate(whitespace_name_filter)
            names.append(name)
            aln_seq = ''.join(seq.seq.split())
            aln_seqs.append(aln_seq)
            aln_lens.add(len(aln_seq))
        if len(aln_lens) != 1:
            raise MolecbioInvalidAlignmentError("cannot save as a Phylip alignment as sequences have different lengths.")
        aln_len = aln_lens.pop()
        max_name = max(len(name) for name in names)
        buff = ['   {}  {}'.format(len(self), aln_len)]
        if kind == 'interleaved':
            pass
        elif kind == 'sequential':
            for name, seq in zip(names, aln_seqs):
                line_buff, line_count = [], 0
                for start in range(0, aln_len, chunk_size):
                    seq_str = seq[start:start+chunk_size]
                    if start == 0:
                        if strict:
                            line_buff = ['{}{}'.format(name, seq_str)]
                            line_count = len(line_buff[0])
                        else:
                            line_buff = ['{} {}'.format(name, seq_str)]
                            line_count = max_name + len(seq_str) + 1
                    else:
                        line_buff.append(seq_str)
                        line_count += len(seq_str)
                    line_total = line_count + len(line_buff) - 1
                    next_chunk_len = min(aln_len - start - chunk_size, chunk_size)
                    if line_total + next_chunk_len + 1 > per_line:
                        buff.append(' '.join(line_buff))
                        line_buff, line_count = [], 0
                if line_buff:
                    buff.append(' '.join(line_buff))
        else:
            raise ValueError("unrecognized `kind` argument '{}'. It must be one of: 'interleaved', 'i', 'sequential', or 's'.".format(kind))
        return '\n'.join(buff)

    # #  Dunder / built-in implementations
    def append(self, seqobj):
        self.data.append(seqobj)
        self._register_seq(seqobj)
    def insert(self, index, seqobj):
        self.data.insert(index, seqobj)
        self._register_seq(seqobj)
    def extend(self, new_seqs):
        for seq in new_seqs:
            self.data.append(seq)
            self._register_seq(seq)
    def pop(self, index=-1):
        seq = self.data.pop(index)
        self._deregister_seq(seq)
        return seq
    def remove(self, seqobj):
        self.data.remove(seqobj)
        self._deregister_seq(seqobj)
    def __setitem__(self, index, seqobj):
        if isinstance(index, slice):
            for ind, newseq in zip(range(*index.indices(len(self))), seqobj):
                self._deregister_seq(self.data[ind])
                self._register_seq(newseq)
                self.data[ind] = newseq
        else:
            self._deregister_seq(self.data[index])
            self._register_seq(seqobj)
            self.data[index] = seqobj
    def __delitem__(self, index):
        if isinstance(index, slice):
            for ind in range(*index.indices(len(self))):
                self._deregister_seq(self.data[ind])
        else:
            self._deregister_seq(self.data[index])
        del self.data[index]
    def __add__(self, seqlist):
        if not isinstance(seqlist, SeqList):
            return NotImplemented
        return SeqList(sequences=self.data+seqlist.data)
    def __radd__(self, *args):
        return NotImplemented
    def __iadd__(self, seqlist):
        if not isinstance(seqlist, SeqList):
            return NotImplemented
        for seq in seqlist:
            self.append(seq)
        return self
    def __mul__(self, num):
        if not isinstance(num, int):
            return NotImplemented
        if num < 1:
            return NotImplemented
        return SeqList(sequences=self.data*num)
    def __rmul__(self, num):
        if not isinstance(num, int):
            return NotImplemented
        if num < 1:
            return NotImplemented
        return SeqList(sequences=self.data*num)
    def __imul__(self, num):
        if not isinstance(num, int):
            return NotImplemented
        if num < 1:
            return NotImplemented
        elif num == 1:
            return self
        else:
            for seq in self.data * (num-1):
                self.append(seq)
            return self
    def __del__(self):
        for seq in self.data:
            self._deregister_seq(seq)
    def __eq__(self, other):
        """Only valid comparison is to other SeqList objects."""
        if not isinstance(other, SeqList):
            return False
        if len(other) != len(self.data):
            return False
        for seq1, seq2 in zip(self.data, other.data):
            if seq1 != seq2:
                return False
        return True
    def __str__(self):
        return 'SeqList with {} sequences'.format(len(self))
    def __repr__(self):
        return '<SeqList at {} of length {}>'.format(hex(id(self)).upper(), len(self))

    # #  Private methods
    def _register_seq(self, seqobj):
        seqobj.parent_refs.append(self.self_ref)
        if self.names_cache:
            self.names_cache = {}
    def _deregister_seq(self, seqobj):
        seqobj.parent_refs.remove(self.self_ref)
        if self.names_cache:
            self.names_cache = {}

    # #  Properties
    @property
    def names(self):
        return [seq.name for seq in self.data]
    @property
    def size(self):
        return sum(len(seq.seq) for seq in self.data)
    @property
    def nongaps(self):
        return sum(seq.nongaps for seq in self.data)
    @property
    def gaps(self):
        return sum(seq.gaps for seq in self.data)
    @property
    def is_alignment(self):
        lens = {len(seq) for seq in self.data}
        return len(lens) == 1
    @property
    def lengths(self):
        return [seq.nongaps for seq in self.data]


class Sequence():
    """

    Implements the following built-in methods: len()"""
    def __init__(self, name='', sequence=''):
        # #  Public attributes
        self.seq = sequence # Also available as Sequence.sequence
        # self.name = str # Property
        # self.nongaps = int # Read-only property; sum of all alpha characters.
        # self.gaps = int # Read-only property; sum of all `-` characters.
        # #  Private attributes
        self._name = name
        self.parent_refs = [] # The SeqList objects that are holding self
    # #  Sequence manipulations
    def strip_gaps(self):
        """Removes all gap characters in the sequence."""
        self.seq = self.seq.replace('-','')
        return self.seq
    def strip_spaces(self):
        """Removes all whitespace characters in the sequences."""
        self.seq = ''.join(self.seq.split())
        return self.seq
    def trim_to(self, start=None, end=None):
        """Keeps only the sequence between indices `start` and `end`."""
        self.seq = self.seq[start:end]
        return self.seq
    def upper(self):
        """Sets the sequence to upper case, returning it."""
        self.seq = self.seq.upper()
        return self.seq
    def lower(self):
        """Sets the sequence to lower case, returning it."""
        self.seq = self.seq.lower()
        return self.seq
    # #  Name manipulations
    def clean(self, allowed={'-'}):
        """Removes all numbers, whitespace, and symbols other than -. Also converts to upper case."""
        new_seq = []
        for c in self.seq:
            if c.isalpha() or c in allowed:
                new_seq.append(c.upper())
        self.seq = ''.join(new_seq)
    def clean_name(self, replacements=phylo_name_filter):
        """Removes or replaces characters in `name`. The default is to be compatible with various phylogenetic software."""
        self.name = self.name.translate(replacements)
    # #  Outputs
    def copy(self):
        return Sequence(name=self.name, sequence=self.seq)
    def to_fasta(self, line=70, spaces=False, numbers=False):
        """`line` is an int indicating the characters per line of output, `spaces` and `numbers` are booleans. If `line` is False-y, no linebreaks or spaces or numbers will be used."""
        spaces_int = 10
        if not line:
            return str(self)
        else:
            if numbers:
                num_max = len(self.seq) // line * line
                num_len = len(str(num_max))
                num_fmt = '{{:>{}}}'.format(num_len)
            buff = []
            start_ind = 0
            for end_ind in range(line, len(self)+line, line):
                seq_line = self.seq[start_ind:end_ind]
                line_buff = []
                if numbers:
                    line_buff.append(num_fmt.format(start_ind+1))
                if spaces:
                    s_i = 0
                    for e_i in range(spaces_int, len(seq_line)+spaces_int, spaces_int):
                        line_buff.append(seq_line[s_i:e_i])
                        s_i = e_i
                else:
                    line_buff.append(seq_line)
                buff.append(' '.join(line_buff))
                start_ind = end_ind
            fasta_seq = '\n'.join(buff)
            return '>{}\n{}'.format(self.name, fasta_seq)
    # #  Dunder implementations
    def __str__(self):
        """Called by str(Sequence) or print(Sequence)."""
        return '>{}\n{}'.format(self.name, self.seq)
    def __repr__(self):
        """Called by repr(Sequence)."""
        return "<Sequence at {} length {} named '{}'>".format(hex(id(self)).upper(), len(self), self.name)
    def __len__(self):
        return len(self.seq)
    def __getitem__(self, index):
        return self.seq[index]
    def __eq__(self, other):
        """Valid comparisons are against other Sequences or against a string representing a sequence."""
        if isinstance(other, str):
            return self.seq == other
        if not isinstance(other, Sequence):
            return False
        return self.name == other.name and self.seq == other.seq
    def __contains__(self, other):
        if isinstance(other, Sequence):
            return other.seq in self.seq
        return other in self.seq
    # #  Properties
    @property
    def sequence(self):
        return self.seq
    @sequence.setter
    def sequence(self, value):
        self.seq = value
    @property
    def name(self):
        return self._name
    @name.setter
    def name(self, new_name):
        # Updates all containing SeqLists that `name` has changed.
        self._name = new_name
        for p_ref in self.parent_refs:
            parent = p_ref() # weakref
            if parent is not None and parent.names_cache:
                parent.names_cache = {}
    @property
    def nongaps(self):
        # Does not validate characters, just counts all alphabetic characters.
        return len(list(filter(str.isalpha, self.seq)))
    @property
    def gaps(self):
        return self.seq.count('-')


# # #  Errors
class MolecbioFileFormatError(ValueError):
    """Raised when loading a file of the incorrect format."""
class MolecbioInvalidAlignmentError(ValueError):
    """Raised when a user attempts to save as an alignment file sequences that do not constitute a valid alignment."""


# # #  TESTING
if __name__ == '__main__':
    seqs = load_phylip('test_dir\\test_sequential.phylip', kind='s')
    #seqs = load('test_dir\confirmed.aln')
    #seqs.strip_gaps().make_unique().clean_names()
    #save_fasta(seqs, 'test_dir\\test.aln')
    print(len(seqs))
    print(len(seqs[0]))
    print(len(seqs[0][0]))
