# -*- coding: utf-8 -*-

from __future__ import division

'''Utilities to manipulate Compact Idiosyncratic Gapped Alignment Report (CIGAR) strings'''



from libc.stdio                cimport sprintf
from libc.stdlib               cimport malloc, realloc, free
from libc.string               cimport memcpy
from cpython.ref               cimport PyObject
from cpython                   cimport PyBytes_FromStringAndSize
from pysam.libchtslib          cimport bam_get_cigar, bam1_t


cdef extern from "Python.h":
    ctypedef struct PySliceObject
    int PySlice_GetIndicesEx(PySliceObject *slice, Py_ssize_t length, Py_ssize_t *start, Py_ssize_t *stop,
                                                   Py_ssize_t *step, Py_ssize_t *slicelen) except -1


DEF BAM_CIGAR_SHIFT=4
DEF BAM_CIGAR_MASK=((1 << BAM_CIGAR_SHIFT) - 1)


cdef class CigarOperator(object):
    # see pxd for attribute defs

    def __cinit__(self, str name, str code, char binary, char consumes_read_bases, char consumes_ref_bases):
        if len(code) != 1:
            raise ValueError('only single character CIGAR codes allowed')

        self.name                = name
        self.code                = ord(code[0])
        self.codestr             = code
        self.binary              = binary
        self.consumes_read_bases = consumes_read_bases
        self.consumes_ref_bases  = consumes_ref_bases
        self.inverse             = None

    def __repr__(self):
        return '<CigarOperator:%s>' % self.name


MATCH        = CigarOperator('MATCH',        'M', BAM_MATCH,        1, 1)
INSERTION    = CigarOperator('INSERTION',    'I', BAM_INSERTION,    1, 0)
DELETION     = CigarOperator('DELETION',     'D', BAM_DELETION,     0, 1)
SKIPPED      = CigarOperator('SKIPPED',      'N', BAM_SKIPPED,      0, 1)
SOFT_CLIP    = CigarOperator('SOFT_CLIP',    'S', BAM_SOFT_CLIP,    1, 0)
HARD_CLIP    = CigarOperator('HARD_CLIP',    'H', BAM_HARD_CLIP,    0, 0)
PADDING      = CigarOperator('PADDING',      'P', BAM_PADDING,      0, 0)
SEQ_MATCH    = CigarOperator('SEQ_MATCH',    '=', BAM_SEQ_MATCH,    1, 1)
SEQ_MISMATCH = CigarOperator('SEQ_MISMATCH', 'X', BAM_SEQ_MISMATCH, 1, 1)

MATCH.inverse         = MATCH
INSERTION.inverse     = DELETION
DELETION.inverse      = INSERTION
SEQ_MATCH.inverse     = SEQ_MATCH
SEQ_MISMATCH.inverse  = SEQ_MISMATCH

cpdef tuple BAM_OPS = (MATCH, INSERTION, DELETION, SKIPPED, SOFT_CLIP, HARD_CLIP, PADDING, SEQ_MATCH, SEQ_MISMATCH)

# Create mapping from ord(char_op) to op
_temp_map = [None]*256
for _op in BAM_OPS:
    _temp_map[_op.code] = _op

cpdef tuple BAM_OP_CHAR_MAP = tuple(_temp_map)
del _temp_map, _op


cdef inline CigarOperator get_op_by_char(char char_op):
    '''
    Obtain the (singleton) CIGAR operator object corresponding to a character.

    Args:
        char_op: Character (e.g. 'M', 'I', 'D') associated with CIGAR operator.

    Returns:
        The CIGAR operator object corresponding to char_op.
    '''
    op = BAM_OP_CHAR_MAP[char_op]
    if op is None:
        raise KeyError('cigar op=%s (%d) not found' % (chr(char_op),char_op))
    return op


cpdef inline CigarOperator get_op_by_str(bytes str_op):
    '''
    Obtain the (singleton) CIGAR operator object corresponding to a character.

    Args:
        char_op: Character (e.g. 'M', 'I', 'D') associated with CIGAR operator.

    Returns:
        The CIGAR operator object corresponding to char_op.
    '''
    if len(str_op) != 1:
        raise KeyError('invalid cigar op=%s' % str_op)
    cdef char_op = ord(str_op)
    op = BAM_OP_CHAR_MAP[char_op]
    if op is None:
        raise KeyError('cigar op=%s (%d) not found' % (chr(char_op),char_op))
    return op


cpdef inline CigarOperator get_op_by_binary(char binary_rep):
    '''
    Obtain the (singleton) CIGAR operator object corresponding to a number.

    Args:
        binary_rep: Number associated with the CIGAR operator

    Returns:
        The CIGAR operator object corresponding to binary_rep.
    '''
    return BAM_OPS[binary_rep]


cdef inline int isdigit(char c):
    return (c >= 48) & (c < 58)  # (c >= '0') & (c <= '9')


cdef inline int integer_width(uint32_t i):
    if i < 10:
        return 1
    elif i < 100:
        return 2
    elif i < 1000:
        return 3
    elif i < 10000:
        return 4
    elif i < 100000:
        return 5
    else:
        return 5 + integer_width(i // 100000)


cdef inline _decode_cigar_to_op(uint32_t c):
    op = get_op_by_binary(c & BAM_CIGAR_MASK)
    n  = c >> BAM_CIGAR_SHIFT
    return op, n


cdef inline _decode_cigar_to_pysam(uint32_t c):
    op = c &  BAM_CIGAR_MASK
    n  = c >> BAM_CIGAR_SHIFT
    return op, n


cdef inline CigarOperator _decode_op(op_obj):
    if isinstance(op_obj, CigarOperator):
        return <CigarOperator>op_obj
    elif isinstance(op_obj, int):
        return get_op_by_binary(op_obj)
    elif isinstance(op_obj, bytes):
        return get_op_by_str(op_obj)
    else:
        raise TypeError(f"_decode_op failed with {op_obj} ({type(op_obj)})")


cdef class CigarSequence(object):
    '''
    Object representing a sequence of CIGAR operators and counts.

    CigarSequence can be instantiated with any of the following:

            1. CIGAR from an aligned read (pysam.AlignedRead)
            2. another CigarSequence
            3. sequence of tuples of CigarOperator and count
            4. sequence of tuples of BAM operator and count (pysam.AlignedRead.cigar)
            5. sequence of tuples of CIGAR operator character and count
            6. CIGAR string

    '''
    # see pxd for attribute defs

    def __cinit__(self, src=None):
        self.cigar_buf     = <uint32_t*>0
        self.cigar_len     = 0
        self.cigar_buf_len = 0

        if src is not None:
            self.extend(src)

    cdef _extend_from_binary(self, uint32_t other_cigar_len, uint32_t* other_cigar_buf):
        cdef uint32_t  last_op, first_op, n

        # Ensure last element of the existing CIGAR and first element of the
        # new CIGAR are properly normalized
        if self.cigar_len and other_cigar_len:
            last_op  = self.cigar_buf[self.cigar_len - 1] & BAM_CIGAR_MASK
            first_op = other_cigar_buf[0] & BAM_CIGAR_MASK
            if last_op == first_op:
                # Compute new length
                n  = self.cigar_buf[self.cigar_len - 1] >> BAM_CIGAR_SHIFT
                n += other_cigar_buf[0] >> BAM_CIGAR_SHIFT

                # Update last element
                self.cigar_buf[self.cigar_len - 1] = (n << BAM_CIGAR_SHIFT) | last_op

                # Do not copy first new element
                other_cigar_len -= 1
                other_cigar_buf += 1

        self._alloc_buf(self.cigar_len + other_cigar_len)
        memcpy(self.cigar_buf + self.cigar_len, other_cigar_buf, other_cigar_len * sizeof(uint32_t))
        self.cigar_len += other_cigar_len

    cdef _extend_from_align(self, AlignedRead align):
        cdef bam1_t* src = align._delegate
        self._extend_from_binary(src.core.n_cigar, bam_get_cigar(src))

    cdef _extend_from_cigar_sequence(self, CigarSequence other):
        self._extend_from_binary(other.cigar_len, other.cigar_buf)

    cdef _extend_from_cigar_slice(self, CigarSequence other, slice sliceobj):
        cdef Py_ssize_t start, stop, step, slicelen
        cdef Py_ssize_t i, j

        PySlice_GetIndicesEx(<PySliceObject*>sliceobj, other.cigar_len, &start, &stop, &step, &slicelen)

        self._alloc_buf(self.cigar_len + slicelen)

        j = start
        for i in range(slicelen):
            self.cigar_buf[i] = other.cigar_buf[j]
            j += step
        self.cigar_len = slicelen

    cdef _extend_from_str(self, bytes cigar_str):
        cdef int  digits = 0
        cdef int  op, n, new_op, new_n
        cdef char c

        if self.cigar_len:
            self.cigar_len -= 1
            op = self.cigar_buf[self.cigar_len] &  BAM_CIGAR_MASK
            n  = self.cigar_buf[self.cigar_len] >> BAM_CIGAR_SHIFT
        else:
            op = -1
            n  = 0

        for c in cigar_str:
            if isdigit(c):
                digits = 10*digits + (c-48)  # ord('0')==48
            else:
                new_op = get_op_by_char(c).binary
                new_n  = digits or 1
                digits = 0

                if new_op != op and n:
                    self._append_op(op, n)
                    n = 0

                op  = new_op
                n  += new_n

        if digits:
            raise ValueError('Invalid CIGAR string: trailing digits with no operation')

        if n:
            self._append_op(op, n)

    cdef _extend_from_sequence(self, cigar_seq):
        cdef int op, new_op, n, new_n

        if self.cigar_len:
            self.cigar_len -= 1
            op = self.cigar_buf[self.cigar_len] &  BAM_CIGAR_MASK
            n  = self.cigar_buf[self.cigar_len] >> BAM_CIGAR_SHIFT
        else:
            op = -1
            n  = 0

        for op_obj,new_n in cigar_seq:
            new_op = _decode_op(op_obj).binary

            if op != new_op and n:
                self._append_op(op, n)
                n = 0

            op  = new_op
            n  += new_n

        if n:
            self._append_op(op, n)

    def __dealloc__(self):
        self.cigar_len = 0
        self.cigar_buf_len = 0
        if self.cigar_buf:
            free(self.cigar_buf)
            self.cigar_buf = <uint32_t*>0

    cdef _alloc_buf(self, uint32_t n):
        if n < 4:
            n = 4
        if self.cigar_buf_len < n:
            if self.cigar_buf:
                self.cigar_buf = <uint32_t*>realloc(self.cigar_buf, n * sizeof(uint32_t))
            elif n:
                self.cigar_buf = <uint32_t*>malloc(n * sizeof(uint32_t))

            if not self.cigar_buf:
                self.cigar_len     = 0
                self.cigar_buf_len = 0
                raise MemoryError

            self.cigar_buf_len = n

    cdef inline _append_op(self, char op, uint32_t n):
        cdef uint32_t last_op

        if n == 0:
            return

        if self.cigar_len:
            last_op = self.cigar_buf[self.cigar_len - 1] & BAM_CIGAR_MASK
            if op == last_op:
                n += self.cigar_buf[self.cigar_len - 1] >> BAM_CIGAR_SHIFT
                self.cigar_buf[self.cigar_len - 1] = (n << BAM_CIGAR_SHIFT) | op
                return

        if self.cigar_buf_len == self.cigar_len:
            self._alloc_buf(self.cigar_buf_len * 2)

        self.cigar_buf[self.cigar_len] = (n << BAM_CIGAR_SHIFT) | op
        self.cigar_len += 1

    def __len__(self):
        '''
        Returns the number of CIGAR operators
        '''
        return self.cigar_len

    def __contains__(self, op_obj):
        '''
        Returns True if op is one of the CIGAR operators present; False otherwise.
        '''
        cdef int i, op = _decode_op(op_obj).binary
        for i in range(self.cigar_len):
            if (self.cigar_buf[i] & BAM_CIGAR_MASK) == op:
                return True
        return False

    cpdef int count(self, op_obj):
        '''
        Returns the total count for the specified CIGAR op.

        >>> cigar = CigarSequence(b'6H5S4M3I5M2D6S11H')
        >>> cigar.count(b'H')
        17
        >>> cigar.count(SOFT_CLIP)
        11
        >>> cigar.count(b'M')
        9
        >>> cigar.count(b'X')
        0
        '''
        cdef int i, n = 0, op = _decode_op(op_obj).binary
        for i in range(self.cigar_len):
            if (self.cigar_buf[i] & BAM_CIGAR_MASK) == op:
                n += self.cigar_buf[i] >> BAM_CIGAR_SHIFT
        return n

    cpdef int gapped_len(self, int include_soft_clip=0):
        '''
        Return the gapped aligment length with our without considering soft clipped bases

        >>> CigarSequence(b'6H5S4M3I5M2D6S11H').gapped_len()
        14
        >>> CigarSequence(b'6H5S4M3I5M2D6S11H').gapped_len(True)
        25
        '''
        cdef int           i
        cdef uint32_t      n = 0
        cdef CigarOperator op

        for i in range(self.cigar_len):
            op = get_op_by_binary(self.cigar_buf[i] & BAM_CIGAR_MASK)
            if not include_soft_clip and op is SOFT_CLIP:
                continue
            if op.consumes_ref_bases or op.consumes_read_bases:
                n += self.cigar_buf[i] >> BAM_CIGAR_SHIFT
        return n

    cpdef int ref_len(self, int query_bases=-1):
        '''
        Return the reference length.  If query_bases is greater than zero, then the query length
        returned is based on that number of reference bases.

        >>> CigarSequence(b'6H5S4M3I5M2D6S11H').ref_len()
        11
        >>> CigarSequence(b'6H5S4M3I5M2D6S11H').ref_len(0)
        0
        '''
        cdef int           i, done = 0
        cdef uint32_t      n, ref = 0, query = 0
        cdef CigarOperator op

        for i in range(self.cigar_len):
            op = get_op_by_binary(self.cigar_buf[i] & BAM_CIGAR_MASK)
            n = self.cigar_buf[i] >> BAM_CIGAR_SHIFT

            if op.consumes_read_bases:
                if query + n >= query_bases >= 0:
                    n = query_bases - query
                    done = 1
                query += n

            if op.consumes_ref_bases:
                ref += n

            if done:
                break

        return ref

    cpdef int query_len(self, int ref_bases=-1, int include_soft_clip=0):
        '''
        Return the query or read length with or without considering soft clipped bases.  If ref_bases
        is greater than zero, then the query length returned is based on that number of reference
        bases.

        >>> CigarSequence(b'6H5S4M3I5M2D6S11H').query_len()
        12
        >>> CigarSequence(b'6H5S4M3I5M2D6S11H').query_len(include_soft_clip=True)
        23
        >>> CigarSequence(b'6H5S4M3I5M2D6S11H').query_len(5)
        8
        >>> CigarSequence(b'6H5S4M3I5M2D6S11H').query_len(5, True)
        13
        >>> CigarSequence(b'6H5S4M3I5M2D6S11H').query_len(0)
        0
        >>> CigarSequence(b'6H5S4M3I5M2D6S11H').query_len(0, True)
        5
        '''
        cdef int           i, done = 0
        cdef uint32_t      n, ref = 0, query = 0
        cdef CigarOperator op

        for i in range(self.cigar_len):
            op = get_op_by_binary(self.cigar_buf[i] & BAM_CIGAR_MASK)
            n = self.cigar_buf[i] >> BAM_CIGAR_SHIFT

            if not include_soft_clip and op is SOFT_CLIP:
                continue

            if op.consumes_ref_bases:
                if ref + n >= ref_bases >= 0:
                    n = ref_bases - ref
                    done = 1
                ref += n

            if op.consumes_read_bases:
                query += n

            if done:
                break

        return query

    def __getitem__(self, index):
        if isinstance(index, slice):
            return self._get_slice(index)

        cdef Py_ssize_t i = index

        if i < 0:
            i += self.cigar_len

        if i < 0 or i >= self.cigar_len:
            raise IndexError

        return _decode_cigar_to_op(self.cigar_buf[i])

    cdef _get_slice(self, sliceobj):
        cdef CigarSequence new_cigar = CigarSequence()
        new_cigar._extend_from_cigar_slice(self, sliceobj)
        return new_cigar

    def __repr__(self):
        '''
        >>> repr(CigarSequence(b'10M10I10D'))
        '[(<CigarOperator:MATCH>, 10), (<CigarOperator:INSERTION>, 10), (<CigarOperator:DELETION>, 10)]'
        '''
        return str(self.to_op_list())

    cpdef append(self, item):
        '''
        Append to a CigarSequence a tuple of a CIGAR operator and a count

        >>> c = CigarSequence()
        >>> c.append( (b'I', 5) )
        >>> c
        [(<CigarOperator:INSERTION>, 5)]
        >>> c.append( (HARD_CLIP, 2) )
        >>> c
        [(<CigarOperator:INSERTION>, 5), (<CigarOperator:HARD_CLIP>, 2)]
        >>> c.append( (0, 5) )
        >>> c.to_string()
        b'5I2H5M'
        '''
        cdef uint32_t op
        op_obj, n = item
        op = _decode_op(op_obj).binary
        self._append_op(op, n)

    cpdef extend(self, items):
        '''
        Extend current CigarSequence with one of the following:

            1. CIGAR from an aligned read (pysam.AlignedRead)
            2. another CigarSequence
            3. sequence of tuples of CigarOperator and count
            4. sequence of tuples of BAM operator and count (pysam.AlignedRead.cigar)
            5. sequence of tuples of CIGAR operator character and count
            6. CIGAR string

        >>> c = CigarSequence(b'I510M')
        >>> c.extend(b'10M5I')
        >>> c.to_string()
        b'1I520M5I'
        >>> c.extend( [(1,5), (0, 10)] )
        >>> c.to_string()
        b'1I520M10I10M'
        >>> c.extend( [(SOFT_CLIP, 5)] )
        >>> c.to_string()
        b'1I520M10I10M5S'
        '''
        if isinstance(items, AlignedRead):
            self._extend_from_align(items)
        elif isinstance(items, CigarSequence):
            self._extend_from_cigar_sequence(items)
        elif isinstance(items, bytes):
            self._extend_from_str(items)
        else:
            self._extend_from_sequence(items)

    def copy(self):
        '''
        Return a new copy of a CigarSequence
        '''
        return CigarSequence(self)

    def __add__(self, other):
        '''
        Return concatination of a CigarSequence and an object that can be
        instantiated as a CigarSequence
        '''
        cdef CigarSequence c = CigarSequence(self)
        c.extend(other)
        return c

    def __iadd__(self, other):
        '''
        In-place concatination (extend) between a CigarSequence and object that can be
        instantiated as a CigarSequence
        '''
        self.extend(other)
        return self

    cpdef pop(self):
        '''
        Pop last CIGAR element, returning a tuple of CigarOperator and count

        >>> c = CigarSequence(b'10M10I10D')
        >>> c.pop()
        (<CigarOperator:DELETION>, 10)
        >>> c.pop()
        (<CigarOperator:INSERTION>, 10)
        >>> c.pop()
        (<CigarOperator:MATCH>, 10)
        >>> c.pop()
        Traceback (most recent call last):
        ...
        IndexError
        '''
        if not self.cigar_len:
            raise IndexError

        self.cigar_len -= 1
        return _decode_cigar_to_op(self.cigar_buf[self.cigar_len])

    cpdef reverse(self):
        '''
        Reverse CigarSequence in-place

        >>> c = CigarSequence(b'10M10I')
        >>> c.reverse()
        >>> c.to_string()
        b'10I10M'
        '''
        cdef int i = 0
        cdef int j = self.cigar_len - 1
        cdef uint32_t temp

        while i < j:
            temp = self.cigar_buf[i]
            self.cigar_buf[i] = self.cigar_buf[j]
            self.cigar_buf[j] = temp
            i += 1
            j -= 1

    cpdef invert(self, left_clip=None, right_clip=None):
        '''
        Inverts the reference and query sequence.
        :return: New CigarSequence with the cigar operators inverted and soft and hard clips stripped.
                 Will add to the cigar the specified number left and right soft clips.
                 Returns also the number of soft clips stripped from the left and right.
        '''
        cdef CigarSequence inv_cigar = CigarSequence()
        cdef CigarOperator op
        cdef uint32_t i, n, first_clip = 1
        cdef uint32_t s_left = 0
        cdef uint32_t s_right = 0
        if left_clip:
            if left_clip < 0:
                raise ValueError("left_clip value must be > 0")
            inv_cigar._append_op(SOFT_CLIP.binary, left_clip)
        for op, n in self:
            if op is HARD_CLIP:
                continue
            if op is SOFT_CLIP:
                if first_clip:
                    s_left = n
                else:
                    s_right = n
                continue
            first_clip = 0
            if op.inverse is None:
                raise KeyError('cigar op=%s has no inverse' % op.codestr)
            inv_cigar._append_op(op.inverse.binary, n)
        if right_clip:
            if right_clip < 0:
                raise ValueError("right_clip value must be > 0")
            inv_cigar._append_op(SOFT_CLIP.binary, right_clip)
        return inv_cigar, s_left, s_right

    cpdef convert_N_to_S(self):
        for i in range(self.cigar_len):
            op_binary = self.cigar_buf[i] & BAM_CIGAR_MASK
            if op_binary == SKIPPED.binary:
                n = self.cigar_buf[i] >> BAM_CIGAR_SHIFT
                self.cigar_buf[i] = (n << BAM_CIGAR_SHIFT) | SOFT_CLIP.binary

    def to_string(self):
        '''
        Convert CigarSequence to a CIGAR string

        >>> CigarSequence(b'150M3S5H').to_string()
        b'150M3S5H'
        >>> CigarSequence(b'150M3I').to_string()
        b'150M3I'
        '''
        cdef bytes    cigar_str
        cdef char*    cigar_ptr
        cdef uint32_t str_len, i

        str_len = self.cigar_len
        for i in range(self.cigar_len):
            str_len += integer_width(self.cigar_buf[i] >> BAM_CIGAR_SHIFT)

        cigar_str = PyBytes_FromStringAndSize(<char*>0, str_len)
        cigar_ptr = cigar_str

        for i in range(self.cigar_len):
            cigar_ptr   += sprintf(cigar_ptr, '%d', self.cigar_buf[i] >> BAM_CIGAR_SHIFT)
            cigar_ptr[0] = get_op_by_binary(self.cigar_buf[i] & BAM_CIGAR_MASK).code
            cigar_ptr   += 1

        #cigar_ptr[0] = 0  # Should not be needed
        return cigar_str

    def to_pysam_list(self):
        '''
        Convert CigarSequence to a pysam CIGAR list of tuplesa

        >>> CigarSequence(b'MID').to_pysam_list()
        [(0, 1), (1, 1), (2, 1)]
        >>> CigarSequence(b'MMMMMIIIIDDDDD').to_pysam_list()
        [(0, 5), (1, 4), (2, 5)]
        >>> CigarSequence(b'10M10I10D').to_pysam_list()
        [(0, 10), (1, 10), (2, 10)]
        '''
        if not self.cigar_len:
            return None

        cdef list cigar = []
        for i in range(self.cigar_len):
            cigar.append( _decode_cigar_to_pysam(self.cigar_buf[i]) )
        return cigar

    def to_op_list(self):
        '''
        Convert CigarSequence to a list of tuples of CigarOperators and counts

        >>> CigarSequence(b'MID').to_op_list()
        [(<CigarOperator:MATCH>, 1), (<CigarOperator:INSERTION>, 1), (<CigarOperator:DELETION>, 1)]
        >>> CigarSequence(b'MMMMMIIIIDDDDD').to_op_list()
        [(<CigarOperator:MATCH>, 5), (<CigarOperator:INSERTION>, 4), (<CigarOperator:DELETION>, 5)]
        >>> CigarSequence(b'10M10I10D').to_op_list()
        [(<CigarOperator:MATCH>, 10), (<CigarOperator:INSERTION>, 10), (<CigarOperator:DELETION>, 10)]
        '''
        cdef list cigar = []
        for i in range(self.cigar_len):
            cigar.append( _decode_cigar_to_op(self.cigar_buf[i]) )
        return cigar


def test_CigarSequence():
    '''
    Various input and output tests

    >>> CigarSequence( [(0, 1), (0, 4), (2, 1), (3, 1), (4, 1)] ).to_string()
    b'5M1D1N1S'
    >>> CigarSequence( (b'150M3S5H') ).to_string()
    b'150M3S5H'
    >>> CigarSequence( (b'150M3II') ).to_string()
    b'150M4I'
    >>> CigarSequence(b'MID').to_pysam_list()
    [(0, 1), (1, 1), (2, 1)]
    >>> CigarSequence(b'MMMMMIIIIDDDDD').to_pysam_list()
    [(0, 5), (1, 4), (2, 5)]
    >>> CigarSequence(b'10M10I10D').to_pysam_list()
    [(0, 10), (1, 10), (2, 10)]

    Test __contains__

    >>> c = CigarSequence(b'10M10I10D')
    >>> b'M' in c
    True
    >>> b'X' in c
    False
    >>> 0 in c
    True
    >>> HARD_CLIP in c
    False

    Test reverse

    >>> c = CigarSequence()
    >>> c.reverse()
    >>> c.to_string()
    b''
    >>> c = CigarSequence(b'1M')
    >>> c.reverse()
    >>> c.to_string()
    b'1M'

    >>> c = CigarSequence(b'1M1I')
    >>> c.reverse()
    >>> c.to_string()
    b'1I1M'

    >>> c = CigarSequence(b'10M10I10D')
    >>> c.to_string()
    b'10M10I10D'
    >>> c.reverse()
    >>> c.to_string()
    b'10D10I10M'

    >>> c = CigarSequence(b'10M10I')
    >>> c.to_string()
    b'10M10I'
    >>> c.reverse()
    >>> c.to_string()
    b'10I10M'

    Test add, extend

    >>> CigarSequence(b'5I') + CigarSequence(b'5M')
    [(<CigarOperator:INSERTION>, 5), (<CigarOperator:MATCH>, 5)]
    >>> CigarSequence(b'5I') + CigarSequence(b'5I')
    [(<CigarOperator:INSERTION>, 10)]
    >>> CigarSequence(b'5I255M') + b'12D'
    [(<CigarOperator:INSERTION>, 5), (<CigarOperator:MATCH>, 255), (<CigarOperator:DELETION>, 12)]

    >>> c  = CigarSequence(b'HSIMD')
    >>> c += b'5D'
    >>> c.to_string()
    b'1H1S1I1M6D'
    >>> c += [(0, 5), (1, 1)]
    >>> c.to_string()
    b'1H1S1I1M6D5M1I'

    >>> c  = CigarSequence()
    >>> c += CigarSequence(b'M')
    >>> c += CigarSequence(b'MM')
    >>> c += CigarSequence(b'MMM')
    >>> c += CigarSequence(b'4M')
    >>> c.to_string()
    b'10M'
    '''


def cigar_to_extended(CigarSequence cigar, bytes ref, bytes query, int ref_start=0, query_start=0):
    cdef char*     sref      = ref
    cdef char*     squery    = query
    cdef int       ref_pos   = ref_start
    cdef int       ref_len   = len(ref)
    cdef int       query_pos = query_start
    cdef int       query_len = len(query)
    cdef int       i, n
    cdef CigarOperator op
    cdef CigarSequence new_cigar = CigarSequence()

    if not isinstance(cigar, CigarSequence):
        raise TypeError('cigar must be instance of CigarSequence')

    for opobj, n in cigar:
        op = <CigarOperator>opobj

        if op.consumes_ref_bases and ref_pos + n > ref_len:
            raise IndexError('reference sequence does not match CIGAR length')
        if op.consumes_read_bases and query_pos + n > query_len:
            raise IndexError('query sequence does not match CIGAR length')

        if op is MATCH:
            for i in range(n):
                op         = SEQ_MATCH if sref[ref_pos] == squery[query_pos] else SEQ_MISMATCH
                ref_pos   += 1
                query_pos += 1
                new_cigar._append_op(op.binary, 1)
        else:
            new_cigar._append_op(op.binary, n)
            if op.consumes_read_bases:
                query_pos += n
            if op.consumes_ref_bases:
                ref_pos += n

    return new_cigar

## <LICENSE>
## Copyright 2014 uta-align Contributors (https://bitbucket.org/biocommons/uta-align)
## 
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
## 
##     http://www.apache.org/licenses/LICENSE-2.0
## 
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
## </LICENSE>

