# -*- coding: utf-8 -*-

from    libc.stdint   cimport uint32_t
from pysam.libcsamfile   cimport AlignedRead

cdef enum:
    BAM_MATCH        = 0
    BAM_INSERTION    = 1
    BAM_DELETION     = 2
    BAM_SKIPPED      = 3
    BAM_SOFT_CLIP    = 4
    BAM_HARD_CLIP    = 5
    BAM_PADDING      = 6
    BAM_SEQ_MATCH    = 7
    BAM_SEQ_MISMATCH = 8


cdef class CigarOperator(object):
    cdef readonly char code, binary, consumes_read_bases, consumes_ref_bases
    cdef readonly str  name, codestr
    cdef public CigarOperator inverse


cdef class CigarSequence(object):
    cdef Py_ssize_t cigar_buf_len, cigar_len
    cdef uint32_t*  cigar_buf

    cdef _extend_from_binary(self, uint32_t other_cigar_len, uint32_t* other_cigar_buf)
    cdef _extend_from_align(self, AlignedRead align)
    cdef _extend_from_cigar_sequence(self, CigarSequence other)
    cdef _extend_from_cigar_slice(self, CigarSequence other, slice sliceobj)
    cdef _extend_from_str(self, bytes cigar_str)
    cdef _extend_from_sequence(self, cigar_seq)
    cdef _alloc_buf(self, uint32_t n)
    cdef inline _append_op(self, char op, uint32_t n)
    cdef _get_slice(self, sliceobj)
    cpdef int count(self, op_obj)
    cpdef int gapped_len(self, int include_soft_clip=?)
    cpdef int ref_len(self, int query_bases=?)
    cpdef int query_len(self, int ref_bases=?, int include_soft_clip=?)
    cpdef append(self, item)
    cpdef extend(self, items)
    cpdef pop(self)
    cpdef reverse(self)
    cpdef invert(self, left_clip=?, right_clip=?)
    cpdef convert_N_to_S(self)


cdef  CigarOperator get_op_by_char(char char_rep)
cpdef CigarOperator get_op_by_str(bytes str_rep)
cpdef CigarOperator get_op_by_binary(char binary_rep)
