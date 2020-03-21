# -*- coding: utf-8 -*-

from  __future__ import division

'''Various alignment algorithms'''


### EXTERNAL CODE ###
#
# Original copyright & software license
#
# Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health &
# Human Services. Funded by NCI under Contract N01-CO-12400.
#
# BioInformed LLC and the U.S. Department of Health & Human Services
# (Copyright Holders) hereby grant to the public a perpetual, irrevocable,
# royalty-free non-exclusive, worldwide license to use, copy, distribute,
# display and prepare derivatives of the software in source or binary forms,
# with or without modification, subject only to the following conditions:
#
#   o Redistributions of source code must retain the above copyright notice,
#     this list of conditions and the following disclaimer.
#
#   o Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
#
#   o Neither the name of BioInformed LLC nor the names of its contributors
#     nor the U.S. Department of Health & Human Services may be used to
#     endorse or promote products derived from this software without specific
#     prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

# C imports
cimport cython

from    libc.stdint                     cimport int32_t, uint32_t
from    libc.limits                     cimport INT_MIN
from    libc.stdlib                     cimport malloc, free
from    libc.string                     cimport memset
from    cpython                         cimport PyErr_NoMemory
from    cpython                         cimport PyBytes_FromStringAndSize

from    uta_align.align.cigar_utils cimport CigarSequence, CigarOperator, get_op_by_char

# Python imports
from    uta_align.align.cigar_utils  import INSERTION, DELETION, SOFT_CLIP, SKIPPED


cdef enum:
    # For Gotoh version
    TRACE_DEL      =  1
    TRACE_INS      =  2
    TRACE_MATCH    =  4
    TRACE_DIR_MASK = TRACE_DEL|TRACE_INS|TRACE_MATCH
    TRACE_NEXT_DEL =  8
    TRACE_NEXT_INS = 16

    # For Altshul-Erikson version
    TRACE_A = 0x01
    TRACE_B = 0x02
    TRACE_C = 0x04
    TRACE_D = 0x08
    TRACE_E = 0x10
    TRACE_F = 0x20
    TRACE_G = 0x40

    # Alignment modes
    MODE_GLOBAL       = 1
    MODE_LOCAL        = 2
    MODE_GLOCAL       = 3
    MODE_LOCAL_GLOBAL = 4


cdef inline int32_t min2(int32_t a, int32_t b) nogil:
    return a if a < b else b


cdef inline int32_t max2(int32_t a, int32_t b) nogil:
    return b if a < b else a


cdef inline int32_t max3(int32_t a, int32_t b, int32_t c) nogil:
    return max2(a, max2(b, c))


cdef inline int32_t max4(int32_t a, int32_t b, int32_t c, int32_t d) nogil:
    return max2(max2(a,b), max2(c, d))


def invert_alignment(cigar, start=None, stop=None, length=None):
    '''
    Inverts the alignment between query sequence and reference.

    :param cigar: Query cigar string
    :param start: Alignment start in reference coordinates
    :param stop: Alignment stop in reference coordinates
    :param length: Length of reference sequence
    :return: Inverted cigar string with the right soft clips,
             alignment start in query coordinates,
             alignment stop in query coordinates
    '''
    if start and stop:
        if stop - start != cigar.ref_len():
            raise ValueError('Reference sequence length between start and stop does not match CIGAR')

    left_clip = start
    right_clip = None

    if stop and length:
        right_clip = length - stop

    inv_cigar, s_left, s_right = cigar.invert(left_clip, right_clip)

    return inv_cigar, s_left, s_left + cigar.query_len()


def cigar_alignment(str1, str2, cigar, hide_match=True):
    cdef CigarOperator op
    cdef int hide, k, n
    cdef char *s1
    cdef char *s2
    cdef char *a1
    cdef char *a2
    cdef char op_char

    if len(str1) != cigar.ref_len():
        raise ValueError('Reference length does not match CIGAR')

    if len(str2) != cigar.query_len(include_soft_clip=True):
        raise ValueError('Query length does not match CIGAR')

    n     = cigar.gapped_len(include_soft_clip=True)
    astr1 = PyBytes_FromStringAndSize(<char*>0, n)
    astr2 = PyBytes_FromStringAndSize(<char*>0, n)
    a1    = astr1
    a2    = astr2
    s1    = str1
    s2    = str2
    hide  = bool(hide_match)

    for opobj,n in cigar:
        op = <CigarOperator>opobj

        if op.consumes_read_bases and op.consumes_ref_bases:
            for k in range(n):
                a1[0] = s1[0]
                if hide and s1[0] == s2[0]:
                    a2[0] = b'.'
                else:
                    a2[0] = s2[0]
                a1 += 1
                a2 += 1
                s1 += 1
                s2 += 1

        elif op.consumes_read_bases:
            if op is SOFT_CLIP:
                op_char = b' '
            else:
                op_char = b'-'

            for k in range(n):
                a1[0] = op_char
                a2[0] = s2[0]
                a1 += 1
                a2 += 1
                s2 += 1

        elif op.consumes_ref_bases:
            if op is SKIPPED:
                op_char = b' '
            else:
                op_char = b'-'

            for k in range(n):
                a1[0] = s1[0]
                a2[0] = op_char
                a1 += 1
                a2 += 1
                s1 += 1

    return astr1, astr2


class Alignment(object):
    '''Object with which to store pairwise alignment results

    ref:         reference sequence
    ref_start:   start index of reference alignment
    ref_stop:    stop index of reference alignment
    query:       query sequence
    query_start: start index of reference alignment
    query_stopy: stop index of reference alignment
    cigar:       alignment operations (CIGARSequence object)
    score:       alignment score
    '''
    __slots__ = ('ref','ref_start','ref_stop','query','query_start','query_stop','cigar','score')

    def __init__(self, ref, ref_start, ref_stop, query, query_start, query_stop, cigar, score):
        self.ref         = ref
        self.ref_start   = ref_start
        self.ref_stop    = ref_stop
        self.query       = query
        self.query_start = query_start
        self.query_stop  = query_stop
        self.score       = score

        if cigar is None or isinstance(cigar, CigarSequence):
            self.cigar = cigar
        else:
            self.cigar = CigarSequence(cigar)

    def gapped_alignment(self, hide_match=True):
        if self.cigar is None:
            raise ValueError('Gapped alignment requires full alignment to be performed')

        ref   = self.ref[self.ref_start:self.ref_stop]
        query = self.query[self.query_start:self.query_stop]

        return cigar_alignment(ref, query, self.cigar, hide_match)


@cython.boundscheck(False)
cdef align_local_full(bytes   s1,
                      bytes   s2,
                      int32_t match_score=10,
                      int32_t mismatch_score=-9,
                      int32_t gap_open_score=-15,
                      int32_t gap_extend_score=-6,
                      char    extended_cigar=0,
                      char    soft_clip=0):
    '''
    see align() docs
    '''
    cdef size_t   n      = len(s1)
    cdef size_t   m      = len(s2)
    cdef char*    ss1    = s1
    cdef char*    ss2    = s2

    cdef char*    edits  = <char*>malloc( (n+1)*(m+1)*sizeof(char))
    cdef int32_t* S      = <int32_t*>malloc( (m+1)*sizeof(int32_t))
    cdef int32_t* coldel = <int32_t*>malloc( (m+1)*sizeof(int32_t))
    cdef int32_t  NINF   = INT_MIN+100000
    cdef int32_t  Sup, Sleft, Sdiag, Sij, Smax
    cdef int32_t  inscost, insccost, delcost, delccost, matcost

    cdef size_t   i,j,start_i,start_j,stop_i,stop_j
    cdef char     c1,c2,back
    cdef char*    editpos = edits+m+1

    # legend of space optimizations for i>0 and j>0:
    #   Sij       = S[i,   j  ]
    #   Sdiag     = S[i-1, j-1]
    #   Sup       = S[i-1, j  ]
    #   Sleft     = S[i,   j-1]
    #   coldel[j] = D[i-1, j  ]
    #   inscost   = I[i,   j-1]

    with nogil:
        Sij = Smax = stop_i = stop_j = 0

        for j in range(m+1):
            S[j]       = 0
            coldel[j]  = NINF
            edits[j]   = 0

        for i in range(1,n+1):
            c1         = ss1[i-1]
            Sdiag      = 0
            Sleft      = 0
            inscost    = NINF
            editpos[0] = 0
            editpos   += 1

            for j in range(1,m+1):
                c2        = ss2[j-1]
                Sup       = S[j]

                matcost   = Sdiag + (match_score if c1==c2 else mismatch_score)
                delccost  = coldel[j]  + gap_extend_score
                delcost   = max2(Sup   + gap_open_score, delccost)
                insccost  = inscost    + gap_extend_score
                inscost   = max2(Sleft + gap_open_score, insccost)
                Sij       = max4(0, delcost, inscost, matcost)

                Sdiag     = Sup
                S[j]      = Sij
                Sleft     = Sij
                coldel[j] = delcost

                if Sij > Smax:
                    Smax   = Sij
                    stop_i = i
                    stop_j = j

                back = 0
                if Sij == matcost:
                    back |= TRACE_MATCH
                if Sij == inscost:
                    back |= TRACE_INS
                if Sij == delcost:
                    back |= TRACE_DEL
                if delcost == delccost:
                    back |= TRACE_NEXT_DEL
                if inscost == insccost:
                    back |= TRACE_NEXT_INS

                editpos[0] = back
                editpos   += 1

    start_i,stop_i,start_j,stop_j,cigar = _roll_cigar_gotoh(s1, s2, stop_i, stop_j, edits, MODE_LOCAL, soft_clip, extended_cigar)

    free(coldel)
    free(S)
    free(edits)

    return Alignment(s1, <int32_t>start_i, <int32_t>stop_i,
                     s2, <int32_t>start_j, <int32_t>stop_j,
                     cigar, <int32_t>Smax)


@cython.boundscheck(False)
cdef align_local_score(bytes   s1,
                       bytes   s2,
                       int32_t match_score=10,
                       int32_t mismatch_score=-9,
                       int32_t gap_open_score=-15,
                       int32_t gap_extend_score=-6):
    '''
    see align() docs
    '''
    cdef size_t   n      = len(s1)
    cdef size_t   m      = len(s2)
    cdef char*    ss1    = s1
    cdef char*    ss2    = s2

    cdef int32_t* S      = <int32_t*>malloc(m*sizeof(int32_t))
    cdef int32_t* coldel = <int32_t*>malloc(m*sizeof(int32_t))
    cdef int32_t  NINF   = INT_MIN+100000
    cdef int32_t  Sup, Sleft, Sdiag, Sij, Smax, inscost, delcost

    cdef size_t   i,j,stop_i,stop_j
    cdef char     c1,c2

    with nogil:
        Sij = Smax = 0
        stop_i = stop_j = 0

        for j in range(m):
            coldel[j] = NINF
            S[j]      = 0

        for i in range(n):
            c1      = ss1[i]
            Sleft   = 0
            Sdiag   = 0
            inscost = NINF

            for j in range(m):
                c2        = ss2[j]
                Sup       = S[j]

                matcost   = Sdiag + (match_score if c1==c2 else mismatch_score)
                delcost   = max2(Sup   + gap_open_score, coldel[j] + gap_extend_score)
                inscost   = max2(Sleft + gap_open_score, inscost   + gap_extend_score)
                Sij       = max4(0, delcost, inscost, matcost)

                S[j]      = Sij
                Sleft     = Sij
                Sdiag     = Sup
                coldel[j] = delcost

                if Sij > Smax:
                    Smax   = Sij
                    stop_i = i+1
                    stop_j = j+1

        free(coldel)
        free(S)

    return Alignment(s1, None, <int32_t>stop_i,
                     s2, None, <int32_t>stop_j,
                     None, <int32_t>Smax)


@cython.boundscheck(False)
cdef align_glocal_full(bytes   s1,
                       bytes   s2,
                       int     mode,
                       int32_t match_score=10,
                       int32_t mismatch_score=-9,
                       int32_t gap_open_score=-15,
                       int32_t gap_extend_score=-6,
                       char    extended_cigar=0,
                       char    soft_clip=0):
    '''
    see align() docs
    '''
    cdef size_t   n      = len(s1)
    cdef size_t   m      = len(s2)
    cdef char*    ss1    = s1
    cdef char*    ss2    = s2

    cdef char*    edits  = <char*>malloc( (n+1)*(m+1)*sizeof(char))
    cdef int32_t* S      = <int32_t*>malloc( (m+1)*sizeof(int32_t))
    cdef int32_t* coldel = <int32_t*>malloc( (m+1)*sizeof(int32_t))
    cdef int32_t  NINF   = INT_MIN+100000
    cdef int32_t  Sup, Sleft, Sdiag, Sij, Smax, Smax_col, Smax_row
    cdef int32_t  inscost, insccost, delcost, delccost, matcost

    cdef size_t   i,j,start_i,start_j,stop_i,stop_j
    cdef size_t   stop_i_row,stop_i_col,stop_j_row,stop_j_col
    cdef char     c1,c2,back
    cdef char*    editpos = edits+m+1

    # legend of space optimizations for i>0 and j>0:
    #   Sij       = S[i,   j  ]
    #   Sdiag     = S[i-1, j-1]
    #   Sup       = S[i-1, j  ]
    #   Sleft     = S[i,   j-1]
    #   coldel[j] = D[i-1, j  ]
    #   inscost   = I[i,   j-1]

    with nogil:
        Sij = 0
        Smax = Smax_row = Smax_col = NINF
        i = j = stop_i_row = stop_i_col = stop_j_row = stop_j_col = 0

        for j in range(m+1):
            S[j]       = 0
            coldel[j]  = NINF
            edits[j]   = 0

        for i in range(1,n+1):
            c1         = ss1[i-1]
            Sdiag      = 0
            Sleft      = 0
            inscost    = NINF
            editpos[0] = 0
            editpos   += 1

            for j in range(1,m+1):
                c2       = ss2[j-1]
                Sup      = S[j]

                matcost   = Sdiag + (match_score if c1==c2 else mismatch_score)
                delccost  = coldel[j]  + gap_extend_score
                delcost   = max2(Sup   + gap_open_score, delccost)
                insccost  = inscost    + gap_extend_score
                inscost   = max2(Sleft + gap_open_score, insccost)
                Sij       = max3(delcost, inscost, matcost)

                Sdiag     = Sup
                S[j]      = Sij
                Sleft     = Sij
                coldel[j] = delcost

                if i == n and Sij > Smax_row:
                    Smax_row   = Sij
                    stop_i_row = i
                    stop_j_row = j

                back = 0
                if Sij == matcost:
                    back |= TRACE_MATCH
                if Sij == inscost:
                    back |= TRACE_INS
                if Sij == delcost:
                    back |= TRACE_DEL
                if delcost == delccost:
                    back |= TRACE_NEXT_DEL
                if inscost == insccost:
                    back |= TRACE_NEXT_INS

                editpos[0] = back
                editpos   += 1

            if Sij > Smax_col:
                Smax_col   = Sij
                stop_i_col = i
                stop_j_col = j

        if mode == MODE_LOCAL_GLOBAL:
            Smax, stop_i, stop_j = Sij, n, m
        elif Smax_col >= Smax_row:
            Smax, stop_i, stop_j = Smax_col, stop_i_col, stop_j_col
        else:
            Smax, stop_i, stop_j = Smax_row, stop_i_row, stop_j_row

    start_i,stop_i,start_j,stop_j,cigar = _roll_cigar_gotoh(s1, s2, stop_i, stop_j, edits, MODE_LOCAL, soft_clip, extended_cigar)

    free(coldel)
    free(S)
    free(edits)

    return Alignment(s1, <int32_t>start_i, <int32_t>stop_i,
                     s2, <int32_t>start_j, <int32_t>stop_j,
                     cigar, <int32_t>Smax)


@cython.boundscheck(False)
cdef align_glocal_score(bytes   s1,
                        bytes   s2,
                        int     mode,
                        int32_t match_score=10,
                        int32_t mismatch_score=-9,
                        int32_t gap_open_score=-15,
                        int32_t gap_extend_score=-6):
    '''
    see align() docs
    '''
    cdef size_t   n      = len(s1)
    cdef size_t   m      = len(s2)
    cdef char*    ss1    = s1
    cdef char*    ss2    = s2

    cdef int32_t* S      = <int32_t*>malloc(m*sizeof(int32_t))
    cdef int32_t* coldel = <int32_t*>malloc(m*sizeof(int32_t))
    cdef int32_t  NINF   = INT_MIN+1000000
    cdef int32_t  Sup, Sleft, Sdiag, Sij, Smax, Smax_col, Smax_row
    cdef int32_t  inscost, delcost, matcost

    cdef size_t   i,j,start_i,start_j,stop_i,stop_j
    cdef size_t   stop_i_row,stop_i_col,stop_j_row,stop_j_col
    cdef char     c1,c2,back

    # legend of space optimizations for i>0 and j>0:
    #   Sij       = S[i,   j  ]
    #   Sdiag     = S[i-1, j-1]
    #   Sup       = S[i-1, j  ]
    #   Sleft     = S[i,   j-1]
    #   coldel[j] = D[i-1, j  ]
    #   inscost   = I[i,   j-1]

    with nogil:
        Sij = 0
        Smax = Smax_row = Smax_col = NINF
        i = j = stop_i_row = stop_i_col = stop_j_row = stop_j_col = 0

        for j in range(m):
            S[j]       = 0
            coldel[j]  = NINF

        for i in range(n):
            c1      = ss1[i]
            Sdiag   = 0
            Sleft   = 0
            inscost = NINF

            for j in range(m):
                c2        = ss2[j]
                Sup       = S[j]

                matcost   = Sdiag + (match_score if c1==c2 else mismatch_score)
                delcost   = max2(Sup   + gap_open_score, coldel[j] + gap_extend_score)
                inscost   = max2(Sleft + gap_open_score, inscost   + gap_extend_score)
                Sij       = max3(delcost, inscost, matcost)

                Sdiag     = Sup
                S[j]      = Sij
                Sleft     = Sij
                coldel[j] = delcost

                if i == n-1 and Sij > Smax_row:
                    Smax_row   = Sij
                    stop_i_row = i+1
                    stop_j_row = j+1

            if Sij > Smax_col:
                Smax_col   = Sij
                stop_i_col = i+1
                stop_j_col = j+1

        if mode == MODE_LOCAL_GLOBAL:
            Smax, stop_i, stop_j = Sij, n, m
        elif Smax_col >= Smax_row:
            Smax, stop_i, stop_j = Smax_col, stop_i_col, stop_j_col
        else:
            Smax, stop_i, stop_j = Smax_row, stop_i_row, stop_j_row

        free(coldel)
        free(S)

    return Alignment(s1, None, <int32_t>stop_i,
                     s2, None, <int32_t>stop_j,
                     None, <int32_t>Smax)


@cython.boundscheck(False)
cdef align_global_full(bytes   s1,
                       bytes   s2,
                       int32_t match_score=10,
                       int32_t mismatch_score=-9,
                       int32_t gap_open_score=-15,
                       int32_t gap_extend_score=-6,
                       char    extended_cigar=0):
    '''
    see align() docs
    '''
    cdef size_t   n      = len(s1)
    cdef size_t   m      = len(s2)
    cdef char*    ss1    = s1
    cdef char*    ss2    = s2

    cdef char*    edits  = <char*>malloc( (n+1)*(m+1)*sizeof(char))
    cdef int32_t* S      = <int32_t*>malloc( (m+1)*sizeof(int32_t))
    cdef int32_t* coldel = <int32_t*>malloc( (m+1)*sizeof(int32_t))
    cdef int32_t  NINF   = INT_MIN+100000
    cdef int32_t  Sup, Sleft, Sdiag, Sij, Smax
    cdef int32_t  inscost, insccost, delcost, delccost, matcost

    cdef size_t   i,j,start_i,start_j
    cdef char     c1,c2,back
    cdef char*    editpos = edits+m+1

    # legend of space optimizations for i>0 and j>0:
    #   Sij       = S[i,   j  ]
    #   Sdiag     = S[i-1, j-1]
    #   Sup       = S[i-1, j  ]
    #   Sleft     = S[i,   j-1]
    #   coldel[j] = D[i-1, j  ]
    #   inscost   = I[i,   j-1]

    with nogil:
        Sij = S[0] = 0
        coldel[0]  = NINF
        edits[0]   = 0

        for j in range(1,m+1):
            S[j]       = gap_open_score + (j-1)*gap_extend_score
            coldel[j]  = NINF
            edits[j]   = TRACE_INS

        for i in range(1,n+1):
            c1         = ss1[i-1]
            Sdiag      = S[0]
            S[0]       = gap_open_score + (i-1)*gap_extend_score
            Sleft      = S[0]
            inscost    = NINF
            editpos[0] = TRACE_DEL
            editpos   += 1

            for j in range(1,m+1):
                c2        = ss2[j-1]
                Sup       = S[j]

                matcost   = Sdiag + (match_score if c1==c2 else mismatch_score)
                delccost  = coldel[j]  + gap_extend_score
                delcost   = max2(Sup   + gap_open_score, delccost)
                insccost  = inscost    + gap_extend_score
                inscost   = max2(Sleft + gap_open_score, insccost)
                Sij       = max3(delcost, inscost, matcost)

                Sdiag     = Sup
                S[j]      = Sij
                Sleft     = Sij
                coldel[j] = delcost

                back = 0
                if Sij == matcost:
                    back |= TRACE_MATCH
                if Sij == inscost:
                    back |= TRACE_INS
                if Sij == delcost:
                    back |= TRACE_DEL
                if delcost == delccost:
                    back |= TRACE_NEXT_DEL
                if inscost == insccost:
                    back |= TRACE_NEXT_INS

                editpos[0] = back
                editpos   += 1

        Smax = S[m]

    _,_,_,_,cigar = _roll_cigar_gotoh(s1, s2, n, m, edits, MODE_GLOBAL, 0, extended_cigar)

    free(coldel)
    free(S)
    free(edits)

    return Alignment(s1, 0, <int32_t>n,
                     s2, 0, <int32_t>m,
                     cigar, <int32_t>Smax)


@cython.boundscheck(False)
cdef align_global_score(bytes   s1,
                        bytes   s2,
                        int32_t match_score=10,
                        int32_t mismatch_score=-9,
                        int32_t gap_open_score=-15,
                        int32_t gap_extend_score=-6):
    '''
    see align() docs
    '''
    cdef size_t   n      = len(s1)
    cdef size_t   m      = len(s2)
    cdef char*    ss1    = s1
    cdef char*    ss2    = s2

    cdef int32_t* S      = <int32_t*>malloc(m*sizeof(int32_t))
    cdef int32_t* coldel = <int32_t*>malloc(m*sizeof(int32_t))
    cdef int32_t  NINF   = INT_MIN+1000000
    cdef int32_t  Sfirst, Sup, Sleft, Sdiag, Sij, Smax
    cdef int32_t  inscost, delcost, matcost

    cdef size_t   i,j,start_i,start_j
    cdef char     c1,c2

    # legend of space optimizations for i>0 and j>0:
    #   Sij       = S[i,   j  ]
    #   Sdiag     = S[i-1, j-1]
    #   Sup       = S[i-1, j  ]
    #   Sleft     = S[i,   j-1]
    #   coldel[j] = D[i-1, j  ]
    #   inscost   = I[i,   j-1]

    with nogil:
        Sij = S[0] = 0
        coldel[0]  = NINF

        for j in range(m):
            S[j]       = gap_open_score + j*gap_extend_score
            coldel[j]  = NINF

        Sdiag  = 0
        Sfirst = Sleft = gap_open_score

        for i in range(n):
            c1           = ss1[i]
            inscost      = NINF

            for j in range(m):
                c2        = ss2[j]
                Sup       = S[j]

                matcost   = Sdiag + (match_score if c1==c2 else mismatch_score)
                delcost   = max2(Sup   + gap_open_score, coldel[j] + gap_extend_score)
                inscost   = max2(Sleft + gap_open_score, inscost   + gap_extend_score)
                Sij       = max3(delcost, inscost, matcost)

                Sdiag     = Sup
                S[j]      = Sij
                Sleft     = Sij
                coldel[j] = delcost

            Sdiag   = Sfirst
            Sfirst += gap_extend_score
            Sleft   = Sfirst

        Smax = S[m-1]

    free(coldel)
    free(S)

    return Alignment(s1, None, <int32_t>n,
                     s2, None, <int32_t>m,
                     None, <int32_t>Smax)


cdef _roll_cigar_gotoh(bytes s1, bytes s2, size_t i, size_t j, char* edits, int mode, int soft_clip, int extended_cigar):
    '''
    Compute the sequence of edits required to transform sequence s1 to s2
    using the operations encoded in the supplied matrix of edit operations.
    '''
    cdef size_t   n          = len(s1)+1
    cdef size_t   m          = len(s2)+1
    cdef char*    ss1        = s1
    cdef char*    ss2        = s2
    cdef size_t   ref_stop   = i
    cdef size_t   query_stop = j
    cdef int32_t  count
    cdef char     op, last_op, back, back_dir
    cdef CigarSequence cigar = CigarSequence()

    if soft_clip and j+1 < m:
        cigar.append( (SOFT_CLIP, m - j - 1) )
        query_stop = m

    if mode == MODE_GLOBAL and i+1 < n:
        cigar.append( (DELETION, n - i - 1) )
        ref_stop = n

    if mode == MODE_GLOBAL and j+1 < m:
        cigar.append( (INSERTION, m - j - 1) )
        query_stop = m

    op        = 0
    count     = 0
    back      = TRACE_MATCH
    back_dir  = back

    while i>=0 or j>=0: # back_dir:
        last_back = back
        last_op   = op
        back      = edits[m*i+j]
        back_dir  = back & TRACE_DIR_MASK

        if op==b'D' and (last_back & TRACE_NEXT_DEL):
            op = b'D'
        elif op==b'I' and (last_back & TRACE_NEXT_INS):
            op = b'I'
        elif back_dir & TRACE_MATCH:
            op = b'M'
        elif back_dir & TRACE_DEL:
            op = b'D'
        elif back_dir & TRACE_INS:
            op = b'I'
        else:
            break

        if op==b'M':
            i -= 1
            j -= 1
        elif op==b'D':
            i -= 1
        elif op==b'I':
            j -= 1
        else:
            raise ValueError('Invalid edit operation')

        if extended_cigar and op==b'M':
            if ss1[i]==ss2[j]:
                op = b'='
            else:
                op = b'X'

        if count and last_op != op:
            cigar.append( (get_op_by_char(last_op),count) )
            count = 1
        else:
            count += 1

    if count:
        cigar.append( (get_op_by_char(op),count) )

    if soft_clip and j > 0:
        cigar.append( (SOFT_CLIP, j) )
        j = 0

    if mode == MODE_GLOBAL and i > 0:
        cigar.append( (DELETION, j) )
        i = 0

    if mode == MODE_GLOBAL and j > 0:
        cigar.append( (INSERTION, j) )
        j = 0

    cigar.reverse()

    return i, ref_stop, j, query_stop, cigar

###################################################################################################

def align(bytes   ref,
          bytes   query,
          bytes   mode,
          char    score_only=0,
          int32_t match_score=10,
          int32_t mismatch_score=-9,
          int32_t gap_open_score=-15,
          int32_t gap_extend_score=-6,
          char    extended_cigar=0,
          char    soft_clip=0):
    '''
    Align ref to query using a global, glocal or hybrid global/local (glocal) gapped alignment with
    affine gap penalties using the scoring parameters provided.  An Alignment object is returned.

    :param              ref: reference sequence  (str/bytes)
    :param            query: query sequence (str/bytes)
    :param             mode: 'global', 'local' or 'glocal'
    :param       score_only: Compute only optimal alignment score and stop coordinates (boolean, default=False)
    :param      match_score: match score (int, default=10)
    :param   mismatch_score: mismatch score (int, default=-9)
    :param   gap_open_score: gap open score (int, default=-15)
    :param gap_extend_score: gap extend score (int, default=-6)
    :param   extended_cigar: produce extended CIGAR (boolean, default=False)
    :param        soft_clip: add soft-clipping to local/glocal alignments (boolean, default=False)
    :return:                 Alignment object

    If mode == 'global' then the Needleman-Wunsch-Gotoh algorithm is used to return the global
    alignment of ref to query.

    If mode == 'local' then the Smith-Waterman-Gotoh algorithm is used to return the local
    alignment of ref to query.

    If mode == 'glocal' then a variant of the Needleman-Wunsch-Gotoh algorithm is used to return the
    alignment of ref to query.  This alignment is similar to global alignment with no penalty for leading
    or trailing gaps, provided that the maximum scoring alignments that spans from the start of the
    reference or the query sequence and read the end of the reference or the query sequence.

    If mode == 'local_global' then a variant of the Needleman-Wunsch-Gotoh algorithm is used to
    return the alignment of ref to query.  This alignment is similar to global alignment with no
    penalty for leading gaps and to glocal alignments except the alignment must span to the end of
    both the reference and query alignments.

    If score_only evaluates to True, then only the optimal alignment score and end coordinates of the
    alignment are returned.  Otherwise, the start coordinates and alignment operations are also returned.
    The operations to align ref to query are returned as a sequence, represented by tuples of extended
    CIGAR (Compact Idiosyncratic Gapped Alignment Report) operations and counts:

      Match or mismatch: (MATCH,        count) if not extended_cigar
      Match:             (SEQ_MATCH,    count) if extended_cigar
      Mismatch:          (SEQ_MISMATCH, count) if extended_cigar
      Insertion:         (INSERTION,    count)
      Deletion:          (DELETION,     count)

    This implementation is based on a standard dynamic programming algorithm, requiring O(n*m) time,
    where n and m are the lengths of the two sequences.  If full alignments are requested (not score_only),
    then O(n*m) space is also required.

    See also:

         http://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm

         http://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm

         Needleman, Saul B.; and Wunsch, Christian D. (1970). "A general method applicable to the
         search for similarities in the amino acid sequence of two proteins".  Journal of Molecular
         Biology 48 (3): 443.53.  doi:10.1016/0022-2836(70)90057-4.  PMID 5420325.

         Smith, Temple F.; and Waterman, Michael S. (1981). "Identification of Common Molecular
         Subsequences".  Journal of Molecular Biology 147: 195–197.
         http://gel.ym.edu.tw/~chc/AB_papers/03.pdf

         Gotoh O (1982). "An improved algorithm for matching biological sequences." J Mol Biol,
         162:705-708.
    '''
    if match_score <= mismatch_score:
        raise ValueError('match_score must be greater than mismatch_score')
    if match_score <= gap_open_score:
        raise ValueError('match_score must be greater than gap_open_score')
    if match_score <= gap_extend_score:
        raise ValueError('match_score must be greater than gap_extend_score')
    if gap_open_score > gap_extend_score:
        raise ValueError('gap_open_score must be less than or equal to gap extend_score')

    mode = mode.lower()

    if mode == b'local':
        if score_only:
            return align_local_score(ref, query, match_score, mismatch_score, gap_open_score, gap_extend_score)
        else:
            return align_local_full(ref, query, match_score, mismatch_score, gap_open_score, gap_extend_score, extended_cigar, soft_clip)
    elif mode == b'global':
        if score_only:
            return align_global_score(ref, query, match_score, mismatch_score, gap_open_score, gap_extend_score)
        else:
            return align_global_full(ref, query, match_score, mismatch_score, gap_open_score, gap_extend_score, extended_cigar)
    elif mode == b'glocal':
        if score_only:
            return align_glocal_score(ref, query, MODE_GLOCAL, match_score, mismatch_score, gap_open_score, gap_extend_score)
        else:
            return align_glocal_full(ref, query, MODE_GLOCAL, match_score, mismatch_score, gap_open_score, gap_extend_score, extended_cigar, soft_clip)
    elif mode == b'local_global':
        if score_only:
            return align_glocal_score(ref, query, MODE_LOCAL_GLOBAL, match_score, mismatch_score, gap_open_score, gap_extend_score)
        else:
            return align_glocal_full(ref, query, MODE_LOCAL_GLOBAL, match_score, mismatch_score, gap_open_score, gap_extend_score, extended_cigar, soft_clip)
    else:
        raise ValueError('Unknown alignment mode ' + mode)


###################################################################################################


def test_align():
    '''
    >>> s1,s2=b'b',b'abc'
    >>> a = align(s1,s2,b'local')
    >>> a.score
    10
    >>> a.cigar.to_string()
    b'1M'
    >>> a1,a2 = a.gapped_alignment()
    >>> a1
    b'b'
    >>> a2
    b'.'

    >>> s1,s2=b'abc',b'b'
    >>> a = align(s1,s2,b'local')
    >>> a.ref_start, a.ref_stop
    (1, 2)
    >>> a.query_start, a.query_stop
    (0, 1)
    >>> a.score
    10
    >>> a.cigar.to_string()
    b'1M'
    >>> a1,a2 = a.gapped_alignment()
    >>> a1
    b'b'
    >>> a2
    b'.'

    >>> s1,s2=b'abbcbbd',b'acd'
    >>> a = align(s1,s2,b'local',match_score=30)
    >>> a.ref_start, a.ref_stop
    (0, 7)
    >>> a.query_start, a.query_stop
    (0, 3)
    >>> a.score
    48
    >>> a.cigar.to_string()
    b'1M2D1M2D1M'
    >>> a1,a2 = a.gapped_alignment()
    >>> a1
    b'abbcbbd'
    >>> a2
    b'.--.--.'

    >>> s1=b'AGACCAAGTCTCTGCTACCGTACATACTCGTACTGAGACTGCCAAGGCACACAGGGGATAG'
    >>> s2=b'GCTGGTGCGACACAT'
    >>> a = align(s1,s2,b'local',mismatch_score=-20)
    >>> a.ref_start, a.ref_stop
    (46, 53)
    >>> a.query_start, a.query_stop
    (6, 14)
    >>> a.score
    55
    >>> a.cigar.to_string()
    b'2M1I5M'
    >>> a1,a2 = a.gapped_alignment()
    >>> a1
    b'GC-ACACA'
    >>> a2
    b'..G.....'

    >>> s1,s2=b'b',b'abc'
    >>> a = align(s1,s2,b'glocal')
    >>> a.score
    10
    >>> a.cigar.to_string()
    b'1M'
    >>> a1,a2 = a.gapped_alignment()
    >>> a1
    b'b'
    >>> a2
    b'.'

    >>> s1,s2=b'abc',b'b'
    >>> a = align(s1,s2,b'glocal')
    >>> a.ref_start, a.ref_stop
    (1, 2)
    >>> a.query_start, a.query_stop
    (0, 1)
    >>> a.score
    10
    >>> a.cigar.to_string()
    b'1M'
    >>> a1,a2 = a.gapped_alignment()
    >>> a1
    b'b'
    >>> a2
    b'.'

    >>> s1,s2=b'abbcbbd',b'acd'
    >>> a = align(s1,s2,b'glocal',match_score=30)
    >>> a.ref_start, a.ref_stop
    (0, 7)
    >>> a.query_start, a.query_stop
    (0, 3)
    >>> a.score
    48
    >>> a.cigar.to_string()
    b'1M2D1M2D1M'
    >>> a1,a2 = a.gapped_alignment()
    >>> a1
    b'abbcbbd'
    >>> a2
    b'.--.--.'

    >>> s1=b'AGACCAAGTCTCTGCTACCGTACATACTCGTACTGAGACTGCCAAGGCACACAGGGGATAG'
    >>> s2=b'GCTGGTGCGACACAT'
    >>> a = align(s1,s2,b'glocal',mismatch_score=-20)
    >>> a.ref_start, a.ref_stop
    (36, 53)
    >>> a.query_start, a.query_stop
    (0, 15)
    >>> a.score
    27
    >>> a.cigar.to_string()
    b'1M1D3M4D1M1I2M1I5M1I'
    >>> a1,a2 = a.gapped_alignment()
    >>> a1
    b'GACTGCCAAG-GC-ACACA-'
    >>> a2
    b'.-...----.T..G.....T'

    >>> s1,s2=b'b',b'abc'
    >>> a = align(s1,s2,b'glocal',score_only=True)
    >>> a.ref_stop, a.query_stop, a.score
    (1, 2, 10)

    >>> s1,s2=b'abc',b'b'
    >>> a = align(s1,s2,b'glocal',score_only=True)
    >>> a.ref_stop, a.query_stop, a.score
    (2, 1, 10)

    >>> s1,s2=b'abbcbbd',b'acd'
    >>> a = align(s1,s2,b'glocal',score_only=True,match_score=30)
    >>> a.ref_stop, a.query_stop, a.score
    (7, 3, 48)

    >>> s1=b'AGACCAAGTCTCTGCTACCGTACATACTCGTACTGAGACTGCCAAGGCACACAGGGGATAG'
    >>> s2=b'GCTGGTGCGACACAT'
    >>> a = align(s1,s2,b'glocal',score_only=True,mismatch_score=-20)
    >>> a.ref_stop, a.query_stop, a.score
    (53, 15, 27)

    >>> s1,s2=b'b',b'abc'
    >>> a = align(s1,s2,b'local',score_only=True)
    >>> a.ref_stop, a.query_stop, a.score
    (1, 2, 10)

    >>> s1,s2=b'abc',b'b'
    >>> a = align(s1,s2,b'local',score_only=True)
    >>> a.ref_stop, a.query_stop, a.score
    (2, 1, 10)

    >>> s1,s2=b'abbcbbd',b'acd'
    >>> a = align(s1,s2,b'local',score_only=True,match_score=30)
    >>> a.ref_stop, a.query_stop, a.score
    (7, 3, 48)

    >>> s1=b'AGACCAAGTCTCTGCTACCGTACATACTCGTACTGAGACTGCCAAGGCACACAGGGGATAG'
    >>> s2=b'GCTGGTGCGACACAT'
    >>> a = align(s1,s2,b'local',score_only=True,mismatch_score=-20)
    >>> a.ref_stop, a.query_stop, a.score
    (53, 14, 55)

    >>> s1,s2=b'b',b'abc'
    >>> a = align(s1,s2,b'global')
    >>> a.score
    -20
    >>> a.cigar.to_string()
    b'1I1M1I'
    >>> a1,a2 = a.gapped_alignment()
    >>> a1
    b'-b-'
    >>> a2
    b'a.c'

    >>> s1,s2=b'abc',b'b'
    >>> a = align(s1,s2,b'global')
    >>> a.score
    -20
    >>> a.cigar.to_string()
    b'1D1M1D'
    >>> a1,a2 = a.gapped_alignment()
    >>> a1
    b'abc'
    >>> a2
    b'-.-'

    >>> s1,s2=b'abbcbbd',b'acd'
    >>> a = align(s1,s2,b'global',match_score=30)
    >>> a.score
    48
    >>> a.cigar.to_string()
    b'1M2D1M2D1M'
    >>> a1,a2 = a.gapped_alignment()
    >>> a1
    b'abbcbbd'
    >>> a2
    b'.--.--.'

    >>> s1=b'AGACCAAGTCTCTGCTACCGTACATACTCGTACTGAGACTGCCAAGGCACACAGGGGATAG'
    >>> s2=b'GCTGGTGCGACACAT'
    >>> a = align(s1,s2,b'global')
    >>> a.score
    -198
    >>> a.cigar.to_string()
    b'1D1M9D3M5D1M19D3M3D1M2D4M5D2M2D'
    >>> a1,a2 = a.gapped_alignment()
    >>> a1
    b'AGACCAAGTCTCTGCTACCGTACATACTCGTACTGAGACTGCCAAGGCACACAGGGGATAG'
    >>> a2
    b'-.---------...-----.-------------------...---.--....-----..--'

    >>> s1,s2=b'b',b'abc'
    >>> a = align(s1,s2,b'global',score_only=True)
    >>> a.score
    -20

    >>> s1,s2=b'abc',b'b'
    >>> a = align(s1,s2,b'global',score_only=True)
    >>> a.score
    -20

    >>> s1,s2=b'abbcbbd',b'acd'
    >>> a = align(s1,s2,b'global',score_only=True,match_score=30)
    >>> a.score
    48

    >>> s1=b'AGACCAAGTCTCTGCTACCGTACATACTCGTACTGAGACTGCCAAGGCACACAGGGGATAG'
    >>> s2=b'GCTGGTGCGACACAT'
    >>> a = align(s1,s2,b'global',score_only=True)
    >>> a.score
    -198
    '''

###################################################################################################


@cython.boundscheck(False)
def needleman_wunsch_altshul_erikson(bytes   s1,
                                     bytes   s2,
                                     int32_t match_score=10,
                                     int32_t mismatch_score=-9,
                                     int32_t gap_open_score=-15,
                                     int32_t gap_extend_score=-6,
                                     char    extended_cigar=0):
    '''
    Align s1 to s2 using the Needleman-Wunsch-Altshul-Erikson algorithm for
    global gapped alignment with affine gap penalties.  The optimal
    alignment score and a corresponding sequence of alignment operations are
    returned.

    See: http://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm

         Needleman, Saul B.; and Wunsch, Christian D. (1970). "A general
         method applicable to the search for similarities in the amino acid
         sequence of two proteins".  Journal of Molecular Biology 48 (3):
         443.53.  doi:10.1016/0022-2836(70)90057-4.  PMID 5420325.

         Altschul SF, Erickson BW.  "Optimal sequence alignment using affine
         gap costs." Bull Math Biol.  1986;48(5-6):603-16.

    >>> s1,s2=b'b',b'abc'
    >>> a = needleman_wunsch_altshul_erikson(s1,s2)
    >>> a.score
    -20
    >>> a.cigar.to_string()
    b'1I1M1I'
    >>> a1,a2 = a.gapped_alignment()
    >>> a1
    b'-b-'
    >>> a2
    b'a.c'

    >>> s1,s2=b'abc',b'b'
    >>> a = needleman_wunsch_altshul_erikson(s1,s2)
    >>> a.score
    -20
    >>> a.cigar.to_string()
    b'1D1M1D'
    >>> a1,a2 = a.gapped_alignment()
    >>> a1
    b'abc'
    >>> a2
    b'-.-'

    >>> s1,s2=b'abbcbbd',b'acd'
    >>> a = needleman_wunsch_altshul_erikson(s1,s2,match_score=30)
    >>> a.score
    48
    >>> a.cigar.to_string()
    b'1M2D1M2D1M'
    >>> a1,a2 = a.gapped_alignment()
    >>> a1
    b'abbcbbd'
    >>> a2
    b'.--.--.'

    >>> s1=b'AGACCAAGTCTCTGCTACCGTACATACTCGTACTGAGACTGCCAAGGCACACAGGGGATAG'
    >>> s2=b'GCTGGTGCGACACAT'
    >>> a = needleman_wunsch_altshul_erikson(s1,s2)
    >>> a.score
    -198
    >>> a.cigar.to_string()
    b'1D1M9D3M5D1M19D3M3D1M2D4M5D2M2D'
    >>> a1,a2 = a.gapped_alignment()
    >>> a1
    b'AGACCAAGTCTCTGCTACCGTACATACTCGTACTGAGACTGCCAAGGCACACAGGGGATAG'
    >>> a2
    b'-.---------...-----.-------------------...---.--....-----..--'
    '''
    cdef size_t   n      = len(s1)
    cdef size_t   m      = len(s2)
    cdef char*    ss1    = s1
    cdef char*    ss2    = s2

    # Effective negative infinity
    cdef int32_t  NINF   = INT_MIN+1000000

    # Score matricies
    cdef char*    trace =    <char*>malloc( (n+2)*(m+2)*sizeof(char) )
    cdef int32_t* R     = <int32_t*>malloc( (m+1)*sizeof(int32_t) )
    cdef int32_t* P     = <int32_t*>malloc( (m+1)*sizeof(int32_t) )
    cdef int32_t  Qij, Qleft
    cdef int32_t  Rij, Rup, Rleft, Rdiag, Rmax
    cdef int32_t  Pij, Pup

    # Temp score variables
    cdef int32_t  delopen, delext, insopen, insext, match

    # Various indices
    cdef int32_t  i,j,ij,k,start_i,start_j

    # Current sequence characters
    cdef char     c1,c2

    # Backtrace flags
    cdef char     a,b,c,d,e,f,g

    R[0] = 0
    P[0] = NINF

    memset(<void*>trace, 0, (n+2)*(m+2))
    trace[0]               = TRACE_E|TRACE_G
    trace[(n+2)*(m+2) - 1] = TRACE_C

    for j in range(1,m+1):
        R[j]     = gap_open_score + (j-1)*gap_extend_score
        P[j]     = NINF
        trace[j] = TRACE_B|TRACE_E
        if j<m:
            trace[j] |= TRACE_F

    for i in range(1,n+1):
        c1         = ss1[i-1]
        Rdiag      = R[0]
        R[0]       = gap_open_score + (i-1)*gap_extend_score
        P[0]       = NINF
        Qij        = NINF
        Rleft      = R[0]
        Qleft      = Qij
        trace[i*(m+2)] |= TRACE_A|TRACE_G
        if i<n:
            trace[i*(m+2)] |= TRACE_D

        for j in range(1,m+1):
            c2       = ss2[j-1]
            Rup      = R[j]
            Pup      = P[j]

            delopen  = Rup + gap_open_score
            delext   = Pup + gap_extend_score
            Pij      = max2(delopen, delext)

            insopen  = Rleft + gap_open_score
            insext   = Qleft + gap_extend_score
            Qij      = max2(insopen, insext)

            match    = Rdiag + (match_score if c1==c2 else mismatch_score)

            Rij      = max3(Pij, Qij, match)

            k = (i-1)*(m+2) + j
            if Pij == delext:
                trace[k] |= TRACE_D
            if Pij == delopen:
                trace[k] |= TRACE_E

            k = i*(m+2) + j - 1
            if Qij == insext:
                trace[k] |= TRACE_F
            if Qij == insopen:
                trace[k] |= TRACE_G

            k = i*(m+2)+j
            if Rij == Pij:
                trace[k] |= TRACE_A
            if Rij == Qij:
                trace[k] |= TRACE_B
            if Rij == match:
                trace[k] |= TRACE_C

            Rdiag        = Rup
            R[j] = Rleft = Rij
            P[j]         = Pij
            Qleft        = Qij

    for i in range(n,-1,-1):
        for j in range(m,-1,-1):
            ij = i*(m+2)+j
            a = trace[ij+m+2] & TRACE_A   # A[i+1,j  ]
            b = trace[ij+1]   & TRACE_B   # B[i,  j+1]
            c = trace[ij+m+3] & TRACE_C   # C[i+1,j+1]
            d = trace[ij]     & TRACE_D   # D[i,  j  ]
            e = trace[ij]     & TRACE_E   # E[i,  j  ]
            f = trace[ij]     & TRACE_F   # F[i,  j  ]
            g = trace[ij]     & TRACE_G   # G[i,  j  ]

            if a or b or c:
                if a and d:
                    trace[ij+m+2]  = (trace[ij+m+2]&(~TRACE_D))|((e^TRACE_E)>>1)
                    trace[ij]      = (trace[ij]&(~TRACE_E))|(((trace[ij]&TRACE_A)^TRACE_A)<<4)|TRACE_A
                else:
                    trace[ij+m+2] &= ~TRACE_D
                    trace[ij]     &= ~TRACE_E

                if b and f:
                    trace[ij+1]    = (trace[ij+1]&(~TRACE_F))|((g^TRACE_G)>>1)
                    trace[ij]      = (trace[ij]&(~TRACE_G))|(((trace[ij]&TRACE_B)^TRACE_B)<<5)|TRACE_B
                else:
                    trace[ij+1]   &= ~TRACE_F
                    trace[ij]     &= ~TRACE_G

    Rmax = R[m]
    start_i,start_j,cigar = _roll_cigar_altshul_erikson(s1,s2,n,m,trace,0,extended_cigar)

    free(P)
    free(R)
    free(trace)

    return Alignment(s1, 0, <int32_t>n,
                     s2, 0, <int32_t>m,
                     cigar, <int32_t>Rmax)


cdef _roll_cigar_altshul_erikson(bytes s1, bytes s2, size_t i, size_t j, char* trace, int soft_clip, int extended_cigar):
    '''
    Compute the sequence of edits required to transform sequence s1 to s2
    using the operations encoded in the supplied matrix of edit operations.
    '''
    cdef size_t   n      = len(s1)
    cdef size_t   m      = len(s2)
    cdef char*    ss1    = s1
    cdef char*    ss2    = s2
    cdef int32_t  count
    cdef char     op, last_op, back
    cdef CigarSequence cigar = CigarSequence()

    if soft_clip and j < m:
        cigar.append( (SOFT_CLIP, m - j) )

    op    = 0
    count = 0

    while i>=0 or j>=0: # back_dir:
        last_op = op
        back    = trace[(m+2)*i+j]

        if op==b'D' and back&(TRACE_A|TRACE_D):
            op = b'D'
        elif op==b'I' and back&(TRACE_B|TRACE_F):
            op = b'I'
        elif back&TRACE_C:
            op = b'M'
        elif back&(TRACE_A|TRACE_E):
            op = b'D'
        elif back&(TRACE_B|TRACE_G):
            op = b'I'
        else:
            break
            #raise ValueError('Invalid traceback at %d,%d = %d' % (i,j,back))

        if op==b'M':
            i -= 1
            j -= 1
        elif op==b'D':
            i -= 1
        elif op==b'I':
            j -= 1
        else:
            raise ValueError('Invalid edit operation')

        if extended_cigar and op==b'M':
            if ss1[i]==ss2[j]:
                op = b'='
            else:
                op = b'X'

        if count and last_op != op:
            cigar.append( (get_op_by_char(last_op),count) )
            count = 1
        else:
            count += 1

    if count:
        cigar.append( (get_op_by_char(op),count) )

    if soft_clip and j > 0:
        cigar.append( (SOFT_CLIP, j) )

    cigar.reverse()

    return i, j, cigar


# Backward compatible function wrappers

def smith_waterman_gotoh_align(bytes   s1,
                               bytes   s2,
                               int32_t match_score=10,
                               int32_t mismatch_score=-9,
                               int32_t gap_open_score=-15,
                               int32_t gap_extend_score=-6,
                               char    extended_cigar=0,
                               char    add_soft_clips=0):

    a = align_local_full(s1, s2, match_score, mismatch_score, gap_open_score, gap_extend_score, extended_cigar, add_soft_clips)

    return (slice(a.ref_start,   a.ref_stop),
            slice(a.query_start, a.query_stop),
            a.score,
            a.cigar)


def smith_waterman_gotoh_score(bytes   s1,
                               bytes   s2,
                               int32_t match_score=10,
                               int32_t mismatch_score=-9,
                               int32_t gap_open_score=-15,
                               int32_t gap_extend_score=-6):

    a = align_local_score(s1, s2, match_score, mismatch_score, gap_open_score, gap_extend_score)

    return a.ref_stop, a.query_stop, a.score


def needleman_wunsch_gotoh_align(bytes s1,
                                 bytes s2,
                                 int32_t match_score=10,
                                 int32_t mismatch_score=-9,
                                 int32_t gap_open_score=-15,
                                 int32_t gap_extend_score=-6,
                                 char    extended_cigar=0):

    a = align_global_full(s1, s2, match_score, mismatch_score, gap_open_score, gap_extend_score, extended_cigar)

    return a.score,a.cigar


def needleman_wunsch_gotoh_score(bytes   s1,
                                 bytes   s2,
                                 int32_t match_score=10,
                                 int32_t mismatch_score=-9,
                                 int32_t gap_open_score=-15,
                                 int32_t gap_extend_score=-6):

    a = align_global_score(s1, s2, match_score, mismatch_score, gap_open_score, gap_extend_score)

    return a.score


def glocal_align(bytes   s1,
                 bytes   s2,
                 int32_t match_score=10,
                 int32_t mismatch_score=-9,
                 int32_t gap_open_score=-15,
                 int32_t gap_extend_score=-6,
                 char    extended_cigar=0,
                 char    add_soft_clips=0):

    a = align_glocal_full(s1, s2, MODE_GLOCAL, match_score, mismatch_score, gap_open_score, gap_extend_score, extended_cigar, add_soft_clips)

    return (slice(a.ref_start,   a.ref_stop),
            slice(a.query_start, a.query_stop),
            a.score,
            a.cigar)

def glocal_score(bytes   s1,
                 bytes   s2,
                 int32_t match_score=10,
                 int32_t mismatch_score=-9,
                 int32_t gap_open_score=-15,
                 int32_t gap_extend_score=-6):

    a = align_glocal_score(s1, s2, MODE_GLOCAL, match_score, mismatch_score, gap_open_score, gap_extend_score)

    return a.ref_stop, a.query_stop, a.score

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

