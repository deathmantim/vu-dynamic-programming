#!/usr/bin/python3

"""Template for dynamic programming assignment.

The code in the template is compatible with both Python 2 and Python 3
When you finish this code, it should be at least compatible with Python 3.
"""

# Packages for commandline options:
import argparse
import sys
import pickle as pk


# Built-in exchange matrices.
with open('substitution_matrices/identity.pkl', 'rb') as f:
    identity = pk.load(f)

with open('substitution_matrices/pam250.pkl', 'rb') as f:
    pam250 = pk.load(f)

with open('substitution_matrices/blosum62.pkl', 'rb') as f:
    blosum62 = pk.load(f)


def get_args():
    """Collect the inputs."""
    parser = argparse.ArgumentParser(
        prog='PROG',
        usage='%(prog)s [options]',
        description='Aligning two sequences',
        epilog='The code was co-opted from Anton Feenstra\'s and'
        'modified by Cico Zhang'
    )
    parser.add_argument('-f', '--fasta', dest='fasta', metavar='FILE',
                        required=True, help='input alignment file (fasta)')
    parser.add_argument('-e,', '--exchange_matrix', dest='exchange_matrix',
                        metavar='SUBSTITUTION MATRIX NAME', help='Substitution '
                        'matrix: pam250, blosum62 or identity',
                        default='pam250')
    parser.add_argument('-l', '--local', dest='align_local',
                        action='store_true', help='Local alignment',
                        default=False)
    parser.add_argument('-g', '--global', dest='align_global',
                        action='store_true', help='Global alignment',
                        default=False)
    parser.add_argument('-s', '--semi_global', dest='align_semiglobal',
                        action='store_true', help='Semi-global alignment',
                        default=False)
    parser.add_argument('-p', '--penalty', dest='gap_penalty', type=int,
                        help='Gap penalty', default=2)
    parser.add_argument('-o', '--output', dest='alignment', metavar='FILE',
                        default='output.align', help='The file to store the alignment')
    parser.add_argument('-m', '--score_matrix', dest='score_matrix',
                        metavar='FILE', default='output.align',
                        help='The file to store the score matrix')
    parser.add_argument('-v', dest='print_on_screen', action='store_true',
                        help='Print the output (alignment(s) and score '
                        'matrix) on the screen', default=False)

    args = parser.parse_args()

    if args.fasta is None:
        sys.exit('Error: no input file (fasta)')

    if not (args.align_local or args.align_global or args.align_semiglobal):
        sys.exit('Error: No alignment strategy is given: global, local or '
                 'semi-global')
    if args.align_local + args.align_global + args.align_semiglobal > 1:
        sys.exit('Error: More than one alignment strategy is given.')

    if args.exchange_matrix not in ['pam250', 'blosum62', 'identity']:
        sys.exit('Unknown exchange matrix ' + args.exchange_matrix)

    return args


class Sequence:
    """Stores a sequence object."""

    def __init__(self, Label="", Sequence=""):
        """Initialize a new Sequence object.

        Label -- identifier of sequence (text)
        Sequence -- sequence string in single-letter alphabet
        """
        self.Label = Label
        self.Sequence = Sequence

    # this makes that you can do 'print sequence' and get nice output:
    def __str__(self):
        """Return string representation of a Sequence object."""
        # newline-delimited values of all the attributes
        return ">%s\n%s" % (self.Label, self.Sequence)


def readSequences(lines):
    """Return Sequences object.

    lines -- list of lines or any object that behaves like it

    This routine parses a fasta file and returns a list of Sequence objects
    containing the sequences with label and sequence data set
    """
    seqs = []
    label = None
    seq_lines = []
    for line in lines:
        line = line.strip()      # strip off white space
        if not line:             # skip empty lines
            continue
        if line.startswith(';'):  # ignore comment lines
            continue
        # check for start of next sequence:
        if line.startswith('>'):  # label line
            # first, store the previous sequence if we had one:
            if seq_lines:
                seqs.append(Sequence(label, ''.join(seq_lines)))
                seq_lines = []
            # get the label (name) for the next sequence
            label = line[1:].strip()
        else:
            # collect all lines with sequence information for this sequence:
            seq_lines.append(line)
    # take care of the last sequence in the file
    seqs.append(Sequence(label, ''.join(seq_lines)))
    return seqs


def do_global_alignment(sequences, matrix, penalty):
    """Do pairwise global alignment using DP."""
    #########################
    # INSERT YOUR CODE HERE #
    #########################
    seq1 = '-' + sequences[0].Sequence
    seq2 = '-' + sequences[1].Sequence

    scoring = [
        [i for i in range(0, (len(seq2)) * -2, -2)]
    ]
    for i in range(-2, (len(seq1)) * -2, -2):
        scoring.append([i])

    aa_start = ord('A')
    for i in range(1, len(seq1)):
        aa_x = seq1[i]
        for j in range(1, len(seq2)):
            aa_y = seq2[j]
            xgap = scoring[i][j-1] - penalty
            ygap = scoring[i-1][j] - penalty
            match = scoring[i-1][j-1] + matrix[ord(aa_x) - aa_start][ord(aa_y) - aa_start]
            scoring[i].append(max([xgap, ygap, match]))

    alignment = global_traceback(scoring, seq1, seq2, penalty, matrix)

    scoring = add_sequences_to_scoring(scoring, seq1, seq2)
    return alignment, scoring
    #########################
    #   END YOUR CODE HERE  #
    #########################


def global_traceback(scoring, seq1, seq2, penalty, matrix):
    i = len(seq1) - 1
    j = len(seq2) - 1
    alignment = ['', '', '']
    alignmentscore = scoring[i][j]
    aa_start = ord('A')

    while True:
        if i == 0 and j == 0:
            break
        aa_x = seq1[i]
        aa_y = seq2[j]
        score = scoring[i][j]

        if score == scoring[i - 1][j] - penalty:  # high road
            alignment[0] += seq1[i]
            alignment[1] += ' '
            alignment[2] += '-'
            i -= 1
            continue
        if score == scoring[i][j - 1] - penalty:  # low road
            alignment[0] += '-'
            alignment[1] += ' '
            alignment[2] += seq2[j]
            j -= 1
            continue
        if score == scoring[i - 1][j - 1] + matrix[ord(aa_x) - aa_start][ord(aa_y) - aa_start]:  # match
            alignment[0] += seq1[i]
            alignment[2] += seq2[j]
            if seq1[i] == seq2[j]:
                alignment[1] += '|'
            else:
                alignment[1] += ' '
            i, j = i-1, j-1
            continue

    alignment = [[x[::-1]] for x in alignment]
    alignment.append(['score = %s' % alignmentscore])

    return alignment


def add_sequences_to_scoring(scoring, seq1, seq2):
    if seq1[0] != '-':
        seq1 = '-' + seq1
    if seq2[0] != '-':
        seq2 = '-' + seq2

    if len(seq1) == len(scoring):
        new_scoring = [list(' ' + seq2)]

        for i in range(len(seq1)):
            new_scoring.append([seq1[i]])
            new_scoring[-1].extend(scoring[i])

    else:
        new_scoring = [list(' ' + seq1)]

        for i in range(len(seq2)):
            new_scoring.append([seq2[i]])
            new_scoring[-1].extend(scoring[i])

    return new_scoring


def do_local_alignment(sequences, matrix, penalty):
    """Do pairwise local alignment using DP."""
    #########################
    # INSERT YOUR CODE HERE #
    #########################

    #########################
    #   END YOUR CODE HERE  #
    #########################


def do_semiglobal_alignment(sequences, matrix, penalty):
    """Do pairwise semi-global alignment using DP."""
    #########################
    # INSERT YOUR CODE HERE #
    #########################

    #########################
    #   END YOUR CODE HERE  #
    #########################


def print_matrix_to_file(matrix, fileName):
    """Write a matrix into file.

    matrix: a list of list in Python, storing a score
    matrix.
    fileName: str, a file name (with a path) to store the matrix.
    It is not recommended to tinker with this function.
    """
    with open(fileName, 'w') as f:
        for row in matrix:
            print('\t'.join(map(str, row)), file=f)


def print_alignment_to_file(alig, fileName):
    """Write a matrix into file.

    alig: a list of list in Python, storing an alignment.
    fileName: str, a file name (with a path) to store the alignment.
    It is not recommended to tinker with this function.
    """
    with open(fileName, 'w') as f:
        for row in alig:
            print(''.join(map(str, row)), file=f)


def print_matrix_on_screen(matrix, width=5):
    """Print a matrix on the screen.

    matrix: a list of list in Python, storing an alignment or a score
    matrix.
    width: that of the space one cell occupies.
    This will facilitate your testing.
    """
    for row in matrix:
        print(''.join(['{0:>{w}}'.format(item, w=width) for item in row]))


def main():
    """Main function.

    Please change it accordingly to make the program work.
    """
    # Get command line options
    args = get_args()

    # Set substitution matrix:
    if args.exchange_matrix == "pam250":
        exchangeMatrix = pam250
    elif args.exchange_matrix == "blosum62":
        exchangeMatrix = blosum62
    else:
        exchangeMatrix = identity

    # Read sequences from fasta file, and catch error reading file
    try:
        sequences = readSequences(open(args.fasta))
    except OSError as e:
        print("ERROR: cannot open or read fasta input file:", e.filename)

    for seq in sequences:
        print(seq)

    # Call alignment routine(s):
    if args.align_global:
        alignment, score_matrix = do_global_alignment(
                sequences, exchangeMatrix, args.gap_penalty)
    elif args.align_local:
        alignment, score_matrix = do_local_alignment(
                sequences, exchangeMatrix, args.gap_penalty)
    elif args.align_semiglobal:
        alignment, score_matrix = do_semiglobal_alignment(
                sequences, exchangeMatrix, args.gap_penalty)
    else:
        sys.exit("BUG! this should not happen.")
    

    # Print the result to files
    if args.alignment: 
        print_alignment_to_file(alignment, args.alignment)
    if args.score_matrix:
        print_matrix_to_file(score_matrix, args.score_matrix)
    
    # Print the result on screen
    if args.print_on_screen:
        print_matrix_on_screen(alignment)
        print_matrix_on_screen(score_matrix)


if __name__ == "__main__":
    main()

# last line
