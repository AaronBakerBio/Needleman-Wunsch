"""Aaron Baker
Class: Cs167
Professor: Donna Slonim
HWK 2-2 Needleman-Wunsch
"""
import numpy as np
from node import node
import sys


class GlobalSeqAlignment:
    """Class GlobalSeqAlignment contains methods that calculate the optimal global alignment of
    two sequences using the Needleman-Wunsch algorithm."""

    def __init__(self, x: str, y: str, M=4, m=-2, g=-2):
        """

        :param x: First sequence to align.
        :param y: Second sequence to align.
        :param M: Match value for a matching score.
        :param m: Mismatch value for a mismatch score.
        :param g: The gap cost.
        """
        self.x, self.y, self.M, self.m, self.g = x, y, M, m, g
        self.matrix = np.empty((len(y) + 1, len(x) + 1), dtype=node)
        self.setup_start()
        self.complete_matrix()

    def complete_matrix(self):
        """
        Method complete matrix will fill and compute the initial alignment scoring matrix
        with pointers for traceback.
        """
        rows, cols = self.matrix.shape
        for i in range(1, rows):
            for j in range(1, cols):
                current_score = self.M if self.x[j - 1] == self.y[i - 1] else self.m
                self.matrix[i, j] = node(current_score)
                left = self.matrix[i, j - 1]
                up = self.matrix[i - 1, j]
                diagonal = self.matrix[i - 1, j - 1]
                self.matrix[i, j].create_pointers(left, up, diagonal, self.g)

    def print_matrix(self):
        """@deprecated
        Method print_matrix returns a printout of the matrix for visualization. It is not
        fully implemented nor-properly space as this was not necessary for the program's goal
        to print the scores, solely used for debugging.
        """
        for row in self.matrix:
            for my_node in row:
                print(str(my_node), end=' ')
            print()

    def update_alignment(self, i, j, move):
        """
        Updates alignment based on move direction.
        :param i: Current row index in the alignment matrix, indicating position in sequence y.
        :param j: Current column index in the alignment matrix, indicating position in sequence x.
        :param move: Direction for the traceback ('left', 'up', 'diagonal'), determining the alignment step.
        :return: Updated (i, j) indices and the new aligned sequence fragments.
        """
        aligned_x, aligned_y = "", ""
        if move == "left":
            j -= 1
            aligned_x = self.x[j] + aligned_x
            aligned_y = "-" + aligned_y
        elif move == "up":
            i -= 1
            aligned_x = "-" + aligned_x
            aligned_y = self.y[i] + aligned_y
        elif move == "diagonal":
            i -= 1
            j -= 1
            aligned_x = self.x[j] + aligned_x
            aligned_y = self.y[i] + aligned_y
        return i, j, aligned_x, aligned_y

    def no_print_traceback(self):
        """
        Performs traceback from the bottom-right of the matrix to find the optimal alignment path.
        :return: Tuple containing the aligned sequences x_prime and y_prime.
        """
        i, j = len(self.y), len(self.x)
        aligned_x, aligned_y = "", ""
        while i > 0 or j > 0:
            current_node = self.matrix[i][j]
            if current_node.diagonal and i > 0 and j > 0:
                move = "diagonal"
            elif current_node.left and j > 0:
                move = "left"
            elif current_node.up and i > 0:
                move = "up"
            else:
                break
            i, j, update_x, update_y = self.update_alignment(i, j, move)
            aligned_x = update_x + aligned_x
            aligned_y = update_y + aligned_y
        return aligned_x, aligned_y

    def setup_start(self):
        """
        Initializes the scoring matrix's first row and column with gap penalties,
        setting up for the alignment calculation.
        """
        self.matrix[0, 0] = node()
        x = self.matrix.shape
        for i in range(1, x[1]):
            self.matrix[0, i] = node(self.m)
            self.matrix[0, i].create_pointers(self.matrix[0, i - 1], None, None, self.g)
        for i in range(1, x[0]):
            self.matrix[i, 0] = node(self.m)
            self.matrix[i, 0].create_pointers(self.matrix[i - 1, 0], None, None, self.g)

    def get_optimal_alignment(self):
        """
        Returns the optimal alignment by using the helper method @no_print_traceback to calculate
        the optimal alignment via traceback.
        :return: Tuple (x,y) containing the optimal alignment.
        """
        return self.no_print_traceback()


def parse(filename: str):
    """
    Method parse receives a file name as input and returns a tuple containing the
    sequences inside the file as a tuple.
    :param filename: Filename for fasta file.
    :return: Sequences x and y from the fasta file as a tuple (x,y) of Strings.
    """
    sequences = []
    sequence = ''
    with open(filename, 'r') as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                if sequence:
                    sequences.append(sequence)
                    sequence = ''
                if len(sequences) == 2:
                    break
            else:
                sequence += line
        if sequence and len(sequences) < 2:
            sequences.append(sequence)
    while len(sequences) < 2:
        sequences.append('')
    return sequences[0], sequences[1]


def parse_FASTA_file(fasta_stream):
    """
    Method pase_FASTA_file will return the sequences in the fasta_file passed from stdin.
    :param fasta_stream: input from stdin.
    :return: The sequences from the analyzed
    """
    sequences = []
    sequence = ''
    for line in fasta_stream:
        line = line.strip()
        if line.startswith(">"):
            if sequence:
                sequences.append(sequence)
                sequence = ''
            if len(sequences) == 2:
                break
        else:
            sequence += line
    if sequence and len(sequences) < 2:
        sequences.append(sequence)
    while len(sequences) < 2:
        sequences.append('')
    return sequences[0], sequences[1]


def write_alignment_to_file(x_prime, y_prime, filename="my.output1"):
    """
    Method write_alignment to file is a deprecated test method for generating a file
    containing the optimal alignments. It will write the alignments passed to the file my.output1
    :param x_prime: First sequence that was aligned.
    :param y_prime: Second sequence aligned.
    :param filename: Filename to write to, by default it will be my.output1.
    """
    with open(filename, 'w') as file:
        file.write(x_prime + '\n')
        file.write(y_prime + '\n')


def main():
    """
    Main method is used by bash to return the optimal alignment.
    """
    x, y = parse_FASTA_file(sys.stdin)
    align = GlobalSeqAlignment(x, y)
    x_prime, y_prime = align.get_optimal_alignment()
    print(x_prime)
    print(y_prime)


"""
Test method for main
def main():
    x, y = parse("aligntest.input1.txt")  # Assuming parse is a function that reads and returns the sequences
    align = GlobalSeqAlignment(x, y)  # Assuming GlobalSeqAlignment is initialized properly
    x_prime, y_prime = align.get_optimal_alignment()  # Assuming this method returns the aligned sequences
    print(x_prime)
    print(y_prime)
    write_alignment_to_file(x_prime, y_prime)
"""

if __name__ == "__main__":
    main()
