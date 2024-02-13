class node:
    def __init__(self, matches=0):
        """
        Initializes node with match score and pointers set to None.
        :param matches: Match score for the node.
        """
        self.match = matches
        self.left = None
        self.up = None
        self.diagonal = None
        self.value = 0

    @staticmethod
    def score_with_gap(neighbor_node, gap):
        """
        Calculates score with gap penalty if node exists.
        :param neighbor_node: The neighboring node to calculate score from.
        :param gap: Gap penalty value.
        :return: Score considering the gap penalty or -inf if node is None.
        """
        return neighbor_node.value + gap if neighbor_node else float('-inf')

    def create_pointers(self, left_node=None, up_node=None, diagonal_node=None, gap=-1):
        """
        Determines pointers based on scores and updates node's value.
        :param left_node: Left node.
        :param up_node: Up node.
        :param diagonal_node: Diagonal node.
        :param gap: Gap penalty value.
        """
        left_score = self.score_with_gap(left_node, gap)
        up_score = self.score_with_gap(up_node, gap)
        diagonal_score = diagonal_node.value + self.match if diagonal_node else float('-inf')
        self.value = max(left_score, up_score, diagonal_score)
        if self.value == left_score and left_node: self.left = left_node
        if self.value == up_score and up_node: self.up = up_node
        if self.value == diagonal_score and diagonal_node: self.diagonal = diagonal_node

    def print_val(self):
        """
        Prints the node's value.
        """
        print(self.value)

    def __str__(self):
        """
        String representation of the node's value.
        :return: String of value
        """
        return str(self.value)
