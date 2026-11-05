from __future__ import absolute_import


class ORFFinder(object):
    """
    Find the longest ORF and return its length, frame and sequence.
    """

    def __init__(self, seq):
        self.seq = seq.upper()
        self.result = ("+", 0, 0, 0, 0)
        self.winner = 0

    def _reverse_comp(self):
        swap = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N", "X": "X"}
        return "".join(swap[base] for base in self.seq)[::-1]

    def codons(self, frame):
        start = frame
        while start + 3 <= len(self.seq):
            yield (self.seq[start:start + 3], start)
            start += 3

    def run_one(self, frame_number, direction, start_coden, stop_coden):
        codon_gen = self.codons(frame_number)
        while True:
            try:
                codon, index = next(codon_gen)
            except StopIteration:
                break
            if codon in start_coden or (not start_coden and codon not in stop_coden):
                orf_start = index
                ended = False
                while True:
                    try:
                        codon, index = next(codon_gen)
                    except StopIteration:
                        ended = True
                    if codon in stop_coden:
                        ended = True
                    if ended:
                        orf_end = index + 3
                        length = orf_end - orf_start
                        if length > self.winner:
                            self.winner = length
                            self.result = (direction, frame_number + 1, orf_start, orf_end, length)
                        if length == self.winner and orf_start < self.result[2]:
                            self.result = (direction, frame_number + 1, orf_start, orf_end, length)
                        break

    def longest_orf(self, direction, start_coden=None, stop_coden=None):
        if start_coden is None:
            start_coden = ["ATG"]
        if stop_coden is None:
            stop_coden = ["TAG", "TAA", "TGA"]

        if direction == "+":
            for frame in range(3):
                self.run_one(frame, direction, start_coden, stop_coden)
        elif direction == "-":
            self.seq = self._reverse_comp()
            for frame in range(3):
                self.run_one(frame, direction, start_coden, stop_coden)
        else:
            raise ValueError("direction must be '+' or '-'")

        return (self.result[4], self.result[1], self.seq[self.result[2]:self.result[3]])


class FindCDS(object):
    """
    Find the longest CDS-like ORF and keep the original integrity flag.
    """

    def __init__(self, seq):
        self.seq = seq
        self.result = (0, 0, 0, 0, 0)
        self.longest = 0
        self.basepair = {"A": "T", "T": "A", "U": "A", "C": "G", "G": "C", "N": "N", "X": "X"}

    def _reverse_complement(self):
        return "".join(self.basepair[base] for base in self.seq)[::-1]

    def get_codons(self, frame_number):
        coordinate = frame_number
        while coordinate + 3 <= len(self.seq):
            yield (self.seq[coordinate:coordinate + 3], coordinate)
            coordinate += 3

    def find_longest_in_one(self, frame_number, direction, start_codon, stop_codon):
        triplet_gen = self.get_codons(frame_number)
        while True:
            try:
                codon, index = next(triplet_gen)
            except StopIteration:
                break
            if codon in start_codon and codon not in stop_codon:
                orf_start = index
                end_extension = False
                while True:
                    try:
                        codon, index = next(triplet_gen)
                    except StopIteration:
                        end_extension = True
                        integrity = -1
                    if codon in stop_codon:
                        integrity = 1
                        end_extension = True
                    if end_extension:
                        orf_end = index + 3
                        length = orf_end - orf_start
                        if length > self.longest:
                            self.longest = length
                            self.result = [direction, orf_start, orf_end, length, integrity]
                        if length == self.longest and orf_start < self.result[1]:
                            self.result = [direction, orf_start, orf_end, length, integrity]
                        break

    def longest_orf(self, direction, start_codon=None, stop_codon=None):
        if start_codon is None:
            start_codon = {"ATG": None}
        if stop_codon is None:
            stop_codon = {"TAG": None, "TAA": None, "TGA": None}

        for frame in range(3):
            self.find_longest_in_one(frame, "+", start_codon, stop_codon)
        return_orf = self.seq[self.result[1]:self.result[2]][:]
        start_coordinate = self.result[1]
        strand_direction = "+"
        orf_integrity = self.result[4]

        if direction == "-":
            self.seq = self._reverse_complement()
            for frame in range(3):
                self.find_longest_in_one(frame, "-", start_codon, stop_codon)
            if self.result[0] == "-":
                return_orf = self.seq[self.result[1]:self.result[2]][:]
                start_coordinate = self.result[1]
                strand_direction = "-"
                orf_integrity = self.result[4]

        return return_orf, start_coordinate, strand_direction, orf_integrity

