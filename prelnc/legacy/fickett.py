from __future__ import absolute_import


class Fickett(object):
    """
    Calculate the Fickett TESTCODE score for a full transcript sequence.
    """

    def __init__(self):
        self.position_parameter = [1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 0.0]
        self.content_parameter = [0.33, 0.31, 0.29, 0.27, 0.25, 0.23, 0.21, 0.19, 0.17, 0]
        self.position_probability = {
            "A": [0.51, 0.55, 0.57, 0.52, 0.48, 0.58, 0.57, 0.54, 0.50, 0.36],
            "C": [0.29, 0.44, 0.55, 0.49, 0.52, 0.60, 0.60, 0.56, 0.51, 0.38],
            "G": [0.62, 0.67, 0.74, 0.65, 0.61, 0.62, 0.52, 0.41, 0.31, 0.17],
            "T": [0.51, 0.60, 0.69, 0.64, 0.62, 0.67, 0.58, 0.48, 0.39, 0.24],
        }
        self.position_weight = {"A": 0.062, "C": 0.093, "G": 0.205, "T": 0.154}
        self.content_probability = {
            "A": [0.40, 0.55, 0.58, 0.58, 0.52, 0.48, 0.45, 0.45, 0.38, 0.19],
            "C": [0.50, 0.63, 0.59, 0.50, 0.46, 0.45, 0.47, 0.56, 0.59, 0.33],
            "G": [0.21, 0.40, 0.47, 0.50, 0.52, 0.56, 0.57, 0.52, 0.44, 0.23],
            "T": [0.30, 0.49, 0.56, 0.53, 0.48, 0.48, 0.52, 0.57, 0.60, 0.51],
        }
        self.content_weight = {"A": 0.084, "C": 0.076, "G": 0.081, "T": 0.055}

    def _lookup_position_probability(self, value, base):
        if float(value) < 0:
            return None
        for idx, threshold in enumerate(self.position_parameter):
            if float(value) >= threshold:
                return float(self.position_probability[base][idx]) * float(self.position_weight[base])
        return 0.0

    def _lookup_content_probability(self, value, base):
        if float(value) < 0:
            return None
        for idx, threshold in enumerate(self.content_parameter):
            if float(value) >= threshold:
                return float(self.content_probability[base][idx]) * float(self.content_weight[base])
        return 0.0

    def value(self, dna):
        if len(dna) < 2:
            return 0.0

        dna = dna.upper()
        total_base = float(len(dna))
        phase_0 = dna[::3]
        phase_1 = dna[1::3]
        phase_2 = dna[2::3]

        phase_counts = {}
        for base in ("A", "C", "G", "T"):
            counts = [phase_0.count(base), phase_1.count(base), phase_2.count(base)]
            phase_counts[base] = counts

        fickett_score = 0.0
        for base in ("A", "C", "G", "T"):
            counts = phase_counts[base]
            content = float(sum(counts)) / total_base
            position = max(counts) / (min(counts) + 1.0)
            fickett_score += self._lookup_content_probability(content, base)
            fickett_score += self._lookup_position_probability(position, base)

        return fickett_score

