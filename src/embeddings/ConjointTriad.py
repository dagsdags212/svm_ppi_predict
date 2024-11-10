from itertools import product
from typing import Iterable
import numpy as np
from AminoAcidClassifier import AminoAcidClassifier


class ConjointTriad:
    """Represent a protein sequence as a feature vector
    comprised of 7x7x7 elements."""

    CONTEXT = 3

    def __init__(
        self,
        seq: str,
        categories: Iterable = range(7),
        classifier=AminoAcidClassifier,
    ) -> None:
        self.seq = seq
        self.categories = categories

    def _create_all_triads(self) -> list[str]:
        cats = list(product(self.categories, repeat=self.CONTEXT))
        for i, cat in enumerate(cats):
            cats[i] = "".join(map(str, cat))
        return cats

    def classify_triad(self, triad: str) -> str:
        assert len(triad) == 3, "triad must have exactly THREE amino acids"
        label = []
        for aa in triad:
            if aa in ("A", "G", "V"):
                label.append(1)
            elif aa in ("I", "L", "F", "P"):
                label.append(2)
            elif aa in ("Y", "M", "T", "S"):
                label.append(3)
            elif aa in ("H", "N", "Q", "W"):
                label.append(4)
            elif aa in ("R", "K"):
                label.append(5)
            elif aa in ("D", "E"):
                label.append(6)
            elif aa == "C":
                label.append(7)
            else:
                raise ValueError("Invalid amino acid")
        return str("".join(map(str, label)))

    def _normalize(self, F):
        return (F - F.min()) / F.max()

    def vectorize(self):
        """Converts a protein sequence into a frequenct vector of size 343."""
        categories = self._create_all_triads()
        F = np.zeros(len(categories))
        catmap = {val: i for i, val in enumerate(categories)}
        for i in range(len(self.seq) - self.CONTEXT + 1):
            triad = self.seq[i : i + self.CONTEXT]
            label = self.classify_triad(triad)
            idx = catmap[label]
            F[idx] += 1
        F = self._normalize(F)
        return F
