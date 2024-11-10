class AminoAcidClassifier:
    """Places an amino acid into one of seven groups
    based on physicochemical properties."""

    AA_MAPPING = {
        "A": "0",
        "G": "0",
        "V": "0",
        "C": "1",
        "I": "2",
        "L": "2",
        "F": "2",
        "P": "2",
        "Y": "3",
        "M": "3",
        "T": "3",
        "S": "3",
        "H": "4",
        "N": "4",
        "Q": "4",
        "W": "4",
        "R": "5",
        "K": "5",
        "D": "6",
        "E": "6",
    }

    def __init__(self, seq: str) -> None:
        self.seq = seq

    def __call__(self) -> str:
        """Map each amino acid into a category."""
        return "".join([self.AA_MAPPING[n] for n in self.seq.upper()])
