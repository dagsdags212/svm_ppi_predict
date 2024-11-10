from enum import StrEnum, auto
import numpy as np
import pandas as pd


def load_aa_properties() -> pd.DataFrame:
    """Load table of physicochemical properties for amino acids."""
    aa_properties = "../data/aa_properties.csv"
    df = pd.read_csv(aa_properties, header=0)
    df.set_index("aa", inplace=True)
    df.astype({col: "float32" for col in df.columns})
    return df


def load_normalized_properties(P: pd.DataFrame) -> pd.DataFrame:
    """Normalize values to have μ=0 and σ=1."""
    return (P - P.mean(axis=0)) / P.std(axis=0)


class AminoAcidProperty(StrEnum):
    """
    Wrapper for amino acid properties to be considered when computing
    the AC matrix.

    Legend:
        H : hydrophobicity
        VSC :
    """

    H = auto()
    VSC = auto()
    P1 = auto()
    P2 = auto()
    SASA = auto()
    NCISC = auto()


class AutoCovariance:
    """
    Converts a protein sequence into a P x N feature vector when
    called where P is the number of amino acid properties and
    N is the sequence length.

    Parameters:
        seq (str) : amino acid sequence
        lag (int) : distance between amino acids for computing variance

    Return
        M (ndarray) : a P x N matrix comprising of float values
    """

    _PROP_TABLE = load_aa_properties()

    def __init__(self, seq: str, lag: int):
        self.seq = seq
        self.lag = lag
        self.prop_table = load_normalized_properties(self._PROP_TABLE)
        self.M = self._populate_aa_matrix()

    def __call__(self):
        N = []
        for j in range(len(AminoAcidProperty)):
            Nj = []
            for lag in range(self.lag):
                Nj.append(self.compute_ac(lag, j))
            N.append(Nj)
        return np.matrix(N)

    def _populate_aa_matrix(self) -> np.ndarray:
        """Computes for the sum of all properties across the sequence."""
        # indexes the property table
        aa2prop = lambda aa, prop: self.prop_table.loc[aa, prop]
        # updates the matrix with property values
        M = np.zeros((len(AminoAcidProperty), len(self.seq)))
        for aa_idx, aa in enumerate(self.seq):
            for prop_idx, prop in enumerate(AminoAcidProperty):
                M[prop_idx, aa_idx] = aa2prop(aa, prop)
        return M

    def compute_ac(self, lag: int, j: int) -> float:
        """Computes for the autocovariance between amino acids
        that are `lag` units apart for the jth property."""
        _, n = self.M.shape
        rowsums = self.M.sum(axis=0)

        ac = 0
        for i in range(n - lag):
            left = self.M[j, i] - (rowsums[j] / n)
            right = self.M[j, i + lag] - (rowsums[j] / n)
            ac += left * right
        ac_norm = ac / (n - lag)
        return ac_norm
