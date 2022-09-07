from dataclasses import dataclass


@dataclass
class PredictedCathDomain:
    """
    Holds data on a PredictedCathDomain (from Gene3D)
    """

    uniprot_acc: str
    sequence_md5: str
    gene3d_domain_id: str
    bitscore: str
    chopping: str
