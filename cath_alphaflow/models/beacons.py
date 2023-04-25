# generated by datamodel-codegen:
#   filename:  https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/openapi.json
#   timestamp: 2023-03-22T14:53:33+00:00

from __future__ import annotations

from enum import Enum
from typing import List, Optional, Union

from pydantic import Field
from pydantic import BaseModel as PyBaseModel

class BaseModel(PyBaseModel):
    class Config:
        use_enum_values = True


class ConfidenceType(Enum):
    pLDDT = "pLDDT"
    QMEANDisCo = "QMEANDisCo"
    ipTM_pTM = "ipTM+pTM"


class EnsembleSampleFormat(Enum):
    PDB = "PDB"
    MMCIF = "MMCIF"
    BCIF = "BCIF"


class EntityPolyType(Enum):
    CYCLIC_PSEUDO_PEPTIDE = "CYCLIC-PSEUDO-PEPTIDE"
    PEPTIDE_NUCLEIC_ACID = "PEPTIDE NUCLEIC ACID"
    POLYDEOXYRIBONUCLEOTIDE = "POLYDEOXYRIBONUCLEOTIDE"
    POLYDEOXYRIBONUCLEOTIDE_POLYRIBONUCLEOTIDE_HYBRID = (
        "POLYDEOXYRIBONUCLEOTIDE/POLYRIBONUCLEOTIDE HYBRID"
    )
    POLYPEPTIDE_D_ = "POLYPEPTIDE(D)"
    POLYPEPTIDE_L_ = "POLYPEPTIDE(L)"
    POLYRIBONUCLEOTIDE = "POLYRIBONUCLEOTIDE"
    OTHER = "OTHER"


class EntityType(Enum):
    BRANCHED = "BRANCHED"
    MACROLIDE = "MACROLIDE"
    NON_POLYMER = "NON-POLYMER"
    POLYMER = "POLYMER"
    WATER = "WATER"


class ExperimentalMethod(Enum):
    ELECTRON_CRYSTALLOGRAPHY = "ELECTRON CRYSTALLOGRAPHY"
    ELECTRON_MICROSCOPY = "ELECTRON MICROSCOPY"
    EPR = "EPR"
    FIBER_DIFFRACTION = "FIBER DIFFRACTION"
    FLUORESCENCE_TRANSFER = "FLUORESCENCE TRANSFER"
    INFRARED_SPECTROSCOPY = "INFRARED SPECTROSCOPY"
    NEUTRON_DIFFRACTION = "NEUTRON DIFFRACTION"
    X_RAY_POWDER_DIFFRACTION = "X-RAY POWDER DIFFRACTION"
    SOLID_STATE_NMR = "SOLID-STATE NMR"
    SOLUTION_NMR = "SOLUTION NMR"
    X_RAY_SOLUTION_SCATTERING = "X-RAY SOLUTION SCATTERING"
    THEORETICAL_MODEL = "THEORETICAL MODEL"
    X_RAY_DIFFRACTION = "X-RAY DIFFRACTION"
    HYBRID = "HYBRID"


class ExperimentalMethod1(Enum):
    ELECTRON_CRYSTALLOGRAPHY = "ELECTRON CRYSTALLOGRAPHY"
    ELECTRON_MICROSCOPY = "ELECTRON MICROSCOPY"
    EPR = "EPR"
    FIBER_DIFFRACTION = "FIBER DIFFRACTION"
    FLUORESCENCE_TRANSFER = "FLUORESCENCE TRANSFER"
    INFRARED_SPECTROSCOPY = "INFRARED SPECTROSCOPY"
    NEUTRON_DIFFRACTION = "NEUTRON DIFFRACTION"
    POWDER_DIFFRACTION = "POWDER DIFFRACTION"
    SOLID_STATE_NMR = "SOLID-STATE NMR"
    SOLUTION_NMR = "SOLUTION NMR"
    SOLUTION_SCATTERING = "SOLUTION SCATTERING"
    THEORETICAL_MODEL = "THEORETICAL MODEL"
    X_RAY_DIFFRACTION = "X-RAY DIFFRACTION"
    HYBRID = "HYBRID"


class IdentifierCategory(Enum):
    UNIPROT = "UNIPROT"
    RFAM = "RFAM"
    CCD = "CCD"
    SMILES = "SMILES"
    INCHI = "INCHI"
    INCHIKEY = "INCHIKEY"


class ModelCategory(Enum):
    EXPERIMENTALLY_DETERMINED = "EXPERIMENTALLY DETERMINED"
    TEMPLATE_BASED = "TEMPLATE-BASED"
    AB_INITIO = "AB-INITIO"
    CONFORMATIONAL_ENSEMBLE = "CONFORMATIONAL ENSEMBLE"


class ModelFormat(Enum):
    PDB = "PDB"
    MMCIF = "MMCIF"
    BCIF = "BCIF"


class ModelType(Enum):
    ATOMIC = "ATOMIC"
    DUMMY = "DUMMY"
    MIX = "MIX"


class OligomericState(Enum):
    MONOMER = "MONOMER"
    HOMODIMER = "HOMODIMER"
    HETERODIMER = "HETERODIMER"
    HOMO_OLIGOMER = "HOMO-OLIGOMER"
    HETERO_OLIGOMER = "HETERO-OLIGOMER"


class Residue(BaseModel):
    confidence: Optional[float] = Field(
        None,
        description="Confidence score in the range of [0,1]",
        example=0.99,
        title="Confidence",
    )
    model_residue_label: int = Field(
        ..., description="Model residue index", example=1, title="Model Residue Label"
    )
    uniprot_residue_number: int = Field(
        ...,
        description="UniProt residue index",
        example=1,
        title="Uniprot Residue Number",
    )


class Seqres(BaseModel):
    aligned_sequence: str = Field(
        ...,
        description="Sequence of the model",
        example="AAGTGHLKKKYT...",
        title="Aligned Sequence",
    )
    from_: int = Field(
        ...,
        alias="from",
        description="1-indexed first residue",
        example=32,
        title="From",
    )
    to: int = Field(..., description="1-indexed last residue", example=976, title="To")


class Template(BaseModel):
    template_id: str = Field(
        ...,
        description="Identifier of the template",
        example="2aqa",
        title="Template Id",
    )
    chain_id: str = Field(
        ...,
        description="Identifier of the chain of the template; this is label_asym_id in mmCIF",
        example="C",
        title="Chain Id",
    )
    template_sequence_identity: float = Field(
        ...,
        description="Sequence identity of the template with the  UniProt accession, in the range of [0,1]\n",
        example=0.97,
        title="Template Sequence Identity",
    )
    last_updated: str = Field(
        ...,
        description="Date of release of the last update in  the format of YYYY-MM-DD\n",
        example="2021-08-06",
        title="Last Updated",
    )
    provider: str = Field(
        ..., description="Provider of the template", example="PDB", title="Provider"
    )
    experimental_method: ExperimentalMethod1 = Field(
        ...,
        description="Experimental method used to determine the template",
        example="HYBRID",
    )
    resolution: float = Field(
        ...,
        description="Resolution of the template, in Angstrom",
        example=2.1,
        title="Resolution",
    )
    preferred_assembly_id: Optional[str] = Field(
        None,
        description="Identifier of the preferred assembly of the template",
        example="1",
        title="Preferred Assembly Id",
    )


class Uniprot(BaseModel):
    aligned_sequence: str = Field(
        ...,
        description="Sequence of the UniProt accession",
        example="AAGTGHLKKKYTAAGTGHLKKKYT...",
        title="Aligned Sequence",
    )
    from_: int = Field(
        ...,
        alias="from",
        description="1-indexed first residue",
        example=23,
        title="From",
    )
    to: int = Field(..., description="1-indexed last residue", example=868, title="To")


class UniprotEntry(BaseModel):
    ac: str = Field(..., description="UniProt accession", example="P00520", title="Ac")
    id: Optional[str] = Field(
        None, description="UniProt identifier", example="ABL1_MOUSE", title="Id"
    )
    uniprot_checksum: Optional[str] = Field(
        None,
        description="CRC64 checksum of the UniProt sequence",
        example="5F9BA1D4C7DE6925",
        title="Uniprot Checksum",
    )
    sequence_length: Optional[int] = Field(
        None,
        description="Length of the UniProt sequence",
        example=76,
        title="Sequence Length",
    )
    segment_start: Optional[int] = Field(
        None,
        description="1-indexed first residue of the UniProt sequence segment",
        example=1,
        title="Segment Start",
    )
    segment_end: Optional[int] = Field(
        None,
        description="1-indexed last residue of the UniProt sequence segment",
        example=86,
        title="Segment End",
    )


class ValidationError(BaseModel):
    loc: List[Union[str, int]] = Field(..., title="Location")
    msg: str = Field(..., title="Message")
    type: str = Field(..., title="Error Type")


class Entity(BaseModel):
    entity_type: EntityType = Field(
        ...,
        description="The type of the molecular entity; similar to _entity.type in mmCIF",
        example="POLYMER",
    )
    entity_poly_type: Optional[EntityPolyType] = Field(
        None,
        description="The type of the molecular entity; similar to _entity_poly.type in mmCIF",
        example="PEPTIDE NUCLEIC ACID",
    )
    identifier: Optional[str] = Field(
        None,
        description="Identifier of the molecule",
        example="Q13033",
        title="Identifier",
    )
    identifier_category: Optional[IdentifierCategory] = Field(
        None, description="Category of the identifier", example="UNIPROT"
    )
    description: str = Field(
        ...,
        description="A textual label of the molecule",
        example="Striatin-3",
        title="Description",
    )
    chain_ids: List[str] = Field(..., title="Chain Ids")


class HTTPValidationError(BaseModel):
    detail: Optional[List[ValidationError]] = Field(None, title="Detail")


class Segment(BaseModel):
    templates: Optional[List[Template]] = Field(
        None,
        description="Information on the template(s) used for the model",
        title="Templates",
    )
    seqres: Seqres = Field(
        ..., description="Information on the sequence of the model", title="Seqres"
    )
    uniprot: Uniprot
    residues: List[Residue] = Field(..., title="Residues")


class SummaryItems(BaseModel):
    model_identifier: str = Field(
        ...,
        description="Identifier of the model, such as PDB id",
        example="8kfa",
        title="Model Identifier",
    )
    model_category: ModelCategory = Field(
        ..., description="Category of the model", example="TEMPLATE-BASED"
    )
    model_url: str = Field(
        ...,
        description="URL of the model coordinates",
        example="https://www.ebi.ac.uk/pdbe/static/entry/1t29_updated.cif",
        title="Model Url",
    )
    model_format: ModelFormat = Field(
        ..., description="File format of the coordinates", example="MMCIF"
    )
    model_type: Optional[ModelType] = Field(
        None,
        description="Defines if the coordinates are atomic-level or contains dummy atoms (e.g. SAXS models), or a mix of both (e.g. hybrid models)\n",
        example="ATOMIC",
    )
    model_page_url: Optional[str] = Field(
        None,
        description="URL of a web page of the data provider that show the model",
        example="https://alphafold.ebi.ac.uk/entry/Q5VSL9",
        title="Model Page Url",
    )
    provider: str = Field(
        ...,
        description="Name of the model provider",
        example="SWISS-MODEL",
        title="Provider",
    )
    number_of_conformers: Optional[float] = Field(
        None,
        description="The number of conformers in a conformational ensemble",
        example=42,
        title="Number Of Conformers",
    )
    ensemble_sample_url: Optional[str] = Field(
        None,
        description="URL of a sample of conformations from a conformational ensemble",
        example="https://proteinensemble.org/api/ensemble_sample/PED00001e001",
        title="Ensemble Sample Url",
    )
    ensemble_sample_format: Optional[EnsembleSampleFormat] = Field(
        None,
        description="File format of the sample coordinates, e.g. PDB",
        example="PDB",
    )
    created: str = Field(
        ...,
        description="Date of release of model generation in the format of YYYY-MM-DD",
        example="2021-12-21",
        title="Created",
    )
    sequence_identity: float = Field(
        ...,
        description="Sequence identity in the range of [0,1] of the model to the UniProt sequence\n",
        example=0.97,
        title="Sequence Identity",
    )
    uniprot_start: int = Field(
        ...,
        description="1-indexed first residue of the model according to UniProt sequence numbering\n",
        example=1,
        title="Uniprot Start",
    )
    uniprot_end: int = Field(
        ...,
        description="1-indexed last residue of the model according to UniProt sequence numbering\n",
        example=142,
        title="Uniprot End",
    )
    coverage: float = Field(
        ...,
        description="Fraction in range of [0, 1] of the UniProt sequence covered by the model.  This is calculated as (uniprot_end - uniprot_start + 1) / uniprot_sequence_length\n",
        example=0.4,
        title="Coverage",
    )
    experimental_method: Optional[ExperimentalMethod] = Field(
        None,
        description="Experimental method used to determine the structure, if applicable",
    )
    resolution: Optional[float] = Field(
        None,
        description="The resolution of the model in Angstrom, if applicable",
        example=1.4,
        title="Resolution",
    )
    confidence_type: Optional[ConfidenceType] = Field(
        None,
        description="Type of the confidence measure. This is required for  theoretical models.\n",
        example="QMEANDisCo",
    )
    confidence_version: Optional[str] = Field(
        None,
        description="Version of confidence measure software used to calculate quality. This is required for theoretical models.\n",
        example="v1.0.2",
        title="Confidence Version",
    )
    confidence_avg_local_score: Optional[float] = Field(
        None,
        description="Average of the confidence measures in the range of [0,1] for QMEANDisCo  and [0,100] for pLDDT. Please contact 3D-Beacons developers if other  estimates are to be added. This is required for theoretical models.\n",
        example=0.95,
        title="Confidence Avg Local Score",
    )
    oligomeric_state: Optional[OligomericState] = Field(
        None, description="Oligomeric state of the model", example="MONOMER"
    )
    preferred_assembly_id: Optional[str] = Field(
        None,
        description="Identifier of the preferred assembly in the model",
        example="1A",
        title="Preferred Assembly Id",
    )
    entities: List[Entity] = Field(
        ..., description="A list of molecular entities in the model", title="Entities"
    )


class Chain(BaseModel):
    chain_id: str = Field(..., title="Chain Id")
    segments: Optional[List[Segment]] = Field(None, title="Segments")


class Chains(BaseModel):
    __root__: List[Chain] = Field(..., title="Chains")


class Detailed(BaseModel):
    summary: SummaryItems
    chains: Chains


class Overview(BaseModel):
    summary: SummaryItems


class UniprotDetails(BaseModel):
    uniprot_entry: Optional[UniprotEntry] = None
    structures: Optional[List[Detailed]] = Field(None, title="Structures")


class UniprotSummary(BaseModel):
    uniprot_entry: Optional[UniprotEntry] = None
    structures: Optional[List[Overview]] = Field(None, title="Structures")
