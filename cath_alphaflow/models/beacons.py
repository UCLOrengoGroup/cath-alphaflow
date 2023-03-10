from enum import Enum

from typing import Optional
from bson import ObjectId
from pydantic import BaseModel, Field


class ModelCategory(Enum):  # pragma: no cover
    EXPERIMENTALLY_DETERMINED = "EXPERIMENTALLY DETERMINED"
    TEMPLATE_BASED = "TEMPLATE-BASED"
    AB_INITIO = "AB-INITIO"
    CONFORMATIONAL_ENSEMBLE = "CONFORMATIONAL ENSEMBLE"
    DEEP_LEARNING = "DEEP-LEARNING"


class ModelType(Enum):
    SINGLE = "single"
    COMPLEX = "complex"


class PyObjectId(ObjectId):  # pragma: no cover
    @classmethod
    def __get_validators__(cls):
        yield cls.validate

    @classmethod
    def validate(cls, v):
        if not ObjectId.is_valid(v):
            raise ValueError("Invalid objectid")
        return ObjectId(v)

    @classmethod
    def __modify_schema__(cls, field_schema):
        field_schema.update(type="string")


class ModelMetadata(BaseModel):  # pragma: no cover
    mappingAccession: str = Field(
        ...,
        description="Protein accession to which this "
        "model has been mapped, e.g. P38398",
    )
    mappingAccessionType: str = Field(
        ..., description="Type of accession code, e.g. uniprot"
    )
    start: int = Field(
        ...,
        description="The index of the first residue of the model according to UniProt"
        " sequence numbering, e.g. 1",
    )
    end: int = Field(
        ...,
        description="The index of the last residue of the model according to UniProt "
        "sequence numbering, e.g. 142",
    )
    modelCategory: ModelCategory = Field(
        ..., description="Category of the model, e.g. EXPERIMENTALLY DETERMINED"
    )
    modelType: ModelType = Field(
        ..., description="Category of the model type, e.g. single"
    )

    confidenceType: str = Field(
        ...,
        description="Type of the confidence measure. This is required for theoretical models.",
    )

    confidenceVersion: Optional[str] = Field(
        None,
        description=(
            "Version of confidence measure software used to calculate quality. This is "
            "required for theoretical models."
        ),
    )

    confidenceAvgLocalScore: str = Field(
        ...,
        description=(
            "Average of the confidence measures in the range of [0,1] for QMEANDisCo "
            "and [0,100] for pLDDT. Please contact 3D-Beacons developers if other estimates "
            "are to be added. This is required for theoretical models."
        ),
    )

    class Config:
        use_enum_values = True


class ModelEntry(BaseModel):  # pragma: no cover
    id: PyObjectId = Field(
        default_factory=PyObjectId, description="Unique ID", alias="_id"
    )
    entryId: str = Field(..., description="Model identifier")
    gene: str = Field(..., description="Gene")
    uniprotAccession: str = Field(..., description="UniProt accession, e.g. P00520")
    uniprotId: str = Field(..., description="UniProt identifier, e.g. ABL1_MOUSE")
    uniprotDescription: str = Field(..., description="Description for the protein")
    taxId: int = Field(
        ...,
        description="The NCBI taxonomy identifier (taxid) that points to a node of the"
        " taxonomy tree",
    )
    organismScientificName: str = Field(
        ..., description="The scientific name of the taxonomy node"
    )
    uniprotStart: int = Field(
        ...,
        description="The index of the first residue of the model according to UniProt"
        " sequence numbering, e.g. 1",
    )
    uniprotEnd: int = Field(
        ...,
        description="The index of the last residue of the model according to UniProt "
        "sequence numbering, e.g. 142",
    )
    modelCategory: ModelCategory = Field(
        ..., description="Category of the model, e.g. EXPERIMENTALLY DETERMINED"
    )
