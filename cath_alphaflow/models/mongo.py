from enum import Enum

from typing import Optional
from bson import ObjectId
from pydantic import ConfigDict, BaseModel, Field

from .beacons import UniprotSummary


class PyObjectId(ObjectId):  # pragma: no cover
    @classmethod
    def validate(cls, v):
        if not ObjectId.is_valid(v):
            raise ValueError(f"Invalid objectid: {v}")
        return ObjectId(v)


class AFFileType(Enum):  # pragma: no cover
    CONFIDENCE = "CONFIDENCE"
    MODEL_CIF = "MODEL_CIF"
    MODEL_PDB = "MODEL_PDB"
    PAE = "PAE"


class AFFile(BaseModel):  # pragma: no cover
    id: PyObjectId = Field(
        default_factory=PyObjectId, description="Unique ID", alias="_id"
    )
    dataset: str = Field(..., description="Dataset that this file belongs to")
    fileType: AFFileType = Field(..., description="AlphaFold file type, e.g. 'PAE'")
    fileName: str = Field(
        ..., description="AlphaFold file name, e.g. 'AF-P00520-F1-model_v4.cif'"
    )
    afVersion: int = Field(..., description="AlphaFold version number, e.g. 4")
    fragNum: int = Field(..., description="AlphaFold model fragment number, e.g. 1")
    uniprotAccession: str = Field(..., description="UniProt accession, e.g. 'P00520'")
    contents: str = Field(..., description="File contents")
    uniprot_summary: UniprotSummary = Field(
        ..., description="UniProt Summary (3D Beacons)"
    )
    # TODO[pydantic]: The following keys were removed: `json_encoders`.
    # Check https://docs.pydantic.dev/dev-v2/migration/#changes-to-config for more information.
    model_config = ConfigDict(
        populate_by_name=True,
        arbitrary_types_allowed=True,
        json_encoders={ObjectId: str},
        use_enum_values=True,
    )

    def __str__(self):
        return f"<AFFile id={str(self.id)} dataset={self.dataset} fileName={self.fileName} fileType={self.fileType} uniprotAccession={self.uniprotAccession}>"
