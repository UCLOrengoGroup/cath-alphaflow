from setuptools import setup
import os

VERSION = "0.1"


def get_long_description():
    with open(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "README.md"),
        encoding="utf8",
    ) as fp:
        return fp.read()


setup(
    name="cath-alphaflow",
    description="Workflow to assign CATH domains to protein structural models predicted by AlphaFold",
    long_description=get_long_description(),
    long_description_content_type="text/markdown",
    author=["Ian Sillitoe", "Nicola Bordin"],
    url="https://github.com/sillitoe/cath-alphaflow",
    project_urls={
        "Issues": "https://github.com/sillitoe/cath-alphaflow/issues",
        "CI": "https://github.com/sillitoe/cath-alphaflow/actions",
        "Changelog": "https://github.com/sillitoe/cath-alphaflow/releases",
    },
    license="Apache License, Version 2.0",
    version=VERSION,
    packages=["cath_alphaflow"],
    entry_points="""
        [console_scripts]
        cath-af-cli=cath_alphaflow.cli:cli
    """,
    install_requires=[
        "click",
        "google-cloud-storage",
    ],
    extras_require={"test": ["pytest"]},
    python_requires=">=3.7",
)
