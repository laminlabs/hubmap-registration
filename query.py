import requests
from typing import NamedTuple, Any, TypeVar, Iterable
import re
import pandas as pd

from rich.progress import (
    Progress,
    SpinnerColumn,
    TimeElapsedColumn,
    BarColumn,
    TextColumn,
)

from datetime import datetime, timezone

DataSetClass = TypeVar("DataSetClass", bound=NamedTuple)


class scRNAseqDataset(NamedTuple):
    raw_expr: str
    expr: str | None = None
    secondary_analysis: str | None = None
    scvelo: str | None = None


class BulkseqDataset(NamedTuple):
    expression_matrices: str


def get_dataset_urls(
    uuid: str, file_types: Iterable[str], dataset_class: DataSetClass
) -> DataSetClass | None:
    """Gets direct URLs to datasets that can be registered as Artifacts."""
    url = "https://search.api.hubmapconsortium.org/v3/param-search/datasets"
    params = {"uuid": uuid}
    response = requests.get(url, params=params)
    response.raise_for_status()
    data = response.json()

    if not data or not isinstance(data, list):
        raise ValueError("Invalid response format from HubMAP API")

    urls = {}
    for descendant in data[0].get("descendant_ids", []):
        for file_type in file_types:
            if file_type not in urls:
                potential_url = (
                    f"https://assets.hubmapconsortium.org/{descendant}/{file_type}"
                )
                if requests.head(potential_url).ok:
                    urls[file_type] = potential_url

    variant_mapping = {
        "raw_expr": {"raw_expr.h5ad", "out.h5ad"},
        "scvelo": {"scvelo.h5ad", "scvelo_annotated.h5ad"},
    }

    kwargs = {}
    for field in dataset_class.__annotations__:
        variants = variant_mapping.get(field, {f"{field}.h5ad"})
        for v in variants:
            if v in urls:
                kwargs[field] = urls[v]
                break
        else:
            kwargs[field] = None

    if all(value is None for value in kwargs.values()):
        # if no valid files found for this uuid — skip it
        return None

    return dataset_class(**kwargs)  # type: ignore


def get_dataset_info(uuid: str) -> Any:
    """Fetches all associated metadata of a dataset UUID.

    Does not fetch detailed Donors or Samples.
    See also https://docs.hubmapconsortium.org/param-search/
    """
    url = "https://search.api.hubmapconsortium.org/v3/param-search/datasets"
    params = {"uuid": uuid}
    response = requests.get(url, params=params)

    if response.status_code == 200:
        return response.json()
    else:
        response.raise_for_status()

    return response


def create_hubmap_metadata_df(
    hubmap_metadata: pd.DataFrame,
    file_types: Iterable[str],
    dataset_class: DataSetClass,
):
    data = []
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TextColumn("[blue]{task.fields[uuid]}"),
        TimeElapsedColumn(),
    ) as progress:
        """Fetches dataset URLs and metadata using the API to collect it in a DataFrame."""
        task = progress.add_task(
            "[cyan]Processing datasets...", total=len(hubmap_metadata), uuid=""
        )
        for uuid in hubmap_metadata["uuid"].values:
            progress.update(task, uuid=uuid)
            dataset_info = get_dataset_info(uuid)[0]
            donor_metadata = (
                dataset_info.get("donor", {})
                .get("metadata", {})
                .get("organ_donor_data", [])
            )

            dataset_urls = get_dataset_urls(
                uuid, file_types=file_types, dataset_class=dataset_class
            )
            urls = {}
            for attr in dir(dataset_urls):
                if not attr.startswith("_") and not callable(
                    getattr(dataset_urls, attr)
                ):
                    value = getattr(dataset_urls, attr)
                    urls[f"{attr}_url"] = value if value is not None else ""

            row = {
                "assay": hubmap_metadata.loc[
                    hubmap_metadata["uuid"] == uuid, "assay_type"
                ].iloc[0],
                "rnaseq_assay_method": hubmap_metadata.loc[
                    hubmap_metadata["uuid"] == uuid, "rnaseq_assay_method"
                ].iloc[0],
                "title": dataset_info.get("title", ""),
                "group_name": dataset_info.get("group_name", ""),
                "consortium": "HuBMAP",
                "doi": dataset_info.get("registered_doi", ""),
                "publication_date": datetime.fromtimestamp(
                    dataset_info.get("published_timestamp", 0) / 1000, tz=timezone.utc
                ).strftime("%Y-%m-%d")
                if dataset_info.get("published_timestamp")
                else "",
                "status": dataset_info.get("data_access_level", ""),
                "dataset_type": dataset_info.get("dataset_type", ""),
                "processing": "raw",
                "organ": next(
                    (
                        sample.get("organ", "")
                        for sample in dataset_info.get("origin_samples", [])
                        if sample.get("organ")
                    ),
                    "",
                ),
                "sample_category": next(
                    (
                        sample.get("sample_category", "")
                        for sample in dataset_info.get("source_samples", [])
                        if sample.get("sample_category")
                    ),
                    "",
                ),
                "analyte_class": next(
                    (
                        a
                        for a in [
                            "RNA",
                            "Protein",
                            "DNA",
                            "Metabolite",
                            "Lipid",
                            "Nucleic acid + protein",
                            "Endogenous fluorophore",
                            "Polysaccharide",
                            "Peptide",
                            "DNA + RNA",
                            "Lipid + metabolite",
                        ]
                        if a in dataset_info.get("dataset_type", "")
                        or any(
                            a in d.get("dataset_type", "")
                            for d in dataset_info.get("descendants", [])
                        )
                    ),
                    "",
                ),
                "bmi": next(
                    (
                        item.get("data_value", "")
                        for item in donor_metadata
                        if item.get("grouping_concept_preferred_term")
                        == "Body Mass Index"
                    ),
                    "",
                ),
                "age": next(
                    (
                        item.get("data_value", "")
                        for item in donor_metadata
                        if item.get("grouping_concept_preferred_term") == "Age"
                    ),
                    "",
                ),
                "ethnicity": next(
                    (
                        item.get("data_value", "")
                        for item in donor_metadata
                        if item.get("grouping_concept_preferred_term") == "Race"
                    ),
                    "",
                ),
                "sex": next(
                    (
                        item.get("data_value", "")
                        for item in donor_metadata
                        if item.get("grouping_concept_preferred_term") == "Sex"
                    ),
                    "",
                ),
                "diseases": [
                    item.get("data_value", "")
                    for item in donor_metadata
                    if item.get("grouping_concept_preferred_term") == "Medical History"
                ]
                or ["normal"],
                "donor_id": dataset_info.get("donor", {}).get("hubmap_id", ""),
                "sample_id": next(
                    (
                        sample.get("hubmap_id", "")
                        for sample in dataset_info.get("source_samples", [])
                        if sample.get("hubmap_id")
                    ),
                    "",
                ),
                "ancestor_id": dataset_info.get("immediate_ancestor_ids", [])[0]
                if dataset_info.get("immediate_ancestor_ids")
                else "",  # Always the first ancestor ID
                **urls,
            }
            data.append(row)
            progress.update(task, advance=1)

    df = pd.DataFrame(data)

    return df


def extract_hubmap_uid(url: str) -> str:
    match = re.search(r"hubmapconsortium\.org/([a-f0-9]{32})", url)
    if not match:
        raise ValueError("Invalid HubMap URL format")
    return match.group(1)
