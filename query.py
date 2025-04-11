import re
import time
from datetime import datetime, timezone
from typing import Any, Iterable, NamedTuple, TypeVar

import pandas as pd
from lamin_utils import logger
import requests
from requests.exceptions import SSLError, RequestException
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
)

DataSetClass = TypeVar("DataSetClass", bound=NamedTuple)


class scRNAseqDataset(NamedTuple):
    raw_expr: str | None = None
    expr: str | None = None
    secondary_analysis: str | None = None
    scvelo: str | None = None


class BulkseqDataset(NamedTuple):
    expression_matrices: str


def safe_head(url: str, retries: int = 3, delay: float = 0.5) -> bool:
    for attempt in range(retries):
        try:
            response = requests.head(url, timeout=5)
            return response.ok
        except SSLError:
            logger.error(f"SSL error for {url}. (retry {attempt + 1}/{retries})")
            time.sleep(delay)
        except RequestException:
            logger.error(f"Request error for {url}. (retry {attempt + 1}/{retries})")
            time.sleep(delay)
    return False


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
                if safe_head(potential_url):
                    urls[file_type] = potential_url

    variant_mapping = {
        "raw_expr": {"raw_expr.h5ad", "out.h5ad"},
        "scvelo": {"scvelo.h5ad", "scvelo_annotated.h5ad"},
    }

    kwargs = {}
    for field in dataset_class.__annotations__:
        variants = variant_mapping.get(field, {f"{field}.h5ad", f"{field}.h5"})
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
) -> pd.DataFrame:
    data = []
    valid_uuids = []

    metadata_lookup = hubmap_metadata.set_index("uuid")

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

        for uuid in metadata_lookup.index:
            progress.update(task, uuid=uuid)

            try:
                dataset_info = get_dataset_info(uuid)[0]
                donor_metadata = (
                    dataset_info.get("donor", {})
                    .get("metadata", {})
                    .get("organ_donor_data", [])
                )

                dataset_urls = get_dataset_urls(
                    uuid,
                    file_types=file_types,
                    dataset_class=dataset_class,
                )

                if dataset_urls is None:
                    logger.warning(
                        f"No usable files for uuid {uuid}, HuBMAP ID {dataset_info.get('hubmap_id')}."
                    )
                    continue

                # Dynamically collect all URLs from dataset class
                urls = {
                    f"{field}_url": getattr(dataset_urls, field) or ""
                    for field in dataset_class.__annotations__
                }

                row_meta = metadata_lookup.loc[uuid]
                row = {
                    "assay": row_meta["assay_type"],
                    "rnaseq_assay_method": row_meta["rnaseq_assay_method"],
                    "title": dataset_info.get("title", ""),
                    "group_name": dataset_info.get("group_name", ""),
                    "consortium": "HuBMAP",
                    "doi": dataset_info.get("registered_doi", ""),
                    "publication_date": datetime.fromtimestamp(
                        dataset_info.get("published_timestamp", 0) / 1000,
                        tz=timezone.utc,
                    ).strftime("%Y-%m-%d")
                    if dataset_info.get("published_timestamp")
                    else "",
                    "status": dataset_info.get("data_access_level", ""),
                    "dataset_type": dataset_info.get("dataset_type", ""),
                    "processing": "raw",
                    "organ": next(
                        (
                            s.get("organ", "")
                            for s in dataset_info.get("origin_samples", [])
                            if s.get("organ")
                        ),
                        "",
                    ),
                    "sample_category": next(
                        (
                            s.get("sample_category", "")
                            for s in dataset_info.get("source_samples", [])
                            if s.get("sample_category")
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
                            i.get("data_value", "")
                            for i in donor_metadata
                            if i.get("grouping_concept_preferred_term")
                            == "Body Mass Index"
                        ),
                        "",
                    ),
                    "age": next(
                        (
                            i.get("data_value", "")
                            for i in donor_metadata
                            if i.get("grouping_concept_preferred_term") == "Age"
                        ),
                        "",
                    ),
                    "ethnicity": next(
                        (
                            i.get("data_value", "")
                            for i in donor_metadata
                            if i.get("grouping_concept_preferred_term") == "Race"
                        ),
                        "",
                    ),
                    "sex": next(
                        (
                            i.get("data_value", "")
                            for i in donor_metadata
                            if i.get("grouping_concept_preferred_term") == "Sex"
                        ),
                        "",
                    ),
                    "diseases": [
                        i.get("data_value", "")
                        for i in donor_metadata
                        if i.get("grouping_concept_preferred_term") == "Medical History"
                    ]
                    or ["normal"],
                    "donor_id": dataset_info.get("donor", {}).get("hubmap_id", ""),
                    "sample_id": next(
                        (
                            s.get("hubmap_id", "")
                            for s in dataset_info.get("source_samples", [])
                            if s.get("hubmap_id")
                        ),
                        "",
                    ),
                    "ancestor_id": dataset_info.get("immediate_ancestor_ids", [])[0]
                    if dataset_info.get("immediate_ancestor_ids")
                    else "",
                    **urls,
                }

                data.append(row)
                valid_uuids.append(uuid)

            except Exception as e:
                logger.error(f"Error processing uuid {uuid}: {e}")

            progress.update(task, advance=1)

    df = pd.DataFrame(data, index=valid_uuids)
    df.index.name = "uuid"
    return df


def extract_hubmap_uid(url: str) -> str:
    match = re.search(r"hubmapconsortium\.org/([a-f0-9]{32})", url)
    if not match:
        raise ValueError("Invalid HubMap URL format")
    return match.group(1)
