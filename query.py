import requests
from typing import NamedTuple, Any


class HubmapSinglecellDataset(NamedTuple):
    raw_expr: str
    expr: str
    secondary_analysis: str
    scvelo: str | None = None


def get_dataset_urls(uuid: str) -> HubmapSinglecellDataset:
    url = "https://search.api.hubmapconsortium.org/v3/param-search/datasets"
    params = {"uuid": uuid}
    response = requests.get(url, params=params)
    response.raise_for_status()
    data = response.json()

    if not data or not isinstance(data, list):
        raise ValueError("Invalid response format from HubMAP API")

    urls = {}
    file_types = [
        "raw_expr.h5ad",
        "expr.h5ad",
        "secondary_analysis.h5ad",
        "scvelo.h5ad",
    ]

    for descendant in data[0].get("descendant_ids", []):
        for file_type in file_types:
            if file_type not in urls:
                potential_url = (
                    f"https://assets.hubmapconsortium.org/{descendant}/{file_type}"
                )
                head_response = requests.head(potential_url)
                if head_response.ok:
                    urls[file_type] = potential_url

    return HubmapSinglecellDataset(
        raw_expr=urls.get("raw_expr.h5ad"),
        expr=urls.get("expr.h5ad"),
        secondary_analysis=urls.get("secondary_analysis.h5ad"),
        scvelo=urls.get("scvelo.h5ad"),
    )


def get_dataset_info(uuid: str) -> Any:
    url = "https://search.api.hubmapconsortium.org/v3/param-search/datasets"
    params = {"uuid": uuid}
    response = requests.get(url, params=params)

    if response.status_code == 200:
        return response.json()
    else:
        response.raise_for_status()

    return response
