#!/usr/bin/env python3
# coding: utf-8

from enum import Enum
from typing import Dict, List, Optional, Tuple
import pandas as pd


import argparse
import yaml


def read_metadata(filepath):
    return pd.read_csv(
        filepath,
        sep="\t",
        usecols=(
            "strain",
            "date",
            "country",
            "division",
            "QC_overall_status",
            "Nextstrain_clade",
            "Nextclade_pango",
        ),
        dtype={
            "country": "category",
            "division": "category",
            "Nextstrain_clade": "category",
            "Nextclade_pango": "category",
        },
        chunksize=100000,
    )


class SpatialLevel(Enum):
    COUNTRY = "country"
    DIVISION = "division"

    def __str__(self):
        return str(self.value)


def filter_location(
    metadata: pd.DataFrame, spatial_level: SpatialLevel, locations: List[str]
) -> pd.DataFrame:
    return metadata[metadata[str(spatial_level)].isin(locations)].copy()


def filter_ambigious_dates(metadata: pd.DataFrame) -> pd.DataFrame:
    # Remove ambigious dates
    unambiguous_dates = (metadata["date"] != "?") & (
        metadata["date"].str.count("-") == 2
    )
    unambiguous_dates = unambiguous_dates & ~(
        metadata["date"].str.contains("X", na=True)
    )
    metadata = metadata[unambiguous_dates].copy()
    return metadata


def filter_date_range(
    metadata: pd.DataFrame, start_date: str, end_date: str
) -> pd.DataFrame:
    # Filter to date range
    metadata["date"] = pd.to_datetime(metadata["date"])
    date_since_start_date = metadata["date"] >= start_date
    date_before_end_date = metadata["date"] <= end_date
    return metadata[(date_since_start_date) & (date_before_end_date)].copy()


def filter_quality_control(metadata: pd.DataFrame) -> pd.DataFrame:
    return metadata[metadata["QC_overall_status"] != "bad"].copy()


def filter_has_clade(metadata: pd.DataFrame) -> pd.DataFrame:
    return metadata[~pd.isnull(metadata["Nextstrain_clade"])].copy()


# Borrowing from John here
def clean_metadata(
    metadata,
    spatial_levels: List[SpatialLevel],
    locations: Dict[SpatialLevel, List[str]],
    start_date: str,
    end_date: str,
    quality_control: bool = True,
):
    # Filter to locations at spatial levels of interest
    for spl in spatial_levels:
        metadata = filter_location(metadata, spl, locations[spl])

    metadata = filter_ambigious_dates(metadata)
    metadata = filter_date_range(metadata, start_date, end_date)

    # Remove low quality sequences
    if quality_control:
        metadata = filter_quality_control(metadata)

    # Filer sequences without a Nextstrain clade
    metadata = filter_has_clade(metadata)
    return metadata


def map_variant(clade: str, lineage: str, variant_map: Dict[str, str]) -> str:
    # If smaller level is present in mapping, use that
    if lineage in variant_map:
        return variant_map[lineage]
    elif clade in variant_map:
        return variant_map[clade]
    else:
        return "other"


def map_variants(
    metadata: pd.DataFrame, variant_map: Dict[str, str]
) -> pd.DataFrame:
    metadata["variant"] = list(
        map(
            lambda clade, lineage: map_variant(clade, lineage, variant_map),
            metadata["Nextstrain_clade"],
            metadata["Nextclade_pango"],
        )
    )

    return metadata


def count_metadata(
    metadata: pd.DataFrame,
    spatial_level: SpatialLevel,
    obs_date: Optional[str] = None,
) -> pd.DataFrame:
    # metadata = metadata[metadata["date_submitted"] < obs_date].copy()

    return (
        metadata.groupby(["date", str(spatial_level), "variant"])["strain"]
        .count()
        .reset_index()
        .rename(
            columns={"strain": "sequences", str(spatial_level): "location"}
        )
        .sort_values(["location", "variant", "date"])
    )


def get_sequence_counts(
    metadata_path: str,
    spatial_levels: List[SpatialLevel],
    locations: Dict[SpatialLevel, List[str]],
    start_date: str,
    end_date: str,
    variant_map: Dict[str, str],
    export_path: Optional[str] = None,
):

    # Load metadata
    raw_metadata_reader = read_metadata(metadata_path)

    metadata_chunks = [
        clean_metadata(
            chunk,
            spatial_levels=spatial_levels,
            locations=locations,
            start_date=start_date,
            end_date=end_date,
        )
        for chunk in raw_metadata_reader
    ]

    metadata = pd.concat(metadata_chunks, ignore_index=True)

    metadata = map_variants(metadata, variant_map)

    seq_counts = count_metadata(metadata, spatial_levels[-1])
    if export_path is not None:
        seq_counts.to_csv(export_path, sep="\t", index=False)
    return seq_counts


class SeqCountsConfig:
    def __init__(self, config_path: str):
        # Load config
        self.config_path = config_path
        self.config = self.read_config(self.config_path)
        self.unpack_config()

    def read_config(self, path):
        with open(path, "r") as file:
            config = yaml.safe_load(file)
        return config

    def get_config_val(self, val: str):
        return self.config[val]

    @staticmethod
    def get_spatial_filters(
        config: Dict,
    ) -> Tuple[List[SpatialLevel], Dict[SpatialLevel, List[str]]]:
        spatial_level_values = tuple(item.value for item in SpatialLevel)
        spatial_levels = []
        locations = dict()
        for filter, values in config["filters"].items():
            if filter in spatial_level_values:
                spatial_levels.append(SpatialLevel(filter))
                locations[SpatialLevel(filter)] = values
        return spatial_levels, locations

    @staticmethod
    def get_date_filters(config: Dict):
        date_range = config["filters"]["date_range"]
        return date_range["start"], date_range["end"]

    def unpack_config(self):
        self.spatial_levels, self.locations = self.get_spatial_filters(
            self.config
        )
        self.start_date, self.end_date = self.get_date_filters(self.config)
        self.variant_map = self.load_variant_map(
            self.config["variant_map_path"]
        )
        return self

    @staticmethod
    def load_variant_map(variant_map_path: str):
        # Opening JSON file
        variant_map = pd.read_csv(variant_map_path, sep="\t")
        return {x["input"]: x["variant"] for _, x in variant_map.iterrows()}

    def get_sequence_counts(
        self, metadata_path: str, export_path: Optional[str] = None
    ) -> pd.DataFrame:
        return get_sequence_counts(
            metadata_path,
            self.spatial_levels,
            self.locations,
            self.start_date,
            self.end_date,
            self.variant_map,
            export_path,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generating variant sequence counts"
    )
    parser.add_argument("--config", help="path to config file")
    parser.add_argument(
        "--metadata-path",
        help="File path to sequence metadata.",
    )
    parser.add_argument(
        "--export-path",
        help="Path to exported file.",
    )

    args = parser.parse_args()

    config = SeqCountsConfig(args.config)

    get_sequence_counts(
        args.metadata_path,
        config.spatial_levels,
        config.locations,
        config.start_date,
        config.end_date,
        config.variant_map,
        args.export_path,
    )
