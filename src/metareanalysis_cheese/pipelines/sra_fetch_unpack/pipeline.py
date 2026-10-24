# src/metareanalysis_cheese/pipelines/sra_fetch_unpack/pipeline.py
from kedro.pipeline import Pipeline, node, pipeline

from .nodes import build_runs_table, fetch_unpack_stage


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=build_runs_table,
                inputs=dict(final_wgs="final_wgs", final_amplicon="final_amplicon"),
                outputs="sra_runs_table",
                name="sra_build_runs_table",
            ),
            node(
                func=fetch_unpack_stage,
                inputs=dict(runs_table="sra_runs_table", params="params:sra_fetch_unpack"),
                outputs="sra_manifest_staged",
                name="sra_fetch_unpack_stage",
            ),
        ]
    )