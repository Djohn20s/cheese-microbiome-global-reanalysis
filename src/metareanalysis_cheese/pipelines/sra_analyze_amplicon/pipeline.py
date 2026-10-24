# src/metareanalysis_cheese/pipelines/sra_analyze_amplicon/pipeline.py
from kedro.pipeline import Pipeline, node, pipeline
from .nodes import analyze_amplicon_from_manifest


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=analyze_amplicon_from_manifest,
                inputs=dict(
                    sra_manifest_staged="sra_manifest_staged",
                    params="params:sra_analyze_amplicon",
                ),
                outputs=dict(
                    sra_profile_amplicon_run_matrix="sra_profile_amplicon_run_matrix",
                    sra_profile_amplicon_tax_matrix="sra_profile_amplicon_tax_matrix",
                ),
                name="sra_analyze_amplicon",
                tags=["sra", "amplicon", "qiime2"],
            )
        ]
    )