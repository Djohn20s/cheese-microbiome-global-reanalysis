from kedro.pipeline import Pipeline, node, pipeline

from .nodes import analyze_wgs_from_manifest


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=analyze_wgs_from_manifest,
                inputs=dict(
                    sra_manifest_staged="sra_manifest_staged",
                    params="params:sra_analyze_wgs",
                ),
                outputs=dict(
                    sra_profile_wgs_run_matrix="sra_profile_wgs_run_matrix",
                    sra_profile_wgs_kraken_matrix="sra_profile_wgs_kraken_matrix",
                    sra_profile_wgs_motus_matrix="sra_profile_wgs_motus_matrix",
                ),
                name="sra_analyze_wgs",
                tags=["sra", "wgs"],
            )
        ]
    )