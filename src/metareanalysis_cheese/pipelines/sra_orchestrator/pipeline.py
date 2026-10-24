from kedro.pipeline import Pipeline, node, pipeline

from .nodes import run_sra_orchestrator


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=run_sra_orchestrator,
                inputs=dict(
                    final_wgs="final_wgs",
                    final_amplicon="final_amplicon",
                    params_orchestrator="params:sra_orchestrator",
                    params_fetch="params:sra_fetch_unpack",
                    params_wgs="params:sra_analyze_wgs",
                    params_amplicon="params:sra_analyze_amplicon",
                ),
                outputs="sra_orchestrator_state",
                name="sra_orchestrator",
                tags=["sra", "orchestrator"],
            )
        ]
    )