from kedro.pipeline import Pipeline, node, pipeline

from .nodes import summarize_stage5_metadata


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=summarize_stage5_metadata,
                inputs=dict(
                    wgs="final_wgs",
                    amplicon="final_amplicon",
                ),
                outputs="stage5_metadata_summary",
                name="stage5_metadata_summary_tables",
            )
        ]
    )