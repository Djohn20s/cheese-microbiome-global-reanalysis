from kedro.pipeline import Pipeline, node, pipeline
from .nodes import fetch_runinfo


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=fetch_runinfo,
                inputs="params:ncbi",
                outputs=["runinfo_raw", "runinfo_download_summary"],
                name="fetch_runinfo",
            )
        ]
    )
