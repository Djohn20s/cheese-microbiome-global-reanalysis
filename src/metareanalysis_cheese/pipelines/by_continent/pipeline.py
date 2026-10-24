# src/cheese_metagenome_meta/pipelines/by_continent/pipeline.py
from kedro.pipeline import Pipeline, node
from .nodes import save_by_continent_strategy_pies_map

def create_pipeline(**kwargs) -> Pipeline:
    return Pipeline([
        node(
            func=save_by_continent_strategy_pies_map,
            inputs=dict(
                by_continent="by_continent",
                output_path="params:by_continent.strategy_map_output",
            ),
            outputs="by_continent_strategy_world_map_png",
            name="by_continent_strategy_world_pies_map",
        )
    ])