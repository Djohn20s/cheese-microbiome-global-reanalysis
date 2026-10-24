"""
Pipeline: process_metadata (Kedro)

0) add_continent_by_centername: join runinfo + center_continent_map -> runinfo_with_continent
1) Stage1 run-level curation -> summary
2) Stage3 bioproject curation -> manual template -> gate for FILLED -> merge manual -> final summary
3) Split final into WGS / AMPLICON / OTHERS
"""

from kedro.pipeline import Pipeline, node, pipeline

from .nodes import (
    add_continent_by_centername,
    stage1_runlevel_curation,
    summarize_runinfo_markdown,
    stage3_bioproject_curation,
    make_manual_bioproject_template,
    require_and_load_manual_bioprojects_filled,   # ✅ NEW
    apply_manual_and_merge_kept,
    split_strategy,
)


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=add_continent_by_centername,
                inputs=dict(
                    runinfo="runinfo_raw",
                    center_map="center_continent_map",
                ),
                outputs="runinfo_with_continent",
                name="add_continent_by_centername",
            ),
            node(
                func=stage1_runlevel_curation,
                inputs=dict(
                    runinfo="runinfo_with_continent",
                    config="params:process_metadata.stage1",
                ),
                outputs=[
                    "stage1_kept",
                    "stage1_dropped",
                    "stage1_drop_reason_counts",
                ],
                name="stage1_runlevel_curation",
            ),
            node(
                func=summarize_runinfo_markdown,
                inputs=dict(
                    runinfo="stage1_kept",
                    params="params:process_metadata.stage1_summary",
                ),
                outputs="stage1_summary_md",
                name="stage1_summary_md",
            ),
            node(
                func=stage3_bioproject_curation,
                inputs=dict(
                    stage1_kept="stage1_kept",
                    params="params:process_metadata.stage3",
                ),
                outputs=[
                    "bioproject_decisions",
                    "stage3_auto_kept",
                    "stage3_dropped",
                    "stage3_manual",
                ],
                name="stage3_bioproject_curation",
            ),

            # 2a) TEMPLATE always generated
            node(
                func=make_manual_bioproject_template,
                inputs="stage3_manual",
                outputs="bioprojects_manual_decisions_template",
                name="stage3_manual_decisions_template",
            ),

            # 2b) Gate: STOP with nice message until FILLED exists; then load FILLED into DF
            node(
                func=require_and_load_manual_bioprojects_filled,
                inputs=dict(
                    manual_template="bioprojects_manual_decisions_template",
                    manual_filled_path="params:process_metadata.stage3.manual_filled_path",
                    manual_template_path="params:process_metadata.stage3.manual_template_path",
                ),
                outputs="bioprojects_manual_decisions_filled_df",
                name="stage3_require_manual_filled",
            ),

            # 2c) Merge decisions (auto_keep + manual OK)
            node(
                func=apply_manual_and_merge_kept,
                inputs=dict(
                    stage1_kept_full="stage1_kept",
                    stage3_auto_kept="stage3_auto_kept",
                    manual_decisions="bioprojects_manual_decisions_filled_df",  # ✅ now DF from gate node
                ),
                outputs=[
                    "stage3_final_kept",
                    "bioprojects_final_keep",
                ],
                name="stage3_apply_manual_and_merge",
            ),

            node(
                func=summarize_runinfo_markdown,
                inputs=dict(
                    runinfo="stage3_final_kept",
                    params="params:process_metadata.final_summary",
                ),
                outputs="final_summary_md",
                name="final_summary_md",
            ),

            node(
                func=split_strategy,
                inputs="stage3_final_kept",
                outputs=[
                    "final_wgs",
                    "final_amplicon",
                    "final_others",
                    "final_strategy_counts",
                ],
                name="stage5_split_strategy",
            ),
        ]
    )