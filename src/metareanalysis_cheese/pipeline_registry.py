from kedro.pipeline import Pipeline

from metareanalysis_cheese.pipelines.fetch_runinfo import create_pipeline as fetch_runinfo
from metareanalysis_cheese.pipelines.process_metadata import create_pipeline as process_metadata
from metareanalysis_cheese.pipelines.summarize_runinfo import create_pipeline as summarize_runinfo
from metareanalysis_cheese.pipelines.by_continent import create_pipeline as by_continent

from metareanalysis_cheese.pipelines.sra_fetch_unpack import create_pipeline as sra_fetch_unpack
from metareanalysis_cheese.pipelines.sra_analyze_wgs import create_pipeline as sra_analyze_wgs
from metareanalysis_cheese.pipelines.sra_analyze_amplicon import create_pipeline as sra_analyze_amplicon

from metareanalysis_cheese.pipelines.sra_orchestrator import create_pipeline as sra_orchestrator


def register_pipelines() -> dict[str, Pipeline]:
    p_fetch = fetch_runinfo()
    p_proc = process_metadata()
    p_sum = summarize_runinfo()
    p_cont = by_continent()

    p_sra_fetch = sra_fetch_unpack()
    p_wgs = sra_analyze_wgs()
    p_amp = sra_analyze_amplicon()

    p_orch = sra_orchestrator()


    p_default = p_fetch + p_proc + p_sum + p_cont + p_sra_fetch + p_wgs + p_amp + p_orch

    return {
        "fetch_runinfo": p_fetch,
        "process_metadata": p_proc,
        "summarize_runinfo": p_sum,
        "by_continent": p_cont,
        "sra_fetch_unpack": p_sra_fetch,
        "sra_analyze_wgs": p_wgs,
        "sra_analyze_amplicon": p_amp,
        "sra_orchestrator": p_orch,
        "__default__": p_default,
    }