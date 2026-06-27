"""
centro_peak_agent_v1.1.py

Single-agent LangGraph workflow for identifying repeat-period peaks
from ULTRA tandem repeat data.

This agent does NOT run ULTRA yet.

It expects an existing ULTRA TSV file and looks for prominent repeat
periods above a user-defined minimum period, defaulting to Period > 60.

Main workflow:

1. Validate input file and parameters.
2. Read ULTRA tandem repeat table.
3. Filter repeats by Period > min_period.
4. Count repeats by period.
5. Detect prominent period peaks.
6. Group adjacent candidate periods, such as 91 and 92 bp.
7. Optionally ask Claude to write a scientific interpretation.
8. Write CSV, JSON, and Markdown report outputs.

Example usage:

python centro_peak_agent_v1.1.py \
    --ultra-tsv ultra.Wm82.gnm6.S97D.tsv \
    --species-name "Glycine max Wm82.gnm6" \
    --output-dir results/glycine_peak_agent \
    --min-period 60 \
    --max-period 300

Required packages:

pip install -U langgraph langchain langchain-anthropic python-dotenv pandas numpy scipy typing-extensions

# Future provider support, when you are ready to enable it:
# pip install -U langchain-openai langchain-google-genai

Optional .env file:

AI_PROVIDER=anthropic
AI_API_KEY=your_key_here
AI_MODEL=claude-sonnet-4-6

If AI_API_KEY is not found, the script still runs and writes a
basic deterministic report without the LLM review.

For now, only AI_PROVIDER=anthropic is supported. OpenAI and Gemini imports
are included below as commented future placeholders.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from dotenv import load_dotenv
from scipy.signal import find_peaks
from typing_extensions import TypedDict

from langchain_anthropic import ChatAnthropic
# Future OpenAI support: uncomment this import when you add OpenAI provider logic.
# from langchain_openai import ChatOpenAI
# Future Gemini support: uncomment this import when you add Gemini provider logic.
# from langchain_google_genai import ChatGoogleGenerativeAI
from langchain_core.messages import HumanMessage, SystemMessage
from langgraph.graph import END, START, StateGraph


load_dotenv()


class PeakAgentState(TypedDict, total=False):
    """
    Shared state passed between LangGraph nodes.

    The early fields are user inputs and analysis parameters.
    The later fields are outputs produced by each workflow step.
    """

    ultra_tsv: str
    species_name: str
    output_dir: str

    min_period: int
    max_period: int | None
    min_count: int
    neighborhood: int
    group_gap: int
    neighbor_support_fraction: float
    min_prominence: float
    min_prominence_fraction: float
    min_prominence_ratio: float
    top_n: int

    input_summary: dict[str, Any]
    period_distribution: list[dict[str, Any]]
    candidate_period_table: list[dict[str, Any]]
    candidate_groups: list[dict[str, Any]]

    scientific_review: str
    llm_used: bool
    output_files: dict[str, str]


def to_jsonable(obj: Any) -> Any:
    """
    Convert numpy/pandas objects into values that can be safely written as JSON.
    """

    if obj is None:
        return None

    if isinstance(obj, dict):
        return {k: to_jsonable(v) for k, v in obj.items()}

    if isinstance(obj, list):
        return [to_jsonable(v) for v in obj]

    if isinstance(obj, tuple):
        return [to_jsonable(v) for v in obj]

    if isinstance(obj, np.ndarray):
        return obj.tolist()

    if isinstance(obj, np.integer):
        return int(obj)

    if isinstance(obj, np.floating):
        value = float(obj)
        if np.isnan(value) or np.isinf(value):
            return None
        return value

    if isinstance(obj, float):
        if np.isnan(obj) or np.isinf(obj):
            return None
        return obj

    if isinstance(obj, pd.Timestamp):
        return obj.isoformat()

    try:
        if pd.isna(obj):
            return None
    except (TypeError, ValueError):
        pass

    return obj


def normalize_column_name(name: str) -> str:
    """
    Normalize ULTRA column names for easier matching.

    Example:
    '#Subrepeats' becomes 'subrepeats'
    'Period' becomes 'period'
    """

    return (
        str(name)
        .strip()
        .lstrip("#")
        .lower()
        .replace(" ", "_")
        .replace("-", "_")
    )


def read_ultra_table(path: str | Path) -> pd.DataFrame:
    """
    Read ULTRA output into a pandas DataFrame.

    Only the Period column is required for this first workflow.
    Other ULTRA columns, such as SeqID, Start, End, Score, Consensus,
    and Subrepeats, are useful later but are not required here.
    """

    path = Path(path)

    try:
        df = pd.read_csv(path, sep="\t", engine="python")
    except Exception:
        df = pd.read_csv(path, sep=r"\s+", engine="python")

    if df.shape[1] == 1:
        df = pd.read_csv(path, sep=r"\s+", engine="python")

    original_columns = list(df.columns)
    normalized_columns = {col: normalize_column_name(col) for col in original_columns}
    df = df.rename(columns=normalized_columns)

    # Accept a few common alternative names for the Period column.
    period_aliases = {
        "period",
        "period_bp",
        "repeat_period",
        "motif_period",
        "monomer_size",
        "monomer_size_bp",
    }

    detected_period_columns = [col for col in df.columns if col in period_aliases]

    if "period" not in df.columns:
        if detected_period_columns:
            df = df.rename(columns={detected_period_columns[0]: "period"})
        else:
            raise ValueError(
                "Could not find a required Period column.\n\n"
                "For this first workflow, the only required column is Period.\n"
                "The #Subrepeats column is optional and is not analyzed yet.\n\n"
                f"Original columns detected: {original_columns}\n"
                f"Normalized columns detected: {list(df.columns)}"
            )

    optional_columns = {
        "seqid",
        "start",
        "end",
        "score",
        "consensus",
        "subrepeats",
    }

    missing_optional_columns = sorted(optional_columns - set(df.columns))

    if missing_optional_columns:
        print(
            "Warning: optional columns are missing and will be ignored for now: "
            + ", ".join(missing_optional_columns)
        )

    df["period"] = pd.to_numeric(df["period"], errors="coerce")
    df = df.dropna(subset=["period"]).copy()

    if df.empty:
        raise ValueError(
            "The Period column was found, but no valid numeric Period values "
            "remained after parsing."
        )

    df["period_int"] = df["period"].round().astype(int)

    return df


def local_baseline(
    distribution: pd.DataFrame,
    period: int,
    neighborhood: int,
    exclude_radius: int = 1,
) -> float:
    """
    Calculate the local median baseline around a period.

    The period itself and nearby periods within exclude_radius are excluded.
    This helps the local baseline represent nearby background instead of the
    peak itself.
    """

    mask = (
        (distribution["period"] >= period - neighborhood)
        & (distribution["period"] <= period + neighborhood)
        & (
            (distribution["period"] < period - exclude_radius)
            | (distribution["period"] > period + exclude_radius)
        )
    )

    nearby_counts = distribution.loc[mask, "count"].to_numpy()

    if len(nearby_counts) == 0:
        return 0.0

    return float(np.median(nearby_counts))


def group_periods(periods: list[int], group_gap: int) -> list[list[int]]:
    """
    Group candidate periods that are close to each other.

    Example:
    periods = [91, 92, 120]
    group_gap = 2

    Returns:
    [[91, 92], [120]]
    """

    if not periods:
        return []

    sorted_periods = sorted(set(periods))
    groups = [[sorted_periods[0]]]

    for period in sorted_periods[1:]:
        if period - groups[-1][-1] <= group_gap:
            groups[-1].append(period)
        else:
            groups.append([period])

    return groups


def analyze_period_peaks(
    df: pd.DataFrame,
    min_period: int = 60,
    max_period: int | None = None,
    min_count: int = 10,
    neighborhood: int = 10,
    group_gap: int = 2,
    neighbor_support_fraction: float = 0.35,
    min_prominence: float = 5.0,
    min_prominence_fraction: float = 0.02,
    min_prominence_ratio: float = 3.0,
    top_n: int = 10,
) -> dict[str, Any]:
    """
    Deterministically identify prominent repeat-period peaks.

    Core idea:
    1. Filter to Period > min_period.
    2. Count repeat records per period.
    3. Detect local peaks.
    4. Rank candidate periods by count and local prominence.
    5. Group nearby candidate periods, such as 91 and 92 bp.
    """

    filtered = df[df["period_int"] > min_period].copy()

    if max_period is not None:
        filtered = filtered[filtered["period_int"] <= max_period].copy()

    if filtered.empty:
        raise ValueError(
            f"No records remained after filtering Period > {min_period}."
        )

    counts_by_period = (
        filtered["period_int"]
        .value_counts()
        .sort_index()
        .rename_axis("period")
        .reset_index(name="count")
    )

    full_range = pd.RangeIndex(
        int(counts_by_period["period"].min()),
        int(counts_by_period["period"].max()) + 1,
    )

    distribution = (
        counts_by_period
        .set_index("period")
        .reindex(full_range, fill_value=0)
        .rename_axis("period")
        .reset_index()
    )

    total_filtered = int(distribution["count"].sum())
    max_count = float(distribution["count"].max())

    distribution["percent_of_filtered"] = (
        distribution["count"] / total_filtered * 100
    )

    distribution["local_baseline"] = distribution["period"].apply(
        lambda p: local_baseline(
            distribution=distribution,
            period=int(p),
            neighborhood=neighborhood,
            exclude_radius=group_gap,
        )
    )

    distribution["prominence_ratio"] = distribution.apply(
        lambda row: float(row["count"]) / max(float(row["local_baseline"]), 1.0),
        axis=1,
    )

    y = distribution["count"].to_numpy(dtype=float)
    period_values = distribution["period"].to_numpy(dtype=int)

    effective_min_prominence = max(
        float(min_prominence),
        float(max_count) * float(min_prominence_fraction),
    )

    peak_indices, peak_properties = find_peaks(
        y,
        prominence=effective_min_prominence,
    )

    seed_periods: set[int] = set()
    peak_prominence_lookup: dict[int, float] = {}

    for idx, prominence in zip(
        peak_indices,
        peak_properties.get("prominences", []),
    ):
        period = int(period_values[idx])
        row = distribution.loc[distribution["period"] == period].iloc[0]

        peak_prominence_lookup[period] = float(prominence)

        if (
            int(row["count"]) >= min_count
            and float(row["prominence_ratio"]) >= min_prominence_ratio
        ):
            seed_periods.add(period)

    # Fallback: if strict local peak detection finds nothing, use top periods.
    if not seed_periods:
        fallback = (
            distribution[distribution["count"] >= min_count]
            .sort_values(["count", "prominence_ratio"], ascending=False)
            .head(top_n)
        )
        seed_periods = set(int(p) for p in fallback["period"].tolist())

    count_lookup = {
        int(row["period"]): int(row["count"])
        for _, row in distribution.iterrows()
    }

    candidate_periods: set[int] = set(seed_periods)

    # Add nearby supporting periods so adjacent peaks like 91/92 can group.
    for seed in seed_periods:
        seed_count = count_lookup.get(seed, 0)

        for period in range(seed - group_gap, seed + group_gap + 1):
            if period not in count_lookup:
                continue

            support_threshold = max(
                min_count,
                int(seed_count * neighbor_support_fraction),
            )

            if count_lookup[period] >= support_threshold:
                candidate_periods.add(period)

    groups = group_periods(list(candidate_periods), group_gap=group_gap)

    group_records: list[dict[str, Any]] = []
    period_to_group: dict[int, str] = {}

    for i, group in enumerate(groups, start=1):
        group_id = f"group_{i:03d}"
        group_df = distribution[distribution["period"].isin(group)].copy()

        peak_row = group_df.sort_values("count", ascending=False).iloc[0]
        peak_period = int(peak_row["period"])
        peak_count = int(peak_row["count"])
        total_group_count = int(group_df["count"].sum())

        left = min(group) - neighborhood
        right = max(group) + neighborhood

        baseline_mask = (
            (distribution["period"] >= left)
            & (distribution["period"] <= right)
            & (~distribution["period"].isin(group))
        )

        baseline_values = distribution.loc[baseline_mask, "count"].to_numpy()
        group_local_baseline = (
            float(np.median(baseline_values)) if len(baseline_values) else 0.0
        )

        group_prominence_ratio = peak_count / max(group_local_baseline, 1.0)

        # Simple ranking score: abundance multiplied by local contrast.
        ranking_score = float(
            np.log10(total_group_count + 1)
            * np.log2(group_prominence_ratio + 1)
        )

        group_record = {
            "group_id": group_id,
            "periods": [int(p) for p in group],
            "period_range": (
                f"{min(group)}"
                if min(group) == max(group)
                else f"{min(group)}-{max(group)}"
            ),
            "peak_period": peak_period,
            "peak_count": peak_count,
            "total_group_count": total_group_count,
            "percent_of_filtered": round(
                total_group_count / total_filtered * 100,
                4,
            ),
            "local_baseline": round(group_local_baseline, 4),
            "prominence_ratio": round(group_prominence_ratio, 4),
            "ranking_score": round(ranking_score, 4),
            "seed_periods": sorted(int(p) for p in seed_periods if p in group),
            "interpretation_hint": (
                "Strong candidate"
                if group_prominence_ratio >= 5 and total_group_count >= min_count
                else "Candidate; inspect further"
            ),
        }

        group_records.append(group_record)

        for period in group:
            period_to_group[int(period)] = group_id

    group_records = sorted(
        group_records,
        key=lambda x: (x["ranking_score"], x["total_group_count"]),
        reverse=True,
    )

    distribution["is_seed_peak"] = distribution["period"].apply(
        lambda p: int(p) in seed_periods
    )

    distribution["is_candidate"] = distribution["period"].apply(
        lambda p: int(p) in candidate_periods
    )

    distribution["candidate_group_id"] = distribution["period"].apply(
        lambda p: period_to_group.get(int(p), "")
    )

    distribution["peak_prominence"] = distribution["period"].apply(
        lambda p: peak_prominence_lookup.get(int(p), 0.0)
    )

    candidate_table = (
        distribution[distribution["is_candidate"]]
        .sort_values(["count", "prominence_ratio"], ascending=False)
        .copy()
    )

    input_summary = {
        "total_ultra_records": int(len(df)),
        "records_after_period_filter": int(len(filtered)),
        "min_period": int(min_period),
        "max_period": None if max_period is None else int(max_period),
        "min_count": int(min_count),
        "neighborhood": int(neighborhood),
        "group_gap": int(group_gap),
        "neighbor_support_fraction": float(neighbor_support_fraction),
        "min_prominence": float(min_prominence),
        "min_prominence_fraction": float(min_prominence_fraction),
        "effective_min_prominence": float(effective_min_prominence),
        "min_prominence_ratio": float(min_prominence_ratio),
        "number_of_seed_peaks": int(len(seed_periods)),
        "number_of_candidate_groups": int(len(group_records)),
    }

    return {
        "input_summary": input_summary,
        "period_distribution": distribution.to_dict(orient="records"),
        "candidate_period_table": candidate_table.to_dict(orient="records"),
        "candidate_groups": group_records[:top_n],
    }


def validate_input_node(state: PeakAgentState) -> dict[str, Any]:
    """
    LangGraph node 1.

    Validate input path and set default parameters.
    """

    ultra_tsv = Path(state["ultra_tsv"])

    if not ultra_tsv.exists():
        raise FileNotFoundError(f"Could not find input file: {ultra_tsv}")

    max_period = state.get("max_period", None)

    return {
        "ultra_tsv": str(ultra_tsv),
        "species_name": state.get("species_name", "Unknown species"),
        "output_dir": state.get("output_dir", "centro_peak_results"),
        "min_period": int(state.get("min_period", 60)),
        "max_period": None if max_period is None else int(max_period),
        "min_count": int(state.get("min_count", 10)),
        "neighborhood": int(state.get("neighborhood", 10)),
        "group_gap": int(state.get("group_gap", 2)),
        "neighbor_support_fraction": float(
            state.get("neighbor_support_fraction", 0.35)
        ),
        "min_prominence": float(state.get("min_prominence", 5.0)),
        "min_prominence_fraction": float(
            state.get("min_prominence_fraction", 0.02)
        ),
        "min_prominence_ratio": float(state.get("min_prominence_ratio", 3.0)),
        "top_n": int(state.get("top_n", 10)),
    }


def peak_analysis_node(state: PeakAgentState) -> dict[str, Any]:
    """
    LangGraph node 2.

    Read the ULTRA table and run deterministic peak analysis.
    """

    df = read_ultra_table(state["ultra_tsv"])

    results = analyze_period_peaks(
        df=df,
        min_period=state["min_period"],
        max_period=state.get("max_period"),
        min_count=state["min_count"],
        neighborhood=state["neighborhood"],
        group_gap=state["group_gap"],
        neighbor_support_fraction=state["neighbor_support_fraction"],
        min_prominence=state["min_prominence"],
        min_prominence_fraction=state["min_prominence_fraction"],
        min_prominence_ratio=state["min_prominence_ratio"],
        top_n=state["top_n"],
    )

    return to_jsonable(results)


def make_basic_review(state: PeakAgentState) -> str:
    """
    Create a deterministic report if no Anthropic API key is available
    or if the LLM call fails.
    """

    groups = state.get("candidate_groups", [])

    if not groups:
        return (
            "No candidate repeat-period groups were detected with the current "
            "parameters. Consider lowering min_count, min_prominence_ratio, or "
            "inspecting the period distribution manually."
        )

    top_group = groups[0]

    return f"""
# Scientific Review

The strongest candidate repeat-period group is **{top_group["period_range"]} bp**.

This group has:

- Peak period: {top_group["peak_period"]} bp
- Peak count: {top_group["peak_count"]}
- Total group count: {top_group["total_group_count"]}
- Percent of filtered repeats: {top_group["percent_of_filtered"]}%
- Local prominence ratio: {top_group["prominence_ratio"]}

This result supports advancing this repeat-period group to the next analysis step:
chromosome-level localization and array-size inspection.

Important limitation: this step only identifies prominent repeat periods. It does
not by itself prove centromeric identity. The next evidence layer should test
whether the selected repeats localize to discrete chromosomal hotspots.
""".strip()


def create_llm(provider: str, model_name: str):
    """
    Create the chat model used for the optional scientific review step.

    This helper makes the rest of the workflow provider-flexible without
    pretending every AI provider uses the same LangChain class.

    Currently supported:
    - anthropic: uses ChatAnthropic

    Future providers are intentionally guarded for now. When you are ready to
    use them, uncomment the relevant imports near the top of this file and add
    provider-specific logic below.
    """

    provider = provider.lower().strip()

    if provider == "anthropic":
        return ChatAnthropic(
            model=model_name,
            temperature=0,
            max_tokens=1400,
        )

    if provider == "openai":
        raise ValueError(
            "AI_PROVIDER='openai' was requested, but OpenAI support is not "
            "enabled yet. To enable it later, uncomment the ChatOpenAI import "
            "near the top of this file and add ChatOpenAI logic inside "
            "create_llm()."
        )

    if provider in {"gemini", "google", "google_genai"}:
        raise ValueError(
            f"AI_PROVIDER='{provider}' was requested, but Gemini support is not "
            "enabled yet. To enable it later, uncomment the "
            "ChatGoogleGenerativeAI import near the top of this file and add "
            "Gemini logic inside create_llm()."
        )

    raise ValueError(
        f"Unsupported AI_PROVIDER: {provider}. "
        "This version currently supports only 'anthropic'."
    )


def scientific_review_node(state: PeakAgentState) -> dict[str, Any]:
    """
    LangGraph node 3.

    Ask the configured AI model to interpret the deterministic peak results.

    The LLM is not allowed to invent candidates. It only reviews the
    candidate groups already produced by the Python analysis.
    """

    provider = os.getenv("AI_PROVIDER", "anthropic").lower().strip()
    api_key = os.getenv("AI_API_KEY")
    model_name = os.getenv("AI_MODEL", "claude-sonnet-4-6")

    if not api_key:
        return {
            "llm_used": False,
            "scientific_review": make_basic_review(state),
        }

    try:
        llm = create_llm(provider=provider, model_name=model_name)
    except ValueError as error:
        fallback = make_basic_review(state)
        fallback += (
            "\n\nLLM review was not completed because of this configuration error:\n"
            f"`{error}`"
        )
        return {
            "llm_used": False,
            "scientific_review": fallback,
        }

    system_message = SystemMessage(
        content=(
            "You are a scientific repeat-analysis assistant for plant "
            "comparative genomics. You do not invent candidate repeats. "
            "You only interpret the deterministic peak-analysis results "
            "provided to you. Clearly distinguish data-supported conclusions "
            "from hypotheses. Recommend the next analysis step."
        )
    )

    top_candidate_rows = state.get("candidate_period_table", [])[:25]

    human_message = HumanMessage(
        content=f"""
Species:
{state.get("species_name", "Unknown species")}

Input summary:
{json.dumps(to_jsonable(state.get("input_summary", {})), indent=2)}

Ranked candidate groups:
{json.dumps(to_jsonable(state.get("candidate_groups", [])), indent=2)}

Top candidate period rows:
{json.dumps(to_jsonable(top_candidate_rows), indent=2)}

Write a concise Markdown report with these sections:

1. Summary conclusion
2. Strongest candidate repeat periods
3. Why these were selected
4. Caveats
5. Recommended next step

Remember: this step only identifies prominent repeat periods from ULTRA output.
It does not prove centromeric identity without chromosome localization evidence.
"""
    )

    try:
        response = llm.invoke([system_message, human_message])
        return {
            "llm_used": True,
            "scientific_review": response.content,
        }

    except Exception as error:
        fallback = make_basic_review(state)
        fallback += (
            "\n\nLLM review was not completed because of this error:\n"
            f"`{error}`"
        )

        return {
            "llm_used": False,
            "scientific_review": fallback,
        }


def write_outputs_node(state: PeakAgentState) -> dict[str, Any]:
    """
    LangGraph node 4.

    Write output files:
    1. Full period distribution CSV
    2. Candidate period rows CSV
    3. Candidate group summary JSON
    4. Markdown scientific report
    """

    output_dir = Path(state["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    input_stem = Path(state["ultra_tsv"]).stem

    period_distribution_path = output_dir / f"{input_stem}.period_distribution.csv"
    candidate_periods_path = output_dir / f"{input_stem}.candidate_periods.csv"
    candidate_groups_path = output_dir / f"{input_stem}.candidate_groups.json"
    report_path = output_dir / f"{input_stem}.candidate_repeat_report.md"

    pd.DataFrame(state["period_distribution"]).to_csv(
        period_distribution_path,
        index=False,
    )

    pd.DataFrame(state["candidate_period_table"]).to_csv(
        candidate_periods_path,
        index=False,
    )

    summary_json = {
        "species_name": state.get("species_name"),
        "input_file": state.get("ultra_tsv"),
        "input_summary": state.get("input_summary"),
        "candidate_groups": state.get("candidate_groups"),
        "llm_used": state.get("llm_used"),
    }

    candidate_groups_path.write_text(
        json.dumps(to_jsonable(summary_json), indent=2),
        encoding="utf-8",
    )

    report = f"""# Candidate Centromeric Repeat Period Report

## Input

- Species: {state.get("species_name")}
- ULTRA file: `{state.get("ultra_tsv")}`

## Analysis Parameters

```json
{json.dumps(to_jsonable(state.get("input_summary", {})), indent=2)}
```

## Scientific Review

{state.get("scientific_review", "")}

## Output Files

- Period distribution CSV: `{period_distribution_path}`
- Candidate periods CSV: `{candidate_periods_path}`
- Candidate groups JSON: `{candidate_groups_path}`
"""

    report_path.write_text(report, encoding="utf-8")

    return {
        "output_files": {
            "period_distribution_csv": str(period_distribution_path),
            "candidate_periods_csv": str(candidate_periods_path),
            "candidate_groups_json": str(candidate_groups_path),
            "report_markdown": str(report_path),
        }
    }


def build_graph():
    """
    Build and compile the LangGraph workflow.

    The graph is intentionally simple:

    START
      -> validate_input
      -> peak_analysis
      -> scientific_review
      -> write_outputs
      -> END
    """

    builder = StateGraph(PeakAgentState)

    builder.add_node("validate_input", validate_input_node)
    builder.add_node("peak_analysis", peak_analysis_node)
    builder.add_node("scientific_review", scientific_review_node)
    builder.add_node("write_outputs", write_outputs_node)

    builder.add_edge(START, "validate_input")
    builder.add_edge("validate_input", "peak_analysis")
    builder.add_edge("peak_analysis", "scientific_review")
    builder.add_edge("scientific_review", "write_outputs")
    builder.add_edge("write_outputs", END)

    return builder.compile()


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.
    """

    parser = argparse.ArgumentParser(
        description=(
            "Identify prominent repeat-period peaks from ULTRA tandem repeat TSV."
        )
    )

    parser.add_argument(
        "--ultra-tsv",
        required=True,
        help="Path to ULTRA tandem repeat TSV file.",
    )

    parser.add_argument(
        "--species-name",
        default="Unknown species",
        help="Species or accession name for the report.",
    )

    parser.add_argument(
        "--output-dir",
        default="centro_peak_results",
        help="Directory where output files will be written.",
    )

    parser.add_argument(
        "--min-period",
        type=int,
        default=60,
        help="Keep repeats with Period > this value. Default: 60.",
    )

    parser.add_argument(
        "--max-period",
        type=int,
        default=None,
        help="Optional upper period cutoff, for example 300.",
    )

    parser.add_argument(
        "--min-count",
        type=int,
        default=10,
        help="Minimum count for a period to be considered. Default: 10.",
    )

    parser.add_argument(
        "--neighborhood",
        type=int,
        default=10,
        help="Window size for local baseline calculation. Default: 10.",
    )

    parser.add_argument(
        "--group-gap",
        type=int,
        default=2,
        help="Group candidate periods within this distance. Default: 2.",
    )

    parser.add_argument(
        "--neighbor-support-fraction",
        type=float,
        default=0.35,
        help=(
            "Neighboring periods are grouped if their count is at least this "
            "fraction of the seed peak. Default: 0.35."
        ),
    )

    parser.add_argument(
        "--min-prominence",
        type=float,
        default=5.0,
        help="Absolute minimum scipy peak prominence. Default: 5.",
    )

    parser.add_argument(
        "--min-prominence-fraction",
        type=float,
        default=0.02,
        help=(
            "Minimum prominence as a fraction of the highest count. "
            "Default: 0.02."
        ),
    )

    parser.add_argument(
        "--min-prominence-ratio",
        type=float,
        default=3.0,
        help=(
            "Minimum count/local-baseline ratio for seed peaks. Default: 3."
        ),
    )

    parser.add_argument(
        "--top-n",
        type=int,
        default=10,
        help="Number of ranked candidate groups to keep. Default: 10.",
    )

    return parser.parse_args()


def main() -> None:
    """
    Main command-line entry point.
    """

    args = parse_args()
    app = build_graph()

    initial_state: PeakAgentState = {
        "ultra_tsv": args.ultra_tsv,
        "species_name": args.species_name,
        "output_dir": args.output_dir,
        "min_period": args.min_period,
        "max_period": args.max_period,
        "min_count": args.min_count,
        "neighborhood": args.neighborhood,
        "group_gap": args.group_gap,
        "neighbor_support_fraction": args.neighbor_support_fraction,
        "min_prominence": args.min_prominence,
        "min_prominence_fraction": args.min_prominence_fraction,
        "min_prominence_ratio": args.min_prominence_ratio,
        "top_n": args.top_n,
    }

    final_state = app.invoke(initial_state)

    print("\nAnalysis complete.\n")
    print("Top candidate groups:")

    candidate_groups = final_state.get("candidate_groups", [])

    if not candidate_groups:
        print("- No candidate groups detected.")
    else:
        for group in candidate_groups:
            print(
                f"- {group['group_id']}: {group['period_range']} bp "
                f"(peak={group['peak_period']} bp, "
                f"count={group['total_group_count']}, "
                f"prominence_ratio={group['prominence_ratio']})"
            )

    print("\nOutput files:")

    for label, path in final_state.get("output_files", {}).items():
        print(f"- {label}: {path}")


if __name__ == "__main__":
    main()
