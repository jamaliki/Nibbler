"""Inspect the frozen benchmark plan before executable adapters are implemented."""

from __future__ import annotations

import argparse

from .baselines import BASELINES
from .workloads import WORKLOADS


def build_parser() -> argparse.ArgumentParser:
    """Create the benchmark-plan command-line parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--list",
        action="store_true",
        dest="list_plan",
        help="list workloads and baseline engines",
    )
    return parser


def render_plan() -> str:
    """Render the stable human-readable benchmark plan."""
    lines = ["Workloads:"]
    lines.extend(f"  {workload.name}: {workload.description}" for workload in WORKLOADS)
    lines.append("Baselines:")
    lines.extend(
        f"  {baseline.name} ({baseline.distribution}): {baseline.purpose}"
        for baseline in BASELINES
    )
    return "\n".join(lines)


def main() -> int:
    """List the plan or reject attempts to report premature timings."""
    arguments = build_parser().parse_args()
    if arguments.list_plan:
        print(render_plan())
        return 0
    build_parser().error(
        "timed adapters are added with the shared parser; use --list during Phase 0"
    )


if __name__ == "__main__":
    raise SystemExit(main())
