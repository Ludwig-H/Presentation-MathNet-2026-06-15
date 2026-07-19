#!/usr/bin/env python3
"""Validate relative links in Markdown files without external dependencies."""

from __future__ import annotations

import re
import sys
from pathlib import Path
from urllib.parse import unquote


ROOT = Path(__file__).resolve().parent.parent
LINK_PATTERN = re.compile(r"(?<!!)\[[^\]]+\]\(([^)]+)\)")
SKIPPED_SCHEMES = (
    "http://",
    "https://",
    "mailto:",
    "file://",
    "data:",
)


def markdown_files() -> list[Path]:
    return sorted(
        path
        for path in ROOT.rglob("*.md")
        if ".git" not in path.parts
    )


def local_target(raw_target: str) -> str | None:
    target = raw_target.strip()
    if not target or target.startswith("#") or target.startswith(SKIPPED_SCHEMES):
        return None
    target = target.split("#", 1)[0]
    if not target or target.startswith("/"):
        return None
    return unquote(target)


def main() -> int:
    failures: list[str] = []
    checked = 0
    for markdown in markdown_files():
        content = markdown.read_text(encoding="utf-8")
        for match in LINK_PATTERN.finditer(content):
            if content.count("```", 0, match.start()) % 2:
                continue
            target = local_target(match.group(1))
            if target is None:
                continue
            checked += 1
            resolved = markdown.parent / target
            if not resolved.exists():
                line = content.count("\n", 0, match.start()) + 1
                relative = markdown.relative_to(ROOT)
                failures.append(f"{relative}:{line}: missing {target}")

    if failures:
        print("Broken relative Markdown links:")
        for failure in failures:
            print(f"- {failure}")
        return 1

    print(f"Markdown links valid: checked {checked} relative links.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
