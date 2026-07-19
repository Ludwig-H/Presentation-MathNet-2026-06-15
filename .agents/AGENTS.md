# Workspace Rules: Presentation-MathNet-2026-06-15

This workspace configuration enforces specific rules for files and math notations.

## Math Notation Standards in Markdown (.md)
To ensure equations render perfectly on GitHub via MathJax:

1. **Inline Math**: Always use `$equation$` with no whitespace immediately following the opening `$` or preceding the closing `$`.
2. **Conflict Prevention**: For equations containing Markdown control characters, use GitHub's exact dollar-backtick syntax: the opening dollar sign and backtick must be adjacent, as must the closing backtick and dollar sign. Never insert spaces inside either delimiter.
3. **Block Math**: Place double-dollar delimiters on their own lines. For complex multi-line math (e.g. using `\\`, `\substack`, or environments), you **must** use a fenced `math` block to prevent GFM parser interference.
4. **Repository Command Compatibility**: GitHub uses MathJax. This repository nevertheless forbids `\operatorname{...}` because its validator and secondary rendering paths do not handle it consistently. Always use `\mathrm{...}` or `\text{...}` instead.
5. **Dollar Signs**: Escape literal dollar signs inside math as `\$` and wrap them in `<span>$</span>` outside math on the same line.
6. **No Nested Math Fences**: Every opening and closing fence of a `math` block must begin at column 1. Never indent or nest it in a list, blockquote, or callout: GitHub otherwise renders it as ordinary code rather than MathJax. Restructure the surrounding Markdown instead.

Refer to [agents.md](../agents.md) for full examples and guides.

## Automatic Verification
To verify that all math formulas in the workspace are valid, run:
```bash
python3 .agents/check_math.py
python3 .agents/check_markdown_links.py
```
The first script checks balanced braces, correct delimiter escaping, and
invalid LaTeX commands. The second checks every relative Markdown link.
These checks are necessary but not sufficient: inspect the actual
GitHub-rendered Markdown after every push and correct any raw delimiter,
unrendered formula, broken table, or misparsed code block immediately on
`main`.
