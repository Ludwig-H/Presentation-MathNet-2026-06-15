# Workspace Rules: Presentation-MathNet-2026-06-15

This workspace configuration enforces specific rules for files and math notations.

## Math Notation Standards in Markdown (.md)
To ensure equations render perfectly on GitHub via MathJax:

1. **Inline Math**: Always use `$equation$` with no whitespace immediately following the opening `$` or preceding the closing `$`.
2. **Conflict Prevention**: For equations containing characters like `_`, `*`, `[`, or `]`, use the backtick syntax `$ `equation` $` to prevent Markdown parsing interference.
3. **Block Math**: Place double dollar signs `$$` on their own lines. For complex multi-line math (e.g. using `\\`, `\substack`, or environments), you **must** use the ` ```math ` code block syntax to prevent GFM parser interference.
4. **Unsupported LaTeX Commands**: Do **not** use `\operatorname{...}` (e.g. `\operatorname{Exp}`, `\operatorname{ov}`) as it is often unsupported or causes rendering errors in GitHub/KaTeX. Always use `\mathrm{...}` or `\text{...}` (e.g. `\mathrm{Exp}`, `\mathrm{ov}`, `\mathrm{RG}`, `\mathrm{tr}`, `\mathrm{diag}`) instead.
5. **Dollar Signs**: Escape literal dollar signs inside math as `\$` and wrap them in `<span>$</span>` outside math on the same line.

Refer to [agents.md](file:///workspaces/Presentation-MathNet-2026-06-15/agents.md) for full examples and guides.

## Automatic Verification
To verify that all math formulas in the workspace are valid, run:
```bash
python3 .agents/check_math.py
python3 .agents/check_markdown_links.py
```
The first script checks balanced braces, correct delimiter escaping, and
invalid LaTeX commands. The second checks every relative Markdown link.

