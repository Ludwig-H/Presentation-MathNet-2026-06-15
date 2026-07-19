# Guidelines for Agents & GitHub-Flavored Markdown Math Standards

Welcome! This repository, **Presentation-MathNet-2026-06-15**, contains presentation materials (LaTeX Beamer) and documents related to community detection on signed graphs, mathematical networks, and thesis manuscripts.

To ensure consistency, readability, and correct rendering across GitHub, all agents and contributors must follow the guidelines detailed below.

---

## 1. General Repository Guidelines

*   **Compilation Artifacts:** LaTeX compilation produces many auxiliary files (e.g., `.aux`, `.log`, `.nav`, `.out`, `.snm`, `.toc`, `.fls`, `.fdb_latexmk`). Never commit these files. Ensure your `.gitignore` is active and correct.
*   **Asset Management:** Store image assets and figures in designated directories (e.g., `theme/imgs/` or `imgs/`). Use relative paths with correct casing.
*   **Documentation:** Maintain all documentation in Markdown format using the standards defined in the next section.
*   **Git Branch Policy:** Work and push directly on `main`. Do not create auxiliary, topic, `agent/*`, or other parasite branches unless the user explicitly requests one. Before every push, verify that the current branch is `main` and target `origin/main` explicitly.

---

## 2. GitHub-Flavored Markdown (GFM) Math Equation Standards

GitHub uses **MathJax** to render LaTeX mathematical expressions. The canonical reference is GitHub's official guide, [Writing mathematical expressions](https://docs.github.com/en/get-started/writing-on-github/working-with-advanced-formatting/writing-mathematical-expressions).

> [!IMPORTANT]
> **Correct rendering on GitHub is a release requirement, not a stylistic preference.** Follow GitHub's delimiters and Markdown nesting rules scrupulously. Never invent delimiter variants, rely on a local renderer, or assume that a syntactically balanced formula will render on GitHub. Work containing raw delimiters, formulas displayed as code, broken tables, or unrendered LaTeX is not complete.

### A. Inline Math Expressions
To render math inline with your text, surround the LaTeX expression with dollar signs (`$`).

*   **Spacing Rule:** Do **not** add spaces immediately after the opening `$` or before the closing `$`.
    *   *Correct:* `$\theta = \alpha + \beta$` renders as $\theta = \alpha + \beta$.
    *   *Incorrect:* `$ \theta = \alpha + \beta $` (may fail to render).
*   **Escaping Conflicts:** If the equation contains Markdown control characters (such as underscores `_` or asterisks `*` or square brackets `[` / `]`), the Markdown parser may process them before MathJax sees them.
    *   *Solution (Backtick Syntax):* Use GitHub's exact dollar-backtick form: the opening delimiter is a dollar sign immediately followed by a backtick, and the closing delimiter is a backtick immediately followed by a dollar sign. Do not insert spaces between either pair of delimiter characters. Copy the canonical form from GitHub's official guide when in doubt.

### B. Block Math Expressions
To render math as a centered block on its own line, use double-dollar delimiters or a fenced code block whose info string is `math`.

*   **Double Dollar Sign (`$$`) Spacing:** The opening `$$` and closing `$$` must be on their own lines, with no other content on those lines.
    *   *Example:*
        ```markdown
        $$
        \mathbf{P}(A \mid B) = \frac{\mathbf{P}(B \mid A)\mathbf{P}(A)}{\mathbf{P}(B)}
        $$
        ```
*   **Math Code Block Syntax (Mandatory for complex/multiline math):** You **must** use fenced `math` blocks for any multi-line equations (equations containing newlines, `\\`, `\substack`, or environments like `align`, `matrix`, etc.). This treats the entire block as raw LaTeX, preventing any Markdown parsing conflicts.
    *   *Example:*
        ````markdown
        ```math
        \left( \sum_{k=1}^n a_k b_k \right)^2 \leq \left( \sum_{k=1}^n a_k^2 \right) \left( \sum_{k=1}^n b_k^2 \right)
        ```
        ````

*   **Block-layout rules:** Leave a blank line before and after every display expression. Never combine a fenced `math` block with double-dollar delimiters, place a display block inside a Markdown table cell, or indent it accidentally as an ordinary code block. Move long formulas out of tables and refer to them from the relevant cell.
*   **Context rules:** Do not put TeX delimiters inside Mermaid labels, URLs, HTML attributes, or ordinary fenced code blocks and expect GitHub to render them as mathematics. In lists and callouts, preserve the blank lines and indentation required for GitHub to recognize the math block.

### C. Escaping Dollar Signs on the Same Line
If you need to use a literal dollar sign in the same line as a math expression:
*   **Inside Math:** Use a backslash: `\$`.
    *   Use the official GitHub dollar-sign example as the canonical reference; do not add spaces to dollar-backtick delimiters.
*   **Outside Math:** Wrap the literal dollar sign in `<span>$</span>`.
    *   *Example:* `To split <span>$</span>100 in half, we calculate $100/2$` renders as "To split <span>$</span>100 in half, we calculate $100/2$".

### D. Repository-Specific Command Compatibility
*   **Avoid `\operatorname`**: GitHub uses MathJax, not KaTeX. Nevertheless, this repository forbids `\operatorname{...}` (e.g. `\operatorname{Exp}`, `\operatorname{ov}`) because its validator and secondary rendering paths do not handle it consistently.
    *   *Solution:* Always use `\mathrm{...}` or `\text{...}` instead (e.g., `\mathrm{Exp}`, `\mathrm{ov}`, `\mathrm{RG}`, `\mathrm{tr}`, `\mathrm{diag}`).

---

## 3. LaTeX Validity Checklist

When adding mathematical equations, ensure:
1.  All parentheses `()`, brackets `[]`, and braces `{}` are properly balanced.
2.  Special characters like `%`, `&`, `_`, `#` are properly escaped when not used in their LaTeX command contexts.
3.  Any newline in multiline equations (like inside `align` or `substack` environments) uses `\\` properly and the equation is wrapped in a fenced `math` block.
4.  No `\operatorname` is used in markdown equations.
5.  Inline dollar and dollar-backtick delimiters use one exact GitHub-supported form, are adjacent to their delimiter characters, and close on the same line.
6.  Display blocks are separated from prose by blank lines, are not embedded in table cells, and are not accidentally parsed as ordinary code.
7.  Every edited formula is inspected in context for collisions with lists, blockquotes, HTML, links, code fences, and Mermaid diagrams.

---

## 4. Automatic and Visual Verification

A Python script is available to automatically validate math syntax in all Markdown files within this workspace:
*   File: [.agents/check_math.py](.agents/check_math.py)
*   **How to run**: `python3 .agents/check_math.py` or `./.agents/check_math.py`

**Mandatory Rule for Agents:** Whenever you edit or add Markdown containing formulas, you must:

1. run `python3 .agents/check_math.py`;
2. run `python3 .agents/check_markdown_links.py`;
3. inspect the changed Markdown on GitHub after pushing, paying particular attention to dollar-backtick expressions, blocks near lists or callouts, formulas near tables, and any raw delimiter, backtick, or TeX text;
4. correct rendering defects immediately on `main`.

The scripts are necessary preflight checks, but they are not proof that GitHub rendered the document correctly. The task is complete only when both the automated checks and the actual GitHub rendering are clean.
