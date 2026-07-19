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

GitHub uses **MathJax** to render LaTeX mathematical expressions. To ensure that math renders correctly and doesn't get corrupted or ignored by the Markdown parser, follow these rules:

### A. Inline Math Expressions
To render math inline with your text, surround the LaTeX expression with dollar signs (`$`).

*   **Spacing Rule:** Do **not** add spaces immediately after the opening `$` or before the closing `$`.
    *   *Correct:* `$\theta = \alpha + \beta$` renders as $\theta = \alpha + \beta$.
    *   *Incorrect:* `$ \theta = \alpha + \beta $` (may fail to render).
*   **Escaping Conflicts:** If the equation contains Markdown control characters (such as underscores `_` or asterisks `*` or square brackets `[` / `]`), the Markdown parser may process them before MathJax sees them.
    *   *Solution (Backtick Syntax):* Start the expression with `$`` ` and end it with ` ``$`.
    *   *Example:* `$ `x_i = y_i \times z_i` $` renders as $`x_i = y_i \times z_i`$.

### B. Block Math Expressions
To render math as a centered block on its own line, use double dollar signs (`$$`) or a ` ```math ` code block.

*   **Double Dollar Sign (`$$`) Spacing:** The opening `$$` and closing `$$` must be on their own lines, with no other content on those lines.
    *   *Example:*
        ```markdown
        $$
        \mathbf{P}(A \mid B) = \frac{\mathbf{P}(B \mid A)\mathbf{P}(A)}{\mathbf{P}(B)}
        $$
        ```
*   **Math Code Block Syntax (Mandatory for complex/multiline math):** You **must** use ` ```math ` code blocks for any multi-line equations (equations containing newlines, `\\`, `\substack`, or environments like `align`, `matrix`, etc.). This treats the entire block as raw LaTeX, preventing any Markdown parsing conflicts.
    *   *Example:*
        ```markdown
        ```math
        \left( \sum_{k=1}^n a_k b_k \right)^2 \leq \left( \sum_{k=1}^n a_k^2 \right) \left( \sum_{k=1}^n b_k^2 \right)
        ```
        ```

### C. Escaping Dollar Signs on the Same Line
If you need to use a literal dollar sign in the same line as a math expression:
*   **Inside Math:** Use a backslash: `\$`.
    *   *Example:* `$ `\sqrt{\$4}` $` renders as $`\sqrt{\$4}`$.
*   **Outside Math:** Wrap the literal dollar sign in `<span>$</span>`.
    *   *Example:* `To split <span>$</span>100 in half, we calculate $100/2$` renders as "To split <span>$</span>100 in half, we calculate $100/2$".

### D. Unsupported LaTeX Commands (GitHub/KaTeX Compatibility)
*   **Avoid `\operatorname`**: Do **not** use `\operatorname{...}` (e.g. `\operatorname{Exp}`, `\operatorname{ov}`). GFM's math renderer (KaTeX) often fails to parse it or throws errors under various environments.
    *   *Solution:* Always use `\mathrm{...}` or `\text{...}` instead (e.g., `\mathrm{Exp}`, `\mathrm{ov}`, `\mathrm{RG}`, `\mathrm{tr}`, `\mathrm{diag}`).

---

## 3. LaTeX Validity Checklist

When adding mathematical equations, ensure:
1.  All parentheses `()`, brackets `[]`, and braces `{}` are properly balanced.
2.  Special characters like `%`, `&`, `_`, `#` are properly escaped when not used in their LaTeX command contexts.
3.  Any newline in multiline equations (like inside `align` or `substack` environments) uses `\\` or double backslashes properly and the equation is wrapped in a ` ```math ` block.
4.  No `\operatorname` is used in markdown equations.

---

## 4. Automatic Verification

A Python script is available to automatically validate math syntax in all Markdown files within this workspace:
*   File: [.agents/check_math.py](file:///workspaces/Presentation-MathNet-2026-06-15/.agents/check_math.py)
*   **How to run**: `python3 .agents/check_math.py` or `./.agents/check_math.py`

**Mandatory Rule for Agents:** Whenever you edit or add markdown files containing math formulas, you **must** run this verification script before finishing/pushing to ensure that no rendering errors are introduced.

