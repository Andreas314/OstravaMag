#!/usr/bin/env bash
set -euo pipefail

OUT="${1:-manual.pdf}"

# Order: title -> theory -> tutorials (basic to advanced) -> settings reference -> references
FILES=(
  title.md
  theory.md
  tutorials/basic.md
  tutorials/minimp.md
  tutorials/qha.md
  settings/inputs.md
  settings/keywords.md
  settings/outputs.md
  references.md
)

# Check that all files exist
for f in "${FILES[@]}"; do
  [[ -f "$f" ]] || { echo "Missing: $f" >&2; exit 1; }
done

# Check for pandoc and a LaTeX engine
command -v pandoc  >/dev/null || { echo "pandoc not found"  >&2; exit 1; }
command -v xelatex >/dev/null || { echo "xelatex not found" >&2; exit 1; }

# Insert a blank line between files so the last line of one never merges with the first of the next.
# --resource-path lets image links like figures/fe_bands.png resolve from the root,
# and also from the subfolders (tutorials/, settings/).
pandoc "${FILES[@]}" \
  --from markdown+tex_math_dollars+pipe_tables+implicit_figures \
  --pdf-engine=xelatex \
  --resource-path=.:figures:tutorials:settings \
  --toc --toc-depth=3 \
  --number-sections \
  -V geometry:margin=2.5cm \
  -V colorlinks=true \
  -V fontsize=11pt \
   -V mainfont="DejaVu Serif" \
  -o "$OUT" \
  --from markdown-blank_before_header-space_in_atx_header+tex_math_dollars+pipe_tables+implicit_figures

echo "Written: $OUT"
