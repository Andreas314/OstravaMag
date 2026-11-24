# Configuration file for the Sphinx documentation builder.
# Place this file in docs/conf.py
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
project = "OstravaMag"      
author = "Ondřej Faiman"
copyright = "2026, Ondřej Faiman"
release = "0.1"                 

# -- General configuration ---------------------------------------------------
extensions = [
    "myst_parser",                 # Markdown support
    "sphinx.ext.mathjax",          # LaTeX math in theory.md
]

source_suffix = {
    ".md": "markdown",
    ".rst": "restructuredtext",
}

# The root document (must contain the toctree)
root_doc = "index"

# MyST: enable $...$ / $$...$$ math, LaTeX environments, and a few conveniences
myst_enable_extensions = [
    "dollarmath",
    "amsmath",
    "colon_fence",
    "deflist",
    "attrs_inline",
]
myst_heading_anchors = 3           # auto-generate anchors for h1-h3 (for [link](#anchor))

# Files that are only used by the pandoc/LaTeX build, not by Sphinx
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "header.tex",
    "make_doc.sh",
    "title.md",                    # pandoc title block
]

# -- Options for HTML output -------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_static_path = []              # add "_static" if you create docs/_static/

# Figures are referenced relatively (e.g. ![](figures/fe_bands.png) or
# ![](../figures/fe_bands.png) from tutorials/), so no extra config is needed.
