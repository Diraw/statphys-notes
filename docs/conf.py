project = "OpenStatPhys Notes"
copyright = "2026, OpenStatPhys"
author = "OpenStatPhys"
release = "0.1.0"

extensions = [
    "myst_parser",
    "sphinx.ext.mathjax",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinxcontrib.bibtex",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

language = "en"

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

master_doc = "index"

html_theme = "sphinx_book_theme"
html_title = "OpenStatPhys Notes"
html_static_path = ["_static"]

html_theme_options = {
    "repository_url": "https://github.com/OpenStatPhys/statphys-notes",
    "use_repository_button": True,
    "use_issues_button": True,
    "home_page_in_toc": True,
}

bibtex_bibfiles = ["refs.bib"]

myst_enable_extensions = [
    "amsmath",
    "dollarmath",
    "colon_fence",
    "deflist",
    "fieldlist",
]

myst_heading_anchors = 3
