# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------

from datetime import datetime
project = "subsettools"
copyright = f"{datetime.now().year} The Trustees of Princeton University and The Arizona Board of Regents on behalf of the University of Arizona, College of Science Hydrology & Atmospheric Sciences"
author = "George Artavanis, Amanda Triplett"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "myst_nb",
    "autoapi.extension",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

autoapi_dirs = ["../src"]
autoapi_template_dir = './sphinx-autoapi-templates'
autoapi_options = ['members',
                   'undoc-members',
                   'show-inheritance', 
                   'show-module-summary', 
                   'imported-members', 
                  ]
autoapi_add_toctree_entry = False

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

nb_execution_mode = 'off'
