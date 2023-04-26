# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'spacipy'
copyright = '2023, kai han'
author = 'kai han'
release = '0.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

extensions = ['myst_parser']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_theme_options = dict(
    navigation_depth=4,
    logo_only=True
)
html_context = dict(
    display_github=True,  # Integrate GitHub
    github_user='lskfs',  # Username
    github_repo='spacipy',  # Repo name
    github_version='master',  # Version
    conf_py_path='/docs/',  # Path in the checkout to the docs root
)
html_static_path = ['_static']
html_show_sphinx = False
