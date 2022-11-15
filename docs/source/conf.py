# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'Multiphase Test-Bench'
copyright = '2022, Aschmoneit'
author = 'Fynn Aschmoneit'

release = '0.1'
version = '0.1.0'

import os
import sys
from unittest.mock import MagicMock

sys.path.insert(0, os.path.abspath('../../src/MultiphaseTestBench'))
sys.path.insert(0, os.path.abspath('../../src/ValidationTests'))


class Mock(MagicMock):
    @classmethod
    def __getattr__(cls, name):
        return MagicMock()


# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'nbsphinx',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'

MOCK_MODULES = ['numpy', 'scipy']
sys.modules.update((mod_name, Mock()) for mod_name in MOCK_MODULES)

