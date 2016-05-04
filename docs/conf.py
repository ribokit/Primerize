from datetime import date
import sys
import os

sys.path.insert(0, os.path.abspath("../"))

from primerize import __version__

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx.ext.doctest',
    'sphinx.ext.todo',
    'sphinx.ext.mathjax',
    'sphinx.ext.githubpages',
]
templates_path = ['_templates']
html_static_path = ['_static']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
source_suffix = '.rst'
master_doc = 'index'

# General information about the project.
project = u'Primerize'
copyright = u'2008-%s The Board of Trustees of the Leland Stanford Junior University. All Rights Reserved' % date.today().year
author = u'Siqi Tian, Rhiju Das'
language = 'en'

version = __version__
release = version


pygments_style = 'sphinx'
add_function_parentheses = True
add_module_names = True
show_authors = False
todo_include_todos = True
man_show_urls = True

modindex_common_prefix = []

html_domain_indices = True
html_use_smartypants = True
html_use_index = True
html_use_modindex = False
html_split_index = False
html_show_sourcelink = False
html_show_sphinx = True
html_show_copyright = True
html_copy_source = False
html_compact_lists = True
html_last_updated_fmt = '%b %d, %Y'
html_search_language = 'en'

html_short_title = None
html_logo = None
html_favicon = None
html_extra_path = []

html_theme = 'sphinx_rtd_theme'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#html_theme_options = {}

# Add any paths that contain custom themes here, relative to this directory.
#html_theme_path = []

# Custom sidebar templates, maps document names to template names.
#html_sidebars = {}

# Additional templates that should be rendered to pages, maps page names to
# template names.
#html_additional_pages = {}

html_context = {}

# Output file base name for HTML help builder.
htmlhelp_basename = 'Primerize_doc'


autodoc_member_order = 'groupwise'
autosummary_generate = True
napoleon_google_docstring = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_rtype = False
