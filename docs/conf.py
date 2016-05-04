from datetime import date
import sys
import os

sys.path.insert(0, os.path.abspath("../"))

from primerize import __version__

extensions = [
    'sphinx.ext.autodoc',
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

html_theme = 'sphinx_rtd_theme'
html_use_smartypants = True
html_domain_indices = True
html_use_index = True
html_split_index = False
html_show_sourcelink = False
html_show_sphinx = True
html_show_copyright = True
html_search_language = 'en'

html_short_title = None
html_logo = None
html_favicon = None
html_extra_path = []

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#html_theme_options = {}

# Add any paths that contain custom themes here, relative to this directory.
#html_theme_path = []

# The name for this set of Sphinx documents.
# "<project> v<release> documentation" by default.
#html_title = u'Primerize v1.2.4'

# If not None, a 'Last updated on:' timestamp is inserted at every page
# bottom, using the given strftime format.
# The empty string is equivalent to '%b %d, %Y'.
#html_last_updated_fmt = None

# Custom sidebar templates, maps document names to template names.
#html_sidebars = {}

# Additional templates that should be rendered to pages, maps page names to
# template names.
#html_additional_pages = {}

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
#html_use_opensearch = ''

# Output file base name for HTML help builder.
htmlhelp_basename = 'Primerize_doc'


# -- Options for manual page output ---------------------------------------
man_pages = [
    (master_doc, 'primerize', u'Primerize Documentation',
     [author], 1)
]
texinfo_documents = [
    (master_doc, 'Primerize', u'Primerize Documentation',
     author, 'Primerize', 'One line description of project.',
     'Miscellaneous'),
]
texinfo_show_urls = 'footnote'

