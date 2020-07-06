# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
#
# Astropy documentation build configuration file.
#
# This file is execfile()d with the current directory set to its containing dir.
#
# Note that not all possible configuration values are present in this file.
#
# All configuration values have a default. Some values are defined in
# the global Astropy configuration which is loaded here before anything else. 
# See astropy.sphinx.conf for which values are set there.

import datetime
#import sys
#import os

import sphinx_rtd_theme

#sys.path.insert(0, os.path.abspath("static"))

# -- General configuration ----------------------------------------------------
# If your documentation needs a minimal Sphinx version, state it here.

needs_sphinx = '1.3'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.inheritance_diagram',
    'sphinx.ext.mathjax',
    'sphinx.ext.extlinks',
#    'sphinx.ext.linkcode',
    'jupyter_sphinx',
    'sphinx.ext.doctest',
    'numpydoc',
]


intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'numpy': ('https://docs.scipy.org/doc/numpy/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/reference/', None),
    'astropy': ('http://docs.astropy.org/en/stable/', None),
    'synphot': ('https://synphot.readthedocs.io/en/latest/', None),
    }

extlinks = {'python': ('https://docs.python.org/3/%s', None),
            'numpy': ('https://docs.scipy.org/doc/numpy/', None),
            'scipy': ('https://docs.scipy.org/doc/scipy/reference/', None),
            'astropy': ('http://docs.astropy.org/en/stable/', None),
            'synphot': ('https://synphot.readthedocs.io/en/latest/', None),
            }


numpydoc_show_class_members = False
autosummary_generate = True
autoclass_content = "class"
autodoc_default_flags = ["members", "inherited-members"]
autodoc_docstring_signature = False

source_suffix = '.rst'
source_encoding = 'utf-8-sig'

master_doc = 'index'


# -- Project information ------------------------------------------------------

# This does not *have* to match the package name, but typically does
project = 'spextra'
author = 'Miguel Verdugo'
current_year = datetime.datetime.now().year
copyright = '2020-{:d}, {}'.format(current_year, author)


# import spextra
# The short X.Y version.
version = "0.1" #spextra.__version__.split('-', 1)[0]
release = version
# The full version, including alpha/beta/rc tags.

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_build']

# The reST default role (used for this markup: `text`) to use for all
# documents.

default_role = 'obj'

# The name of the Pygments (syntax highlighting) style to use.

pygments_style = 'default' #'sphinx'

# -- Options for HTML output ---------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = [sphinx_gallery.glr_path_static()]


# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
#html_last_updated_fmt = ''

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
html_title = '{0} v{1}'.format(project, release)

# Output file base name for HTML help builder.
htmlhelp_basename = project + 'doc'

# Add local templates path to modify autosummary templates
#templates_path = ['_templates']

# Static files to copy after template files
html_static_path = ['_static']


# -- Options for LaTeX output --------------------------------------------------

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, documentclass [howto/manual]).

latex_documents = [('index', project + '.tex', project + u' Documentation',
                    author, 'manual')]


# -- Options for manual page output --------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).

man_pages = [('index', project.lower(), project + u' Documentation',
              [author], 1)]

