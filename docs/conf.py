project = 'VDJ-AQP'
extensions = ['sphinx.ext.githubpages']
html_theme = 'pydata_sphinx_theme'
exclude_patterns = ['_build', 'superpowers/**']

author = 'ISALGO laboratory'
html_static_path = ['_static']
html_theme_options = {
    'logo': {'image_light': '_static/vdjaqp_light.svg', 'image_dark': '_static/vdjaqp_dark.svg'},
    'github_url': 'https://github.com/antigenomics/vdjaqp',
    'show_nav_level': 2,
    'header_links_before_dropdown': 4,
    'navigation_depth': 2,
    'collapse_navigation': False,
    'navbar_align': 'left',
}
