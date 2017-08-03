#!/usr/bin/env python
#
# report ->
#
# code to generate HTML summary report of SHELXC/D/E results.

from jinja2.environment import Environment
from jinja2.loaders import PackageLoader


def render_html_report(params):
    '''Render Jinja2 template of fast_ep report using SHELXC/D/E results.'''

    env = Environment(loader=PackageLoader('lib.report', 'templates/html'))
    tmpl = env.get_template('fastep_report.html')
    html_page = tmpl.render(**params)
    with open('fastep_report.html', 'wt') as f:
        f.write(html_page)
