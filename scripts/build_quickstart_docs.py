"""Build the public quickstart site with local search and checked download links."""
import argparse
from html import escape
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import shutil
from urllib.parse import unquote, urlsplit

import mistune

ROOT = Path(__file__).resolve().parents[1]


def build(output):
    output.mkdir(parents=True, exist_ok=True)
    sources = {
        'index': ROOT / 'docs/QUICKSTART.md',
        'python-guide': ROOT / 'docs/PYTHON-USER-GUIDE.md',
        'r-guide': ROOT / 'docs/R-USER-GUIDE.md',
        'installation': ROOT / 'docs/installation.md',
        'modeling': ROOT / 'docs/modeling.md',
        'fitting': ROOT / 'docs/fitting.md',
    }
    routes = {path.resolve(): name + '.html' for name, path in sources.items()}
    pages = {name: path.read_text() for name, path in sources.items()}
    navigation = list(pages)

    class Renderer(mistune.HTMLRenderer):
        def local_url(self, url):
            parts = urlsplit(url)
            if parts.scheme or parts.netloc or not parts.path:
                return url
            source = (self.source.parent / unquote(parts.path)).resolve()
            if source.is_relative_to(ROOT) and source.is_dir():
                return 'https://github.com/unibio-intelligence/pharmflux/tree/main/' + str(source.relative_to(ROOT))
            if not source.is_relative_to(ROOT) or not source.is_file():
                raise ValueError(f'Broken documentation link: {self.source.relative_to(ROOT)}: {url}')
            relative = source.relative_to(ROOT)
            if relative.parts[0] in ('internal', '.public-release') or source.name in ('AGENTS.md', 'README.internal.md'):
                raise ValueError(f'Private file linked from documentation: {relative}')
            if source not in routes:
                name = 'reference-' + re.sub(r'[^a-zA-Z0-9_-]', '-', str(relative.with_suffix('')))
                if source.suffix == '.md':
                    routes[source] = name + '.html'
                    sources[name] = source
                    pages[name] = source.read_text()
                else:
                    routes[source] = name + source.suffix
                    shutil.copyfile(source, output / routes[source])
            return routes[source] + ('?' + parts.query if parts.query else '') + ('#' + parts.fragment if parts.fragment else '')

        def heading(self, text, level, **attributes):
            plain = re.sub('<[^>]+>', '', text)
            slug = re.sub(r'[^\w\s-]', '', plain).lower().replace(' ', '-')
            count = self.headings.get(slug, 0)
            self.headings[slug] = count + 1
            anchor = slug + (f'-{count}' if count else '')
            return f'<h{level} id="{escape(anchor)}">{text}</h{level}>\n'

        def link(self, text, url, title=None):
            return super().link(text, self.local_url(url), title)

        def image(self, text, url, title=None):
            return super().image(text, self.local_url(url), title)

    renderer = Renderer(escape=False)
    markdown = mistune.create_markdown(renderer=renderer, plugins=['table'])
    bodies = {}
    while len(bodies) < len(pages):
        name = next(name for name in pages if name not in bodies)
        renderer.source, renderer.headings = sources[name], {}
        bodies[name] = markdown(pages[name])
    nav = ''.join(f'<a href="{name}.html">{escape(pages[name].splitlines()[0].lstrip("# "))}</a>' for name in navigation)
    index = [{'title': text.splitlines()[0].lstrip('# '), 'url': name + '.html', 'text': re.sub(r'[#*`|]', '', text)} for name, text in pages.items()]
    data = json.dumps(index).replace('<', '\\u003c').replace('>', '\\u003e').replace('&', '\\u0026')
    search = '''const input=document.querySelector('#search'), results=document.querySelector('#results');
input.addEventListener('input',()=>{results.replaceChildren();const terms=input.value.toLowerCase().trim().split(/\\s+/).filter(Boolean);if(!terms.length)return;
for(const p of window.searchIndex){const s=p.text.toLowerCase();if(!terms.every(t=>s.includes(t)))continue;
const item=document.createElement('li'),a=document.createElement('a');a.href=p.url;a.textContent=p.title;item.append(a);results.append(item);}});'''
    (output / 'search.js').write_text('window.searchIndex=' + data + ';\n' + search)
    (output / 'style.css').write_text('''body{font:16px system-ui;line-height:1.65;color:#172033;margin:0;background:#fafbfd}header{background:#123047;color:white;padding:18px 24px}header a{color:white}a{color:#165d9c}aside{width:235px;position:fixed;top:88px;bottom:0;overflow:auto;padding:20px}aside a{display:block;margin:8px 0}main{max-width:940px;margin:30px 30px 80px 300px}input{padding:10px;border:1px solid #9aafbd;border-radius:5px;width:min(80%,600px)}pre{overflow:auto;background:#edf2f7;padding:16px;border-radius:6px}code{font-size:14px}table{border-collapse:collapse;display:block;overflow:auto}td,th{padding:8px 12px;border:1px solid #cdd8e1;text-align:left}h1,h2,h3{line-height:1.3}@media(max-width:800px){aside{position:static;width:auto;display:flex;flex-wrap:wrap;gap:10px}main{margin:20px}}''')
    for name, body in bodies.items():
        title = escape(pages[name].splitlines()[0].lstrip('# '))
        document = f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{title} — PharmFlux 0.1.0</title><link rel="stylesheet" href="style.css"><script src="search.js" defer></script></head><body><header><a href="index.html">PharmFlux 0.1.0 documentation</a></header><aside aria-label="Documentation navigation">{nav}</aside><main><label for="search">Search documentation</label><br><input id="search" type="search" placeholder="Search tasks, API names or scientific limits"><ul id="results" aria-live="polite"></ul>{body}</main></body></html>'
        (output / (name + '.html')).write_text(document)
    # Executed HTML uses the same guide-relative links as its Markdown source.
    for guide in ('PYTHON-USER-GUIDE', 'R-USER-GUIDE'):
        source = ROOT / 'docs' / (guide + '.html')
        target = output / routes[source.resolve()]
        renderer.source = ROOT / 'docs' / (guide + '.md')
        content = source.read_text()
        class GuideLinks(HTMLParser):
            def handle_starttag(self, tag, attributes):
                url = dict(attributes).get('href')
                if url is None or urlsplit(url).scheme or url.startswith(('#', '//')):
                    return
                line, column = self.getpos()
                offset = sum(len(value) for value in content.splitlines(keepends=True)[:line - 1]) + column
                original = self.get_starttag_text()
                replacement = re.sub(r'(href=[\"\'])(.*?)([\"\'])',
                                     lambda match: match.group(1) + renderer.local_url(url) + match.group(3), original)
                edits.append((offset, offset + len(original), replacement))
        edits = []
        GuideLinks().feed(content)
        for start, end, replacement in reversed(edits):
            content = content[:start] + replacement + content[end:]
        target.write_text(content)
    (output / '.nojekyll').touch()
    # Parse actual HTML attributes; embedded JavaScript may contain href strings.
    class SiteLinks(HTMLParser):
        def handle_starttag(self, tag, attributes):
            for attribute, url in attributes:
                if attribute not in ('href', 'src') or url is None:
                    continue
                parts = urlsplit(url)
                if parts.scheme or parts.netloc or not parts.path:
                    continue
                if not (page.parent / unquote(parts.path)).is_file():
                    raise ValueError(f'Broken site link: {page.name}: {url}')
    for page in output.glob('*.html'):
        SiteLinks().feed(page.read_text())
    return output


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT / 'docs/_site')
    print(build(parser.parse_args().out))
