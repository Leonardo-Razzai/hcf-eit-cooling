"""Keep only the bibliography entries cited in a LaTeX .aux file and drop DOIs.

Nature-style reference lists carry no DOIs. Run after a first pdflatex pass:
    python paper/tools/trim_bib.py paper/manuscript/main.aux paper/manuscript/refs.bib
The .bib file is rewritten in place; new entries must be added before running.
"""
import re
import sys


def main(aux_path, bib_path):
    aux = open(aux_path, encoding='utf-8').read()
    keys = set()
    for m in re.findall(r'\\citation\{([^}]*)\}', aux):
        keys.update(k.strip() for k in m.split(','))
    bib = open(bib_path, encoding='utf-8').read()
    out = []
    for e in re.split(r'\n(?=@)', bib):
        m = re.match(r'\s*@\w+\{([^,]+),', e)
        if m and m.group(1).strip() in keys:
            out.append(re.sub(r',?\s*\n\s*doi\s*=\s*\{[^}]*\}', '', e).strip())
    found = {re.match(r'@\w+\{([^,]+),', e).group(1).strip() for e in out}
    print(len(keys), 'cited;', len(out), 'kept; missing:', keys - found or 'none')
    open(bib_path, 'w', encoding='utf-8').write('\n\n'.join(out) + '\n')


if __name__ == '__main__':
    main(*sys.argv[1:3])
