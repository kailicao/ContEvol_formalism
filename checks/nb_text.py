"""Print a Mathematica notebook's cells (Input, Text, sections; outputs cut short) as plain text.

Usage, from the repository's root: python3 -I checks/nb_text.py <notebook.nb> [output length]
Box structure is flattened to its strings, so expressions read approximately.
"""
import re
import sys

STYLES = r'(Input|Text|Section|Subsection|Subsubsection|Title|Chapter|Output|Print)'


def match_bracket(s, start):
    depth, j, instr = 0, start, False
    while j < len(s):
        c = s[j]
        if instr:
            if c == '\\':
                j += 2
                continue
            if c == '"':
                instr = False
        else:
            if c == '"':
                instr = True
            elif c == '[':
                depth += 1
            elif c == ']':
                depth -= 1
                if depth == 0:
                    return j
        j += 1
    return -1


def main(path, outlen):
    src = open(path, encoding='utf-8', errors='replace').read()
    pos, out = 0, []
    while True:
        k = src.find('Cell[', pos)
        if k < 0:
            break
        if src.startswith('Cell[CellGroupData', k):
            pos = k + 5
            continue
        e = match_bracket(src, k + 4)
        body = src[k:e + 1]
        sm = re.search(r',\s*"' + STYLES + r'"', body)
        if sm is None:
            pos = k + 5
            continue
        style = sm.group(1)
        strs = re.findall(r'"((?:[^"\\]|\\.)*)"', body[:sm.start()])
        txt = ' '.join(strs)
        txt = (txt.replace('\\[LeftDoubleBracket]', '[[').replace('\\[RightDoubleBracket]', ']]')
               .replace('\\n', '\n').replace('\\<', '').replace('\\>', '').replace('\\"', '"'))
        txt = re.sub(r'\\\[(\w+)\]', r'<\1>', txt)
        if style in ('Output', 'Print'):
            txt = txt[:outlen]
        out.append(f'--- {style}\n{txt}')
        pos = e + 1
    print('\n'.join(out))


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 300)
