import copy

from .util import *
import string
import re
from ruamel.yaml.comments import CommentedMap
import collections

# Color dictionary for xgrow colors...
import pkg_resources
import os.path

rgbv = pkg_resources.resource_stream(__name__, os.path.join('data', 'rgb.txt'))
xcolors = {" ".join(y[3:]): "rgb({},{},{})".format(y[0], y[1], y[2])
           for y in [x.decode().split() for x in rgbv]}
del (rgbv)

edp_closetoopen = {
    x: y
    for x, y in zip(string.ascii_lowercase, string.ascii_uppercase)
}
edp_closetoopen.update({')': '(', ']': '[', '}': '{'})


def check_edotparen_consistency(expr):
    expr = expand_compact_edotparen(expr)
    expr = re.sub("\s+", "", expr)
    counts = collections.Counter()
    strand = 0
    strandloc = 0
    for s in expr:
        if s in edp_closetoopen.values():
            counts[s] += 1
        elif s in edp_closetoopen.keys():
            try:
                counts[edp_closetoopen[s]] -= 1
            except KeyError:
                raise ValueError("Opening not found", s, strand, strandloc)
        elif s == ".":
            pass
        elif s == "+":
            strand += 1
            strandloc = 0
            continue
        else:
            raise ValueError("Unknown char", s, strand, strandloc)
        strandloc += 1
    if max(counts.values()) > 0:
        raise ValueError(counts)


def check_edotparen_sequence(edotparen, sequence):
    expr = re.sub("\s+", "", expand_compact_edotparen(edotparen))
    seq = re.sub("\s+", "", sequence).lower()
    if len(expr) != len(seq):
        raise ValueError("Unequal lengths")
    stacks = {}
    strand = 0
    strandloc = 0
    for s, v in zip(expr, seq):
        if s in edp_closetoopen.values():
            if s not in stacks.keys():
                stacks[s] = []
            stacks[s].append(v)
        elif s in edp_closetoopen.keys():
            ss = edp_closetoopen[s]
            if ss not in stacks.keys():
                raise ValueError("Opening not found", s, strand, strandloc)
            vv = stacks[ss].pop()
            if v != wc[vv]:
                raise ValueError(
                    "{} != WC({}) at strand {} loc {} (both from 0)".format(
                        v, vv, strand, strandloc), v, vv, strand, strandloc)
        elif s == ".":
            assert v in wc.keys()
        elif s == "+":
            assert v == "+"
            strand += 1
            strandloc = 0
            continue
        else:
            raise ValueError("Unknown char", s, strand, strandloc)
        strandloc += 1
    if max(len(stack) for stack in stacks.values()) > 0:
        raise ValueError(stacks)


def expand_compact_edotparen(expr):
    return re.sub(r"(\d+)([\[\]\(\)\{\}A-Za-z\.])", lambda m: int(m.group(1)) * m.group(2),
                  expr)


def prettify_edotparen(expr):
    # This is evil:
    return re.sub(r"(([\[\]\(\)\{\}A-Za-z\.])\2+)",
                  lambda m: "{}{}".format(len(m.group(1)), m.group(2)), expr)


class tile_daoe(object):
    def __init__(self, defdict):
        self._defdict = defdict
        pass

    def sequence_diagram(self):
        if 'extra' in self._defdict:
            ttype = self['type'] + '_' + self['extra']
        else:
            ttype = self['type']
        from lxml import etree
        base_svg = etree.parse(
            pkg_resources.resource_stream(__name__, os.path.join(
                'seqdiagrambases', ttype + '.svg')))

        strings = self._seqdiagseqstrings + [
            e for e, t in self.tile_ends
            if not (t in ('hairpin', 'blunt', 'inert'))
        ] + [self['name']]

        texts = base_svg.findall("//{http://www.w3.org/2000/svg}text")

        for text, string in zip(texts, strings):
            if text.getchildren():
                text.getchildren()[0].text = string
            else:
                text.text = string

        return base_svg.xpath(
            '/svg:svg/svg:g',
            namespaces={'svg': 'http://www.w3.org/2000/svg'})[0]

    @property
    def tile_ends(self):
        return zip(self['ends'], self._endtypes)

    def __getitem__(self, *args):
        return self._defdict.__getitem__(*args)

    def check_sequence(self):
        try:
            check_edotparen_sequence(self.edotparen,
                                     "+".join(self['fullseqs']))
        except ValueError as e:
            raise ValueError("{} core is inconsistent.".format(self['name']),
                             self['name']) from e

    def get_end_defs(self):
        es = list()
        for (strand, start, end), endtype, endname in zip(
                self._endlocs, self._endtypes, self['ends']):

            if endtype == 'DT':
                sl = slice(start - 1, end)
            elif endtype == 'TD':
                sl = slice(start, end + 1)
            else:
                sl = None

            if sl and 'fullseqs' in self._defdict.keys():
                seq = self['fullseqs'][strand][sl]

                if endtype == 'DT':
                    seq = (seq + 'n').lower()
                elif endtype == 'TD':
                    seq = ('n' + seq).lower()
            else:
                seq = None

            if endname[-1] == '/':
                endname = endname[:-1]
                # FIXME: this seems wrong. Check.
                # if endtype == 'DT': seq=util.wc[seq[0]]+seq[1:]
                # if endtype == 'TD': seq=seq[:-1]+util.wc[seq[-1]]
                # simple revcomp
                if seq:
                    seq = "".join(reversed([wc[x] for x in seq]))

            e = CommentedMap({'name': endname, 'type': endtype})

            if seq:
                e['fseq'] = seq

            es.append(e)
        return es

    @property
    def orderableseqs(self):
        seqs = copy.deepcopy(self['fullseqs'])
        assert ('label' not in self._defdict.keys())
        return [("{}-{}".format(self['name'], i+1), seq)
                for i, seq in enumerate(seqs)]
    

class tile_daoe_single(tile_daoe):
    """Base class for single DAO-E tiles."""

    def __init__(self, defdict, orient=None):
        """TODO: to be defined1.

        :defdict: TODO

        """
        tile_daoe.__init__(self, defdict)
        self._orient = orient
        self._endtypes = None

    def abstract_diagram(self, drawing):
        tilediag = drawing.g()

        if 'color' in self._defdict.keys():
            fill = xcolors[self['color']]
        else:
            fill = "rgb(255,255,255)"

        # Tile Box
        tilediag.add(
            drawing.rect(
                (0, 0), (100, 100), stroke="black", fill=fill))

        # End Names
        tilediag.add(
            drawing.text(
                self['ends'][0], insert=(50, 10), text_anchor='middle',
                dominant_baseline="mathematical", font_size="14pt"))
        tilediag.add(
            drawing.text(
                self['ends'][1], insert=(90, 50), text_anchor='middle',
                dominant_baseline="mathematical", transform="rotate(90,90,50)",
                font_size="14pt"))
        tilediag.add(
            drawing.text(
                self['ends'][2], insert=(50, 90), text_anchor='middle',
                dominant_baseline="mathematical", font_size="14pt"))
        tilediag.add(
            drawing.text(
                self['ends'][3], insert=(10, 50), text_anchor='middle',
                dominant_baseline="mathematical",
                transform="rotate(-90,10,50)", font_size="14pt"))

        # Tile Name
        tilediag.add(
            drawing.text(
                self['name'],
                insert=(50, 50),
                text_anchor='middle',
                dominant_baseline="mathematical",
                font_size="14pt"))

        if self._orient:
            tilediag.add(
                drawing.text(
                    self._orient[0],
                    insert=(92, 8),
                    text_anchor='middle',
                    dominant_baseline="mathematical",
                    font_size="9pt"))
            tilediag.add(
                drawing.text(
                    self._orient[1],
                    insert=(8, 92),
                    text_anchor='middle',
                    dominant_baseline="mathematical",
                    font_size="9pt"))

        return tilediag

    @property
    def _seqdiagseqstrings(self):
        s = self['fullseqs']
        return [
            s[0][:5] + "--" + s[0][5:13],
            s[0][:-6:-1] + "--" + s[0][-6:-14:-1],
            s[1][:8] + "--" + s[1][8:16], s[1][16:24], s[1][24:32][::-1],
            (s[1][32:40] + "--" + s[1][40:48])[::-1],
            (s[2][:8] + "--" + s[2][8:16])[::-1], s[2][16:24][::-1],
            s[2][24:32], (s[2][32:40] + "--" + s[2][40:48]),
            (s[3][:5] + "--" + s[3][5:13])[::-1],
            (s[3][:-6:-1] + "--" + s[3][-6:-14:-1])[::-1]
        ]

    @property
    def _short_bound_full(self):
        s = self['fullseqs']
        return [s[0][5:-5], s[3][5:-5]]

    @property
    def _side_bound_regions(self):
        s = self['fullseqs']
        return [s[0][5:5 + 8], s[0][-5 - 8:-5], s[3][5:5 + 8], s[3][-5 - 8:-5]]


class tile_daoe_5up(tile_daoe_single):
    def __init__(self, defdict):
        tile_daoe_single.__init__(self, defdict, orient=('5', '3'))
        self._endtypes = ['TD', 'TD', 'DT', 'DT']
        # endlocs is strand, loc, length
        self._endlocs = [(0, 0, 5), (3, 0, 5), (3, 21, None), (0, 21, None)]
    # valid edotparen for both 3up and 5up
    edotparen = "5.16(5.+8)16[16{8)+8(16]16}8(+5.16)5."

    @property
    def orderableseqs(self):
        seqs = copy.deepcopy(self['fullseqs'])
        if 'label' in self._defdict.keys() and self['label'] == 'both':
            assert seqs[1][16] == 'T'
            assert seqs[2][16] == 'T'
            seqs[1] = seqs[1][:16]+'/iBiodT/'+seqs[1][17:]
            seqs[2] = seqs[2][:16]+'/iBiodT/'+seqs[2][17:]
            assert re.sub('/iBiodT/', 'T', seqs[1]) == self['fullseqs'][1]
            assert re.sub('/iBiodT/', 'T', seqs[2]) == self['fullseqs'][2]
        return [("{}-{}".format(self['name'], i+1), seq)
                for i, seq in enumerate(seqs)]


class tile_daoe_3up(tile_daoe_single):
    def __init__(self, defdict):
        tile_daoe_single.__init__(self, defdict, orient=('3', '5'))
        self._endtypes = ['DT', 'DT', 'TD', 'TD']
        self._endlocs = [(0, 21, None), (3, 21, None), (3, 0, 5), (0, 0, 5)]
    # valid edotparen for both 3up and 5up
    edotparen = "5.16(5.+8)16[16{8)+8(16]16}8(+5.16)5."

    @property
    def orderableseqs(self):
        seqs = copy.deepcopy(self['fullseqs'])
        if 'label' in self._defdict.keys() and self['label'] == 'both':
            assert seqs[1][31] == 'T'
            assert seqs[2][31] == 'T'
            seqs[1] = seqs[1][:31]+'/iBiodT/'+seqs[1][32:]
            seqs[2] = seqs[2][:31]+'/iBiodT/'+seqs[2][32:]
            assert re.sub('/iBiodT/', 'T', seqs[1]) == self['fullseqs'][1]
            assert re.sub('/iBiodT/', 'T', seqs[2]) == self['fullseqs'][2]
        return [("{}-{}".format(self['name'], i+1), seq)
                for i, seq in enumerate(seqs)]


class tile_daoe_5up_2h(tile_daoe_single):
    def __init__(self, defdict):
        tile_daoe_single.__init__(self, defdict, orient=('5', '3'))
        self._endtypes = ['TD', 'hairpin', 'DT', 'DT']
        self._endlocs = [(0, 0, 5), (3, 0, 18), (3, 21, None), (0, 21, None)]

    edotparen = "5.16(5.+8)16[16{8)+8(16]16}8(+5.16)7(4.7)"

    @property
    def _short_bound_full(self):
        s = self['fullseqs']
        return [s[0][5:-5], s[3][18:-5]]

    @property
    def _side_bound_regions(self):
        s = self['fullseqs']
        return [s[0][5:5 + 8], s[0][-5 - 8:-5], s[3][18:18 + 8],
                s[3][-5 - 8:-5]]

    @property
    def _seqdiagseqstrings(self):
        s = copy.copy(self['fullseqs'])
        hp = s[3][0:13]
        s[3] = s[3][13:]
        return [
            s[0][:5] + "--" + s[0][5:13],
            s[0][:-6:-1] + "--" + s[0][-6:-14:-1],
            s[1][:8] + "--" + s[1][8:16], s[1][16:24], s[1][24:32][::-1],
            (s[1][32:40] + "--" + s[1][40:48])[::-1],
            (s[2][:8] + "--" + s[2][8:16])[::-1], s[2][16:24][::-1],
            s[2][24:32], (s[2][32:40] + "--" + s[2][40:48]),
            hp[0:5] + '-' + hp[5:7] + '-' + hp[7:9],
            (hp[9:11] + '-' + hp[11:13] + '-')[::-1],
            (s[3][:5] + "--" + s[3][5:13])[::-1],
            (s[3][:-6:-1] + "--" + s[3][-6:-14:-1])[::-1]
        ]


class tile_daoe_3up_2h(tile_daoe_single):
    @property
    def _short_bound_full(self):
        s = self['fullseqs']
        return [s[0][5:-5], s[3][5:-18]]

    @property
    def _side_bound_regions(self):
        s = self['fullseqs']
        return [s[0][5:5 + 8], s[0][-5 - 8:-5], s[3][5:5 + 8],
                s[3][-18 - 8:-18]]

    def __init__(self, defdict):
        tile_daoe_single.__init__(self, defdict, orient=('3', '5'))
        self._endtypes = ['DT', 'hairpin', 'TD', 'TD']
        self._endlocs = [(0, 21, None), (3, 21, None), (3, 0, 5), (0, 0, 5)]

    @property
    def _seqdiagseqstrings(self):
        s = copy.copy(self['fullseqs'])
        hp = s[3][26:]
        s[3] = s[3][:26]
        return [
            s[0][:5] + "--" + s[0][5:13],
            s[0][:-6:-1] + "--" + s[0][-6:-14:-1],
            s[1][:8] + "--" + s[1][8:16], s[1][16:24], s[1][24:32][::-1],
            (s[1][32:40] + "--" + s[1][40:48])[::-1],
            (s[2][:8] + "--" + s[2][8:16])[::-1], s[2][16:24][::-1],
            s[2][24:32], (s[2][32:40] + "--" + s[2][40:48]),
            (s[3][:5] + "--" + s[3][5:13])[::-1],
            (s[3][:-6:-1] + "--" + s[3][-6:-14:-1])[::-1],
            '-' + hp[0:2] + '-' + hp[2:4],
            (hp[4:6] + '-' + hp[6:8] + '-' + hp[8:13])[::-1]
        ]


class tile_daoe_doublehoriz(tile_daoe):
    def __init__(self, defdict):
        tile_daoe.__init__(self, defdict)
        self._orient = None

    def abstract_diagram(self, drawing):
        tilediag = drawing.g()

        if 'color' in self._defdict.keys():
            fill = xcolors[self['color']]
        else:
            fill = "rgb(255,255,255)"

        # Tile Box
        tilediag.add(
            drawing.rect(
                (0, 0), (200, 100), stroke="black", fill=fill))

        # End Names
        tilediag.add(
            drawing.text(
                self['ends'][0], insert=(50, 10), text_anchor='middle',
                dominant_baseline="mathematical", font_size="14pt"))
        tilediag.add(
            drawing.text(
                self['ends'][1], insert=(150, 10), text_anchor='middle',
                dominant_baseline="mathematical", font_size="14pt"))
        tilediag.add(
            drawing.text(
                self['ends'][2], insert=(190, 50), text_anchor='middle',
                dominant_baseline="mathematical",
                transform="rotate(90,190,50)", font_size="14pt"))
        tilediag.add(
            drawing.text(
                self['ends'][3], insert=(150, 90), text_anchor='middle',
                dominant_baseline="mathematical", font_size="14pt"))
        tilediag.add(
            drawing.text(
                self['ends'][4], insert=(50, 90), text_anchor='middle',
                dominant_baseline="mathematical", font_size="14pt"))
        tilediag.add(
            drawing.text(
                self['ends'][5], insert=(10, 50), text_anchor='middle',
                dominant_baseline="mathematical",
                transform="rotate(-90,10,50)", font_size="14pt"))

        # Tile Name
        tilediag.add(
            drawing.text(
                self['name'],
                insert=(100, 50),
                text_anchor='middle',
                dominant_baseline="mathematical",
                font_size="14pt"))

        if self._orient:
            # Orientation
            tilediag.add(
                drawing.text(
                    self._orient[0],
                    insert=(92, 8),
                    text_anchor='middle',
                    dominant_baseline="mathematical",
                    font_size="9pt"))
            tilediag.add(
                drawing.text(
                    self._orient[1],
                    insert=(192, 8),
                    text_anchor='middle',
                    dominant_baseline="mathematical",
                    font_size="9pt"))
            tilediag.add(
                drawing.text(
                    self._orient[1],
                    insert=(8, 92),
                    text_anchor='middle',
                    dominant_baseline="mathematical",
                    font_size="9pt"))
            tilediag.add(
                drawing.text(
                    self._orient[0],
                    insert=(108, 92),
                    text_anchor='middle',
                    dominant_baseline="mathematical",
                    font_size="9pt"))

        return tilediag


class tile_daoe_doublevert(tile_daoe):
    def __init__(self, defdict):
        tile_daoe.__init__(self, defdict)
        self._orient = None

    def abstract_diagram(self, drawing):
        tilediag = drawing.g()

        if 'color' in self._defdict.keys():
            fill = xcolors[self['color']]
        else:
            fill = "rgb(255,255,255)"

        # Tile Box
        tilediag.add(
            drawing.rect(
                (0, 0), (100, 200), stroke="black", fill=fill))

        # End Names
        tilediag.add(
            drawing.text(
                self['ends'][0], insert=(50, 10), text_anchor='middle',
                dominant_baseline="mathematical", font_size="14pt"))
        tilediag.add(
            drawing.text(
                self['ends'][1], insert=(90, 50), text_anchor='middle',
                dominant_baseline="mathematical", transform="rotate(90,90,50)",
                font_size="14pt"))
        tilediag.add(
            drawing.text(
                self['ends'][2], insert=(90, 150), text_anchor='middle',
                dominant_baseline="mathematical",
                transform="rotate(90,90,150)", font_size="14pt"))
        tilediag.add(
            drawing.text(
                self['ends'][3], insert=(50, 190), text_anchor='middle',
                dominant_baseline="mathematical", font_size="14pt"))
        tilediag.add(
            drawing.text(
                self['ends'][4], insert=(10, 150), text_anchor='middle',
                dominant_baseline="mathematical",
                transform="rotate(-90,10,150)", font_size="14pt"))
        tilediag.add(
            drawing.text(
                self['ends'][5], insert=(10, 50), text_anchor='middle',
                dominant_baseline="mathematical",
                transform="rotate(-90,10,50)", font_size="14pt"))

        # Tile Name
        tilediag.add(
            drawing.text(
                self['name'],
                insert=(50, 100),
                text_anchor='middle',
                dominant_baseline="mathematical",
                font_size="14pt"))

        if self._orient:
            # Orientation
            tilediag.add(
                drawing.text(
                    self._orient[0],
                    insert=(92, 8),
                    text_anchor='middle',
                    dominant_baseline="mathematical",
                    font_size="9pt"))
            tilediag.add(
                drawing.text(
                    self._orient[0],
                    insert=(8, 192),
                    text_anchor='middle',
                    dominant_baseline="mathematical",
                    font_size="9pt"))
            tilediag.add(
                drawing.text(
                    self._orient[1],
                    insert=(92, 108),
                    text_anchor='middle',
                    dominant_baseline="mathematical",
                    font_size="9pt"))
            tilediag.add(
                drawing.text(
                    self._orient[1],
                    insert=(8, 92),
                    text_anchor='middle',
                    dominant_baseline="mathematical",
                    font_size="9pt"))

        return tilediag


class tile_daoe_doublehoriz_35up(tile_daoe_doublehoriz):
    def __init__(self, defdict):
        tile_daoe_doublehoriz.__init__(self, defdict)
        self._endtypes = ['DT', 'TD', 'TD', 'DT', 'TD', 'TD']
        self._orient = ('3', '5')
        self._endlocs = [(0, -5, None), (2, 0, 5), (5, 0, 5), (5, -5, None),
                         (3, 0, 5), (0, 0, 5)]

    @property
    def _short_bound_full(self):
        s = self['fullseqs']
        el = self._endlocs
        return [s[0][el[5][2]:el[0][1]], s[5][el[2][2]:el[3][1]]]

    @property
    def _side_bound_regions(self):
        s = self['fullseqs']
        el = self._endlocs

        def g(e):
            if e[2] is None or e[2] == len(s[e[0]]):
                return s[e[0]][-8 + e[1]:e[1]]
            elif e[1] == 0:
                return s[e[0]][e[2]:e[2] + 8]
            else:
                raise ValueError()

        return [g(e) for e in el]

    @property
    def _seqdiagseqstrings(self):
        s = self['fullseqs']
        return [s[0][:5] + "--" + s[0][5:13],
                s[0][:-6:-1] + "--" + s[0][-6:-14:-1],
                s[1][:8] + "--" + s[1][8:16], s[1][16:24], s[1][24:32][::-1],
                (s[1][32:40] + "--" + s[1][40:48])[::-1],
                s[2][:5] + "--" + s[2][5:13],
                (s[2][13:21] + "--" + s[2][21:26] + "--" + s[2][26:34] + "--" +
                 s[2][34:42])[::-1], (s[2][42:50])[::-1], s[2][50:58],
                s[2][58:66] + "--" + s[2][66:74],
                (s[3][0:5] + "--" + s[3][5:13])[::-1], s[3][13:21] + "--" +
                s[3][21:26] + "--" + s[3][26:34] + "--" + s[3][34:42],
                s[3][42:50], s[3][50:58][::-1],
                (s[3][58:66] + "--" + s[3][66:74])[::-1],
                (s[4][0:8] + "--" + s[4][8:16])[::-1], (s[4][16:24])[::-1],
                s[4][24:32], s[4][32:40] + "--" + s[4][40:48],
                (s[5][0:5] + "--" + s[5][5:13])[::-1],
                s[5][13:21] + "--" + s[5][21:26]]


class tile_daoe_doublevert_35up(tile_daoe_doublevert):
    @property
    def _short_bound_full(self):
        s = self['fullseqs']
        el = self._endlocs
        return [s[0][el[5][2]:el[0][1]], s[5][el[2][2]:el[3][1]]]

    @property
    def _side_bound_regions(self):
        s = self['fullseqs']
        el = self._endlocs

        def g(e):
            if e[2] is None or e[2] == len(s[e[0]]):
                return s[e[0]][-8 + e[1]:e[1]]
            elif e[1] == 0:
                return s[e[0]][e[2]:e[2] + 8]
            else:
                raise ValueError()

        return [g(e) for e in el]

    def __init__(self, defdict):
        tile_daoe_doublevert.__init__(self, defdict)
        self._endtypes = ['DT', 'DT', 'TD', 'DT', 'DT', 'TD']
        self._orient = ('3', '5')
        self._endlocs = [(0, -5, None), (3, -5, None), (5, 0, 5),
                         (5, -5, None), (2, -5, None), (0, 0, 5)]


class tile_daoe_doublehoriz_35up_1h2i(tile_daoe_doublehoriz_35up):
    edotparen = '5.16(7(4.7)+8)16[16{8)+5(29(16]16}8(+5.29)16[16{8)5)+8(16]16}8(+5.16)5.'

    def __init__(self, defdict):
        tile_daoe_doublehoriz_35up.__init__(self, defdict)
        self._endtypes[0] = 'hairpin'
        self._endtypes[1] = 'blunt'

    @property
    def _seqdiagseqstrings(self):
        s = self['fullseqs']
        return [s[0][:5] + "--" + s[0][5:13],
                (s[0][13:21] + "--" + s[0][21:26] + "-" + s[0][26:28] + "-" +
                 s[0][28:30])[::-1], s[1][:8] + "--" + s[1][8:16], s[1][16:24],
                s[1][24:32][::-1], (s[1][32:40] + "--" + s[1][40:48])[::-1],
                s[2][:5] + "--" + s[2][5:13],
                (s[2][13:21] + "--" + s[2][21:26] + "--" + s[2][26:34] + "--" +
                 s[2][34:42])[::-1], (s[2][42:50])[::-1], s[2][50:58],
                s[2][58:66] + "--" + s[2][66:74],
                (s[3][0:5] + "--" + s[3][5:13])[::-1], s[3][13:21] + "--" +
                s[3][21:26] + "--" + s[3][26:34] + "--" + s[3][34:42],
                s[3][42:50], s[3][50:58][::-1],
                (s[3][58:66] + "--" + s[3][66:74] + "--" + s[3][74:79])[::-1],
                (s[4][0:8] + "--" + s[4][8:16])[::-1], (s[4][16:24])[::-1],
                s[4][24:32], s[4][32:40] + "--" + s[4][40:48],
                (s[5][0:5] + "--" + s[5][5:13])[::-1],
                s[5][13:21] + "--" + s[5][21:26],
                s[0][30:32] + "-" + s[0][32:34] + "-" + s[0][34:39]]


class tile_daoe_doublehoriz_35up_4h5i(tile_daoe_doublehoriz_35up):
    edotparen = '5.16(5.+8)16[16{8)+5.29(16]16}8(5(+5)29)16[16{8)+8(16]16}8(+5.16)7(4.7)'

    def __init__(self, defdict):
        tile_daoe_doublehoriz_35up.__init__(self, defdict)
        self._endtypes[3] = 'hairpin'
        self._endtypes[4] = 'blunt'

    @property
    def _seqdiagseqstrings(self):
        s = self['fullseqs']
        return [s[0][:5] + "--" + s[0][5:13],
                s[0][:-6:-1] + "--" + s[0][-6:-14:-1],
                s[1][:8] + "--" + s[1][8:16], s[1][16:24], s[1][24:32][::-1],
                (s[1][32:40] + "--" + s[1][40:48])[::-1],
                s[2][:5] + "--" + s[2][5:13],
                (s[2][13:21] + "--" + s[2][21:26] + "--" + s[2][26:34] + "--" +
                 s[2][34:42])[::-1], (s[2][42:50])[::-1], s[2][50:58],
                s[2][58:66] + "--" + s[2][66:74] + "--" + s[2][74:79],
                (s[3][0:5] + "--" + s[3][5:13])[::-1], s[3][13:21] + "--" +
                s[3][21:26] + "--" + s[3][26:34] + "--" + s[3][34:42],
                s[3][42:50], s[3][50:58][::-1],
                (s[3][58:66] + "--" + s[3][66:74])[::-1],
                (s[4][0:8] + "--" + s[4][8:16])[::-1], (s[4][16:24])[::-1],
                s[4][24:32], s[4][32:40] + "--" + s[4][40:48],
                (s[5][0:5] + "--" + s[5][5:13])[::-1], s[5][13:21] + "--" +
                s[5][21:26] + "-" + s[5][26:28] + "-" + s[5][28:30],
                (s[5][30:32] + "-" + s[5][32:34] + "-" + s[5][34:39])[::-1]]


class tile_daoe_doublehoriz_35up_2h3h(tile_daoe_doublehoriz_35up):
    def __init__(self, defdict):
        tile_daoe_doublehoriz_35up.__init__(self, defdict)
        self._endtypes[1] = 'hairpin'
        self._endtypes[2] = 'hairpin'
        self._endlocs = [(0, -5, None), (2, 0, 18), (5, 0, 18), (5, -5, None),
                         (3, 0, 5), (0, 0, 5)]

    edotparen = '5.16(5.+8)16[16{8)+7(4.7)29(16]16}8(+5.29)16[16{8)+8(16]16}8(+7(4.23)5.'

    @property
    def _seqdiagseqstrings(self):
        s = copy.copy(self['fullseqs'])
        hp2 = s[2][0:13]
        s[2] = s[2][13:]
        hp3 = s[5][0:13]
        s[5] = s[5][13:]
        return [
            s[0][:5] + "--" + s[0][5:13],
            s[0][:-6:-1] + "--" + s[0][-6:-14:-1],
            s[1][:8] + "--" + s[1][8:16], s[1][16:24], s[1][24:32][::-1],
            (s[1][32:40] + "--" + s[1][40:48])[::-1],
            (hp2[:5] + '-' + hp2[5:7] + '-' + hp2[7:9])[::-1],
            hp2[9:11] + '-' + hp2[11:13] + '-', s[2][:5] + "--" + s[2][5:13],
            (s[2][13:21] + "--" + s[2][21:26] + "--" + s[2][26:34] + "--" +
             s[2][34:42])[::-1], (s[2][42:50])[::-1], s[2][50:58],
            s[2][58:66] + "--" + s[2][66:74],
            (s[3][0:5] + "--" + s[3][5:13])[::-1], s[3][13:21] + "--" +
            s[3][21:26] + "--" + s[3][26:34] + "--" + s[3][34:42], s[3][42:50],
            s[3][50:58][::-1], (s[3][58:66] + "--" + s[3][66:74])[::-1],
            (s[4][0:8] + "--" + s[4][8:16])[::-1], (s[4][16:24])[::-1],
            s[4][24:32], s[4][32:40] + "--" + s[4][40:48],
            hp3[:5] + '-' + hp3[5:7] + '-' + hp3[7:9],
            (hp3[9:11] + '-' + hp3[11:13] + '-')[::-1],
            (s[5][0:5] + "--" + s[5][5:13])[::-1],
            s[5][13:21] + "--" + s[5][21:26]
        ]


class tile_daoe_doublevert_35up_4h5h(tile_daoe_doublevert_35up):
    def __init__(self, defdict):
        tile_daoe_doublevert_35up.__init__(self, defdict)
        self._endtypes[3] = 'hairpin'
        self._endtypes[4] = 'hairpin'
        self._endlocs = [(0, -5, None), (3, -5, None), (5, 0, 5),
                         (5, -18, None), (2, -18, None), (0, 0, 5)]

    edotparen = "5.16(5.+8)16[16{8)+8(16]16}8(5(16(7(4.7)+8)16[16{8)5)16)5.+8(16]16}8(+5.16)7(4.7)"

    @property
    def _seqdiagseqstrings(self):
        s = copy.copy(self['fullseqs'])
        return [
            s[0][:5] + "--" + s[0][5:13],
            s[0][:-6:-1] + "--" + s[0][-6:-14:-1],
            s[1][:8] + "--" + s[1][8:16],
            s[1][16:24],
            s[1][24:32][::-1],
            (s[1][32:40] + '--' + s[1][40:48])[::-1],
            (s[2][0:8] + "--" + s[2][8:16])[::-1],
            s[2][16:24][::-1],
            s[2][24:32],
            s[2][32:40] + "--" + s[2][40:48] + "--" + s[2][48:53] + "--" +
            s[2][53:61],
            (s[2][61:69] + "--" + s[2][69:74])[::-1],
            ('-' + s[2][74:76] + '-' + s[2][76:78])[::-1],
            s[2][78:80] + '-' + s[2][80:82] + '-' + s[2][82:87],
            s[3][0:8] + "--" + s[3][8:16],
            s[3][16:24],
            s[3][24:32][::-1],
            (s[3][32:40] + "--" + s[3][40:48] + "--" + s[3][48:53] + "--" +
             s[3][53:61])[::-1],
            (s[3][61:69] + "--" + s[3][69:74]),
            (s[4][0:8] + "--" + s[4][8:16])[::-1],
            (s[4][16:24])[::-1],
            s[4][24:32],
            s[4][32:40] + "--" + s[4][40:48],
            (s[5][0:5] + "--" + s[5][5:13])[::-1],
            s[5][13:21] + "--" + s[5][21:26],
            "-" + s[5][26:28] + '-' + s[5][28:30],  # hp5
            (s[5][30:32] + '-' + s[5][32:34] + '-' + s[5][34:39])[::-1]  # hp5
        ]


_tiletypes = {
    'tile_daoe_5up': tile_daoe_5up,
    'tile_daoe_3up': tile_daoe_3up,
    'tile_daoe_5up_2h': tile_daoe_5up_2h,
    'tile_daoe_3up_2h': tile_daoe_3up_2h,
    'tile_daoe_doublehoriz_35up': tile_daoe_doublehoriz_35up,
    'tile_daoe_doublehoriz_35up_2h3h': tile_daoe_doublehoriz_35up_2h3h,
    'tile_daoe_doublehoriz_35up_1h2i': tile_daoe_doublehoriz_35up_1h2i,
    'tile_daoe_doublehoriz_35up_4h5i': tile_daoe_doublehoriz_35up_4h5i,
    'tile_daoe_doublehoriz_35up_4h5b': tile_daoe_doublehoriz_35up_4h5i,
    'tile_daoe_doublevert_35up_4h5h': tile_daoe_doublevert_35up_4h5h
}


class TileFactory(object):
    def __init__(self, tiletypes=_tiletypes):
        self.tiletypes = tiletypes

    def parse(self, defdict):
        if 'extra' in defdict:
            ttype = defdict['type'] + '_' + defdict['extra']
        else:
            ttype = defdict['type']
        return self.tiletypes[ttype](defdict)

    # Utility Functions


def get_tile_ends(tile):
    r = list()
    for x in zip(tile['ends'], tileendtypes[tile['type']]):
        if x[0] == 'hp':
            x = (x[0], 'hairpin')
        r.append(x)
    return r


def compname(endname):
    if endname[-1] == '/':
        return endname[:-1]
    else:
        return endname + '/'


def order_pepper_strands(strandlist):
    # We're going to assume, for now, that they're already ordered by Pepper.
    # We'll see.
    return [strandseq for (strandname, strandseq) in strandlist]


def gettile(tset, tname):
    foundtiles = [x for x in tset['tiles'] if x['name'] == tname]
    assert len(foundtiles) == 1
    return foundtiles[0]

# default tilefactory
tfactory = TileFactory()


def endlist_from_tilelist(tilelist, fail_immediate=True, tilefactory=tfactory):
    """\
Given a named_list of tiles (or just a list, for now), extract the sticky ends 
from each tile, and merge these (using merge_ends) into a named_list of sticky
ends.

Parameters
----------

tilelist: either a list of tile description dicts, or a list of tile instances.
    If 

fail_immediate: (default True) if True, immediately fail on a failure,
    with ValueError( tilename, exception ) from exception  If False, collect 
    exceptions, then raise ValueError( "message", [(tilename, exception)] , 
    output ).

tilefactory: (default tiletypes.tfactory) if tilelist is a list of tile 
    descriptions, use this TileFactory to parse them.
"""
    endlist = named_list()
    errors = []

    for tile in tilelist:
        try:
            if not isinstance(tile, tile_daoe):
                tile = tilefactory.parse(tile)
            ends = tile.get_end_defs()
            endlist = merge_endlists(endlist, ends, in_place=True)
        except BaseException as e:
            if fail_immediate:
                raise ValueError(tile['name']) from e
            else:
                errors.append((tile['name'], e))

    if errors:
        raise ValueError("End list generation failed on:", errors, endlist)

    return endlist
