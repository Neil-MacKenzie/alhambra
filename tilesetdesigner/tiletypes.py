import copy

# Color dictionary for xgrow colors...
import pkg_resources
import os.path
rgbv = pkg_resources.resource_stream(__name__, os.path.join('data','rgb.txt'))
xcolors={ " ".join(y[3:]): "rgb({},{},{})".format(y[0],y[1],y[2]) for y in [str(x.split()) for x in rgbv] }
del(rgbv)

class tile_daoe(object):
    def __init__(self, defdict):
        self._defdict = defdict
        pass

    def sequence_diagram(self):
        if 'extra' in self._defdict:
            ttype = self['type']+'_'+self['extra']
        else:
            ttype = self['type']
        from lxml import etree
        base_svg = etree.parse( pkg_resources.resource_stream(__name__, os.path.join('seqdiagrambases',ttype+'.svg') ) )

        strings = self._seqdiagseqstrings + \
                [ e for e,t in self.tile_ends if not (t in ('hairpin','blunt','inert')) ] + \
                [ self['name'] ]

        texts = base_svg.findall( "//{http://www.w3.org/2000/svg}text" )

        for text, string in zip(texts,strings):
            text.text = string

        return base_svg.xpath('/svg:svg/svg:g',namespaces={'svg':'http://www.w3.org/2000/svg'})[0]

    @property
    def tile_ends(self):
        return zip(self['ends'],self._endtypes)


    def __getitem__(self, *args):
        return self._defdict.__getitem__(*args)


class tile_daoe_single(tile_daoe):

    """Base class for single DAO-E tiles."""

    def __init__(self, defdict, orient=None):
        """TODO: to be defined1.

        :defdict: TODO

        """
        tile_daoe.__init__(self,defdict)
        self._orient = orient
        self._endtypes = None
    
    def abstract_diagram(self, drawing):
        tilediag = drawing.g()
        
        if 'color' in self._defdict.keys():
            fill = xcolors[self['color']]
        else:
            fill = "rgb(255,255,255)"

        # Tile Box
        tilediag.add( drawing.rect((0,0),(100,100),stroke="black",fill=fill) )

        # End Names
        tilediag.add( drawing.text( self['ends'][0], insert=(50,10), text_anchor='middle', alignment_baseline="middle", font_size="14pt" ) )
        tilediag.add( drawing.text( self['ends'][1], insert=(90,50), text_anchor='middle', alignment_baseline="middle", transform="rotate(90,90,50)", font_size="14pt") )
        tilediag.add( drawing.text( self['ends'][2], insert=(50,90), text_anchor='middle', alignment_baseline="middle", font_size="14pt") )
        tilediag.add( drawing.text( self['ends'][3], insert=(10,50), text_anchor='middle', alignment_baseline="middle", transform="rotate(-90,10,50)", font_size="14pt") ) 

        # Tile Name
        tilediag.add( drawing.text( self['name'], insert=(50,50), text_anchor='middle', alignment_baseline="middle", font_size="14pt" ) )
        
        if self._orient:
            tilediag.add( drawing.text( self._orient[0], insert=(92,8), text_anchor='middle', alignment_baseline="middle", font_size="9pt"))
            tilediag.add( drawing.text( self._orient[1], insert=(8,92), text_anchor='middle', alignment_baseline="middle", font_size="9pt"))

        return tilediag 
    @property
    def _seqdiagseqstrings(self):
        s = self['fullseqs']
        return [ \
             s[0][:5]+"--"+s[0][5:13], \
             s[0][:-6:-1]+"--"+s[0][-6:-14:-1], \
             s[1][:8]+"--"+s[1][8:16], \
             s[1][16:24], \
             s[1][24:32][::-1],
             (s[1][32:40]+"--"+s[1][40:48])[::-1],
             (s[2][:8]+"--"+s[2][8:16])[::-1], \
             s[2][16:24][::-1], \
             s[2][24:32],
             (s[2][32:40]+"--"+s[2][40:48]),
             (s[3][:5]+"--"+s[3][5:13])[::-1], \
             (s[3][:-6:-1]+"--"+s[3][-6:-14:-1])[::-1]         
             ]


class tile_daoe_5up(tile_daoe_single):
    def __init__(self, defdict):
        tile_daoe_single.__init__(self, defdict, orient = ('5','3'))
        self._endtypes = ['TD','TD','DT','DT']

class tile_daoe_3up(tile_daoe_single):
    def __init__(self, defdict):
        tile_daoe_single.__init__(self, defdict, orient = ('3','5'))
        self._endtypes = ['DT','DT','TD','TD']

class tile_daoe_5up_2h(tile_daoe_single):
    def __init__(self, defdict):
        tile_daoe_single.__init__(self, defdict, orient = ('5','3'))
        self._endtypes = ['TD','hairpin','DT','DT']

    @property
    def _seqdiagseqstrings(self):
        s = copy.copy(self['fullseqs'])
        hp = s[3][0:13]; s[3]=s[3][13:]
        return [ \
             s[0][:5]+"--"+s[0][5:13], \
             s[0][:-6:-1]+"--"+s[0][-6:-14:-1], \
             s[1][:8]+"--"+s[1][8:16], \
             s[1][16:24], \
             s[1][24:32][::-1],
             (s[1][32:40]+"--"+s[1][40:48])[::-1],
             (s[2][:8]+"--"+s[2][8:16])[::-1], \
             s[2][16:24][::-1], \
             s[2][24:32],
             (s[2][32:40]+"--"+s[2][40:48]),
             hp[0:5]+'-'+hp[5:7]+'-'+hp[7:9],
             (hp[9:11]+'-'+hp[11:13]+'-')[::-1],
             (s[3][:5]+"--"+s[3][5:13])[::-1], \
             (s[3][:-6:-1]+"--"+s[3][-6:-14:-1])[::-1]         
             ]

class tile_daoe_3up_2h(tile_daoe_single):
    def __init__(self, defdict):
        tile_daoe_single.__init__(self, defdict, orient = ('3','5'))
        self._endtypes = ['DT','hairpin','TD','TD']

    @property
    def _seqdiagseqstrings(self):
        s = copy.copy(self['fullseqs'])
        hp = s[3][26:]; s[3]=s[3][:26]
        return [ \
             s[0][:5]+"--"+s[0][5:13], \
             s[0][:-6:-1]+"--"+s[0][-6:-14:-1], \
             s[1][:8]+"--"+s[1][8:16], \
             s[1][16:24], \
             s[1][24:32][::-1],
             (s[1][32:40]+"--"+s[1][40:48])[::-1],
             (s[2][:8]+"--"+s[2][8:16])[::-1], \
             s[2][16:24][::-1], \
             s[2][24:32],
             (s[2][32:40]+"--"+s[2][40:48]),
             (s[3][:5]+"--"+s[3][5:13])[::-1], \
             (s[3][:-6:-1]+"--"+s[3][-6:-14:-1])[::-1],
             '-'+hp[0:2]+'-'+hp[2:4],
             (hp[4:6]+'-'+hp[6:8]+'-'+hp[8:13])[::-1]
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
        tilediag.add( drawing.rect((0,0),(200,100),stroke="black",fill=fill) )

        # End Names
        tilediag.add( drawing.text( self['ends'][0], insert=(50,10), text_anchor='middle', alignment_baseline="middle", font_size="14pt" ) )
        tilediag.add( drawing.text( self['ends'][1], insert=(150,10), text_anchor='middle', alignment_baseline="middle", font_size="14pt" ) )
        tilediag.add( drawing.text( self['ends'][2], insert=(190,50), text_anchor='middle', alignment_baseline="middle", transform="rotate(90,190,50)", font_size="14pt") )
        tilediag.add( drawing.text( self['ends'][3], insert=(150,90), text_anchor='middle', alignment_baseline="middle", font_size="14pt") )
        tilediag.add( drawing.text( self['ends'][4], insert=(50,90), text_anchor='middle', alignment_baseline="middle", font_size="14pt") )
        tilediag.add( drawing.text( self['ends'][5], insert=(10,50), text_anchor='middle', alignment_baseline="middle", transform="rotate(-90,10,50)", font_size="14pt") ) 

        # Tile Name
        tilediag.add( drawing.text( self['name'], insert=(100,50), text_anchor='middle', alignment_baseline="middle", font_size="14pt" ) )
        
        if self._orient:
            # Orientation
            tilediag.add( drawing.text( self._orient[0], insert=(92,8), text_anchor='middle', alignment_baseline="middle", font_size="9pt"))
            tilediag.add( drawing.text( self._orient[1], insert=(192,8), text_anchor='middle', alignment_baseline="middle", font_size="9pt"))
            tilediag.add( drawing.text( self._orient[1], insert=(8,92), text_anchor='middle', alignment_baseline="middle", font_size="9pt"))
            tilediag.add( drawing.text( self._orient[0], insert=(108,92), text_anchor='middle', alignment_baseline="middle", font_size="9pt"))

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
        tilediag.add( drawing.rect((0,0),(100,200),stroke="black",fill=fill) )

        # End Names
        tilediag.add( drawing.text( self['ends'][0], insert=(50,10), text_anchor='middle', alignment_baseline="middle", font_size="14pt" ) )
        tilediag.add( drawing.text( self['ends'][1], insert=(90,50), text_anchor='middle', alignment_baseline="middle", transform="rotate(90,90,50)", font_size="14pt") )
        tilediag.add( drawing.text( self['ends'][2], insert=(90,150), text_anchor='middle', alignment_baseline="middle", transform="rotate(90,90,150)", font_size="14pt") )
        tilediag.add( drawing.text( self['ends'][3], insert=(50,190), text_anchor='middle', alignment_baseline="middle", font_size="14pt") )
        tilediag.add( drawing.text( self['ends'][4], insert=(10,150), text_anchor='middle', alignment_baseline="middle", transform="rotate(-90,10,150)", font_size="14pt") ) 
        tilediag.add( drawing.text( self['ends'][5], insert=(10,50), text_anchor='middle', alignment_baseline="middle", transform="rotate(-90,10,50)", font_size="14pt") ) 
        
        # Tile Name
        tilediag.add( drawing.text( self['name'], insert=(50,100), text_anchor='middle', alignment_baseline="middle", font_size="14pt" ) )
 
        if self._orient:
            # Orientation
            tilediag.add( drawing.text( self._orient[0], insert=(92,8), text_anchor='middle', alignment_baseline="middle", font_size="9pt"))
            tilediag.add( drawing.text( self._orient[0], insert=(8,192), text_anchor='middle', alignment_baseline="middle", font_size="9pt"))
            tilediag.add( drawing.text( self._orient[1], insert=(92,108), text_anchor='middle', alignment_baseline="middle", font_size="9pt"))
            tilediag.add( drawing.text( self._orient[1], insert=(8,92), text_anchor='middle', alignment_baseline="middle", font_size="9pt")) 
        
        return tilediag

class tile_daoe_doublehoriz_35up(tile_daoe_doublehoriz):
    def __init__(self, defdict):
        tile_daoe_doublehoriz.__init__(self, defdict)
        self._endtypes = ['DT', 'TD', 'TD', 'DT', 'TD', 'TD']
        self._orient = ('3','5')
    
    @property
    def _seqdiagseqstrings(self):
        s = self['fullseqs']
        return [ s[0][:5]+"--"+s[0][5:13], \
              s[0][:-6:-1]+"--"+s[0][-6:-14:-1],
              s[1][:8]+"--"+s[1][8:16], \
              s[1][16:24], \
              s[1][24:32][::-1], \
              (s[1][32:40]+"--"+s[1][40:48])[::-1],
              s[2][:5]+"--"+s[2][5:13], \
              (s[2][13:21]+"--"+s[2][21:26]+"--"+s[2][26:34]+"--"+s[2][34:42])[::-1],
              (s[2][42:50])[::-1],
              s[2][50:58],
              s[2][58:66]+"--"+s[2][66:74],
              (s[3][0:5]+"--"+s[3][5:13])[::-1],
              s[3][13:21]+"--"+s[3][21:26]+"--"+s[3][26:34]+"--"+s[3][34:42],
              s[3][42:50],
              s[3][50:58][::-1],
              (s[3][58:66]+"--"+s[3][66:74])[::-1],
              (s[4][0:8]+"--"+s[4][8:16])[::-1],
              (s[4][16:24])[::-1],
              s[4][24:32],
              s[4][32:40]+"--"+s[4][40:48],
              (s[5][0:5]+"--"+s[5][5:13])[::-1],
              s[5][13:21]+"--"+s[5][21:26]
             ] 

class tile_daoe_doublevert_35up(tile_daoe_doublevert):
    def __init__(self, defdict):
        tile_daoe_doublevert.__init__(self, defdict)
        self._endtypes = ['DT','DT','TD','DT','DT','TD']
        self._orient = ('3','5')


class tile_daoe_doublehoriz_35up_1h2i(tile_daoe_doublehoriz_35up):
    def __init__(self, defdict):
        tile_daoe_doublehoriz_35up.__init__(self, defdict)
        self._endtypes[0]='hairpin'
        self._endtypes[1]='blunt'

    @property
    def _seqdiagseqstrings(self):
        s = self['fullseqs']
        return [ s[0][:5]+"--"+s[0][5:13], \
                (s[0][13:21]+"--"+s[0][21:26]+"-"+s[0][26:28]+"-"+s[0][28:30])[::-1],
              s[1][:8]+"--"+s[1][8:16], \
              s[1][16:24], \
              s[1][24:32][::-1], \
              (s[1][32:40]+"--"+s[1][40:48])[::-1],
              s[2][:5]+"--"+s[2][5:13], \
              (s[2][13:21]+"--"+s[2][21:26]+"--"+s[2][26:34]+"--"+s[2][34:42])[::-1],
              (s[2][42:50])[::-1],
              s[2][50:58],
              s[2][58:66]+"--"+s[2][66:74],
              (s[3][0:5]+"--"+s[3][5:13])[::-1],
              s[3][13:21]+"--"+s[3][21:26]+"--"+s[3][26:34]+"--"+s[3][34:42],
              s[3][42:50],
              s[3][50:58][::-1],
              (s[3][58:66]+"--"+s[3][66:74]+"--"+s[3][74:79])[::-1],
              (s[4][0:8]+"--"+s[4][8:16])[::-1],
              (s[4][16:24])[::-1],
              s[4][24:32],
              s[4][32:40]+"--"+s[4][40:48],
              (s[5][0:5]+"--"+s[5][5:13])[::-1],
              s[5][13:21]+"--"+s[5][21:26],
              s[0][30:32]+"-"+s[0][32:34]+"-"+s[0][34:39]
             ] 

class tile_daoe_doublehoriz_35up_4h5i(tile_daoe_doublehoriz_35up):
    def __init__(self, defdict):
        tile_daoe_doublehoriz_35up.__init__(self, defdict)
        self._endtypes[3]='hairpin'
        self._endtypes[4]='blunt'

    @property
    def _seqdiagseqstrings(self):
        s = self['fullseqs']
        return [ s[0][:5]+"--"+s[0][5:13], \
              s[0][:-6:-1]+"--"+s[0][-6:-14:-1],
              s[1][:8]+"--"+s[1][8:16], \
              s[1][16:24], \
              s[1][24:32][::-1], \
              (s[1][32:40]+"--"+s[1][40:48])[::-1],
              s[2][:5]+"--"+s[2][5:13], \
              (s[2][13:21]+"--"+s[2][21:26]+"--"+s[2][26:34]+"--"+s[2][34:42])[::-1],
              (s[2][42:50])[::-1],
              s[2][50:58],
              s[2][58:66]+"--"+s[2][66:74]+"--"+s[2][74:79],
              (s[3][0:5]+"--"+s[3][5:13])[::-1],
              s[3][13:21]+"--"+s[3][21:26]+"--"+s[3][26:34]+"--"+s[3][34:42],
              s[3][42:50],
              s[3][50:58][::-1],
              (s[3][58:66]+"--"+s[3][66:74])[::-1],
              (s[4][0:8]+"--"+s[4][8:16])[::-1],
              (s[4][16:24])[::-1],
              s[4][24:32],
              s[4][32:40]+"--"+s[4][40:48],
              (s[5][0:5]+"--"+s[5][5:13])[::-1],
              s[5][13:21]+"--"+s[5][21:26]+"-"+s[5][26:28]+"-"+s[5][28:30],
              (s[5][30:32]+"-"+s[5][32:34]+"-"+s[5][34:39])[::-1]
             ] 

class tile_daoe_doublehoriz_35up_2h3h(tile_daoe_doublehoriz_35up):
    def __init__(self, defdict):
        tile_daoe_doublehoriz_35up.__init__(self, defdict)
        self._endtypes[1]='hairpin'; self._endtypes[2]='hairpin'

    @property
    def _seqdiagseqstrings(self):
        s = copy.copy(self['fullseqs'])
        hp2 = s[2][0:13]; s[2]=s[2][13:]
        hp3 = s[5][0:13]; s[5]=s[5][13:]
        return [ s[0][:5]+"--"+s[0][5:13], \
              s[0][:-6:-1]+"--"+s[0][-6:-14:-1],
              s[1][:8]+"--"+s[1][8:16], \
              s[1][16:24], \
              s[1][24:32][::-1], \
              (s[1][32:40]+"--"+s[1][40:48])[::-1],
              (hp2[:5]+'-'+hp2[5:7]+'-'+hp2[7:9])[::-1],
              hp2[9:11]+'-'+hp2[11:13]+'-',
              s[2][:5]+"--"+s[2][5:13], \
              (s[2][13:21]+"--"+s[2][21:26]+"--"+s[2][26:34]+"--"+s[2][34:42])[::-1],
              (s[2][42:50])[::-1],
              s[2][50:58],
              s[2][58:66]+"--"+s[2][66:74],
              (s[3][0:5]+"--"+s[3][5:13])[::-1],
              s[3][13:21]+"--"+s[3][21:26]+"--"+s[3][26:34]+"--"+s[3][34:42],
              s[3][42:50],
              s[3][50:58][::-1],
              (s[3][58:66]+"--"+s[3][66:74])[::-1],
              (s[4][0:8]+"--"+s[4][8:16])[::-1],
              (s[4][16:24])[::-1],
              s[4][24:32],
              s[4][32:40]+"--"+s[4][40:48],
              hp3[:5]+'-'+hp3[5:7]+'-'+hp3[7:9],
              (hp3[9:11]+'-'+hp3[11:13]+'-')[::-1],
              (s[5][0:5]+"--"+s[5][5:13])[::-1],
              s[5][13:21]+"--"+s[5][21:26]
             ] 

class tile_daoe_doublevert_35up_4h5h(tile_daoe_doublevert_35up):
    def __init__(self, defdict):
        tile_daoe_doublevert_35up.__init__(self, defdict)
        self._endtypes[3]='hairpin'; self._endtypes[4]='hairpin'

    @property
    def _seqdiagseqstrings(self):
        s = copy.copy(self['fullseqs'])
        return [ s[0][:5]+"--"+s[0][5:13], 
            s[0][:-6:-1]+"--"+s[0][-6:-14:-1],
            s[1][:8]+"--"+s[1][8:16], 
            s[1][16:24], 
            s[1][24:32][::-1], 
            (s[1][32:40]+'--'+s[1][40:48])[::-1], 
            (s[2][0:8]+"--"+s[2][8:16])[::-1],
            s[2][16:24][::-1],
            s[2][24:32],
            s[2][32:40]+"--"+s[2][40:48]+"--"+s[2][48:53]+"--"+s[2][53:61],
            (s[2][61:69]+"--"+s[2][69:74])[::-1],
            ('-'+s[2][74:76]+'-'+s[2][76:78])[::-1],
            s[2][78:80]+'-'+s[2][80:82]+'-'+s[2][82:87],
            s[3][0:8]+"--"+s[3][8:16],
            s[3][16:24],
            s[3][24:32][::-1],
            (s[3][32:40]+"--"+s[3][40:48]+"--"+s[3][48:53]+"--"+s[3][53:61])[::-1],
            (s[3][61:69]+"--"+s[3][69:74]),
            (s[4][0:8]+"--"+s[4][8:16])[::-1],
            (s[4][16:24])[::-1],
            s[4][24:32],
            s[4][32:40]+"--"+s[4][40:48],
            (s[5][0:5]+"--"+s[5][5:13])[::-1],
            s[5][13:21]+"--"+s[5][21:26],
            "-"+s[5][26:28]+'-'+s[5][28:30], #hp5
            (s[5][30:32]+'-'+s[5][32:34]+'-'+s[5][34:39])[::-1] #hp5
            ]

_tiletypes = { \
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
    def __init__(self, tiletypes = _tiletypes):
        self.tiletypes = tiletypes

    def parse( self, defdict ):
        if 'extra' in defdict:
            ttype = defdict['type']+'_'+defdict['extra']
        else:
            ttype = defdict['type']
        return self.tiletypes[ ttype ](defdict)
       

        
