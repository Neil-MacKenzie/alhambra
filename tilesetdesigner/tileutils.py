tileendtypes = {
        'tile_daoe_5up': ['TD','TD','DT','DT'],
        'tile_daoe_3up': ['DT','DT','TD','TD'],
        'tile_daoe_doublehoriz_35up': ['DT', 'TD', 'TD', 'DT', 'TD', 'TD'],
        'tile_daoe_doublevert_35up': ['DT','DT','TD','DT','DT','TD']
        }

def get_tile_ends(tile):
    r = list()
    for x in zip(tile['ends'],tileendtypes[tile['type']]):
        if x[0]=='hp':
            x = (x[0],'hairpin')
        r.append(x)
    return r

def compname( endname ):
    if endname[-1] == '/':
        return endname[:-1]
    else:
        return endname+'/'

def order_pepper_strands( strandlist ):
    # We're going to assume, for now, that they're already ordered by Pepper.
    # We'll see.
    return [ strandseq for (strandname, strandseq) in strandlist ]

def gettile(tset,tname):
    foundtiles = [x for x in tset['tiles'] if x['name']==tname]
    assert len(foundtiles)==1
    return foundtiles[0]
