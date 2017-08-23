import stickydesign as sd
import numpy as np
import re
import copy
from . import stickyends


def generate_xgrow_dict(ts,
                        perfect=False,
                        rotate=False,
                        energetics=None):
    
    # Combine ends and tile-specified adjacents
    newtiles = []
    newends = []
    doubleends = []
    doubles = []
    vdoubleends = []
    vdoubles = []
    ts = copy.deepcopy(ts)
    newtiles.append({'name': 'origami',
                     'edges': ['origami', 'origami', 'origami', 'origami'],
                     'stoic': 0,
                     'color': 'white'})

    atiles = [None]*16
    for tilename in ts['seed']['use_adapters']:
        try:
            tile = [ x for x in ts['seed']['adapters'] if x['name']==tilename ][0]
        except IndexError:
            raise Exception("Can't find {}".format(tilename))
        newtile = {}
        newtile['edges'] = [ 'origami' ] +  [ re.sub('/','_c',x) for x in tile['ends'] ] + [ 'origami' ]
        newtile['name'] = tile['name']
        newtile['stoic'] = 0
        newtile['color'] = 'white'
        atiles[tile['loc']-1] = newtile
    for tile in atiles:
        if tile:
            newtiles.append(tile)
        else:
            newtiles.append( { 'name': 'emptyadapt', 'edges': ['origami',0,0,'origami'], 'stoic': 0, 'color': 'white'} )


    if rotate:
        rotatedtiles = []
        for tile in ts['tiles']:
            if tile['type'] == 'tile_daoe_3up' or tile['type'] == 'tile_daoe_5up':
                newtile = copy.deepcopy(tile)
                newtile['name']+='_lrf'
                newtile['ends']=[tile['ends'][x] for x in (1,0,3,2)]
                rotatedtiles.append(newtile)
                newtile = copy.deepcopy(tile)
                newtile['name']+='_udf'
                newtile['type']='tile_daoe_'+{'5up':'3up','3up':'5up'}[tile['type'][-3:]]
                newtile['ends']=[tile['ends'][x] for x in (3,2,1,0)]
                rotatedtiles.append(newtile)
                newtile = copy.deepcopy(tile)
                newtile['name']+='_bf'
                newtile['type']='tile_daoe_'+{'5up':'3up','3up':'5up'}[tile['type'][-3:]]
                newtile['ends']=[tile['ends'][x] for x in (2,3,0,1)]
                rotatedtiles.append(newtile)
            elif tile['type'] == 'tile_daoe_doublehoriz_35up':
                newtile = copy.deepcopy(tile)
                newtile['name']+='_lrf'
                newtile['type']='tile_daoe_doublevert_53up'
                newtile['ends']=[tile['ends'][x] for x in (2,1,0,5,4,3)]
                rotatedtiles.append(newtile)
                newtile = copy.deepcopy(tile)
                newtile['name']+='_udf'
                newtile['type']='tile_daoe_doublevert_53up'
                newtile['ends']=[tile['ends'][x] for x in (5,4,3,2,1,0)]
                rotatedtiles.append(newtile)
                newtile = copy.deepcopy(tile)
                newtile['name']+='_bf'
                newtile['ends']=[tile['ends'][x] for x in (3,4,5,0,1,2)]
                rotatedtiles.append(newtile)
            elif tile['type'] == 'tile_daoe_doublevert_35up':
                newtile = copy.deepcopy(tile)
                newtile['name']+='_lrf'
                newtile['type']='tile_daoe_doublehoriz_53up'
                newtile['ends']=[tile['ends'][x] for x in (2,1,0,5,4,3)]
                rotatedtiles.append(newtile)
                newtile = copy.deepcopy(tile)
                newtile['name']+='_udf'
                newtile['type']='tile_daoe_doublehoriz_53up'
                newtile['ends']=[tile['ends'][x] for x in (5,4,3,2,1,0)]
                rotatedtiles.append(newtile)
                newtile = copy.deepcopy(tile)
                newtile['name']+='_bf'
                newtile['ends']=[tile['ends'][x] for x in (3,4,5,0,1,2)]
                rotatedtiles.append(newtile)

        ts['tiles'] += rotatedtiles
    
    for tile in ts['tiles']:
        if tile['type'] == 'tile_daoe_3up' or tile['type'] == 'tile_daoe_5up':
            newtile = {}
            newtile['edges'] = [ re.sub('/','_c',x) for x in tile['ends'] ]
            if 'name' in tile: newtile['name'] = tile['name']
            if 'conc' in tile: newtile['stoic'] = tile['conc']
            if 'color' in tile: newtile['color'] = tile['color']
            newtiles.append(newtile)

        if tile['type'] == 'tile_daoe_doublehoriz_35up' or tile['type'] == 'tile_daoe_doublehoriz_53up':
            newtile1 = {}
            newtile2 = {}
            newtile1['edges'] = [ re.sub('/','_c',x) for x in tile['ends'][0:1] ] \
                + [ tile['name']+'_db' ] \
                + [ re.sub('/','_c',x) for x in tile['ends'][4:] ]
            newtile2['edges'] = [ re.sub('/','_c',x) for x in tile['ends'][1:4] ] \
                + [ tile['name']+'_db' ]            
            newtile1['name'] = tile['name']+'_left'
            newtile2['name'] = tile['name']+'_right'
                        
            doubleends.append( tile['name']+'_db' )
            doubles.append( (newtile1['name'], newtile2['name']) )
            
            if 'conc' in tile: 
                newtile1['stoic'] = tile['conc']
                newtile2['stoic'] = tile['conc']
                
            if 'color' in tile: 
                newtile1['color'] = tile['color']
                newtile2['color'] = tile['color']
                
            newtiles.append(newtile1)
            newtiles.append(newtile2)
        if tile['type'] == 'tile_daoe_doublevert_35up' or tile['type'] == 'tile_daoe_doublevert_53up':
            newtile1 = {}
            newtile2 = {}
            newtile1['edges'] = [ re.sub('/','_c',x) for x in tile['ends'][0:2] ] \
                + [ tile['name']+'_db' ] \
                + [ re.sub('/','_c',x) for x in tile['ends'][5:] ]
            newtile2['edges'] = [ tile['name']+'_db' ] + [ re.sub('/','_c',x) for x in tile['ends'][2:5] ] 
            newtile1['name'] = tile['name']+'_top'
            newtile2['name'] = tile['name']+'_bottom'
                        
            vdoubleends.append( tile['name']+'_db' )
            vdoubles.append( (newtile1['name'], newtile2['name']) )
            
            if 'conc' in tile: 
                newtile1['stoic'] = tile['conc']
                newtile2['stoic'] = tile['conc']
                
            if 'color' in tile: 
                newtile1['color'] = tile['color']
                newtile2['color'] = tile['color']
                
            newtiles.append(newtile1)
            newtiles.append(newtile2)
    
    newends.append( { 'name': 'origami', 'strength': 100 } )
 
    for end in doubleends:
        newends.append( { 'name': end, 'strength': 10 } )
    for end in vdoubleends:
        newends.append( { 'name': end, 'strength': 10 } )

    gluelist = []
    if not perfect: 
        glueends = {'DT': [], 'TD': []}
        for end in ts['ends']:
            newends.append( { 'name': end['name'], 'strength': 0 } )
            newends.append( { 'name': end['name']+'_c', 'strength': 0 } )
            if (end['type'] == 'TD') or (end['type'] == 'DT'):
                glueends[end['type']].append((end['name'],end['fseq']))
        
        if energetics:
            ef = energetics
        else:
            ef = stickyends.DEFAULT_ENERGETICS
       
        eavg = {}
        for t in ['DT', 'TD']:
            names, fseqs = zip(*glueends[t])
            ea = sd.endarray(fseqs, t)
            eavg[t] = np.average( ef.matching_uniform( ea ) )
        eavg_combined = ( eavg['DT'] + eavg['TD'] ) / 2.0
        
        for t in ['DT','TD']:
            names, fseqs = zip(*glueends[t])
            allnames = names + tuple( x+'_c' for x in names )
            ea = sd.endarray(fseqs, t)
            ar = sd.energy_array_uniform(ea,ef) / eavg_combined
            for i1,n1 in enumerate(names):
                for i2,n2 in enumerate(allnames):
                    gluelist.append([n1,n2,max(float(ar[i1,i2]),0.0)])

         
    else:
        if 'ends' not in ts.keys():
            ts['ends']=[]
        endsinlist = set( e['name'] for e in ts['ends'] )
        endsintiles = set()
        for tile in ts['tiles']:
            endsintiles.update( re.sub('/','',e) for e in tile['ends'] if e != 'hp')
        for end in ts['ends'] + list({'name': e} for e in endsintiles):
            newends.append( { 'name': end['name'], 'strength': 0 } )
            newends.append( { 'name': end['name']+'_c', 'strength': 0 } )
            gluelist.append([end['name'],end['name']+'_c',1.0]) 
            

    
    newends.append( {'name': 'hp', 'strength': 0} )

    xga = {}
    xga['doubletiles'] = [ list(x) for x in doubles ]
    xga['vdoubletiles'] = [ list(x) for x in vdoubles ]
    xga.update( ts['xgrow_options'] )
    xga.update( ts['xgrow_options'] )
    #if not perfect:
    #    xga['gse_calc_avg'] = eavg_combined
     
        
    sts = { 'tiles': newtiles, 'bonds': newends, 'xgrowargs': xga, 'glues': gluelist }
    
    return sts
