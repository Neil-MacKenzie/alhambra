import ruamel.yaml as yaml
from ruamel.yaml.comments import CommentedMap
from ruamel.yaml.representer import RoundTripRepresenter
from collections import Counter 
from .util import named_list, merge_endlists
from . import tiletypes
from functools import reduce
import operator
import copy


class TileSet(CommentedMap):
    def __init__(self, val={}):
        CommentedMap.__init__(self, val)
        if 'ends' in self.keys():
            self['ends'] = named_list(self['ends'])
        else:
            self['tiles'] = named_list()
        if 'tiles' in self.keys():
            self['tiles'] = named_list(self['tiles'])
        else:
            self['tiles'] = named_list()


    def check_consistent(self):
        # * END LIST The end list itself must be consistent.
        # ** Each end must be of understood type
        # ** Each end must have a valid sequence or no sequence
        # ** There must be no more than one instance of each name
        # ** WARN if there are ends with no namecounts
        # * TILE LIST
        # ** each tile must be of understood type (must parse)
        for tile in self['tiles']:
            parsed = tiletypes.tfactory.parse(tile)
            # ** the tile type edotparen must be consistent, if it has one
            if parsed.edotparen:
                tiletypes.check_edotparen_consistency(parsed.edotparen)
            else:
                log.warning("tile type {} has no edotparen".format(tile['type']))
            # ** each tile must have no sequence, or a valid sequence
            if 'fullseqs' in tile.keys():
                parsed.check_sequence()
            # ** each tile must have the right number of ends
            if 'ends' in tile.keys():
                assert len(parsed._endtypes) == len(tile['ends'])
        # ** ends in the tile list must be consistent (must merge)
        endsfromtiles = tiletypes.endlist_from_tilelist(self['tiles'])
        # ** there must be no more than one tile with each name
        self['tiles'].check_consistent()
        # ** WARN if any end that appears does not have a complement used or vice versa
        # ** WARN if there are tiles with no name
        # * TILE + END
        # ** The tile and end lists must merge validly
        # (checks sequences, adjacents, types, complements)
        merge_endlists(self['ends'], endsfromtiles)
        
        # ** WARN if tilelist has end references not in ends
        # ** WARN if merge is not equal to the endlist
        # ** WARN if endlist has ends not used in tilelist
        # * ADAPTERS / SEEDS
        # ** seeds must be of understood type
        # ** adapter locations must be valid
        # ** each adapter must have no sequence or a consistent sequence
        # *** the RH strand must match the associated tile
        # *** the ends in the sequence must match the ends in the endlist
        # *** the LH sequence must be validly binding to both RH and origami
        # ** each adapter must have valid definition, which means for us:
        # *** if both tile mimic and ends are specified, they must match
    
    def summary(self):
        self.check_consistent()
        info = {'ntiles': len(self['tiles']),
                'nends':  len(self['ends']),
                'ntends': len(tiletypes.endlist_from_tilelist(self['tiles'])),
                'tns':    " ".join(x['name'] for x in self['tiles'] if 'name' in x.keys()),
                'ens':    " ".join(x['name'] for x in self['ends'] if 'name' in x.keys()),
                'name':   " {}".format(self['info']['name']) if \
                                       ('info' in self.keys() and \
                                        'name' in self['info'].keys()) else ""}
        tun = sum( 1 for x in self['tiles'] if 'name' not in x.keys() )
        if tun > 0:
            info['tns'] += " ({} unnamed)".format(tun)
        eun = sum( 1 for x in self['ends'] if 'name' not in x.keys() )
        if eun > 0:
            info['ens'] += " ({} unnamed)".format(eun)
        return "TileSet{name}: {ntiles} tiles, {nends} ends, {ntends} ends in tiles.\nTiles: {tns}\nEnds:  {ens}".format(**info)

    def __str__( self ):
        return self.summary()

    def copy(self):
        return copy.deepcopy(self)

    def __deepcopy__(self, memo):
        # type: (Any) -> Any
        res = self.__class__()
        memo[id(self)] = res
        for k in self:
            res[k] = copy.deepcopy(self[k])
            self.copy_attributes(res, deep=True)
        return res

    def dump(self, stream):
        return yaml.round_trip_dump(self, stream)
    
RoundTripRepresenter.add_representer(TileSet,
                                     RoundTripRepresenter.represent_dict)


def load_tileset_dict(*args, **kwargs):
    return TileSet(yaml.round_trip_load(*args, **kwargs))
