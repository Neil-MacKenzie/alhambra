from ruamel.yaml.comments import CommentedMap
from .util import NamedList
from .ends import EndList
from .tilestructures import TileStructure, getstructure


class Tile(CommentedMap):
    def __init__(self, val={}):
        CommentedMap.__init__(self, val)

        if 'type' in self.keys():
            self['structure'] = self['type']
            del(self['type'])

        self._structure = getstructure(self.get('structure', None),
                                       extra=self.get('extra', None))

    def structure():
        doc = """Doc string"""
        
        def fget(self):
            return self._structure
    
        def fset(self, value):
            if isinstance(value, TileStructure):
                self._structure = value
                self['structure'] = value.name
            else:
                self._structure = getstructure(self.get('structure', None),
                                               extra=self.get('extra', None))
                self['structure'] = self._structure.name
    
        def fdel(self):
            self._structure = None
            del self['structure']
        
        return locals()
    structure = property(**structure())

    @property
    def strands(self):
        return self.get('fullseqs', None)

    def check_strands(self):
        self.structure.check_strands(self)

    def sequence_diagram(self):
        return self.structure.sequence_diagram(self)

    def abstract_diagram(self, tileset=None):
        return self.structure.abstract_diagram(self, tileset=None)
        
    def name():
        doc = """Doc string"""
        def fget(self):
            return self['name']
    
        def fset(self, value):
            self['name'] = value
    
        return locals()
    name = property(**name())
    
    def ends():
        doc = """Doc string"""
        
        def fget(self):
            return self['ends']
    
        def fset(self, value):
            # Make sure there are the right number of ends.
            # FIXME: do a better job checking, raise a better error
            if self.structure:
                assert len(self.structure._endlocs) == len(value)
            self._ends = value
    
        def fdel(self):
            del self['ends']
        return locals()
    ends = property(**ends())

    @property
    def endlist(self):
        return self.structure.get_endlist(self)

    @property
    def orderableseqs(self):
        return self.structure.orderableseqs(self)
    
    def check_consistent(self):
        if self.structure:
            self.structure.check_consistent()
            if self.ends:
                assert(self.structure.numends == len(self.ends))
            if self.strands:
                self.check_strands()

    def check_sequences(self):
        self.structure.check_strands(self.strands)


class TileList(NamedList):
    def __init__(self, val=[]):
        NamedList.__init__(self, val)
        for i, tile in enumerate(self):
            self[i] = Tile(tile)

    def check_consistent(self):
        NamedList.check_consistent(self)
        for tile in self:
            tile.check_consistent()

    def endlist(self, fail_immediate=True):
        """\
    Given a NamedList of tiles (or just a list, for now), extract the sticky ends 
    from each tile, and merge these (using merge_ends) into a NamedList of sticky
    ends.

    Parameters
    ----------

    tilelist: either a list of tile description dicts, or a list of tile instances.
        If 

    fail_immediate: (default True) if True, immediately fail on a failure,
        with ValueError( tilename, exception ) from exception  If False, collect 
        exceptions, then raise ValueError( "message", [(tilename, exception)] , 
        output ).

        """
        endlist = EndList()
        errors = []

        for tile in self:
            try:
                endlist = endlist.merge(tile.endlist, in_place=True)
            except BaseException as e:
                if fail_immediate:
                    raise ValueError(tile.name) from e
                else:
                    errors.append((tile.name, e))

        if errors:
            raise ValueError("End list generation failed on:", errors, endlist)

        return endlist
