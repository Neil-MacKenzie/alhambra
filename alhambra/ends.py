import warnings
from ruamel.yaml.comments import CommentedMap
from .util import NamedList, merge_ends
import copy
from peppercompiler.DNA_classes import wc

class End(CommentedMap):
    def __str__(self):
        if self.fseq:
            if self.etype == 'DT':
                s = self.seq[0]+'-'+self.seq[1:]
                c = self.comp[0]+'-'+self.comp[1:]
            elif self.etype == 'TD':
                s = self.seq[:-1]+'-'+self.seq[-1]
                c = self.comp[:-1]+'-'+self.comp[-1]
            return "<end {} ({}{}): {} | {}>".format(
                self['name'], self['type'], len(self.seq), s, c)
        else:
            return "<end {} ({})>".format(
                self['name'], self.get('type', '?'))

    def fseq():
        def fget(self):
            return self.get('fseq', None)
    
        def fset(self, value):
            if self.get('fseq', None) and len(value) != len(self['fseq']):
                warnings.warn("Changing end length")
            self['fseq'] = value
    
        def fdel(self):
            del self['fseq']
        return locals()
    fseq = property(**fseq())

    def name():

        def fget(self):
            return self['name']
    
        def fset(self, value):
            self['name'] = value
    
        return locals()
    name = property(**name())
    
    @property
    def seq(self):
        if not self.fseq:
            return None
        if self.etype == 'TD':
            return self.fseq[1:]
        elif self.etype == 'DT':
            return self.fseq[:-1]

    @property
    def comp(self):
        if not self.fseq:
            return None
        if self.etype == 'TD':
            return wc(self.fseq[:-1].upper()).lower()
        elif self.etype == 'DT':
            return wc(self.fseq[1:].upper()).lower()

    @property
    def etype(self):
        return self['type']


class EndList(NamedList):
    def __init__(self, val=[]):
        NamedList.__init__(self, val)
        for i, end in enumerate(self):
            self[i] = End(end)

    def merge(endlist1, endlist2, fail_immediate=False,
              in_place=False):
        """\
        Given end lists `endlist1` and `endlist2`, merge the two lists, using
        `merge_ends` to merge any named ends that are present in both.

        Parameters
        ----------

        endlist1: NamedList of sticky ends.

        endlist2: NamedList OR just a list of sticky ends.  If it just a list,
            which may have multiple copies of the same named sticky end, each
            sticky end in order is merged.

        fail_immediate: (default False) if True, fail immediately on a
        failed merge, passing through the ValueError from merge_ends.
        If False, finish merging the two lists, then raise a
        ValueError listing *all* ends that failed to merge.

        in_place: (default False) if True, do merging in place in endlist1.
        Note the merged and added ends from endlist2 will be copies regardless.

        output: a merged NamedList of sticky ends.


        Exceptions
        ----------

        ValueError: In the event of a failed merge, when named ends
        cannot be merged.  If fail_immediate is True, then this is passed
        through from merge_ends.  If fail_immediate is False, then the
        ValueError has the following args: ("message", [exceptions],
        [failed_name,failed_end1,failed_end2) ...],out)

        """
        # Check consistency of each NamedList
        endlist1.check_consistent()
        try:
            endlist2.check_consistent()
        except AttributeError:
            pass

        if not in_place:
            out = copy.deepcopy(endlist1)
        else:
            out = endlist1

        exceptions = []
        errors = []

        for end in endlist2:
            if end['name'] in out.keys():
                try:
                    out[end['name']] = merge_ends(out[end['name']], end)
                except ValueError as e:
                    if fail_immediate:
                        raise e
                    else:
                        exceptions.append(e)
                        errors.append((end['name'], out[end['name']], end))
            else:
                out.append(copy.deepcopy(end))

        if errors:
            errorstring = " ".join([e[0] for e in errors])
            raise ValueError("Errors merging {}".format(errorstring),
                             exceptions, errors, out)

        return out

