import copy
import time

from ruamel.yaml.comments import CommentedSeq
from ruamel.yaml.representer import RoundTripRepresenter

from . import seq


class NamedList(CommentedSeq):
    """A class for a list of dicts, where some dicts have a 'name' item
which should be unique, but others might not.  Indexing works with
either number or name.  Note that updating dicts may make the list
inconsistent.

    """

    def __init__(self, x=[]):
        CommentedSeq.__init__(self, x)

    def __getitem__(self, i):
        if isinstance(i, str):
            r = [x for x in self if x.get('name', None) == i]
            if len(r) > 1:
                raise KeyError(
                    "There are {} elements named {}.".format(len(r), i))
            elif len(r) == 0:
                raise KeyError("No element named {} found.".format(i))
            else:
                return r[0]
        else:
            return CommentedSeq.__getitem__(self, i)

    def check_consistent(self):
        """Checks that each name appears only once.  On failure, returns a
ValueError with (message, {failed_name: count}).  Otherwise, return
with no output.

        """
        names = [v['name'] for v in self if 'name' in v.keys()]

        # if we have *no* names, the tests will fail, but we are obviously
        # consistent, so:
        if not names:
            return

        from collections import Counter
        namecounts = Counter(names)

        if max(namecounts.values()) > 1:
            badcounts = {n: v for n, v in namecounts.items() if v > 1}
            raise ValueError("Inconsistent NamedList.", badcounts)

    def __setitem__(self, i, v):
        if isinstance(i, str):
            r = [(ii, x) for ii, x in enumerate(self)
                 if x.get('name', None) == i]
            if len(r) > 1:
                raise KeyError(
                    "There are {} elements named {}.".format(len(r), i))
            elif len(r) == 0:
                self.append(v)
            else:
                CommentedSeq.__setitem__(self, r[0][0], v)
        else:
            CommentedSeq.__setitem__(self, i, v)

    def keys(self):
        return [x['name'] for x in self if 'name' in x.keys()]

    def __deepcopy__(self, memo):
        # FIXME: this is here waiting on ruamel.yaml bugfix.
        res = self.__class__()
        memo[id(self)] = res
        for k in self:
            res.append(copy.deepcopy(k))
            self.copy_attributes(res, deep=True)
        return res


RoundTripRepresenter.add_representer(NamedList,
                                     RoundTripRepresenter.represent_list)


class ProgressLogger(object):
    def __init__(self, logger, N, seconds_interval=60):
        self.logger = logger
        stime = time.perf_counter()
        self.stime = stime
        self.ltime = stime
        self.li = 0
        self
        self.seconds_interval = seconds_interval
        self.N = N
        self.logger.info("starting {} tasks".format(self.N))

    def update(self, i):
        ctime = time.perf_counter()
        if ctime - self.ltime > self.seconds_interval:
            self.logger.info(
                "finished {}/{}, {} s elapsed, {} s est remaining".format(
                    i, self.N,
                    int(ctime - self.stime),
                    int((self.N - i) * (ctime - self.stime) / i)))
            self.ltime = ctime
            self.li = i
