#!/usr/bin/env python3

# take a sequence, contig, scaffold of DNA, and split into n windows of size k, overlapping or not.
# stolen directly from https://stackoverflow.com/questions/6822725/rolling-or-sliding-window-iterator

from collections import deque
from itertools import islice

def slidingWindow(iterable, size=2, step=1, fillvalue=None):
    if size < 0 or step < 1:
        raise ValueError
    it = iter(iterable)
    q = deque(islice(it, size), maxlen=size)
    if not q:
        return  # empty iterable or size == 0
    q.extend(fillvalue for _ in range(size - len(q)))  # pad to size
    while True:
        yield iter(q)  # iter() to avoid accidental outside modifications
        try:
            q.append(next(it))
        except StopIteration: # Python 3.5 pep 479 support
            return
        q.extend(next(it, fillvalue) for _ in range(step - 1))
