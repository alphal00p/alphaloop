#!/bin/sh
""":" .

exec python "$0" "$@"
"""

from __future__ import print_function

import argparse
import copy
import re
import sys

__doc__ = """\
FORM log profiler.

Example
-------
$ formprof.py myformprogram.log
$ formprof.py -u myformprogram.log
$ formprof.py -t myformprogram.log
$ formprof.py -m myformprogram.log
$ formprof.py -e myformprogram.log

Python versions
---------------
2.7, 3.2, 3.3, 3.4, 3.5

"""

if sys.version_info[0] >= 3:
    string_types = str,
else:
    string_types = basestring,


class Stat(object):
    """Statistics object."""

    __slots__ = (
        'name',             # module name (str)
        'expr',             # expression name (str)
        'start',            # starting time (float)
        'end',              # end time (float)
        'elapsed',          # elapsed time (float)
        'count',            # number of occurrence (int)
        'generated_terms',  # number of generated terms (int)
        'terms_in_output',  # number of terms in output (int)
        'bytes_used',       # number of bytes used in output (int)
        'parent',           # parent node (Stat)
        'children',         # child nodes (List[Stat])
    )

    def is_leaf(self):
        """Return True if the object is a leaf."""
        # A node has one or more module statistics information, so `start` and
        # `end` are useless.
        return hasattr(self, 'start')

    def add_child(self, other):
        """Add a child item."""
        if not isinstance(other, Stat):
            raise TypeError('the object is not "Stat"')

        if getattr(other, 'parent', None):
            raise ValueError('the object already has a parent.')

        n = getattr(self, 'parent', None)
        while n:
            if n == other:
                raise ValueError('tried to add an ancestor as a child')
            n = n.parent

        # NOTE: objects in a tree must have the attributes `parent`,
        #       `children`, `elapsed` and `count`, and most likely have `name`.

        if not hasattr(self, 'children'):
            if self.is_leaf():
                # Convert this object from a leaf to a node.
                s = copy.copy(self)
                self.children = []
                self.children.append(s)
                s.parent = self
                s.children = []
                for a in self.__slots__:
                    if a not in ('name', 'elapsed', 'count', 'parent',
                                 'children'):
                        if hasattr(self, a):
                            delattr(self, a)
            else:
                self.children = []

        if not hasattr(self, 'parent'):
            self.parent = None

        self.children.append(other)
        other.parent = self
        if not hasattr(other, 'children'):
            other.children = []

        if not hasattr(self, 'elapsed'):
            self.elapsed = 0.0
        if not hasattr(self, 'count'):
            self.count = 0
        if not hasattr(other, 'elapsed'):
            other.elapsed = 0.0
        if not hasattr(other, 'count'):
            other.count = 0

        # Add the elapsed time and count to this node and ancestors.
        n = self
        while n:
            n.elapsed += other.elapsed
            n.count += other.count
            n = n.parent

    def __str__(self):
        """Return the string representation."""
        items = []
        for a in self.__slots__:
            if hasattr(self, a):
                items.append('{0}={1}'.format(a, getattr(self, a)))
        return '[Stat{0}{1}]'.format(
            ' ' if items else '',
            ', '.join(items),
        )


def analyze_logfile(file):
    """Generator of Stat objects from a file."""
    if isinstance(file, string_types):
        with open(file, 'r') as f:
            for s in analyze_logfile(f):
                yield s
            return

    SEARCH_TIME = 0  # noqa: N806
    SKIP_TIME = 1    # noqa: N806
    SEARCH_EXPR = 2  # noqa: N806
    SEARCH_NAME = 3  # noqa: N806

    state = SEARCH_TIME
    name = None
    expr = None
    time = 0.0
    old_time = 0.0
    generated_terms = 0
    terms_in_output = 0
    bytes_used = 0

    for line in file:
        line = line.rstrip()

        if state == SEARCH_TIME:
            m = re.match('\s*(?:Thread|Process) \d+ reporting', line)
            if m:
                state = SKIP_TIME
                continue

            m = re.match('W?Time =\s*([0-9.]+) sec \s*'
                         'Generated terms =\s*(\d+)', line)
            if m:
                time = float(m.group(1))
                generated_terms = int(m.group(2))
                state = SEARCH_EXPR
        elif state == SKIP_TIME:
            m = re.match('W?Time =\s*[0-9.]+ sec \s*'
                         'Generated terms =\s*\d+', line)
            if m:
                state = SEARCH_TIME
        elif state == SEARCH_EXPR:
            m = re.search('[0-9]+ Terms (?:left|active)', line)
            if m:
                state = SEARCH_TIME

            m = re.match('\s*(\S+?)\s* Terms in output =\s*(\d+)', line)
            if m:
                expr = m.group(1)
                terms_in_output = int(m.group(2))
                state = SEARCH_NAME
        elif state == SEARCH_NAME:
            m = re.match('(.+?)Bytes used \s*=\s*(\d+)', line)
            if m:
                name = m.group(1).strip()
                if not name:
                    name = '<unnamed>'
                bytes_used = int(m.group(2))
                state = SEARCH_TIME

                stat = Stat()
                stat.name = name
                stat.expr = expr
                stat.start = old_time
                stat.end = time
                stat.elapsed = time - old_time
                stat.count = 1
                stat.generated_terms = generated_terms
                stat.terms_in_output = terms_in_output
                stat.bytes_used = bytes_used
                old_time = time
                yield stat


def print_normal(stats, sort):
    """Print statistics in the normal mode."""
    total_time = stats[-1].end

    # Sort.
    if sort:
        stats.sort(key=lambda s: -s.elapsed)

    # Stringification.
    stats = [(
        s.name,
        s.expr,
        '{0:.2f}'.format(s.elapsed),
        '{0:.2%}'.format(s.elapsed / total_time),
        '{0:.2f}'.format(s.start),
        '{0:.2f}'.format(s.end),
        str(s.generated_terms),
        str(s.terms_in_output),
        str(s.bytes_used),
    ) for s in stats]

    # Construct the format.
    columns = [
        'module  ',
        'expr    ',
        'time',
        ' ',
        'start',
        'end',
        'genterms',
        'outterms',
        'bytes',
    ]

    column_widths = [
        max(max(len(s[i]) for s in stats), len(columns[i]))
        for i in range(len(columns))
    ]

    fmt = (
        '{{0:<{0}}}  {{1:<{1}}}  {{2:>{2}}}  {{3:>{3}}}  {{4:>{4}}}  '
        '{{5:>{5}}}  {{6:>{6}}}  {{7:>{7}}}  {{8:>{8}}}'
    ).format(*column_widths)

    # Print the result.
    print(fmt.format(*columns))
    for s in stats:
        print(fmt.format(*s))


def print_module(stats, sort):
    """Print statistics combined for each module."""
    total_time = stats[-1].end

    # Combine.
    new_stats = {}  # str -> Stat
    for s in stats:
        if s.name not in new_stats:
            t = Stat()
            t.name = s.name
            t.count = 1
            t.elapsed = s.elapsed
            new_stats[s.name] = t
        else:
            t = new_stats[s.name]
            t.count += 1
            t.elapsed += s.elapsed
    stats = list(new_stats.values())

    # Sort.
    if sort:
        stats.sort(key=lambda s: -s.elapsed)

    # Stringification.
    stats = [(
        s.name,
        str(s.count),
        '{0:.2f}'.format(s.elapsed),
        '{0:.2%}'.format(s.elapsed / total_time),
    ) for s in stats]

    # Construct the format.
    columns = [
        'module  ',
        'count',
        'time',
        ' ',
    ]

    column_widths = [
        max(max(len(s[i]) for s in stats), len(columns[i]))
        for i in range(len(columns))
    ]

    fmt = (
        '{{0:<{0}}}  {{1:>{1}}}  {{2:>{2}}}  {{3:>{3}}}'
    ).format(*column_widths)

    # Print the result.
    print(fmt.format(*columns).rstrip())
    for s in stats:
        print(fmt.format(*s))


def print_expr(stats, sort):
    """Print statistics combined for each expression."""
    total_time = stats[-1].end

    # Combine.
    new_stats = {}  # str -> Stat
    for s in stats:
        if s.expr not in new_stats:
            t = Stat()
            t.expr = s.expr
            t.count = 1
            t.elapsed = s.elapsed
            new_stats[s.expr] = t
        else:
            t = new_stats[s.expr]
            t.count += 1
            t.elapsed += s.elapsed
    stats = list(new_stats.values())

    # Sort.
    if sort:
        stats.sort(key=lambda s: -s.elapsed)

    # Stringification.
    stats = [(
        s.expr,
        str(s.count),
        '{0:.2f}'.format(s.elapsed),
        '{0:.2%}'.format(s.elapsed / total_time),
    ) for s in stats]

    # Construct the format.
    columns = [
        'expr    ',
        'count',
        'time',
        ' ',
    ]

    column_widths = [
        max(max(len(s[i]) for s in stats), len(columns[i]))
        for i in range(len(columns))
    ]

    fmt = (
        '{{0:<{0}}}  {{1:>{1}}}  {{2:>{2}}}  {{3:>{3}}}'
    ).format(*column_widths)

    # Print the result.
    print(fmt.format(*columns).rstrip())
    for s in stats:
        print(fmt.format(*s))


def print_tree(stats, sort):
    """Print statistics combined for each module (tree-like)."""
    total_time = stats[-1].end

    # Construct the tree. Each node is determined from the module name,
    # separated by '-'. For example,
    #   MyModule-1
    #   MyModule-2-a
    #   MyModule-2-b
    #   MyModule-3
    #   MyModule-3-a
    # are converted as
    #   MyModule*
    #   +- MyModule-1
    #   +- MyModule-2*
    #   |  +- MyModule-2-a
    #   |  +- MyModule-2-b
    #   +- MyModule-3*
    #      +- MyModule-3
    #      +- MyModule-3-a
    root = Stat()
    root.name = '*'
    nodes = {}  # str -> Stat
    for s in stats:
        names = s.name.split('-')
        names = ['-'.join(names[:i + 1]) for i in range(len(names))]
        n = root
        for name in names:
            if name not in nodes:
                new_node = Stat()
                new_node.name = name + '*'
                nodes[name] = new_node
                n.add_child(new_node)
            n = nodes[name]
        n.add_child(s)

    # Sort.
    def sort_tree(node):
        node.children.sort(key=lambda s: -s.elapsed)
        for n in node.children:
            sort_tree(n)

    if sort:
        sort_tree(root)

    # Stringification.
    stats = []

    def walk_tree(node, indent_str, last):
        if len(node.children) == 1:
            walk_tree(node.children[0], indent_str, last)
            return

        s = node
        ss = indent_str + ('+- ' if node != root else '') + s.name

        if node.is_leaf():
            stats.append([
                ss,
                s.expr,
                '',
                '{0:.2f}'.format(s.elapsed),
                '{0:.2%}'.format(s.elapsed / total_time),
                '{0:.2f}'.format(s.start),
                '{0:.2f}'.format(s.end),
                str(s.generated_terms),
                str(s.terms_in_output),
                str(s.bytes_used),
            ])
        else:
            stats.append([
                ss,
                '',
                str(s.count),
                '{0:.2f}'.format(s.elapsed),
                '{0:.2%}'.format(s.elapsed / total_time),
                '',
                '',
                '',
                '',
                ''
            ])

        if node != root:
            if (last):
                indent_str += '  '
            else:
                indent_str += '| '

        for n in node.children:
            walk_tree(n, indent_str, n == node.children[-1])

    walk_tree(root, '', True)

    # Construct the format.
    columns = [
        'module  ',
        'expr    ',
        'count',
        'time',
        ' ',
        'start',
        'end',
        'genterms',
        'outterms',
        'bytes',
    ]

    column_widths = [
        max(max(len(s[i]) for s in stats), len(columns[i]))
        for i in range(len(columns))
    ]

    fmt = (
        '{{0:<{0}}}  {{1:<{1}}}  {{2:>{2}}}  {{3:>{3}}}  {{4:>{4}}}  '
        '{{5:>{5}}}  {{6:>{6}}}  {{7:>{7}}}  {{8:>{8}}}  {{9:>{9}}}'
    ).format(*column_widths)

    # Print the result.
    print(fmt.format(*columns))
    for s in stats:
        print(fmt.format(*s).rstrip())


def main():
    """Entry point."""
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE, SIG_DFL)

    # Parse the command line arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument('logfile',
                        type=str,
                        metavar='LOGFILE',
                        help='log file to be analyzed')
    parser.add_argument('-t',
                        '--tree',
                        action='store_const',
                        const='tree',
                        dest='mode',
                        help='combine statistics for each module (tree-like)')
    parser.add_argument('-m',
                        '--module',
                        action='store_const',
                        const='module',
                        dest='mode',
                        help='combine statistics for each module')
    parser.add_argument('-e',
                        '--expr',
                        action='store_const',
                        const='expr',
                        dest='mode',
                        help='combine statistics for each expression')
    parser.add_argument('-s',
                        '--sort',
                        action='store_const',
                        const=True,
                        default=True,
                        help='sort by time (default)')
    parser.add_argument('-u',
                        '--unsort',
                        action='store_const',
                        const=False,
                        dest='sort',
                        help='do not sort')
    args = parser.parse_args()

    # Parse the log file.
    stats = list(analyze_logfile(args.logfile))
    if not stats:
        print('empty log', file=sys.stderr)
        return

    # Print statistics.
    if args.mode == 'tree':
        print_tree(stats, args.sort)
    elif args.mode == 'module':
        print_module(stats, args.sort)
    elif args.mode == 'expr':
        print_expr(stats, args.sort)
    else:
        print_normal(stats, args.sort)


if __name__ == '__main__':
    main()
