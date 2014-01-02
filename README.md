# rbtree - a C library implementing red/black trees

  This code implements a clean API for working with
mutable red/black trees.  It's unique in working with the
user's own data structure, which only requires space
to store L, R child pointers and one bit to hold
black (off) or red (on) node coloring.

  No parent pointers are required, as the tree traversals
used by the library are written in a creative, recursive
way to store that information on the call stack.

# Example:

  There is an example implementation using structs holding
integers named 'test' in the top-level dir, along
with some great display code using dotty from [www.graphviz.org](GraphViz).

