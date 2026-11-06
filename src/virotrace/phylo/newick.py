"""A small Newick reader, enough for IQ-TREE output (tip names, branch lengths, support labels)."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..errors import ContractError


@dataclass
class Tree:
    parent: list[int]  # -1 for the root
    length: list[float]  # length of the edge above each node (0 for the root)
    children: list[list[int]]
    name: list[str]  # tip names; "" for internal nodes
    label: list[str]  # support label on the edge above each internal node ("" if none)
    tips: list[int] = field(default_factory=list)  # node ids of tips, in preorder

    @property
    def n_nodes(self) -> int:
        return len(self.parent)


def parse_newick(text: str) -> Tree:
    """Parse one Newick tree. Internal labels (support values) are ignored."""
    s = "".join(text.split())
    if not s.endswith(";"):
        raise ContractError("newick: missing terminating ';'")
    s = s[:-1]
    parent, length, children, name, label = [-1], [0.0], [[]], [""], [""]
    stack: list[int] = []
    cur = 0
    i = 0

    def new_child(p: int) -> int:
        parent.append(p)
        length.append(0.0)
        children.append([])
        name.append("")
        label.append("")
        children[p].append(len(parent) - 1)
        return len(parent) - 1

    while i < len(s):
        ch = s[i]
        if ch == "(":
            stack.append(cur)
            cur = new_child(cur)
            i += 1
        elif ch == ",":
            cur = new_child(stack[-1])
            i += 1
        elif ch == ")":
            cur = stack.pop()
            i += 1
        elif ch == ":":
            j = i + 1
            while j < len(s) and s[j] not in ",()":
                j += 1
            try:
                length[cur] = float(s[i + 1 : j])
            except ValueError as exc:
                raise ContractError(f"newick: bad branch length {s[i + 1 : j]!r}") from exc
            i = j
        else:
            j = i
            while j < len(s) and s[j] not in ",():":
                j += 1
            if not children[cur]:
                name[cur] = s[i:j]
            else:  # a label after ')' is a support value
                label[cur] = s[i:j]
            i = j
    if stack:
        raise ContractError("newick: unbalanced parentheses")

    # renumber in preorder so each subtree's tips are contiguous
    order: list[int] = []
    todo = [0]
    while todo:
        node = todo.pop()
        order.append(node)
        todo.extend(reversed(children[node]))
    new_id = {old: k for k, old in enumerate(order)}
    tree = Tree(
        parent=[-1 if parent[o] < 0 else new_id[parent[o]] for o in order],
        length=[length[o] for o in order],
        children=[[new_id[c] for c in children[o]] for o in order],
        name=[name[o] for o in order],
        label=[label[o] for o in order],
    )
    tree.tips = [k for k in range(tree.n_nodes) if not tree.children[k]]
    if len(tree.tips) < 3:
        raise ContractError("newick: fewer than three tips")
    return tree


def tip_distances(tree: Tree) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(D, lo, hi)``: D[node, tip] is the path length; tips of node k are lo[k]:hi[k]."""
    n, t = tree.n_nodes, len(tree.tips)
    tip_index = {node: k for k, node in enumerate(tree.tips)}
    lo, hi = np.zeros(n, dtype=int), np.zeros(n, dtype=int)
    for node in range(n - 1, -1, -1):  # reverse preorder: children before parents
        if node in tip_index:
            lo[node], hi[node] = tip_index[node], tip_index[node] + 1
        else:
            kids = tree.children[node]
            lo[node], hi[node] = lo[kids[0]], hi[kids[-1]]
    depth = np.zeros(n)
    for node in range(1, n):
        depth[node] = depth[tree.parent[node]] + tree.length[node]
    dist = np.zeros((n, t))
    dist[0] = depth[tree.tips]
    for node in range(1, n):
        row = dist[tree.parent[node]] + tree.length[node]
        row[lo[node] : hi[node]] -= 2 * tree.length[node]
        dist[node] = row
    return dist, lo, hi


def midpoint_edge(tree: Tree) -> tuple[int, float]:
    """``(node, x)``: the midpoint of the longest tip-to-tip path lies ``x`` above ``node``
    on the edge between ``node`` and its parent."""
    dist, _, _ = tip_distances(tree)
    tip_rows = dist[tree.tips]
    a, b = np.unravel_index(int(np.argmax(tip_rows)), tip_rows.shape)
    half = float(tip_rows[a, b]) / 2

    def ancestors(node: int) -> list[int]:
        out = []
        while node != -1:
            out.append(node)
            node = tree.parent[node]
        return out

    up_a, up_b = ancestors(tree.tips[a]), ancestors(tree.tips[b])
    shared = set(up_a) & set(up_b)
    climb_a = [n for n in up_a if n not in shared]
    climb_b = [n for n in up_b if n not in shared]
    acc = 0.0
    for node in climb_a:  # tip a upwards
        if acc + tree.length[node] >= half:
            return node, half - acc
        acc += tree.length[node]
    for node in reversed(climb_b):  # then down to tip b
        if acc + tree.length[node] >= half:
            return node, tree.length[node] - (half - acc)
        acc += tree.length[node]
    return climb_a[-1] if climb_a else climb_b[0], 0.0


def to_newick(tree: Tree, root_edge: tuple[int, float] | None = None) -> str:
    """Write the tree, optionally re-rooted ``x`` above ``node`` on its parent edge."""
    nbrs: dict[int, list[tuple[int, float, str]]] = {k: [] for k in range(tree.n_nodes)}
    for c in range(1, tree.n_nodes):
        p = tree.parent[c]
        nbrs[p].append((c, tree.length[c], tree.label[c]))
        nbrs[c].append((p, tree.length[c], tree.label[c]))
    if root_edge is None:
        stack = [(0, -1, 0.0, "")]
    else:
        u, x = root_edge
        p = tree.parent[u]
        # a virtual root joins u (x below the cut) and p (the rest of the edge)
        stack = [(p, u, tree.length[u] - x, ""), (u, p, x, tree.label[u])]
    order: list[tuple[int, int, float, str]] = []  # (node, came_from, edge_length, label)
    while stack:
        node, frm, ln, lab = stack.pop()
        order.append((node, frm, ln, lab))
        stack.extend(reversed([(n, node, e, lb) for n, e, lb in nbrs[node] if n != frm]))
    text: dict[int, str] = {}
    kids: dict[int, list[int]] = {}
    for node, frm, _, _ in order:
        kids.setdefault(frm, []).append(node)  # children of the virtual root hang off from=-1 ids
    for node, frm, ln, lab in reversed(order):
        children = [n for n in kids.get(node, []) if n != frm]
        body = tree.name[node] if not tree.children[node] else ""
        if children:
            body = "(" + ",".join(text[c] for c in children) + ")" + lab
        text[node] = body + (f":{ln:.8g}" if frm != -1 or ln else "")
    if root_edge is None:
        return text[0].removesuffix(":0") + ";"
    return "(" + text[root_edge[0]] + "," + text[tree.parent[root_edge[0]]] + ");"
