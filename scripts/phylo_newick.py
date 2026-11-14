#!/usr/bin/env python3
"""Minimal Newick utilities for Figure 2 gene-tree processing."""

import re


class Node(object):
    def __init__(self, name="", length="", children=None):
        self.name = name
        self.length = length
        self.children = children or []


def tokenize(text):
    text = re.sub(r"\[[^\]]*\]", "", text.strip())
    return re.findall(r"\(|\)|,|:|;|'(?:[^']|'')*'|[^\s\(\),:;]+", text)


def parse_newick(text):
    tokens = tokenize(text)
    position = [0]

    def parse_label_length():
        name = ""
        length = ""
        if position[0] < len(tokens) and tokens[position[0]] not in [
            "(",
            ")",
            ",",
            ":",
            ";",
        ]:
            name = tokens[position[0]].strip("'")
            position[0] += 1
        if position[0] < len(tokens) and tokens[position[0]] == ":":
            position[0] += 1
            if position[0] >= len(tokens):
                raise ValueError("missing branch length")
            length = tokens[position[0]]
            position[0] += 1
        return name, length

    def parse_node():
        if position[0] < len(tokens) and tokens[position[0]] == "(":
            position[0] += 1
            children = [parse_node()]
            while position[0] < len(tokens) and tokens[position[0]] == ",":
                position[0] += 1
                children.append(parse_node())
            if position[0] >= len(tokens) or tokens[position[0]] != ")":
                raise ValueError("unbalanced Newick parentheses")
            position[0] += 1
            name, length = parse_label_length()
            return Node(name, length, children)
        name, length = parse_label_length()
        if not name:
            raise ValueError("unnamed Newick leaf")
        return Node(name, length, [])

    root = parse_node()
    if position[0] < len(tokens) and tokens[position[0]] == ";":
        position[0] += 1
    if position[0] != len(tokens):
        raise ValueError("unexpected Newick tokens")
    return root


def serialize(node):
    prefix = (
        "(" + ",".join(serialize(child) for child in node.children) + ")"
        if node.children
        else ""
    )
    suffix = node.name
    if node.length:
        suffix += ":" + node.length
    return prefix + suffix


def ufboot_support(label):
    if not label:
        return None
    raw = label.split("/")[-1]
    try:
        return float(raw)
    except ValueError:
        return None


def contract_low_support(node, threshold, is_root=True):
    contracted = []
    for child in node.children:
        child = contract_low_support(child, threshold, False)
        support = ufboot_support(child.name) if child.children else None
        if child.children and support is not None and support < threshold:
            contracted.extend(child.children)
        else:
            contracted.append(child)
    node.children = contracted
    return node


def graph(root):
    adjacency = {}
    leaves = {}
    counter = [0]

    def visit(node, parent=None):
        ident = counter[0]
        counter[0] += 1
        adjacency[ident] = []
        if parent is not None:
            adjacency[ident].append(parent)
            adjacency[parent].append(ident)
        if not node.children:
            leaves[node.name] = ident
        for child in node.children:
            visit(child, ident)

    visit(root)
    return adjacency, leaves


def distance(adjacency, start, target):
    queue = [(start, 0)]
    seen = set([start])
    for node, depth in queue:
        if node == target:
            return depth
        for neighbor in adjacency[node]:
            if neighbor not in seen:
                seen.add(neighbor)
                queue.append((neighbor, depth + 1))
    raise ValueError("disconnected tree")


def quartet_topology(root, taxa):
    adjacency, leaves = graph(root)
    if any(taxon not in leaves for taxon in taxa):
        return "missing"
    a, b, c, d = taxa
    scores = [
        (
            distance(adjacency, leaves[a], leaves[b])
            + distance(adjacency, leaves[c], leaves[d]),
            a + "+" + b + "|" + c + "+" + d,
        ),
        (
            distance(adjacency, leaves[a], leaves[c])
            + distance(adjacency, leaves[b], leaves[d]),
            a + "+" + c + "|" + b + "+" + d,
        ),
        (
            distance(adjacency, leaves[a], leaves[d])
            + distance(adjacency, leaves[b], leaves[c]),
            a + "+" + d + "|" + b + "+" + c,
        ),
    ]
    scores.sort(key=lambda item: item[0])
    if scores[0][0] == scores[1][0]:
        return "unresolved"
    return scores[0][1]
