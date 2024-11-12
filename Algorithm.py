import Circuit

#ralph_position = [7.101560294445499, 6.796879411109002]
ralph_position = [5, 10]
vertices = [[1, 2], [7, 14], [3, 5], [6, 9], [6.1, 5], [8, 4]]
#vertices = [[1, 2], [7, 14], [3, 5], [6, 9], [8, 5], [9, 10]]
lista = []


class TreeNode:
    def __init__(self, value):
        self.value = value
        self.children = []
        self.parent = []


def build_tree(n):

    root = TreeNode(0)
    remaining_elements = list(range(1, n + 1))

    print(remaining_elements)

    build_tree_recursive(root, remaining_elements, n)

    return root


def build_tree_recursive(node, remaining_elements, n):
    if not remaining_elements:
        return

    for value in remaining_elements:
        child = TreeNode(value)
        node.children.append(child)
        child.parent.append(node)

        if(value == n):
            lista.append([])
        if value < n:
            child_remaining = list(range(value + 1, n + 1))
            build_tree_recursive(child, child_remaining, n)



def print_tree(node, depth=0):
    if not node:
        return

    print("  " * depth + str(node.value))
    for child in node.children:
        print("child value")
        print(child.value)
        print("parent")
        print(child.parent[0].value)
        print_tree(child, depth + 1)


def create_parent_list(node, parent_list, value):
    if(node.value == 0):
        return parent_list
    if(node.value == value):
        parent_list.append(node.value)
    parent_list.append(node.parent[0].value)
    return create_parent_list(node.parent[0], parent_list, value)

def solve1(node, solution, slopes, lengths, parent_list, final_list):
    if not node:
        return

    for child in node.children:

        if (node.value == 0):

            result = Circuit.find_interception_interval_forward(vertices, ralph_position, node.value, child.value, slopes, lengths)


            if(result != []):
                final_list.append([child.value, node.value])
                solve1(child, result, slopes, lengths, parent_list, final_list)
        else:
            result = Circuit.find_next_interception_interval(solution, vertices, node.value, child.value, slopes, lengths)
            if(result != []):
                parent_list = []
                parent_list = create_parent_list(child, parent_list, child.value)
                final_list.append(parent_list)
                solve1(child, result, slopes, lengths, parent_list, final_list)


def solve():
    slopes = Circuit.calculate_slope(vertices)
    lengths = Circuit.calculate_lengths(vertices)
    max = 0
    path = []
    parent_list = []
    final_list = []
    final_results = []

    n = len(vertices) - 2
    tree_root = build_tree(n)
    solve1(tree_root, [], slopes, lengths, parent_list, final_list)

    for element in final_list:
        if(len(element) > max):
            max = len(element)
            path = element
    for i in range(len(path) - 1, 0, -1):
        if(path[i] == 0):
            result = Circuit.find_interception_interval_forward(vertices, ralph_position, path[i], path[i-1], slopes, lengths)
            final_results.append(result)
        else:
            result = Circuit.find_next_interception_interval(result, vertices, path[i], path[i-1], slopes, lengths)
            final_results.append(result)


    Circuit.plot_circuit(vertices, vertices[0], final_results)
