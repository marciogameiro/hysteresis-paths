### HysteresisPathsBlowup.py
### MIT LICENSE 2025 Marcio Gameiro

import DSGRN
import DSGRN_utils
import itertools
import functools
import sys

class NetworkAnalyzer:
    def __init__(self, network, P_gene, level):
        self.network = network
        # Signal output gene index
        self.P_index = network.index(P_gene)
        self.level = level
        self.dim = self.network.size()
        self.parameter_graph = DSGRN.ParameterGraph(network)

    def classify_FP(self, morse_node, graded_complex):
        """Classify FP as quiescent or proliferative"""
        # Get cells in the Morse set. The Morse graph vertex (morse_node) is the cell grading
        morse_set = [c for c in graded_complex.complex()(self.dim) if graded_complex.value(c) == morse_node]
        if len(morse_set) != 1:
            return 'O'
        # Get cell coordinates in the cubical blowup complex
        cell_coords = graded_complex.complex().coordinates(morse_set[0])
        # Cells with 0 coordinate are boundary cells
        # Quiescent if lowest non-boundary coordinate
        if cell_coords[self.P_index] == 1:
            return 'q' # Quiescent
        # Proliferative if not lowest non-boundary coordinate
        if cell_coords[self.P_index] > 1:
            return 'p' # Proliferative
        return 'O'

    # def quiescent_FP(self, morse_node, graded_complex):
    #     # Get cells in the Morse set. The Morse graph vertex (morse_node) is the cell grading
    #     morse_set = [c for c in graded_complex.complex()(self.dim) if graded_complex.value(c) == morse_node]
    #     if len(morse_set) != 1:
    #         return False
    #     # Get cell coordinates in the cubical blowup complex
    #     cell_coords = graded_complex.complex().coordinates(morse_set[0])
    #     # Cells with 0 coordinate are boundary cells
    #     # Check if lowest non-boundary coordinate
    #     if cell_coords[self.P_index] == 1:
    #         return True
    #     return False

    # def proliferative_FP(self, morse_node, graded_complex):
    #     # Get cells in the Morse set. The Morse graph vertex (morse_node) is the cell grading
    #     morse_set = [c for c in graded_complex.complex()(self.dim) if graded_complex.value(c) == morse_node]
    #     if len(morse_set) != 1:
    #         return False
    #     # Get cell coordinates in the cubical blowup complex
    #     cell_coords = graded_complex.complex().coordinates(morse_set[0])
    #     # Cells with 0 coordinate are boundary cells
    #     # Check if not lowest non-boundary coordinate
    #     if cell_coords[self.P_index] > 1:
    #         return True
    #     return False

    def ClassifyParameter(self, par_index):
        parameter = self.parameter_graph.parameter(par_index)
        morse_graph, stg, graded_complex = DSGRN_utils.ConleyMorseGraph(parameter, level=self.level)
        conley_index = lambda node: morse_graph.vertex_label(node).split(':')[1].strip()
        # conley_index = lambda node: eval(morse_graph.vertex_label(node).split(':')[1].strip())
        FP_conley_index = lambda node: conley_index(node).count('1') == 1
        stable_nodes = [node for node in morse_graph.vertices() if len(morse_graph.adjacencies(node)) == 0]
        stable_FPs = [node for node in stable_nodes if FP_conley_index(node)]
        stable_FP_types = [self.classify_FP(node, graded_complex) for node in stable_FPs]
        quiescent = any(type == 'q' for type in stable_FP_types)
        proliferative = any(type == 'p' for type in stable_FP_types)
        if len(stable_FPs) == 1 and quiescent:
            # Quiescent and monostable
            return 'Q'
        if len(stable_FPs) == 1 and proliferative:
            # Proliferative and monostable
            return 'P'
        if quiescent and proliferative:
            # Quiescent and proliferative
            return 'B'
        return 'O'

    def is_hysteretic_path(self, path, strict=True, hyst_type='ascending'):
        # The labels are:
        # 'Q' : Quiescent and monostable
        # 'B' : Quiescent and proliferative (at least bi-stable)
        # 'P' : Proliferative and monostable
        #
        # A path exhibits (ascending) hysteresis if:
        # 1) The path starts with one or more 'Q'
        # 2) Followed by one or more 'B'
        # 3) The path ends with one or more 'P'
        #
        # A more strict notion of (ascending) hysteresis is:
        # 1) The path starts with a single 'Q'
        # 2) Followed by one or more 'B'
        # 3) The path ends with a single 'P'
        #
        # Get classification labels for each parameter in path
        labels = [self.ClassifyParameter(par_index) for par_index in path]
        # If strict make sure a single Q and a single P
        if strict and labels.count('Q') != 1:
            return False
        if strict and labels.count('P') != 1:
            return False
        # Remove repeated consecutive labels and make a string
        path_str = ''.join(key for key, group in itertools.groupby(labels))
        if hyst_type == 'ascending':
            # Path is hysteretic iff path_str == 'QBP'
            return path_str == 'QBP'
        # Else assume type == 'descending'
        # Path is hysteretic iff path_str == 'PBQ'
        return path_str == 'PBQ'

def reverse_topological_sort(graph):
    """Return list of vertices in reverse topologically sorted order."""
    result = []
    explored = set()
    dfs_stack = [(v,0) for v in graph.vertices]
    while dfs_stack:
        (v,i) = dfs_stack.pop()
        if (v,i) in explored: continue
        explored.add((v,i))
        if i == 0: # preordering visit
            dfs_stack.extend([(v,1)] + [(u,0) for u in graph.adjacencies(v)])
        elif i == 1: # postordering visit
            result.append(v)
    return result

def get_all_paths(graph, source=None, target=None, allowed=None):
    """Returns a list of all allowed paths with at least two edges
    from a source vertex u to a target vertex v. An allowed path
    is a path where allowed returns true for every vertex.
    """
    if source == None: source = lambda v : True
    if target == None: target = lambda v : True
    if allowed == None: allowed = lambda x : True
    top_sorted_verts = reverse_topological_sort(graph)
    paths = {} # Allowed paths from a vertex v to a taget vertex
    result = []
    num_verts = len(top_sorted_verts)
    for k, v in enumerate(top_sorted_verts):
        if not allowed(v):
            continue
        paths[v] = [[v] + path for u in graph.adjacencies(v) for path in paths[u] if allowed(u)]
        if target(v):
            paths[v].append([v])
        if source(v):
            result.extend([path for path in paths[v] if len(path) > 2])
        sys.stdout.write(f'\rComputing list of factor graph paths: {(k + 1) / num_verts:0.1%} complete.')
    sys.stdout.write(f'\rComputing list of factor graph paths: {num_verts / num_verts:0.1%} complete.\n')
    return result

def compute_hysteresis_paths(network, S_gene, P_gene, rpi, strict=True, hyst_type='ascending', level=3):
    """Given a network, a input gene name (S_gene), a output
    gene name (P_gene) and a relative parameter index (rpi)
    compute all hysteretic paths.
    """
    parameter_graph = DSGRN.ParameterGraph(network)
    # Get number of network nodes
    D = parameter_graph.dimension()
    # DSGRN uses an indexing scheme to refer to parameters. It is based on
    # a mixed-radix number scheme where the place value of each digit varies
    # according to the number of logic parameters for each node and the number
    # of order parameter for each node. Specifically, the ordering of the digits
    # is (from least significant) the sizes of each factor graph, followed by the
    # number of permutations of the out-edges for each node. We call these "bases"
    # (as in number base) and we compute the place value for each digit.
    indexing_place_bases = [parameter_graph.logicsize(d) for d in range(D)] + [parameter_graph.ordersize(d) for d in range(D)]
    indexing_place_values = functools.reduce(lambda x, y : x + [x[-1]*y], indexing_place_bases[:-1], [1])
    # Get index of input gene
    S_index = network.index(S_gene)
    # Get size of the factor graph associated with S_gene
    num_gene_param = indexing_place_bases[S_index]
    # Get the product of the sizes of all remaining factor graphs,
    # and reorderings of all genes (including S_gene)
    num_reduced_param = parameter_graph.size() / num_gene_param
    # Create factor graph
    hexcodes = parameter_graph.factorgraph(S_index)
    vertices = set(range(len(hexcodes)))
    edges = [(u, v) for u in vertices for v in vertices if DSGRN.isAdjacentHexcode(hexcodes[u], hexcodes[v])]
    factor_graph = DSGRN.Graph(vertices, edges)
    # Giving a reduced parameter index (rpi) and a gene parameter index (gpi),
    # corresponding to gene S_index, get the corresponding full parameter index
    full_parameter_index = lambda rpi, gpi, S_index : rpi % indexing_place_values[S_index] + \
                                                      gpi * indexing_place_values[S_index] + \
                                                      (rpi // indexing_place_values[S_index]) * \
                                                      indexing_place_values[S_index + 1]
    # Give a full parameter index (pi) get the reduced parameter index
    reduced_parameter_index = lambda pi, S_index : (pi % indexing_place_values[S_index] + \
                                                   (pi // indexing_place_values[S_index + 1]) * \
                                                   indexing_place_values[S_index], \
                                                   (pi // indexing_place_values[S_index]) % \
                                                   indexing_place_bases[S_index])
    # Set network analyzer
    network_analyzer = NetworkAnalyzer(network, P_gene, level)
    # Get all paths in S_gene factor graph
    factor_paths = get_all_paths(factor_graph)
    # Get hysteresis paths
    hysteresis_paths = []
    num_factor_paths = len(factor_paths)
    for k, factor_path in enumerate(factor_paths):
        # Get path with indices in full indices
        path = [full_parameter_index(rpi, gpi, S_index) for gpi in factor_path]
        if network_analyzer.is_hysteretic_path(path, strict=strict, hyst_type=hyst_type):
            hysteresis_paths.append(path)
        sys.stdout.write(f'\rAnalyzing factor graph paths for hysteresis: {(k + 1) / num_factor_paths:0.1%} complete.')
    sys.stdout.write(f'\rAnalyzing factor graph paths for hysteresis: {num_factor_paths / num_factor_paths:0.1%} complete.\n')
    # Get the total number of paths
    num_paths = len(factor_paths)
    return num_paths, hysteresis_paths
