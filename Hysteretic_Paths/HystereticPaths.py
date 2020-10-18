# HystereticPaths.py
# Marcio Gameiro
# 2020-10-18
# MIT LICENSE

import DSGRN
import itertools
import functools

class NetworkAnalyzer:
    def __init__(self, network, P_gene):
        self.network = network
        # Signal output gene index
        self.P_index = network.index(P_gene)
        self.parameter_graph = DSGRN.ParameterGraph(network)

    def is_FP(self, annotation):
        return annotation.startswith("FP")

    def is_quiescent_FP(self, annotation):
        if self.is_FP(annotation):
            digits = [int(s) for s in annotation.replace(",", "").split() if s.isdigit()]
            if digits[self.P_index] == 0:
                return True
        return False

    def is_proliferative_FP(self, annotation):
        if self.is_FP(annotation):
            digits = [int(s) for s in annotation.replace(",", "").split() if s.isdigit()]
            if digits[self.P_index] > 0:
                return True
        return False

    def ClassifyParameter(self, par_index):
        parameter = self.parameter_graph.parameter(par_index)
        domain_graph = DSGRN.DomainGraph(parameter)
        morse_decomp = DSGRN.MorseDecomposition(domain_graph.digraph())
        morse_graph = DSGRN.MorseGraph(domain_graph, morse_decomp)
        morse_graph_poset = morse_graph.poset()
        stable_annotations = [morse_graph.annotation(n)[0] for n in range(morse_graph.poset().size()) if len(morse_graph.poset().children(n)) == 0]
        monostable = len(stable_annotations) == 1
        quiescent = any(self.is_quiescent_FP(annotation) for annotation in stable_annotations)
        proliferative = any(self.is_proliferative_FP(annotation) for annotation in stable_annotations)
        if monostable and quiescent:
            return 'Q'
        if monostable and proliferative:
            return 'P'
        if quiescent and proliferative:
            return 'B'
        if quiescent:
            return 'q'
        if proliferative:
            return 'p'
        return 'O'

    def is_hysteretic_path(self, path, hyst_type='ascending'):
        # Recalling that the labels are:
        # 'Q' : Quiescent and monostable
        # 'q' : Quiescent but not monostable
        # 'B' : Quiescent and proliferative (at least bi-stable)
        # 'P' : Proliferative and monostable
        # 'p' : Proliferative but not monostable
        # A path is defined as hysteretic if the following are true:
        # 1) The path starts with a 'Q'
        # 2) Followed by one or more 'Q' or 'q'
        # 3) Followed by one or more 'B'
        # 4) Followed by one or more 'P' or 'p'
        # 5) The path ends with a 'P'

        # Get classification label for each parameter in path
        labels = [self.ClassifyParameter(par_index) for par_index in path]
        # Transform labels into string
        path_str = ''.join(labels)
        if hyst_type == 'ascending':
            # Return false if does not start with 'Q'
            if not path_str.startswith('Q'):
                return False
            # Return false if does not end with 'P'
            if not path_str.endswith('P'):
                return False
        else: # Assume type == 'descending'
            # Return false if does not start with 'P'
            if not path_str.startswith('P'):
                return False
            # Return false if does not end with 'Q'
            if not path_str.endswith('Q'):
                return False
        # Replace 'q' by 'Q' and 'p' by 'P'
        path_str = path_str.replace('q', 'Q')
        path_str = path_str.replace('p', 'P')
        # Remove repeated consecutive labels from string
        path_str = ''.join(key for key, group in itertools.groupby(list(path_str)))
        if hyst_type == 'ascending':
            # Path is hysteretic iff path_str == 'QBP'
            return path_str == 'QBP'
        else: # Assume type == 'descending'
            # Path is hysteretic iff path_str == 'PBQ'
            return path_str == 'PBQ'

    def is_hysteretic_path_strict(self, path, hyst_type='ascending'):
        # A more strict notion of hysteresis is to require:
        # 1) The path starts with one or more 'Q'
        # 2) Followed by one or more 'B'
        # 3) The path ends with one or more 'P'

        # Get classification label for each parameter in path
        labels = [self.ClassifyParameter(par_index) for par_index in path]
        # Remove repeated consecutive labels and make a string
        path_str = ''.join(key for key, group in itertools.groupby(labels))
        if hyst_type == 'ascending':
            # Path is hysteretic iff path_str == 'QBP'
            return path_str == 'QBP'
        else: # Assume type == 'descending'
            # Path is hysteretic iff path_str == 'PBQ'
            return path_str == 'PBQ'

def reverse_topological_sort(graph):
    """Return list of vertices in reverse topologically sorted order.
    """
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
    for v in top_sorted_verts:
        if not allowed(v):
            continue
        paths[v] = [[v] + path for u in graph.adjacencies(v) for path in paths[u] if allowed(u)]
        if target(v):
            paths[v].append([v])
        if source(v):
            result.extend([path for path in paths[v] if len(path) > 2])
    return result

def compute_hysteresis_paths(network, S_gene, P_gene, rpi, strict=True, hyst_type='ascending'):
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
    network_analyzer = NetworkAnalyzer(network, P_gene) 
    # Get all paths in S_gene factor graph
    factor_paths = get_all_paths(factor_graph)
    # Get hysteresis paths
    hysteretic_paths = []
    for factor_path in factor_paths:
        # Get path with indices in full indices
        path = [full_parameter_index(rpi, gpi, S_index) for gpi in factor_path]
        if strict:
            if network_analyzer.is_hysteretic_path_strict(path, hyst_type):
                hysteretic_paths.append(path)
        else:
            if network_analyzer.is_hysteretic_path(path, hyst_type):
                hysteretic_paths.append(path)
    # Get the total number of paths
    num_paths = len(factor_paths)
    return num_paths, hysteretic_paths
