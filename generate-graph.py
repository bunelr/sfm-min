#!/usr/bin/env python
import argparse
import random

def write_graph(out_file, graph):
    nb_nodes = len(graph)
    # Write all the middle nodes, then the source node, then the sink
    for from_node_idx in range(1, nb_nodes-1) + [0, nb_nodes-1]:
        for to_node_idx in range(1, nb_nodes-1):
            out_file.write(str(graph[from_node_idx][to_node_idx]) + ' ')
        out_file.write(str(graph[from_node_idx][0]) + ' ')
        out_file.write(str(graph[from_node_idx][nb_nodes-1]) + '\n')


def generate_graph(number_nodes, edge_proba):
    # The graph is going to be represented as a list of arrays
    # Each array representing the outgoing arcs from a node
    graph = []

    for node_idx in range(number_nodes):
        outgoing_arcs = [0 for _ in range(number_nodes)]
        graph.append(outgoing_arcs)

        if node_idx == 0:
            continue

        # Ensure that we get at least one edge going to this
        reachable = False
        while not reachable:
            for prev_node_idx in range(node_idx):
                # Only add nodes from the previous edges
                if random.random() < edge_proba:
                    # Need to add this edge.
                    # The value is
                    edge_val = random.randint(1, 100)
                    graph[prev_node_idx][node_idx] = edge_val
                    reachable = True

    # Find all the nodes that can lead to the sink
    can_reach_sink = set([number_nodes-1])
    to_process = [number_nodes-1]
    while to_process:
        processing = to_process.pop()
        for node_idx in range(number_nodes):
            if graph[node_idx][processing] > 0:
                if node_idx not in can_reach_sink:
                    to_process.append(node_idx)
                    can_reach_sink.add(node_idx)

    can_reach_sink.discard(0)  # don't want possible to lead to source
    for node_idx in range(number_nodes):
        if node_idx not in can_reach_sink:
            edge_val = random.randint(1, 100)
            target_node = random.choice(list(can_reach_sink))
            graph[node_idx][target_node] = edge_val

    return graph

def main():
    parser = argparse.ArgumentParser(description='Generate textfile representing graphs')
    parser.add_argument('-nb', default=20, type=int,
                        help='how many nodes')
    parser.add_argument('-out', default="graph.txt",
                        type=argparse.FileType('w'),
                        help="Which file to write the graph to")
    parser.add_argument('-edge_proba', default=0.2,
                        type=float,
                        help="Lower bound on the probability of creating an edge")
    args = parser.parse_args()
    graph = generate_graph(args.nb, args.edge_proba)
    write_graph(args.out, graph)

if __name__=='__main__':
    main()
