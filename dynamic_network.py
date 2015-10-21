import networkx as nx
import numpy as np
nodes = dict()
edges = dict()

for iteration in np.arange(0, 300):
    f = nx.read_dot("cp_{}_cp_network.dot".format(iteration))
    for node in f.nodes(data=True):
        if node[0] in nodes:
            for at in node[1]:
                if nodes[node[0]][at][-1][2] == node[1][at]:
                    continue
                else:
                    nodes[node[0]][at][-1][1] = iteration
                    nodes[node[0]][at].append([iteration, 'Infinity', node[1][at]])
            if not nodes[node[0]]['Time Interval'][-1][1] == 'Infinity':
                nodes[node[0]]['Time Interval'].append([iteration, 'Infinity'])
            nodes[node[0]]['present'] = 1    
        else:
            nodes[node[0]] = dict()
            for at in node[1]:
                nodes[node[0]][at] = [[iteration, 'Infinity', node[1][at]]]
            nodes[node[0]]['Time Interval'] = [[iteration, 'Infinity']]
            nodes[node[0]]['present'] = 1
    for node in nodes:
        if nodes[node]['present'] == 1:
            nodes[node]['present'] = 0
            continue
        if nodes[node]['Time Interval'][-1][1] == 'Infinity':
            nodes[node]['Time Interval'][-1][1] == iteration
        else:
            continue
            
            
    for edge in f.edges(data=True):
        if (edge[0], edge[1]) in edges:
            for at in edge[2]:
                if edges[(edge[0], edge[1])][at][-1][2] == edge[2][at]:
                    continue
                else:
                    edges[(edge[0], edge[1])][at][-1][1] = iteration
                    edges[(edge[0], edge[1])][at].append([iteration, 'Infinity', edge[2][at]])
            if not edges[(edge[0], edge[1])]['Time Interval'][-1][1] == 'Infinity':
                edges[(edge[0], edge[1])]['Time Interval'].append([iteration, 'Infinity'])
            edges[(edge[0], edge[1])]['present'] = 1    
        else:
            edges[(edge[0], edge[1])] = dict()
            for at in edge[2]:
                edges[(edge[0], edge[1])][at] = [[iteration, 'Infinity', edge[2][at]]]
            edges[(edge[0], edge[1])]['Time Interval'] = [[iteration, 'Infinity']]
            edges[(edge[0], edge[1])]['present'] = 1
    for edge in edges:
        if edges[edge]['present'] == 1:
            edges[edge]['present'] = 0
            continue
        if edges[edge]['Time Interval'][-1][1] == 'Infinity':
            edges[edge]['Time Interval'][-1][1] == iteration
        else:
            continue

for node in nodes:
    del nodes[node]['present']
for edge in edges:
    del edges[edge]['present']
    
g = nx.Graph()
for node in nodes:
    g.add_node(
    
#g = open('dyn_net.gexf', 'w')
#
#g.write('<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n')
#g.write('<gexf>\n')
#g.write('\t<graph mode="dynamic" defaultedgetype="undirected" timeformat=\"time\">')
#idnum = 0
#g.write('\t\t<attributes class=\"node\" mode=\"dynamic\">\n')
#for at in nodes.values()[0]:
#    g.write('\t\t\t<attribute id="{}" title="{}" type="{}"/>'.format(idnum, at, type(nodes.values()[0][at][0][-1])))
#g.write('\t\t</attributes>')
#g.write('\t\t<attributes class=\"edge\" mode=\"dynamic\">\n')
#for at in edges.values()[0]:
#    g.write('\t\t\t<attribute id="{}" title="{}" type="{}"/>'.format(idnum, at, type(edges.values()[0][at][0][-1])))
#g.write('\t\t</attributes>')
#g.write('\t\t<nodes>')
#for node in nodes:
#    g.write('\t\t\t<node id="{}" label="{}" start="{}"


    
            
#n = open('nodes.csv', 'w')
#e = open('edges.csv', 'w')
#
#n.write("Nodes\tId\tLabel\t{}\n".format('\t'.join(nodes.values()[0].keys())))
#for node in nodes:
#    ats = []
#    for at in nodes[node]:
#        gs = []
#        for group in nodes[node][at]:
#            gs.append("[{})".format(", ".join(map(str, group))))
#        ats.append("\"<{}>\"".format("; ".join(gs)[:-1] + "]"))
#    n.write("{}\t{}\t{}\t{}\n".format(node, node, node, "\t".join(ats)))
#n.close()
#    
#e.write("Source\tTarget\tType\tId\tLabel\t{}\n".format('\t'.join(edges.values()[0].keys())))
#for edge in edges:
#    num = 1
#    ats = []
#    for at in edges[edge]:
#        gs = []
#        for group in edges[edge][at]:
#            gs.append("[{})".format(", ".join(map(str, group))))
#        ats.append("\"<{}>\"".format("; ".join(gs)[:-1] + "]"))
#    e.write("{}\t{}\tUndirected\t{}\t\t{}\n".format(edge[0], edge[1], num, "\t".join(ats)))
#e.close()