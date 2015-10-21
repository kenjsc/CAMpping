import networkx as nx
import numpy as np

def parse_timeline(num=300):
    nodes = dict()
    edges = dict()

    for iteration in np.arange(0, num):
        f = nx.read_dot("cp_{}_cp_network.dot".format(iteration))
        for node in f.nodes(data=True):
            if node[0] in nodes:
                for at in node[1]:
                    if nodes[node[0]][at][-1][0] == node[1][at]:
                        continue
                    else:
                        nodes[node[0]][at][-1][2] = iteration
                        nodes[node[0]][at].append([node[1][at], iteration, None])
                if not nodes[node[0]]['spells'][-1][1] == None:
                    nodes[node[0]]['spells'].append([iteration, None])
                nodes[node[0]]['present'] = 1    
            else:
                nodes[node[0]] = dict()
                for at in node[1]:
	                nodes[node[0]][at] = [[node[1][at], iteration, None]]
                nodes[node[0]]['spells'] = [[iteration, None]]
                nodes[node[0]]['present'] = 1
        for node in nodes:
            if nodes[node]['present'] == 1:
                nodes[node]['present'] = 0
                continue
            if nodes[node]['spells'][-1][1] == None:
                nodes[node]['spells'][-1][1] == iteration
            else:
                continue
            
            
        for edge in f.edges(data=True):
            if (edge[0], edge[1]) in edges:
                for at in edge[2]:
                    if edges[(edge[0], edge[1])][at][-1][0] == edge[2][at]:
 	                   continue
                    else:
                        edges[(edge[0], edge[1])][at][-1][2] = iteration
                        edges[(edge[0], edge[1])][at].append([edge[2][at], iteration, None])
                if not edges[(edge[0], edge[1])]['spells'][-1][1] == None:
                    edges[(edge[0], edge[1])]['spells'].append([iteration, None])
                edges[(edge[0], edge[1])]['present'] = 1    
            else:
                edges[(edge[0], edge[1])] = dict()
                for at in edge[2]:
                    edges[(edge[0], edge[1])][at] = [[edge[2][at], iteration, None]]
                edges[(edge[0], edge[1])]['spells'] = [[iteration, None]]
                edges[(edge[0], edge[1])]['present'] = 1
        for edge in edges:
            if edges[edge]['present'] == 1:
                edges[edge]['present'] = 0
                continue
            if edges[edge]['spells'][-1][1] == None:
                edges[edge]['spells'][-1][1] = iteration
            else:
                continue
    
    for node in nodes:
        del nodes[node]['present']
    for edge in edges:
        del edges[edge]['present']
    return nodes, edges
    
def gen_graph(nodes, edges):
    
    g = nx.Graph()
    g.add_nodes_from([(node, ats) for node, ats in nodes.items()])
    g.add_edges_from([(edge[0], edge[1], ats) for edge, ats in edges.items()])
    return g

def write_graph(g):
    nx.write_gexf(g, 'dynamic_gexf.gexf')

def write_direct(nodes, edges):    
    g = open('dyn_net.gexf', 'w')
    
    g.write('<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n')
    g.write('<gexf>\n')
    g.write('  <graph mode="dynamic" defaultedgetype="undirected" timeformat=\"Integer\">\n')
    g.write('    <attributes class=\"node\" mode=\"dynamic\">\n')
    g.write('      <attribute id="0" title="weight" type="float"/>\n')
    g.write('    </attributes>\n')
    g.write('    <attributes class=\"edge\" mode=\"dynamic\">\n')
    g.write('      <attribute id="1" title="weight" type="float"/>\n')
    g.write('    </attributes>\n')
    g.write('    <nodes>\n')
    for node in nodes:
        g.write('      <node id="{}" label="{}">\n'.format(node, node))
        g.write('        <spells>\n')
        for start, end in nodes[node]['spells']:
            if end == None:
                g.write('          <spell start="{}" />\n'.format(start))
            else:
                g.write('          <spell start="{}" end="{}" />\n'.format(start, end))
        g.write('        </spells>\n')
        g.write('        <attvalues>\n')
        for val, start, end in nodes[node]['weight']:
            if end == None:
                g.write('          <attvalue for="0" value="{}" start="{}" />\n'.format(float(val), start))
            else:
                g.write('          <attvalue for="0" value="{}" start="{}" endopen="{}" />\n'.format(float(val), start, end))
        g.write('        </attvalues>\n')
        g.write('      </node>\n')
    g.write('    </nodes>\n')
    g.write('    <edges>\n')
    n = 0
    for edge in edges:
        g.write('      <edge id="{}" source="{}" target="{}">\n'.format(n, edge[0], edge[1]))
        g.write('        <spells>\n')
        for start, end in edges[edge]['spells']:
            if end == None:
                g.write('          <spell start="{}" />\n'.format(start))
            else:
                g.write('          <spell start="{}" end="{}" />\n'.format(start, end))
        g.write('        </spells>\n')
        g.write('        <attvalues>\n')
        for val, start, end in edges[edge]['weight']:
            if end == None:
                g.write('          <attvalue for="1" value="{}" start="{}" />\n'.format(float(val), start))
            else:
                g.write('          <attvalue for="1" value="{}" start="{}" endopen="{}" />\n'.format(float(val), start, end))
        g.write('        </attvalues>\n')
        g.write('      </edge>\n')
        n += 1
    g.write('    </edges>\n')
    g.write('  </graph>\n')
    g.write('</gexf>\n')