import numpy as np
import networkx as nx
from progressbar import ProgressBar, Percentage, Bar, RotatingMarker, ETA, FileTransferSpeed
import itertools

def co_express_network(g):
    h = nx.Graph()
    edges = dict()
    for fin in g.ad_heat.baskets():
        h.add_node(fin, weight=fin.activity * fin.cluster_score, mass=fin.m, rt=fin.rt, bid=fin.name, num=len(fin.keys()))
    widgets = ['VectorMove: ', Percentage(), ' ', Bar(marker=RotatingMarker()), ' ',\
	ETA(), ' ', FileTransferSpeed()]
	
    pbar = ProgressBar(widgets=widgets, maxval=len(g.ad_heat.idruns())).start()
    
    extracts = set([idrun.name[:-2] for idrun in g.ad_heat.idruns()])
	
    for num, extract in enumerate(extracts):
        pbar.update(num)
        basks = []
        for idrun in [extract + let + "_" for let in 'ABCDEF']:
            try:
                basks.extend(g.ad_heat[idrun].keys())
            except KeyError:
                continue
        basks = list(set(basks))
        for pair in itertools.combinations(basks, 2):
            edges.setdefault((g.ad_heat[pair[0]], g.ad_heat[pair[1]]), [])
#            edges[(g.ad_heat[pair[0]], g.ad_heat[pair[1]])].append(np.log(g.bask_prob(g.ad_heat[pair[0]], g.cp[idrun.name[:-1]]) * g.bask_prob(g.ad_heat[pair[1]], g.cp[idrun.name[:-1]])))
            edges[(g.ad_heat[pair[0]], g.ad_heat[pair[1]])].append(1)
    pbar.finish()
    for edge in edges:
        h.add_edge(*edge, weight=len(edges[edge]))
    return h