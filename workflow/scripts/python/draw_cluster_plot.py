from itertools import combinations
#   import matplotlib.pyplot as plt
import math
import pydot
import pandas as pd
import pydot
from PIL import Image

def set_node_label(row):
    return f'{row.SampleAlias}'

class FullGraph:
    # initialize object with manifest data and cdata
    def __init__(self, manifest_df, cdata_df):
        self.df = manifest_df
        self.dfc = cdata_df
        self.node_mappings = { 
            "ObservedGender": { "F": "circle", "M": "square", "U": "star", "nan": "star" },
            "SampleLabel":    { "Number": "square", "NWD": "diamond", "TOE": "triangle" },
            "SampleAlias":    { "Number": "square", "NWD": "diamond", "TOE": "triangle" },
            "SampleTypeCode": { "d": "square", "m": "circle", "r": "diamond" },
            "dataset_id":     { "tmp0": "square", "tmp1": "diamond", "tmp2": "triangle", "tmp3": "circle" },
            "S_SAMPLEID":    { "tmp0": "square", "tmp1": "diamond", "tmp2": "triangle", "tmp3": "circle" }
        }
        self.mapping_funcs = {
            "ObservedGender": lambda x: x,
            "SampleLabel": None,
            "SampleAlias": None,
            "SampleTypeCode": lambda x: x,
            "dataset_id": None,
            "S_SAMPLEID": lambda x: x[0:3] + x[4]
        }
    
    # method to call to generate graph
    def graph(self, shape_by):
        # filtering by group_id
        df_filter_groupid = self.df
        df_filter_groupid['NodeLabel'] = df_filter_groupid.apply(set_node_label, axis=1)
        
        # get nodes
        nodes = df_filter_groupid["SampleLabel"].values
        
        # group nodes by shape mappings
        #shape_nodes = self.map_nodes(df_filter_groupid, shape_by)
        shape_nodes = []

        # getting all combinations of samplelabels
        #edges, labels = self.generate_edges(nodes)
        edges = []
        labels = []
        group_id = None

        return self.generate_graph(shape_nodes, edges, labels, group_id, shape_by, df_filter_groupid)
        
    def generate_graph(self, shape_nodes, edges, labels, groupid, shape_by, df_filter_groupid):
        
        ### Get the list of subject ids and set colors for subjects
        subjects = list(set(df_filter_groupid["S_SUBJECTID"].values))
        
        # create graph
        G = pydot.Dot('tmp/graph_path.png', graph_type="graph", compound='true', fontname='Verdana')
        
        # create legend
        # print(f'{subjects}')
        ## node = pydot.Node(shape_nodes[node_shape]["attribute"], shape=node_shape, style="filled", fillcolor="lightblue")
                
        for (s, subject) in enumerate(subjects):
            cluster = pydot.Cluster(graph_name=subject, label=f"{subject}", 
                                   fontsize="10", style="filled", fillcolor="white")

            # add nodes
            samples = df_filter_groupid[df_filter_groupid.S_SUBJECTID == subject]
            fill = s + 1
            for (s_idx, sample) in samples.iterrows():
                node_label = f'{sample.NodeLabel}'
                sample_label = sample.SampleLabel
                shape = self.node_mappings['SampleTypeCode'][sample.SampleTypeCode]
                node = pydot.Node(sample_label, label=node_label, shape=shape, style="filled", colorscheme="accent8", color=fill)
                cluster.add_node(node)

            G.add_subgraph(cluster)
        
        # add edges
        for edge in edges:
            edge1 = edge[0]
            edge2 = edge[1]
            # if edge1[0:3] == "TOE":
            #     edge1 = edge1[0:9]
            # if edge2[0:3] == "TOE":
            #     edge2 = edge2[0:9]
            
            new_edge = pydot.Edge(edge1, edge2, label=labels[edge])
            G.add_edge(new_edge)
        
        #G.create_png(graph_path)
        G.write_png('tmp/huge_graph.png')
        return G.create_png()
        
    # returns a list of edges and a list of edge labels
    # given a list of nodes
    def generate_edges(self, n):
        all_edges = list(combinations(n, 2))
        filtered_edges = []
        labels = {}
        
        for edge in all_edges:
            try:
                query_string = f'ID1 == "{edge[0]}" and ID2 == "{edge[1]}"'
                concord = self.dfc.query(query_string)["Concord"].values[0]

                # only create edge for pairs that have concord values
                if not math.isnan(concord):
                    labels[edge] = concord
                    filtered_edges.append(edge)
            except:
                pass
                
        return filtered_edges, labels
    
    # returns a dictionary with nodes grouped by shape such as { "shape": { "nodes": [nodes], "attribute": "attr" }... }
    # given a list of nodes
    def map_nodes(self, n, by):
        shape_nodes = {}
        
        try:
            curr_map = self.node_mappings[by]
        except Exception:
            raise Exception(f"shape_by is not a valid option! ({self.node_mappings.keys()})")
        
        for (node_index, node) in n.iterrows():
            # get attribute for node
            node_label = node['SampleLabel']
            query_string = f'SampleLabel == "{node}"'
            
            # cast to String for nan case
            #attr = str(self.df.query(query_string)[by].values[0])
            attr = node_label

            # get final attribute String
            # will not throw exception because of try except above
            # final_attr = self.mapping_funcs[by](attr)
            final_attr = 'SampleLabel'
            
            # get mapped shape
            # shape = curr_map[final_attr]
            shape = 'circle'

            shape_nodes.setdefault(shape, {"nodes": [], "attribute": final_attr})
            shape_nodes[shape]["nodes"].append(node)
                
        return shape_nodes
    

def main():
    multiomics_manifest_extracted = pd.read_csv(snakemake.input.manifest)
    c_data = pd.read_csv(snakemake.input.c_data)
    full_graph_generator = FullGraph(multiomics_manifest_extracted, c_data)
    pdot = full_graph_generator.graph(shape_by='SampleTypeCode')
    # plt = Image(pdot)
    display(plt)

if __name__=="__main__":
    main()