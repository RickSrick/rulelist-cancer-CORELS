import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import NullFormatter, FixedLocator
import csv

colors = ['xkcd:coral', 'xkcd:beige', 'xkcd:teal', 'xkcd:goldenrod', 'xkcd:lavender']

matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})

def prune_effect():
    data = "prune.csv"

    table = pd.read_csv(data)
    table = table[["testname", "percentage", "rule_score"]]
    table_1 = table[table["testname"].isin(["OV o COADREAD"])]
    table_2 = table[table["testname"].isin(["LUAD o LUSC"])]

    table_1 = table_1[["percentage", "rule_score"]].groupby(["percentage"],as_index=False).mean()
    table_2 = table_2[["percentage", "rule_score"]].groupby(["percentage"],as_index=False).mean()

    var_1 = table_1["rule_score"]

    print(table_1)
    print(table_2)

def plot_node_precision(str, locp):
    data = "../data1.csv"

    table = pd.read_csv(data)
    table = table[table["testname"].isin([str])]
    table = table.groupby(["testname","percentage","nodes","policy","min_support"], as_index=False).mean()
    table = table[["nodes", "min_support", "rule_score"]].groupby(["nodes", "min_support"], as_index=False).mean()
    
    table001 = table[table["min_support"].isin([0.01])]
    table01 = table[table["min_support"].isin([0.1])]
    table03 = table[table["min_support"].isin([0.3])]
    table04 = table[table["min_support"].isin([0.4])]

    #table001,table01,table03,table04
    p = [table001,table01,table03,table04]
    #"0.01", "0.1", "0.3","0.4"
    x = ["0.01", "0.1", "0.3","0.4"]
    k = 0
    for v in p:
        plt.plot(v["nodes"], v["rule_score"], "-o", color=colors[k], label=x[k])
        k+= 1
    
    plt.ylabel('Score del modello')
    plt.xlabel('nodi massimi')
    plt.xscale('log')
    plt.title('andamento precisione con incremento dei nodi')
    plt.legend(loc=locp, title="min support")
    plt.savefig(str+'_nodes.pgf')
    plt.close()
    return 0

def mean_score(str):
    data = "../data1.csv"

    table = pd.read_csv(data)
    table = table[table["testname"].isin([str])]
    v = table["rule_score"].mean()
    print(v)

def max_score(str):
    data = "../data1.csv"

    table = pd.read_csv(data)
    table = table[table["testname"].isin([str])]
    table = table[table["rule_score"].isin([table["rule_score"].max()])]
    print(table)  

def plot_node_precision_algorithm(str, num):
    data = "../data1.csv"

    table = pd.read_csv(data)
    table = table[table["testname"].isin([str])]
    table = table.groupby(["testname","percentage","nodes","policy","min_support"], as_index=False).mean()
    table = table[["nodes", "policy", "min_support", "rule_score"]].groupby(["nodes", "policy", "min_support"], as_index=False).mean()
    
    table001 = table[table["min_support"].isin([num]) & table["policy"].isin(["objective"])]
    table01  = table[table["min_support"].isin([num]) & table["policy"].isin(["lower_bound"])]
    table02  = table[table["min_support"].isin([num]) & table["policy"].isin(["curious"])]
    table03  = table[table["min_support"].isin([num]) & table["policy"].isin(["bfs"])]

    p = [table001,table01,table02,table03]
    x = ["objective", "lower_bound", "curious", "bfs"]
    k = 0
    for v in p:
        plt.plot(v["nodes"], v["rule_score"], "-o", color=colors[k], label=x[k])
        k+= 1
    
    plt.ylabel('Score del modello')
    plt.xlabel('nodi massimi')
    plt.title('andamento precisione con incremento dei nodi')
    plt.legend(loc='upper right')
    plt.savefig(str+"_"+num+'_nodes.pgf')
    plt.show()
    return 0

def draw_plot_performance():

    data = "../data1.csv"
    table = pd.read_csv(data)
    table["p1"] = 1

    t1 = table.groupby(["min_support"], as_index=False).sum()
    print(t1["p1"])
    rule_list = table[table["winner"] == 1].groupby(["min_support"], as_index=False).sum()
    equal = table[table["winner"] == 2].groupby(["min_support"], as_index=False).sum()
    tree_list = table[table["winner"] == 0].groupby(["min_support"], as_index=False).sum()

    rule_list = rule_list["p1"]
    equal = equal["p1"]
    tree_list = tree_list["p1"]

    print(rule_list)
    print(equal)
    print(tree_list)

    penguin_means = {
        'rule list': [rule_list, colors[0]],
        'parità': [equal, colors[2]],
        'decision tree': [tree_list, colors[1]]
    }

    names = tuple(table["min_support"].astype(str).unique())

    print(names)

    x = np.arange(len(names))  # the label locations
    width = 0.25  # the width of the bars
    multiplier = 0.4

    fig, ax = plt.subplots(layout='constrained')

    for attribute, measurement in penguin_means.items():
        offset = width * multiplier
        rects = ax.bar(x + offset, measurement[0], width, label=attribute, color=measurement[1])
        ax.bar_label(rects, padding=3)
        multiplier += 1

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('numero di test')
    ax.set_xlabel('min support')
    ax.set_title('comparazione tra decision e tree rulelist in base al min support')
    ax.set_xticks(x + width, names)
    ax.legend(loc='upper left', ncols=4)
    ax.set_ylim(0, 140)
    plt.savefig('tree_rule.pgf')
