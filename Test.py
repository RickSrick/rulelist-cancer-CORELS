from sklearn.model_selection import train_test_split
import matplotlib
import matplotlib.pyplot as plt
import NotificationManager
from sklearn import tree
from corels import *
import pandas as pd
import numpy as np
import time
import csv
import os
import gc

matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})


MAX_CARD_SUPPORT = 2
TEST_SIZE = 0.2         #percentage of test size 80-20

class Test:

    def __init__(self, files:list, cancer_type_one:str, cancer_type_two:str)-> None:
        self.hotnet2, self.snvs_strictly_filtered = files
        self.cancer_type_one = cancer_type_one
        self.cancer_type_two = cancer_type_two
        self.hotnet2_table = self.__get_hotnet2_table()
        self.mutations, self.snvs_table = self.__get_snvs_strictly_table()
        NotificationManager.secret_set_env()

    def launch_test(self, percentage:float, number_node:int, policy:str, min_support:float)-> float: 
        
        print("test launched")
        start_time = time.perf_counter()
        score_corels, model_corels = self.__CORELS_body(percentage, number_node, policy, min_support)
        end_time = time.perf_counter()

        score_tree, model_tree = self.__tree_body(percentage)
        execution_time = end_time - start_time

        #START PRODUCING OUTPUT
        self.__output_test(percentage, number_node, policy, min_support, model_corels, model_tree, score_corels, score_tree, execution_time)

        del model_corels, model_tree
        gc.collect()

    def __output_test(self, percentage, number_node, policy, min_support, model_corels, model_tree, score_corels, score_tree, execution_time):
        filename = self.__set_enviroment(number_node, policy, min_support)
        rule_file = open(filename + ".txt", "a")
        
        if os.path.isfile("data1.csv") == False:
            data_file = open("data1.csv", "a")
            data_file.write("testname,percentage,nodes,policy,min_support,time,rule_score,tree_score,winner\n")
        else:
            data_file = open("data1.csv","a")

        #0 = tree WIN
        #1 = rule WIN
        #2 = equal WIN
        win = 0
        if score_tree < score_corels:
            win = 1
        elif score_corels == score_tree:
            win = 2
            

        data_file.write(f'''{self.cancer_type_one} o {self.cancer_type_two},{percentage},{number_node},{policy},{min_support},{round(execution_time, 3)},{round(score_corels, 3)},{round(score_tree,3)},{win}\n''')
        rule_file.write(str(model_corels.rl())+
            "\nSCORE:"+str(score_corels)+
            "\n"+("="*30)+"\n")
        
        plt.figure(figsize=(50,50))
        tree.plot_tree(model_tree, feature_names=self.__prune_snvs_table(percentage)[0])
        plt.savefig(filename + '.pgf')
        
        data_file.close()
        rule_file.close()
        plt.close()

    def __set_enviroment(self, number_node, policy, min_support)-> str:
        name = f"{self.cancer_type_one}-{self.cancer_type_two}_n:{str(number_node)}_p:{policy}_m:{str(min_support)}"
        path = r'tests/'+name 

        if not os.path.exists(path):
            os.makedirs(path)
        
        return path+"/"+name
    
    def __get_hotnet2_table(self)-> pd.DataFrame:
        hotnet2_table = pd.read_csv(self.hotnet2, sep=" ", names=["Patient", "CancerType"])
        hotnet2_table = hotnet2_table[hotnet2_table["CancerType"].isin([self.cancer_type_one, self.cancer_type_two])]
        hotnet2_table.CancerType[hotnet2_table.CancerType == self.cancer_type_one] = 1
        hotnet2_table.CancerType[hotnet2_table.CancerType == self.cancer_type_two] = 0
        
        return hotnet2_table

    def __get_snvs_strictly_table(self):
        mutations = {}
        patients = []
        line_num = 0

        """
        Reading TSV file, saving 
            - Key : mutation_name (str)
            - Value : list of patients with mutation (list(int))
        """
        with open(self.snvs_strictly_filtered) as file:
            tsv_file = csv.reader(file, delimiter="\t")

            for line in filter(lambda x: x[0] in self.hotnet2_table['Patient'].values,tsv_file):
                patients.append(line[0])
                for x in line[1:]:
                    if x in mutations:
                        mutations[x].append(line_num)
                    else:
                        mutations[x] = [line_num]
                line_num += 1

        quantity_mutation = dict(sorted(mutations.items(), key=lambda item: len(item[1]), reverse=True))
        binary_data = np.zeros(shape=(len(self.hotnet2_table.index), len(mutations)), dtype=np.uint8)
        mutation_num = 0

        for mut in quantity_mutation:
            for ind in quantity_mutation[mut]:
                binary_data[ind, mutation_num] = 1
            mutation_num+= 1

        bin_mutation = pd.DataFrame(binary_data, columns=quantity_mutation).assign(Patient=patients)
        bin_mutation = self.hotnet2_table.merge(bin_mutation, on="Patient",how="left")

        return quantity_mutation, bin_mutation

    def __prune_snvs_table(self, percentage):
        number_of_people = int (len(self.hotnet2_table.index) * percentage)
        min_mutation = min(self.mutations, key=lambda mutation: abs(len(self.mutations[mutation]) - number_of_people))
        min_value = abs(len(self.mutations[min_mutation]) - number_of_people)
        res = [mutation for mutation in self.mutations if abs(len(self.mutations[mutation]) - number_of_people) == min_value]
        index_mutation = self.snvs_table.columns.get_loc(res[-1])

        #output
        features_name = list(self.snvs_table.columns[2: index_mutation+1])
        mutation_table = self.snvs_table[features_name].to_numpy()
        cancer_type = self.snvs_table["CancerType"].to_numpy()
        
        return features_name, mutation_table, cancer_type    
    
    def __CORELS_body(self, percentage:float, number_node:int, policy:str, min_support:float):
        features, X, y = self.__prune_snvs_table(percentage)

        card = len(features) if len(features) < MAX_CARD_SUPPORT else MAX_CARD_SUPPORT
        
        #"rulelist","label","progress", "minor"
        c = CorelsClassifier(max_card=card, n_iter=number_node, verbosity=[], policy=policy, min_support=min_support, c=0.005)
    
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=TEST_SIZE)
        c.fit(X_train, y_train, features=features, prediction_name="It's "+self.cancer_type_one+" cancer")
        
        return float(c.score(X_test, y_test)), c
    
    def __tree_body(self, percentage:float):
        features, X, y = self.__prune_snvs_table(percentage)

        clf = tree.DecisionTreeClassifier()
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=TEST_SIZE)
        
        #Conversion type
        y_train = y_train.astype(np.bool_)
        y_test = y_test.astype(np.bool_)

        clf.fit(X_train,y_train)

        return clf.score(X_test, y_test), clf