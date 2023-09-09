from sklearn.model_selection import train_test_split
from sklearn import tree
import matplotlib.pyplot as plt
from corels import *
import pandas as pd
import numpy as np
import csv
import NotificationManager
import os
import time

MAX_CARD_SUPPORT = 2
TEST_SIZE = 0.2

class Test:

    def __init__(self, files:list, cancer_type_one:str, cancer_type_two:str)-> None:
        self.hotnet2, self.snvs_strictly_filtered = files
        self.cancer_type_one = cancer_type_one
        self.cancer_type_two = cancer_type_two
        self.hotnet2_table = self.__get_hotnet2_table()
        self.mutations = self.__save_unique_mutation()

    def launch_test(self, percentage:float, number_node:int, policy:str, min_support:float)-> float: 
        #set output enviroment
        filename = self.__set_enviroment(number_node, policy, min_support)
        file = open(filename + ".txt", "a")

        
        
        mutation_list = sorted(set(self.mutations.values()), reverse=True)[1:]
        a = min(mutation_list, key=lambda x:abs(x-(int(len(self.hotnet2_table.values)*(percentage)))))
        patient_bound = mutation_list.index(a)
        scores = []

        #check valid patient bound 
        if patient_bound <= 1:
            patient_bound = 2

        if os.path.isfile("data.csv") == False:
            data_file = open("data.csv", "a")
            data_file.write("testname,percentage,nodes,policy,min_support,time,rule_score,tree_score\n")
        else:
            data_file = open("data.csv","a")
        #for v in range(patient_bound):
        #    print(v)
        start_time = time.perf_counter()
        score, model = self.__CORELS_body(mutation_list[patient_bound],int(number_node), policy, float(min_support))
        end_time = time.perf_counter()
        execution_time = end_time - start_time
        b = self.__trie_body(mutation_list[patient_bound])
        print("Score",score,"\n")
        data_file.write(f'''{self.cancer_type_one} o {self.cancer_type_two},{percentage},{number_node},{policy},{min_support},{round(execution_time, 3)},{round(score, 3)},{round(b,3)}\n''')

        data_file.close()   
        file.write(str(model.rl())+
            "\nSCORE:"+str(score)+
            "\n"+("="*30)+"\n")
        file.close()
        
        NotificationManager.telegram_notify(cancer_type_one=self.cancer_type_one, cancer_type_two=self.cancer_type_two, number_node=number_node, policy=policy, percentage=percentage, min_support=min_support, filename=filename)
        return 0
    
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
    
    def __save_unique_mutation(self)-> dict:
        mutations = {}
        with open(self.snvs_strictly_filtered) as file:
            tsv_file = csv.reader(file, delimiter="\t")

            for line in filter(lambda x: x[0] in self.hotnet2_table['Patient'].values, tsv_file):
                    for x in line[1:]:
                        mutations[x] = mutations.get(x,0) + 1
        return mutations

    def __get_snvs_strictly_filtered_table(self, mutation: list):
        binary_data = np.zeros(shape=(len(self.hotnet2_table.index), len(mutation)), dtype=np.uint8)
        patient_list_order = []
        nline = 0

        with open(self.snvs_strictly_filtered) as file:
            tsv_file = csv.reader(file, delimiter="\t")

            for line in filter(lambda x: x[0] in self.hotnet2_table['Patient'].values,tsv_file):
                patient_list_order.append(line[0])
                for x in line[1:]:
                    if x in mutation:
                        binary_data[nline, mutation.index(x)] = 1
                nline += 1

        bin_mutation = pd.DataFrame(binary_data, columns=mutation).assign(Patient=patient_list_order)
        bin_mutation = self.hotnet2_table.merge(bin_mutation, on="Patient",how="left")

        return bin_mutation[mutation].to_numpy(), bin_mutation["CancerType"].to_numpy()

    def __CORELS_body(self, bound:int, number_node:int, policy:str, min_support:float):
        features = [i for i in self.mutations if (self.mutations[i] > bound)]
        X, y = self.__get_snvs_strictly_filtered_table(features)

        card = len(features) if len(features) < MAX_CARD_SUPPORT else MAX_CARD_SUPPORT
        #"rulelist","label","progress", "minor"
        c = CorelsClassifier(max_card=card, n_iter=number_node, verbosity=[], policy=policy, min_support=min_support, c=0.005)
    
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=TEST_SIZE)
        c.fit(X_train, y_train, features=features, prediction_name="It's "+self.cancer_type_one+" cancer")
        
        return float(c.score(X_test, y_test)), c
    
    def __trie_body(self, bound:int):
        features = [i for i in self.mutations if (self.mutations[i] > bound)]
        X, y = self.__get_snvs_strictly_filtered_table(features)
        clf = tree.DecisionTreeClassifier()
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=TEST_SIZE)
        y_train = y_train.astype(np.bool_)
        y_test = y_test.astype(np.bool_)
        X_test = X_test.astype(np.bool_)
        X_test = X_test.astype(np.bool_)
        clf.fit(X_train,y_train)
        plt.figure(figsize=(25,25))
        a = clf.score(X_test, y_test)
        tree.plot_tree(clf, feature_names=features)
        plt.savefig('tree.png',format='png')
        NotificationManager.send_photo("tree.png")
        return a
        