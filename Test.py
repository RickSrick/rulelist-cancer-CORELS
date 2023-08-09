from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
from corels import *
import pandas as pd
import numpy as np
import csv
import NotificationManager
import os

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
        filename = self.__set_enviroment(number_node, policy, min_support)
        
        mutation_list = sorted(set(self.mutations.values()), reverse=True)[1:]
        a = min(mutation_list, key=lambda x:abs(x-(int(len(self.hotnet2_table.values)*(percentage)))))
        patient_bound = mutation_list.index(a)
        file = open(filename + ".txt", "a")
        scores = []

        if patient_bound <= 1:
            patient_bound = 2

        for v in range(patient_bound):
            print(v)
            score, model = self.__CORELS_body(mutation_list[patient_bound],int(number_node), policy, float(min_support))
            scores.append(score)
            file.write(str(model.rl())+
                "\nSCORE:"+str(score)+
                "\n"+("="*30)+"\n")
            print("Score",score,"\n")
        file.write("AVG SCORE:" + str(np.average(scores)) + "\n"+("="*30)+"\n")
        plt.plot(np.arange(0,patient_bound),scores)
        plt.savefig(filename + "_plot.png")
        plt.close()
        file.close()
        NotificationManager.telegram_notify(cancer_type_one=self.cancer_type_one, cancer_type_two=self.cancer_type_two, number_node=number_node, policy=policy, percentage=percentage, min_support=min_support, filename=filename)
        return np.average(scores)
    
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

            for line in filter(lambda x: x[0] in self.hotnet2_table['Patient'].values,tsv_file):
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
        c = CorelsClassifier(max_card=card, n_iter=number_node, verbosity=["rulelist","label","progress", "minor"], policy=policy, min_support=min_support, c=0.005)
    
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=TEST_SIZE)

        #train_split = int((1-TEST_SIZE) * X.shape[0])
        #X_train, y_train = X[:train_split],y[:train_split]
        #X_test,y_test = X[train_split:], y[train_split:]

        # Fit the model. Features is a list of the feature names
        c.fit(X_train, y_train, features=features, prediction_name="It's "+self.cancer_type_one+" cancer")
        
        return float(c.score(X_test, y_test)), c