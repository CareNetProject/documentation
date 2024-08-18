import sys
import os
import csv
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import re

#takes a string and returns whole words (length > 1) without special characters
#sorted and separated by '-'
#'John K Smith, MD' would return 'john-md-smith'
def extractName(text):
    #strip special characters
    text  = re.sub('[^a-zA-Z0-9 \n]', '', text)
    #set to lowercase, split into a list, and sort
    wordList = sorted(text.lower().split())
    #strip out length 1 strings (initials)
    wordList = [x for x in wordList if len(x) > 1]

    #join with '-'
    return '-'.join(wordList)

#takes the first part of a filename as an argument
#example: use no_re_ward_events to catch no_re_ward_events1.csv, no_re_ward_events2.csv, etc
if len(sys.argv) > 1:
    fileStarts = sys.argv[1]  
else:
    print("No filename provided.")
    exit()

df_full = pd.DataFrame()
for filename in os.listdir('input/'):
    if filename.startswith(fileStarts) and 'output' not in filename:
        print(f'Adding: {filename}')
        df2 = pd.read_csv('input/' +filename)
        df_full = pd.concat([df_full, df2])
df_full = df_full.reset_index(drop=True)

#drop null rows for accessor/accessee
df_full = df_full.dropna(subset=['DEID_ACCESS_USER_NAME'])
df_full = df_full.dropna(subset=['DEID_NOTE_AUTHOR_NAME'])

#remove any non-human names listed in exclude.txt
with open('exclude.txt', 'r') as file:
    for line in file:
        line = line.rstrip('\n')
        df_full = df_full[~df_full['DEID_ACCESS_USER_NAME'].str.contains(line)]
        df_full = df_full[~df_full['DEID_NOTE_AUTHOR_NAME'].str.contains(line)]


#group by patientID
grouped_df = df_full.groupby('DEID_PAT_MRN_ID')

for id_value, df in grouped_df:

    #create a graph for each patient
    g = nx.DiGraph()
    df['DEID_ACCESS_TIME'] = pd.to_datetime(df['DEID_ACCESS_TIME'], format='%m/%d/%y %H:%M')
    df = df.sort_values(by=['DEID_ACCESS_TIME'], ascending=True)

    associatedTimes = {}
    for i in df.index:
        fromName = extractName(df['DEID_ACCESS_USER_NAME'][i])
        toName = extractName(df['DEID_NOTE_AUTHOR_NAME'][i])

        #exclude when from and to are same person
        if (fromName != toName):
            #if this edge already exists, see if we need to increment the weight
            if g.has_edge(str(fromName), str(toName)):
                comboName = str(fromName) + str(toName)

                if comboName in associatedTimes:
                    timeDiff = (df['DEID_ACCESS_TIME'][i] - associatedTimes[comboName]).total_seconds()
                    # 30 minutes
                    if (timeDiff > 1800):
                        g[fromName][toName]['weight'] += 1
                        associatedTimes[comboName] = df['DEID_ACCESS_TIME'][i]
                else:
                    g[fromName][toName]['weight'] += 1
                    associatedTimes[comboName] = df['DEID_ACCESS_TIME'][i]

            #if edge doesn't already exist, make the edge with weight 1
            else:
                g.add_weighted_edges_from([(str(fromName), str(toName), int(1))], weight='weight')
                comboName = str(fromName) + str(toName)
                associatedTimes[comboName] = df['DEID_ACCESS_TIME'][i]


    #save location depending on name (ward vs icu)
    if ('ward' in fileStarts):
        outFile = str(id_value) + '_ward.csv'
    else:
        outFile = str(id_value) + '_icu.csv'
    file = open('graphs/' + outFile, 'w+', newline ='')
    with file:
        write = csv.writer(file)
        for edge in g.edges(data=True):
            source, target, data = edge
            weight = data['weight']
            write.writerow([source, target, weight])
    file.close()
