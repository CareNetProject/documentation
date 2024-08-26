# Steps

1. Install requirements. `pip3 install -r requirements.txt`
2. Generate graphs from the input files. Graphs will be saved in `/graphs`.
   ```
   python3 createGraphs.py no_re_icu_events
   python3 createGraphs.py no_re_ward_events
   ```
3. Calculate statistics from the graphs. `python3 calc.py`
   - Results will be saved in `/results`.
   - Each patient will have an individual file, along with two aggregate files.  One using NetworkX calculations, the other using manually coded formulas for verification.

## Creating a graph visualization

1. Steps 1 and 2 from above must be complete.  Visualizations are created from the .csv files saved in `/graphs`.
2. Run `python3 viewGraph.py graphs/<filename>.csv` where the filename is the graph you want to visualize.
   -  The type of graph can be configured by altering `viewGraph.py`.  Other layouts have been left in the code.  Leave the layout you want to use uncommented and comment out the rest.  The default layout is `pos = nx.fruchterman_reingold_layout(g,  scale = 1, k=1 )`.
