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
