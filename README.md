# 1. Ladderpath JSON Standard Format
## V1.0.0.20240910\_Alpha Version

```yaml
{
"info": "V1.0.0.20240910_Alpha"

"ladderpath-index": 23,
"order-index": 10, # order-index Ω
"size-index": 33, 
"eta": 0.67, // order-rate (η)
"omega_max": 15, // Ω of the most ordered sequenece with the same size-index, used to compute η
"omega_min_list": [3,5,5,9,3], // List of Ω's from randomized sequences; all outcomes are recorded here. The minimum is used to compute η.
      
"ladderons": {  // Dictionary in the format ID (int): List[COMP (List), LEN (int), STR (str), POS (Dict)]
    0: [ //ID  0, 1, 2, ...
     ["ACE", 1, "C", 5, 6, 6, 3, "EE"], // COMP: Basic building blocks + ladderon IDs arranged by the target sequence; uniquely reconstructible
     23, // LEN: Length of the ladderon
     "", // STR: Placeholder string, used to store the reconstructed ladderon
     {-1: [1, 20, 35], 1: [0, 6], 4: [2]} // POS: {int: [int]}, indicates which higher-level ladderons (by ID) use this ladderon, and the positions (indexes) at which it appears in their COMP lists.
     ],
   
    1: [
     [6, 2],
     13, 
     "",
     {3: [21, 30]}
     ],
   
     ...
    }


"basic_building_blocks": List(str)

  
"targets": {  // Dictionary in format ID (int): [COMP (List), LEN (int), STR (str), REP (int)]
    -1: [ // Target identifiers: -1, -2, -3, ...
      ["ACE", 1, "C", 5, 6, 6, 3, "EE"],
      98,
      "",
      1  // REP: Number of times this target appears in the original targets
      ],
   
    -2: [ 
      [1, 3, "A", 6],
      50,
      "",
      2
      ],
   
    ...
    }
  
"duplications_info": Dict(int:[int]) // Stores information about duplicated sequences in targets. Key: target ID; Value: list of its positions in the original targets (0-based indexing)
}
```

Notes:

* `COMP`: ["ACE", 1, "C", 5, 6, 6, 3, "EE"] indicates that this ladderon corresponds to ACEABCRRDADADDDEE (when segmented, it reads: ACE AB C RR DA DA DDD EE). Here, AB has an assigned ID of 1, so it is represented as 1 in COMP; RR has an ID of 5, hence 5 in COMP; The same logic applies to all other components.

* Targets must be sequentially numbered as -1, -2, -3, etc.; Ladderons are numbered 0, 1, 2, etc.; The basic building blocks (bbb) are numbered b1, b2, b3, etc.



## V1.0.1.20240928\_Alpha Version

All else remains unchanged, except for modifications to the handling of eta. Specifically: (1) The separate entries `omega_max` and `omega_min_list` have been removed. (2) A new consolidated information object called eta_info has been added.

```JSON
{
"eta": 0.67, // Order-rate (η), representing an approximate value. Detailed information is stored in `eta_info`


"eta_info": {
    "omega_max_AllIdentical": 20, // Equivalent to generating a sequence of the same length as the original, consisting entirely of 'A's (no longer segmented by target lengths).
    "omega_max_Sorted": 15, // Calculated omega from the fully concatenated and sorted version of the original sequence (also not segmented by target lengths). This is often the maximum omega among all permutations.
  
    "omega_min_Shuffle_list": [3,5,5,9,3], // Omega values obtained by random shuffling of the target sequence. The minimum of these is used in computing eta.
    "omega_min_LocalDist_list": [3,5,8] // Omega values from randomly generated sequences that follow the original letter distribution. These sequences have the same length as the target.
    "omega_min_EvenDist_list": [5,9,2,6] // Omega values from sequences of the same length as the original, with all characters evenly distributed.
    }
}
```

## V1.0.2.20250404\_Alpha Version

An additional input format specification for targets was introduced:

```yaml
{
"input_type": "list" # or "dict"
}
```






# 2. Usage
## 2.1 Main functions: Calculate Ladderpath

This is the main function for computing ladderpaths (the most time-consuming part):

```Python
import ladderpath as lp

# strs = ['ABCDBCDBCDCDEFEF']
strs = ['ABCDBCDBCDCDEFEF', 'AA', 'DBCDBCA', 'BCDCDAEF']
# strs = ['ABCDBCDBCDCDEFEF', 'AA', 'AA', 'DBCDBCA', 'BCDCDAEF','AA', 'DBCDBCA', 'AAAA']
# strs = {'ABCDBCDBCDCDEFEF':1, 'AA':3, 'DBCDBCA':2, 'BCDCDAEF':1, 'AAAA':1}

lpjson = lp.get_ladderpath(strs,
           info='V1.0.0.20240910_Alpha', 
           fill_ladderons_STR = False,
           estimate_eta = False, 
           save_file_name=None, 
           show_version=True) # Output is in the standard JSON format for ladderpaths
```

* `strs`: `List[str]`. A list of target sequences to compute the ladderpath for. It can contain a single string (i.e., compute the ladderpath for just one sequence). It can also contain multiple strings. Duplicate strings are permitted and handled correctly. Alternatively, a dict format can be used, where keys are strings and values indicate frequency.
* `info`: `str`. Metadata to be recorded in the output JSON (e.g., program version, data source, author, etc.).
* `fill_ladderons_STR`: `bool` (default: `False`). Whether to fully populate the string content of ladderons in the output.
* `estimate_eta`: `bool` (default: `False`). Whether to estimate the order-rate eta.
* `estimate_eta_para`: Only required if `estimate_eta=True`. Must be explicitly set. See the next section for details.
* `save_file_name`: `str`. Path and filename for saving the output JSON file.
* `show_version`: `bool` (default: `True`). Whether to print version and metadata (contained in the info field).



The output is a JSON file in the standard ladderpath format (see Section 1 for specifications).

```Python
{'info': 'V1.0.1.20240928_Alpha',

 'ladderpath-index': 20,
 'order-index': 13,
 'size-index': 33,
 
 'eta': None,
 
 'ladderons': {0: [[1, 1], 6, '', {-1: [3], -3: [0]}],
               1: [['D', 4], 3, '', {0: [0, 3]}],
               2: [['DCD'], 3, '', {-4: [2], -1: [9]}],
               3: [['EF'], 2, '', {-1: [12, 14], -4: [6]}],
               4: [['BC'], 2, '', {-1: [1], -4: [0], 1: [1]}]},
 'basic_building_blocks': ['A', 'D', 'C', 'E', 'F', 'B'],
 
 'targets': {-1: [['A', 4, 0, 2, 3, 3], 16, '', 1],
             -2: [['AA'], 2, '', 1],
             -3: [[0, 'A'], 7, '', 1],
             -4: [[4, 2, 'A', 3], 8, '', 1]},
           
 'duplications_info': {}
 
 "eta_info": {
    "omega_max_AllIdentical": None,
    "omega_max_Sorted": None,
    "omega_min_Shuffle_list": [],
    "omega_min_LocalDist_list": [],
    "omega_min_EvenDist_list": []
    }
}
```

## 2.2 Auxiliary Functions of the Main
### Load standard Ladderpath JSON file

```Python
lpjson = lp.load_ladderpath_json(save_file_name)
```

* `save_file_name` is the file path and name of the JSON file, e.g., `'new_lpjson.json'`.

* This function will automatically convert all dictionary keys in the JSON file to integers.



### Clear and Populate the Cache

```Python
lp.clear_lpjson_STR(lpjson) # Clears the explicit string representation (STR) of each ladderon in the JSON (clears the cache)

lp.fill_lpjson_STR(lpjson) # Fills in the explicit string representation (STR) of each ladderon in the JSON (populates the cache)
```

* These functions modify lpjson in place—no reassignment is necessary.



### Display the Partially Ordered Multiset (POM)

```Python
# Display three metrics: ladderpath-index, order-index, size-index
index3 = lp.disp3index(lpjson) 

# Obtain the partially ordered multiset (POM) representation of the ladderpath
pom, pom_str = lp.POM_from_JSON(lpjson, display_str=False)
```

* `pom` is the dictionary representation of the ladderpath's partially ordered multiset (POM):

```Python
{0: {'A': 3, 'D': 3, 'C': 2, 'B': 1, 'E': 1, 'F': 1}, # Level 0
 1: {'EF': 2, 'BC': 2, 'DCD': 1}, # Level 1
 2: {'DBC': 1}, # Level 2
 3: {'DBCDBC': 1}} # Level 3
```

* `pom_str` is the string-based representation of the ladderpath POM, expressed in a format like: `{ A, B // CD(2), // ACD }`. This provides a compact view of the structural hierarchy within the ladderpath.




### Compute order-rate eta (η)

```Python
estimate_eta_para={'n': 10, 'max_method': 'AllIdentical', \
                   'min_method': 'EvenDist', \
                   'min_method_nBase': 26}
eta = lp.get_eta(lpjson, estimate_eta_para=estimate_eta_para)
```

* `estimate_eta_para` is a required dictionary of parameters for computing `eta`. These must be explicitly set.
(Note: Currently, when there are multiple targets, the sequences generated for computing omega_min or omega_max are not split according to individual targets—the sequences are merged into one long sequence before computation.)
  - `n`: Number of randomized sequences to generate when computing `omega_min`
  
  - `max_method`: Method used to define the theoretical maximum for omega: `"AllIdentical"`- assumes that a sequence of repeated 'A's has the highest order (most structured); `"Sorted"`- assumes that sorting the characters of the original sequence yields the most ordered form (less accurate, generally not recommended).
  
  - `min_method`: `"Shuffle"`- randomly shuffles the original sequence to generate a distribution; `"EvenDist"`- generates random sequences using an even distribution over nBase distinct letters; `"LocalDist"`- intended to sample based on the empirical letter distribution of the target sequence (not supported in current version).

  * `min_method_nBase`: Required when `min_method='EvenDist'`; indicates the number of distinct basic letters to distribute evenly.

If an eta value already exists in lpjson (and was computed using the same specified method), it will be returned directly. Otherwise, it is computed anew.

```Python
lp.update_eta(lpjson, estimate_eta_para=estimate_eta_para)
# lp.get_eta(lpjson, estimate_eta_para=estimate_eta_para)
```

This function updates the eta value in-place within lpjson, using the parameters specified in `estimate_eta_para`. It performs the following steps: 1. Generates `n` randomized sequences based on the `min_method`. 2. Computes `omega_min` for each, and appends the results to the corresponding `omega_min_list`. 3. Takes the minimum value from that list as `omega_min` and uses it to compute `eta`. 4. Updates the `eta` field and relevant subfields in `lpjson`.

Note: `update_eta()` does not return a value. To retrieve the computed eta, use `lp.get_eta()` afterward.




### Draw Laddergraph

```Python
g = lp.draw_laddergraph(lpjson, 
                        show_longer_than = 0, 
                        style = "ellipse",
                        TargetNames_UserDefined = {-1:'species1', -2:'new2', -3:'species3'},
                        warning_n_ladderons_to_show = 500, 
                        rankdir = "BT", 
                        color = "grey",
                        figsize = "8,8",
                        save_fig_name = None, 
                        figformat = "pdf", 
                        cleanGVfile = True)
```

* `show_longer_than`: `int`. Display only ladderons whose sequence length exceeds this value. If set to 0, all ladderons are displayed, including the most basic units.

* `style`: `str`. Determines the visual style of the laddergraph. Options include: `"ellipse"`: The sequence is not shown; the ellipse size reflects sequence length. `"ellipse-OnlyShowTargetID"`: Same as "ellipse", but limits ellipse size for target sequences. `"box"`: Displays the sequence within each node (box). `"box-OnlyShowTargetID"`: Displays only the target's ID, not the full sequence.

* `TargetNames_UserDefined`: `dict`. A user-defined mapping of target IDs to human-readable names.

* `warning_n_ladderons_to_show`: `int`. Maximum number of ladderons allowed in the visual output to prevent rendering issues.

* `rankdir`: `str`. Node layout direction. Acceptable values: BT (from bottom to top), TB, LR, RL.

* `color`: `str`. Accepts color names (e.g., `"grey"`, `"red"`) or hex codes (e.g., `"#808080"`).

* `figsize`: `str`. Size of the output figure. Default: None (automatically adjusted). Suggested: `"8,8"` (8 inches × 8 inches).

* `save_fig_name`: `str`. File name (with path) for saving the figure, e.g., `"figs/New\_Laddergraph"`.

* `figformat`: `str`. File format for the saved figure. Common options: `"pdf"`, `"png"`, etc.

* `cleanGVfile`: `True/False`. Whether to delete the intermediate `.gv` file generated during rendering.



### Additional Auxiliary Functions
**Deduplicate Target Sequences**

```Python
strs_unique, duplications_info = uniquenize(strs)
```

**Reconstruct the Target Sequence**
(This is useful when there are duplicate sequences among the targets, as the target IDs may not match their order of appearance.)

```Python
reconstruct_targets_list(lpjson) -> List[int]
```

**Count basic_building_blocks**, e.g., `{'A': 5, 'C': 3, 'B':1}`

```Python
bbb_dict = bbb_multiplicity(lpjson) -> Dict[str:int]
```




## 2.3 Additional Modules and Applications
### Compression

```Python
import ladderpath_tools.compress as lp_c
```

**Compress a standard ladderpath JSON into a compact list**

```Python
# strs = ['ABCDBCDBCDCDEFEF', 'AA', 'AA', 'DBCDBCA', 'BCDCDAEF','AA', 'DBCDBCA', 'AAAA']
# lpjson = lp.get_ladderpath(strs)

compressed_list = lp_c.compress(lpjson, display=False, SEP='@')
```

* The `compressed_list` contains the compressed representation of the ladderpath. Explanation of the structure:

```JSON
# [8, // Number of target sequences
    '@', 'A', 1, 6, 4, 2, 2,  // Target 1
    '@', 3,  // Target 2
    '@', 3, 
    '@', 0, 
    '@', 1, 4, 'A', 2, 
    '@', 3, 
    '@', 0, 
    '@', 3, 3,  // Target 8
  
    '@', 6, 'A',  // Ladderon 0 (note: ordering may differ from lpjson)
    '@', 'BC',   // Ladderon 1
    '@', 'EF',   // Ladderon 2
    '@', 'AA',   // Ladderon 3
    '@', 'DCD',  // Ladderon 4
    '@', 'D', 1, // Ladderon 5
    '@', 5, 5]   // Ladderon 6
```

**Decompress the `compressed_list` back into the original targets**
This correctly restores the original list of sequences, including any duplicates.

```Python
decompressed_str = lp_c.decompressed(compressed_list, display=False)
```


### Compute Ladderpath Distance from a Laddergraph

```Python
import ladderpath_tools.lp_dis_laddergraph2points as lp_dis_laddergraph2points

two_leaf = ['ABCDBCDBCDCDEFEF', 'DBCDBCA']
# two_leaf = [-1, -3]
lp_dis_laddergraph2points.lp_dis_laddergraph2points(two_leaf, lpjson)
```

* `two_leaf` can either be the two target sequences themselves or their IDs. The `lpjson` must be the ladderpath JSON that contains both targets.



### Phylogeny & Distance

```Python
import ladderpath_tools.phylo as lp_p
```

**Generate a Distance Matrix** (Creates an instance of `LADDERPATH_DISTANCE_MATRIX`)

```Python
# Deduplicate target sequences first
strs_unique, duplications_info = lp.uniquenize(strs)

lp_dis_mat = lp_p.LADDERPATH_DISTANCE_MATRIX(strs_unique, strsNames=None， method='lp_dis_laddergraph2points-absolute', lpjson=lpjson)
```

* Deduplication is required: The list of input sequences `strs` must have duplicates removed using `lp.uniquenize(strs)`.

* The `LADDERPATH_DISTANCE_MATRIX` object has only three attributes: `.strs`: the deduplicated list of sequences; `.strsNames`: user-defined or automatically assigned names; `.mat`: the n × n upper triangular matrix of ladderpath distances, where `mat[i, j]` is the distance between target i and j.

* Parameters:
  - `strs_unique`: `List[str]`. The deduplicated list of target sequences.
  - `strsNames`: `List[str]` or `None`. Optional list of names (e.g., species names) corresponding to strs_unique. If `None`, names will be assigned in the order of appearance (e.g., `"0"`, `"1"`, ...).
  - `method`: `str`. Must be one of the following: `'lp_dis_laddergraph2points-absolute'`, `'lp_dis_laddergraph2points-relative'`,`'independent2points-absolute'`, `'independent2points-relative'`. These represent different strategies: `"laddergraph2points"` methods compute distances directly from a common laddergraph; `"independent2points"` methods compute distances by constructing ladderpaths independently for each pair (differences are usually minor).
  - `lpjson`: The ladderpath JSON object computed from `strs_unique`.



**Plot Phylogenetic Trees: UPGMA & NJ**

```Python
# Generate a UPGMA phylogenetic tree and export in Newick format
lp_dis_mat.phylotree_upgma_newick(outgroupName=None, save_newick_file=None, save_fig_pdf_name=None) -> None

# Generate a Neighbor-Joining (NJ) tree and export in Newick format
lp_dis_mat.phylotree_nj_newick(outgroupName=None, save_newick_file=None, save_fig_pdf_name=None) -> None
```

* `outgroupName`: Optional; specify the name of an outgroup if desired.



**Plot a Minimum Spanning Tree**

```Python
# Plot a minimum spanning tree using the computed distance matrix
lp_dis_mat.minimum_spanning_tree(save_fig_pdf_name=None, show_weights=False) -> None
```

* `show_weights`: `False/True`. Whether to display edge weights (distances) on the MST graph.
