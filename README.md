# StableTBN
A Python tool to generate the stable configurations of a Thermodynamic Binding Network.  Also supports user-specified constraints for more specific investigation of a network.  


#### System Requirements
+ Python 3.6 (or higher)
+ [Google OR-tools](https://developers.google.com/optimization/install)

In order to use functions that compute Hilbert bases, the package [4ti2](https://4ti2.github.io/) must also be installed and accessible via the system path.

Note: This software is tested on Linux-based operating systems, but we expect that most functionality will be retained on Windows and macOS-based systems.  Please let us know if you find any compatibility issues!
  
# Command line usage
    
    $ python3 stable-tbn.py <tbn_filename>

## Optional Command line flags 
    
    -h, --help            show help info
    -i, --ip              use IP for the optimization step (instead of CP)
    -1                    report only one stable configuration
    -t, --timed           print elapsed time
    -f, --full            print full configuration (includes singletons)
    -w <bond_weight>      relative worth of bonds vs polymers formed, e.g. 0.5
    -c <constraints_file> filename for constraints text file
    --formulation <name>  specify a specific solution formulation, one of:
                             BOND_AWARE_NETWORK,
                             BOND_OBLIVIOUS_NETWORK,
                             POLYMER_BINARY_MATRIX,
                             POLYMER_INTEGER_MATRIX,
                             POLYMER_UNBOUNDED_MATRIX,
                             VARIABLE_BOND_WEIGHT,
                             GRAVER_BASIS
    -v, --verbose         display solver output
    --benchmark           do not display the stable configuration(s)



## Input file format (TBN specification)

The input file that specifies a TBN must have the following format:

#### Binding Sites
- Binding sites can be any alphanumeric string (underscores are allowed): "a"
- Any binding site that ends with an asterisk (\*) indicates a complement binding site: "a\*"
- Multiple binding sites of the same type can be specified by a number followed by the binding site in parenthesis: "3(a)"

#### Monomers
- Each line of the input file represents a monomer (with whitespace separating the binding sites): "a b c d\*"
- Ending the line with a right-angle bracket (>) followed by a string specifies an optional name: "a b c d\* >m1" (alphanumeric with optional underscores)
- To specify multiple monomers of the same type, specify the quantity and follow with the monomer (and its optional name) in brackets: "2\[a b c d\* >m1\]"
- To specify an excess of monomers of the same type, use the string "inf" as the quantity: "inf\[a b c d* >m1\]" 

Note: This format is backwards-compatible with the input format of [StableGen](http://stablegen.net/).  However, because StableTBN does not implement site-level bonds, any site-level names will be ignored.  

## Constraints file format

In the optional constraints file, you can provide constraints to specify additional properties.  Each line of the file specifies an additional constraint.

#### MAX/MIN POLYMERS

Specifying **MAX/MIN POLYMERS** instructs the solver to only give the configurations with number of components bounded by the provided values.  Optionally, these constraints can be used together to specify both bounds.

    MAX POLYMERS 5
    MIN POLYMERS 3
    
#### MAX/MIN MERGES

Specifying **MAX/MIN MERGES** instructs the solver to only give the configurations that perform the specified number of merges.

    MAX MERGES 10
    MIN MERGES 9
    
#### MAX/MIN ENERGY

Specifying **MAX/MIN ENERGY** instructs the solver to only give the configurations that are bounded in energy by the specified values.  To use this, must provide a relative weighting factor using BOND_WEIGHT or by using the "-w" flag.

    MAX ENERGY 4.6
    MIN ENERGY 2.0
    
#### BOND WEIGHT

Specifying **BOND WEIGHT** forces the solver to use LOW_W_FORMULATION and specifies the relative energetic benefit of forming site-level bonds over establishes separate components.  The energy of a configuration is calculated by (weighting_factor * bond_deficit + number_of_merges), where bond_deficit is the difference between the number of bonds in a saturated configuration and the present configuration.

    MAX ENERGY 4.6
    MIN ENERGY 2.0

#### NO OPTIMIZE

Specifying **NO OPTIMIZE** instructs the solver to bypass the optimization step.  It will report all saturated configurations of the network that satisfy the constraints.  CAUTION: This may be a great many configurations, so be sure to specify more constraints!

    NO OPTIMIZE
    
#### NO SORT

Specifying **NO SORT** instructs the solver to not sort the column order of the matrix received from the black box solver.  It is possible that this may provide a small speedup if only one configuration is required.  It is not recommended to use this when solving for all configurations.

    NO SORT

#### TOGETHER (planned, NYI)

Specifying **TOGETHER** forces a multiset of monomers to group as a motif into the same component.  A global value can also be specified which determines how many times this motif must be repeated.  

    #[TOGETHER #[m1] #[m2] #[m3] ...]

#### NOTTOGETHER (planned, NYI)

Specifying **NOTTOGETHER** prevents a multiset of monomers from grouping into the same component.

    NOTTOGETHER #[m1] #[m2] #[m3] ...

#### FREE (planned, NYI)

Specifying **FREE** forces the specified monomer to be in an explicit singleton component (or if a value is specified, a lower bound on the number of this monomer that must be free)

    FREE #[m1]

#### NOTFREE (planned, NYI)

Specifying **NOTFREE** forces specifies a lower bound on the number of a monomer type that must be present in non-singleton components.

    NOTFREE #[m1]
    
Note: StableTBN does not implement site-level bonds explicitly, and so StableGen-style **PAIRED** constraints are not supported. 

# Examples

#### Example 1 input

    examples/tbn.txt

    a b*
    a* b
    a
    a* b
    b
    b


#### Example 1 output
    
    $ python3 stable_tbn.py examples/tbn.txt

    Configuration 1:
    {[a b*], [a* b]}; {[a* b], [a]}
    
    $ python3 stable_tbn.py examples/tbn.txt --full

    Configuration 1:
    {[a b*], [a* b]}; {[a* b], [a]}; 2{[b]}


#### Example 2 input

    examples/tbn_gg.txt

    x00 x01 x02 >H0
    x10 x11 x12 >H1
    x20 x21 x22 >H2
    x00 x10 x20 >V0
    x01 x11 x21 >V1
    x02 x12 x22 >V2
    x00* x01* x02* x10* x11* x12* x20* x21* x22* >G
    
#### Example 2 output
    
    $ python3 stable_tbn.py examples/tbn_gg.txt

    Configuration 1:
    {G, H0, H1, H2}
    Configuration 2:
    {G, V0, V1, V2}

#### Example 3 input

    examples/tbn_gg_2g_with_cat_and_excess.txt

    inf[x00 x01 x02 >H0]
    inf[x10 x11 x12 >H1]
    inf[x20 x21 x22 >H2]
    inf[x00 x10 x20 >V0]
    inf[x01 x11 x21 >V1]
    inf[x02 x12 x22 >V2]
    2[x00* x01* x02* x10* x11* x12* x20* x21* x22* >G]
    1[x00 x10 x11 x20 x21 x22 >C]

#### Example 3 output

    $ python3 stable_tbn.py examples/tbn_gg_2g_with_cat_and_excess.txt

    Configuration 1:
    {G, H0, H1, H2}; {G, V0, V1, V2}
    Configuration 2:
    2{G, H0, H1, H2}
    Configuration 3:
    {C, G, H0, H1}; {G, V0, V1, V2}
    Configuration 4:
    {C, G, H0, H1}; {G, H0, H1, H2}
    Configuration 5:
    {C, G, H0, V2}; {G, V0, V1, V2}
    Configuration 6:
    {C, G, H0, V2}; {G, H0, H1, H2}
    Configuration 7:
    {C, G, V1, V2}; {G, H0, H1, H2}
    Configuration 8:
    {C, G, V1, V2}; {G, V0, V1, V2}
    Configuration 9:
    2{G, V0, V1, V2}

# Citation

#### Author
David Haley

#### References
David Haley, David Doty, "Computing Properties of Thermodynamic Binding Networks: An Integer Programming Approach". [https://arxiv.org/abs/2011.10677](https://arxiv.org/abs/2011.10677) 

Laurent Perron and Vincent Furnon. OR-tools. [https://developers.google.com/optimization](https://developers.google.com/optimization)

4ti2 team. 4ti2 – a software package for algebraic, geometric and combinatorial problems on linear spaces. [https://4ti2.github.io](https://4ti2.github.io)
