# EndoRex

ENDOsymbiotic REConciliation Software

Takes as input a species tree and a gene tree and outputs a reconciled version of the gene tree.
The leaves of the gene tree must provide the species containing the gene, and the location of the gene (0 = mitochondrial, 1 = nuclear).
A reconciled gene tree is output. 

Example usage : 

> python endorex.py -s "(A, B);" -g "((((A__0__a1, A__1__a2), A__0__a3), A__0__a4), B__1__b1);"

By default, gene losses are not shown.  Thy can be displayed using the --lossmode argument.  

> python endorex.py -s "((A, B)x, C)y;" -g "((A__0__ga1, C__0__gc1)w, B__1__gb1)r;" --lossmode=full

# Arguments

The arguments are as follows : 

>--species_tree SNEWICK, -s SNEWICK \
>                        Species tree in newick format.  Each leaf must have a distinct name. \
> \
>--gene_tree GNEWICK, -g GNEWICK \
>                        Gene tree in newick format. Each leaf must be labeled with 3 parameters separated by a double \
>                        underscore __, as follows : [name] __ [species] __ [Location] where [name] is an arbitrary \
>                        identifier of the gene [species] is the name of the species containing the gene and must be a \
>                        leaf of the species tree, [location] is either 0 or 1 \
> \
>--dupcost DUPCOST, -d DUPCOST \
>                        Duplication cost \
> \
>--losscost LOSSCOST, -l LOSSCOST \
>                        Loss cost \
> \
>--transfercost01 TRANSFERCOST01, -f0 TRANSFERCOST01 \
>                        EGT Transfer cost from 0 to 1 \
> \
>--transfercost10 TRANSFERCOST10, -f1 TRANSFERCOST10 \
>                        EGT Transfer cost from 1 to 0 \
> \
>--transpocost01 TRANSPOCOST01, -p0 TRANSPOCOST01 \
>                        EGT Transposition cost from 0 to 1 \
> \
>--transpocost10 TRANSPOCOST10, -p1 TRANSPOCOST10 \
>                        EGT Transposition cost from 1 to 0
>--lossmode LOSSMODE, -lm LOSSMODE
>                        How losses should be displayed, must be one of [none, partial, full]. none : no info on losses
>                        is shown. partial : internal nodes that appear as losses are shown, but not the lost leaves.
>                        full : leaves are also present.
                 
                        
# Output format 

The gene tree is output is newick format, with an additional label at the internal nodes.
The label has the form 

>[event]__[location]

where event is in ['Spe', 'Dup', 'Egtf'] and location is in [0, 1].  Here Egtf is an EGT Transfer event.
The output may also contain unary nodes labeled Egtp that correspond to EGT Transpositions.

If the loss mode is partial or full, internal nodes labeled L are nodes that needed to be inserted owing to a loss.
If the loss mode is full, a loss leaf is labeled loss_[species_name].  

Note that some losses might be in ancestral species.  If the input species tree does not assign a species name to ancestral nodes, species_name will be an empty string here.

