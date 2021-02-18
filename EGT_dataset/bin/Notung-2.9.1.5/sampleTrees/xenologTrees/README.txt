This folder includes trees for all the figures used in the manuscript
"Xenolog Classification."

Reconcile all gene and species trees using the following settings, unless noted:
  - Infer Transfers box ticked
  - Default costs: 
      . Duplications: 1.5 
      . Transfers: 3.0
      . Losses: 1.0
  - Species tree: species.nwk
  - Species Labels: post-fix 
Note: Reference genes in the figures are referred to as gHat

Each gene tree has an associated annotation file.  This annotation
color codes the xenolog classification species sets and the genes that
descended from the transfer recipient.  When imported (File -> Import
Annotations) for a gene tree, all four categories will be color
annotated.  When imported for a species trees, only the D, R, and O
sets will be color annotated.
  - D: blue
  - R: red
  - O: grey
  - gHat, the transferred genes: purple


GENE TREES:

geneFig1and3.nwk
  - Gene tree used in Figures 1 and 3
  - One duplication (in ancestral species of Y and Z)  
  - One transfer (from species Y to species X)

  - Annotation file: geneFig1and3-annotations.txt  
  - Species tree: species.nwk


geneFig2.nwk
  - Gene tree used in Figure 2
  - One transfer (from species Y to species X)

  - Annotation file: geneFig2-annotations.txt  
  - Species tree: species.nwk


geneFig4.nwk
  - Gene tree used in Figure 4
  - Two comparable transfers 
      from species Y to species (X,U)
      from species X to species W
  - One loss in ancestral species of (X,U)

  - Annotation file: geneFig4-annotations.txt  
      This annotation file color codes sets for both transfer t2 and
      the supertransfer.  The displayed annotations can be changed
      from the Annotations tab, using the "Show/Hide" button
  - Species tree: speciesFig4.nwk


geneFig5.nwk
  - Gene tree used in Figure 5
  - Two incomparable transfers 
      from species (Y,Z) to species V
      from species Y to species X

  - Annotation file: geneFig5-annotations.txt  
      This annotation file color codes sets for both transfers t2 and
      t2.  The displayed annotations can be changed from the
      Annotations tab, using the "Show/Hide" button
  - Species tree: species.nwk


geneFigS1.nwk
  - Gene tree used in Figure S1
  - One duplication (in ancestral species of W, X, Y, and Z)  
  - One transfer (from species Y to species X)

  - Annotation file: geneFigS1-annotations.txt  
  - Species tree: species.nwk


geneFigS2.nwk
  - For reconciliation, use Loss cost of 2.0 
  - Gene tree based on Figure S2
  - Two comparable transfers, form a loop in the species tree
      from species (W,X) to species (Y,Z)
      from species Y to species X

  - Annotation file: geneFigS2-annotations.txt  
      This annotation file color codes sets for both transfer t2 and
      the supertransfer.  The displayed annotations can be changed
      from the Annotations tab, using the "Show/Hide" button
  - Species tree: speciesFig.nwk


BIO4 EXAMPLE:
Figures 7, S6, S7, and S8

Gene Tree: bio4g.nwk
Species Tree: bio4s-annotations.ntg
* These trees were transcribed from:
  Hall and Dietrich. (2007). The reacquisition of biotin prototrophy
    in Saccharomyces cerevisiae involved horizontal gene transfer,
    gene duplication and gene clustering.  Genetics, 177(4),
    2293–2307.
Annotations: bio4Annotations.txt
Reconciled Gene Tree: bio4g-reconciledAnnotated.ntg

For reconciliation, use the following settings:
  - Infer Transfers box ticked
  - Costs: 
      . Duplications: 1.5 
      . Transfers: 5.0
      . Losses: 1.0
  - Species tree: bio4s.nwk
  - Species Labels: post-fix 
  - To hide losses, uncheck the box "Display Loss Nodes"
Note: Reference gene was bio4_sacce

When reconciled, there will be 2 equally optimal solutions; the
difference between these two solutions is the source of the transfer
into gammaproteobacteria and the associated loss.  To cycle through
the solutions, click on the green circle.

