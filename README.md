EasyBadiRate
============

Using badirate to calculate gene families with significant expansion and contraction



How it works
------------
BadiRate is a program that estimates family turnover rates by likelihood-based methods. This EasyBadiRate can help you to investigate expansions and contractions of gene families with the gain-and-death (GD) model in BadiRate. 

Specifically, to calculate which lineages were significantly expanded or contracted, we used the free model in BadiRate to estimate the sizes of the ancestral gene families in all clades of the species tree. For the branches whose gene families did not experience family size changes, they were set to be the background branches, which have the same family turnover rate; thereafter, based on this model, we re-estimated the likelihood, and the model was regarded as the null hypothesis. For branches that experienced size changes, an alternative hypothesis for each branch was built by forcing the given branch to follow the same turnover rate with the background branches. A branch that experienced size changes was considered to be significant, if AIC (alternative hypothesis) - AIC (null hypothesis) > 2 (Akaike’s information criterion (AIC) was computed from the likelihood and numbers of parameters in each model). Hereafter, the significantly expanded and contracted gene families in all lineages were obtained.

Install
-------
The project relies on <a href="https://biopython.org/">Biopython</a>, then you just need to clone this tool.

        pip install biopython
        git clone https://github.com/SouthernCD/EasyBadiRate.git

Usage
-----

**input file**:

`tree_file`: a species tree file in newick

        (Osa:0.1344315,((Xvi:0.412211,(Tze:0.276903,Dal:0.206769):0.09238):0.033643,(Aof:0.319951,(Ash:0.242502,(Gel:0.221933,(Dca:0.087169,Peq:0.121297):0.038069):0.090386):0.150286):0.032129):0.268863);

`size_tsv_file`: a gene family size file in tsv file (Only one family can be counted at a time.)

        FAM_ID  Osa     Xvi     Tze     Dal     Aof     Ash     Gel     Dca     Peq
        OG1     1       1       1       1       1       1       4       1       1

**command**:

        python EasyBadiRate.py -l label_tree.tre OG1 tree_file size_tsv_file

**output**

stdout

        Tag     Gain    Loss    Likelihood
        OG1 ['13->9'] [] -2.47669280793592

label_tree.tre

        (Osa_1:0.1344315,((Xvi_2:0.412211,(Tze_3:0.276903,Dal_4:0.206769)5:0.09238)6:0.033643,(Aof_7:0.319951,(Ash_8:0.242502,(Gel_9:0.221933,(Dca_10:0.087169,Peq_11:0.121297)12:0.038069)13:0.090386)14:0.150286)15:0.032129)16:0.268863)17;