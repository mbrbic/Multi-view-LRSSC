## MultiViewLRSSC

MATLAB implementation of [Multi-view Low-rank Sparse Subspace Clustering Algorithm](https://arxiv.org/abs/1708.08732).

## Datasets

The datasets used in the paper can be found in the 'datasets' directory. [UCI digit](http://archive.ics.uci.edu/ml/datasets/Multiple+Features) and [Reuters](https://archive.ics.uci.edu/ml/datasets/Reuters+RCV1+RCV2+Multilingual,+Multiview+Text+Categorization+Test+collection#) datasets are from the UCI Machine Learning Repository. For Reuters dataset random sample of 100 documents per class is generated, resulting in a dataset consisting of 600 documents. [3-sources](http://mlg.ucd.ie/datasets/3sources.html) dataset is from the University College Dublin. Of 948 articles, we used 169 articles available in all three sources.  
[Prokaryotic phyla](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5137458/) dataset contains 551 prokaryotic species described with heterogeneous multi-view data including textual data and different genomic representations. Textual data consists of bag-of-words representation of documents describing prokaryotic species, while genomic representations include the proteome composition and gene repertoire representations. Proteome composition is encoded as the relative frequencies of amino acids and gene repertoire is encoded as the presence/absence indicators of gene families in a genome. We applied PCA on each of the three views separately and retain principal components explaining 90% of the variance. Each species in the dataset is labeled with the phylum it belongs to (cluster assignments).

## Citing

When using the code in your research work, please cite "Multi-view Low-rank Sparse Subspace Clustering" by Maria Brbic and Ivica Kopriva.

    @article{brbic2018,
    title={Multi-view Low-rank Sparse Subspace Clustering},
    author={Brbi\'c, Maria and Kopriva, Ivica},
    journal={Pattern Recognition},
    volume={73},
    pages={247--258},
    year={2018},
    doi = {https://doi.org/10.1016/j.patcog.2017.08.024}
    }

## Acknowledgements

This research project is supported by the Croatian Science Foundation grant IP-2016-06-5235 (Structured Decompositions of Empirical Data for 
Computationally-Assisted Diagnoses of Disease) and by the Croatian Science Foundation grant HRZZ-9623 (Descriptive Induction).
