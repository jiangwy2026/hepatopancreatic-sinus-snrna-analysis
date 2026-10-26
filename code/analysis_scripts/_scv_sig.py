import inspect
import scvelo as scv
print("filter_and_normalize:", inspect.signature(scv.pp.filter_and_normalize))
print("normalize_per_cell:", inspect.signature(scv.pp.normalize_per_cell))
print("filter_genes:", inspect.signature(scv.pp.filter_genes))
print("filter_genes_dispersion:", inspect.signature(scv.pp.filter_genes_dispersion))
