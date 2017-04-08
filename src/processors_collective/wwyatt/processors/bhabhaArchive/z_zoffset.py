from plot import plot
files=[
    'z_zoffset.root',
    'a_zoffset.root',
    'aa_zoffset.root',
]
plots=[
'in',
'out',
]
plot(files, plots, "Outliers", folder='outliers',log=True, horizontal="outliers")
