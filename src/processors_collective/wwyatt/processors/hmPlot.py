import plot

hitmiss_files = [
    "aa_positron_hitmiss.root",
    "aa_electron_hitmiss.root",
    "a_positron_hitmiss.root",
    "a_electron_hitmiss.root",
    "z_electron_hitmiss.root",
    "z_positron_hitmiss.root",
    ]
names = [
    "hh",
    "mm",
    "hm",
    ]

plot.scatterPlots(hitmiss_files, names)
