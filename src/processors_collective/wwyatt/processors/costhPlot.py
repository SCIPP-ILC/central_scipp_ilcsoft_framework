from plot import plot
c_plt=[
"cos_cm_2.root",
]
lab_plt=[
"cos_lab_2.root",
]
names = [
    ("ETheta", 2),
    ("PTheta", 4),
    ("Theta", 3),
    ]
plot(c_plt, names, "Theta Distribution Center of Mass", True)
plot(lab_plt, names, "Theta Distribution Lab", True)

