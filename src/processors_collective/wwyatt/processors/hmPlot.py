import plot
import sys
import Image

hitmiss_files = [
    "hitmiss_cm/aa.root",
    "hitmiss_cm/a.root",
    "hitmiss_cm/z.root",
    "hitmiss_lab/z.root",
    "hitmiss_lab/a.root",
    "hitmiss_lab/aa.root",
    ]
names = [
    ("hh", 3),
    ("hm", 2),
    ("mm", 4),
    ]
names = [
    ("mm", 4),
    ("hm", 2),
    ("hh", 3),
    ]

filenames = plot.scatterPlots(hitmiss_files, names)
f1 = filenames[:len(filenames)/2]
f2 = filenames[len(filenames)/2:]
fs = [f1,f2]
fs_n = {f1: "file1.png", f2: "file2.png"}
for f in fs:
    images = map(Image.open, f)
    widths, heights = zip(*(i.size for i in images))
    
    total_width = sum(widths)
    max_height = max(heights)
    
    new_im = Image.new('RGB', (total_width, max_height))
    
    x_offset = 0
    for im in images:
        new_im.paste(im, (x_offset,0))
        x_offset += im.size[0]

    new_im.save(fs_n[f])
