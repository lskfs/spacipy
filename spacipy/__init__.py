
colors = [
          '#d70000', '#8c3cff', '#028800', '#00acc7', '#98ff00', 
          '#ff7fd1', '#6c004f', '#ffa530', '#583b00', '#005759', 
          '#0000dd', '#00fdcf', '#a1756a', '#bcb7ff', '#95b578', 
          '#c004b9', '#645474', '#790000', '#0774d8', '#fef590', 
          '#004b00', '#8f7a00', '#ff7266', '#eeb9b9', '#5e7e66', 
          '#9be4ff', '#ec0077', '#a67bb9', '#5a00a4', '#04c600', 
          '#9e4b00', '#9c3b50', '#cbc400', '#718298', '#00af8a', 
          '#8388ff', '#5d373b', '#390000', '#fdc0ff', '#bee7c0', 
          '#db6d01', '#93b8b6', '#e452ff', '#2f5382', '#c46690', 
          '#556220', '#c59f72', '#048287', '#69e780', '#802790'
          ]

def get_cmap(cmap=None, N=None, category=False, theme='light'):

    import re
    import matplotlib as mpl

    if theme == 'light':
        facecolor = 'w'
    elif theme == 'dark':
        facecolor = 'k'

    is_hex_color = lambda x: re.search(r'^#(?:[0-9a-fA-F]{3}){1,2}$', x)
    if isinstance(cmap, str) and is_hex_color(cmap):
        cmap = [facecolor, cmap]

    if isinstance(cmap, list):
        if category:
            if N is None:
                N = len(cmap)
            cmap = mpl.colors.ListedColormap(cmap[:N], N=N)
        else:
            cmap = mpl.colors.LinearSegmentedColormap.from_list('cmap', 
                    cmap, N=256)
    elif isinstance(cmap, str):
        if not category:
            N = 256
        else:
            assert N is not None
        cmap = mpl.cm.get_cmap(cmap, N)
        if category and isinstance(cmap, mpl.colors.LinearSegmentedColormap):
            cmap = [mpl.colors.rgb2hex(cmap(i)) for i in range(cmap.N)]
            cmap = mpl.colors.ListedColormap(cmap, N=N)
    return cmap

_cmap = get_cmap(cmap=colors, N=len(colors), category=True)

from .gem import Gem
from .gemplot import gem_spatial
from .mask import Stoarr
from .segplot import seg_spatial



