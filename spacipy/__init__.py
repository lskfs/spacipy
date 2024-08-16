
_colors_dark = [
        '#d60000', '#e2afaf', '#018700', '#a17569', '#e6a500', '#004b00',
        '#6b004f', '#573b00', '#005659', '#5e7b87', '#0000dd', '#00acc6',
        '#bcb6ff', '#bf03b8', '#645472', '#790000', '#0774d8', '#729a7c',
        '#8287ff', '#ff7ed1', '#8e7b01', '#9e4b00', '#8eba00', '#a57bb8',
        '#5901a3', '#8c3bff', '#a03a52', '#a1c8c8', '#f2007b', '#ff7752',
        '#bac389', '#15e18c', '#60383b', '#546744', '#380000', '#e252ff',
        ]

_colors_light = [
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

def get_cmap(cmap=None, N=None, category=False): #, theme='light'):

    import re
    import matplotlib as mpl

    #if theme == 'light':
    #    facecolor = 'w'
    #elif theme == 'dark':
    #    facecolor = 'k'

    is_hex_color = lambda x: re.search(r'^#(?:[0-9a-fA-F]{3}){1,2}$', x)
    #if isinstance(cmap, str) and is_hex_color(cmap):
    #    cmap = [facecolor, cmap]

    if isinstance(cmap, (list, tuple)):
        if category:
            if N is None:
                N = len(cmap)
            cmap = mpl.colors.ListedColormap(cmap[:N], N=N)
        else:
            if N is None:
                N = 256
            cmap = mpl.colors.LinearSegmentedColormap.from_list('cmap', 
                    cmap, N=N)
    elif isinstance(cmap, str):
        if not category:
            if N is None:
                N = 256
        else:
            assert N is not None
        cmap = mpl.cm.get_cmap(cmap, N)
        if category and isinstance(cmap, mpl.colors.LinearSegmentedColormap):
            cmap = [mpl.colors.rgb2hex(cmap(i)) for i in range(cmap.N)]
            cmap = mpl.colors.ListedColormap(cmap, N=N)
    return cmap

from .gem import Gem
from .stoarr import Stoarr
from .segplot import seg_spatial
from .coords3d import Coords
from .mipplot import mip_spatial

