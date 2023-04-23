
from spacipy import Gem

def gem_spatial(gem, color, 
        as_category=True, 
        groups=None, group_col=None, order=None, 
        label_col='label', uns_key='seg_spatial', img_key=None, 
        dx=0.715, scalebar_length=500,
        alpha_shape=False, alpha=0.05,
        outfile=None,
        ):
    
    obj = gem.copy(deep=True)

    obj.plot(
            on=color,
            outfile=outfile,
            )

    return


