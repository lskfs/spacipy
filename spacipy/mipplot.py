
import numpy as np
from spacipy import Coords

def mip_spatial(adata, uns_key='3d_spatial', 
        color=None, highlight=None, keep_other=False, 
        groups=None, group_col=None, order=None, group_normalize=True, 
        ax=None, cmap=None, theme='dark', outfile=None, 
        vmin=None, vmax=None, percentile=False, 
        ps=1, coords_scale=0.5, 
        angles=None,
        body_loc=None, rotation=None, interval=100,
        figsize=None, density_radius=30,
        ):
    """
    Visualization of 3d values with maximum intensity projection

    Parameters
    ----------
    adata: anndata.AnnData object
        with '3d_spatial' key in adata.uns saved 3d spatial coordinates
    color: str
        column name in adata.obs or gene name in adata.var
    highlight: str, list
        if colored with category dtype, only highligh selected ones,
    keep_other: bool
        weather or not keep other non-highlight groups when using highlight
    group_col: str
        column name in adata.obs contains the batch information, if multiple groups 
        included, the program will try to rerrange the location of each group to make
        sure to plot all groups in a single panel.
    groups: str, list
        which group to use when input adata contains multiple batches
    order: list
        order of the groups in the plot
    group_normalize: bool
        weather or not normalize for each group when doing MIP
    ax
        ax to use
    cmap: str|list|dict of hex color code
        color used to plot
    theme: dark or light
    outfile: str
        output image name to save in local disk
    
    ps: float
        point size in plot
    coords_scale: float
        the coordinates will scaled by timing this number before mip transform
    
    """

    from pandas.api.types import is_numeric_dtype

    obj = adata.copy()

    if isinstance(color, str):
        color = [color]
    if isinstance(highlight, str):
        highlight = [highlight]
    
    def flatten(nested_list):
        for x in nested_list:
            if isinstance(x, list):
                yield from flatten(x)
            else:
                yield x

    if any(isinstance(i, list) for i in color):
        _flatten_color = list(flatten(color))
    else:
        _flatten_color = color
    _numeric_color = []
    _category_color = []
    for c in color:
        if c in obj.obs_keys():
            if is_numeric_dtype(obj.obs[c]):
                _numeric_color.append(c)
            else:
                _category_color.append(c)
        elif c in obj.var_names:
            _numeric_color.append(c)
    
    if _numeric_color and _category_color:
        raise 'only single panel supportted'
    if len(_category_color) > 1:
        raise 'only single category color type supportted'

    for c in _numeric_color:
        if c in obj.obs_keys():
            continue
        # detect feature as gene name
        arr = obj[:, c].X.toarray()
        obj.obs[c] = arr.reshape(arr.shape[0])
        #if exp_cutoff is not None:
        #    obj.obs.loc[obj.obs[c] < exp_cutoff, c] = 0
    
    td_uns = adata.uns.get(uns_key, {})

    cols = _numeric_color + _category_color
    if group_col is not None:
        cols += [group_col]
    data = obj.obs[cols].copy(deep=True)
    data[['x', 'y', 'z']] = td_uns

    if groups is not None:
        assert group_col is not None
        if isinstance(groups, str):
            groups = [groups]
        data = data[data[group_col].isin(groups)]
        if data[group_col].dtype == 'category':
            data[group_col] = data[group_col].cat.remove_unused_categories()

    coords3d = Coords(data)
    if angles:
        if group_col is not None:
            coords3d = coords3d.batch_rotation3d(
                    group_col,
                    angles=angles,
                    )
        else:
            coords3d = coords3d.rotation3d(
                    angles=angles,
                    )

    if group_col is not None:
        coords3d = coords3d.rearrange(
                group_col, order=order, 
                preset_mode=body_loc,
                rotation=rotation,
                interval=interval,
                )
    
    if _category_color:
        if highlight and keep_other:
            for _cc, _hl in zip(_category_color, highlight):
                if not _hl:
                    continue
                if isinstance(_hl, str):
                    _hl = [_hl]
                coords3d.loc[~coords3d[_cc].isin(_hl), _cc] = 'other'
        _category_values = np.unique(coords3d[_category_color].values).tolist()
        coords3d = coords3d.calculate_density(
                _category_color, 
                group_col=group_col, 
                radius=density_radius,
                )
        if highlight:
            if any(isinstance(i, list) for i in highlight):
                _flatten_highlight = list(flatten(highlight))
            else:
                _flatten_highlight = highlight
            if keep_other:
                _flatten_highlight += ['other']
            _del_category_values = [x for x in _category_values 
                    if x not in _flatten_highlight]
            cols = [x for x in coords3d.columns if x not in _del_category_values]
            coords3d = coords3d[cols]
            _flatten_color = _flatten_highlight
            color = highlight
        else:
            _flatten_color = _category_values
    
    coords3d = coords3d.mip_transform(
            _flatten_color, 
            group_col=group_col, 
            group_normalize=group_normalize,
            coords_scale=coords_scale,
            )
    
    if cmap is None:
        from spacipy import _colors_dark
        cmap = _colors_dark[:len(_flatten_color)]

    ax = coords3d.mip_scatter(
            color=_flatten_color, 
            cmap=cmap, 
            vmin=vmin,
            vmax=vmax,
            percentile=percentile, 
            outfile=outfile, 
            theme=theme,
            ax=ax,
            ps=ps, 
            figsize=figsize,
            )
    del coords3d
    return


