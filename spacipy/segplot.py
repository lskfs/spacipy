
import pandas as pd
from spacipy import Stoarr

def seg_spatial(adata, color=None, highlight=None,
        as_category=True, vmin=None, vmax=None,
        groups=None, group_col=None, order=None, dis=10,
        label_col='label', uns_key='seg_spatial', img_key=None, 
        scalebar=True, dx=0.715, scalebar_length=500,
        alpha_shape=False, alpha=0.05,
        cmap=None, outfile=None, cols=6,
        cbar=True, border=True, 
        bg='black', highlight_bg='#808080',
        ):

    """imshow plot in segmentation matrix

    format of adata.uns[uns_key]:
        {
          'batch1': 
            {
              'seg': 2d labeled arr,
              'img': multi-dims arr
            },
          'batch2':
            {
              'seg': 2d labeled arr,
              'img': multi-dims arr
            }
        }
    """

    obj = adata.copy()

    if color not in obj.obs_keys():
        # detect feature as gene name
        arr = obj[:, [color]].X.toarray()
        obj.obs[color] = arr.reshape(arr.shape[0])
        as_category = False
    assert color in obj.obs_keys()

    seg_uns = adata.uns.get(uns_key, {})
    seg_groups = list(seg_uns.keys())

    if group_col is None:
        group_col = 'Batch'
    if groups is not None:
        if isinstance(groups, str):
            groups = [groups]
        obj = obj[obj.obs[group_col].isin(groups)]
    else:
        if group_col not in obj.obs_keys():
            # assert as single batch adata object
            obj.obs[group_col] = seg_groups[0]
        groups = obj.obs[group_col].unique().tolist()
    
    data = obj.obs[[color, group_col, label_col]].copy(deep=True)
    
    base = 0
    seg_arr = []
    relabeled_data = []
    for index, name in enumerate(groups):
        data_group = data[data[group_col] == name]
        data_group[label_col] = data_group[label_col].astype(int).astype(str)
        data_group = data_group[[label_col, color]].set_index(label_col)
        data_group.index.name = None
        seg = seg_uns[name]['seg']
        seg = seg.view(Stoarr)
        seg, label_map = seg.relabel(base=base, return_map=True)
        base = seg.max()
        data_group = data_group.merge(label_map, how='left', 
                left_index=True, right_index=True)
        data_group = data_group.set_index('label')

        relabeled_data.append(data_group)
        seg_arr.append(seg)

    if len(groups) > 1:
        concat_seg = seg_arr[0].concat(seg_arr[1:], cols=cols, dis=dis)
    else:
        concat_seg = seg_arr[0]
    relabeled_data = pd.concat(relabeled_data)
    print(relabeled_data)

    ax = concat_seg.plot(
            annotation=relabeled_data, 
            highlight=highlight, 
            dx=dx, 
            scalebar=scalebar,
            scalebar_length=scalebar_length,
            alpha_shape=alpha_shape, 
            alpha=alpha, 
            cmap=cmap, 
            category=as_category, 
            outfile=outfile, 
            border=border,
            vmin=vmin,
            vmax=vmax,
            cbar=cbar,
            bg=bg,
            highlight_bg=highlight_bg,
            )

    return ax


