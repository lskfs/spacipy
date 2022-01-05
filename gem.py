#!/usr/bin/env python3

import os
import sys
import numpy as np
import pandas as pd

#from pathlib import Path

class GemSeries(pd.Series):
    pass

class Gem(pd.DataFrame):
    @property
    def _constructor(self):
        return Gem

    _constructor_sliced = GemSeries
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.offset_x = None
        self.offset_y = None

    @classmethod
    def readin(cls, gem, **kwargs):
        """
        sample name should be SN or variant
        """

        df = pd.read_csv(gem, sep='\t', header=0, 
                compression='infer', comment='#', 
                **kwargs)

        if 'UMICount' in df.columns:
            count = 'UMICount'
        elif 'UMICounts' in df.columns:
            count = 'UMICounts'
        elif 'MIDCount' in df.columns:
            count = 'MIDCount'
        else:
            count = None
        if count:
            df = df.rename(columns={count: 'MIDCounts', })

        if 'cell' in df.columns:
            cell = 'cell'
        else:
            cell = None

        if cell:
            df = df.rename(columns={cell: 'label'})

        #if 'label' in df.columns:
        #    df['label'] = df['label'].astype(str)

        return cls(df, copy=True)

    def binning(self, bin_size=None, inplace=False, outfile=None, 
            add_cols=False):

        obj = self.copy(deep=True)

        obj.bin_size = bin_size
        
        if not add_cols:
            obj['x'] = (obj['x'] / bin_size).astype(int) * bin_size
            obj['y'] = (obj['y'] / bin_size).astype(int) * bin_size
        else:
            obj['xbin'] = (obj['x'] / bin_size).astype(int) * bin_size
            obj['ybin'] = (obj['y'] / bin_size).astype(int) * bin_size
        
        if outfile:
            obj.to_csv(outfile, sep='\t', index=False, compression='infer')

        if inplace:
            self._update_inplace(obj)
        else:
            return obj

    def to_img(self, on=None, bin_size=None, outfile=None, reset_coords=True, pad=0):

        if not on:
            on = 'nCount'

        obj = self.copy(deep=True)
        
        if reset_coords:
            obj['x'] = obj['x'] - obj['x'].min()
            obj['y'] = obj['y'] - obj['y'].min()

        if bin_size:
            if hasattr(obj, 'bin_size'):
                sys.stderr.write('WARN: your data is now in bin size {}, '
                        'DNB binning should always processed on DNB!'.format(
                                obj.bin_size))
            else:
                obj = obj.binning(bin_size=bin_size)

        if hasattr(obj, 'bin_size'):
            obj['x'] = (obj['x'] / obj.bin_size).astype(int)
            obj['y'] = (obj['y'] / obj.bin_size).astype(int)

        max_x = int(obj['x'].max()) + 1
        max_y = int(obj['y'].max()) + 1
        
        if on in ['nCount', 'nFeature']:
            obj = obj.stat()
        else:
            assert on in obj.columns

            obj = obj[['x', 'y', on]]
            obj = obj.drop_duplicates()
        
        image = np.zeros(shape=(max_y, max_x), dtype=int)
        image[obj['y'], obj['x']] = obj[on]

        if pad:
            image = np.pad(image, pad_width=pad)

        if outfile:
            import cv2
            cv2.imwrite(outfile, image)

        return image

    @property
    def img_shape(self):
        width = self['x'].max() - self['x'].min() + 1
        height = self['y'].max() - self['y'].min() + 1

        if isinstance(width, float):
            width = int(width) + 1
        if isinstance(height, float):
            height = int(height) + 1
        return (height, width)
    
    def to_matrix(self, bin_size=None, outfile=None):
        obj = self.copy(deep=True)

        adata = obj.to_anndata(bin_size=bin_size)
        
        mtx = pd.DataFrame.sparse.from_spmatrix(
                adata.X, 
                index=adata.obs, 
                columns=adata.var
                )
        if outfile:
            mtx.to_csv(outfile, sep=',', index=True, columns=True)
        return mtx

    def lasso(self, coords=None, mask=None, bin_size=None):

        obj = self.copy(deep=True)
        
        coords = pd.read_csv()
        
        pass

    def mask(self, matrix, outfile=None, label_object=False):
        
        obj = self.copy(deep=True)
        
        assert obj.img_shape == matrix.shape
        
        obj['x_mtx'] = obj['y'] - obj['y'].min()
        obj['y_mtx'] = obj['x'] - obj['x'].min()
        
        import scipy.sparse
        
        mtx= scipy.sparse.csc_matrix(matrix)
        mtx = mtx.tocoo()
        tmp = []
        for x, y, mask in zip(mtx.row, mtx.col, mtx.data):
            tmp.append([x, y, int(mask)])
        
        if label_object:
            columns = ['x_mtx', 'y_mtx', 'label']
            if 'label' in obj.columns:
                # relabeling gem
                sys.stdout.write('relabeling gem ... \n')
                obj = obj.drop(columns=['label'])
        else:
            columns = ['x_mtx', 'y_mtx', 'mask']
        triplet = pd.DataFrame(tmp, columns=columns)
        
        obj = obj.merge(triplet, how='inner', on=['x_mtx', 'y_mtx'])
        
        offset_x = obj.x.min() - obj.y_mtx.min()
        offset_y = obj.y.min() - obj.x_mtx.min()
        
        if (not label_object) or ('label' not in obj.columns):
            header = ['geneID', 'x', 'y', 'MIDCounts']
        else:
            header = ['geneID', 'x', 'y', 'MIDCounts', 'label']
        obj = obj[header]
        
        if outfile:
            obj.to_csv(outfile, sep='\t', index=False, compression='infer')
        
        return obj, offset_x, offset_y
    
    def stat(self, batch=None, outfile=None):
        obj = self.copy(deep=True)

        if batch:
            assert batch in obj.columns

        if 'label' in obj.columns:
            cols = ['x', 'y', 'label']
        else:
            cols = ['x', 'y']

        #if batch:
        #    cols.append(batch)

        obj = obj.groupby(cols).agg(
                nCount=('MIDCounts', 'sum'), 
                nFeature=('geneID', 'nunique')
                )
        
        obj = obj.reset_index()

        if outfile:
            obj.to_csv(outfile, sep='\t', index=False, compression='infer')

        return obj

    def to_anndata(self, bin_size=None, cell_loc=None):

        from scipy import sparse
        import anndata

        obj = self.copy(deep=True)
        
        if bin_size:
            if 'label' in obj.columns:
                sys.stderr.write('guess you have labeled Gem, but treat '
                        'as binning data\n'
                        )
            obj = obj.binning(bin_size=bin_size)
            obj['label'] = obj['x'].astype(str) + '_' + obj['y'].astype(str)
        else:
            if 'label' not in obj.columns:
                sys.stderr.write('guess your data is in binning size\n')
            else:
                pass

        label_obj = obj.groupby(['label', 'geneID']).agg(csum=('MIDCounts', 'sum'))
        label_obj = label_obj.reset_index()[['label', 'geneID', 'csum']]
        
        labels = set(label_obj.label.values)
        genes = set(label_obj.geneID.values)

        label2loc = dict(zip(labels, range(0, len(labels))))
        gene2loc = dict(zip(genes, range(0, len(genes))))

        rows = [label2loc[x] for x in label_obj.label]
        cols = [gene2loc[x] for x in label_obj.geneID]

        mtx = sparse.csr_matrix((label_obj.csum, (rows, cols)))
        
        obs = pd.DataFrame(index = labels)
        var = pd.DataFrame(index = genes)

        if cell_loc is None:
            obsm = obj[['label', 'x', 'y']].drop_duplicates()
            obsm = cell_loc.set_index('label')
        else:
            if 'cell' in cell_loc.columns:
                cell_loc = cell_loc.rename(columns={'cell':'label'})
            cell_loc['label'] = cell_loc['label'].astype(str)
            obsm = cell_loc.set_index('label')
        obsm = obs.merge(obsm, how='left', left_index=True, right_index=True)
        obsm = {"spatial": obsm.loc[:, ['x', "y"]].values}

        adata = anndata.AnnData(
                X=mtx, 
                obs=obs, 
                var=var, 
                obsm=obsm, 
                )

        return adata

    def add_layer(self, image, mask=None):
        pass

    def add_metadata(self, metadata, batch=None, bin_size=None, right_on=None,
            inplace=False):
        """ add metadata info for each DNB
        metadata: csv file or pandas DataFrame
        batch: colname in metadata which used as batch info
        left_on: colnames to join on in gem, default label
        right_on: colnames to join on in metadata, default label
        """
        import pandas as pd
        obj = self.copy(deep=True)

        if not isinstance(metadata, pd.DataFrame):
            metadata = pd.read_csv(metadata, sep='\t', header=0)
            #if 'seurat_clusters' in metadata.columns:
            #    metadata['seurat_clusters'] = metadata['seurat_clusters'] + 1
        if bin_size:
            obj = obj.binning(bin_size=bin_size, add_cols=True)
            obj['label'] = obj['xbin'].astype(str) + '_' + obj['ybin'].astype(str)

        left_on = ['label']
        if not right_on:
            if 'cell' in metadata.columns:
                right_on = ['cell']
            elif 'label' in metadata.columns:
                right_on = ['label']
        else:
            if isinstance(right_on, str):
                right_on = [right_on]

        if batch:
            assert batch in obj.columns, f'column {batch} not find in gem'
            left_on.append(batch)
            right_on.append(batch)

        obj = obj.merge(metadata, how='left', 
                left_on=left_on, 
                right_on=right_on
                )

        if inplace:
            self._update_inplace(obj)
        else:
            return obj

    def plot(self, 
            on=None, bin_size=None,
            metadata=None, batch=None,
            cmap=None, robust=False, 
            image=None, mask=None, border=None,
            offset_x=0, offset_y=0, 
            i_alpha=0.5, h_alpha=0.5, 
            outfile=None, ax=None,
            dx=715, units='nm', scalebar_length=50,
            ps=0.1, alpha=1, cbar=True,
            atx=None, aty=None, rotation=None,
            category=False, font_scale=1,
            ):
        
        obj = self.copy(deep=True)

        if atx or aty or rotation:
            obj, affine = obj.relocation(x=atx, y=aty, rotation=rotation, 
                    return_affine=True)
            if image is not None:
                image = scipy.ndimage.affine_transform(image.T, affine, order=0).T
            if mask is not None:
                mask = scipy.ndimage.affine_transform(mask.T, affine, order=0).T
            if border is not None:
                border = scipy.ndimage.affine_transform(border.T, affine, order=0).T
        
        if metadata is not None:
            obj = obj.add_metadata(metadata, batch=batch, bin_size=bin_size, 
                    right_on='cell')

        if offset_x or offset_y:
            obj['x'] = obj['x'] - offset_x
            obj['y'] = obj['y'] - offset_y

        if on in ['nCount', 'nFeature']:
            obj = obj.stat(batch=batch)

        assert on in obj.columns

        obj = obj[['x', 'y', on]]
        obj = obj.drop_duplicates()

        import matplotlib as mpl
        import matplotlib.pyplot as plt
        #from mpl_toolkits.axes_grid1 import axes_size, make_axes_locatable
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        from matplotlib_scalebar.scalebar import ScaleBar

        #plt.rc('savefig', dpi=100, bbox='tight', pad_inches=0)
        mpl.rcParams.update({"scalebar.color": 'white'})
        mpl.rcParams.update({"scalebar.width_fraction": 0.025})
        mpl.rcParams.update({"scalebar.length_fraction": 0.3})
        mpl.rcParams.update({"scalebar.scale_loc": 'top'})
        mpl.rcParams.update({"scalebar.location": 'lower right'})
        mpl.rcParams.update({"scalebar.box_alpha": 0})
        mpl.rcParams.update({'font.size': font_scale * 18})
        
        if image is None:
            shape = [i / 100 for i in obj.img_shape[::-1]]
        else:
            shape = [i / 100 for i in image.shape[::-1]]
        
        plt.style.use('dark_background')
        if ax is None:
            fig, ax = plt.subplots(figsize=shape, dpi=100)
        else:
            fig = plt.gcf()
            ax = ax
        
        x = obj.x.values
        y = obj.y.values
        value = obj[on].values

        if robust:
            vmin = np.nanpercentile(value, 0.02)
            vmax = np.nanpercentile(value, 0.98)
        else:
            vmin, vmax = None, None
        
        if isinstance(cmap, list):
            if category:
                N = len(set(value))
            else:
                N = len(cmap)
            cmap = mpl.colors.ListedColormap(cmap, N=N)
            if not category:
                cmap = cmap.colors
                cmap = mpl.colors.LinearSegmentedColormap.from_list('cmap', cmap)
        elif isinstance(cmap, str):
            if category:
                N = len(set(value))
            else:
                N = 256
            cmap = plt.get_cmap(cmap, N)
            if isinstance(cmap, mpl.colors.LinearSegmentedColormap):
                import matplotlib.colors as mc
                cmap = [mc.rgb2hex(cmap(i)) for i in range(cmap.N)]
                cmap = mpl.colors.ListedColormap(cmap, N=N)

        scatter = ax.scatter(x, y, 
                s=ps,
                c=value, 
                cmap=cmap,
                vmin=vmin, 
                vmax=vmax,
                alpha=alpha,
                zorder=1,
                )

        if image is not None:
            assert isinstance(image, np.ndarray)
            
            if mask is not None:
                assert isinstance(mask, np.ndarray)
                mask = np.isin(mask, [0])
                image[mask] = 0
            if border is not None:
                assert isinstance(border, np.ndarray)

                mask = np.isin(border, [1])
                image[mask] = 0
            
            imap = ax.imshow(image, cmap='gray', zorder=0, interpolation='none')

        scalebar = ScaleBar(
                dx=dx, 
                units=units, 
                fixed_value=scalebar_length, 
                fixed_units='um'
                )
        ax.add_artist(scalebar)

        ax.set_yticks([])
        ax.set_xticks([])
        ax.set_ylabel('')
        ax.set_xlabel('')
        
        _all = ["top", "bottom", "right", "left"]
        for spine in _all:
            ax.spines[spine].set_visible(False)
        
        if cbar:
            cax = inset_axes(ax, width='2%', height='30%', loc='upper left')
            if not category:
                cb = fig.colorbar(scatter, cax=cax, label=on,)
            else:
                if isinstance(cmap, mpl.colors.ListedColormap):
                    norm = mpl.colors.NoNorm()
                else:
                    bounds = np.arange(-1, len(set(value)), 1)
                    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
                cb = fig.colorbar(
                        mpl.cm.ScalarMappable(cmap=cmap, norm=norm),
                        cax=cax,
                        #label=on,
                        )
            #cb.ax.tick_params(color='white', labelcolor='white', labelsize=fontsize)
            #cb.ax.set_ylabel(on, fontsize=fontsize, fontweight='bold')
            cb.ax.tick_params(color='white', labelcolor='white')
            cb.ax.set_ylabel(on, fontweight='bold')
        
        plt.subplots_adjust(left=0, bottom=0, right=1, top=1, 
                hspace=0, wspace=0)
        if outfile:
            fig.savefig(outfile)
        return ax

    def center(self):
        
        x = int((self.x.max() - self.x.min()) / 2)
        gx = x + self.x.min()
        y = int((self.y.max() - self.y.min()) / 2)
        gy = y + self.y.min()
        
        return x, y, gx, gy

    def relocation(self, x=None, y=None, canvas=None, width=None, height=None, 
            rotation=0, return_affine=False):
        """ given a canvas (a matrix), relocate the gem in the middow
        x, y: target center
        canvas: target canvas
        width, height: target canvas shape
        rotation: clock wise rotation degree
        """

        import math

        obj = self.copy(deep=True)
        
        if canvas:
            width, height = canvas.shape
        if width:
            x = int(width / 2)
        if height:
            y = int(height / 2)
        
        orig_x, orig_y, orig_gx, orig_gy = obj.center()

        if x:
            offsetx = x - orig_x
        else:
            offsetx = 0

        if y:
            offsety = y - orig_y
        else:
            offsety = 0

        obj['x'] = obj['x'] - obj['x'].min() + offsetx
        obj['y'] = obj['y'] - obj['y'].min() + offsety
        
        def srotate(angle, xs, ys, px, py):
            angle = math.radians(angle)
            x = (xs - px) * math.cos(angle) + \
                    (ys - py) * math.sin(angle) + px
            y = (ys - py) * math.cos(angle) - \
                    (xs - px) * math.sin(angle) + py
            return x, y
        if rotation:
            #gx = orig_gx + offsetx
            #gy = orig_gy + offsety

            xs, ys = srotate(
                    rotation, 
                    obj.x.values,
                    obj.y.values,
                    x,
                    y
                    )
            
            obj['x'] = xs
            obj['y'] = ys
        
        if return_affine:
            affine = np.array([
                    [math.cos(rotation), math.sin(rotation), offsetx], 
                    [-math.sin(rotation), math.cos(rotation), offsety], 
                    [0, 0, 1]
                    ])
            return obj, affine
        else:
            return obj

    def rearrange(self, batch, order=None, rotation=None, ncols=6, pad=True):
        """auto rearrange the location for batch gem
        batch: which column should be detected as batch
        order: values from batch column for ordering control
        rotation: degree for clockwise rotation
        ncols: col number in the figure
        pad: whether append pad around each gem
        """

        obj = self.copy(deep=True)
        assert batch in obj.columns
        
        max_width, max_height = 0, 0
        batches = []
        for sample, group in obj.groupby(batch):
            batches.append(sample)
            w, h = group.img_shape
            if w > max_width:
                max_width = w
            if h > max_height:
                max_height = h
        if pad:
            max_width += int(max_width / 10)
            max_height += int(max_height / 10)

        if order:
            batches = order[:]
        if rotation:
            batches = list(zip(batches, rotation))
        else:
            batches = list(zip(batches, [0] * len(batches)))
        
        merged_gem = []
        for index, (sample, degree) in enumerate(batches):
            ncol = index % ncols
            nrow = index // ncols

            x = int(max_width / 2) + ncol * max_width
            y = int(max_height / 2) + nrow * max_height
            
            gem = obj[obj[batch] == sample]
            gem, affine = gem.relocation(x=x, y=y, rotation=degree,
                    return_affine=True)
            merged_gem.append(gem)
        
        merged_gem = pd.concat(merged_gem)
        return merged_gem

GemSeries._constructor_expanddim = Gem


