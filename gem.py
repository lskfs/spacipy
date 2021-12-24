#!/usr/bin/env python3

import os
import sys
import numpy as np
import pandas as pd

#from pathlib import Path

class Gem(pd.DataFrame):
    @property
    def _constructor(self):
        return Gem
    
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

        if 'label' in df.columns:
            df['label'] = df['label'].astype(str)

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

        max_x = obj['x'].max() + 1
        max_y = obj['y'].max() + 1
        
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
    def shape(self):
        img = self.to_img()
        return img.shape
    
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
        
        assert obj.to_img().shape == matrix.shape
        
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
    
    def stat(self, outfile=None):
        obj = self.copy(deep=True)

        if 'label' in obj.columns:
            cols = ['x', 'y', 'label']
        else:
            cols = ['x', 'y']
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

    def plot(self, on=None, metadata=None, 
            bin_size=None, cmap=None, 
            robust=False, image=None, 
            mask=None, border=None,
            offset_x=0, offset_y=0, 
            i_alpha=0.5, h_alpha=0.5, 
            outfile=None, ax=None,
            dx=715, units='nm',
            ps=0.1, alpha=1, cbar=True,
            atx=None, aty=None, rotation=None,
            category=False,
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
            import pandas as pd

            if not isinstance(metadata, pd.DataFrame):
                metadata = pd.read_csv(metadata, sep='\t', header=0)
            #if 'seurat_clusters' in metadata.columns:
            #    metadata['seurat_clusters'] = metadata['seurat_clusters'] + 1
            if bin_size is None:
                assert ('cell' in metadata.columns) or ('label' in metadata.columns)
                if 'cell' in metadata.columns:
                    metadata = metadata.rename(columns={'cell':'label'})
                metadata['label'] = metadata['label'].astype(str)
                obj = obj.merge(metadata, how='left', on='label')
            else:
                assert ('xbin' in metadata.columns) and ('ybin' in metadata.columns)
                obj = obj.binning(bin_size=bin_size, add_cols=True)
                obj = obj.merge(metadata, how='left', on=['xbin', 'ybin'])

        if offset_x or offset_y:
            obj['x'] = obj['x'] - offset_x
            obj['y'] = obj['y'] - offset_y

        if on in ['nCount', 'nFeature']:
            obj = obj.stat()

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
        
        if image is None:
            if offset_x or offset_y:
                heatmap = obj.to_img(reset_coords=False)
            else:
                heatmap = obj.to_img()
            shape = [i / 100 for i in heatmap.shape[::-1]]
        else:
            shape = [i / 100 for i in image.shape[::-1]]
        
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
                cmap = cmap[:max(value) + 1]
            cmap = mpl.colors.ListedColormap(cmap)

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
                fixed_value=50, 
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
                norm = mpl.colors.NoNorm()
                cb = fig.colorbar(
                        mpl.cm.ScalarMappable(cmap=cmap, norm=norm),
                        cax=cax,
                        label=on,
                        )
            cb.ax.tick_params(color='white', labelcolor='white')

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

        obj['x'] = obj['x'] + offsetx
        obj['y'] = obj['y'] + offsety
        
        def srotate(angle, xs, ys, px, py):
            import math
            xs = (xs - px) * math.cos(angle) + \
                    (ys - py) * math.sin(angle) + px
            ys = (ys - py) * math.cos(angle) - \
                    (xs - px) * math.sin(angle) + py
            return xs, ys
        if rotation:
            gx = orig_gx + offsetx
            gy = orig_gy + offsety

            xs, ys = srotate(
                    rotation, 
                    obj.x.values,
                    obj.y.values,
                    gx,
                    gy
                    )
            
            obj['x'] = xs
            obj['y'] = ys
        
        if return_affine:
            affine = np.ndarray([
                    [math.cos(rotation), math.sin(rotation), offsetx], 
                    [-math.sin(rotation), math.cos(rotation), offsety], 
                    [0, 0, 1]
                    ])
            return obj, affine
        else:
            return obj

