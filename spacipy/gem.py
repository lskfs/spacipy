#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd

class GemSeries(pd.Series):
    pass

class Gem(pd.DataFrame):
    """
    A Gem object is a pandas.DataFrame with adapted method for stereo-seq GEM file.
    """
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
        Alternate constructor to create a Gem from a file.
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

    def binning(self, bin_size=None, batch=None, inplace=False, outfile=None, 
            add_cols=False):
        """
        Group gem into square bins with fixed width and height.
        
        Parameters
        --------------------
        bin_size: int
            size of width and height of each bin
        batch: str
            column name for Batch information if contain
        add_cols: bool
            add bin coords values in new columns other than change the 
            original x, y coords
        outfile: str
            save grouped gem into file
        inplace: bool
            Whether the binning should modify the data in place or return
            a modified copy.
        """

        obj = self.copy(deep=True)

        obj.bin_size = bin_size
        
        if not add_cols:
            obj['x'] = (obj['x'] / bin_size).astype(int) * bin_size
            obj['y'] = (obj['y'] / bin_size).astype(int) * bin_size

            if 'label' in obj.columns:
                cols = ['geneID', 'x', 'y', 'label']
            else:
                cols = ['geneID', 'x', 'y']

            if batch and batch in obj.columns:
                cols.append(batch)

            obj = obj.groupby(cols).agg(
                    MIDCounts=('MIDCounts', 'sum'), 
                    )
            obj = obj.reset_index()
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
        """
        Convert gem into image using nCount or nFeature as pixel value.

        Parameters
        ----------
        on: str, default nCount
            how to represent each pixel of the image, 'nCount' or 'nFeature'
        bin_size: int
            group into bins before convert with width and height of bin_size
        outfile: str
            save the image
        reset_coords: bool, default True
            if True, the DNB/Bin with minimum (x,y) will be the origin 
            point in the output image; if False, the image will start at
            (0, 0), which counld cause huge image
        pad: int
            how many pixels will be padded around the image
        """

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
        """
        the shape of representative image.
        """
        width = self['x'].max() - self['x'].min() + 1
        height = self['y'].max() - self['y'].min() + 1

        if isinstance(width, float):
            width = int(width) + 1
        if isinstance(height, float):
            height = int(height) + 1
        return (height, width)
    
    def to_matrix(self, bin_size=None, outfile=None):
        """
        convert the tabular structure to sparse matrix.

        Parameters
        ----------
        bin_size: int
            group into bins before convert with width and height of bin_size
        outfile: str
            save the matrix
        """
        obj = self.copy(deep=True)

        adata = obj.to_anndata(bin_size=bin_size)
        
        mtx = pd.DataFrame.sparse.from_spmatrix(
                adata.X, 
                index=adata.obs.index.values, 
                columns=adata.var.index.values
                )
        if outfile:
            mtx.to_csv(outfile, sep=',', index=True, header=True)
        return mtx

    def lasso(self, x=None, y=None, width=None, height=None, 
            coords=None, mask=None, bin_size=None):

        obj = self.copy(deep=True)
        
        pass

    def roi(self, x=None, y=None, width=None, height=None, 
            inplace=False, outfile=None,):
        """
        Subset the gem object with specified coordinate, width and height

        Parameters
        ----------
        x, y: int
            bottom-left point coordinate
        width, height: int
            width and height from the bottom-left point
        outfile: str
            save the subsetted gem
        inplace: bool, default False
            Whether should modify the data in place or return
            a modified copy.
        """

        obj = self.copy(deep=True)

        obj = obj.loc[(
            (obj.x >= x) & 
            (obj.x < (x + width)) & 
            (obj.y >= y) & 
            (obj.y < (y + height))
            )]

        if outfile:
            obj.to_csv(outfile, sep='\t', index=False, compression='infer')

        if inplace:
            self._update_inplace(obj)
        else:
            return obj

    def mask(self, matrix, outfile=None, label_object=False, 
            return_offset=True, return_cropped_mask=False):
        """
        Subset the gem object with a numpy array mask

        Parameters
        ----------
        matrix: np.ndarray(int)
            mask matrix in int dtype, coordinate with overlapped 
            0 or nan will be dropped
        outfile: str
            save the subsetted gem
        label_object: bool, default False
            Whether or not to label the masked object, which will 
            generate a new 'label' column in the result object
        return_offset: bool, default True
            Whether or not to return the offset
        return_cropped_mask: bool, default False
            Whether or not to return the cropped mask matrix, which
            will only keep the kept region in result object
        """
        
        obj = self.copy(deep=True)
        
        assert obj.img_shape == matrix.shape, 'input with different shape'
        
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
        
        if return_offset and not return_cropped_mask:
            return obj, (offset_x, offset_y)
        elif return_cropped_mask and not return_offset:
            w, h = obj.img_shape
            matrix = matrix[offset_x: offset_x + w, offset_y: offset_y + h]
            return obj, matrix
        else:
            return obj
    
    def stat(self, genes=None, batch=None, bycell=False, outfile=None):
        """
        Compute some basic statistic for the object

        Parameters
        ----------
        genes: str, list
            Also consider certain genes for computing beside by
            default total statistic.
        batch: str
            column name for Batch information if contain
        bycell: bool, defaul False
            statistic at cell level other than at spot level
        outfile: str
            save the results
        """

        obj = self.copy(deep=True)
        if genes:
            obj_copy = obj.copy(deep=True)

        if batch:
            assert batch in obj.columns

        if bycell:
            assert 'label' in obj.columns
            cols = ['label']
            if batch:
                cols.append(batch)

            obj['coord'] = obj.x.astype(str) + '_' + obj.y.astype(str)
            cell_stat = obj.groupby(cols).agg(
                    #nCount_Spatial=('MIDCounts', 'sum'),
                    #nFeature_Spatial=('geneID', 'nunique'),
                    nCount=('MIDCounts', 'sum'),
                    nFeature=('geneID', 'nunique'),
                    nDNB=('coord', 'nunique')
                    ).reset_index()
            cell_stat = cell_stat[['label', 'nCount', 'nFeature', 'nDNB']]
            obj = obj.merge(cell_stat, how='left', on=['label'])
        else:
            cols = ['x', 'y']
            if batch:
                cols.append(batch)
        
            obj = obj.groupby(cols).agg(
                nCount=('MIDCounts', 'sum'), 
                nFeature=('geneID', 'nunique')
                )
        obj = obj.reset_index()

        if genes:
            if isinstance(genes, str):
                genes = [genes]
            obj_copy = obj_copy[obj_copy['geneID'].isin(genes)]
            obj_copy = obj_copy.groupby(['x', 'y', 'geneID']).agg(
                    MIDCounts=('MIDCounts', 'sum')
                    )
            obj_copy = obj_copy.pivot(
                    index=['x', 'y'], 
                    columns='geneID', 
                    values='MIDCounts').reset_index()
            obj = obj.merge(obj_copy, how='left', on=['x', 'y'])
        obj = obj.fillna(0)

        if outfile:
            obj.to_csv(outfile, sep='\t', index=False, compression='infer')

        return obj

    def to_anndata(self, bin_size=None, cell_loc=None, batch=None):
        """
        Convert the binning or cell labeling object to anndata

        Parameters
        ----------
        bin_size: int
            group into bins before convert with width and height of bin_size
        cell_loc: pandas.DataFrame
            dataframe with x and y coordinates for each cell label
        batch: str
            column name for Batch information if contain
        """

        from scipy import sparse
        import anndata

        obj = self.copy(deep=True)
        
        if bin_size:
            if 'label' in obj.columns:
                sys.stderr.write('guess you have labeled Gem, but treat '
                        'as binning data\n'
                        )
            obj = obj.binning(bin_size=bin_size, batch=batch)
            obj['label'] = obj['x'].astype(str) + '_' + obj['y'].astype(str)
        else:
            if 'label' not in obj.columns:
                sys.stderr.write('guess your data is in binning size\n')
            else:
                pass
        
        if batch:
            assert batch in obj.columns
            obj['label'] = obj['label'].astype(str) + '.' + obj[batch].astype(str)

            label2batch = obj[['label', batch]].drop_duplicates()
            label2batch = label2batch.set_index('label')

        label_obj = obj.groupby(['label', 'geneID']).agg(csum=('MIDCounts', 'sum'))
        label_obj = label_obj.reset_index()[['label', 'geneID', 'csum']]
        label_obj['label'] = label_obj['label'].astype(str)
        
        labels = list(set(label_obj.label.values))
        genes = list(set(label_obj.geneID.values))

        label2loc = dict(zip(labels, range(0, len(labels))))
        gene2loc = dict(zip(genes, range(0, len(genes))))

        rows = [label2loc[x] for x in label_obj.label]
        cols = [gene2loc[x] for x in label_obj.geneID]

        mtx = sparse.csr_matrix((label_obj.csum, (rows, cols)))
        
        obs = pd.DataFrame(index = labels)
        if batch:
            obs = obs.merge(label2batch, how='left', left_index=True, right_index=True)
        var = pd.DataFrame(index = genes)

        if cell_loc is not None and bin_size is None:
            if 'cell' in cell_loc.columns:
                cell_loc = cell_loc.rename(columns={'cell':'label'})
            if batch:
                assert batch in cell_loc.columns
                cell_loc['label'] = cell_loc['label'].astype(str) + '.' + cell_loc[batch].astype(str)
                cell_loc = cell_loc[['label', 'x', 'y']]
        elif cell_loc is None and bin_size is not None:
            cell_loc = obj[['label', 'x', 'y']].drop_duplicates()

        if cell_loc is not None:
            cell_loc['label'] = cell_loc['label'].astype(str)
            cell_loc = cell_loc.set_index('label')
            obsm = obs.merge(cell_loc, how='left', left_index=True, right_index=True)
            obsm = {'spatial': obsm.loc[:, ['x', 'y']].values}
        else:
            obsm = None

        adata = anndata.AnnData(
                X=mtx, 
                obs=obs, 
                var=var, 
                obsm=obsm, 
                )

        return adata

    def add_layer(self, image, mask=None):
        pass

    def _add_metadata(self, metadata, batch=None, bin_size=None, right_on=None,
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
                if right_on == 'index':
                    right_index = True
                    right_on = []
                else:
                    right_index = False
                    right_on = [right_on]

        if batch:
            assert batch in obj.columns, f'column {batch} not find in gem'
            left_on.append(batch)
            right_on.append(batch)

        obj = obj.merge(metadata, how='left', 
                left_on=left_on, 
                right_on=right_on,
                right_index = right_index
                )

        if inplace:
            self._update_inplace(obj)
        else:
            return obj

    def plot(self, 
            color=None, cmap=None, category=False, ncols=6,
            ps=0.1, alpha=1, robust=False, vmin=None, vmax=None,
            cbar=True, font_scale=1, dark_bg=True,
            dx=715, units='nm', scalebar_length=50,
            bin_size=None, bycell=False, batch=None, 
            metadata=None, meta_on='cell', image=None, 
            ax=None, outfile=None, return_fig=False, 
            ):
        """
        Plot the object in spatial coordinate

        Parameters
        ----------
        color: str, list
            Which information to plot, can be nCount, nFeature or any gene name
        cmap: str, list
            color map used for plot, can be matplotlib colormap or hex code
        category: bool, default False
            Whether or not the color is category
        ncols: int
            column number for multi-panel plot
        ps: float
            point size for scatter
        alpha: float, between 0-1
            alpha value for transparency
        robust: bool, default False
            color bar robust with 0.02 ~ 0.98
        vmin, vmax: float
            minimal and maximal value for color
        cbar: bool, default True
            Whether or not plot the color bar
        font_scale: float
            scale value for font size
        dark_bg: bool, default True
            Whether or not to plot in black background
        dx, units, scalebar_length: float
            deprecated scale bar parameters
        bin_size: int
            group into bins before convert with width and height of bin_size
        bycell: bool, default False
            Whether or not plot at cell level
        batch: str
            Column name for Batch information
        metadata: pandas.DataFrame
            deprecated
        meta_on: str
            deprecated
        image: numpy.ndarray
            numpy ndarray format image matrix, will be plot as background
        ax: matplotlib.Axis, list
            matplotlib.Axis object(s)
        outfile: str
            save the plot in file
        return_fig: bool, default False
            Whether or not to return the matplotlib figure
        """
        
        obj = self.copy(deep=True)

        if metadata is not None:
            obj = obj._add_metadata(metadata, batch=batch, bin_size=bin_size, 
                    right_on=meta_on)
        
        if isinstance(color, str):
            color = [color]
        genes = [x for x in color if x in obj.geneID.unique()]
        obj = obj.stat(genes=genes, batch=batch, bycell=bycell)

        assert all(x in obj.columns for x in color), f'{color} not find'
        
        obj = obj[['x', 'y'] + color]
        obj = obj.drop_duplicates()
        
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        #from matplotlib_scalebar.scalebar import ScaleBar
        
        if max(obj.img_shape) > 1000:
            shape_scale = 300
        else:
            shape_scale = 10

        if image is None:
            shape = [i / shape_scale for i in obj.img_shape[::-1]]
        else:
            shape = [i / shape_scale for i in image.shape[::-1]]
        
        if dark_bg:
            plt.style.use('dark_background')

        panels = len(color)
        if ax is None:
            if panels == 1:
                nrows = 1
                ncols = 1
            elif panels > 1 and panels <= ncols:
                nrows = 1
                ncols = panels
            elif panels > ncols:
                nrows = panels // ncols + 1
                ncols = ncols
            fig, ax = plt.subplots(nrows, ncols, figsize=shape, dpi=shape_scale)
            if nrows == 1 or ncols == 1:
                axes = [ax[index] for index in range(panels)]
            else:
                axes = [ax[index // ncols, index % ncols] for index in range(panels)]
        else:
            if panels != 1:
                assert isinstance(ax, list), 'must provide list of ax ' \
                                             'for multi-panel plot'
                assert len(ax) < panels, 'no enough panel for plot'
            fig = plt.gcf()
            axes = ax
        
        x = obj.x.values
        y = obj.y.values
        for c, axis in zip(color, axes):
            value = obj[c].values

            if robust:
                vmin = np.nanpercentile(value, 0.02)
                vmax = np.nanpercentile(value, 0.98)

            if vmin:
                vmin = vmin
            if vmax:
                vmax = vmax
        
            if category:
                N = len(set(value))
            else:
                N = 256
            cmap = get_cmap(cmap=cmap, N=N, category=category, theme='dark')
        
            if isinstance(value[0], str):
                value_map = dict((n, v + 1) for v, n in enumerate(set(value)))
                for n, v in value_map.items():
                    print(n, v)
                value = [value_map[n] for n in value]
            #cmap = cmap(np.linspace(0, 1, len(set(value))))

            if not ps:
                ps = (72. / fig.dpi) ** 2
                marker = ','
            else:
                ps = ps
                marker = '.'
            scatter = ax.scatter(x, y, 
                    s=ps,
                    c=value, 
                    cmap=cmap,
                    marker=marker,
                    linewidths=0,
                    edgecolors='none',
                    vmin=vmin, 
                    vmax=vmax,
                    alpha=alpha,
                    zorder=5,
                    )

            if image is not None:
                assert isinstance(image, np.ndarray)
            
                if bin_size:
                    import cv2
                    width = int(image.shape[0] / bin_size)
                    height = int(image.shape[1] / bin_size)
                    image = cv2.resize(image, (width, height), 
                            interpolation=cv2.INTER_NEAREST)
            
                imap = ax.imshow(image, cmap='gray', zorder=0, 
                        interpolation='none')


            """
            if bin_size:
                dx = dx * bin_size
            scalebar = ScaleBar(
                    dx=dx, 
                    units=units, 
                    fixed_value=scalebar_length, 
                    fixed_units='um',
                    border_pad=0.5,
                    color='white',
                    width_fraction=0.025,
                    length_fraction=0.3,
                    scale_loc='top',
                    location='lower right',
                    box_alpha=0,
                    label_loc='top',
                    sep=5*font_scale,
                    )
            ax.add_artist(scalebar)
            """

            ax.set_yticks([])
            ax.set_xticks([])
            ax.set_ylabel('')
            ax.set_xlabel('')
        
            _all = ["top", "bottom", "right", "left"]
            for spine in _all:
                ax.spines[spine].set_visible(False)
        
            if cbar:
                cax = inset_axes(ax, width='2%', height='30%', loc='upper right')
                if not category:
                    mappable = scatter
                    #cb = fig.colorbar(scatter, cax=cax, label=on,)
                else:
                    if isinstance(cmap, mpl.colors.ListedColormap):
                        norm = mpl.colors.NoNorm()
                    else:
                        bounds = np.arange(-1, len(set(value)), 1)
                        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
                    mappable = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
    
                if vmax and vmin:
                    extend = 'both'
                elif vmax:
                    extend = 'max'
                elif vmin:
                    extend = 'min'
                else:
                    extend = 'neither'

                cb = fig.colorbar(mappable, cax=cax, extend=extend) #label=on,)

                cb.outline.set_visible(False)
                cb.ax.tick_params(color='white', labelcolor='white', 
                        width=2*font_scale, length=10*font_scale)
                cb.ax.set_ylabel(color, fontweight='bold')
                cb.ax.yaxis.set_ticks_position('left')
                cb.ax.yaxis.set_label_position('left')
        
        plt.subplots_adjust(left=0, bottom=0, right=1, top=1, 
                hspace=0, wspace=0)
        mpl.rcParams.update({'font.size': font_scale * 18})

        if outfile:
            fig.savefig(outfile)
        if return_fig:
            return fig
        return ax

    @property
    def center(self):
        """
        Center of the object in both axis coordinate and data coordinate
            (axis_x, axis_y, data_x, data_y)
        """
        
        x = int((self.x.max() - self.x.min()) / 2)
        gx = x + self.x.min()
        y = int((self.y.max() - self.y.min()) / 2)
        gy = y + self.y.min()
        
        return (x, y, gx, gy)

    def relocation(self, x=None, y=None, canvas=None, width=None, height=None, 
            rotation=0, return_affine=False):
        """
        Re-locate the object coordinates with 
            - a target center coordinate
            - a target canvas (numpy ndarray matrix)
            - a target shape

        Parameters
        ----------
        x, y: int
            coordinates of the target center
        canvas: target canvas
        width, height: int
            pass
        rotation: int
            clock wise rotation degree
        return_affine: bool, default False
            Whether or not to return the affine matrix dict
        """

        import math

        obj = self.copy(deep=True)
        
        if canvas:
            width, height = canvas.shape
        if width:
            x = int(width / 2)
        if height:
            y = int(height / 2)
        
        orig_x, orig_y, orig_gx, orig_gy = obj.center

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

    #@profile
    def rearrange(self, batch, order=None, rotation=None, ncols=6, pad=True, 
            return_affine=False):
        """
        Auto rearrange the location for batch gem

        Parameters
        ----------
        bin_size: int
            Group into bins before convert with width and height of bin_size
        cell_loc: pandas.DataFrame
            Dataframe with x and y coordinates for each cell label
        batch: str
            Column name for Batch information if contain
        order: list
            Values from batch column for ordering control
        rotation: int
            Degree for clockwise rotation
        ncols: int
            Column number in the figure
        pad: bool, default True
            Whether append pad around each object
        return_affine: bool, default False
            Whether or not to return the affine matrix dict
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
        merged_affine = {}
        for index, (sample, degree) in enumerate(batches):
            ncol = index % ncols
            nrow = index // ncols

            x = int(max_width / 2) + ncol * max_width
            y = int(max_height / 2) + nrow * max_height
            
            gem = obj[obj[batch] == sample]
            gem, affine = gem.relocation(x=x, y=y, rotation=degree,
                    return_affine=True)
            merged_gem.append(gem)
            merged_affine[sample] = affine
        
        merged_gem = pd.concat(merged_gem)
        if return_affine:
            return merged_gem, merged_affine
        else:
            return merged_gem

GemSeries._constructor_expanddim = Gem


