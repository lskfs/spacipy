
import numpy as np
from spacipy import get_cmap

class Stoarr(np.ndarray):
    """A Stoarr object is a numpy.ndarray with adapted method for spatial matrix"""

    @classmethod
    def readin(cls, fname, img_mode=-1):
        """
        Alternate constructor to create a Stoarr from a file.
        """
        if not isinstance(fname, str):
            fname = str(fname)

        if fname.endswith(('.txt', '.txt.gz')):
            obj = np.loadtxt(fname, dtype=np.int32)
        elif fname.endswith('.npy'):
            obj = np.load(fname)
        elif fname.endswith(('.tif', '.png')):
            import cv2
            obj = cv2.imread(fname, img_mode)

        obj = obj.view(Stoarr)

        return obj

    def write(self, outfile, fmt='.txt'):
        """
        save the object into txt
        """

        if fmt == '.txt':
            np.savetxt(outfile, self, fmt='%d')
        elif fmt == '.npy':
            np.save(outfile, self)
        elif fmt == '.tif':
            import cv2
            cv2.imwrite(outfile, self)

        return

    def relabel(self, base=0, return_map=False):
        
        obj = self.copy()
        
        unique_labels, labels = np.unique(obj, return_inverse=True)
        if base:
            labels = np.add(labels, base, 
                    where=np.isin(labels, [0], invert=True)
                    )
        relabeled_obj = labels.reshape(obj.shape).view(Stoarr)
        
        if return_map:
            import pandas as pd

            orig_arr = obj.reshape(-1)
            relabeled_arr = relabeled_obj.reshape(-1)
            
            label_map = pd.DataFrame(
                    data=list(zip(orig_arr, relabeled_arr)), 
                    columns=['orig', 'label']
                    )
            label_map = label_map.drop_duplicates()
            label_map['orig'] = label_map['orig'].astype(int).astype(str)
            label_map['label'] = label_map['label'].astype(int).astype(str)
            label_map = label_map.set_index('orig')
            label_map.index.name = None
            return relabeled_obj, label_map
        else:
            return relabeled_obj

    def rotation(self, angle=None, contour_slice=True):
        """
        rotate the object clockwise
        
        Parameters
        ------------
        angle: int
            The rotation angle in degrees.
        contour_slice: bool, default True
            Whether or not slice the object with contour before rotate
        """

        from scipy.ndimage import rotate

        obj = self.copy()
        if contour_slice:
            obj = obj.contour_slice()

        obj = rotate(obj, angle=angle, reshape=True, order=0)
        """
        import numpy as np
        import scipy.ndimage
        w = len(obj.diagonal())

        affineR = np.array([
                [np.cos(angle), np.sin(angle), w],
                [-np.sin(angle), np.cos(angle), w],
                [0, 0, 1]
                ])
        affined_obj = scipy.ndimage.affine_transform(
                obj.T,
                affineR, 
                output_shape=(3*w, 3*w),
                order=0
                )
        obj = affined_obj.T
        """
        return obj

    def affine_transform(self, affineR, 
            reference=None, shape=None, order=0, cval=0,  
            outfile=None):
        """
        transform the objec with affine matrix

        Parameters
        ----------------
        affineR: array
            transform matrix
        reference: array
            reference matrix
        order: int
            pass
        outfile: str
            out file name
        """
        import scipy.ndimage
        
        obj = self.copy()

        if reference is not None:
            shape = reference.T.shape
        else:
            shape = None

        if shape is not None:
            shape = shape

        affined_matrix = scipy.ndimage.affine_transform(
                obj.T, 
                affineR, 
                output_shape=shape, 
                order=order, 
                cval=cval,
                )
        obj = affined_matrix.T.view(Stoarr)
        
        if outfile:
            obj.write(outfile, fmt='.tif')
        return obj

    def to_triplet(self, name='mask'):
        """
        Convert the object to triplet structure

        Parameters
        -----------------
        name: str
            column name of the value in triplet table
        """
        import pandas as pd
        import scipy.sparse

        obj = self.copy()
        mtx= scipy.sparse.csc_matrix(obj)
        mtx = mtx.tocoo()
        tmp = []
        for x, y, mask in zip(mtx.row, mtx.col, mtx.data):
            tmp.append([x, y, str(int(mask))])
        triplet = pd.DataFrame(tmp, columns=['x', 'y', name])
        return triplet

    def padding(self, pad=None):
        """
        Padding empty pixels around the object

        Parameters
        --------------
        pad: int
            number of pixels to pad
        """

        obj = self.copy()
        """
        if isinstance(pad, (list, tuple)) and len(pad) == 2:
            pad_width, pad_height = pad
        elif isinstance(pad, int):
            pad_width = pad
            pad_height = pad
        elif isinstance(pad, bool) and pad:
            w, h = obj.shape
            pad_width = int(w / 20)
            pad_height = int(h / 20)
        elif isinstance(pad, float) and pad <= 1:
            w, h = obj.shape
            pad_width = int(w * pad)
            pad_height = int(h * pad)

        obj = np.pad(obj, pad_width=(pad_width, pad_height))
        """
        obj = np.pad(obj, pad_width=pad, mode='constant', constant_values=0)
        return obj

    def contour_slice(self, transfer_to=None):
        """
        Slice the object with nearest contour box
        
        Parameters
        -------------
        transfer_to: np.ndarry or Stoarr object
            also slice other provided object with the same contour,
            will return a new sliced object
        """
        
        obj = self.copy()

        mask = np.isin(obj, [0], invert=True)
        index = np.nonzero(mask)

        min_x = np.min(index[0])
        max_x = np.max(index[0])
        min_y = np.min(index[1])
        max_y = np.max(index[1])

        obj = obj[min_x:max_x, min_y:max_y].view(Stoarr)

        if transfer_to is not None:
            if not isinstance(transfer_to, list):
                transfer_to = [transfer_to]
            out_mtx = []
            for mtx in transfer_to:
                mtx = mtx[min_x:max_x, min_y:max_y]
                out_mtx.append(mtx)
            if len(out_mtx) == 1:
                out_mtx = out_mtx[0]
            return obj, out_mtx
        else:
            return obj

    def _concat(self, other, axis=1, dis=10):

        obj = self.copy()
        if isinstance(other, np.ndarray):
            other = [other]
        
        contour_ = []
        max_ = 0
        chan_ = None
        for m in [obj] + other:
            if hasattr(m, 'hull'):
                m, h = m.contour_slice(transfer_to=[m.hull])
                m.hull = h
            else:
                m = m.contour_slice()
            if not chan_:
                chan_ = len(m.shape)
            if axis == 1 and m.shape[0] > max_:
                    max_ = m.shape[0]
            elif axis == 0 and m.shape[1] > max_:
                max_ = m.shape[1]
            contour_.append(m)
        
        padding_ = []
        hull_padding_ = []
        for m in contour_:
            if axis == 1:
                _pad = max_ - m.shape[0]
                pad_width = ((_pad // 2, _pad // 2 + _pad % 2), (dis, dis))
            elif axis == 0:
                _pad = max_ - m.shape[1]
                pad_width = ((dis, dis), (_pad // 2, _pad // 2 + _pad % 2))
            if chan_ == 3:
                pad_width += ((0, 0),)
            m = np.pad(m, pad_width=pad_width, 
                    mode='constant', constant_values=0)
            if hasattr(m, 'hull'):
                hull = np.pad(m.hull, pad_width=pad_width,
                        mode='constant', constant_values=0)
                hull_padding_.append(hull)
            padding_.append(m)
        
        obj = np.concatenate(padding_, axis=axis).view(Stoarr)
        if hull_padding_:
            obj.hull = np.concatenate(hull_padding_, axis=axis)
        
        return obj

    def concat(self, other, cols=None, dis=10):
        """
        Concat multiple objects into a single one
        
        Parameters
        -----------------
        other: np.ndarray/stoarr, list
            np.ndarray/stoarr object or list of objects
        cols: int
            column numbers 
        dis: int
            pixel number to pad between the objects
        """

        obj = self.copy()
        if isinstance(other, np.ndarray):
            other = [other]

        objs = [obj] + other
        if len(objs) <= cols:
            obj = objs[0]._concat(objs[1:], axis=1, dis=dis)
        else:
            rows_ = []
            for index in range(1, len(objs), cols):
                row_ = objs[index: index + cols]
                row_ = row_[0]._concat(row_[1:], axis=1, dis=dis)
                rows_.append(row_)
            
            obj = rows_[0]._concat(rows_[1:], axis=0, dis=dis)
        return obj

    def alpha_shape(self, alpha=0.05):
        """
        Compute the alpha shape (concave hull) around the object
        
        Will return a new object with attribute 'hull' which saving 
        the xy coordinates of the alpha shape

        Parameters
        ---------------
        alpha: float
            alpha value
        """
        import alphashape

        obj = self.copy()

        triplet = obj.to_triplet(name='label')

        xy = triplet[['x', 'y']].values
        
        hull = alphashape.alphashape(xy, alpha)
        xy = np.array(list(zip(*hull.exterior.xy)))
        #patch = Polygon(xy, fill=False, color='w', linestyle='--', 
        #        linewidth=2, zorder=10)
        obj.hull = xy

        return obj

    def find_boundaries(self, mode='inner'):
        """
        Compute the boundaries of the object.
        Return bool array where boundaries between labeled regions are True.

        Parameters
        ----------------
        mode: str, default inner
            mode in ['thick', 'inner', 'outer', 'subpixel'], 
            see skimage.segmentation.find_boundaries for detail
        """

        import skimage.segmentation as ss
        
        obj = self.copy()

        boundaries = ss.find_boundaries(
                obj, mode=mode).astype(bool)

        return boundaries

    def find_centroid(self):
        """
        Compute the centroid of the labelled object
        """

        import pandas as pd
        from skimage.measure import regionprops_table

        obj = self.copy()

        zero_mask = np.isin(obj, [0])
        obj[zero_mask] = obj.max() + 1

        unique, unique_indices, unique_inverse = np.unique(
                obj, 
                return_index=True, 
                return_inverse=True
                )
        relabed_obj = unique_inverse.reshape(obj.shape)
        min_mask = np.isin(relabed_obj, [0])
        relabed_obj[min_mask] = relabed_obj.max() + 1
        relabed_obj[zero_mask] = 0

        unique[-1] = 0
        label_map = pd.DataFrame(
                dict(orig_label=unique, 
                    label=unique_inverse[unique_indices]
                    )
            )
        props = regionprops_table(
                relabed_obj,
                properties=('label', 'centroid', 'area')
                )
        df = pd.DataFrame(props)
        df = df.merge(label_map, how='left', on='label')
        df = df[['orig_label', 'centroid-1', 'centroid-0', 'area']]
        df = df.rename(
                columns={
                    'orig_label':'label',
                    'centroid-0':'y',
                    'centroid-1':'x',
                    }
                )
        df['label'] = df['label'].astype(int).astype(str)
        return df
    
    @staticmethod
    def color_mapping(triplet, color=None, cmap='deep', category=True):
        
        import matplotlib as mpl
        hex2rgb = lambda x: tuple(int(x[i: i+2], 16) for i in (1, 3, 5))

        if category:
            values = triplet[color].unique()
            N = len(values)
            
            if isinstance(cmap, (str, list)):
                cmap = get_cmap(cmap, N=N, category=True)
                color_dict = dict((v, c) for v, c in zip(values, cmap.colors))
            elif isinstance(cmap, dict):
                color_dict = cmap
                label, colors = zip(*tuple(color_dict.items()))
                cmap = get_cmap(colors, category=True)
            
            triplet['hex'] = [color_dict[i] for i in triplet[color].values]
            triplet[['r', 'g', 'b']] = [hex2rgb(i) for i in triplet['hex'].values]
        else:
            cmap = get_cmap(cmap, category=False)
        
            norm = mpl.colors.Normalize(vmin=min(triplet[color].values), vmax=max(triplet[color].values))
            rgba = [cmap(i) for i in norm(triplet[color].values)]
            triplet[['r', 'g', 'b']] = [(round(r*255), round(g*255), round(b*255)) for r, g, b, a in rgba]
            #triplet[['r', 'g', 'b']] = [hex2rgb(i) for i in triplet['hex'].values]
        
        return triplet, cmap
    
    def plot(self, annotation=None, highlight=None, 
            bg='black', highlight_bg='#808080',
            scalebar=True, dx=0.715, scalebar_length=500, 
            alpha_shape=True, alpha=0.05, border=True, border_color='#000000', 
            cmap=None, category=True, cbar=False, 
            ax=None, outfile=None, vmin=None, vmax=None):
        """
        Imaging the object
        
        Parameters
        ---------------------
        annotation: pandas.DataFrame
            one column pandas DataFrame containing annotation,
            with object labels as index
        highlight: list
            which value in the annotation table to highlight
        highlight_bg: str
            hex code for the non-highlight values
        bg: str
            background color
        scalebar: bool, default True
            whether or not to plot the scale bar
        dx: float
            pixel scale for scale bar calculation
        scalebar_length: float
            length of the scale bar
        alpha_shape: bool
            pass
        alpha: float
            pass
        border: bool, default True
            whether or not to plot the border for cell segment
        border_color: str, default $000000
            cell border color
        cmap: str, list
            color map used for plot, can be matplotlib colormap or hex code
        category: bool, default True
            whether the value is category or continous type
        cbar: bool, default False
            whether or not to plot the color bar
        ax: matplotlib.Axis
            matplotlib.Axis
        outfile: str
            save the image
        vmin, vmax: int
            minimal and maximal value range for color        
        """

        if isinstance(cmap, dict):
            color_dict = cmap

        import matplotlib as mpl
        import matplotlib.pyplot as plt
        #from matplotlib_scalebar.scalebar import ScaleBar
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        from matplotlib.patches import Polygon, Rectangle
        
        obj = self.copy()

        if border:
            _border_mask = obj.find_boundaries(mode='inner')
        else:
            _border_mask = None

        if annotation is not None:
            import pandas as pd
            from pandas.api.types import is_numeric_dtype
            from pandas.api.types import is_string_dtype

            triplet = obj.to_triplet(name='label')
            triplet = triplet.set_index('label')

            if isinstance(annotation, str):
                annotation = pd.read_csv(annotation, sep='\t', header=0, index_col=0)

            triplet = triplet.merge(annotation, how='inner', left_index=True, right_index=True)
            col = list(annotation.columns)[0]
            if is_string_dtype(triplet[col]) or category:
                category = True
            else:
                category = False
            triplet, cmap = self.color_mapping(triplet, col, cmap=cmap, category=category)
            
            w, h = obj.shape[:2]
            if bg in ['black', 'k']:
                obj = np.full((w, h, 3), 0, dtype=int)
            elif bg in ['white', 'w']:
                obj = np.full((w, h, 3), 255, dtype=int)
            obj[triplet['x'], triplet['y']] = triplet[['r', 'g', 'b']]
        
        if _border_mask is not None:
            obj[_border_mask] = (0, 0, 0)

        if max(obj.shape) > 1000:
            shape_scale = 300
        else:
            shape_scale = 100

        shape = [i / shape_scale for i in obj.shape[:2][::-1]]
        plt.style.use('dark_background')
        if ax is None:
            fig, ax = plt.subplots(figsize=shape, dpi=shape_scale)
        else:
            fig = plt.gcf()
            ax = ax
        
        im = ax.imshow(obj, zorder=0, interpolation='none', 
                vmin=vmin, vmax=vmax, norm=None)
        
        from matplotlib.lines import Line2D
        def _scatter_legend(cmap, ax=None, loc='upper center', ncol=4,):
            legend_elements = []
            for label, color in cmap.items():
                element = Line2D(
                        [0], [0], marker='o', 
                        markersize=5, 
                        markeredgecolor='none',
                        linewidth=0,
                        color=color, 
                        label=label
                        )
                legend_elements.append(element)

            ax.legend(handles=legend_elements, loc=loc, ncol=ncol, 
                    markerscale=0.5, fontsize='small')
            return ax

        if cbar:
            #cax = inset_axes(ax, width='2%', height='30%', loc='upper right')
            if category:
                _scatter_legend(color_dict, ax=ax, loc='upper right', ncol=1,)
            #else:

            """
            if not category:
                if vmin is None:
                    vmin = 0
                if vmax is None:
                    vmax = obj.max()
                norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
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

            try:
                cb = fig.colorbar(mappable, cax=cax, ticks=list(range(len(label))))
                cb.ax.set_yticklabels(label)
            except:
                cb = fig.colorbar(mappable, cax=cax)

            cb.outline.set_visible(False)
            cb.ax.tick_params(color='white', labelcolor='white', )
            #cb.ax.set_ylabel(on, fontweight='bold')
            cb.ax.yaxis.set_ticks_position('left')
            cb.ax.yaxis.set_label_position('left')
            """
        
        if alpha_shape:
            xyhull = obj.alpha_shape(alpha=alpha)
            patch = Polygon(xyhull, fill=False, color='w', linestyle='--', 
                    linewidth=2, zorder=10)
            ax.add_patch(patch)
        
        if scalebar:
            # add scale bar
            x_start = int(obj.shape[1] / 50)
            x_length = scalebar_length / dx
            y_start = int(obj.shape[0] / 20)
            y_length = x_length / 50
            if y_length < 1:
                y_length = 1
            rect = Rectangle(
                    (x_start, y_start),
                    width=x_length,
                    height=y_length,
                    ec=None, fc='w'
                    )
            ax.add_patch(rect)
            ax.text(x_start + x_length + 5, y_start + y_length / 2, f'{scalebar_length} um', 
                    size=5, color='w', ha='left', va='center')
        #ax.invert_yaxis()

        ax.set_yticks([])
        ax.set_xticks([])
        ax.set_ylabel('')
        ax.set_xlabel('')
        
        _all = ["top", "bottom", "right", "left"]
        for spine in _all:
            ax.spines[spine].set_visible(False)
        
        plt.subplots_adjust(left=0, bottom=0, right=1, top=1, 
                hspace=0, wspace=0)

        if outfile:
            fig.savefig(outfile, dpi=600)
        return ax


