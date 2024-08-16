
import numpy as np
import pandas as pd

from spacipy import get_cmap

class CoordsSeries(pd.Series):
    pass

class Coords(pd.DataFrame):
    """
    x   y   z
    """
    @property
    def _constructor(self):
        return Coords

    _constructor_sliced = CoordsSeries

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @property
    def center(self):
        x = int((self['x'].max() - self['x'].min()) / 2)
        px = x + self['x'].min()

        y = int((self['y'].max() - self['y'].min()) / 2)
        py = y + self['y'].min()
        return (px, py)

    def center3d(self):
        x = int((self['x'].max() - self['x'].min()) / 2)
        px = x + self['x'].min()

        y = int((self['y'].max() - self['y'].min()) / 2)
        py = y + self['y'].min()
        
        z = int((self['z'].max() - self['z'].min()) / 2)
        pz = z + self['z'].min()
        return (px, py, pz)

    def batch(self, group_col=None, func=None, **keywords):
        
        obj = self.copy(deep=True)
        
        obj = obj.groupby(group_col).apply(
                lambda x: getattr(x, func)(**keywords)
                )
        
        return obj
        
    def rotation3d(self, angles):
        
        assert len(angles) == 3
        from scipy.spatial.transform import Rotation as R

        obj = self.copy(deep=True)
        array = obj[['x', 'y', 'z']].values

        r = R.from_euler('XYZ', angles, degrees=True)
        rotated_array = r.apply(array)
        obj[['x', 'y', 'z']] = rotated_array

        return obj

    def batch_rotation3d(self, group_col, angles=None):

        if angles == 'left':
            angles = [90, 0, 90]
        elif angles == 'front':
            angles = [0, 90, 90]
        
        obj = self.copy(deep=True)
        
        obj = obj.batch(
                group_col=group_col,
                func='rotation3d',
                angles=angles
                )
        
        return obj
    
    def rotation(self, angle=90):
        import math

        obj = self.copy(deep=True)

        px, py = obj.center

        def srotate(angle, xs, ys, px, py):
            angle = math.radians(angle)
            x = (xs - px) * math.cos(angle) + \
                    (ys - py) * math.sin(angle) + px
            y = (ys - py) * math.cos(angle) - \
                    (xs - px) * math.sin(angle) + py

            return x, y
        xs, ys = srotate(
                angle,
                obj['x'].values,
                obj['y'].values,
                px,
                py
                )
        obj['x'] = xs
        obj['y'] = ys
        return obj

    def relocation(self, rotation=0, cx=0, cy=0):
        """relocated the coords by first rotating and then offset
        rotation: degree in clock wise
        offsetx: shift alone the x+, toward the right after rotation
        offsety: shift alone the y+, toward the upper after rotation
        """
        import math

        obj = self.copy(deep=True)

        if rotation:
            obj = obj.rotation(angle=rotation)

        px, py = obj.center
        if cx:
            offsetx = int(cx - px)
        else:
            offsetx = None
        if cy:
            offsety = int(cy - py)
        else:
            offsety = None

        if offsetx:
            obj['x'] = obj['x'] + offsetx
        if offsety:
            obj['y'] = obj['y'] + offsety

        return obj

    def location_stitch(self, group_col, order=None, interval=100):

        obj = self.copy(deep=True)
        
        groups = {}
        height = 0
        for name, group in obj.groupby(group_col):
            groups[name] = group

            h = group['y'].max() - group['y'].min() + 1
            if h > height:
                height = h
        cy = height // 2
        
        if not order:
            order = list(groups.keys())
        
        width = 0
        prev_x = 0
        location = {}
        for index, name in enumerate(order):
            group = groups[name]
            
            w = group['x'].max() - group['x'].min() + 1
            cx = w // 2 + prev_x
            #location[name] = {'cx': cx, 'cy': cy}
            location[name] = [0, cx, cy]
            
            prev_x += (w + interval)
        
        return location

    def rearrange(self, group_col, 
            order=None, 
            rotation=None,
            cx=None, cy=None,
            interval=100,
            preset_mode='auto'):

        obj = self.copy(deep=True)

        if preset_mode == 'auto':
            if rotation:
                obj = obj.batch(
                        group_col=group_col,
                        func='rotation',
                        angle=rotation
                        )
            preset_mode = obj.location_stitch(
                    group_col=group_col,
                    order=order, 
                    interval=interval,
                    )
        elif isinstance(preset_mode, dict):
            pass
        
        objs = []
        for name, group in obj.groupby(group_col):
            rotation, cx, cy = preset_mode[name]

            group = group.relocation(rotation=rotation, 
                    cx=cx, cy=cy)
            objs.append(group)
        obj = pd.concat(objs)
        
        return Coords(obj)

    def binning(self, bin_size=5):
        obj = self.copy(deep=True)
        obj['x'] = (obj['x'] / bin_size).astype(int)
        obj['y'] = (obj['y'] / bin_size).astype(int)
        return obj

    @property
    def img_shape(self):
        width = self['x'].max() - self['x'].min() + 1
        height = self['y'].max() - self['y'].min() + 1
        if isinstance(width, float):
            width = int(width) + 1
        if isinstance(height, float):
            height = int(height) + 1
        return (width, height)

    def alpha_shape(self, alpha=0.05):
        import alphashape
        from shapely.geometry import MultiPolygon

        obj = self.copy(deep=True)

        obj = obj[['x', 'y']].drop_duplicates()
        xy = obj[['x', 'y']].values
        
        hull = alphashape.alphashape(xy, alpha)
        if isinstance(hull, MultiPolygon):
            xy = [np.array(list(zip(*h.exterior.xy))) for h in hull.geoms]
        else:
            xy = [np.array(list(zip(*hull.exterior.xy)))]
        #patch = Polygon(xy, fill=False, color='w', linestyle='--', 
        #        linewidth=2, zorder=10)
        return xy

    def calculate_density(self, cluster_col, group_col=None, radius=20):
        from sklearn.neighbors import NearestNeighbors

        obj = self.copy(deep=True)

        if isinstance(cluster_col, str):
            cluster_col = [cluster_col]
        by = cluster_col[:]
        if group_col is not None:
            by.append(group_col)
        
        def func(group, radius=20):
            array = group[['x', 'y', 'z']].values
        
            nbrs = NearestNeighbors(n_neighbors=2, radius=5, algorithm='ball_tree').fit(array)
            #distances, indices = nbrs.radius_neighbors(array, radius = radius)
            indices = nbrs.radius_neighbors(array, radius=radius, return_distance=False)
        
            group['density'] = np.array([len(x) for x in indices])
            return group
        
        obj = obj.groupby(by=by).apply(func, radius=radius) #.reset_index()
        index = [x for x in obj.columns if x not in cluster_col + ['density']]
        obj = obj.pivot(index=index, columns=cluster_col, values='density')
        obj = obj.fillna(0).reset_index()

        return Coords(obj)

    def mip_transform(self, value, group_col=None, coords_scale=None, 
            group_normalize=True):
        from sklearn.preprocessing import QuantileTransformer

        obj = self.copy(deep=True)
        if not coords_scale:
            coords_scale = 1
        obj['x'] = (obj['x'] * coords_scale).astype(int)
        obj['y'] = (obj['y'] * coords_scale).astype(int)

        if isinstance(value, str):
            value = [value]
        if group_col and group_col in obj.columns:
            obj = obj[['x', 'y', group_col] + value]
            func = dict((v, 'max') for v in value)
            obj = obj.groupby(['x', 'y', group_col]).agg(func).dropna().reset_index()
        else:
            obj = obj[['x', 'y'] + value]
            func = dict((v, 'max') for v in value)
            obj = obj.groupby(['x', 'y']).agg(func).reset_index()
        
        def func(arr):
            if isinstance(arr, pd.Series):
                index = arr.index
                arr = arr.to_numpy()
            else:
                index = None
            qt = QuantileTransformer(output_distribution='uniform', random_state=0)
            scale_value = qt.fit_transform(arr.reshape(-1, 1))
            
            scale_value = scale_value.reshape(-1)
            scale_value = scale_value * 8
            scale_value = np.power((np.ones(len(scale_value)) * 2), scale_value)
            scale_value = scale_value - 1

            if len(arr.shape) > 1:
                scale_value = scale_value.reshape(arr.shape)
            
            if index is not None:
                scale_value = pd.Series(scale_value, index=index)
            return scale_value
        
        for v in value:
            if group_normalize and group_col:
                obj[v] = obj.groupby(group_col)[v].apply(func)
            else:
                obj[v] = func(obj[v].values)
        
        return Coords(obj)
    
    def mip_scatter(self,
            color=None, 
            cmap=None, ps=0.1, 
            vmin=None, vmax=None, percentile=False, 
            outfile=None, ax=None,
            theme='light',
            figsize=None,
            ):
        
        obj = self.copy(deep=True)

        if isinstance(color, str):
            color = [color]
        if isinstance(cmap, dict):
            cmap = [cmap[x] for x in color]

        from pathlib import Path
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from matplotlib.patches import Polygon, Rectangle

        if theme == 'dark':
            plt.style.use('dark_background')
        
        if not figsize:
            w, h = [i / 600 for i in obj.img_shape]
        else:
            w, h = [i / 600 for i in figsize]
        if ax is None:
            fig, ax = plt.subplots(figsize=(w, h), dpi=600)
        else:
            fig = plt.gcf()
            ax = ax

        import io
        import cv2
            
        obj['x'] = obj['x'].astype(int)
        obj['y'] = obj['y'].astype(int)
        x, y = obj.x.values, obj.y.values
        
        channels = []
        extent = []
        for name, c in zip(color, cmap):
            if name is None:
                continue
            if c == '#cccccc':
                ps = ps * 0.6
            layer_cmap = get_cmap(['#000000', c])
            v = obj[name].values
            zorder = np.argsort(v)

            if percentile:
                if (vmin is None) and (vmax is None):
                    vmin = np.percentile(v, 2)
                    vmax = np.percentile(v, 98)
                
            layer_fig, layer_ax = plt.subplots(figsize=(w, h), dpi=600)
            layer_ax.scatter(x[zorder], y[zorder], s=ps*10, c=v[zorder], lw=0, marker='.',  
                    edgecolors='none', cmap=layer_cmap, vmin=vmin, vmax=vmax)
            layer_ax.axis('off')
            layer_fig.subplots_adjust(left=0, bottom=0, right=1, top=1,
                    hspace=0, wspace=0)
            layer_fig.tight_layout(pad=0)
            if not extent:
                left, right = layer_ax.get_xlim()
                bottom, top = layer_ax.get_ylim()
                extent = [left, right, bottom, top]
            buf = io.BytesIO()
            layer_fig.savefig(buf, format='png', dpi=600) #, transparent=True)
            #plt.close(layer_fig)
            layer_fig.clear()
            buf.seek(0)
            img_arr = np.frombuffer(buf.getvalue(), dtype=np.uint8)
            buf.close()
            img = cv2.imdecode(img_arr, 1)
            img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
            channels.append(img)
            
        if theme == 'dark':
            image_stacks = np.flip(np.max(np.stack(channels, axis=3), axis=3), axis=0)
        elif theme == 'light':
            image_stacks = np.flip(np.min(np.stack(channels, axis=3), axis=3), axis=0)
        ax.imshow(image_stacks, norm=None, zorder=0, interpolation='none', 
                origin='lower', extent=extent) #, aspect='auto')
        ax.axis('off')
        
        if not outfile:
            plt.show()
        else:
            fig.subplots_adjust(left=0, bottom=0, right=1, top=1,
                    hspace=0, wspace=0)
            fig.tight_layout(pad=0)
            fig.savefig(outfile, dpi=600)
        return ax

CoordsSeries._constructor_expanddim = Coords


