"""usage: {} --cell <cell_mask> --tissue <tissue_mask> --blur <blur_mask> --image <orig_img.tif>
"""

import os
import sys
import getopt
import copy

from collections import Counter

import numpy as np
import pandas as pd
import cv2

from pathlib import Path

HOME = Path(__file__).parent

def getfootprint(struc, a, b=None):
    from skimage.morphology import (
            square, 
            rectangle, 
            diamond, 
            disk, 
            octagon, 
            star)

    struc_lib = {
            'square': square, 
            'rectangle': rectangle, 
            'diamond': diamond, 
            'disk': disk, 
            'octagon': octagon, 
            'star': star
            }

    morph = struc_lib[struc]

    if struc in ['rectangle', 'octagon']:
        if b is None:
            sys.stderr.write('two args required\n')
            sys.exit()
        return morph(a, b)
    else:
        if b is not None:
            sys.stderr.write('only one arg required\n')
            sys.exit()
        return morph(a)

class Mask:
    def __init__(self, matrix):
        if isinstance(matrix, str):
            if matrix.endswith('.txt'):
                matrix = np.loadtxt(matrix, dtype=int)
            elif matrix.endswith(('.tif', '.png')):
                matrix = cv2.imread(matrix, cv2.IMREAD_UNCHANGED)
        self.matrix = matrix.astype(int)
        
    def to_triplet(self, name='mask'):
        import scipy.sparse

        mtx= scipy.sparse.csc_matrix(self.matrix)
        mtx = mtx.tocoo()
        tmp = []
        for x, y, mask in zip(mtx.row, mtx.col, mtx.data):
            tmp.append([x, y, int(mask)])
        triplet = pd.DataFrame(tmp, columns=['x', 'y', name])
        return triplet

    def binning(self, bin_size):

        sys.stdout.write('binning ... ')
        sys.stdout.flush()

        triplet = self.to_triplet()
    
        triplet['xbin'] = (triplet.x / bin_size).astype(int) * bin_size
        triplet['ybin'] = (triplet.y / bin_size).astype(int) * bin_size
        triplet['bin'] = triplet.xbin.astype(str) + '_' + triplet.ybin.astype(str)

        index = [(-i, x) for i, x in enumerate(triplet['bin'].unique())]
        index = pd.DataFrame(index, columns=['N', 'bin'])

        triplet = triplet.merge(index, how='left', on='bin')
    
        matrix = np.zeros(shape=self.matrix.shape, dtype=int)
        matrix[triplet['x'], triplet['y']] = triplet['N']

        sys.stdout.write('done\n')
        return Mask(matrix)

    def to_binary(self):
        obj = copy.deepcopy(self)
        mask = np.isin(obj.matrix, [0], invert=True)
        obj.matrix[mask] = 1
        return obj
    
    def subtract(self, other):

        sys.stdout.write('subtracting ... ')
        sys.stdout.flush()

        obj = copy.deepcopy(self)

        obj = obj.to_binary()
        other = other.to_binary()

        obj.matrix = obj.matrix - other.matrix

        sys.stdout.write('done\n')
        return obj
    
    def intersection(self, other, label_area_cutoff=0.3):
        """intersection of label mask and binary mask
        * mask: binary matrix
        * label_area_cutoff: labels with greater area will be dropped
        """
        
        sys.stdout.write('intersection ... ')
        sys.stdout.flush()

        obj = copy.deepcopy(self)

        if isinstance(other, Mask):
            other = other.to_binary()

        values = np.unique(obj.matrix)
        if len(values) == 2:
            mask = cv2.bitwise_and(obj.matrix, other.matrix)
            mask = np.invert(mask.astype(bool))
        else:
            binary = self.to_binary()
            
            mask = cv2.bitwise_and(binary.matrix, other.matrix)
            mask = np.invert(mask.astype(bool))

            orig_counter = Counter(obj.matrix.flatten())

            filter_part = obj.matrix[mask]
            filter_counter = Counter(filter_part.flatten())

            filter_labels = []
            for label, pixels in filter_counter.items():
                if label == 0:
                    continue
                ratio = pixels / orig_counter[label]
                if ratio < label_area_cutoff:
                    continue
                filter_labels.append(label)

            filter_labels = list(set(filter_labels))
            mask = np.isin(obj.matrix, filter_labels)

        obj.matrix[mask] = 0

        sys.stdout.write('{} labels removed\n'.format(len(filter_labels)))
        return obj

    def relabel(self, label_map=None, coord_map=None):
        if label_map is None and coord_map is None:
            unique_labels, labels = np.unique(self.matrix, return_inverse=True)
            matrix = labels.reshape(self.matrix.shape)

            #obj = Mask(matrix)
            #obj.unique_labels = unique_labels
            #obj.labels = labels
            return Mask(matrix)
        else:
            triplet = self.to_triplet()
            
            if label_map is not None:
                triplet = triplet.merge(label_map, how='left', 
                        left_on='mask', right_index=True)
            elif coord_map is not None:
                min_x, min_y = coord_map.x.min(), coord_map.y.min()
                
                if min_x != 0 or min_y != 0:
                    if min_x != 0:
                        min_x = min_x
                    if min_y != 0:
                        min_y = min_y
                else:
                    min_x = 0
                    min_y = 0
                
                coord_map['x'] = coord_map['x'] - min_x
                coord_map['y'] = coord_map['y'] - min_y
                
                triplet = triplet.merge(
                        coord_map, 
                        how='left', 
                        left_on=['x', 'y'], 
                        right_on=['y', 'x']
                        )
            triplet = triplet.fillna(0)
    
            matrix = np.zeros(shape=self.matrix.shape, dtype=int)
            if label_map is not None:
                matrix[triplet['x'], triplet['y']] = triplet['mask_y']
            elif coord_map is not None:
                matrix[triplet['x_x'], triplet['y_x']] = triplet['mask_y']
            return Mask(matrix)
    
    def retrieve(self):
        if not self.unique_labels and not self.labels:
            return

        matrix = self.unique_labels[self.labels]
        matrix = matrix.reshape(self.shape)
        obj = Mask(matrix)

        return obj

    def minimum_filter(self, footprint='octagon', ksize=(4, 4), iterations=2):

        sys.stdout.write('minimum filter ... ')
        sys.stdout.flush()

        obj = copy.deepcopy(self)
        obj.matrix = obj.matrix.astype(np.uint8)
        #obj.matrix = cv2.applyColorMap(
        #        obj.matrix,
        #        cv2.COLORMAP_JET
        #        )
        
        try:
            n, m = ksize
        except:
            n = ksize
            m = None
        footprint = getfootprint(footprint, n, m)
        obj.matrix = cv2.erode(
                obj.matrix, 
                kernel=footprint, 
                iterations=iterations
                )
        #cv2.imwrite('blur.png', obj.matrix)

        sys.stdout.write('done\n')
        return obj
    
    def filter_by_matrix(self, on=None, min_value=None, max_value=None, 
            draw=False, prefix=None):
        """label mask method
        * on: filter by minimum value of the input matrix
        """
        
        sys.stdout.write('filter by matrix ... ')
        sys.stdout.flush()

        obj = copy.deepcopy(self)

        triplet = obj.to_triplet()
        ref = on.to_triplet()
        triplet = triplet.merge(ref, how='left', on=('x', 'y'))
        triplet = triplet.fillna(0)

        medians = triplet.groupby('mask_x')['mask_y'].median()
        medians = medians.to_frame()

        if draw:
            fig = self.relabel(medians)
            cv2.imwrite(f'{prefix}.median.png', fig.matrix)

        if min_value:
            filter_labels = medians[medians['mask_y'] < min_value].index.values
        if max_value:
            filter_labels = medians[medians['mask_y'] > max_value].index.values
        
        mask = np.isin(obj.matrix, filter_labels)
        obj.matrix[mask] = 0

        sys.stdout.write('{} labels removed\n'.format(len(filter_labels)))
        return obj

    def filter_by_diameter(self, min_size=1, max_size=None):
        """label mask method
        * min_size: max circo radius
        """

        sys.stdout.write('filter by diameter ... ')
        sys.stdout.flush()

        from skimage.measure import regionprops

        obj = copy.deepcopy(self)
        #obj.matrix = obj.matrix.astype(np.uint8)
        
        filter_labels = []
        regions = regionprops(obj.matrix)
        for index, props in enumerate(regions):
            if props.minor_axis_length <= 8 and (props.minor_axis_length * 5 
                    <= props.major_axis_length):
                # abnormity cell with large aspect ratio
                filter_labels.append(props.label)
                continue
            if props.area > 1000 or props.area < 6:
                # extreme large cell caused by non-detected blur region
                # extreme small cell original segmentation fault
                filter_labels.append(props.label)
                continue
            if props.extent < 0.3:
                filter_labels.append(props.label)
                continue
            if props.minor_axis_length < min_size:
                # extreme thin cell
                filter_labels.append(props.label)
                continue
            if max_size and props.major_axis_length > max_size:
                # extreme fat cell
                filter_labels.append(props.label)
                continue

        mask = np.isin(obj.matrix, filter_labels)
        obj.matrix[mask] = 0

        sys.stdout.write('{} labels removed\n'.format(len(filter_labels)))
        return obj

    def merge(self, other, how='left'):
        
        sys.stdout.write('merge mix labels ... ')
        sys.stdout.flush()

        if how == 'left':
            obj = copy.deepcopy(self)
            mask1 = obj.to_binary()
            mask2 = copy.deepcopy(other)
        elif how == 'right':
            obj = copy.deepcopy(other)
            mask1 = obj.to_binary()
            mask2 = copy.deepcopy(self)
        else:
            pass

        intersection = cv2.bitwise_and(mask1.matrix, mask2.matrix)

        mask2.matrix[intersection] = 0

        obj.matrix += mask2.matrix

        sys.stdout.write('done\n')
        return obj
    
    def save(self, outfile='out.txt'):
        if outfile.endswith(('.tif', '.png')):
            cv2.imwrite(outfile, self.matrix)
        else:
            np.savetxt(outfile, self.matrix, fmt='%d')
        return 

    def overlayoutlines(self, image=None, prefix=None):

        sys.stdout.write('draw outlines ... ')
        sys.stdout.flush()
        
        import skimage.io
        import skimage.segmentation

        if isinstance(image, str):
            image = skimage.io.imread(image)

        outlines = skimage.segmentation.mark_boundaries(
                image, 
                self.matrix, 
                color=(1, 0, 0),
                mode='inner',
                )
        b, g, r = cv2.split(outlines)

        sys.stdout.write('{} labels\n'.format(len(np.unique(self.matrix))))

        mask = np.isin(b, [1])
        image[mask] = 255
        
        if prefix:
            np.savetxt(f'{prefix}.outlines.txt', b, fmt='%d')
            cv2.imwrite(f'{prefix}.outlines.png', image)
        return b, image

    def affine_transform(self, affineR, reference=None, order=0, outfile=None):
        import scipy.ndimage
        
        obj = copy.deepcopy(self)

        if reference is not None:
            shape = reference.T.shape
        else:
            shape = None

        affined_matrix = scipy.ndimage.affine_transform(
                obj.matrix.T, 
                affineR, 
                output_shape=shape, 
                order=order
                )
        obj.matrix = affined_matrix.T
        
        if outfile:
            obj.save(outfile)

        return obj

    def crop(self, mask=None, outfile=None):
        obj = copy.deepcopy(self)
        
        if mask is None:
            mask = obj.matrix != 0
        else:
            if isinstance(mask, Mask):
                mask = mask.matrix
            assert obj.matrix.shape == mask.shape
            mask = mask != 0
        coords = np.argwhere(mask)
        xs, ys = zip(*coords)
        x0, y0 = coords.min(axis=0)
        x1, y1 = coords.max(axis=0) + 1
        
        obj.matrix = obj.matrix[x0:x1, y0:y1]
        
        if outfile:
            obj.save(outfile)

        return obj

    def plot(self, 
            label_map=None, 
            coord_map=None, 
            outlines=None, 
            cmap=None, 
            backgroud=None, 
            outfile=None, 
            dx=715, units='nm', 
            robust=False,
            perc=(2,98),
            interpolation='none',
            ax=None,
            ):

        background_mask = np.isin(self.matrix, [0])
        
        if label_map is not None:
            if not isinstance(label_map, pd.DataFrame):
                label_map = pd.read_csv(label_map, sep='\t', header=None, 
                        names=['cell', 'mask'], index_col='cell')
            self = self.relabel(label_map=label_map)
        elif coord_map is not None:
            self = self.relabel(coord_map=coord_map)
        
        if outlines:
            if isinstance(outlines, str):
                outlines = np.loadtxt(outlines)
            elif isinstance(outlines, Mask):
                outlines = outlines.matrix
            mask = np.where(outlines == 1)
            self.matrix[mask] = 0
        
        import matplotlib.pyplot as plt
        from matplotlib_scalebar.scalebar import ScaleBar
        import seaborn_image as isns
        #from skimage.exposure import adjust_gamma

        if ax is None:
            fig, ax = plt.subplots()
            
        if not cmap:
            cmap = 'deep'

        ax = isns.imgplot(
                self.matrix, 
                ax=ax,
                robust=robust, 
                perc=perc, 
                cmap=cmap, 
                despine=True,
                #map_func=adjust_gamma, 
                #gamma=0.5
                origin='upper',
                interpolation=interpolation,
                )

        scalebar = ScaleBar(
                dx=dx, 
                units=units, 
                fixed_value=50, 
                fixed_units='um'
                )
        ax.add_artist(scalebar)
        #ax.invert_yaxis()

        if outfile:
            fig.savefig(outfile, dpi=600)
        return ax

    @property
    def centroid(self):

        from skimage.measure import regionprops_table

        label_mask = copy.deepcopy(self.matrxi)
        assert len(np.unique(label_mask)) > 2, 'must be label mask to find centroid'

        zero_mask = np.isin(label_mask, [0])
        label_mask[zero_mask] = label_mask.max() + 1

        unique, unique_indices, unique_inverse = np.unique(
                label_mask, 
                return_index=True, 
                return_inverse=True
                )

        relabel_mask = unique_inverse.reshape(label_mask.shape)
        min_mask = np.isin(relabel_mask, [0])
        relabel_mask[min_mask] = relabel_mask.max() + 1
        relabel_mask[zero_mask] = 0
    
        unique[-1] = 0
        label_map = pd.DataFrame(dict(
                orig_label=unique, 
                label=unique_inverse[unique_indices]
                ))
    
        props = regionprops_table(
                relabel_mask,
                properties=('label', 'centroid')
                )

        df = pd.DataFrame(props)
        df = df.merge(label_map, how='left', on='label')
        df = df[['orig_label', 'centroid-1', 'centroid-0']]
        df = df.rename(
                columns={
                    'orig_label':'label', 
                    'centroid-0':'y',
                    'centroid-1':'x'
                    }
                )
        df['label'] = df['label'].astype(int)
        return df

