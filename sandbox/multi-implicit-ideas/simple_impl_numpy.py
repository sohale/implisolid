import numpy as np
# from tkinter import Y


def signed_pow(x, p):
    absx = np.abs(x)
    sgnx = np.sign(x)
    return np.power(absx, p) * sgnx

def sdf1(X,Y):
    return signed_pow((X - 0.5)**2 + (Y-0.5)**2, 0.5) - 0.5


def sdf2(X,Y):
    return signed_pow((X - 0.0)**2 + (Y-0.5)**2, 0.5) - 0.5

def mesh():
    x = np.linspace(-2, 2, num=900)
    y = np.linspace(-2, 2, num=1100)
    X, Y = np.meshgrid(x, y)
    return X,Y

def pl1(img_vals):
    from matplotlib import pyplot
    pyplot.imshow(img_vals,)
    pyplot.colorbar()
    pyplot.show()

def pl3(img_vals):
    # cmap_name = 'Wistia'  # Wistia: Sun
    # see: https://matplotlib.org/3.5.0/tutorials/colors/colormaps.html
    cmap_name = 'copper'  # Wistia: Sun
    from matplotlib import pyplot
    cmap=pyplot.get_cmap(cmap_name)
    pyplot.imshow(img_vals, cmap=cmap)
    pyplot.colorbar()
    pyplot.show()

def pl4(img_vals, more_info):
    # cmap_name = 'Wistia'  # Wistia: Sun
    # see: https://matplotlib.org/3.5.0/tutorials/colors/colormaps.html
    cmap_name = 'copper'  # Wistia: Sun
    from matplotlib import pyplot
    cmap=pyplot.get_cmap(cmap_name)

    X,Y = more_info
    # extent = (-3, 4, -4, 3)
    extent = (np.min(X.ravel()), np.max(X.ravel()), np.min(Y.ravel()), np.max(Y.ravel()))

    levels = np.linspace(-1, 5, num=20)
    pyplot.imshow(img_vals, cmap=cmap, extent=extent, origin='lower')
    #origin='lower'
    pyplot.colorbar()
    pyplot.contour(X,Y, img_vals, levels=levels, colors='r' , origin='image') #, extent=extent) # cmap=cmap)
    # origin='upper' deafult: lower?
    # https://matplotlib.org/stable/gallery/images_contours_and_fields/contour_image.html
    pyplot.show()

def pl2(img_vals):
    from matplotlib import pyplot
    import matplotlib as mpl

    #cmap = mpl.cm.cool
    #norm = mpl.colors.Normalize(vmin=5, vmax=10)
    #fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
    #         cax=ax, orientation='horizontal', label='Some Units')

    #fig, ax = pyplot.subplots(figsize=(6, 1))
    #fig.subplots_adjust(bottom=0.5)
    # https://matplotlib.org/3.5.0/tutorials/colors/colorbar_only.html
    cmap = (mpl.colors.ListedColormap(['royalblue', 'cyan', 'yellow', 'orange'])
        .with_extremes(over='red', under='blue'))
    bounds = [-1.0, -0.5, 0.0, 0.5, 1.0]
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    pyplot.colorbar( mpl.cm.ScalarMappable(cmap=cmap, norm=norm),
        ticks=bounds,
        boundaries=[-10] + bounds + [10],
        #cax=ax,
        )

    pyplot.plot([1,2], [1,4])
    ih = pyplot.imshow(img_vals,
    #cmap=pyplot.get_cmap(name)
    cmap=cmap
    )

    # pyplot.colorbar()
    # pyplot.show()

    pyplot.show()

def main():
    X,Y = mesh()
    v1 = sdf1(X,Y)
    v2 = sdf2(X,Y)
    v = np.minimum(v1,v2)
    # v = np.maximum(v1,v1*0)
    v[v<0] = -1
    pl4(v, (X,Y))
    from matplotlib import pyplot




main()
