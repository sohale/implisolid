import numpy as np
# from tkinter import Y


def signed_pow(x, p):
    absx = np.abs(x)
    sgnx = np.sign(x)
    #return np.power(absx, p) * sgnx, sgnx
    return np.power(x, p), None

def sdf1(X,Y):
    v = signed_pow((X - 0.5)**2 + (Y-0.5)**2, 0.5)[0] - 0.2
    (gx,gy) = 2*(X - 0.5),  2*(Y-0.5)
    return v, (gx,gy)


def sdf2(X,Y):
    v = signed_pow((X - 0.0)**2 + (Y-0.5)**2, 0.5)[0] - 0.5
    #, (gx,gy)
    (gx,gy) = 2*(X - 0.0),  2*(Y-0.5)
    return v, (gx,gy)

from math import prod
def apply_third_index(min12, v12_list):
    # todo: take and take_along_axis
    axis_ = 2
    print('min12.shape', min12.shape)  # min12.shape (1100, 900)
    # vr = v12.reshape((axis_, prod(list(min12.shape))))
    assert len(v12_list) > 0

    v120 = v12_list[0]
    ne = prod(list(v120.shape))
    num_choices = v120.shape[axis_]
    aea = int(ne / num_choices)
    min12r = min12.reshape((aea, ))

    for v12a in v12_list:
        assert tuple(list(min12.shape) + [num_choices]) == v12a.shape
        print('>', v12a.shape)

    main_shape = min12.shape

    print('outlist')
    outlist = []
    for v12a in v12_list:
        # vr = v12a.reshape((axis_, aea))  # flatten all dims except axis_
        vr = v12a.reshape((aea, axis_))  # flatten all dims except axis_
        print('vr.shape',vr.shape)
        print('min12r.shape',min12r.shape)
        print('slow')
        # vr.shape (990000, 2)
        # min12r.shape (990000,)
        # voutr = vr[:, min12r]  # the core operation

        def little_test():
            a = np.array([[1,2,3,4],[10,20,30,40]]).T
            print(a)
            print(a.shape)  # (4,2)
            #print(a[[0,1,2]])
            #print(a[:,[0,0,1]])
            #print(a[[0,0,1]])
            #print(a[0,1], a[1,0])
            ind = np.array([[0,1,0,1]]).T
            print(ind.shape)  # 4x1
            print(np.take_along_axis(a,ind, axis=1))  # shape: (4,1)
        #little_test()
        #exit()
        print(min12r.shape)
        # take
        # put_along_axis
        # ravel_multi_index
        #voutr = np.take_along_axis(vr, min12r, axis=1)  # the core operation
        #voutr = np.take(vr, min12r, axis=1)  # the core operation
        voutr = np.take_along_axis(vr, min12r[:,None], axis=1)  # the core operation

        print('slow done')
        print('voutr.shape', voutr.shape)
        print('main_shape', main_shape)
        voutu = voutr.reshape(main_shape)
        outlist.append(voutu)
    for voutu in outlist:
        print('voutu.shape', voutu.shape)
    return tuple(outlist)


def im_sdf(X,Y):
    v1,(g1x,g1y) = sdf1(X,Y)
    v2,(g2x,g2y) = sdf2(X,Y)
    print('00')
    #v = np.minimum(v1,v2)
    v12 = np.concatenate((v1[:,:,None],v2[:,:,None],),axis=2)
    print('11')
    min12 = np.argmin(v12, axis=2)
    print('22', min12.shape)
    print('22', X.shape)
    print('22', min12)

    print('g1.shape', g1x.shape)
    g12x = np.concatenate((g1x[:,:,None],g2x[:,:,None],),axis=2)
    g12y = np.concatenate((g1y[:,:,None],g2y[:,:,None],),axis=2)

    print('g1.apply_third_index')
    v, gx, gy = apply_third_index(min12, (v12, g12x, g12y))

    #v,gx,gy = v1,g1x,g1y
    #v,gx,gy = v2,g2x,g2y

    # #v = v12[:,:,min12]
    # v = v12[:,:],min12
    # gx,gy = None,None

    yield v, (gx,gy)

    near_balance = np.abs(v1-v2) < 0.1 # 0.02
    #v[~near_balance] = np.nan
    #mgx # multiple gx
    mv = v12.copy()
    mgx = g12x.copy()
    mgy = g12y.copy()
    mv[~near_balance] = np.nan
    mgx[~near_balance] = np.nan
    mgy[~near_balance] = np.nan
    yield mv, (mgx,mgy)

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

    X,Y,gx,gy = more_info
    # extent = (-3, 4, -4, 3)
    extent = (np.min(X.ravel()), np.max(X.ravel()), np.min(Y.ravel()), np.max(Y.ravel()))

    #levels = np.linspace(-1, 5, num=20)
    levels = np.linspace(-1, 5, num=10)
    pyplot.imshow(img_vals, cmap=cmap, extent=extent, origin='lower')
    #origin='lower'
    pyplot.colorbar()
    pyplot.contour(X,Y, img_vals, levels=levels, colors='r' , origin='image') #, extent=extent) # cmap=cmap)
    # origin='upper' deafult: lower?
    # https://matplotlib.org/stable/gallery/images_contours_and_fields/contour_image.html
    #pyplot.show()

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

def plot_grads(X,Y, gx,gy):
    from matplotlib import pyplot

    gn = np.power(gx**2 + gy**2, 0.5)
    gx = gx/gn * 0.1
    gy = gy/gn * 0.1
    # rotate!
    gx,gy = gy, -gx

    pyplot.streamplot(X,Y,gx,gy, \
        density=3, color='k', \
        maxlength=0.4, linewidth=0.2, \
        arrowsize=1, arrowstyle='->'
    )
    # https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.FancyArrowPatch.html#matplotlib.patches.FancyArrowPatch

def main():
    X,Y = mesh()

    # v = np.maximum(v1,v1*0)
    #v = np.maximum(v1,v1*0)
    g = im_sdf(X,Y)

    v, (gx,gy) = next(g)
    v[v<0] = -1
    print('vv')
    pl4(v, (X,Y, gx, gy))

    plot_grads(X,Y, gx,gy)


    from matplotlib import pyplot
    pyplot.show()

    pyplot.figure()
    mv, (mgx,mgy) = next(g)
    print(mv.shape, mgx.shape,mgy.shape)
    pl4(mv[:,:,0], (X,Y, mgx[:,:,0], mgy[:,:,0]))
    for i in range(mgx.shape[2]):
        v,gx,gy = mv[:,:,i], mgx[:,:,i], mgy[:,:,i]
        plot_grads(X,Y, gx,gy)
    pyplot.show()

main()
