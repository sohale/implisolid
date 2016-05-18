import numpy as np

import sys
print"Mayavi: importing..", ; sys.stdout.flush()
from mayavi import mlab
print "Imported."; sys.stdout.flush()

def noisy(v, ampl):
    v = v + (np.random.rand(v.shape[0],v.shape[1])*2.-1.) * ampl
    return v

from basic_types import normalize_vector4_vectorized
from basic_types import check_vector4_vectorized

def visualise_gradients(mlab, pos, iobj, arrow_size):
    lm = arrow_size  # 1.  # STEPSIZE
    pos4 = np.concatenate((pos, np.ones((pos.shape[0],1))),axis=1)
    pnormals = - iobj.implicitGradient(pos4)
    pnormals = normalize_vector4_vectorized(pnormals)
    check_vector4_vectorized(pos4)
    xyz = pos
    uvw = pnormals [:,0:3] / 2.
    xx, yy, zz = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    uu, vv, ww = uvw[:, 0], uvw[:, 1], uvw[:, 2]
    #ax.quiver
    #ax.quiver(xx, yy, zz,   uu, vv, ww,  length=np.abs(lm), arrow_length_ratio=0.3, alpha=0.3, pivot="tail")
    #arrow_length_ratio=   length=np.abs(lm)
    #pivot: tail | middle | tip
    #mlab.quiver3d(x_verts,y_verts,z_verts, UVW_normals[:,0],UVW_normals[:,1],UVW_normals[:,2],color=(0,0,0))
    mlab.quiver3d(xx, yy, zz, uu, vv, ww, color=(0, 0, 0), scale_factor=np.abs(lm), line_width=0.5)


def visualise_displacements(mlab, verts_from, verts_to):
    n = verts_from.shape[0]
    assert verts_from.shape == (n, 3)
    assert verts_to.shape == (n, 3)

    xyz1 = verts_from
    xyz2 = verts_to
    set_trace()
    ddd = xyz2 - xyz1  #displacement
    arrow_size = np.mean(np.linalg.norm(ddd, axis=1)) * 0.3

    xx, yy, zz = xyz1[:, 0], xyz1[:, 1], xyz1[:, 2]
    uu, vv, ww = ddd[:, 0], ddd[:, 1], ddd[:, 2]
    mlab.quiver3d(xx, yy, zz, uu, vv, ww, color=(0, 0, 0), scale_factor=arrow_size)
    #line_width=0.5


def display_simple_using_mayavi_2(vf_list, pointcloud_list, minmax=(-1,1), mayavi_wireframe=False, opacity=1.0,
        separate_panels=True, gradients_at=None, gradients_from_iobj=None, pointsizes=None, pointcloud_opacity=1.,
        add_noise=[], noise_added_before_broadcast=False, fromto=None, random_interior_point=False, labels=None, grad_arrow_len=None):
    """Two separate panels"""

    #print"Mayavi: importing..", ; sys.stdout.flush()

    #from mayavi import mlab
    #print "Imported."; sys.stdout.flush()

    if pointsizes is None:
        pointsizes = [0.2]*10

    if type(opacity) is list:
        opacities = opacity  # 1.0
    else:
        opacities = [opacity] + [0.2]*(len(vf_list)-1)  # 1.0, 0.2 #0.1

    assert len(opacities) == len(vf_list)

    for fi in range(len(vf_list)):
        if separate_panels:
            mlab.figure()

        #allpoints are plottedon all panels?
        color_list = [(1, 0, 0), (0, 0, 0), (1, 1, 0), (0, 0, 1), (0,1,0)]
        i = 0
        for c in pointcloud_list:
            #if separate_panels:
            #    if i != fi:
            #        continue
            mlab.points3d(c[:, 0], c[:, 1], c[:, 2], color=color_list[i], scale_factor=pointsizes[i], opacity=pointcloud_opacity)
            i += 1
        del i

        vf = vf_list[fi]
        verts, faces = vf

        if verts is None:
            continue
        if verts.size == 0:
            print("Warning: empty vertices")
            continue
        if faces.size == 0:
            print("Warning: no faces")
            continue

        assert verts.ndim == 2
        assert faces.ndim == 2
        assert verts.shape == (verts.shape[0], 3), str(verts.shape)
        assert faces.shape == (faces.shape[0], 3), str(faces.shape)
        if type(mayavi_wireframe) is list:
            wire_frame1 = mayavi_wireframe[fi]
            assert len(mayavi_wireframe) == len(vf_list)
        else:
            wire_frame1 = mayavi_wireframe
        M = 0.1*0
        if len(add_noise)>0:
            assert len(add_noise) == len(vf_list)
            M = float(add_noise[fi])

        ignore_noise = False
        if ignore_noise:
            _v = verts
            _f = faces
        else:
            if noise_added_before_broadcast:
                pre_expand_noise = True
                post_expand_noise = False
            else:
                pre_expand_noise = False
                post_expand_noise = True

            if pre_expand_noise:
                print "noise", M/2.
                verts = noisy(verts, M/2.)

            _v = verts[faces, :]
            #print _v.shape

            _nv = np.prod(faces.shape)
            _v = _v.reshape( (_nv, 3) )
            assert _nv % 3 == 0
            _f = np.arange(_nv).reshape( (_nv/3, 3) )

            qq = _v[_f,:]

            #print np.nonzero( np.isnan(_f.ravel()) )
            #print np.nonzero( np.isinf(_f.ravel()) )
            vv=_f.ravel().copy(); vv.sort()
            #print np.nonzero(vv != np.arange(vv.size))

            #_v = np.concatenate( (_v, np.zeros( (10000, 3) )), axis=0)
            if post_expand_noise:
                #_v = _v + (np.random.rand( _v.shape[0], _v.shape[1] ) -0.5) * M
                _v = noisy(_v, M/2.)
                print "noise2 ", M/2.

            #qq = _v[_f, :]
            #assert  np.all( _f.ravel() == np.arange( _f.size ) )
            #print qq
            #print _v
            #assert np.all(qq.ravel() == _v.ravel())

        #v1 = [_v[0]+M*np.random.rand() for vert in verts],
        #v2 = [_v[1]+M*np.random.rand() for vert in verts],
        #v3 = [_v[2]+M*np.random.rand() for vert in verts],
        v1 = [v[0] for v in _v]
        v2 = [v[1] for v in _v]
        v3 = [v[2] for v in _v]
        mlab.triangular_mesh(
                        v1, v2, v3, _f,
                        representation="surface" if not wire_frame1 else "wireframe",
                        opacity=opacities[fi], scale_factor = 100.0)
        #opacity = 0.2 #0.1

    #del fi
    del vf_list
    if True:
        if minmax is not None:
            (RANGE_MIN, RANGE_MAX) = minmax
            x = np.linspace(RANGE_MIN,RANGE_MAX,2).reshape(2,1)
            y = np.zeros((2,1))
            z = np.zeros((2,1))

            mlab.plot3d(x, y, z,line_width=3,name="x-axis")
            mlab.plot3d(y, x, z,line_width=3,name="y-axis")
            mlab.plot3d(z, y, x,line_width=3,name="z-axis")

            mlab.text3d(RANGE_MAX,0,0, "x", scale=0.3)
            mlab.text3d(0,RANGE_MAX,0, "y", scale=0.3)
            mlab.text3d(0,0,RANGE_MAX, "z", scale=0.3)
            mlab.text3d(RANGE_MIN,0,0, "-x", scale=0.3)
            mlab.text3d(0,RANGE_MIN,0, "-y", scale=0.3)
            mlab.text3d(0,0,RANGE_MIN, "-z", scale=0.3)

        visualise_centroid_ids = False
        if visualise_centroid_ids:
            for fii in range(faces.shape[0]):
                if fii in [65, 274, 310, 362]:
                    f3 = faces[fii, :]
                    v123 = verts[f3, :]
                    centroid = np.mean(v123, axis=0)
                    mlab.text3d(centroid[0],centroid[1],centroid[2], str(fii), scale=0.1 * 0.2)

        if not labels == None:
            (labels_xyz, labels_ids) = labels
            if labels_ids.size > 400:
                print "**** SKIPPING LABELS *****"
                vii = labels_ids[0]
                v123 = labels_xyz[vii, :]
                text = "** SKIPPING LABELS (>400) **"
                mlab.text3d(v123[0], v123[1], v123[2], text, scale=0.03*5)
            else:
                for vii in range(labels_xyz.shape[0]):
                    if vii in labels_ids:
                        v123 = labels_xyz[vii, :]
                        mlab.text3d(v123[0], v123[1], v123[2], str(vii), scale=0.03)

    def add_random_interior_points(ax, iobj, avg_edge_len):
        """ Adding random points """
        n = 10000
        import basic_types
        # ******
        #print avg_edge_len, "WHY USED BEFORE DEFINED?"
        ampl = avg_edge_len
        #ampl = 2
        x = basic_types.make_random_vector_vectorized(n, ampl, 1, type="rand", normalize=False)
        v = iobj.implicitFunction(x)
        x_sel =  x[ v >= 0 , :]
        if x_sel.size ==0:
            print("No points")
            return
        ax.points3d(x_sel[:,0], x_sel[:,1], x_sel[:,2], color=(0,0,0), scale_factor=0.2)

    if gradients_at is not None:
        if grad_arrow_len is None:
            verts1, faces1 = vf_list[0]
            avg_edge_len = compute_average_edge_length(verts, faces)
            print "avg_edge_len", avg_edge_len
            grad_arrow_len = avg_edge_len / 20. *10.
        visualise_gradients(mlab, gradients_at, gradients_from_iobj, grad_arrow_len)
    if gradients_from_iobj is not None:
        if random_interior_point:
            add_random_interior_points(mlab, gradients_from_iobj, avg_edge_len)

    if not fromto == None:
        (verts_from, verts_to) = fromto
        assert verts_from.shape[1] == 3
        assert verts_to.shape[1] == 3
        assert verts_to.shape[0] < 10000
        visualise_displacements(mlab, verts_from, verts_to)
    mlab.show()
    return


def quick_vis(verts, facets, face_idx):
    #print facets.shape
    #print facets[face_idx, :].shape
    display_simple_using_mayavi_2( [(verts, facets), (verts, facets[face_idx, :]), (verts, facets[face_idx, :])],
       pointcloud_list=[ ], pointcloud_opacity=0.2,
       mayavi_wireframe=[True, True, True], opacity=[0.1, 1, 0.2],
       gradients_at=None, separate_panels=False, gradients_from_iobj=None, add_noise=[0, 0, 0.01],
       noise_added_before_broadcast=True
       )


# for visualisation only
"""
def compute_average_edge_length(verts, faces):
    nfaces = faces.shape[0]
    expand = verts[faces, :]
    assert expand.shape == (nfaces, 3, 3)
    assert expand[:, 2, :].shape == (nfaces, 3)
    ea_sum = 0.
    for i in range(3):
        i1 = i
        i2 = (i+1) % 3
        e1 = np.linalg.norm(expand[:, i1, :] - expand[:, i2, :], axis=1)  # bug fixed! axis=1 was missing.
        assert e1.shape == (nfaces,)
        ea_sum += np.mean(e1)

    #set_trace()
    return ea_sum / 3.
"""