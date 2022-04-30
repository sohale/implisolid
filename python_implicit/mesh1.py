import mesh_utils
import numpy as np
from basic_types import check_vector4_vectorized, normalize_vector4_vectorized

SUPPRESS_DEBUG_PRINT = True


class Mesh_1(object):
    def __init__(self, faces, verts):
        self.faces = faces
        self.verts = verts
        self.centroids = None
        self.vertex_neighbours_list = None
        self.centroid_gradients = None
        self.facet_areas = None

    def build_centroids(self):
        self.calculate_face_areas()

        expand = self.verts[self.faces, :]
        nfaces = self.faces.shape[0]
        print(expand.shape)
        assert expand.shape == (nfaces, 3, 3)
        centroids3 = np.mean(self.verts[self.faces[:], :], axis=1)
        print(centroids3)
        if False:
            """ Check if there are repeats in the centroids """
            repeats = 0
            for i in range(centroids3.shape[0]):
                for j in range(centroids3.shape[0]):
                    if i == j:
                        continue
                    dif = centroids3[i, :] - centroids3[j, :]
                    if np.max(np.abs(dif)) < 0.0001:
                        # print("eq ", centroids3[i,:], centroids3[j,:], "faces: ", i,j)
                        # print(self.verts[self.faces[i,:],:], self.verts[self.faces[i,:],:] )
                        repeats += 1
            if repeats > 0:
                print("repeats: ", repeats)

        n = self.faces.shape[0]
        self.centroids = np.concatenate((centroids3, np.ones((n, 1))), axis=1)
        # print(self.centroids.shape)
        self.centroid_gradients = None

    def evaluate_centroid_gradients(self, iobj):
        assert self.centroids is not None
        n = self.centroids.shape[0]
        # centroids4 = np.concatenate( (self.centroids, np.ones((n,1))), axis=1)
        centroids4 = self.centroids
        check_vector4_vectorized(centroids4)
        self.centroid_gradients = iobj.implicitGradient(centroids4)
        assert not np.any( np.isnan(self.centroid_gradients) )
        assert not np.any(np.isinf(self.centroid_gradients))
        self.centroid_normals = normalize_vector4_vectorized(self.centroid_gradients)
        # self.centroid_normals used only in (Laplacian) weighted resampling

    def calculate_face_areas(self):
        nfaces = self.faces.shape[0]
        expand = self.verts[self.faces, :]
        assert expand.shape == (nfaces, 3, 3)
        assert expand[:, 2, :].shape == (nfaces, 3)
        a = np.cross(expand[:, 1, :] - expand[:, 0, :], expand[:, 2, :] - expand[:, 0, :], axis=1)
        self.facet_areas = np.linalg.norm(a, axis=1, ord=2) / 2.0
        # print(list(-np.sort(-self.facet_areas)))
        # print(self.facet_areas[self.facet_areas < 0.00001])
        degenerates_count = len(self.facet_areas[self.facet_areas < 0.00001])   # 9
        print(degenerates_count)

    def build_neighbours(self):
        self.vertex_neighbours_list = mesh_utils.make_neighbour_faces_of_vertex(self.faces)

    # high level functions:
    def correct_centroids(self):
        # calls optimize_centroid
        # todo
        pass

    def get_A_b(self, vertex_id):
        nlist = self.vertex_neighbours_list[vertex_id]
        nai = np.array(nlist)
        center_array = self.centroids[nai, :]
        # print(center_array)
        # self.evaluate_centroid_gradients(iobj)
        # print("nai", nai)
        # print(self.centroid_gradients.shape)

        # note some centers may not be projected successfully in the previous step
        not_projected_successfully = np.isnan(center_array[0])
        if np.any(not_projected_successfully):
            pass

        # todo:
        # normals = self.centroid_normals[nai,:]
        normals = self.centroid_gradients[nai, :]  # why do we have repeats??
        # note : not normalised. But it works.

        norms = np.linalg.norm(normals, ord=2, axis=1)
        # can be either 0, 1 or Nan
        if np.any(norms < 0.000001):    # can be exactly 0.0
            print("Error: bad normal", normals)

        # can be 0,0,0, inf, nonsharp, degenerate, ...
        degenerate_normals = np.logical_or(np.isnan(np.sum(normals, axis=1)), norms < 0.0000001)
        # simpler: degenerate_normals = np.logical_or(np.isnan(norms), norms < 0.0000001 )
        # todo:



        #print(normals)
        assert not np.any( np.isnan(normals) )
        assert not np.any( np.isinf(normals) )

        #normals = normalize_vector4_vectorized( normals ) #todo: move it to evaluate_centroid_gradients or self.centroid_normals

        #print("normals", normals) # can be all 0,0,0

        pa = center_array
        na = normals   # (4)x4
        #grad = Ax+b
        A = np.zeros((3,3))
        b = np.zeros((3,1))
        #assert len(pa) == len(na)
        assert na.shape == pa.shape
        for i in range(na.shape[0]):
            #assert na[i]***
            n_i = na[i,0:3, np.newaxis]
            assert n_i.shape == (3,1)
            #print("n_i: ", n_i, n_i.shape)
            nnt = np.dot(n_i, np.transpose(n_i))
            #print("nnt",nnt)
            assert nnt.shape == (3,3)
            A += nnt
            #A can contain equal rows. It's fine. In this case, we have faces that are parallel or on the same plane (e.g. on the same side of a cube)
            p_i = pa[i, 0:3, np.newaxis]
            assert p_i.shape == (3,1)
            b += -np.dot( nnt, p_i)

        return A, b

    def update_centroids_and_gradients(self, iobj, lambda_=0.5, maxdist=np.infty, method="v0.1_inplace"):
        tolerance = 0.00001
        #todo: lambda = <edge length> / 2
        #lambda_ = 0.5
        #self.centroids =
        #mesh_utils.project_points2(iobj, self.centroids, tolerance = tolerance, lambda_=lambda_, maxdist=maxdist)
        #xa =
        if method=="v0.1_inplace":
            xa = mesh_utils.optimize_points1_inplace_oldworking(iobj, self.centroids, tolerance = tolerance, lambda_=lambda_, maxdist=maxdist,
                inplace=True)
        elif method=="v0.1_not_inplace":
            xa = mesh_utils.optimize_points1_inplace_oldworking(iobj, self.centroids, tolerance = tolerance, lambda_=lambda_, maxdist=maxdist,
                inplace=False)
        else:
            raise NotImplementedError()
        self.unsuccessfully_projected_centroids = iobj.implicitFunction(self.centroids) <= tolerance
        self.evaluate_centroid_gradients(iobj)  # rendudant calculation

    #identical to the old version
    def quadratic_optimise_vertices(self, alpha=1.0):
        print(("="*10+"\n")*7)
        assert not self.centroids is None
        assert not self.vertex_neighbours_list is None
        assert not self.centroid_gradients is None

        #qem_utils.quadratic_optimise_vertices(iobj,)

        self.new_verts = self.verts*0 -1

        nvert = len(self.vertex_neighbours_list)

        self.verts_ranks = np.zeros((nvert,))

        assert nvert == self.verts.shape[0]
        for vi in range(nvert):
            A, b = self.get_A_b(vi)

            #print A
            #print(vi)
            #print(self.facet_areas[vi])
            u, s, v = np.linalg.svd(A)
            assert np.allclose(A, np.dot(u, np.dot(np.diag(s), v)))
            #print(s)  # [  1.48148148e+01   1.67928330e-15   1.01592270e-50]
            assert s[0] == np.max(s)
            #print( s / s[0] )  # [  1.00000000e+00   1.13351623e-16   6.85747820e-52]

            tau = 10. ** 3.
            s[s / s[0] < 1.0/tau] = 0
            #print(s , s[0] , tau)
            rank = np.sum(s / s[0] > 1.0/tau)
            #print(s)
            #print("rank = ", rank)

            #rank will never be 0: s[0]/s[0] is always 1, even when s[0] is too small.
            #assert s[0] > 0.000001

            if not  s[0] > 0.000001:
                print("Warning! sigma_1 == 0" )
                print(s)
                print("A", A)

                #not tested
                self.verts_ranks[vi] = 0
                self.new_verts[vi,0:3] = new_x[:,0]

            assert np.all(s[:rank]/s[0] >= 1.0/tau)


            #x
            #x= self.centroids[fi] [0:3, np.newaxis]
            x= self.verts[vi, 0:3, np.newaxis]
            assert x.shape == (3, 1)

            y = np.dot(v,x).copy()
            utb = np.dot(-np.transpose(u), b)
            #print("rank", rank, "1/tau=", 1./tau)
            #print s
            for i in range(rank):
                #print(np.dot(-np.transpose(u), b), "scalar")
                assert np.dot(-np.transpose(u), b).shape == (3,1)
                #print s[i] , 1.0/tau
                #assert s[i] >= 1.0/tau #fails when s[0] is small
                assert s[i]/s[0] >= 1.0/tau
                y[i] = utb[i] / s[i]
            new_x = np.dot(np.transpose(v), y)
            #print(x.ravel(), " -> ", new_x.ravel())
            #print("    delta=", (new_x - x).ravel())

            self.new_verts[vi,0:3] = new_x[:,0]
            #self.new_verts[vi,3] = 1

            assert x.shape == (3, 1)
            #alpha = 1
            self.new_verts[vi,0:3] = new_x[:,0] * alpha + x[:,0] * (1.0-alpha)

            if not np.all(np.abs(utb.ravel()[rank:] ) < 0.0001):
                #print("s", s.ravel()/s[0], "   utb", utb.ravel()/s[0])
                pass
            self.verts_ranks[vi] = rank

            #exit()
        print("max rank = ", np.max(self.verts_ranks))
        print("min rank = ", np.min(self.verts_ranks))
        if not np.min(self.verts_ranks) >= 1:
            print("Warning: assertion: np.min(self.verts_ranks) >= 1 failed." )

        if False:
            assert np.min(self.verts_ranks) >= 1

    # Some tests failed. so discontinued for now:
    def quadratic_optimise_vertices_new(self, alpha=1.0):
        assert False
        assert not self.centroids is None
        assert not self.vertex_neighbours_list is None
        assert not self.centroid_gradients is None
        self.new_verts = self.verts*0 -1

        nvert = len(self.vertex_neighbours_list)

        self.verts_ranks = np.zeros((nvert,))

        assert nvert == self.verts.shape[0]
        for vi in range(nvert):
            #old_x = self.centroids[vi, :]
            #old_x = self.verts[vi, 0:3, np.newaxis]

            A, b = self.get_A_b(vi)

            #print A
            #print(vi)
            #print(self.facet_areas[vi])
            u, s, v = np.linalg.svd(A)
            assert np.allclose(A, np.dot(u, np.dot(np.diag(s), v)))
            #print(s)  # [  1.48148148e+01   1.67928330e-15   1.01592270e-50]
            assert s[0] == np.max(s)
            #print( s / s[0] )  # [  1.00000000e+00   1.13351623e-16   6.85747820e-52]

            tau = 10. ** 3.

            old_x= self.verts[vi, 0:3, np.newaxis]
            #old_x
            #old_x= self.centroids[fi] [0:3, np.newaxis]
            old_x= self.verts[vi, 0:3, np.newaxis]
            assert old_x.shape == (3, 1)

            nihl0 = s[0]<0.0000000001
            if nihl0:
                assert np.all(s < 0.00001)
                rank = 0
                s = s * 0
                new_x = old_x

            else:

                if np.any(np.isnan(s/s[0])): # or np.any(np.isinf(s/s[0])):
                    #print("NAN")
                    print('nan', s)
                if not s[0] < tau*1000:
                    if not  SUPPRESS_DEBUG_PRINT:
                        print('s[0] not< big', s) #i.e. too big

                s[s / s[0] < 1.0/tau] = 0
                #print(s , s[0] , tau)
                rank = np.sum(s / s[0] > 1.0/tau)
                #print(s)
                #print("rank = ", rank)


                old_y = np.dot(v,old_x)
                #y = old_y
                y = old_y.copy()
                utb = np.dot(-np.transpose(u), b)
                for i in range(rank):
                    #print(np.dot(-np.transpose(u), b), "scalar")
                    assert np.dot(-np.transpose(u), b).shape == (3,1)
                    assert s[i] > 1.0/tau
                    y[i] = utb[i] / s[i]
                    #for other dimensions: keep the old y

                #for i in range(rank,3):
                #    assert np.abs(utb[i]-0.0)
                should_bezero = utb[rank:]
                incorrect = np.any(np.abs(should_bezero - 0.0) > 0.001*2)
                if incorrect:
                    if not SUPPRESS_DEBUG_PRINT:
                        print "should have been zero: ", should_bezero

                new_x = np.dot(np.transpose(v), y)
                #print(old_x.ravel(), " -> ", new_x.ravel())
                #print("    delta=", (new_x - old_x).ravel())
            del(s)
            new_x
            self.new_verts[vi,0:3] = new_x[:,0]
            #self.new_verts[vi,3] = 1

            dd = np.linalg.norm( old_x - new_x)
            if dd > 1.0:
                print("moved far: += ", dd, (new_x-old_x).ravel() )

            assert old_x.shape == (3, 1)
            #alpha = 1
            self.new_verts[vi,0:3] = new_x[:,0] * alpha + old_x[:,0] * (1.0-alpha)

            if not np.all(np.abs(utb.ravel()[rank:] ) < 0.0001):
                #print("s", s.ravel()/s[0], "   utb", utb.ravel()/s[0])
                pass
            self.verts_ranks[vi] = rank

        print("max rank = ", np.max(self.verts_ranks), ",  min rank = ", np.min(self.verts_ranks))
        if not nihl0:
            if not np.min(self.verts_ranks) >= 1:
                print("Warning: assertion: np.min(self.verts_ranks) >= 1 failed." )
        if nihl0:
            print("degenerate: (0,0,0) triangle")

        if False:
            assert np.min(self.verts_ranks) >= 1

    def adaptive_subdivision(self):
        pass

    #def vertex_resampling
    def vertex_resampling(self):
        #assert normals (gradients) are normalised
        c = 2.0
        def kij(i,j):
            assert i != j
            #i,j are centroids
            pi, pj = (self.centroids[i,0:3], self.centroids[j,0:3])

            mi, mj = (self.centroid_normals[i,0:3], self.centroid_normals[j,0:3])
            assert mi.shape == (3,)
            assert mj.shape == (3,)
            assert np.abs(np.linalg.norm(mi) - 1.0) < 0.0000001
            assert np.abs(np.linalg.norm(mj) - 1.0) < 0.0000001
            mimj = np.dot(np.transpose(mi), mj)
            #mimj[mimj>1.0] = 1.0
            if mimj>1.0:
                mimj= 1.0
            pipj = np.linalg.norm(pi - pj)
            print "pipj ", pipj, "  arccos= pi*",np.arccos(mimj)/np.pi #why is it zero??
            assert pipj == np.linalg.norm(pi - pj, ord=2)

            kij = np.arccos(mimj) / pipj  # radians?
            return kij

        def wi(i, ja):
            """ ja = list of centroid indices (face index). i is a face index. """
            assert not i in ja
            ki = 0
            for j in ja:
                ki += kij(i,j)
            wi = 1.0 + c*ki
            print("w_i=", wi)
            return wi
        vi = 1
        ca = self.vertex_neighbours_list[vi]
        print("ca: ", ca)
        for i in ca:
            print("i",i)
            ca_ = filter( lambda idx: idx != i, ca )
            print(i, ca_)
            print( wi(i, ca_) )
        #exit()


if __name__ == '__main__':
    pass
