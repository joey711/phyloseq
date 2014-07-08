# python/pycogent code to create reference UniFrac results for unit testing.
# Obviously, assumes that pycogent has generated correct values as well.
#
from cogent.maths.unifrac.fast_unifrac import fast_unifrac_file
import numpy
envs_in = open("inst/extdata/gp500test.env.txt")
tree_in = open("inst/extdata/gp500test.tree")
res_uuf = fast_unifrac_file(tree_in, envs_in, weighted=False)
numpy.savetxt("inst/extdata/gp500-uuf.csv", res_uuf['distance_matrix'][0], delimiter=",")
# Now wuf-unnormalized
envs_in = open("inst/extdata/gp500test.env.txt")
tree_in = open("inst/extdata/gp500test.tree")
res_wufu = fast_unifrac_file(tree_in, envs_in, weighted=True)
numpy.savetxt("inst/extdata/gp500-wufu.csv", res_wufu['distance_matrix'][0], delimiter=",")
# Now wuf (normalized... weighted='correct')
envs_in = open("inst/extdata/gp500test.env.txt")
tree_in = open("inst/extdata/gp500test.tree")
res_wuf = fast_unifrac_file(tree_in, envs_in, weighted='correct')
numpy.savetxt("inst/extdata/gp500-wuf.csv", res_wuf['distance_matrix'][0], delimiter=",")