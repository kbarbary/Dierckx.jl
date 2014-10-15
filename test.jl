import Dierckx

x = [1., 2., 3., 4., 5.]
y = [1., 2., 3., 4., 5.]
z = ones(5, 5)
tx, ty, c, fp = Dierckx.regrid(x, y, z)

println(tx)
