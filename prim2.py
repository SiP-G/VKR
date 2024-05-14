from dolfin import *

mesh1 = Mesh('prim11.xml')
subdomains1 = MeshFunction('size_t', mesh1, 'prim11_physical_region.xml')
boundaries1 = MeshFunction ('size_t', mesh1 , 'prim11_facet_region.xml')


mesh2 = Mesh('prim12.xml')
subdomains2 = MeshFunction('size_t', mesh2, 'prim12_physical_region.xml')
boundaries2 = MeshFunction ('size_t', mesh2 , 'prim12_facet_region.xml')


mesh3 = Mesh('prim13.xml')
subdomains3 = MeshFunction('size_t', mesh3, 'prim13_physical_region.xml')
boundaries3 = MeshFunction ('size_t', mesh3 , 'prim13_facet_region.xml')

 
mesh4 = Mesh('prim14.xml')
subdomains4 = MeshFunction('size_t', mesh4, 'prim14_physical_region.xml')
boundaries4 = MeshFunction ('size_t', mesh4 , 'prim14_facet_region.xml')


meshEt = Mesh('prim15.xml')
subdomainsEt = MeshFunction('size_t', meshEt, 'prim15_physical_region.xml')
boundariesEt = MeshFunction ('size_t', meshEt , 'prim15_facet_region.xml')

mu = 0.8*1.0e11
lmbda = 1.25*1.0e11
f = Constant((0.0,-1.0e8))
g = Constant(0.0)
def epsilon(u):
	return 0.5*( grad(u) + grad(u).T)
def sigma(u):
	return lmbda*div(u)*Identity(2) + 2*mu*epsilon(u)

VEt = VectorFunctionSpace (meshEt , "CG" , 1)

dxEt = Measure('dx' , domain=meshEt , subdomain_data = subdomainsEt)
dsEt = Measure('ds' , domain=meshEt , subdomain_data = boundariesEt)

V1 = VectorFunctionSpace (mesh1 , "CG" , 1)

dx1 = Measure('dx' , domain=mesh1 , subdomain_data = subdomains1)
ds1 = Measure('ds' , domain=mesh1 , subdomain_data = boundaries1)

V2 = VectorFunctionSpace (mesh2 , "CG" , 1)

dx2 = Measure('dx' , domain=mesh2 , subdomain_data = subdomains2)
ds2 = Measure('ds' , domain=mesh2 , subdomain_data = boundaries2)

V3 = VectorFunctionSpace (mesh3 , "CG" , 1)

dx3 = Measure('dx' , domain=mesh3 , subdomain_data = subdomains3)
ds3 = Measure('ds' , domain=mesh3 , subdomain_data = boundaries3)

V4 = VectorFunctionSpace (mesh4 , "CG" , 1)

dx4 = Measure('dx' , domain=mesh4 , subdomain_data = subdomains4)
ds4 = Measure('ds' , domain=mesh4 , subdomain_data = boundaries4)


u0_val = Constant(0.1)
p0 = Constant(1000.0)

u0Et = interpolate(u0_val, VEt)
u01 = interpolate(u0_val, V1)
u02 = interpolate(u0_val, V2)
u03 = interpolate(u0_val, V3)
u04 = interpolate(u0_val, V4)

bc1Et = DirichletBC(VEt.sub (0), g , boundariesEt , 8)
bc11 = DirichletBC(V1.sub (0), g , boundaries1 , 8)
bc12 = DirichletBC(V2.sub (0), g , boundaries2 , 8)
bc13 = DirichletBC(V3.sub (0), g , boundaries3 , 8)
bc14 = DirichletBC(V4.sub (0), g , boundaries4 , 8)

bc2Et = DirichletBC(VEt.sub (1), g , boundariesEt , 11)
bc21 = DirichletBC(V1.sub (1), g , boundaries1 , 11)
bc22 = DirichletBC(V2.sub (1), g , boundaries2 , 11)
bc23 = DirichletBC(V3.sub (1), g , boundaries3 , 11)
bc24 = DirichletBC(V4.sub (1), g , boundaries4 , 11)

bcsEt = [bc1Et, bc2Et]
bcs1 = [bc11, bc21]
bcs2 = [bc12, bc22]
bcs3 = [bc13, bc23]
bcs4 = [bc14, bc24]

k = Constant(10)

T=6.0
N = 100
tau = T/N

info("Solving MeshEt")
uEt = TrialFunction(VEt)
vEt = TestFunction(VEt)
a = (1/tau)*uEt*vEt*dxEt + k*inner(sigma(uEt), epsilon(vEt)) * dxEt
L = f * vEt * dxEt +(1/tau)*u0Et*vEt*dxEt + inner(f, vEt)*dsEt(10) + inner(f, vEt)*dsEt(9)

uEt = Function(VEt)

print(uEt.vector().get_local())
file = File('./results3/time_depEt.pvd')
t=0
while t<T:
	t+=tau
	solve(a == L, uEt, bcsEt)
	file << uEt
	u0Et.assign(uEt)


info("Solving mesh1")
u1 = TrialFunction(V1)
v1 = TestFunction(V1)
a1 = (1/tau)*u1*v1*dx1 + k*inner(sigma(u1), epsilon(v1)) * dx1
L1 = f * v1 * dx1 + (1/tau)*u01*v1*dx1 + inner(f, v1)*ds1(10) + inner(f, v1)*ds1(9)
u1 = Function(V1)
file1 = File('./results3/time_dep1.pvd')
t=0
while t<T:
        t+=tau
        solve(a1 == L1, u1, bcs1)
        file1 << u1
        u01.assign(u1)

info("Solving mesh2")
u2 = TrialFunction(V2)
v2 = TestFunction(V2)
a2 = (1/tau)*u2*v2*dx2 + k*inner(sigma(u2), epsilon(v2)) * dx2
L2 = f * v2 * dx2+ (1/tau)*u02*v2*dx2 + inner(f, v2)*ds2(10) + inner(f, v2)*ds2(9)
u2 = Function(V2)
file2 = File('./results3/time_dep2.pvd')
t=0
while t<T:
        t+=tau
        solve(a2 == L2, u2, bcs2)
        file2 << u2
        u02.assign(u2)

info("Solving mesh3")
u3 = TrialFunction(V3)
v3 = TestFunction(V3)
a3 = (1/tau)*u3*v3*dx3 + k*inner(sigma(u3), epsilon(v3))*dx3
L3 = f * v3 * dx3 + (1/tau)*u03*v3*dx3 + inner(f, v3)*ds3(10) + inner(f, v3)*ds3(9)
u3 = Function(V3)
file3 = File('./results3/time_dep3.pvd')
t=0
while t<T:
        t+=tau
	solve(a3 == L3, u3, bcs3)
	file3 << u3
        u03.assign(u3)

info("Solving mesh4")
u4 = TrialFunction(V4)
v4 = TestFunction(V4)
a4 = (1/tau)*u4*v4*dx4 + k*inner(sigma(u4), epsilon(v4)) * dx4
L4 = f * v4 * dx4 + (1/tau)*u04*v4*dx4 + inner(f, v4)*ds4(10) + inner(f, v4)*ds4(9)
u4 = Function(V4)
file4 = File('./results3/time_dep4.pvd')
t=0
while t<T:
        t+=tau
        solve(a4 == L4, u4, bcs4)
        file4 << u4
        u04.assign(u4)
        
u_e = interpolate(u2, VEt)
E1 = (u_e - uEt)*dx
E1_error = abs(assemble(E1))*100
E2_t = inner(uEt-u_e,uEt-u_e)*dx
E2_b = inner(uEt,uEt)*dx
E2 = sqrt(abs(assemble(E2_t))/abs(assemble(E2_b)))*100
E3_t = inner(grad(uEt-u_e),grad(uEt-u_e))*dx
E3_b = inner(grad(uEt),grad(uEt))*dx
E3 = sqrt(abs(assemble(E3_t))/abs(assemble(E3_b)))*100
print("abs error = " + str(E1_error) + "%\nerror L2 = " + str(E2) + "%\nerror H1 = " + str(E3) + "%")
